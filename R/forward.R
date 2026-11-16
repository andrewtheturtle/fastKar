# In this file I define some functions used in the prediction of 3D structure at rearrangements.

#' Simulate the Hi-C data of a set of gWalks in reference coordinates using the fastKar framework
#' 
#' @param walks either a gGnome::gWalk or a list containing objects "graph" (a gGnome gGraph), "snode.id" (a list of arrays of node.ids inside "graph" in their walk order), and "circular" (an array of booleans indicating the circularity of the walks)
#' @param target_region a gRange corresponding to the region of interest for the simulation. If left blank, assume the footprint of the walks.
#' @param pix.size resolution of the simulation tiles, default is 100kb. If set to 0, each tile corresponds to half a node in the gGraph.
#' @param if.comps boolean setting whether or not to simulate compartments (defaults to FALSE)
#' @param mc.cores accepted for compatibility (defaults to 1). The walk kernel runs single-threaded in C++, where forking over walks costs more than it saves; results do not depend on it. Parallelize over candidates instead.
#' @param if.sum boolean deciding whether to sum the simulated Hi-C signal from all walks or whether to return a list of Hi-C maps (defaults to TRUE)
#' @param depth Hi-C depth as returned by estimate.depthratio() (defaults to 1); halved here to the per-haploid-copy depth that simulate_walks() and simulate_diploid_map() take
#' @param purity tumor purity/fraction of the sample explained by the supplied walks (defaults to 1). If less than 1, the remaining fraction is simulated as a diploid reference genome over the same context and added to the tumor simulation.
#' @param model sets whether to simulate Hi-C (model = 0) or long-read data (model = length of reads in bp)
#' @param if.interchr boolean, sets whether to simulate interchromosomal contacts (default is TRUE). If if.sum = FALSE, then this is set to FALSE. 
#' @param gm.out boolean, chooses whether to return a gMatrix or a data.table (defaults to TRUE, returning a gMatrtix)
#' @param prepped optional result of prep_for_sim() for this graph, target_region and pix.size, reused instead of tiling again (target_region, pix.size and if.comps are then not used for tiling)
#' @param conserve use the contact-conserving trans term (R/conserve.R): every copy-bp has the same total contacts, so copies that lose cis (circles, short walks, walk ends) get more trans. Default getOption('fastKar.conserve', FALSE)
#' @return either a gMatrix (gm.out = if.sum = T), a list of gMatrices (gm.out = T, if.sum = F), a data.table (gm.out = F, if.sum = T), or a list of data.tables (gm.out = if.sum = F)
forward_simulate <- function(walks,target_region = NULL,pix.size=1e5,if.comps=FALSE,mc.cores=1,if.sum=TRUE,depth=1,purity=1,model=0,if.interchr=T,gm.out=T,prepped=NULL,conserve=getOption('fastKar.conserve',FALSE)){
    if (!is.numeric(purity) || length(purity) != 1 || is.na(purity) || purity < 0 || purity > 1){
        stop('purity must be a single numeric value between 0 and 1')
    }
    if (purity < 1 && !if.sum){
        stop('purity < 1 is currently supported only when if.sum = TRUE')
    }
    haploid.depth = depth/2
    if (is.null(prepped)){
        if (is.null(target_region)){target_region = walks$footprint}
        prepped.data = prep_for_sim(walks,target_region,pix.size,if.comps)
    } else{
        prepped.data = prepped
    }
    if (purity == 1){
        return(simulate_walks(walks,prepped.data$tiled.target,prepped.data$widthdt,if.comps,mc.cores,if.sum,haploid.depth,model,if.interchr,gm.out,conserve=conserve))
    }
    tumor.sim = simulate_walks(walks,prepped.data$tiled.target,prepped.data$widthdt,if.comps,mc.cores,TRUE,haploid.depth*purity,model,if.interchr,F,conserve=conserve)
    normal.sim = simulate_diploid_map(prepped.data$tiled.target,if.comps=if.comps,depth=haploid.depth*(1-purity),model=model,if.interchr=if.interchr,gm.out=F,
                                      conserve=conserve,seqlengths=GenomeInfoDb::seqlengths(walks$graph$gr))
    if (identical(normal.sim$dat$id, tumor.sim$dat$id) && !is.unsorted(tumor.sim$dat$id, strictly = TRUE)){
        # both maps on the same template, in pixel id order: add the value vectors
        combined.dat = normal.sim$dat[,.(i,j,id)]
        combined.dat[,value:=normal.sim$dat$value + tumor.sim$dat$value]
        setattr(combined.dat,'sorted','id')
    } else{
        combined.dat = merge.data.table(normal.sim$dat[,.(i,j,id,value)],
                                        tumor.sim$dat[,.(id,value)],
                                        by='id',all=TRUE,suffixes=c('.normal','.tumor'))
        combined.dat[is.na(value.normal),normal.value:=0]
        combined.dat[is.na(value.tumor),value.tumor:=0]
        combined.dat = combined.dat[,.(i,j,id,value=value.normal+value.tumor)]
    }
    if (gm.out){
        return(gM(gr=dt2gr(prepped.data$tiled.target),dat=combined.dat))
    } else{
        return(list(gr=prepped.data$tiled.target,dat=combined.dat))
    }
}

#' Tile target region and prepare data necessary for Hi-C simulation. 
#'
#' @param walks either a gGnome::gWalk or a list containing objects "graph" (a gGnome gGraph), "snode.id" (a list of arrays of node.ids inside "graph" in their walk order), and "circular" (an array of booleans indicating the circularity of the walks)
#' @param target_region a gRange corresponding to the region of interest for the simulation. If left blank, assume the footprint of the walks.
#' @param pix.size resolution of the simulation tiles, default is 100kb. If set to 0, each tile corresponds to half a node in the gGraph.
#' @param if.comps boolean setting whether or not to simulate compartments (defaults to FALSE)
#' @return a list with two elements: "tiled.target", a data.table containing the tiles, their coordinates, and compartment identity, and "widthdt", a data.table containing the widths of all nodes from the input gGraph
prep_for_sim = function(walks,target_region=NULL,pix.size=1e5,if.comps=F){
    #first, get nodes from graph
    if (is.null(walks$graph) | is.null(walks$graph$gr$cn)){
        stop('Please supply gGraph with CN in input walks object')
    }
    nodesgr = walks$graph$gr[,c('node.id','cn')] 
    nodesgr = nodesgr[strand(nodesgr)=='+']
    widthdt = gr2dt(nodesgr)[,.(nodewidth=width,node.id,start,end)] #start/end: walk_rows() needs them for nodes the target only partly covers
    #then, subset to target region or make target region if not provided
    if (is.null(target_region)){
        if (is.null(walks$footprint)){
            allnodes = unique(abs(do.call(c,walks$snode.id)))
            target_region = nodesgr[allnodes]
        } else{
            target_region = walks$footprint
        }
    }
    #now, tile target region, merging with node boundaries where the tumor graph is defined.
    #If target_region extends outside the supplied graph, keep those extra tiles with
    #node.id = NA and tumor cn = 0 so a diploid normal component can still be simulated
    #on the full requested region without rebinning.
    # parts of the target not covered by any node, within the graph's seqlengths (the reference took
    # these from the no-query rows of gr.merge(streduce(nodes %&% target), target, all = TRUE))
    bare = function(g) GenomicRanges::GRanges(as.character(GenomicRanges::seqnames(g)), IRanges::ranges(g))
    target.ns = bare(target_region)
    nodes.ns = bare(nodesgr)
    any.in.target = any(IRanges::overlapsAny(nodes.ns, target.ns))
    gaps.ns = suppressWarnings(GenomicRanges::setdiff(target.ns, nodes.ns))
    sl = GenomeInfoDb::seqlengths(nodesgr)[as.character(GenomicRanges::seqnames(gaps.ns))]
    sl[is.na(sl)] = Inf
    keep = GenomicRanges::start(gaps.ns) <= sl
    gaps.ns = gaps.ns[keep]
    sl = sl[keep]
    over = GenomicRanges::end(gaps.ns) > sl
    if (any(over)) GenomicRanges::end(gaps.ns)[over] = as.integer(sl[over])
    outside.target = streduce(gaps.ns)
    tiled.parts = list()
    if (any.in.target){
        tiled.parts[[length(tiled.parts)+1]] = gr2dt(gr.merge(gr.tile(target_region,pix.size),nodesgr))[,.(start,end,seqnames,width,node.id,cn)]
    }
    if (length(outside.target)){
        outside.dt = gr2dt(gr.tile(outside.target,pix.size))[,.(start,end,seqnames,width)]
        outside.dt[,`:=`(node.id=NA_integer_,cn=0)]
        tiled.parts[[length(tiled.parts)+1]] = outside.dt[,.(start,end,seqnames,width,node.id,cn)]
    }
    tiled.target = rbindlist(tiled.parts,use.names=TRUE,fill=TRUE)
    setorder(tiled.target,seqnames,start,end,node.id)
    tiled.target[,tile.id:=.I]
    tiled.target = tiled.target[,.(start,end,seqnames,tile.id,width,node.id,cn)]
    #check for compartment data
    if (if.comps==FALSE){
        comps.gr = NULL
        tiled.target$compartment = 'A'
    }else{
        stop('Compartments not supported in this version')
        comps.gr = gr.nochr(fastKar::compartment_lookup)
        tiled.target = eval_comps(gr.tile(target_region,pix.size),comps.gr) #need to update eval_comps to dt mode
    }
    return(list(tiled.target=tiled.target,widthdt=widthdt))
}

#' Double the pixels between two distinct tiles at the same locus (in place)
#'
#' On allelic graphs the h1/h2 copies of a node share coordinates, so the tiling has one tile per haplotype
#' at each locus. The kernel stores contacts between tiles p != q once, on the upper-triangle pixel (p, q),
#' and contacts within one tile on (p, p) with both orders counted. When p and q are summed into one bin
#' (gM() disjoining, or .tiles_to_bins()), that bin's diagonal needs (p, q) twice, as (p, q) and (q, p).
#' Graphs without shared tile coordinates are left untouched.
#'
#' @param dat data.table with columns i, j (tile indices) and value, modified by reference
#' @param tiled.target tiling from prep_for_sim()
#' @return dat, invisibly
double_samelocus_pixels <- function(dat, tiled.target){
    key = paste(as.character(tiled.target$seqnames), tiled.target$start, tiled.target$end)
    if (!anyDuplicated(key)) return(invisible(dat))
    locus = match(key, key)
    idx = which(dat$i != dat$j & locus[dat$i] == locus[dat$j])
    if (length(idx)) set(dat, i = idx, j = 'value', value = 2 * dat$value[idx])
    invisible(dat)
}

#' Downstream of prep_for_sim, calculate Hi-C data in reference coordinates given sets of walks on a graph and a tiled target region
#' 
#' @param walks either a gGnome::gWalk or a list containing objects "graph" (a gGnome gGraph), "snode.id" (a list of arrays of node.ids inside "graph" in their walk order), and "circular" (an array of booleans indicating the circularity of the walks)
#' @param tiled.target a data.table produced by prep_for_sim(), containing target tile coordinates, tile.ids, and corresponding node.ids, as well as CN and width
#' @param widthdt a data.table produced by prep_for_sim(), containing width data for all the nodes in the gGraph simulated here.
#' @param if.comps boolean setting whether or not to simulate compartments (defaults to FALSE)
#' @param mc.cores accepted for compatibility (defaults to 1). The walk kernel runs single-threaded in C++, where forking over walks costs more than it saves; results do not depend on it. Parallelize over candidates instead.
#' @param if.sum boolean deciding whether to sum the simulated Hi-C signal from all walks or whether to return a list of Hi-C maps (defaults to TRUE)
#' @param depth per-haploid-copy depth, i.e. estimate.depthratio() / 2 (forward_simulate() halves its depth before calling); defaults to 1
#' @param model sets whether to simulate Hi-C (model = 0) or long-read data (model = length of reads in bp)
#' @param if.interchr boolean, sets whether to simulate interchromosomal contacts (default is TRUE). If if.sum = FALSE, then this is set to FALSE. 
#' @param gm.out boolean, chooses whether to return a gMatrix or a data.table (defaults to TRUE, returning a gMatrtix)
#' @param conserve contact-conserving trans term (trans_weights() replaces the tile CN in the trans term); default getOption('fastKar.conserve', FALSE)
#' @return either a gMatrix (gm.out = if.sum = T), a list of gMatrices (gm.out = T, if.sum = F), a data.table (gm.out = F, if.sum = T), or a list of data.tables (gm.out = if.sum = F)
simulate_walks <- function(walks,tiled.target,widthdt,if.comps=F,mc.cores=1,if.sum=T,depth=1,model=0,if.interchr=T,gm.out=T,conserve=getOption('fastKar.conserve',FALSE)){
    # walks is either a gWalk object or a list with elements: snode.id, circular
    # tiled.target is a gr2dt with metadata columns of tile.id, node.id, and cn
    # mc.cores is accepted for compatibility; the C++ kernel is single-threaded and results do not depend on it.
    if (if.comps){
        stop('Compartments not supported in this version')
    }
    snode.id = walks$snode.id
    # subset to walks that hit target region
    hit_target = walks_hit_target(snode.id,tiled.target$node.id)
    if (!any(hit_target)){
        stop('none of the walks hits the target region')
    }
    circular = walks$circular[hit_target]
    #then, check for walk CN and if not provided, assume all 1s
    walk.cn.all = walks$dt$cn
    if(is.null(walk.cn.all)){
        walk.cn = rep(1,sum(hit_target))
    } else{
        walk.cn = walk.cn.all[hit_target]
    }
    # every walk as one run of local rows (tiles in walk order, flipped on the minus strand)
    rows = walk_rows(snode.id[hit_target],tiled.target,widthdt)
    kern = forward_kernel_args(model)
    n.tiles = nrow(tiled.target)
    walk_acc = function(w, cn){
        # local rows of walks w (consecutive), walk CNs cn
        r = seq.int(rows$walk_ptr[w[1]] + 1, length.out = rows$walk_ptr[w[length(w)] + 1] - rows$walk_ptr[w[1]])
        ptr = as.integer(rows$walk_ptr[c(w, w[length(w)] + 1)] - rows$walk_ptr[w[1]])
        forward_walks_cpp(ptr, rows$orig[r], rows$width[r], rows$int_width, as.logical(circular[w]), as.numeric(cn), n.tiles,
                          kern$z, kern$lmin, kern$slope, kern$maxval, kern$maxval_density, kern$model, kern$lr_denom, kern$diagf)
    }
    if(if.sum){
        acc = walk_acc(seq_along(walk.cn), walk.cn)
        interchr = if.interchr & model==0
        trans.cn = if (interchr && isTRUE(conserve)) trans_weights(walks, tiled.target, widthdt)$B else tiled.target$cn
        total.counts = setDT(forward_assemble_cpp(as.numeric(tiled.target$width), as.numeric(trans.cn), acc, depth,
                                                  fastKar::interchrom_density, interchr))
        if (interchr){
            setattr(total.counts,'sorted','id') # downstream code expects the summed map keyed by id
        }
        double_samelocus_pixels(total.counts, tiled.target)
        if (gm.out){
            return(gM(gr = dt2gr(tiled.target),dat=total.counts))
        } else{
            return(list(gr = tiled.target, dat=total.counts))
        }
    }else{
        target.dat = as.data.table(make_template_dat_cpp(tiled.target))
        setattr(target.dat,'sorted',c('i','j'))
        if (gm.out){
            target.gr = dt2gr(tiled.target)
        }
        scaled_out = lapply(seq_along(walk.cn),function(ii){
            x = data.table::copy(target.dat)
            set(x,j='value',value=depth*(walk.cn[ii]*walk_acc(ii, 1)))
            double_samelocus_pixels(x, tiled.target)
            if (gm.out){
                return(gM(gr=target.gr,dat=x))
            } else{
                return(x)
            }
        })
        return(scaled_out)
    }
}

#' Which walks visit at least one node of the target tiling
#'
#' @param snode.id list of signed node id vectors
#' @param target.nodes node.id column of the tiled target (may contain NA)
#' @return logical vector, one per walk
walks_hit_target = function(snode.id,target.nodes){
    lens = base::lengths(snode.id)
    if (!length(lens)) return(logical(0))
    flat = abs(unlist(snode.id,use.names=FALSE))
    hit = flat %in% target.nodes[!is.na(target.nodes)]
    tabulate(rep.int(seq_along(lens), lens)[hit], nbins = length(lens)) > 0
}

#' Local rows (tiles in walk order) of a set of walks, as simulate_walks() lays them out
#'
#' Each node visit expands into the node's rows in genomic order (reversed on the minus strand): its tiles in
#' the target tiling plus one untiled row for every part of the node the target does not cover (before, between
#' and after its tiles), or one untiled row with the node's full width if the node has no tiles. Without the
#' untiled parts, tiles on either side of a gap in the target would sit next to each other along the walk.
#' @param widthdt node.id, nodewidth, start, end of every node (prep_for_sim())
#' @return list(walk_ptr = 0-based row offsets per walk (length W+1), orig = tile.id or NA,
#'   width = row width (numeric), int_width = whether the reference width column is integer)
walk_rows = function(snode.id,tiled.target,widthdt){
    if (!all(c('start','end') %in% names(widthdt))){
        stop('widthdt needs the node start and end columns (see prep_for_sim())')
    }
    lens = base::lengths(snode.id)
    flat = unlist(snode.id,use.names=FALSE)
    v = abs(flat)
    # per tiled node, its rows in genomic order: tiles (orig = tile.id) and the uncovered parts (orig = NA)
    tt = data.table(node.id = tiled.target$node.id, orig = as.integer(tiled.target$tile.id), start = as.numeric(tiled.target$start),
                    end = as.numeric(tiled.target$end), width = as.numeric(tiled.target$width), r = seq_len(nrow(tiled.target)))[!is.na(node.id)]
    setorder(tt, node.id, start, r)
    nd = match(tt$node.id, widthdt$node.id)
    tt[, `:=`(nstart = as.numeric(widthdt$start[nd]), nend = as.numeric(widthdt$end[nd]))]
    tt[, prev_end := shift(end, fill = nstart[1] - 1), by = node.id]
    gaps = rbind(tt[start - prev_end > 1, .(node.id, pos = start - 0.5, orig = NA_integer_, width = start - prev_end - 1)],
                 tt[tt[, .I[.N], by = node.id]$V1][nend > end, .(node.id, pos = end + 0.5, orig = NA_integer_, width = nend - end)])
    tpl = rbind(tt[, .(node.id, pos = start, orig, width)], gaps)
    setorder(tpl, node.id, pos)
    rl = rle(tpl$node.id)
    rl_first = cumsum(c(1L, rl$lengths))[seq_along(rl$values)]
    mt = match(v, rl$values)
    ntile = ifelse(is.na(mt), 0L, rl$lengths[mt])
    nrow_visit = pmax(ntile, 1L)
    vis = rep.int(seq_along(v), nrow_visit)
    k = sequence(nrow_visit)
    minus = flat[vis] < 0
    k[minus] = nrow_visit[vis[minus]] - k[minus] + 1L
    tiled = ntile[vis] > 0
    row = rep(NA_integer_, length(vis))
    row[tiled] = rl_first[mt[vis[tiled]]] + k[tiled] - 1L
    orig = tpl$orig[row]
    width = tpl$width[row]
    untiled = !tiled
    width[untiled] = as.numeric(widthdt$nodewidth[match(v[vis[untiled]], widthdt$node.id)])
    walk_ptr = c(0, cumsum(as.numeric(nrow_visit)))[cumsum(c(0, lens)) + 1]
    list(walk_ptr = walk_ptr,
         orig = orig,
         width = width,
         int_width = is.integer(tiled.target$width))
}

#' Model constants for the C++ walk kernel, read from the package data at call time
forward_kernel_args = function(model){
    f = fastKar::splineobj$distance_spline
    z = environment(f)$z
    lmin = z$x[1]
    list(z = z, lmin = lmin, slope = f(lmin, deriv = 1), maxval = 1e8,
         maxval_density = distance_density(1e8), model = model,
         lr_denom = if (model > 0) 10*min(fastKar::small_lookup[d>0]$d)^2 else 1,
         diagf = fastKar::splineobj$diag_spline)
}

#' Pixel id of reference pixel (i, j), i <= j, on an n-tile template: its row in make_template_dat()
#'
#' @param i,j integer tile indices (1-based) with i <= j
#' @param n number of tiles
#' @return integer pixel ids
pixel_id = function(i, j, n){
    i = as.numeric(i); j = as.numeric(j)
    as.integer((i - 1)*n - (i - 1)*(i - 2)/2 + (j - i + 1))
}

#' Hi-C contact density per bp^2 at genomic distance dist (bp), from the Rao 2014 distance spline
#'
#' The spline is trained on distances of 10 kb to 249 Mb. Below 10 kb its cubic (fmm) extrapolation
#' diverges as dist -> 0 (~1e10 at 2 bp), which walks that revisit tiny nodes back to back can reach.
#' Below the first knot we instead extend it as a power law with the slope at that knot (~ -1).
#' @param dist numeric vector of distances in bp (> 0)
#' @return numeric vector of contact densities
distance_density = function(dist){
    ldist = log(dist)
    lmin = environment(fastKar::splineobj$distance_spline)$z$x[1]
    out = fastKar::splineobj$distance_spline(pmax(ldist, lmin))
    below = ldist < lmin
    out[below] = out[below] + fastKar::splineobj$distance_spline(lmin, deriv = 1) * (ldist[below] - lmin)
    exp(out)
}

#' Simulate Hi-C data using spline trained on Rao 2014 data
#' 
#' @param walk.dt data.table prepared by simulate_walks(), each row is a tile to simulate and columns are start,end,midpoints of the tile in local walk coordinates.
#' @param is.circular boolean indicating whether the walk is circular or not (defaults to FALSE)
#' @param model sets whether to simulate Hi-C (model = 0) or long-read data (model = length of reads in bp)
#' @return data.table of long-format gMatrix with columns i,j,id,value
simulate_map<-  function(walk.dt,is.circular=FALSE,model=0){
    dat.new = as.data.table(make_template_dat_cpp(walk.dt))
    totl = 0 # circle length for the folded pixel integral (0 = linear)
    if (is.circular){
        totl = walk.dt$end[nrow(walk.dt)]
        dat.new[,dist_temp:= abs(walk.dt$mid[j]-walk.dt$mid[i])]
        dat.new[,dist:=pmin(dist_temp,totl-dist_temp)][,dist_temp:=NULL]
    }else{
        dat.new[,dist:= abs(walk.dt$mid[j]-walk.dt$mid[i])]
    }
    if(model==0){
        ka = forward_kernel_args(0) # pixel-integrated kernel (as the C++ kernels), 1/distance beyond 1e8
        dat.new[dist>0,value:=forward_pixel_cpp(walk.dt$width[i], walk.dt$width[j], dist, ka$z, ka$lmin, ka$slope, ka$maxval, ka$maxval_density, totl)]
        dat.new[dist==0,value:=exp(fastKar::splineobj$diag_spline(log(widthprod)))]
   }else {
        min.res = min(fastKar::small_lookup[d>0]$d)
        dat.new[,value:=widthprod/(10*min.res^2)*(exp(-dist/model))]
        #some analytical model for short/long reads with readlength as an input parameter
   }
   dat.new[is.na(value),value:=0]
#   return(dat.new[,.(i,j,value,id)])
   return(dat.new[,.(i,j,value,id,dist)])
}

#' Evaluate compartment identity of all tiles
#' 
#' TO BE UPDATED, CURRENTLY DEPRECATED
eval_comps <- function(tiled.walk,comps.gr=NULL){
    nbin = length(tiled.walk)
    if (is.null(comps.gr)){
        afrac = rep(1,nbin)
        bfrac = rep(0,nbin)
    } else if(class(comps.gr)=='GRanges'){
        afrac = tiled.walk %O% (comps.gr %Q% (!is.na(compartment) & compartment=='A'))
        bfrac = tiled.walk %O% (comps.gr %Q% (!is.na(compartment) & compartment=='B'))
    } else if(comps.gr=='rand'){
        afrac = sample(c(0,1),nbin,replace=TRUE)
        bfrac = 1-afrac
    }
    compartment = ifelse(afrac > 0.5 & bfrac < 0.5, 'A', ifelse(afrac < 0.5 & bfrac > 0.5, 'B', NA))
    tiled.walk$compartment = compartment
    return(tiled.walk)
}

#' Given a set of bins, make the gMatrix data.table corresponding to all pairwise contacts for those bins
#' 
#' @param target.bins either a GRanges, or a data.table. If DT, must have a column "width", giving the width of all the bins
#' @return data.table with columns i,j,id,widthprod corresponding to all pairwise contacts, widthprod is the product of widths of tiles i and j
make_template_dat = function(target.bins,if.comps=FALSE){
	if (inherits(target.bins,'GRanges')){
		return(as.data.table(make_template_dat_cpp(gr2dt(target.bins))))
	}
	else if (inherits(target.bins,'data.table')){
		return(as.data.table(make_template_dat_cpp(target.bins)))
	}
}

# make_template_dat_cpp() is compiled with the package (src/forward_kernel.cpp).

#' Symmetrize a matrix
#' 
#' @param input.mat some input matrix, must be square
#' @return matrix of the same shape as input.mat, with same diagonal, and non-diagonal being sum of matrix and its transpose
#' NOTE: could be sped up with cpp
symmetrize <- function(input.mat){
    output.mat = input.mat + Matrix::t(input.mat)
    Matrix::diag(output.mat) = Matrix::diag(input.mat)
    return(output.mat)
}

#' Directly simulate a diploid reference Hi-C map on an existing tiling
#'
#' This is the fast/reference-specialized equivalent of simulating two linear
#' reference homologs through make_diploid_walks() + simulate_walks(), but it
#' takes the already-prepared tile table from prep_for_sim() and computes the
#' diploid signal directly in reference coordinates.  Using the tumor simulation
#' tiling makes the normal and tumor data.tables directly combinable by pixel id.
#'
#' @param tiled.target data.table produced by prep_for_sim(), containing at least start, end, seqnames, width, and tile.id columns
#' @param if.comps boolean setting whether or not to simulate compartments (currently unsupported)
#' @param depth per-haploid-copy depth, i.e. estimate.depthratio() / 2, as in simulate_walks()
#' @param model 0 (Hi-C); long-read models (model > 0) are not supported here, so forward_simulate() with purity < 1 needs model = 0
#' @param if.interchr boolean, whether to add the same interchromosomal/background term as simulate_walks()
#' @param gm.out boolean, chooses whether to return a gMatrix or a list with gr/dat data.tables
#' @param conserve contact-conserving trans term (reference_trans_weights() in place of CN 2); default getOption('fastKar.conserve', FALSE)
#' @param seqlengths optional chromosome lengths for the reference positions under conserve (else hg38 calibration lengths)
#' @return either a gMatrix or a list with tiled target and contact data.table
simulate_diploid_map = function(tiled.target,if.comps=FALSE,depth=1,model=0,if.interchr=TRUE,gm.out=TRUE,conserve=getOption('fastKar.conserve',FALSE),seqlengths=NULL){
    if (missing(tiled.target) || is.null(tiled.target)){
        stop('tiled.target must be supplied for direct diploid simulation')
    }
    if (if.comps){
        stop('Compartments not supported in direct diploid simulation')
    }
    tiled.target = data.table::copy(as.data.table(tiled.target))
    required.cols = c('start','end','seqnames','width')
    if (!all(required.cols %in% colnames(tiled.target))){
        stop('tiled.target must contain start, end, seqnames, and width columns')
    }
    if (!'tile.id' %in% colnames(tiled.target)){
        tiled.target[,tile.id:=.I]
    }
    tiled.target[,cn:=2]
    tiled.target[,compartment:='A']

    if (model == 0){
        # C++: same-chromosome Hi-C model x 2 (cis), then the trans term and depth as in simulate_walks(): per-tile
        # trans weight 2 (CN), or the reference weights of the conserved trans term
        chrom = match(as.character(tiled.target$seqnames), unique(as.character(tiled.target$seqnames)))
        chrom[is.na(tiled.target$seqnames)] = NA_integer_
        kern = forward_kernel_args(0)
        cis = forward_diploid_cpp(as.integer(chrom), as.numeric((tiled.target$start+tiled.target$end)/2), as.numeric(tiled.target$width),
                                  kern$z, kern$lmin, kern$slope, kern$maxval, kern$maxval_density, kern$diagf)
        trans.cn = if (if.interchr && isTRUE(conserve)) reference_trans_weights(tiled.target, seqlengths) else tiled.target$cn
        out.dat = setDT(forward_assemble_cpp(as.numeric(tiled.target$width), as.numeric(trans.cn), cis, depth,
                                             fastKar::interchrom_density, if.interchr))
    } else{
        # long-read diploid path removed 2026-11-16 (user; it errored at the same_chrom lookup); restore if needed
        stop('simulate_diploid_map: only the Hi-C model (model = 0) is supported; long-read diploid simulation (purity < 1 with model > 0) is not implemented')
    }

    if (gm.out){
        return(gM(gr=dt2gr(tiled.target),dat=out.dat))
    } else{
        return(list(gr=tiled.target,dat=out.dat))
    }
}
