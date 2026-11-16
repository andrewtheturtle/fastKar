## gimme_bfb() tiles chromosome 1 with parse.gr('1'), which takes the chromosome length from
## DEFAULT_GENOME / DEFAULT_BSGENOME and stops when neither is set. For the calling test only,
## point DEFAULT_GENOME at a temporary chrom.sizes file (chromosome 1, hg19 length). A genome
## that is already set is left as it is.
local_default_genome <- function(env = parent.frame()) {
  if (nzchar(Sys.getenv("DEFAULT_GENOME")) || nzchar(Sys.getenv("DEFAULT_BSGENOME"))) {
    return(invisible(NULL))
  }
  f <- withr::local_tempfile(fileext = ".chrom.sizes", .local_envir = env)
  writeLines("1\t249250621", f)
  withr::local_envvar(DEFAULT_GENOME = f, .local_envir = env)
  invisible(f)
}
