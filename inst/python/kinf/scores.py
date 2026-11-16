"""Plain scores f(kar) -> scalar (jittable; static arrays precomputed outside JAX) and energy helpers for optimize.py."""
import numpy as np
import jax
import jax.numpy as jnp


def walk_sum(kar, per_end_values):
    """sum of per-end values over the ends of each walk: array of length 2*nu indexed by walk id (tail state; unused ids are 0)"""
    v = jnp.asarray(per_end_values)
    return jax.ops.segment_sum(v, kar.wid, num_segments=len(v))

def n_walks(g):
    N = 2 * g.nu
    return lambda kar: jnp.zeros(N).at[kar.wid].set(1.0).sum()

def n_circular(g):
    N = 2 * g.nu
    return lambda kar: jnp.zeros(N).at[kar.wid].max(kar.iscyc.astype(float)).sum()

def amp_entropy(g, amp_nodes):
    """entropy of the distribution of amplicon node visits over walks (a visit = left + right half copy, 1/2 each)"""
    N = 2 * g.nu
    w = np.zeros(N); w[0::2] = np.where(np.isin(g.node_of_uid, np.asarray(amp_nodes)), 0.5, 0.0)  # one end per half copy
    w = jnp.asarray(w)
    def f(kar):
        c = walk_sum(kar, w)
        p = c / c.sum()
        return -jnp.sum(jnp.where(p > 0, p * jnp.log(jnp.where(p > 0, p, 1.0)), 0.0))
    return f

def amp_site_fraction(g, amp_nodes, site_nodes):
    """share of the amplicon on the same walk as an integration site (blog amp_site_fraction): copy-weighted amplicon bp on
    walks that visit a site node / amplicon bp on all walks. Site nodes that are also amplicon nodes are dropped. A half
    copy counts its own width, so a node visit counts the node width. site_nodes: graph node ids, e.g. the nodes overlapping
    a locus called by fastKar_integration_tiles()."""
    N = 2 * g.nu
    amp = np.isin(g.node_of_uid, np.asarray(amp_nodes))
    site = np.isin(g.node_of_uid, np.setdiff1d(np.asarray(site_nodes), np.asarray(amp_nodes)))
    assert amp.any(), 'no amplicon node in the graph'
    assert site.any(), 'no site node (outside the amplicon) in the graph'
    wa = np.zeros(N); wa[0::2] = np.where(amp, g.width_uid, 0.0)  # one end per half copy
    ws = np.zeros(N); ws[0::2] = site.astype(float)
    wa, ws = jnp.asarray(wa), jnp.asarray(ws)
    def f(kar):
        a = walk_sum(kar, wa); s = walk_sum(kar, ws)
        return jnp.sum(jnp.where(s > 0, a, 0.0)) / a.sum()
    return f

def schedule(ts, vs):
    """piecewise-linear t -> value, constant outside [ts[0], ts[-1]]"""
    ts, vs = jnp.asarray(ts, float), jnp.asarray(vs, float)
    return lambda t: jnp.interp(t, ts, vs)

def target(k, sched):
    """energy(values, t) = |values[k] - sched(t)|"""
    return lambda values, t: jnp.abs(values[k] - sched(t))
