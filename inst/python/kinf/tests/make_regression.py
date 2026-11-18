"""Freeze the selftest numbers -> regression.json (read by test_problem.py and test_optimize.py).

History: first frozen from the OLD modules (model, delta, walks, optimize, infer) at fastKar commit 3798ff7; problem.py and
optimize.py were then shown to reproduce those numbers exactly (rtol 1e-9 / 1e-10), last checked 2026-11-13 with the old
interchrom_density. Refrozen 2026-11-13 from the current modules after interchrom_density became the genome-wide trans density
(3.01e-13, was chr1 x chr2 2.80e-13) and the selftest was regenerated (scripts/make_kinf_selftest.R).
Refrozen 2026-11-16 after the pixel-integrated kernel (with the circle fold) and the 2026-11-16 parameter refit; selftest
regenerated (now with conserve.csv). Frozen with the worktree optimize.py (site-mixture changes, uncommitted then). Same wirings, moves and
search settings as the tests.
Refrozen 2026-11-18 after the retraining on the reference mask (depth = valid-pair coverage; D 4.918653e-13) and the switch
to diag = False (diagonal pixels not scored) by default; selftest regenerated. Run from a clean checkout of the kinf modules (uncommitted optimize.py changes would be frozen in).
Usage: python -I tests/make_regression.py"""
import os, sys, json
import numpy as np
import jax
jax.config.update('jax_enable_x64', True)
here = os.path.dirname(os.path.abspath(__file__)); kin = os.path.dirname(here)
sys.path.insert(0, kin)
import problem as pr, optimize as op

folder = os.path.join(kin, 'selftest')
old = json.load(open(os.path.join(here, 'regression.json')))
g = pr.Graph.load(folder); h = pr.HiC.load(folder, g)
sn, circ = pr.read_gwalks(folder + '/truth_gwalks.json')
out = dict(truth_nll=float(pr.nll(g, h, pr.b_from_perms(g, pr.walks_to_perms(g, sn, circ)))))
bs = pr.random_b(g, jax.random.PRNGKey(0), 10)
kars = []
for b in bs:
    w, c = pr.perms_to_walks(g, pr.perms_from_b(g, b))
    kars.append(sorted([list(map(int, k[1])), bool(k[0])] for k in (pr.canonical(a, cc) for a, cc in zip(w, c))))
out.update(b=[[int(x) for x in b] for b in bs], nll=[float(pr.nll(g, h, b)) for b in bs], canonical=kars)
# moves: the old list (np.random.RandomState(1).randint(0, len(moves), 300)); incremental NLL along it from wiring 0
mv = [int(x) for x in np.random.RandomState(1).randint(0, len(g.moves), 300)]
T = pr.nll_term(g, h); upd = jax.jit(T.update)
b = bs[0]; ko = pr.kar(g, b); aux = T.init(ko); init_nll = float(T.value(aux)); seq = []
for k in mv:
    i, j = g.moves[k]; b = b.at[i].set(b[j]).at[j].set(b[i]); kn = pr.kar(g, b)
    aux = upd(aux, ko, kn, i, j); ko = kn; seq.append(float(T.value(aux)))
out.update(moves=mv, move_nll=seq, init_nll=init_nll)
out['rescale_ratio'] = float(pr.rescale_depth(g, pr.HiC.load(folder, g), 64, 0))
# search: as test_optimize.py section 1
r = op.optimize(g, [T], T=np.geomspace(1, 1000, 4), chains=4, steps=200, block=50, seed=1, topk=3, log=lambda s: None)
out['search_selftest'] = dict(note='optimize.py on selftest (refrozen 2026-11-16, pixel-integrated kernel + refit parameters): T geomspace(1, 1000, 4), '
                                   'chains 4, steps 200, block 50, seed 1, topk 3',
                              best_nll=float(r['best']['E']), best_b=[int(x) for x in r['best']['b']],
                              topk_nll=[float(t['E']) for t in pr.distinct_karyotypes(g, r['top'], 3)])
assert out['moves'] == old['moves'] and out['b'] == old['b'], 'wirings or moves differ from the previous freeze'
json.dump(out, open(os.path.join(here, 'regression.json'), 'w'))
print('truth', out['truth_nll'], '(was', old['truth_nll'], ') ratio', out['rescale_ratio'], '(was', old['rescale_ratio'], ')',
      'search best', out['search_selftest']['best_nll'], 'best_b same as before', out['search_selftest']['best_b'] == old['search_selftest']['best_b'])
