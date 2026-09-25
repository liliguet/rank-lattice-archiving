# Rank-Lattice Pareto Archiving

Implementation of the rank lattice, of the partition-relaxed algorithm and of
the baselines from

> **Rank-Lattice Pareto Archiving for Exact Invariance under Monotone Recalibration of the Objectives**
> Submitted to *Evolutionary Computation*.


## Installation

Python 3.10 or later is required.

```bash
git clone
cd rank-lattice-archiving
pip install -r requirements.txt
```

## Quick start

```python
import numpy as np
from pdom.archive import QuotientArchive
from pdom.partitions import RankGridPartition

F = np.random.default_rng(0).random((1000, 5))        # objective vectors of a stream
archive = QuotientArchive(RankGridPartition(b=3).fit(F), mode="S")
for i, f in enumerate(F):
    archive.submit(f, data=i, ident=i)                # ident resolves exact duplicates
print(len(archive.reps), "points retained")
```

The identifier passed as `ident` decides between two points with the same
objective vector (rule U0 of Definition 3), so the archive does not depend on
the order in which parallel evaluations finish. Any mutually comparable values
work, such as the position in the stream.

The collector attaches the archive to a pymoo algorithm as its callback.

```python
from pymoo.optimize import minimize
from pdom.problems import get_problem
from pdom.workflow import RankLatticeCollector, as_pymoo_problem, make_baseline

prob = get_problem("maf7", 5)
collector = RankLatticeCollector(b=3)
minimize(as_pymoo_problem(prob), make_baseline("nsga3", 5, 212, seed=0),
         ("n_evals", 20000), seed=0, callback=collector)
A = collector.members()                                # retained objective vectors
F_show, X_show = collector.reduce(212)                 # display set of 212 points
```

## Default configuration

The configuration recommended in the paper is the collector with resolution
`b=3`, cell-order pruning on, and the reference set refreshed with the host's
population every 20 generations, which are the defaults of
`RankLatticeCollector`. Every point the host reported keeps a member within
one cell of the lattice, δ = 2⁻³ in quantile mass per objective, measured
against the reference set that was current when the point was decided.
Proposition 3 of the paper accounts for the later refreshes pair by pair.
`reduce(k)` selects a display set of `k` members by distance-based subset
selection in quantile coordinates, which is invariant along with the archive.

With `prune=False` the archive also keeps strict elitism, so that every point
the host reported is weakly dominated by a member, and it is correspondingly
larger. A finer resolution such as `b=5` or `b=6` tightens the coverage level
to 2⁻⁵ or 2⁻⁶ at the cost of a larger archive.

## Recalibrations evaluated in floating point

The invariance holds for maps that are strictly increasing on the values that
actually occur. A map evaluated in binary64 can fail this, because a
compressive map such as `log` or `log1p` can send two distinct values to the
same number, and so can the rounding of a reciprocal. The module
`pdom.recalibration` checks a recorded stream for this and repairs it.
In the example below the two values are the ones of Remark 7 of the paper.

```python
import numpy as np
from pdom.recalibration import order_violations, order_preserving

F = np.array([[6575303.1262349002], [6575303.1262349030]])
G = np.log(F)
print(order_violations(F, G))        # [1], the two logarithms coincide
H = order_preserving(F, G)           # strictly increasing again, same order as (log v, v)
print(order_violations(F, H))        # [0]
```

## Examples

The folder `examples/` contains three complete scripts, which run in a few
minutes each.

| Script | What it shows |
|---|---|
| `invariance_demo.py` | The rank lattice keeps the same set after a cube and `log(1+v)` recalibration, an affine grid does not |
| `collector_demo.py` | The rank lattice attached to NSGA-III as a passive collector, with its display set, at the budget of the paper |
| `run_prd_and_baselines.py` | The partition-relaxed algorithm and two baselines on MaF2 with five objectives |

Run them from the repository root, for example `PYTHONPATH=. python examples/invariance_demo.py`.

## Package contents

| Module | Content |
|---|---|
| `pdom/partitions.py` | Rank lattice (`RankGridPartition`), affine grid, quantile Mondrian, nearest-anchor and Hilbert-curve partitions |
| `pdom/archive.py` | Conservative quotient update of Definition 3 and Algorithm 1, in Mode S and Mode τ, with cell-order pruning |
| `pdom/workflow.py` | `RankLatticeCollector`, `run_prd`, `run_baseline`, `make_baseline` and helpers |
| `pdom/display.py` | Distance-based subset selection in quantile or normalized coordinates |
| `pdom/recalibration.py` | Arithmetic check and order-preserving repair of a recalibrated stream |
| `pdom/algorithm.py` | Partition-relaxed algorithm (PRD) with separate archive and selection partitions |
| `pdom/calibrate.py` | Calibration of every partition to a target number of occupied cells |
| `pdom/problems.py` | MaF1, MaF2, MaF3, MaF6 and MaF7, and the WFG and DTLZ problems of pymoo with analytic reference fronts |
| `pdom/metrics.py` | IGD+, normalized IGD+, hypervolume and the statistical tests used in the paper |
| `pdom/baseline_*.py` | Ports of GrEA, θ-DEA, SDR, LEO and MaOEA-HAP from their reference implementations |

The configurations of the paper correspond to

```python
run_prd(prob, archive="rank_grid", selection="anchor")      # rank-archive PRD
run_prd(prob, archive="rank_grid", selection="rank_grid")   # ordinal PRD
run_prd(prob, archive="uniform_grid", selection="anchor")   # affine-archive PRD
run_prd(prob, archive="qlmp", selection="anchor")           # Mondrian-archive PRD
```

and the ten baselines are available through `run_baseline(name, prob)` with
`name` in `moead`, `nsga3`, `rvea`, `ctaea`, `grea`, `theta_dea`, `sdr`,
`agemoea2`, `maoea_hap` and `leo`. MOEA/D, NSGA-III, RVEA, CTAEA and
AGE-MOEA2 come from pymoo.

## Reproducibility

Every algorithm draws its randomness from the seed it is given, so a run is
reproducible bit for bit with the same library versions. The MaOEA-HAP and
GrEA runs reported in the paper were made with an earlier version in which
part of their randomness was not seeded, so those runs can be repeated in
distribution but not bit for bit.

