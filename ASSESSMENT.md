# csNelson — Assessment (2026-11-16)

Scope: library code review (~6.7k LOC excluding tests), a Release build and run of `samples/BASchurPerfomances.cpp`,
and a comparison with Ceres, g2o and GTSAM. Companion document: `csBlockMatrix/ASSESSMENT.md`.

**Verdict:** the core design is sound and contains ideas worth building on, two of which suit GPUs
unusually well. The claim "faster than Ceres" is not verifiable today (no comparison harness, synthetic data),
and the measured CPU bottleneck sits exactly in the part that also has to be redesigned for GPU.

## What it is

g2o-inspired nonlinear least squares:

- **Sections**: `SingleSection` (one variable group) and `DoubleSection` (U/V/W partition, Schur-ready).
- **Edges**: unary / binary / n-ary, CRTP; the user writes the residual, the Jacobians, and the Hessian block
  updates (`updateH11Block`, `updateH12Block`, `updateH22Block`).
- **Algorithms**: Gauss-Newton, Levenberg-Marquardt.
- **Linear solvers**: dense / sparse Cholesky (Eigen `SimplicialLDLT`, natural / AMD / COLAMD ordering), Schur
  variants. The pipeline `SolverDiagonalBlocksInverseWWtMultSchur` chains V⁻¹ → W·V⁻¹ → bS →
  S = U − W·V⁻¹·Wᵀ → solve S → Wᵀ·x → back-substitute V.
- **Parallelism**: OpenMP with configurable schedule per stage; per-stage timing statistics.

## Strengths

1. **Lock-free, deterministic Hessian assembly (each output block gathers its own contributions).** At
   `structureReady()` every Hessian block becomes a computation unit with the list of contributions it
   receives; the b-segment is owned by the diagonal block's unit. No locks or atomics, and the summation order
   is fixed, so results are identical for any thread count. Ceres and g2o assemble under locks and don't
   guarantee this. It is also exactly the pattern a GPU wants (segmented reduction).
2. **Plan/execute Schur product.** `MatrixWWtMultiplier::prepare()` computes, once, which pairs of block UIDs
   contribute to each block of S; `multiply()` just walks flat lists. Together with choosing the storage
   layout for U, V, W and S at compile time, this is the real "customizable multiplier" innovation.
   Compile-time block sizes alone are not unique: Ceres (`SchurEliminator<2,3,6>`) and g2o (`BlockSolver_6_3`)
   have them too.
3. **Two-phase setup** (`parametersReady()`, `structureReady()`): the structure is frozen and there is no
   allocation per iteration.
4. **The Schur stages scale well with threads** (see below).

## Measured baseline

Setup: `BASchurPerfomances` built with `g++ -O3 -march=native -fopenmp`. Geometry is synthetic, the sparsity is
BAL ladybug-49 (48 variable cameras, ~7.7k points, 31.8k observations), cameras are SE3 + pinhole (6 dof, not
BAL's 9-dof cameras). Times are per Gauss-Newton iteration on a 16-core laptop.

| Stage | 1 thread | 8 threads |
|---|---|---|
| Edge evaluation + Hessian assembly | 12.5 ms | 7.1 ms (only 1.8×) |
| Schur stages (V⁻¹, W·V⁻¹, bS, S, Wᵀx) | ~3.5 ms | ~0.9 ms |
| Solve S (288×288, Eigen `SimplicialLDLT`) | 1.7 ms | 2.5 ms (no threading) |
| **Total** | **~18–21 ms** | **~11 ms** |

**Edge evaluation and assembly is the bottleneck.** Causes, from reading the code:

- every edge is a separately `new`-allocated polymorphic object;
- every Hessian contribution is a heap `EdgeHessianUpdater` called virtually, so each update follows the
  pointer chain updater → edge → section → block;
- U diagonal blocks in BA have ~650 contributions each (48 units), which balances badly with OpenMP's static
  schedule;
- in BA each W block has exactly one contributor, so the separate assembly pass over W is pure overhead;
- `clearAll()` runs on a single thread.

The S solve uses a simplicial sparse factorization on a nearly dense 288×288 matrix; dense LLT would be
several times faster, and for large problems (thousands of cameras) a supernodal factorization is mandatory.

### How to run the benchmark correctly

`BASchurPerfomances` **must** be given `ladybug_49_7776-sparsityPattern.txt` as its first argument. The built-in
default problem has zero edges and returns `SuccessAtStart` in 0.3 ms. Any old numbers produced without the
argument are meaningless.

## Bugs found

| Where | Issue |
|---|---|
| `SolverDiagonalBlocksInverseWWtMultSchur_impl.hpp:52` | `c = std::max(...)` assigns to the loop counter instead of `_uv_maxAbsHDiag`. This corrupts the loop and leaves V out of LM's initial damping value, in exactly the BA solver. |
| `DoubleSection_impl.hpp:295, :333` | `EdgeUnaryAdapterT::edgeSetterComputer[...]` is missing the `(section)` call, so an n-ary edge on a `DoubleSection` with two or more variable nodes won't compile. |
| `GaussNewton_impl.hpp` | No step rejection and no stopping test on relative cost change (TODO in code). Observed: cost 0.006 at iteration 0, then growing to ~10¹⁵ by iteration 50 → `MaxIterationReached`. |
| `samples/BASchurPerfomances.cpp:82-83` | `Eigen::Vector3d() * 0.001` is uninitialized memory, not random noise; pose noise is effectively absent. |
| `DoubleSection_impl.hpp` (n-ary / binary `addEdge`) | `assert(i.id() < j.id())` forces callers to order node ids; should be normalized internally. |

## Code health

- **About 63 copies of OpenMP loops.** Every kernel exists about 10 times, once per scheduling variant
  (`BaseSection.cpp`, all multipliers, `MatrixDiagInv`). One `parallel_for(range, fn, settings)` helper would
  remove a large share of the code, and it is the natural place to plug in a CUDA backend.
- **Template parameter explosion**: `DoubleSection` takes 11 template parameters, and solver selection goes
  through nested traits. This needs a facade of type aliases or preset configurations
  (e.g. `BAProblem<CamDim, PointDim>`).
- **About 80 test files differ only in storage layout** (`testDoublesection{F,D,V}{F,D}{F,D}{F,D}.cpp`); use
  `TEMPLATE_TEST_CASE`.
- No robust losses, no automatic differentiation, no covariance, only two variable groups (U and V).

## Compared with existing libraries

| | Ceres | g2o | GTSAM | csNelson today |
|---|---|---|---|---|
| Main focus | General nonlinear least squares | Graph SLAM / BA | Factor graphs, incremental (iSAM2) | Structured batch problems |
| Jacobians | Automatic differentiation | Hand-written or numeric | Hand-written or expression-based | Hand-written only |
| Robust loss functions | Yes | Yes | Yes | **No** |
| Schur elimination | Highly optimized, sizes fixed at compile time | Block solver | Multifrontal elimination | Plan/execute, storage layout chosen per matrix |
| Sparse direct solver | SuiteSparse supernodal, Accelerate, Eigen | CHOLMOD, CSparse, Eigen | Own multifrontal | **Eigen simplicial only** |
| Iterative solver (PCG) | Yes (`ITERATIVE_SCHUR`) | Yes | Yes | No |
| GPU | Partial (CUDA dense / CG) | No | No | No |
| Covariance / incremental | Yes / No | Partial / No | Yes / Yes | No / No |
| More than 2 variable groups | n/a | n/a | n/a | No |

**Where csNelson can honestly win today:** small and medium BA with hand-written Jacobians, where Ceres pays
for its generality (automatic differentiation, residual- and parameter-block bookkeeping, locked assembly).
But a test of hand-written Jacobians against Ceres's automatic ones compares Jacobian strategies, not solvers;
a fair benchmark must give every library the same Jacobians.

**Where it loses today:** large BA (simplicial vs supernodal factorization, no PCG), and anything needing
robust losses.

**Positioning:** don't out-generalize Ceres, and don't compete with GTSAM on incremental solving. A credible niche
is *the fastest batch solver for problems whose structure is known up front (BA, pose–landmark SLAM,
calibration), with results that don't depend on thread count, and the same edge code running on CPU and GPU.*

## GPU feasibility

The architecture is more GPU-ready than most libraries:

- the gather-style assembly lists and the precomputed Schur pair lists translate directly into GPU kernels
  (one warp or thread block per output block);
- the Schur stages are batched small dense operations (6×6, 6×3, 3×3); V⁻¹ is a trivially parallel batch of
  3×3 inversions.

Obstacles:

1. **Edges are virtual C++ objects running on the host.** They need to become plain data structs plus a static
   `evaluate(params, J, r)` that compiles for both host and device (Eigen's fixed-size types work in CUDA).
   Store them grouped by edge type, field by field ("structure of arrays"). **This is also the fix for the
   CPU bottleneck above.**
2. **Block storage must be flat** (see the csBlockMatrix assessment, BSR-like buffers).
3. **Precision.** Consumer and laptop Ampere GPUs (the dev machine has an RTX A3000 Laptop, 6 GB) run FP64 at
   about 1/64 of FP32 speed, so a useful GPU path must work in float or mixed precision. References showing
   single precision can be stable for BA: Root BA (Demmel et al., CVPR 2021) and Power Bundle Adjustment
   (Weber et al., CVPR 2023).
4. **Solving S on the GPU.** Options:
   - cuSOLVER dense `potrf` for S up to roughly 2k cameras;
   - cuDSS or BaSpaCho (Meta's supernodal Cholesky, designed for BA and SLAM) for sparse S;
   - PCG with a Schur-Jacobi or power-series preconditioner for very large problems.

   Other GPU BA work to compare against: MegBA, PBA (Wu 2011), DeepLM, Theseus.

## Roadmap

### Phase 0: get the facts (weeks)
- Fix the bugs above. Gauss-Newton gets step acceptance and a relative cost-change stopping test.
- Benchmark harness on **real** BAL datasets (several sizes). Ceres, g2o and GTSAM all use the **same
  hand-written Jacobians** and the same stopping rules. Report time per iteration **and** time to reach a
  given cost; fix the thread count; collect stats with the existing timing hooks.
- Outcome: a factual answer to "does it beat Ceres, and where?".

### Phase 1: rework the CPU core
- Edges stored by type as structure of arrays; Jacobians in contiguous buffers.
- Fuse edge evaluation and assembly where a block has a single contributor (V and W in BA); keep gather-style
  assembly only where needed (U).
- Contribution lists as flat CSR-style index arrays instead of `unique_ptr<EdgeHessianUpdater>` objects.
- `parallel_for` backend helper replacing the duplicated OpenMP loops; dynamic schedule or weight-balanced
  partitions.
- Choose the S solver by size and density: dense LLT, CHOLMOD supernodal, Eigen simplicial.
- Robust losses (applied when weighting the Jacobian and residual, not as an afterthought).
- Target: under 5 ms per iteration on ladybug-49 at 8 threads.

### Phase 2: features
- N variable groups (generalizing U and V), e.g. BA with shared intrinsics, VIO.
- PCG on the reduced system with a Schur-Jacobi preconditioner.
- Optional lightweight automatic differentiation (small dual-number type) for prototyping.
- Covariance and marginals of selected blocks.
- BAL and g2o file readers; facade types to hide the template parameters.

### Phase 3: GPU
- Step 1: move only the Schur pipeline to the GPU, keeping edge evaluation on the CPU; measure how much the
  transfers cost.
- Step 2: device-callable edges, so evaluation and assembly happen on the GPU.
- Step 3: float or mixed precision; cuDSS or BaSpaCho for S, or PCG for very large problems.

## Related work for the linear-algebra layer

Full list in `csBlockMatrix/ASSESSMENT.md` → "Related work". Most relevant to csNelson:

| Library | Relevance |
|---|---|
| **SLAM++** (Brno) | Closest relative: block matrices for SLAM/BA, fixed-size kernels picked at runtime from a compile-time list of block sizes; also an incremental solver. Benchmark candidate. |
| **BaSpaCho** (Meta) | Supernodal Cholesky on variable-size parameter blocks, CPU + CUDA. Prime candidate for the S solve (Phase 1 on CPU, Phase 3 on GPU). |
| **CHOLMOD** (SuiteSparse) | Supernodal CPU factorization used by Ceres and g2o; needed to compete on large BA. |
| **cuDSS** (NVIDIA) | GPU direct sparse solver (CSR input); alternative to BaSpaCho for Phase 3. |
| **cuSPARSE BSR / Kokkos Kernels `BsrMatrix`** | Matrix–vector products for PCG on GPU when block sizes are uniform. |
| **MAGMA / cuBLAS batched, libxsmm** | Batched small dense kernels for V⁻¹, W·V⁻¹ and the S pair products. |
| **SymForce** (Skydio) | Code generation of residuals and Jacobians from symbolic expressions; a model for device-callable edges without hand-written Jacobians. |
| **Theseus** (Meta), **MegBA**, **PBA**, **DeepLM** | GPU nonlinear least squares and BA references for the Phase 3 design and for benchmarking. |

Verify each project's current state before depending on it (list compiled from memory, 2026-11-16).

---

# csNelson as a general SLAM and mapping solver

Goal: csNelson becomes the solver for all applications. The focus is SLAM and mapping, including
continuous-time (Gaussian-process priors, splines) and discrete components, but it is not limited to them.

## Is GTSAM "just Schur + marginalization"?

Algebraically, yes: Schur elimination, marginalization and iSAM2 are all **variable elimination** (block Gaussian
elimination). What differs is four choices:

1. **Which variables get eliminated, and in what order.** Schur elimination is "eliminate the landmarks first";
   GTSAM eliminates on a general ordering (COLAMD, METIS, or constrained versions).
2. **What you keep afterwards.**
   - Marginalization keeps a dense prior on the remaining variables; the eliminated ones are gone for good.
   - iSAM2 keeps the factored form (a Bayes tree), so nothing is thrown away.
3. **When you relinearize.**
   - Marginalization freezes the linearization point of the prior, so it needs first-estimate Jacobians (FEJ)
     to stay consistent.
   - iSAM2 relinearizes only the variables whose estimate moved more than a threshold ("fluid
     relinearization"), and re-eliminates only the affected part of the tree.
4. **How much work an update costs.**
   - Sliding-window cost is bounded.
   - iSAM2 is cheap for odometry updates, but a loop closure can force re-elimination of a large part of the
     tree.

Three regimes are needed for a general solver:

| Regime | Used by | Mechanism |
|---|---|---|
| **Batch** | Mapping, offline SfM/BA, calibration | Full solve with Schur elimination or a supernodal factorization |
| **Fixed-lag / sliding window** | VIO and LIO front ends (OKVIS, VINS-Mono, Basalt) | Schur marginalization into a prior factor, with FEJ |
| **Incremental smoothing** | Full SLAM back ends | iSAM2 Bayes tree, or an incremental Cholesky update as in SLAM++ |

## What GTSAM offers beyond the algebra

- **Factor graph modelling**: heterogeneous variables, noise models, robust kernels, and the `Expression<T>`
  system for building factors with automatic Jacobians.
- **Lie groups out of the box**: `Pose3`, `Rot3`, `NavState`, and so on.
- **SLAM building blocks**: IMU preintegration (`CombinedImuFactor`), smart projection factors that eliminate
  landmarks inside the factor, and fixed-lag smoothers.
- **Marginal covariances** recovered efficiently from the Bayes tree, needed for data association and for
  gating loop closures.
- **Discrete and hybrid inference**: `DiscreteFactorGraph`, and hybrid discrete–continuous graphs for data
  association, mode switching and outlier hypotheses.
- **Ecosystem**: Python bindings, and a lot of papers and code built on it (LIO-SAM, Kimera, GPMP2).

## Is GTSAM the de facto standard? Has anyone beaten it?

It is the standard for **incremental SLAM back ends in research**, but not universal:

- Ceres dominates batch problems and is used by VINS-Mono, OKVIS and Cartographer.
- g2o is used by ORB-SLAM.
- Several of the best systems use their own solver: Basalt, DSO, ICE-BA.

On performance, GTSAM has been beaten in specific setups (figures from memory; verify before quoting):

- **SLAM++** (Ila, Polok et al., IJRR 2017) reported faster incremental solving than iSAM2 on standard SLAM
  datasets.
- **ICE-BA** (Liu et al., CVPR 2018) beat iSAM2 on visual-inertial bundle adjustment.
- **SymForce** (Skydio, RSS 2022) reported large speedups over both GTSAM and Ceres on some robotics problems, by
  generating flattened residual and Jacobian code.
- **Basalt**'s custom VIO solver is much faster than general-purpose ones.

GTSAM's weaknesses are the same as csNelson's current eval path, only more so:

- `shared_ptr` factors with virtual calls everywhere;
- dynamic-size Eigen matrices inside factors;
- heavy allocation during each update;
- limited multithreading and no GPU.

Its breadth, theory and ecosystem are much harder to beat. Realistic strategy: win on speed, deterministic
results and GPU support in the core SLAM regimes, and read the common graph file formats so you can compare
directly.

## What csNelson is missing as a general solver

### Tier 1: the core that must change

1. **Dynamic structure.** This is the biggest gap. `structureReady()` freezes the problem once, but SLAM adds
   variables and factors every frame and sometimes removes them. Needed: incremental structure updates with plan
   reuse (add variables and factors, recompute only the affected parts of the computation plan).
2. **Any number of heterogeneous variable groups**, instead of only Single or Double sections: poses,
   landmarks, velocities, biases, calibration, spline control points, GP states.
3. **A general elimination engine instead of hard-coded Schur.**
   - Express the elimination order as a plan: "eliminate group L, then factor the rest" is Schur; "eliminate
     the oldest state into a prior" is marginalization.
   - This keeps the plan/execute idea and makes it the unifying core.
   - Schur elimination, marginalization and incremental updates then become three uses of one engine.
4. **Manifold variables as first-class objects.** Today the user writes one `oplus` for the whole problem. Each
   variable type should carry its own retraction and tangent dimension; csLie already exists.
5. **Robust estimation**: M-estimators (Huber, Cauchy, and so on), plus graduated non-convexity (GNC) for
   outlier-heavy loop closures.
6. **Trust-region methods**: LM with Marquardt diagonal scaling, Dogleg, and proper stopping criteria.
7. **Linear solvers**: supernodal Cholesky, QR or square-root forms (numerically stable, and they enable float
   on GPU), and PCG.
8. **Covariance and marginal recovery**: a sparse inverse of selected blocks via the Takahashi recursion.

### Tier 2: SLAM and mapping specifics

- **Marginalization prior factor**: a generic dense linear factor, FEJ support, and fixed-lag smoothing on top.
- **Incremental smoothing**: either an iSAM2-style Bayes tree or an incremental Cholesky update in the SLAM++
  style. Do this after the fixed-lag smoother: it is harder and fewer real systems use it in the front end.
- **Continuous time**:
  - **GP motion priors** (Barfoot's STEAM, GPMP2) connect consecutive states, so the matrix is
    block-tridiagonal. That is exactly the `Periodic<6,3,6>`-style pattern proposed in
    `csBlockMatrix/ASSESSMENT.md`: constant structure, ideal for compile-time patterns.
  - **B-splines on Lie groups** (Sommer et al., CVPR 2020; Basalt): each residual touches k consecutive control
    points (4 for cubic splines), so the matrix is banded. These are n-ary edges with fixed arity; csNelson's
    existing n-ary edges are a good base.
- **IMU preintegration** and **smart factors**, where landmarks are Schur-eliminated inside the factor.

### Tier 3: discrete and hybrid inference

- Start with **robust and switchable constraints plus GNC**. They cover most real outlier and data-association
  needs within continuous optimization.
- Full **hybrid discrete–continuous inference** (GTSAM hybrid graphs, multi-hypothesis iSAM) is a
  research-grade feature. Design the factor interface so it is possible later, but don't build it first.

### Tier 4: usability

Automatic differentiation or code generation (SymForce-style), Python bindings, readers and writers for g2o,
BAL and TORO files, and a facade over the long template lists.

## Direction

Keep what differentiates csNelson, and generalize around it:

1. **Core**: a block elimination engine with plan/execute, N variable groups and dynamic structure. Lock-free
   assembly, deterministic results, and edges stored by type as plain data that can also run on the GPU.
2. **Regimes, in order**: batch, then fixed-lag with marginalization, then incremental.
3. **Factors**: robust kernels, IMU preintegration, GP priors, B-splines, marginalization priors.
4. **Backends**: CPU first with OpenMP behind a `parallel_for` helper; GPU later for batch and large windows.

Pitch: *one solver for batch, sliding-window and incremental SLAM, with continuous-time priors as first-class
structure, faster than GTSAM and Ceres in those regimes, with deterministic results and the same factor code on
CPU and GPU.* Beating GTSAM on breadth is not the target.

## Revised roadmap

This supersedes the phases in the "Roadmap" section above; the earlier Phase 0 and Phase 1 items are kept and
folded in.

### Phase 0: get the facts
- Fix the bugs listed above; Gauss-Newton step acceptance and stopping criteria.
- Benchmark harness against Ceres, g2o and GTSAM with the same Jacobians and stopping rules, on:
  - BAL (batch BA);
  - pose-graph datasets in g2o/TORO format (Intel, Manhattan, sphere, garage) for batch and incremental;
  - one VIO sequence (e.g. EuRoC) for fixed-lag.

### Phase 1: the new core (batch)
- Edges stored by type as structure of arrays; flat CSR-style contribution lists; fused evaluation and
  assembly; `parallel_for` backend helper.
- N heterogeneous variable groups with per-type manifold retraction (csLie).
- A general elimination plan, with Schur elimination as one instance; supernodal / dense / PCG S solvers.
- Robust kernels and GNC; LM with Marquardt scaling and Dogleg.
- Covariance recovery for selected blocks.

### Phase 2: dynamic structure and fixed-lag
- Add and remove variables and factors with incremental plan updates.
- Marginalization into a dense prior factor; FEJ; fixed-lag smoother.
- IMU preintegration, GP motion priors, B-spline factors (using csBlockMatrix `Periodic<...>` patterns).

### Phase 3: incremental smoothing
- Incremental factorization: a Bayes-tree-style or incremental Cholesky update, with selective relinearization.
- Benchmark against iSAM2 and SLAM++ on the pose-graph datasets.

### Phase 4: GPU and breadth
- GPU backend (as in the earlier Phase 3): Schur pipeline first, then device-callable edges, then float or
  mixed precision.
- Code generation or automatic differentiation for factors, Python bindings.
- Hybrid discrete–continuous inference, if the applications need it.
