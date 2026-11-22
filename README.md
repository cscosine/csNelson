# csNelson

Nonlinear least squares in C++ on Eigen, inspired by g2o

<!-- cs-template:begin (this block is owned by cs-project-template: do not edit; the rest of the file is the project's) -->
## Getting started

csNelson is built by the project that consumes it: no setup or build CI of its own. Formatting and checks:
`pre-commit install` once, then `pre-commit run --all-files` (also run by `.github/workflows/checks.yml`).

The shared files (setup, checks, CI, agent rules) come from
[cs-project-template](https://github.com/cscosine/cs-project-template); to bring them to its latest version, on a
branch: `copier update --defaults`.
<!-- cs-template:end -->

Nonlinear least squares in C++17 on Eigen, inspired by g2o: variable groups as sections (`SingleSection`, and
`DoubleSection` with the U / V / W partition for the Schur complement), unary / binary / n-ary edges (CRTP; the user
writes the residual, the Jacobians and the Hessian block updates), Gauss-Newton and Levenberg-Marquardt, on the
block-sparse matrices of [csBlockMatrix](https://github.com/cscosine/csBlockMatrix).

Built and released as part of [csoptimization](https://github.com/cscosine/csoptimization), with
[csCMake](https://github.com/cscosine/csCMake). Early stage: the review of the design and its next steps is
[ASSESSMENT.md](ASSESSMENT.md).
