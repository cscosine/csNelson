# csNelson

Nonlinear least squares in C++17 on Eigen, inspired by g2o: variable groups as sections (`SingleSection`, and
`DoubleSection` with the U / V / W partition for the Schur complement), unary / binary / n-ary edges (CRTP; the user
writes the residual, the Jacobians and the Hessian block updates), Gauss-Newton and Levenberg-Marquardt, on the
block-sparse matrices of [csBlockMatrix](https://github.com/cscosine/csBlockMatrix).

Built and released as part of [csoptimization](https://github.com/cscosine/csoptimization), with
[csCMake](https://github.com/cscosine/csCMake). Early stage: the review of the design and its next steps is
[ASSESSMENT.md](ASSESSMENT.md).
