# Parallel Stabilized Column Generation (PSCG) Algorithm


Author: Brian C. Dandurand (2017-2019)  
Postdoctoral Appointee  
Division of Mathematics and Computer Science  
Argonne National Laboratory  
Lemont, Illinois, USA 

PSCG is a C++ implementation of the algorithm developed in

Boland, N., Christiansen, J., Dandurand, B., Eberhard, A., Oliveira, F., 
"A parallelizable augmented Lagrangian method applied to large-scale nonconvex-constrained optimization problems," *Mathematical Programming A*, 1–34, 2018.

It is based on C++ code originally developed by Jeffrey Christiansen, Brian Dandurand, and Fabricio Oliveira
at RMIT in Melbourne Australia. 

The developments in the above citation were funded under the Australian Research Council Project ARC DP 140100985 during 2015-2016.
CIs of that projects were Prof. Andrew Eberhard (RMIT) and Prof. Natashia Boland (Georgia Tech); and PI Prof. Jeffrey Linderoth (U. Wisconsin-Madison).

Functionally, PSCG is similar to the alternating direction method of multipliers (ADMM) with the following two exceptions:

1) Like ADMM, primal updates are performed via a two-block Gauss-Seidel iteration; 
unlike ADMM, one of those block updates correspond to an iteration of the Frank-Wolf method applied to the block specific subproblem rather than computing an exact minimization;
2) Dual updates are taken conditionally based on a serious step condition common with proximal bundle methods, rather than unconditionally in ADMM.

The proof of optimal convergence under mild assumptions is developed in the above citation.
PSCG has been applied to finding the optimal solution to the Lagrangian dual of a two-stage stochastic
integer program. 

Currently, the problem needs to be in SMPS format. 
Furthermore, the current implementation uses a third party library DSP 
developed by Kibaek Kim at Argonne National Laboratory to read the necessary subproblems from SMPS format.
