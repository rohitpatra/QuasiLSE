# QuasiconvexLSE
This is an R package to compute the multivariate quasiconvex/quasiconcave nonparametric LSE with or without additional monotonicity constraints as described in "Least Squares Estimation of a Monotone Quasiconvex Regression Function" by Somabha Mukherjee, Rohit K. Patra, Andrew L. Johnson, and Hiroshi Morita, which can be found at the following link:

https://arxiv.org/abs/2003.04433

To download the R package use the following in R:
```
library(devtools)
install_github(repo = "rohitpatra/QuasiLSE")
library(QuasiLSE)
```

Note that the above package requires an installation CPLEX or Gurobi and the R-package RCPLEX/gurobi 
Files in the folder titled "ReplicationCode" replicate Figure 5, 6, and 7 of https://arxiv.org/abs/2003.04433. 

**References**

Somabha Mukherjee, Rohit K. Patra, Andrew L. Johnson, and Hiroshi Morita. **Least Squares Estimation of a Monotone Quasiconvex Regression Function**. 2020. arXiv:2003.04433

