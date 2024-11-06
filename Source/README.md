# NLCH with source

We solve the NLCH system with a source. 

---
## Main files

* [`Source_NL_CH_Constant.mlx`](Source_NL_CH_Constant.mlx): Solves the equation where the source is computed using the exact solution $\rho \equiv \frac{1}{2}$.

* [`Source_NL_CH.mlx`](Source_NL_CH.mlx): Solves the equation where the source is computed using the exact solution $\rho(x;t) = e^{-t} \sin(2\pi x_1) \cos(2\pi x_2)$.

* [`Source_NLCH_Tester_Plots.m`](Source_NLCH_Tester_Plots.m): Replicates the plots for the `N=20` tests.

	
* [`NLCH_Tester_20.m`](NLCH_Tester_20.m): Runs a test for the `N=20` kernels without DAE. Was used to create [`Animations 20`](Animations%2020).