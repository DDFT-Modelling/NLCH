# Approximating the singular kernel

We compute the action of the Newtonian kernel as

$$
\begin{align}
	[K \ast \rho] (x) &= \int\limits_{ \Omega } K(x-y) \rho(y) \ \mathrm{d}y
	\approx \int\limits_{ \Omega \setminus B_\infty (x;\varepsilon) } K(x-y) \rho(y) \ \mathrm{d}y + \rho(x) \int\limits_{ B_\infty (x;\varepsilon) } K(x-y) \ \mathrm{d}y  = \mathbb{M} \rho + \rho \circ \mathtt{G}_\varepsilon.
\end{align}
$$

---
## Main files

* [`Newtonian.m`](Newtonian.m) yields the necessary handles to compute $\mathtt{G}_\varepsilon$.

* [`SingularKernelIntergralMS_Conv_Maximal_Q.m`](SingularKernelIntergralMS_Conv_Maximal_Q.m) computes $\mathbb{M}$ and $\mathtt{G}_\varepsilon$ using a maximal partition of the domain $[0,1]^2$ at each pseudospectral collocation point based on a Chebyshev grid. The function takes (at most) four arguments: the neighbourhood radius $\varepsilon$, number of collocation points in the first direction $N_x$, number of collocation points in the second direction $N_y$, and resolution factor $\alpha$.

* [`Errors_SingularConv_MORE_EPS_grid_Maximal.m`](Errors_SingularConv_MORE_EPS_grid_Maximal.m) performs several error tests against many choices of $(\varepsilon,N_x,N_y,\alpha)$. It includes subroutines to:
	1. Compute and plot the maximum error of the approximation of $K \ast 1$, see examples [`Errors_Max_Factor[10].pdf`](Errors_Resolution/Errors_Max_Factor[10].pdf) and [`Errors_Max_Factor[20].pdf`](Errors_Resolution/Errors_Max_Factor[20].pdf).
	2. Plot swarm plots comparing the approximation errors against different resolution factors $\alpha$, see examples [`Test_Singular_Convolution_Fixed_Swarm[10,2,10,16,3].pdf`](Errors_Resolution/Test_Singular_Convolution_Fixed_Swarm[10,2,10,16,3].pdf) and [`Test_Singular_Convolution_Fixed_Swarm[20,1,5,8,3].pdf`](Errors_Resolution/Test_Singular_Convolution_Fixed_Swarm[20,1,5,8,3].pdf).
	3. Store a kernel for a given pair of dimensions $(N_x, N_y)$, a fixed radius $\varepsilon > 0$, and several factors $\alpha$ in a structured object. For example, [`Singular_Kernels_Subs_10_epsA.mat`](Singular_Kernels_Subs_10_epsA.mat) contains the fields
		a. `I`, `J`, and `G` which are used to evaluate $\mathtt{G}_\varepsilon$.
		b. `NI`, `NJ`, and `NG` which are the numerical evaluations of the functions above at the collocation ponts.
		c. `eps` is the value of $\varepsilon$.
		d. `box` is a `2DChebClass` box object with Chebyshev collocation points.
		e. `N`$=(N_x,N_y)$ is the number of collocation points for each direction.
		f. `Level` is a structure that provides, for each factor $\alpha$, the approximation of $\mathbb{M}$.
	4. Plot errors over domain as intensity values at each collocation point, see for example [`Test_Singular_Convolution_Domain[10,1,0.01].pdf`](Errors_Domain/Test_Singular_Convolution_Domain[10,1,0.01].pdf) and [`Test_Singular_Convolution_Domain[10,10,0.01].pdf`](Errors_Domain/Test_Singular_Convolution_Domain[10,10,0.01].pdf).

	

* [`Maximal_Grid_40x40_epsB.m`](Maximal_Grid_40x40_epsB.m) a very condensed version of the previous code that only computes and stores the kernel for $\varepsilon = 10^{-5}$, $N=40$, and three factors.

---

## Supplementary files

* [`Newtonian_Grad_IJ.m`](Newtonian_Grad_IJ.m) computes the gradient of $\mathtt{G}_\varepsilon$ (not used in construction).


---

## Outputs

* [`Plot_Errors_Singular_Conv_eps.m`](Plot_Errors_Singular_Conv_eps.m) plots the error of approximating $K \ast 1$ in three different formats:
	1. Swarm plot of the approximation error as $\varepsilon$ changes
	2. Absolute error (per node) in $[0,1]^2$
	3. Absolute error (per node) in $\Omega \setminus B_\infty (x;\varepsilon)$.

	A selection of these plots is included in the folder [`Errors_All_Mesh`](Errors_All_Mesh).

* [`Plot_Errors_Singular_Conv_Facts.m`](Plot_Errors_Singular_Conv_Facts.m) plots a swarm plot comparing the errors for different resolution factors $\alpha$.

* [`Plot_Errors_Singular_Conv_Domain.m`](Plot_Errors_Singular_Conv_Domain.m) plot errors as intensity values over the domain at each collocation point.
