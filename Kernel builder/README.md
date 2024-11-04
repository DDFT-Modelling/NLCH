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
	1. Compute and plot the maximum error of the approximation of $K \ast 1$, see examples [`Errors_Max_Factor[10].pdf`](Errors_Max_Factor[10].pdf) and [`Errors_Max_Factor[20].pdf`](Errors_Max_Factor[20].pdf).
	2. 
	
It also stores the kernels for 


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
