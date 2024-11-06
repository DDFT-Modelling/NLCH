# NLCH
Extension for the Nonlocal Cahn–Hilliard Equation with singular potentials

---

### How do I run tests? ###
 
1. Please run [`AddPaths.m`](AddPaths.m) in the main folder. 
2. Run any Matlab file in this folder or any of the following:

	* [`Kernel builder`](Kernel%20builder): Code for building spectral element approximation of the Newtonian Kernel.
	
	* [`Compressed Kernels`](Compressed%20Kernels): Kernels for $N\in \\{40,50\\}$ that are split into several files due to GitHub's file-size limit.

	* [`Source`](Source): Solves the NLCH with a source. Two examples are included: a constant source and a decaying wave.

	* [`DAE Tests`](DAE%20Tests): Examples of the DAE solver for a 1D and 2D heat equations.


3. (Optional) If the code requires a kernel that has been split into several files (e.g., `Singular_Kernels_Subs_40_epsB_[split_1].mat`), then rebuild it using the code available in [File_Splitter.m](Kernel%20builder/File_Splitter.m).

## Files in this folder

A gallery of examples is presented and based on the following scripts:

* [`NL_CH_Integrator_Simple.m`](NL_CH_Integrator_Simple.m): Solves the NLCH system using Matlab's own DAE solver for Neumann data.
* [`NL_CH_Integrator_DAE.m`](NL_CH_Integrator_DAE.m): Solves the NLCH system with a DAE solver for Neumann data.

Since the kernel has already been created in the other folders, we explore the scalings

$$
\begin{align*}
	K_\eta : \mathbb{R}^2 \setminus\{0\} \ni x  \longmapsto  \frac{\eta}{2\pi} \log \|x\| = \frac{\eta}{4\pi} \log ( x_1^2 + x_2^2 ) \in \mathbb{R},
\end{align*}
$$

for $\eta > 0$.

A minimal working example is presented in [`Run_Example.m`](Run_Example.m) which solves the system for a wave of compact support and the convolution kernel is scaled by $\eta = 300$, see [`MWE.gif`](MWE.gif) for an animation.

Additional examples:

A. [`Periodic_Wave.m`](Periodic_Wave.m): Solves the NLCH system with initial condition $\rho_0(x) = \sin(2\pi x_1) \cos(2\pi x_2)$. Energy plots are presented in [`NLCH_Energy_PW.pdf`](NLCH_Energy_PW.pdf), an example is animated in [`A_eta=-1.gif`](A_eta=-1.gif), example density plots are presented in [`A_Solution_PW.pdf`](A_Solution_PW.pdf). More figures for each case are presented in the folder [`Gallery_PW`](Gallery_PW).

B. [`Constant.m`](Constant.m): The initial condition is a constant state $\rho_0 \equiv -\frac{1}{2}$. Figures associated with this experiment are presented in the folder [`Gallery_C`](Gallery_C).

C. [`Compact.m`](Compact.m): The initial condition is a wave with compact support. Figures associated with this experiment are presented in the folder [`Gallery_CW`](Gallery_CW).

D. [`Regular.m`](Regular.m): The initial condition is the same as above, but we use a regular potential to approximate the action of the singular potential. Outputs are presented in [`Gallery_Regular`](`Gallery_Regular`).




### Requirements (Matlab versions)

The code was tested using Matlab R2023b and R2024a. The fastest kernel creation times (to date) can be achieved in R2023b.


### Additional files

* [`plots_in_box.m`](plots_in_box.m) returns an animation of a spatio–temporal density.
* [`plot_panels.m`](plot_panels.m) generates density plots.


### Additional Python code

Two additional notebooks are presented:

* [`Exact singular convolutions.ipynb`](Exact%20singular%20convolutions.ipynb) contains a lengthly computation of the Newtonian potential against a constant function for a small neighbourhood: $[K \ast \iota_{ B_\infty (\cdot; \varepsilon) \cap [0,1]^2 }] (x)$. Many constants related to this evaluation are computed analytically here and the result is plotted against several values of $\varepsilon$.

* [`Exact singular convolutions - Numerical test.ipynb`](Exact%20singular%20convolutions%20-%20Numerical%20test.ipynb) computationally verifies the analytical formulas obtained in the previous notebook using simple quadrature methods that circumvent the singularity.

---

### About `MultiShape`

`MultiShape` is a library developed by Ben Goddard, Rory Mills-Williams, John Pearson, and Jonna Roden. It builds spectral element methods based on the `2DChebClass` library. The baseline package can be found in the public repository [`MultiShape`](https://bitbucket.org/bdgoddard/multishapepublic/src/master/). A detailed presentation can be found in **[1]**.


---

### About `2DChebClass`


`2DChebClass` is a library developed by Andreas Nold and Ben Goddard. It can be used to solve a wide range of (integro-)differential systems in various 1D and 2D geometries. The baseline package can be found in the public repository [`2DChebClass`](https://github.com/NoldAndreas/2DChebClass). A detailed presentation can be found in **[2]**.

---

**[1]** Jonna C. Roden, Rory D. Mills-Williams, John W. Pearson, and Benjamin D. Goddard, 2024 "MultiShape: a spectral element method, with applications to Dynamical Density Functional Theory and PDE-constrained optimization." _IMA Journal of Numerical Analysis,_ drae066. Links: [ArXiv](https://arxiv.org/abs/2207.05589), [IMA JNA](https://doi.org/10.1093/imanum/drae066)


**[2]** Andreas Nold, Benjamin D. Goddard, Peter Yatsyshin, Nikos Savva & Serafim Kalliadasis, 2017 "Pseudospectral methods for density functional theory in bounded and unbounded domains." _Journal of Computational Physics, Elsevier, 334, 639–664. Links: [ArXiv](https://arxiv.org/abs/1701.06182), [J Comp Phys](https://doi.org/10.1016/j.jcp.2016.12.023)
