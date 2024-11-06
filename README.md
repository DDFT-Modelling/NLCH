# NLCH
Extension for the Nonlocal Cahn–Hilliard Equation with singular potentials

---

### How do I run tests? ###
 
1. Please run [`AddPaths.m`](AddPaths.m) in the main folder. 
2. Run any Matlab file in this folder or any of the following:

	* [`Kernel builder`](Kernel%20builder): Code for building spectral element approximation of the Newtonian Kernel.

3. (Optional) If the code requires a kernel that has been split into several files (e.g., `Singular_Kernels_Subs_40_epsB_[split_1].mat`), then rebuild it using the code available in [File_Splitter.m](Kernel%20builder/File_Splitter.m).





### Requirements (Matlab Versions)

The code was tested using Matlab R2023b and R2024a. The fastest kernel creation times (to date) can be achieved in R2023b.


### Additional Python code

Two additional notebooks are presented:

* [`Exact singular convolutions.ipynb`](Exact%20singular%20convolutions.ipynb) contains a lengthly computation of the Newtonian potential against a constant function for a small neighbourhood: $[K \ast \iota_{ B_\infty (\cdot; \varepsilon) \cap [0,1]^2 }] (x)$. Many constants related to this evaluation are computed analytically here and the result is plotted against several values of $\varepsilon$.

* [`Exact singular convolutions - Numerical test.ipynb`](Exact%20singular%20convolutions%20-%20Numerical%20test.ipynb) computationally verifies the analytical formulas obtained in the previous notebook using simple quadrature methods that circumvent the singularity.

---

## About `MultiShape`

`MultiShape` is a library developed by Ben Goddard, Rory Mills-Williams, John Pearson, and Jonna Roden. It builds spectral element methods based on the `2DChebClass` library. The baseline package can be found in the public repository [`MultiShape`](https://bitbucket.org/bdgoddard/multishapepublic/src/master/). A detailed presentation can be found in **[1]**.


---

## About `2DChebClass`


`2DChebClass` is a library developed by Andreas Nold and Ben Goddard. It can be used to solve a wide range of (integro-)differential systems in various 1D and 2D geometries. The baseline package can be found in the public repository [`2DChebClass`](https://github.com/NoldAndreas/2DChebClass). A detailed presentation can be found in **[2]**.

---

**[1]** Jonna C. Roden, Rory D. Mills-Williams, John W. Pearson, and Benjamin D. Goddard, 2024 "MultiShape: a spectral element method, with applications to Dynamical Density Functional Theory and PDE-constrained optimization." _IMA Journal of Numerical Analysis,_ drae066. Links: [ArXiv](https://arxiv.org/abs/2207.05589), [IMA JNA](https://doi.org/10.1093/imanum/drae066)


**[2]** Andreas Nold, Benjamin D. Goddard, Peter Yatsyshin, Nikos Savva & Serafim Kalliadasis, 2017 "Pseudospectral methods for density functional theory in bounded and unbounded domains." _Journal of Computational Physics, Elsevier, 334, 639–664. Links: [ArXiv](https://arxiv.org/abs/1701.06182), [J Comp Phys](https://doi.org/10.1016/j.jcp.2016.12.023)
