# Multicellular Phase Field Model
## Model
$\phi_i(r,t)$ represents i-th cell, $\phi_i=1$ for interior of the cell and $\phi_i=0$ for exterior, $i=1,2,\cdots,N$. $\phi_e(r)$ repersents the eggshell.

**Surface tension** 

$$\mathbf{F}_{\text{ten},i}=-\gamma\left[\Delta\phi_i-cW'(\phi_i)\right]\frac{\mathbf\nabla\phi_i}{|\mathbf\nabla\phi_i|^2},$$
where $W(\phi_i)=\phi_i^2(1-\phi_i)^2$.

**Repulsion**

$$\mathbf {F}_{\text{rep},i}=\left(g_{e}\phi_i\phi_e^2+g\phi_i\sum_{j\neq{i}}^N\phi_j^2\right)\frac{\mathbf\nabla\phi_i}{|\mathbf\nabla\phi_i|^2}.$$

**Volume constraint**
$$\mathbf{F}_{\text{vol},i}=M\left[\frac{\int_\Omega\phi_idr}{V_i(t)}-1\right]{\mathbf{n}},$$
where $V_i(t)$ is the prescribed volume for the i-th cell. ${\mathbf{n}}=\frac{\nabla\phi_i}{|\nabla\phi_i|}$ is the inward normal vector.

**Cell-cell adhesion**
$$\mathbf{F}_{\text{adh},i}=-\sum^N_{j\neq{i}}\sigma_{ij}\mathbf\nabla\phi_j.$$

**Evolution equation of $\phi_i$**

$$\tau\frac{\partial\phi_i}{\partial{t}}=-(\mathbf{F}_{\text{ten},i}+\mathbf{F}_{\text{rep},i}+\mathbf{F}_{\text{vol},i}+\mathbf{F}_{\text{adh},i})\cdot\mathbf\nabla\phi_i.$$
## Instruction
Sim8to12.m implements the simulations

`mainBDF2VC()` is main funcion for simulation, calling `ShapeInit()`, `GridInit()`, `CellDivision()`, `SimParaModify()`.

CellParameters.mat file storing biological information of cells, including cell names, volumes, lineage, division axes. The biological information is obtained from *[Cao 2020](https://doi.org/10.1038/s41467-020-19863-x)*. `CellParameters` is a struct