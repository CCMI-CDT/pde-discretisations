We can write the 1D Advection equation as:
$$\frac{\partial u}{\partial t} = -c \frac{\partial u}{\partial x}$$

One can show that $\frac{1}{2} \int u^2 dx$ is conserved with periodic boundary conditions.
Now consider the semi-discretised equation in space; i.e. a (finite dimensional) ODE in continuous time:

$$\dot{u}=-\frac{c}{2\Delta x} (u_{n+1}-u_{n-1})=Au$$

where $A$ is the circulant matrix:

```math
A=\begin{pmatrix}
0 & -\frac{c}{2} & 0 & \dots & 0& \frac{c}{2} \\
\frac{c}{2} & \ddots & \ddots & \ddots && 0 \\
0 & \ddots & \ddots & \ddots& \ddots& \vdots \\
\vdots  &\ddots &\ddots&\ddots&\ddots&0\\\
0&&\ddots&\ddots&\ddots& -\frac{c}{2} \\
-\frac{c}{2} & 0 & \dots & 0 & \frac{c}{2} & 0
\end{pmatrix}
```

Note that $A$ is skew-symmetric, i.e.  it satisfies $A^T=-A$.
We can consider the discrete $L^2$ norm:

$$\lVert u \rVert^2_{2}=\Delta x u^Tu$$

$$\frac{d}{dt}\lVert u \rVert^2_{2}=2u^T \dot{u}=u^TAu-u^TA^Tu=0$$

i.e. it is also a conserved quantity.

Now consider the fully discretised system in time using the leapfrog integrator. Let $\mu=\frac{c\Delta t}{\Delta x}$.

This preserves the symplectic structure above (as long as the CFL condition is satisfied):

```math
\begin{aligned}
u^{n+1} - u^{n-1} &= -\mu Au^n \\
(u^{n+1} + u^{n-1})^T(u^{n+1} - u^{n-1})&= -\mu(u^{n+1} + u^{n-1})^TAu^n \\
\lVert u^{n+1} \rVert ^2_{2} - \lVert u^{n-1} \rVert ^2_{2} + \mu(u^{n+1} + u^{n-1})^TAu^n &= 0
\end{aligned}
```

One can show that when $\mu = 1$, we have the exact preservation of the discrete 2-norm of our solution in time increments of $2 \Delta t$.
$$\lVert u^{n+1} \rVert ^2_{2} = \lVert u^{n-1} \rVert ^2_{2}$$
If our first step also respects this, then we have preservation for all time. 


However we note that Leapfrog integration conserves the energy of a perturbed system; when $\mu < 1$ we have oscillatory but bounded energy, and for $\mu > 1$, the energy grows unboundedly.


