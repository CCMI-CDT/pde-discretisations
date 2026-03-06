## Kuramoto–Sivashinsky equation

The Kuramoto–Sivashinsky equation takes the following form:

$$
u_t + u u_x + u_{xx} + u_{xxxx} = 0
$$

Make the substitution $w = u_{xx}$, which gives a coupled system of equations:

$$
\begin{aligned}
u_t + u u_x + w + w_{xx} &= 0 \\
w - u_{xx} &= 0
\end{aligned}
$$

Assuming $L$-periodic boundary conditions and using the integration by parts formula, with test functions $v, s \in \mathcal{H}^1(0,L)$,

$$
\int_\Omega v\, u_{xx}
= \int_{\partial\Omega} v\,u - \int_\Omega v_x u_x
= -\int_\Omega v_x u_x
$$

This gives the variational form

$$
\begin{aligned}
\int_\Omega v\, u_t
+ \int_\Omega v\, u u_x
- \int_\Omega v_x u_x
- \int_\Omega v_x w_x &= 0 \\
\int_\Omega s\, w
+ \int_\Omega s_x u_x &= 0
\end{aligned}
$$

For time stepping, we use a backward Euler discretisation for $u_t$:

$$
u_t \approx \frac{u^{n+1} - u^n}{\Delta t}.
$$

This gives

$$
\begin{aligned}
\int_\Omega v\, \frac{u^{n+1} - u^n}{\Delta t}
+ \int_\Omega v\, u^n u_x^n
- \int_\Omega v_x u_x^n
- \int_\Omega v_x w_x^n &= 0 \\
\int_\Omega s\, w^n
+ \int_\Omega s_x u_x^n &= 0
\end{aligned}
$$
