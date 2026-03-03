# Backward Error Analysis

We consider the linear advection equation

$$
u_t + au_x = 0
$$

We aim to conduct a Wilkinson's backward error analysis by observing for which $a$ and Courant number $c$, the Upwind Backward Euler method is stable.

We are using a backward Euler, hence

$$
u_t \approx \frac{u_k^{n+1} - u_k^n}{\Delta t}
$$

and an upwind method with the backward Euler, hence

$$
u_x \approx \frac{u_k^{n+1} - u_{k-1}^{n+1}}{\Delta x}.
$$

Substituting into the linear advection equation,

$$
\frac{u_k^{n+1} - u_k^n}{\Delta t} + a\frac{u_k^{n+1} - u_{k-1}^{n+1}}{\Delta x} = 0
$$

If we take, for simplicity, $u := u_k^{n+1} $, then we can Taylor expand the terms

$$
u_k^{n} = u - \Delta t u_t + \frac{\Delta t^2}{2}u_{tt} + \mathcal{O}(\Delta t^3)
$$

and

$$
u_{k-1}^{n+1} = u - \Delta x u_x + \frac{\Delta x^2}{2}u_{xx} + \mathcal{O}(\Delta x^3).
$$

