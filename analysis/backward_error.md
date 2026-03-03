# Backward Error Analysis

We consider the linear advection equation

$$
u_t + au_x = 0
$$

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

Now, substituting in to the discretised equation, we find

$$
\frac{u - \left( u - \Delta t u_t + \frac{\Delta t^2}{2}u_{tt} + \mathcal{O}(\Delta t^3) \right)}{\Delta t} + a \frac{u - \left( u - \Delta x u_x + \frac{\Delta x^2}{2}u_{xx} + \mathcal{O}(\Delta x^3) \right)}{\Delta x} = 0.
$$

Simplifying, we have that 

$$
u_t + au_x - \frac{\Delta t}{2}u_{tt} - a\frac{\Delta t}{2}u_{xx} = \mathcal{O}(\Delta t^2, \Delta x^2)
$$

We now aim to eliminate the time derivative $u_{tt}$.

$$
\partial_t : \ \ u_{tt} = -au_{xt} + \frac{\Delta t}{2}u_{ttt} + a\frac{\Delta t}{2}u_{xxt} + \mathcal{O}(\Delta t^2, \Delta x^2) 
$$

$$
\partial_x : \ \ u_{tx} = -au_{xx} + \frac{\Delta t}{2}u_{ttx} + a\frac{\Delta t}{2}u_{xxx} + \mathcal{O}(\Delta t^2, \Delta x^2)
$$

so, we can find

$$
u_{tt} = a^2u_{xx} + \mathcal{O}(\Delta t^2, \Delta x^2)
$$

and therefore,

$$
u_t + au_x = \frac{\Delta t}{2](a^2 + a)u_{xx} + \mathcal{O}(\Delta t^2, \Delta x^2)
$$

as a positive coefficient in front of the even derivative means that we have diffusion, we can see that this is always satisfied as the method is upwind for positive $a$. Therefore, the backward upwind Euler method is unconditionally stable.
