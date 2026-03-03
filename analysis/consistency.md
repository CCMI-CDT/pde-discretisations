## 1D Advection Equation with Periodic B.C.s

We consider the PDE

$$
\begin{cases}
u_t + a u_x = 0, \\
u(x,0) = \varphi(x), \\
u(L,t) = u(0,t).
\end{cases}
$$

One can readily check that

$$
u(x,t) = \varphi(x - at)
$$

is a solution, provided $\varphi$ is sufficiently regular.

---

## Backward Euler–Centered Method

This is a finite difference method. We approximate

$$
u_m^n \approx u(m\Delta x,n\Delta t),
$$

with spatial step

$$
\Delta x = \frac{L}{M}, \quad m = 0,\dots, M-1,
$$

and time levels $n = 0,1,2,\dots$.

Note that periodic boundary conditions are enforced by identifying indices modulo $M$:

$$
u_M^n = u_0^n.
$$

Since we only store indices up to $0, \dots, M-1$, this is already dealt with in the discretized equation.

---

### Discretisation

We use:

**Time discretisation:**

$$
u_t \approx \frac{u_m^{n+1} - u_m^n}{\Delta t}.
$$

**Centered difference (space one timestep ahead, implicit):**

$$
u_x \approx \frac{u_{m+1}^{n+1} - u_{m-1}^{n+1}}{2\Delta x}.
$$

Substituting into the PDE $u_t + a u_x = 0$ gives

$$
\frac{u_m^{n+1} - u_m^n}{\Delta t} + a \frac{u_{m+1}^{n+1} - u_{m-1}^{n+1}}{2\Delta x} = 0.
$$

Multiplying through by $\Delta t$:

$$u_m^{n+1} - u_m^n+ \frac{a\Delta t}{2\Delta x}\left(u_{m+1}^{n+1} - u_{m-1}^{n+1}\right)= 0.$$

Let the Courant number be

$$
c = \frac{a\Delta t}{\Delta x}.
$$

We arrive at the implicit scheme

$$
u_m^{n+1}+ \frac{c}{2}\left(u_{m+1}^{n+1} - u_{m-1}^{\,n+1}\right)= u_m^n.
$$

---

### Consistency Analysis

Using Taylor expansion in space (at fixed time), we write

$$
u(x + \Delta x, \cdot)= u(x, \cdot)+ \Delta xu_x+ \frac{(\Delta x)^2}{2} u_{xx}+ \frac{(\Delta x)^3}{6} u_{xxx}+ \cdots,
$$

and similarly

$$u(x - \Delta x, \cdot)= u(x, \cdot)- \Delta xu_x+ \frac{(\Delta x)^2}{2} u_{xx}- \frac{(\Delta x)^3}{6} u_{xxx}+ \cdots.
$$

Some rearranging gives

$$
2\Delta xu_x = u(x+\Delta x,\cdot) - u(x-\Delta x,\cdot) + \frac{2(\Delta x)^3}{6} u_{xxx}+ \cdots 
$$

Dividing by $2\Delta x$,

$$
u_x  = \frac{u(x+\Delta x,\cdot) - u(x-\Delta x,\cdot)}{2\Delta x}+  \mathcal{O}((\Delta x)^2).
$$

We see the leading truncation error is of order $(\Delta x)^2$, an order better than standard forward-backward methods.

Let $u(x, t)$ solve $u_t + au_x=0$ and set $u_m^n = u(m \Delta x, n \Delta t)$. Define the residual as

$$
\tau_m^n := \left| u_t(m \Delta x, n \Delta t) + au_x(m \Delta x, n \Delta t) - \left( \frac{u_m^{n+1} - u_m^n}{\Delta t} + a\frac{u_{m+1}^{n+1} - u_{m-1}^{n+1}}{2 \Delta x} \right) \right|
$$

By the triangle inequality, we can bound

$$
\tau_m^n \le \left| u_t(m \Delta x, n \Delta t) - \frac{u_m^{n+1} - u_m^n}{\Delta t} \right| + a \left| u_x(m \Delta x, n \Delta t) - \frac{u_{m+1}^{n+1} - u_{m-1}^{n+1}}{2 \Delta x} \right|
$$

Taylor in time and our above centered difference approximation:

$$
\frac{u_m^{n+1} - u_m^n}{\Delta t} = u_t(m \Delta x, n \Delta t) + \mathcal{O}(\Delta t)
$$
$$
\frac{u_{m+1}^{n+1} - u_{m-1}^{n+1}}{2 \Delta x} = u_x(m \Delta x, n \Delta t) + \mathcal{O}((\Delta x)^2) + \mathcal{O}(\Delta t)
$$

and therefore the residual satisfies

$$
\tau_m^n = \mathcal{O}(\Delta t) + \mathcal{O}((\Delta x)^2).
$$

Hence, as $\Delta t, \Delta x \to 0$, the **Backward Euler–Centered scheme is consistent** with the 1D advection equation, with first-order accuracy in time and second-order accuracy in space.
