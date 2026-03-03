```math
\begin{aligned}
\frac{1}{2\Delta t}(u^{n+1}_k - u^{n-1}_k)
+ \frac{a}{2\Delta x}(u_{k+1}^{n} - u_{k-1}^n)
&= \frac{1}{2\Delta t}\Bigl(
u^{n}_k + \Delta t\, u_{t}\big|_k^{n}
+ \frac{\Delta t^2}{2} u_{tt}\big|_k^{n}
+ \mathcal{O}((\Delta t)^3) \\
&\qquad\quad
- u^{n}_k + \Delta t\, u_{t}\big|_k^{n}
- \frac{\Delta t^2}{2} u_{tt}\big|_k^{n}
+ \mathcal{O}((\Delta t)^3)
\Bigr) \\
&\quad + \frac{a}{2\Delta x}\Bigl(
u^{n}_k + \Delta x\, u_{x}\big|_k^{n}
+ \frac{\Delta x^2}{2} u_{xx}\big|_k^{n}
+ \mathcal{O}((\Delta x)^3) \\
&\qquad\quad
- u^{n}_k + \Delta x\, u_{x}\big|_k^{n}
- \frac{\Delta x^2}{2} u_{xx}\big|_k^{n}
+ \mathcal{O}((\Delta x)^3)
\Bigr) \\
&= \frac{1}{2\Delta t}\Bigl(
2\Delta t\, u_t\big|_k^{n}
+ \mathcal{O}((\Delta t)^3)
\Bigr)
+ \frac{a}{2\Delta x}\Bigl(
2\Delta x\, u_x\big|_k^{n}
+ \mathcal{O}((\Delta x)^3)
\Bigr) \\
&= \bigl(u_t\big|_k^{n} + a\, u_x\big|_k^{n}\bigr)
+ \mathcal{O}((\Delta t)^2,\, (\Delta x)^2)
\end{aligned}
```

This residual $\to 0$ as $\Delta t, \Delta x \to 0$, so the scheme is consistent.
