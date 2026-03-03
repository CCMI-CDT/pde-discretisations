$$
\begin{aligned}
\frac{1}{2\Delta t}(u^{n+1}_k - u^{n-1}_k)
+ \frac{a}{2\Delta x}(u_{k+1}^{n} - u_{k-1}^n)
&= \frac{1}{2\Delta t}\Big(
u^{n}_k + \Delta t\, u_{t}\big|_k^{n}
+ \frac{\Delta t^2}{2} u_{tt}\big|_k^{n}
+ \mathcal{O}((\Delta t)^3) \\
&\qquad
- u^{n}_k + \Delta t\, u_{t}\big|_k^{n}
- \frac{\Delta t^2}{2} u_{tt}\big|_k^{n}
+ \mathcal{O}((\Delta t)^3)
\Big) \\
&\quad + \frac{a}{2\Delta x}\Big(
u^{n}_k + \Delta x\, u_{x}\big|_k^{n}
+ \frac{\Delta x^2}{2} u_{xx}\big|_k^{n}
+ \mathcal{O}((\Delta x)^3) \\
&\qquad
- u^{n}_k + \Delta x\, u_{x}\big|_k^{n}
- \frac{\Delta x^2}{2} u_{xx}\big|_k^{n}
+ \mathcal{O}((\Delta x)^3)
\Big) \\
&= \frac{1}{2\Delta t}\Big(
2 \Delta t\, u_t\big|_k^{n}
+ \mathcal{O}((\Delta t)^3)
\Big)
+ \frac{a}{2\Delta x}\Big(
2 \Delta x\, u_x\big|_k^{n}
+ \mathcal{O}((\Delta x)^3)
\Big) \\
&=\big(u_t\big|_k^{n} + a u_x\big|_k^{n}\big)
+ \mathcal{O}((\Delta t)^2, (\Delta x)^2)
\end{aligned}
$$

This residual $\to 0$ as $\Delta t, \Delta x \to 0$, so the scheme is consistent.