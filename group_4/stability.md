# Von Neumann Stability

  
Write the leapfrog scheme for the advection equation as
$$

\begin{align*}
u_{k}^{n+1} &= u_{k}^{n-1} - c(u_{k+1}^n - u_{k-1}^n)
\end{align*}

$$
Use the ansatz $u_{k}^n = A^ne^{i \xi k}$ with $\xi = \pi\Delta x$ and substitute into the above:
$$

\begin{align*}

A^{n+1} e^{i \xi k} &= A^{n-1}e^{i\xi k} - c(A^n e^{i\xi(k+1)} - A^n e^{i`\xi (k-1)}) \\

A^2 &= 1 - cA(e^{i\xi} - e^{i\xi}) \\

A^2 + 2ic\sin(\xi)A -1 &= 0 \\

\implies A &= \frac{-2ic\sin(\xi) \pm \sqrt{-4c^2\sin^2(\xi) + 4} }{2} \\

A &= -ic\sin(\xi) \pm \sqrt{1- c^2 \sin^2(\xi)}

\end{align*}
$$

In the case $\lvert c \rvert \leq1$, the square root is real and $|A| = c^2\sin^2(\xi) + 1 - c^2\sin^2(\xi) = 1$. Hence the solution is stable.
In the case $|c| > 1$, then the roots are purely imaginary. Let $A_{+} = i\left(-c\sin(\xi) + \sqrt{ c^2 \sin^2(\xi) - 1 }\right)$ and $A_{+} = i\left(-c\sin(\xi) - \sqrt{ c^2 \sin^2(\xi) - 1 }\right)$, then

$$

\begin{align*}

A_{+}A_{-} &= - (c^2 \sin^2(\xi) - c^2\sin^2(\xi) - 1) = 1 \\

&= 1 \\

\implies |A_{+}||A_{-}| &= 1

\end{align*}

$$

so either $|A_{+}| > 1$ or $|A_{-}| > 1$ in which case the solution is unstable.

In summary we require $|c| \leq 1$ for stability.