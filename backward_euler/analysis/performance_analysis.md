## Performance Analysis

We now examine the computational cost of solving the linear system

$$
A \mathbf{u}^{n+1} = \mathbf{u}^n
$$

at each time step. Recall that $A$ is an $M \times M$ circulant matrix that remains constant throughout the simulation. The total number of time steps is

$$
N_t = \frac{t_{\max}}{\Delta t} = \frac{t_{\max} \cdot a}{\Delta x \cdot c} = \frac{t_{\max} \cdot a \cdot M}{L \cdot c},
$$

With our defaults ($M = 1000$, $c = 1$, $a = 1$), we have $N_t = 1000$.

---

### Direct Solve (LU Decomposition)

The call `np.linalg.solve(A, u)` performs a full LU factorisation of $A$ at every time step. The cost of LU factorisation for a dense $M \times M$ matrix is

$$
\frac{2}{3}M^3 + \mathcal{O}(M^2)
$$

flops, and the subsequent forward–backward substitution costs $2M^2$ flops. Since we repeat this $N_t$ times, the total work is

$$
W_{\text{direct}} = N_t \left(\frac{2}{3}M^3 + 2M^2\right) = \mathcal{O}(N_t M^3).
$$

<!-- Each step also reads and writes the full $M \times M$ matrix, incurring $\sim 24M^2$ bytes of memory traffic per step. The arithmetic intensity is therefore

$$
I_{\text{direct}} = \frac{(2/3)M^3}{24M^2} \approx \frac{M}{36} \approx 28 \;\text{flops/byte},
$$

which is compute-bound, but all that compute is wasted re-factorising the same matrix. -->

**Measured time: 22.7 s.**

<!-- ---

### Pre-factored LU

The natural improvement: factorise $A$ once at the start, then use only the forward–backward substitution at each step.

The one-time factorisation costs $(2/3)M^3$ flops. Each subsequent step requires only the triangular solves at $2M^2$ flops, giving total work

$$
W_{\text{LU}} = \frac{2}{3}M^3 + N_t \cdot 2M^2 = \mathcal{O}(M^3 + N_t M^2).
$$

Per step, we read the two triangular factors ($\sim 8M^2$ bytes each) and the right-hand side vector ($8M$ bytes), so the memory traffic is $\sim 16M^2$ bytes. The arithmetic intensity drops to

$$
I_{\text{LU}} = \frac{2M^2}{16M^2} \approx 0.13 \;\text{flops/byte},
$$

which is memory-bound. The compute cost dropped by a factor of $M/3$, but the memory access pattern barely changed.

In practice, the theoretical speedup of $\sim M/3 \approx 333\times$ is not fully realised because LAPACK's LU uses highly optimised BLAS-3 routines (cache-friendly block operations), so the per-step factorisation in the direct solver is faster than raw flop counts suggest.

**Measured time: 0.34 s — a 68× speedup over direct.** -->

---

### FFT Solver (Circulant Diagonalisation)

Since $A$ is circulant, its eigenvectors are the discrete Fourier modes. Writing $A = F^{-1} \Lambda F$, where $F$ is the DFT matrix and $\Lambda = \text{diag}(\hat{a}_0, \dots, \hat{a}_{M-1})$ contains the eigenvalues, the system $A \mathbf{u}^{n+1} = \mathbf{u}^n$ becomes

$$
\mathbf{u}^{n+1} = F^{-1} \Lambda^{-1} F \mathbf{u}^n.
$$

#### Complexity of the FFT (Master Method)

The Cooley–Tukey algorithm computes the DFT of a length-$M$ sequence (assuming $M = 2^k$) by splitting it into its even- and odd-indexed entries. The DFT of length $M$ is

$$
X_j = \sum_{m=0}^{M-1} x_m \, \omega^{jm}, \quad \omega = e^{-2\pi i / M}.
$$

Separating even and odd indices $x_m = x_{2r}$ and $x_{2r+1}$:

$$
X_j = \underbrace{\sum_{r=0}^{M/2 - 1} x_{2r} \, \omega^{2rj}}_{\text{DFT of even entries}} + \omega^j \underbrace{\sum_{r=0}^{M/2 - 1} x_{2r+1} \, \omega^{2rj}}_{\text{DFT of odd entries}}.
$$

Each inner sum is itself a DFT of length $M/2$. The $M$ multiplications by the "twiddle factors" $\omega^j$ and the $M$ additions to combine the two half-transforms constitute $\mathcal{O}(M)$ work. This gives the recurrence

$$
T(M) = 2\,T\!\left(\frac{M}{2}\right) + \mathcal{O}(M).
$$

We apply the **Master Theorem**. For a recurrence of the form $T(M) = a\,T(M/b) + \mathcal{O}(M^d)$, the three cases are:

$$
T(M) = \begin{cases}
\mathcal{O}(M^d \log M) & \text{if } a = b^d, \\
\mathcal{O}(M^d) & \text{if } a < b^d, \\
\mathcal{O}(M^{\log_b a}) & \text{if } a > b^d.
\end{cases}
$$

Here we have $a = 2$, $b = 2$, and $d = 1$ (since the combine step is $\mathcal{O}(M^1)$). Checking:

$$
b^d = 2^1 = 2 = a,
$$

so we are in the first case ($a = b^d$), giving

$$
\boxed{T(M) = \mathcal{O}(M^d \log M) = \mathcal{O}(M \log M).}
$$

---

The eigenvalues $\hat{a}_k = \text{FFT}(A_{:,0})$ are computed once. Each time step then requires:

1. **FFT** of $\mathbf{u}^n$: $5M \log_2 M$ flops (Cooley–Tukey, $\sim 5$ flops per butterfly),
2. **Pointwise division** $\hat{u}_k / \hat{a}_k$: $6M$ flops (complex division),
3. **IFFT**: another $5M \log_2 M$ flops,
4. Extracting the real part: $M$ flops.

The total per-step cost is

$$
W_{\text{FFT}}^{\text{step}} = 10 M \log_2 M + 7M \approx 1.07 \times 10^5 \;\text{flops},
$$

and the total over all steps is

$$
W_{\text{FFT}} = N_t (10 M \log_2 M + 7M) = \mathcal{O}(N_t M \log M).
$$

<!-- The memory traffic per step involves a few length-$M$ complex vectors, totalling $\sim 80M$ bytes. The arithmetic intensity is

$$
I_{\text{FFT}} = \frac{10 M \log_2 M}{80 M} = \frac{\log_2 M}{8} \approx 1.25 \;\text{flops/byte}.
$$

This is memory-bound on most hardware, but the total data volume ($\sim 80$ KB) comfortably fits in L2 cache, so the memory bottleneck is much less severe than it appears from the intensity alone. -->

**Measured time: 0.025 s, a 907× speedup over direct**

---

### Summary

As $M$ grows, the asymptotic scaling is:

- **Direct:** $\mathcal{O}(N_t M^3)$ — impractical for large grids.
<!-- - **Pre-factored LU:** $\mathcal{O}(M^3 + N_t M^2)$ — the one-time factorisation dominates for moderate $N_t$, per-step cost dominates for large $N_t$. -->
- **FFT:** $\mathcal{O}(N_t M \log M)$ — optimal for circulant systems, but requires the circulant structure.

<!-- The FFT approach also has the best memory scaling at $\mathcal{O}(M)$ per step. -->