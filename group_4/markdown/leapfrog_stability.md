# Von Neumann Stability

Write the leapfrog scheme for the advection equation as

$$
u_{k}^{n+1} = u_{k}^{n-1} - c(u_{k+1}^n - u_{k-1}^n)
$$

Use the ansatz $u_k^n = A^n e^{i \xi k}$ where $\xi = \pi \Delta x$ and substitute:

$$
A^{n+1} e^{i \xi k} 
= A^{n-1} e^{i \xi k} - c\left(A^n e^{i\xi(k+1)} - A^n e^{i\xi(k-1)}\right)
$$

Cancel $A^{n-1}e^{i\xi k}$:

$$
A^2 = 1 - cA(e^{i\xi} - e^{-i\xi})
$$

Using $e^{i\xi} - e^{-i\xi} = 2i\sin\xi$:

$$
A^2 + 2ic\sin(\xi)A - 1 = 0
$$

Solving:

$$
A = -ic\sin(\xi) \pm \sqrt{1 - c^2\sin^2(\xi)}
$$

---

## Case 1: $|c| \le 1$

Then the square root is real, and

$$
|A|^2 = c^2\sin^2(\xi) + (1 - c^2\sin^2(\xi)) = 1
$$

Hence $|A| = 1$ and the scheme is stable.

---

## Case 2: $|c| > 1$

For some $\xi$, we have $c^2\sin^2(\xi) > 1$, and the roots become

$$
A_\pm = i\left(-c\sin\xi \pm \sqrt{c^2\sin^2\xi - 1}\right).
$$

Their product is

$$
A_+ A_- = 1.
$$

Therefore

$$
|A_+||A_-| = 1.
$$

Since $|A_+| \ne |A_-|$, one must satisfy $|A|>1$.

Thus the scheme is unstable.

---

## Stability condition

$$
|c| \le 1.
$$
