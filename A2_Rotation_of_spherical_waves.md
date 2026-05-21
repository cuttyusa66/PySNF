# Appendix A2 — Rotation of Spherical Waves

**Source:** Hansen, *Spherical Near-Field Antenna Measurements*, Appendix A2 (pp. 343–354).

**Cross-references:**
- Symbols: [A5_List_of_principal_symbols_and_uses.md](A5_List_of_principal_symbols_and_uses.md)
- Chapter 2 (Dipole rotation example uses these): [02_Scattering_matrix_description_of_an_antenna.md](02_Scattering_matrix_description_of_an_antenna.md)
- Appendix A1 (wave functions, $\bar{P}_n^{|m|}$, $\vec{K}_{smn}$): [A1_Spherical_wave_functions_notation_and_properties.md](A1_Spherical_wave_functions_notation_and_properties.md)

> **Purpose.** Computational reference for spherical-wave rotation. Every formula is tagged with its book number `(A2.N)`. Includes the explicit low-order $\Delta^n_{m'm}$ tables for $n = 0\dots 5$ — drop-in unit-test targets.

---

## 1. Euler-Angle Convention (§A2.1)

Two right-handed Cartesian frames, initially coincident. The **primed** frame $(x', y', z')$ is rotated relative to the fixed **unprimed** frame $(x, y, z)$ through three successive rotations *about its own axes*:

| Step | Rotation axis | Angle | Frame after step |
|---|---|---|---|
| 1 | $z$-axis | $\phi_o$ | $(x_1, y_1, z_1)$ |
| 2 | $y_1$-axis | $\theta_o$ | $(x_2, y_2, z_2)$ |
| 3 | $z_2$-axis | $\chi_o$ | $(x', y', z')$ |

**The Euler-angle triple is written in reversed order:** $(\chi_o,\,\theta_o,\,\phi_o)$.

Right-hand rule: a positive $\pi/2$ rotation about $\hat{z}$ takes $\hat{x} \to \hat{y}$.

> **Bug-bait.** The order of rotation operations is $\phi_o \to \theta_o \to \chi_o$ (z-then-y-then-z), but the **printed argument tuple** is $(\chi_o, \theta_o, \phi_o)$. Mismatched conventions are the most common source of rotation errors. **Validate at every entry point.**

> **vs. scipy.** `scipy.spatial.transform.Rotation.from_euler("ZYZ", [α, β, γ])` defaults to **intrinsic** rotations in the listed order `α, β, γ`. To match Hansen, call with `[phi_o, theta_o, chi_o]` (i.e., reverse the printed tuple) and use intrinsic mode.

---

## 2. Rotation of a Spherical Wave Function (§A2.2)

A wave function defined in the *unprimed* frame is expanded in the *primed* frame's basis:

$$
\boxed{\;\vec{F}^{(c)}_{smn}(r,\theta,\phi) = \sum_{\mu=-n}^{n} D^{n}_{\mu m}(\chi_o, \theta_o, \phi_o)\,\vec{F}^{(c)}_{s\mu n}(r', \theta', \phi')\;} \tag{A2.1}
$$

with the rotation function decomposing as:

$$
\boxed{\;D^{n}_{\mu m}(\chi_o, \theta_o, \phi_o) = e^{im\phi_o}\,d^{n}_{\mu m}(\theta_o)\,e^{i\mu\chi_o}\;} \tag{A2.2}
$$

### Invariants under rotation

| Quantity | Behaviour |
|---|---|
| $s$ | **Invariant.** TE stays TE, TM stays TM. |
| $n$ | **Invariant.** Modes of different $n$ don't mix. |
| $c$ | Invariant. Radial dependence is preserved. |
| $m$ | **Varies** (sum over $\mu$). |
| $r$ | Invariant ($r' = r$). |

### Phase-only special cases

- $\phi_o$-rotation alone ($\theta_o = \chi_o = 0$): multiplies coefficients by $e^{im\phi_o}$.
- $\chi_o$-rotation alone ($\theta_o = \phi_o = 0$): multiplies by $e^{i\mu\chi_o}$.
- The interesting work is in $d^n_{\mu m}(\theta_o)$.

### Rotation of mode coefficients

If the unrotated antenna has coefficients $Q^{(c)}_{smn}$, then after rotating the antenna by Euler angles $(\chi_o, \theta_o, \phi_o)$ the new coefficients (in the original frame) are:

$$
Q'^{(c)}_{s\mu n} = \sum_{m=-n}^{n} D^{n}_{\mu m}(\chi_o, \theta_o, \phi_o)\,Q^{(c)}_{smn}
$$

(Derived by transforming the expansion through (A2.1) and reading off coefficients of $\vec{F}^{(c)}_{s\mu n}$ in the new frame.)

> **Quick sanity check** (Chapter 2 dipole example, Eq. 2.123): rotating a $\hat{z}$-electric dipole into a $\hat{x}$-electric dipole uses $(\chi_o, \theta_o, \phi_o) = (0, -\pi/2, 0)$. With $m = 0$ (only $Q_{201}$ is non-zero for the z-dipole), the sum collapses to $\mu = \pm 1$:
> $$\vec{F}^{(3)}_{201}(r,\theta,\phi) = d^1_{-1,0}(-\pi/2)\,\vec{F}^{(3)}_{2,-1,1}(r',\theta',\phi') + d^1_{0,0}(-\pi/2)\,\vec{F}^{(3)}_{201}(r',\theta',\phi') + d^1_{1,0}(-\pi/2)\,\vec{F}^{(3)}_{211}(r',\theta',\phi')$$
> $$= \tfrac{\sqrt{2}}{2}\,\vec{F}^{(3)}_{2,-1,1} + 0 - \tfrac{\sqrt{2}}{2}\,\vec{F}^{(3)}_{211}$$
> matching Eq. (2.123). Use this as a regression test.

---

## 3. The Rotation Coefficient $d^n_{\mu m}(\theta)$ (§A2.3)

### General expression (Edmonds [1], Eq. 4.1.15)

$$
\boxed{\;d^{n}_{\mu m}(\theta) = \sqrt{\frac{(n+\mu)!\,(n-\mu)!}{(n+m)!\,(n-m)!}}\,\sum_{\sigma}\binom{n+m}{n-\mu-\sigma}\binom{n-m}{\sigma}\,(-1)^{n-\mu-\sigma}\,\Bigl(\cos\tfrac{\theta}{2}\Bigr)^{2\sigma+\mu+m}\,\Bigl(\sin\tfrac{\theta}{2}\Bigr)^{2n-2\sigma-\mu-m}\;} \tag{A2.3}
$$

with the binomial coefficient $\binom{i}{j} = i!/((i-j)!\,j!)$ (A2.4).

Summation range for $\sigma$: only those terms where every factorial argument is non-negative (i.e. $\sigma \ge 0$, $\sigma \le n-m$, $\sigma \le n-\mu$, $\sigma \ge -\mu-m$). $d^n_{\mu m}(\theta)$ is **real**.

### Jacobi-polynomial form (Edmonds [1], Eq. 4.1.23)

$$
d^{n}_{\mu m}(\theta) = \sqrt{\frac{(n+\mu)!\,(n-\mu)!}{(n+m)!\,(n-m)!}}\,\Bigl(\cos\tfrac{\theta}{2}\Bigr)^{\mu+m}\,\Bigl(\sin\tfrac{\theta}{2}\Bigr)^{\mu-m}\,P^{(\mu-m,\mu+m)}_{n-\mu}(\cos\theta) \tag{A2.5}
$$

where $P^{(\alpha,\beta)}_n$ is the Jacobi polynomial.

### Normalization

$$
\sum_{\mu=-n}^{n}\bigl(d^{n}_{\mu m}(\theta)\bigr)^2 = 1 \quad\text{for all }(m,n) \tag{A2.6}
$$

### Symmetries

$$
d^{n}_{\mu m}(\theta) = d^{n}_{m\mu}(-\theta) \tag{A2.7}
$$

$$
d^{n}_{\mu m}(\theta) = (-1)^{\mu+m}\,d^{n}_{m\mu}(\theta) \tag{A2.8}
$$

$$
d^{n}_{\mu m}(\theta) = (-1)^{\mu+m}\,d^{n}_{-\mu,-m}(\theta) \tag{A2.9}
$$

Combining (A2.8) and (A2.9): $d^n_{\mu m}(\theta)$ has **parity $(\mu + m)$**. It is $2\pi$-periodic in $\theta$.

### Orthogonality (Edmonds [1], Eq. 4.6.1)

$$
\int_{0}^{\pi} d^{n}_{\mu m}(\theta)\,d^{n'}_{\mu m}(\theta)\,\sin\theta\,d\theta = \frac{2}{2n+1}\,\delta_{nn'} \tag{A2.10}
$$

### Fourier expansion (used in Wacker algorithm)

$$
d^{n}_{\mu m}(\theta) = i^{\mu - m}\sum_{m'=-n}^{n}\Delta^{n}_{m'\mu}\,\Delta^{n}_{m'm}\,e^{-im'\theta} \tag{A2.11}
$$

$$
= i^{m - \mu}\sum_{m'=-n}^{n}\Delta^{n}_{m'\mu}\,\Delta^{n}_{m'm}\,e^{+im'\theta} \tag{A2.12}
$$

The Fourier coefficient is a product of "deltas":

$$
\boxed{\;\Delta^{n}_{m'm} \equiv d^{n}_{m'm}(\pi/2)\;} \tag{A2.13}
$$

### Three-term recurrence in $\mu$ (Fano & Racah [3])

$$
\sqrt{(n+\mu+1)(n-\mu)}\,\sin\theta\,d^{n}_{\mu+1, m}(\theta) + \sqrt{(n+\mu)(n-\mu+1)}\,\sin\theta\,d^{n}_{\mu-1, m}(\theta) + (2m - 2\mu\cos\theta)\,d^{n}_{\mu m}(\theta) = 0 \tag{A2.14}
$$

### Special values

$$
d^{n}_{\mu m}(0) = \delta_{\mu m},\quad\text{all }n \tag{A2.15}
$$

$$
d^{n}_{\mu m}(\pi) = (-1)^{n+m}\,\delta_{\mu,-m} \tag{A2.16}
$$

### Special cases (Larsen [4])

Relation to normalized associated Legendre functions:

$$
\boxed{\;d^{n}_{0m}(\theta) = \Bigl(-\frac{m}{|m|}\Bigr)^m\,\sqrt{\frac{2}{2n+1}}\,\bar{P}_n^{|m|}(\cos\theta)\;} \tag{A2.17}
$$

(Convention: $(-m/|m|)^m = 1$ when $m=0$.)

Combinations relevant for $m = \pm 1$:

$$
d^{n}_{1m}(\theta) + d^{n}_{-1,m}(\theta) = -\frac{2}{\sqrt{n(n+1)}}\,\frac{m\,d^{n}_{0m}(\theta)}{\sin\theta} \tag{A2.18}
$$

$$
d^{n}_{1m}(\theta) - d^{n}_{-1,m}(\theta) = -\frac{2}{\sqrt{n(n+1)}}\,\frac{d}{d\theta}\{d^{n}_{0m}(\theta)\} \tag{A2.19}
$$

### Direct link to far-field pattern components $\vec{K}_{smn}$

Inserting (A2.17) into (A2.18)/(A2.19) and using the $\vec{K}_{smn}$ definitions of (A1.59)/(A1.60):

$$
d^{n}_{1m}(\theta) + d^{n}_{-1,m}(\theta) = \frac{-2i^{n}}{\sqrt{2n+1}}\,e^{-im\phi}\,\{\vec{K}_{1mn}(\theta,\phi)\}_\theta \tag{A2.20}
$$

$$
= \frac{2i^{n+1}}{\sqrt{2n+1}}\,e^{-im\phi}\,\{\vec{K}_{2mn}(\theta,\phi)\}_\phi \tag{A2.21}
$$

$$
d^{n}_{1m}(\theta) - d^{n}_{-1,m}(\theta) = \frac{2i^{n+1}}{\sqrt{2n+1}}\,e^{-im\phi}\,\{\vec{K}_{1mn}(\theta,\phi)\}_\phi \tag{A2.22}
$$

$$
= \frac{-2i^{n}}{\sqrt{2n+1}}\,e^{-im\phi}\,\{\vec{K}_{2mn}(\theta,\phi)\}_\theta \tag{A2.23}
$$

> **Consistency check.** (A2.20) and (A2.21) equating yields $\{\vec{K}_{2mn}\}_\phi = i\,\{\vec{K}_{1mn}\}_\theta$, which is one component of the helicity identity $\vec{K}_{smn} = i\hat{r}\times\vec{K}_{3-s,m,n}$ (Eq. 2.178). Use as a unit test.

---

## 4. The Deltas $\Delta^{n}_{m'm}$ (§A2.4)

$\Delta^{n}_{m'm} = d^{n}_{m'm}(\pi/2)$ for $-n \le m', m \le n$. Indexed inside a "delta pyramid" with $0 \le n \le N$.

### General expressions

From (A2.3) with $\theta = \pi/2$:

$$
\Delta^{n}_{m'm} = \sqrt{\frac{(n+m')!\,(n-m')!}{(n+m)!\,(n-m)!}}\,\frac{1}{2^n}\sum_{\sigma}\binom{n+m}{n-m'-\sigma}\binom{n-m}{\sigma}\,(-1)^{n-m'-\sigma} \tag{A2.24}
$$

From the Jacobi form (A2.5):

$$
\Delta^{n}_{m'm} = \sqrt{\frac{(n+m')!\,(n-m')!}{(n+m)!\,(n-m)!}}\,\frac{1}{2^{m'}}\,P^{(m'-m,\,m'+m)}_{n-m'}(0) \tag{A2.25}
$$

### Symmetries (full triangle from one octant)

| Eq. | Symmetry |
|---|---|
| (A2.26) | $\Delta^{n}_{m'm} = (-1)^{m'+m}\,\Delta^{n}_{mm'}$ |
| (A2.27) | $\Delta^{n}_{m'm} = (-1)^{n+m'}\,\Delta^{n}_{m,-m'}$ |
| (A2.28) | $\Delta^{n}_{m'm} = (-1)^{n+m}\,\Delta^{n}_{-m',m}$ |
| (A2.29) | $\Delta^{n}_{m'm} = (-1)^{m'+m}\,\Delta^{n}_{-m',-m}$ |
| (A2.30) | $\Delta^{n}_{m'm} = \Delta^{n}_{-m,-m'}$ |
| (A2.31) | $\Delta^{n}_{m'm} = (-1)^{n+m}\,\Delta^{n}_{-mm'}$ |
| (A2.32) | $\Delta^{n}_{m'm} = (-1)^{n+m'}\,\Delta^{n}_{m',-m}$ |

These reduce storage to the $(m \ge 0,\,m' \ge 0,\,m' \ge m)$ triangle.

### Recurrence relations

**Horizontal in $m'$** (from A2.14 at $\theta = \pi/2$):

$$
\sqrt{(n+m'+1)(n-m')}\,\Delta^{n}_{m'+1,m} + \sqrt{(n+m')(n-m'+1)}\,\Delta^{n}_{m'-1,m} + 2m\,\Delta^{n}_{m'm} = 0 \tag{A2.33}
$$

**Horizontal in $m$** (from A2.33 + symmetries):

$$
\sqrt{(n+m+1)(n-m)}\,\Delta^{n}_{m',m+1} + \sqrt{(n+m)(n-m+1)}\,\Delta^{n}_{m',m-1} - 2m'\,\Delta^{n}_{m'm} = 0 \tag{A2.34}
$$

**Lewis' diagonal recurrence** [6]:

$$
\frac{\sqrt{(n+m'+1)(n-m')}\,\sqrt{(n+m+1)(n-m)}}{m' + m + 1}\,\Delta^{n}_{m'+1,m+1} + \frac{\sqrt{(n+m')(n-m'+1)}\,\sqrt{(n+m)(n-m+1)}}{m' + m - 1}\,\Delta^{n}_{m'-1,m-1} = \frac{2(m'+m)}{(m'+m)^2 - 1}\bigl\{n(n+1) - (m'+m)^2 + m'm + 1\bigr\}\,\Delta^{n}_{m'm} \tag{A2.35}
$$

**Vertical (in $n$)** from Jacobi-polynomial recurrence:

$$
\sqrt{(n+m'+1)(n-m'+1)(n+m+1)(n-m+1)}\,n\,\Delta^{n+1}_{m'm} + \sqrt{(n+m')(n-m')(n+m)(n-m)}\,(n+1)\,\Delta^{n-1}_{m'm} + (2n+1)\,m'm\,\Delta^{n}_{m'm} = 0 \tag{A2.36}
$$

**Wave recurrence for Fourier coefficients of $d^n_{0m}(\theta)$, $\bar{P}_n^m$, $P_n^m$:**

$$
(n+m'+2)(n-m'-1)\,\Delta^{n}_{m'+2,m}\,\Delta^{n}_{m'+2,0} + (n+m'-1)(n-m'+2)\,\Delta^{n}_{m'-2,m}\,\Delta^{n}_{m'-2,0} - 2(n^2 - m'^2 + n - 2m^2)\,\Delta^{n}_{m'm}\,\Delta^{n}_{m'0} = 0 \tag{A2.37}
$$

### Special cases (from A2.17, A2.24, [2] eqn 8.6.1)

**Column at $m = 0$ or $m' = 0$:**

$$
\Delta^{n}_{m'0} = \begin{cases} 0, & (n+m')\text{ odd} \\ (-1)^{(n-m')/2}\,\dfrac{1}{2^n}\,\sqrt{\dbinom{n+m'}{(n+m')/2}\dbinom{n-m'}{(n-m')/2}}, & (n+m')\text{ even}\end{cases} \tag{A2.38}
$$

$$
\Delta^{n}_{0m} = \begin{cases} 0, & (n+m)\text{ odd} \\ (-1)^{(n+m)/2}\,\dfrac{1}{2^n}\,\sqrt{\dbinom{n+m}{(n+m)/2}\dbinom{n-m}{(n-m)/2}}, & (n+m)\text{ even}\end{cases} \tag{A2.39}
$$

**Diagonal at $m = n$ or $m' = n$ (corner of the layer):**

$$
\Delta^{n}_{m'n} = (-1)^{n+m'}\,\frac{1}{2^n}\,\sqrt{\dbinom{2n}{n-m'}} \tag{A2.40}
$$

$$
\Delta^{n}_{nm} = \frac{1}{2^n}\,\sqrt{\dbinom{2n}{n-m}} \tag{A2.41}
$$

$$
\boxed{\;\Delta^{n}_{nn} = \frac{1}{2^n}\;} \tag{A2.42}
$$

---

## 5. Computational Stability (§A2.5)

### Three-term recurrence theory (§A2.5.1)

For $y_{n+2} + a\,y_{n+1} + b\,y_n = 0$, the characteristic roots $K_{1,2}$ from $K^2 + aK + b = 0$ determine stability:

- **Two distinct real roots** $|K_1| < |K_2|$: the $K_1^n$ solution decreases relative to $K_2^n$. **Recurrence is unstable** for the smaller solution in the direction of increasing index; must run in the **opposite** direction (Miller-style).
- **Two complex conjugate roots**: equal-envelope oscillating real solutions — stable in **both** directions.

### Delta recurrence stability (§A2.5.2)

For the $\Delta$ recurrence (A2.33), the approximate characteristic roots for large $(n+m'), (n-m')$ are:

$$
K_{\pm} \approx -\sqrt{\frac{m^2}{n^2 - m'^2}} \pm \sqrt{\frac{m^2}{n^2 - m'^2} - 1} \tag{A2.45}
$$

**Behaviour by region:**

| Region | Character of $K_\pm$ | Stable direction |
|---|---|---|
| $m'^2 + m^2 < n^2$ (inside circle) | Complex conjugate, $|K| = 1$ | Both directions |
| $m'^2 + m^2 > n^2$ (outside circle) | Two negative real, $|K_1| < 1 < |K_2|$ | **From face inwards** (decreasing $|m'|$) |

> **Implementation rules of thumb:**
> - **(A2.33), (A2.34), (A2.35):** start at the pyramid face ($|m'| = n$ or $|m| = n$) and recurse **inwards** toward $m' = 0$ or $m = 0$.
> - **(A2.36):** start at face and recurse **vertically downward** (decreasing $n$).
> - The corner element $\Delta^n_{nn} = 1/2^n$ (A2.42) is the smallest entry; for large $n$ it **underflows**. Workaround: start the recursion **inside the pyramid but outside the $m'^2 + m^2 = n^2$ circle** with one zero and one arbitrary non-zero seed, then **renormalize** afterward against a known $\Delta^n_{1m}$ value.
> - **Do not normalize against $\Delta^n_{0m}$**: half of those entries are zero (parity rule).
> - **Do normalize against $\Delta^n_{1m}$**: none are zero for $m \ne 0$.

---

## 6. Low-Order Delta Tables (§A2.6)

Tables show the $(m \ge 0,\,m' \ge 0)$ block. Use the symmetries in §4 to extend to negative indices.

### $n = 0$

$$
\Delta^{0}_{00} = 1
$$

### $n = 1$

| $\Delta^{1}_{m'm}$ | $m=0$ | $m=1$ |
|---|---|---|
| $m' = 0$ | $0$ | $-\dfrac{\sqrt{2}}{2}$ |
| $m' = 1$ | $\dfrac{\sqrt{2}}{2}$ | $\dfrac{1}{2}$ |

### $n = 2$

| $\Delta^{2}_{m'm}$ | $m=0$ | $m=1$ | $m=2$ |
|---|---|---|---|
| $m' = 0$ | $-\dfrac{1}{2}$ | $0$ | $\dfrac{\sqrt{6}}{4}$ |
| $m' = 1$ | $0$ | $-\dfrac{1}{2}$ | $-\dfrac{1}{2}$ |
| $m' = 2$ | $\dfrac{\sqrt{6}}{4}$ | $\dfrac{1}{2}$ | $\dfrac{1}{4}$ |

### $n = 3$

| $\Delta^{3}_{m'm}$ | $m=0$ | $m=1$ | $m=2$ | $m=3$ |
|---|---|---|---|---|
| $m' = 0$ | $0$ | $\dfrac{\sqrt{3}}{4}$ | $0$ | $-\dfrac{\sqrt{5}}{4}$ |
| $m' = 1$ | $-\dfrac{\sqrt{3}}{4}$ | $-\dfrac{1}{8}$ | $\dfrac{\sqrt{10}}{8}$ | $\dfrac{\sqrt{15}}{8}$ |
| $m' = 2$ | $0$ | $-\dfrac{\sqrt{10}}{8}$ | $-\dfrac{1}{2}$ | $-\dfrac{\sqrt{6}}{8}$ |
| $m' = 3$ | $\dfrac{\sqrt{5}}{4}$ | $\dfrac{\sqrt{15}}{8}$ | $\dfrac{\sqrt{6}}{8}$ | $\dfrac{1}{8}$ |

### $n = 4$

| $\Delta^{4}_{m'm}$ | $m=0$ | $m=1$ | $m=2$ | $m=3$ | $m=4$ |
|---|---|---|---|---|---|
| $m' = 0$ | $\dfrac{3}{8}$ | $0$ | $-\dfrac{\sqrt{10}}{8}$ | $0$ | $\dfrac{\sqrt{70}}{16}$ |
| $m' = 1$ | $0$ | $\dfrac{3}{8}$ | $\dfrac{\sqrt{2}}{8}$ | $-\dfrac{\sqrt{7}}{8}$ | $-\dfrac{\sqrt{14}}{8}$ |
| $m' = 2$ | $-\dfrac{\sqrt{10}}{8}$ | $-\dfrac{\sqrt{2}}{8}$ | $\dfrac{1}{4}$ | $\dfrac{\sqrt{14}}{8}$ | $\dfrac{\sqrt{7}}{8}$ |
| $m' = 3$ | $0$ | $-\dfrac{\sqrt{7}}{8}$ | $-\dfrac{\sqrt{14}}{8}$ | $-\dfrac{3}{8}$ | $-\dfrac{\sqrt{2}}{8}$ |
| $m' = 4$ | $\dfrac{\sqrt{70}}{16}$ | $\dfrac{\sqrt{14}}{8}$ | $\dfrac{\sqrt{7}}{8}$ | $\dfrac{\sqrt{2}}{8}$ | $\dfrac{1}{16}$ |

### $n = 5$

| $\Delta^{5}_{m'm}$ | $m=0$ | $m=1$ | $m=2$ | $m=3$ | $m=4$ | $m=5$ |
|---|---|---|---|---|---|---|
| $m' = 0$ | $0$ | $-\dfrac{\sqrt{30}}{16}$ | $0$ | $\dfrac{\sqrt{35}}{16}$ | $0$ | $-\dfrac{3\sqrt{7}}{16}$ |
| $m' = 1$ | $\dfrac{\sqrt{30}}{16}$ | $\dfrac{1}{16}$ | $-\dfrac{\sqrt{7}}{8}$ | $-\dfrac{\sqrt{42}}{32}$ | $\dfrac{\sqrt{21}}{16}$ | $\dfrac{\sqrt{210}}{32}$ |
| $m' = 2$ | $0$ | $\dfrac{\sqrt{7}}{8}$ | $\dfrac{1}{4}$ | $-\dfrac{\sqrt{6}}{16}$ | $-\dfrac{\sqrt{3}}{4}$ | $-\dfrac{\sqrt{30}}{16}$ |
| $m' = 3$ | $-\dfrac{\sqrt{35}}{16}$ | $-\dfrac{\sqrt{42}}{32}$ | $\dfrac{\sqrt{6}}{16}$ | $\dfrac{13}{32}$ | $\dfrac{9\sqrt{2}}{32}$ | $\dfrac{3\sqrt{5}}{32}$ |
| $m' = 4$ | $0$ | $-\dfrac{\sqrt{21}}{16}$ | $-\dfrac{\sqrt{3}}{4}$ | $-\dfrac{9\sqrt{2}}{32}$ | $-\dfrac{1}{4}$ | $-\dfrac{\sqrt{10}}{32}$ |
| $m' = 5$ | $\dfrac{3\sqrt{7}}{16}$ | $\dfrac{\sqrt{210}}{32}$ | $\dfrac{\sqrt{30}}{16}$ | $\dfrac{3\sqrt{5}}{32}$ | $\dfrac{\sqrt{10}}{32}$ | $\dfrac{1}{32}$ |

---

## 7. Useful Identities Summary

| | |
|---|---|
| $\Delta^{n}_{nn}$ | $1/2^n$ — smallest, used to detect underflow |
| $d^n_{0,0}(\pi/2)$ | $\Delta^n_{00}$; zero for odd $n$ (parity rule) |
| $d^n_{\mu m}(-\theta)$ | $d^n_{m\mu}(\theta)$ (A2.7) |
| Parity in $\theta$ | $(\mu + m)$ |
| $d^n_{\mu m}(0)$ | $\delta_{\mu m}$ — rotation by 0 is identity |
| $d^n_{\mu m}(\pi)$ | $(-1)^{n+m}\delta_{\mu,-m}$ — rotation by $\pi$ flips $m$ |
| $d^n_{0m}$ vs $\bar P_n^{|m|}$ | $d^n_{0m}(\theta) = (-m/|m|)^m\,\sqrt{2/(2n+1)}\,\bar P_n^{|m|}(\cos\theta)$ |
| $\sum_\mu (d^n_{\mu m})^2$ | $1$ (normalization, A2.6) |

---

## 8. Implementation Checklist for Python Code

### Convention
1. **Euler-angle tuple order.** The printed triple is $(\chi_o, \theta_o, \phi_o)$, but the rotations apply $\phi_o$ first (about $z$), then $\theta_o$ (about $y_1$), then $\chi_o$ (about $z_2$). When using `scipy.spatial.transform.Rotation.from_euler("ZYZ", [phi_o, theta_o, chi_o])`, **reverse the printed tuple** to match scipy's listed order.
2. **Right-hand rule.** Positive $\phi$-rotation about $\hat{z}$: $\hat{x} \to \hat{y}$. Verify by rotating $\hat{x}$ by $\pi/2$ about $\hat{z}$ and asserting the result is $\hat{y}$.

### $d^n_{\mu m}(\theta)$
3. **Edmonds formula (A2.3).** Implement and verify against the Jacobi form (A2.5). Spot-check at $\theta = \pi/2$ against the §6 tables for $n \le 5$.
4. **Realness.** $d^n_{\mu m}(\theta) \in \mathbb{R}$ for real $\theta$. Imaginary part must be zero to machine precision.
5. **Parity.** $d^n_{\mu m}(-\theta) = d^n_{m\mu}(\theta)$ — Eq. (A2.7). Spot-check with random $(\mu, m, n, \theta)$.
6. **$\theta = 0$ identity matrix.** $d^n_{\mu m}(0) = \delta_{\mu m}$ — Eq. (A2.15). Should be exact, not approximate.
7. **$\theta = \pi$ off-diagonal.** $d^n_{\mu m}(\pi) = (-1)^{n+m}\delta_{\mu,-m}$ — Eq. (A2.16).
8. **Normalization.** $\sum_\mu (d^n_{\mu m}(\theta))^2 = 1$ — Eq. (A2.6). Test for $n = 1\dots 20$, random $\theta$.
9. **Legendre link.** $d^n_{0m}(\theta) = (-m/|m|)^m\,\sqrt{2/(2n+1)}\,\bar P_n^{|m|}(\cos\theta)$ — Eq. (A2.17). Numerically compare against your $\bar P$ implementation (Appendix A1) at, e.g., $\theta = \pi/3$.
10. **Chapter-2 dipole regression.** $d^1_{-1,0}(-\pi/2) = \sqrt{2}/2$, $d^1_{0,0}(-\pi/2) = 0$, $d^1_{1,0}(-\pi/2) = -\sqrt{2}/2$ (Eqs. 2.120–2.122).
11. **Orthogonality integral.** Eq. (A2.10): $\int_0^\pi d^n_{\mu m} d^{n'}_{\mu m} \sin\theta\,d\theta = 2/(2n+1)\cdot\delta_{nn'}$ for fixed $\mu, m$. Use Gauss–Legendre quadrature with order $\ge n + n' + 1$.

### Deltas $\Delta^n_{m'm}$
12. **Tables.** For $n = 0\dots 5$, build all $\Delta^n_{m'm}$ and compare to §6 tables exactly (rationals + surds).
13. **Symmetries.** Eqs. (A2.26)–(A2.32). All seven must hold for arbitrary $(n, m', m)$.
14. **Parity in $\Delta^n_{m'0}$, $\Delta^n_{0m}$.** Eqs. (A2.38)/(A2.39): zero when $(n + m')$ (or $(n+m)$) is odd. Inspect the n=3 and n=5 columns/rows in the tables.
15. **Corner identity.** $\Delta^n_{nn} = 1/2^n$ — Eq. (A2.42). Trivial but useful smoke test.
16. **Edge identity.** $\Delta^n_{nm} = \tfrac{1}{2^n}\sqrt{\binom{2n}{n-m}}$ — Eq. (A2.41). Test for $n = 5$, all $m$.
17. **Recurrence consistency.** Implement (A2.33) and verify against table values for $n = 5$. **Run the recurrence inwards from the pyramid face** ($m' = \pm n$) to avoid instability outside the circle $m^2 + m'^2 = n^2$.
18. **Vertical recurrence.** Eq. (A2.36) — run from large $n$ downward (table-comparison test for $n = 0 \to 5$).
19. **Underflow handling.** For large $n$ (say $n \ge 60$), $\Delta^n_{nn} = 2^{-n}$ underflows in IEEE double. Implement seeded recurrence with renormalization against a non-zero $\Delta^n_{1m}$.

### Wave-function rotation
20. **Single-mode rotation (A2.1).** For a pure $\vec{F}^{(c)}_{smn}$ at $(r,\theta,\phi)$ in the unprimed frame, the expansion (A2.1) should reproduce it when the primed frame is identical (all Euler angles = 0). Verify: $D^n_{\mu m}(0,0,0) = \delta_{\mu m}$.
21. **Decomposition (A2.2).** $D^n_{\mu m}(\chi_o, \theta_o, \phi_o) = e^{im\phi_o}\,d^n_{\mu m}(\theta_o)\,e^{i\mu\chi_o}$. Test all three pure-phase cases ($\theta_o = 0$ with random $\phi_o, \chi_o$ — should give unimodular complex numbers).
22. **Chapter-2 dipole rotation.** Use the example from §2 — rotating $\vec{F}^{(3)}_{201}$ by $(0, -\pi/2, 0)$ should reproduce Eq. (2.123) of Chapter 2: equal-and-opposite $\pm\sqrt{2}/2$ on $\vec{F}^{(3)}_{2,\mp 1,1}$ and zero elsewhere.
23. **Rotation invariants.** Rotating a single-mode coefficient set yields only modes of the **same** $n$ and **same** $s$. Diagonal blocks in $n$ and $s$ — non-zero in other blocks would be a bug.
24. **Composition.** Two successive rotations $R_1(\chi_1,\theta_1,\phi_1)$ then $R_2(\chi_2,\theta_2,\phi_2)$ should yield the same result as a single rotation by the composed Euler angles. Generate random rotations, compose with scipy's ZYZ extrinsic/intrinsic convention you chose, and compare to applying $D^n$ twice.
25. **Inverse rotation.** The inverse of Euler angles $(\chi_o, \theta_o, \phi_o)$ in this convention is $(-\phi_o, -\theta_o, -\chi_o)$ (ZYZ inverse). Applying $D$ for $(\chi_o,\theta_o,\phi_o)$ then $D$ for $(-\phi_o,-\theta_o,-\chi_o)$ should give the identity on mode coefficients.

### Fourier expansion (Wacker algorithm support)
26. **Series (A2.11)/(A2.12).** Verify the two equivalent forms reconstruct $d^n_{\mu m}(\theta)$ for random $\theta$. Note the $i^{\mu-m}$ vs $i^{m-\mu}$ phases and the sign of the exponential.
27. **Discrete computation.** Once all $\Delta^n_{m'm}$ for a fixed $n$ are tabulated, $d^n_{\mu m}(\theta_i)$ at a list of equispaced $\theta_i$ values is **one FFT per $(\mu, m)$ pair** — see §A2.5.2 final paragraph and Appendix 4 of the book.

---

*End of Appendix A2 reference.*
