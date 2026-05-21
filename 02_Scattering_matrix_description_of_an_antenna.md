# Chapter 2 — Scattering Matrix Description of an Antenna

**Source:** Hansen, *Spherical Near-Field Antenna Measurements*, Chapter 2 (pp. 8–60).
**Cross-reference:** Symbols defined in [A5_List_of_principal_symbols_and_uses.md](A5_List_of_principal_symbols_and_uses.md).

> **Purpose of this document:** structured reference for validating Python implementations of spherical-wave / scattering-matrix antenna code. Every equation is tagged with its book number `(2.N)` so that test cases and code comments can cite them directly.

---

## 0. Global Conventions (must be honoured by any implementation)

| Convention | Value | Comment |
|---|---|---|
| Time factor | $\exp(-i\omega t)$ | Suppressed throughout. **Sign flip = wrong sign of $i$ in every radial-function expression.** |
| Imaginary unit | $i$ | Not $j$. Outgoing wave $\propto e^{+ikr}$ under this convention. |
| Spherical coords | $r, \theta, \phi$ | $0 \le r < \infty$, $0 \le \theta \le \pi$, $0 \le \phi < 2\pi$ |
| Wavenumber | $k = \omega\sqrt{\mu\varepsilon} = 2\pi/\lambda$ | Real for loss-free; complex allowed (except power & directivity formulas). |
| Specific admittance | $\eta = \sqrt{\varepsilon/\mu}$ | Used in field-coefficient prefactors. |
| Specific impedance | $\zeta = \sqrt{\mu/\varepsilon}$ | Equals $1/\eta$. |
| Index $s$ | $s \in \{1, 2\}$ | $s=1$ ↔ TE in coefficient $Q$; $s=2$ ↔ TM in coefficient $Q$. **In the wave-function $\vec{F}^{(c)}_{smn}$ the meaning is reversed**: $\vec{F}^{(c)}_{1mn}$ is transverse (no radial component); $\vec{F}^{(c)}_{2mn}$ carries the radial component. |
| Index $n$ | $n = 1, 2, 3, \dots$ | No $n=0$. |
| Index $m$ | $m = -n, -n+1, \dots, n-1, n$ | $\|m\| \le n$. |
| Index $c$ | $c \in \{1, 2, 3, 4\}$ | $c=1$: $j_n$ (finite at origin). $c=2$: $n_n$ (singular at origin). $c=3$: $h_n^{(1)}$ (outgoing). $c=4$: $h_n^{(2)}$ (incoming). |
| Power normalization | A single $c=3$ mode with $\|Q^{(3)}_{smn}\| = 1$ radiates $\tfrac{1}{2}$ W. | The factor $k/\sqrt{\eta}$ in $\vec{E}$ is what makes this hold. |
| Wave-coefficient units | $[Q^{(c)}_{smn}] = \text{W}^{1/2}$ | The functions $\vec{F}^{(c)}_{smn}$ are dimensionless. |
| Waveguide-port units | $[v] = [w] = \text{W}^{1/2}$ | Incident power at local port is $\tfrac{1}{2}\|v\|^2$. |
| Convention $(-m/\|m\|)^m$ | $= 1$ when $m = 0$ (Eq. 2.19). | Ensures the Edmonds phase. |

---

## 1. Maxwell's Equations and the Vector Wave Equation (§2.2.1)

In a linear, isotropic, homogeneous medium with assumed sources $\vec{J}$ (electric), $\vec{M}$ (magnetic):

$$
\nabla \times \vec{H} = -i\omega\varepsilon \vec{E} + \vec{J} \tag{2.1}
$$

$$
\nabla \times \vec{E} = i\omega\mu \vec{H} - \vec{M} \tag{2.2}
$$

In a source-free region both $\vec{E}$ and $\vec{H}$ satisfy:

$$
\nabla \times (\nabla \times \vec{C}) - k^2 \vec{C} = 0 \tag{2.3}
$$

### Generating function and Hansen vector functions

The scalar Helmholtz equation $(\nabla^2 + k^2) f = 0$ (Eq. 2.4) generates:

$$
\vec{m} = \nabla f \times \vec{r} \tag{2.5}
$$

$$
\vec{n} = k^{-1} \nabla \times \vec{m} \tag{2.6}
$$

Reciprocal curl relations:

$$
\vec{m} = k^{-2} \nabla \times (\nabla \times \vec{m}) \tag{2.7}
$$

$$
\vec{m} = k^{-1} \nabla \times \vec{n} \tag{2.8}
$$

> The third (irrotational) Hansen function $\vec{l} = \nabla f$ is not needed here.

### Stratton generating function

$$
f^{(c)}_{\sigma mn}(r,\theta,\phi) = z^{(c)}_n(kr)\, P_n^m(\cos\theta)\, \begin{cases}\cos m\phi\\\sin m\phi\end{cases} \tag{2.9}
$$

(With $\sigma \in \{e, o\}$ for even/odd trig.) $P_n^m$ is the unnormalized associated Legendre function.

### Radial functions $z_n^{(c)}(kr)$

| $c$ | $z_n^{(c)}$ | Wave type |
|---|---|---|
| 1 | $j_n(kr)$ — spherical Bessel | Standing, finite at origin (2.10a) |
| 2 | $n_n(kr)$ — spherical Neumann | Standing, singular at origin (2.10b) |
| 3 | $h_n^{(1)}(kr) = j_n(kr) + i\, n_n(kr)$ | **Outward travelling** (2.10c) |
| 4 | $h_n^{(2)}(kr) = j_n(kr) - i\, n_n(kr)$ | **Inward travelling** (2.10d) |

**Bessel identity (used for matching at origin):**

$$
j_n(kr) = \tfrac{1}{2}\bigl(h_n^{(1)}(kr) + h_n^{(2)}(kr)\bigr) \tag{2.30}
$$

### Large-argument (far-field) asymptotics for $kr \to \infty$, $kr \gg n$

$$
z_n^{(3)}(kr) \to (-i)^{n+1} \frac{e^{ikr}}{kr} \tag{2.13}
$$

$$
z_n^{(4)}(kr) \to i^{n+1} \frac{e^{-ikr}}{kr} \tag{2.14}
$$

$$
\frac{1}{kr}\frac{d}{d(kr)}\bigl\{kr\, z_n^{(3)}(kr)\bigr\} \to (-i)^{n}\frac{e^{ikr}}{kr} \tag{2.15}
$$

$$
\frac{1}{kr}\frac{d}{d(kr)}\bigl\{kr\, z_n^{(4)}(kr)\bigr\} \to i^{n}\frac{e^{-ikr}}{kr} \tag{2.16}
$$

---

## 2. Power-Normalized Spherical Wave Functions (§2.2.2)

### Modified scalar generating function (Jensen [13])

$$
F^{(c)}_{mn}(r,\theta,\phi) = \frac{1}{\sqrt{2\pi}}\,\frac{1}{\sqrt{n(n+1)}}\Bigl(-\frac{m}{|m|}\Bigr)^{m} z_n^{(c)}(kr)\, \bar{P}_n^{|m|}(\cos\theta)\, e^{im\phi} \tag{2.18}
$$

with $\bar{P}_n^{|m|}$ the **normalized** associated Legendre function (Belousov [14]; see Appendix A1).

### Vector wave functions

**$s=1$ (transverse — no radial component):**

$$
\vec{F}^{(c)}_{1mn}(r,\theta,\phi) = \nabla F^{(c)}_{mn} \times \vec{r}
$$

$$= \frac{1}{\sqrt{2\pi}}\frac{1}{\sqrt{n(n+1)}}\Bigl(-\frac{m}{|m|}\Bigr)^{m}\left\{
\,z_n^{(c)}(kr)\,\frac{im\,\bar{P}_n^{|m|}(\cos\theta)}{\sin\theta}\,e^{im\phi}\,\hat{\theta}
-z_n^{(c)}(kr)\,\frac{d\bar{P}_n^{|m|}(\cos\theta)}{d\theta}\,e^{im\phi}\,\hat{\phi}
\right\} \tag{2.20}$$

**$s=2$ (carries radial component):**

$$
\vec{F}^{(c)}_{2mn}(r,\theta,\phi) = k^{-1}\nabla \times \vec{F}^{(c)}_{1mn}
$$

$$= \frac{1}{\sqrt{2\pi}}\frac{1}{\sqrt{n(n+1)}}\Bigl(-\frac{m}{|m|}\Bigr)^{m}\Biggl\{
\frac{n(n+1)}{kr}\,z_n^{(c)}(kr)\,\bar{P}_n^{|m|}(\cos\theta)\,e^{im\phi}\,\hat{r}
$$
$$
+\,\frac{1}{kr}\frac{d}{d(kr)}\bigl\{kr z_n^{(c)}(kr)\bigr\}\,\frac{d\bar{P}_n^{|m|}(\cos\theta)}{d\theta}\,e^{im\phi}\,\hat{\theta}
$$
$$
+\,\frac{1}{kr}\frac{d}{d(kr)}\bigl\{kr z_n^{(c)}(kr)\bigr\}\,\frac{im\,\bar{P}_n^{|m|}(\cos\theta)}{\sin\theta}\,e^{im\phi}\,\hat{\phi}
\Biggr\} \tag{2.21}$$

> Both $\vec{F}^{(c)}_{1mn}$ and $\vec{F}^{(c)}_{2mn}$ are **dimensionless**.

### Field expansion (in a source-free region such as region 2 of Fig. 2.2)

$$
\vec{E}(r,\theta,\phi) = \frac{k}{\sqrt{\eta}}\sum_{c\,s\,m\,n} Q^{(c)}_{smn}\,\vec{F}^{(c)}_{smn}(r,\theta,\phi) \tag{2.22}
$$

$$
\vec{H}(r,\theta,\phi) = (i\omega\mu)^{-1}\nabla\times\vec{E} = -ik\sqrt{\eta}\sum_{c\,s\,m\,n} Q^{(c)}_{smn}\,\vec{F}^{(c)}_{3-s,m,n}(r,\theta,\phi) \tag{2.23}
$$

**Note the index swap $s \to 3-s$ in the magnetic field.** A correct implementation must mirror this.

### Summation conventions

$$
\sum_{smn} \;=\; \sum_{s=1}^{2}\sum_{n=1}^{\infty}\sum_{m=-n}^{n} \tag{2.25}
$$

$$
\sum_{csmn} \;=\; \sum_{c=3}^{4}\sum_{s=1}^{2}\sum_{n=1}^{\infty}\sum_{m=-n}^{n} \tag{2.26}
$$

### Single-index transformation (when $M = N$)

$$
j = 2\{n(n+1) + m - 1\} + s \tag{2.27}
$$

$$
\sum_{s=1}^{2}\sum_{n=1}^{N}\sum_{m=-n}^{n} = \sum_{j=1}^{J} \tag{2.28}
$$

$$
J = 2N(N+2) \tag{2.29}
$$

Substitutions: $Q^{(c)}_{smn} = Q^{(c)}_j$, $\vec{F}^{(c)}_{smn} = \vec{F}^{(c)}_j$.

> **Validation table (n=1 modes):**
>
> | $s$ | $m$ | $n$ | $j$ |
> |---|---|---|---|
> | 1 | -1 | 1 | 1 |
> | 2 | -1 | 1 | 2 |
> | 1 | 0  | 1 | 3 |
> | 2 | 0  | 1 | 4 |
> | 1 | +1 | 1 | 5 |
> | 2 | +1 | 1 | 6 |
> | 1 | -2 | 2 | 7 |
> | 2 | -2 | 2 | 8 |

### Power radiated (outgoing-wave field)

$$
P = \frac{1}{2}\sum_{smn}|Q^{(3)}_{smn}|^2 \quad\text{[watts]} \tag{2.24, 2.55}
$$

### Mode-type identification

- For a TM wave: $\vec{H} \propto \vec{F}^{(c)}_{1mn}$ (no radial), $\vec{E} \propto \vec{F}^{(c)}_{2mn}$. Coefficient $Q^{(c)}_{2mn}$ ($s=2$).
- For a TE wave: $\vec{E} \propto \vec{F}^{(c)}_{1mn}$, $\vec{H} \propto \vec{F}^{(c)}_{2mn}$. Coefficient $Q^{(c)}_{1mn}$ ($s=1$).
- Sum of wave-function $s$ indices for $\vec{E}$ and $\vec{H}$ of the same mode = 3.

### Truncation rule

Convergence is obtained for:

$$
N = \lfloor kr_0 \rfloor + n_1 \tag{2.31}
$$

with $r_0$ = minimum-sphere radius, $n_1$ depending on accuracy/geometry. **Empirical default: $n_1 = 10$** when the field point is more than a few wavelengths from the minimum sphere and 4 correct digits suffice. $n_1$ grows roughly as $(kr_0)^{1/3}$ for fixed accuracy.

For an interior (point closer than the source-free inner radius $r_i$):

$$
N = \lfloor kr \rfloor + n_1 \tag{2.32}
$$

### Region definitions for a multi-mode antenna ($N$ = truncation)

| Region | Range |
|---|---|
| Evanescent | $r_0 \lesssim r \lesssim N/k$ |
| Fresnel (near-field) | $N/k \lesssim r \lesssim 4N^2/(\pi k)$ |
| Fraunhofer (far-field) | $4N^2/(\pi k) \lesssim r < \infty$ |

### Rayleigh distance

$$
R = \frac{2D^2}{\lambda} = \frac{2}{\lambda}(2 r_0)^2 \tag{2.33}
$$

$$
\approx \frac{2}{\lambda}\Bigl(2\frac{N}{k}\Bigr)^2 = \frac{4}{\pi}\frac{N^2}{k} \tag{2.34}
$$

### On-axis behaviour (singular at poles, used in dipole derivations)

$$
[\vec{F}^{(c)}_{smn}(r,0,\phi)]_r = [\vec{F}^{(c)}_{smn}(r,\pi,\phi)]_r = 0 \quad\text{for } m \ne 0 \tag{2.35}
$$

$$
[\vec{F}^{(c)}_{smn}(r,0,\phi)]_{\theta,\phi} = [\vec{F}^{(c)}_{smn}(r,\pi,\phi)]_{\theta,\phi} = 0 \quad\text{for } m \ne \pm 1 \tag{2.36}
$$

---

## 3. Power Flow (§2.2.4)

### Complex Poynting vector

$$
\vec{S} = \tfrac{1}{2}\,\vec{E}\times\vec{H}^* \tag{2.37}
$$

### Power radiated through a sphere of radius $r$

$$
P = \int_0^{2\pi}\!\!\int_0^{\pi} \mathrm{Re}(\hat{r}\cdot\vec{S})\, r^2\sin\theta\, d\theta\, d\phi \tag{2.38}
$$

### Reactive-energy relation

$$
2\omega(W_e - W_m) = \int_0^{2\pi}\!\!\int_0^{\pi} \mathrm{Im}(\hat{r}\cdot\vec{S})\, r^2\sin\theta\, d\theta\, d\phi \tag{2.39}
$$

For a TE mode ($s=1$): $W_m > W_e$. For a TM mode ($s=2$): $W_e > W_m$.

### Far-field plane-wave relation

$$
\vec{H} = \eta\,\hat{r}\times\vec{E} \tag{2.40}
$$

$$
\vec{S} = \tfrac{1}{2}\eta|\vec{E}|^2\,\hat{r} \quad\text{[W/m²]} \tag{2.41}
$$

$$
P_1(\theta,\phi) = \tfrac{1}{2}\eta|\vec{E}|^2 r^2 \quad\text{[W/sr]} \tag{2.42}
$$

### Complex-conjugate identity for $\vec{F}^{(c)}_{smn}$

$$
\vec{F}^{(3)*}_{smn}(r,\theta,\phi) = (-1)^m\,\vec{F}^{(4)}_{s,-m,n}(r,\theta,\phi) \tag{2.45}
$$

### Orthogonality integral (used for power and reciprocity)

$$
\int_0^{2\pi}\!\!\int_0^{\pi} \bigl\{\vec{F}^{(c)}_{smn}\times\vec{F}^{(\gamma)}_{\sigma\mu\nu}\bigr\}\cdot\hat{r}\,\sin\theta\,d\theta\,d\phi
$$

$$
= \delta_{3-s,\sigma}\,\delta_{m,-\mu}\,\delta_{n\nu}\,(-1)^{3-s}\,(-1)^m\,R^{(c)}_{sn}(kr)\,R^{(\gamma)}_{3-s,n}(kr) \tag{2.46}
$$

with the radial-function abbreviation

$$
R^{(c)}_{sn}(kr) = \begin{cases} z_n^{(c)}(kr), & s = 1 \\ \dfrac{1}{kr}\dfrac{d}{d(kr)}\{kr\,z_n^{(c)}(kr)\}, & s = 2 \end{cases} \tag{2.47}
$$

### Kronecker delta and Wronskian

$$
\delta_{ij} = \begin{cases}0 & i \ne j \\ 1 & i = j\end{cases} \tag{2.48}
$$

$$
R^{(3)}_{sn}(kr) = R^{(1)}_{sn}(kr) + i\,R^{(2)}_{sn}(kr) \tag{2.50}
$$

$$
R^{(4)}_{sn}(kr) = R^{(1)}_{sn}(kr) - i\,R^{(2)}_{sn}(kr) \tag{2.51}
$$

$$
R^{(1)}_{1n}(kr)R^{(2)}_{2n}(kr) - R^{(1)}_{2n}(kr)R^{(2)}_{1n}(kr) = (kr)^{-2} \tag{2.52}
$$

### Complex-power-flux result (outgoing $c=3$ field)

$$\tfrac{1}{2}\int(\vec{E}\times\vec{H}^*)\cdot\hat{r}\,r^2\sin\theta\,d\theta\,d\phi
= \sum_{smn}\Bigl\{\tfrac{1}{2} + \tfrac{1}{2}i(-1)^{3-s}(kr)^2 V_n(kr)\Bigr\}\bigl|Q^{(3)}_{smn}\bigr|^2 \tag{2.53}$$

with cross-product

$$
V_n(kr) = R^{(1)}_{1n}(kr)R^{(2)}_{1n}(kr) + R^{(1)}_{2n}(kr)R^{(2)}_{2n}(kr) = \Bigl(\frac{1}{kr} + \frac{1}{2}\frac{d}{d(kr)}\Bigr)|h_n^{(1)}(kr)|^2 \tag{2.54}
$$

$V_n(kr)$ is always negative (Abramowitz & Stegun 10.1.27).

**Real power:** independent of $r$, as in (2.55) above.

---

## 4. The Antenna Scattering Matrix (§2.3.1)

### Field outside the minimum sphere

$$
\vec{E}(r,\theta,\phi) = \frac{k}{\sqrt{\eta}}\sum_{j=1}^{J}\bigl\{a_j\,\vec{F}^{(4)}_j(r,\theta,\phi) + b_j\,\vec{F}^{(3)}_j(r,\theta,\phi)\bigr\},\quad r > r_0 \tag{2.56}
$$

Index relation: $Q^{(4)}_{smn} = a_j = a_{smn}$ (incoming), $Q^{(3)}_{smn} = b_j = b_{smn}$ (outgoing), with $j$ from (2.27).

### Antenna scattering equation (order $J+1$)

$$
\begin{bmatrix}\Gamma & \mathbf{R}\\ \mathbf{T} & \mathbf{S}\end{bmatrix}\begin{bmatrix}v\\ \mathbf{a}\end{bmatrix} = \begin{bmatrix}w\\ \mathbf{b}\end{bmatrix} \tag{2.57}
$$

Expanded:

$$
\Gamma v + \sum_{j=1}^{J} R_j a_j = w \tag{2.58}
$$

$$
T_i v + \sum_{j=1}^{J} S_{ij} a_j = b_i,\quad i = 1, 2, \dots, J \tag{2.59}
$$

Condensed: $\hat{\mathbf{S}}\hat{\mathbf{a}} = \hat{\mathbf{b}}$ (2.60), where $\hat{\mathbf{S}}$ is the **total** scattering matrix.

| Symbol | Shape | Meaning |
|---|---|---|
| $\Gamma$ | scalar | Antenna reflection coefficient |
| $\mathbf{R}$ | $1 \times J$ row | Receiving coefficients $R_j$ |
| $\mathbf{T}$ | $J \times 1$ column | Transmitting coefficients $T_i$ |
| $\mathbf{S}$ | $J \times J$ | Scattering coefficients $S_{ij}$ |

All elements of $\hat{\mathbf{S}}$ are dimensionless.

### Lossless antenna — unitarity

$$\hat{\mathbf{S}}^{+}\hat{\mathbf{S}} = \hat{\mathbf{I}}$$ (unit matrix of order $J+1$).

In particular, for the first column:

$$
|\Gamma|^2 + |\mathbf{T}|^2 = 1 \tag{2.61}
$$

### Lossy antenna

$$
\tfrac{1}{2}\sum_{i=1}^{J}|b_i|^2 = \tfrac{1}{2}|v|^2 - \bigl(\tfrac{1}{2}|w|^2 + P_{\text{loss}}\bigr) \tag{2.62}
$$

$$
|\Gamma|^2 + |\mathbf{T}|^2 = 1 - \frac{P_{\text{loss}}}{P_{\text{inc}}} \tag{2.63}
$$

In the lossy case, every row and column of $\hat{\mathbf{S}}$ has norm $\le 1$.

### Empty space

$\hat{\mathbf{S}}$ undefined, but $\mathbf{S} = \mathbf{I}$ (unit matrix of infinite order). Outgoing equals incoming: $\mathbf{b} = \mathbf{a}$.

### Generator coupling — radiated waves

Generator: $v = v_g + \Gamma_g w$ (2.64). With $\mathbf{a} = 0$:

$$
\Gamma v = w \tag{2.65}
$$

$$
\mathbf{T} v = \mathbf{b} \tag{2.66}
$$

$$
\boxed{\;\mathbf{b} = \frac{v_g}{1 - \Gamma_g \Gamma}\,\mathbf{T}\;} \tag{2.67}
$$

$$
\vec{E}(r,\theta,\phi) = \frac{k}{\sqrt{\eta}}\sum_{i=1}^{J} b_i\,\vec{F}^{(3)}_i(r,\theta,\phi) \tag{2.68}
$$

### Load coupling — received wave & scattering

Load: $v = \Gamma_l w$ (2.69). With $\Gamma v + \mathbf{R}\mathbf{a} = w$ (2.70):

$$
\boxed{\;w = \frac{1}{1-\Gamma_l\Gamma}\,\mathbf{R}\mathbf{a}\;} \tag{2.71}
$$

Power accepted by the load:

$$
P = \tfrac{1}{2}(1-|\Gamma_l|^2)|w|^2 = \tfrac{1}{2}(1-|\Gamma_l|^2)\,\frac{|\mathbf{R}\mathbf{a}|^2}{|1-\Gamma_l\Gamma|^2} \tag{2.72}
$$

**Matched load** ($\Gamma_l = 0$): $P = P' = \tfrac{1}{2}|\mathbf{R}\mathbf{a}|^2$ (2.73).

**Conjugate-matched load** ($\Gamma_l = \Gamma^*$): maximum, "available" power:

$$
P_a = \frac{1}{2}\,\frac{|\mathbf{R}\mathbf{a}|^2}{1-|\Gamma|^2} \tag{2.74}
$$

### Scattered field

From $\mathbf{T} v + \mathbf{S}\mathbf{a} = \mathbf{b}$ (2.75) and (2.69), (2.71):

$$
\mathbf{b} = \{\mathbf{T}\Gamma_l(1-\Gamma\Gamma_l)^{-1}\mathbf{R} + \mathbf{S}\}\mathbf{a} \tag{2.76}
$$

**Subtract empty-space scattering** to get the field actually scattered by the antenna:

$$
\mathbf{b}' = \{\mathbf{T}\Gamma_l(1-\Gamma\Gamma_l)^{-1}\mathbf{R} + (\mathbf{S} - \mathbf{I})\}\mathbf{a} \tag{2.77}
$$

Matched load ($\Gamma_l = 0$): $\mathbf{b}' = (\mathbf{S} - \mathbf{I})\mathbf{a}$ (2.78).

$$
\vec{E}'(r,\theta,\phi) = \frac{k}{\sqrt{\eta}}\sum_{i=1}^{J} b'_i\,\vec{F}^{(3)}_i(r,\theta,\phi) \tag{2.79}
$$

---

## 5. Reciprocity (§2.3.2)

### Generalized Lorentz form (unprimed and primed = adjoint device)

$$
vw' - v'w = \sum_{smn}(-1)^m\bigl(b_{smn}\,a'_{s,-m,n} - a_{smn}\,b'_{s,-m,n}\bigr) \tag{2.93}
$$

### Adjoint-antenna relations (with $\mathbf{R}'$, $\mathbf{T}$, etc.)

$$
R'_{smn} = (-1)^m\,T_{s,-m,n} \tag{2.104}
$$

$$
S'^{\sigma\mu\nu}_{smn} = (-1)^{m+\mu}\,S^{s,-m,n}_{\sigma,-\mu,\nu} \tag{2.105}
$$

$$
\Gamma' = \Gamma \tag{2.106}
$$

### Reciprocal antenna (special case — antenna is its own adjoint)

$$
\boxed{\;R_{smn} = (-1)^m\,T_{s,-m,n}\;} \tag{2.107}
$$

$$
S^{\sigma\mu\nu}_{smn} = (-1)^{m+\mu}\,S^{s,-m,n}_{\sigma,-\mu,\nu} \tag{2.108}
$$

---

## 6. Fields of Electric and Magnetic Dipoles (§2.3.3)

### z-directed electric dipole at origin

Dipole moment $d_e = I\ell$.

$$
\vec{E}^{z}_e = -\frac{\zeta k^2}{2\pi}d_e\,\frac{h_1^{(1)}(kr)}{kr}\cos\theta\,\hat{r} + \frac{\zeta k^2}{4\pi}d_e\,\frac{1}{kr}\frac{d}{d(kr)}\{kr\,h_1^{(1)}(kr)\}\sin\theta\,\hat{\theta} \tag{2.109}
$$

$$
\vec{H}^{z}_e = \frac{ik^2}{4\pi}d_e\,h_1^{(1)}(kr)\sin\theta\,\hat{\phi} \tag{2.110}
$$

with

$$
h_1^{(1)}(kr) = -\frac{e^{ikr}}{kr}\Bigl(1 + \frac{i}{kr}\Bigr) \tag{2.111}
$$

$$
\frac{1}{kr}\frac{d}{d(kr)}\{kr\,h_1^{(1)}(kr)\} = \frac{e^{ikr}}{kr}\Bigl\{-i + \frac{1}{kr} + \frac{i}{(kr)^2}\Bigr\} \tag{2.112}
$$

### Spherical-wave form

$$
\vec{E}^{z}_e = \frac{k}{\sqrt{\eta}}\,Q_{201}\,\vec{F}^{(3)}_{201}(r,\theta,\phi) \tag{2.115}
$$

$$
\vec{H}^{z}_e = -ik\sqrt{\eta}\,Q_{201}\,\vec{F}^{(3)}_{101}(r,\theta,\phi) \tag{2.116}
$$

$$
\boxed{\;Q_{201} = -\frac{1}{\sqrt{6\pi}}\,\frac{k}{\sqrt{\eta}}\,d_e\;} \tag{2.117}
$$

### z-directed magnetic dipole

Moment $d_m = I_m \ell = -i\omega\mu S I'$ (2.130, 2.131).

$$
\vec{E}^{z}_m = \frac{k}{\sqrt{\eta}}\,Q_{101}\,\vec{F}^{(3)}_{101}(r,\theta,\phi) \tag{2.132}
$$

$$
\vec{H}^{z}_m = -ik\sqrt{\eta}\,Q_{101}\,\vec{F}^{(3)}_{201}(r,\theta,\phi) \tag{2.133}
$$

$$
\boxed{\;Q_{101} = -\frac{i}{\sqrt{6\pi}}\,k\sqrt{\eta}\,d_m\;} \tag{2.138}
$$

### Duality

If $d_m = -\zeta d_e$ (2.139):

$$
\vec{E}_m = \zeta\vec{H}_e,\quad \vec{H}_m = -\eta\vec{E}_e \tag{2.140, 2.141}
$$

$$
Q_{101} = -i\,Q_{201} \tag{2.142}
$$

### Rotation result for x-directed dipole (Euler angles $(0,-\pi/2,0)$)

Rotation coefficients used:

$$
d^{1}_{-1,0}(-\pi/2) = \tfrac{\sqrt{2}}{2},\quad d^{1}_{0,0}(-\pi/2) = 0,\quad d^{1}_{1,0}(-\pi/2) = -\tfrac{\sqrt{2}}{2} \tag{2.120-2.122}
$$

Yielding:

$$
\vec{F}^{(3)}_{201}(r,\theta,\phi) = \tfrac{\sqrt{2}}{2}\vec{F}^{(3)}_{2,-1,1}(r',\theta',\phi') - \tfrac{\sqrt{2}}{2}\vec{F}^{(3)}_{211}(r',\theta',\phi') \tag{2.123}
$$

### Fields of x- and y-directed dipoles

**x-directed electric dipole:**

$$
\vec{E}^{x}_e = \frac{k}{\sqrt{\eta}}Q_{201}\,\tfrac{\sqrt{2}}{2}\{\vec{F}^{(3)}_{2,-1,1} - \vec{F}^{(3)}_{211}\} \tag{2.124}
$$

$$
\vec{H}^{x}_e = -ik\sqrt{\eta}\,Q_{201}\,\tfrac{\sqrt{2}}{2}\{\vec{F}^{(3)}_{1,-1,1} - \vec{F}^{(3)}_{111}\} \tag{2.125}
$$

**y-directed electric dipole:**

$$
\vec{E}^{y}_e = \frac{k}{\sqrt{\eta}}Q_{201}\,\tfrac{i\sqrt{2}}{2}\{\vec{F}^{(3)}_{2,-1,1} + \vec{F}^{(3)}_{211}\} \tag{2.128}
$$

$$
\vec{H}^{y}_e = -ik\sqrt{\eta}\,Q_{201}\,\tfrac{i\sqrt{2}}{2}\{\vec{F}^{(3)}_{1,-1,1} + \vec{F}^{(3)}_{111}\} \tag{2.129}
$$

**x-directed magnetic dipole:**

$$
\vec{E}^{x}_m = \frac{k}{\sqrt{\eta}}Q_{101}\,\tfrac{\sqrt{2}}{2}\{\vec{F}^{(3)}_{1,-1,1} - \vec{F}^{(3)}_{111}\} \tag{2.134}
$$

$$
\vec{H}^{x}_m = -ik\sqrt{\eta}\,Q_{101}\,\tfrac{\sqrt{2}}{2}\{\vec{F}^{(3)}_{2,-1,1} - \vec{F}^{(3)}_{211}\} \tag{2.135}
$$

**y-directed magnetic dipole:**

$$
\vec{E}^{y}_m = \frac{k}{\sqrt{\eta}}Q_{101}\,\tfrac{i\sqrt{2}}{2}\{\vec{F}^{(3)}_{1,-1,1} + \vec{F}^{(3)}_{111}\} \tag{2.136}
$$

$$
\vec{H}^{y}_m = -ik\sqrt{\eta}\,Q_{101}\,\tfrac{i\sqrt{2}}{2}\{\vec{F}^{(3)}_{2,-1,1} + \vec{F}^{(3)}_{211}\} \tag{2.137}
$$

---

## 7. Scattering Matrices for Electric and Magnetic Dipoles (§2.3.4)

> All dipoles assumed lossless and matched ($\Gamma = 0$, $|\mathbf{T}|^2 = 1$).
> Matrices indexed with the single-index convention (2.27). Layout: full $(J+1)\times(J+1)$ block of $\hat{\mathbf{S}}$, columns/rows $0\dots 7$ (where row/column 0 corresponds to the local-port reflection / receiving / transmitting; rows/cols 1..7 are the first seven spherical-mode ports = all $n=1$ modes plus the first $n=2$ mode).
>
> **Index reminder for $n=1$:** $j=1: (1,-1,1)$, $j=2: (2,-1,1)$, $j=3: (1,0,1)$, $j=4: (2,0,1)$, $j=5: (1,1,1)$, $j=6: (2,1,1)$, $j=7: (1,-2,2)$.

### z-directed electric dipole — Eq. (2.148)

Only mode $j=4$ ($s=2, m=0, n=1$) couples to the port.
- $T_4 = 1$, all other $T_i = 0$.
- $R_4 = 1$, all other $R_j = 0$ (consistent with reciprocity (2.107)).
- $S_{ij} = \delta_{ij}$ for $i,j \ne 4$; $S_{4j} = S_{i4} = 0$.

$$\hat{\mathbf{S}}^{z}_e = \begin{bmatrix}
0 & 0 & 0 & 0 & 1 & 0 & 0 & \cdots\\
0 & 1 & 0 & 0 & 0 & 0 & 0\\
0 & 0 & 1 & 0 & 0 & 0 & 0\\
0 & 0 & 0 & 1 & 0 & 0 & 0\\
1 & 0 & 0 & 0 & 0 & 0 & 0\\
0 & 0 & 0 & 0 & 0 & 1 & 0\\
0 & 0 & 0 & 0 & 0 & 0 & 1\\
\vdots & & & & & & & \ddots
\end{bmatrix}$$

### x-directed electric dipole — Eq. (2.154)

Non-zero $\mathbf{T}$: $T_2 = \sqrt{2}/2$, $T_6 = -\sqrt{2}/2$.
Non-zero $\mathbf{R}$: $R_2 = \sqrt{2}/2$, $R_6 = -\sqrt{2}/2$.
Non-zero $\mathbf{S}$ entries (within $n=1$ block): $S_{22} = S_{26} = S_{62} = S_{66} = 1/2$.
All other entries: $S_{ii} = 1$ outside the $\{2,6\}$ subspace; rest zero.

$$\hat{\mathbf{S}}^{x}_e = \begin{bmatrix}
0 & 0 & \tfrac{\sqrt 2}{2} & 0 & 0 & 0 & -\tfrac{\sqrt 2}{2} & 0 & \cdots\\
0 & 1 & 0 & 0 & 0 & 0 & 0 & 0\\
\tfrac{\sqrt 2}{2} & 0 & \tfrac{1}{2} & 0 & 0 & 0 & \tfrac{1}{2} & 0\\
0 & 0 & 0 & 1 & 0 & 0 & 0 & 0\\
0 & 0 & 0 & 0 & 1 & 0 & 0 & 0\\
0 & 0 & 0 & 0 & 0 & 1 & 0 & 0\\
-\tfrac{\sqrt 2}{2} & 0 & \tfrac{1}{2} & 0 & 0 & 0 & \tfrac{1}{2} & 0\\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 1\\
\vdots & & & & & & & & \ddots
\end{bmatrix}$$

### y-directed electric dipole — Eq. (2.155)

Non-zero $\mathbf{T}$: $T_2 = i\sqrt{2}/2$, $T_6 = i\sqrt{2}/2$.
Non-zero $\mathbf{R}$: $R_2 = -i\sqrt{2}/2$, $R_6 = -i\sqrt{2}/2$.
Non-zero $\mathbf{S}$ in $\{2,6\}$ block: $S_{22} = 1/2$, $S_{26} = -1/2$, $S_{62} = -1/2$, $S_{66} = 1/2$.

$$\hat{\mathbf{S}}^{y}_e = \begin{bmatrix}
0 & 0 & -\tfrac{i\sqrt 2}{2} & 0 & 0 & 0 & -\tfrac{i\sqrt 2}{2} & 0 & \cdots\\
0 & 1 & 0 & 0 & 0 & 0 & 0 & 0\\
\tfrac{i\sqrt 2}{2} & 0 & \tfrac{1}{2} & 0 & 0 & 0 & -\tfrac{1}{2} & 0\\
0 & 0 & 0 & 1 & 0 & 0 & 0 & 0\\
0 & 0 & 0 & 0 & 1 & 0 & 0 & 0\\
0 & 0 & 0 & 0 & 0 & 1 & 0 & 0\\
\tfrac{i\sqrt 2}{2} & 0 & -\tfrac{1}{2} & 0 & 0 & 0 & \tfrac{1}{2} & 0\\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 1\\
\vdots & & & & & & & & \ddots
\end{bmatrix}$$

### z-directed magnetic dipole — Eq. (2.156)

Only mode $j=3$ ($s=1, m=0, n=1$) couples. $T_3 = -i$, $R_3 = -i$, $S_{33} = 0$.

$$\hat{\mathbf{S}}^{z}_m = \begin{bmatrix}
0 & 0 & 0 & -i & 0 & 0 & \cdots\\
0 & 1 & 0 & 0 & 0 & 0\\
0 & 0 & 1 & 0 & 0 & 0\\
-i & 0 & 0 & 0 & 0 & 0\\
0 & 0 & 0 & 0 & 1 & 0\\
0 & 0 & 0 & 0 & 0 & 1\\
\vdots & & & & & & \ddots
\end{bmatrix}$$

### x-directed magnetic dipole — Eq. (2.157)

Non-zero $\mathbf{T}$: $T_1 = -i\sqrt{2}/2$, $T_5 = i\sqrt{2}/2$.
$\{1,5\}$ block: $S_{11} = S_{55} = 1/2$, $S_{15} = S_{51} = 1/2$.

$$\hat{\mathbf{S}}^{x}_m = \begin{bmatrix}
0 & -\tfrac{i\sqrt 2}{2} & 0 & 0 & 0 & \tfrac{i\sqrt 2}{2} & 0 & \cdots\\
-\tfrac{i\sqrt 2}{2} & \tfrac{1}{2} & 0 & 0 & 0 & \tfrac{1}{2} & 0\\
0 & 0 & 1 & 0 & 0 & 0 & 0\\
0 & 0 & 0 & 1 & 0 & 0 & 0\\
0 & 0 & 0 & 0 & 1 & 0 & 0\\
\tfrac{i\sqrt 2}{2} & \tfrac{1}{2} & 0 & 0 & 0 & \tfrac{1}{2} & 0\\
0 & 0 & 0 & 0 & 0 & 0 & 1\\
\vdots & & & & & & & \ddots
\end{bmatrix}$$

### y-directed magnetic dipole — Eq. (2.158)

Non-zero $\mathbf{T}$: $T_1 = \sqrt{2}/2$, $T_5 = \sqrt{2}/2$.
$\{1,5\}$ block: $S_{11} = 1/2$, $S_{15} = -1/2$, $S_{51} = -1/2$, $S_{55} = 1/2$.

$$\hat{\mathbf{S}}^{y}_m = \begin{bmatrix}
0 & -\tfrac{\sqrt 2}{2} & 0 & 0 & 0 & -\tfrac{\sqrt 2}{2} & 0 & \cdots\\
\tfrac{\sqrt 2}{2} & \tfrac{1}{2} & 0 & 0 & 0 & -\tfrac{1}{2} & 0\\
0 & 0 & 1 & 0 & 0 & 0 & 0\\
0 & 0 & 0 & 1 & 0 & 0 & 0\\
0 & 0 & 0 & 0 & 1 & 0 & 0\\
\tfrac{\sqrt 2}{2} & -\tfrac{1}{2} & 0 & 0 & 0 & \tfrac{1}{2} & 0\\
0 & 0 & 0 & 0 & 0 & 0 & 1\\
\vdots & & & & & & & \ddots
\end{bmatrix}$$

### Combined elements (Huygens, turnstile)

**z-directed Huygens source** = $\hat{x}$-electric + $\hat{y}$-magnetic-dual dipoles (the magnetic dipole is the dual source of the electric one). Total matrix (Eq. 2.159):

$$\hat{\mathbf{S}}^{z}_H = \begin{bmatrix}
0 & -\tfrac{1}{2} & \tfrac{1}{2} & 0 & 0 & -\tfrac{1}{2} & -\tfrac{1}{2} & 0 & \cdots\\
\tfrac{1}{2} & \tfrac{3}{4} & \tfrac{1}{4} & 0 & 0 & -\tfrac{1}{4} & -\tfrac{1}{4} & 0\\
\tfrac{1}{2} & -\tfrac{1}{4} & \tfrac{1}{4} & 0 & 0 & -\tfrac{1}{4} & \tfrac{3}{4} & 0\\
0 & 0 & 0 & 1 & 0 & 0 & 0 & 0\\
0 & 0 & 0 & 0 & 1 & 0 & 0 & 0\\
\tfrac{1}{2} & -\tfrac{1}{4} & \tfrac{1}{4} & 0 & 0 & \tfrac{3}{4} & -\tfrac{1}{4} & 0\\
-\tfrac{1}{2} & \tfrac{1}{4} & \tfrac{3}{4} & 0 & 0 & \tfrac{1}{4} & \tfrac{1}{4} & 0\\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 1\\
\vdots & & & & & & & & \ddots
\end{bmatrix}$$

**z-directed turnstile** = $\hat{x}$-electric + $\hat{y}$-electric dipoles in phase quadrature. Total matrix (Eq. 2.160):

$$\hat{\mathbf{S}}^{z}_T = \begin{bmatrix}
0 & 0 & 1 & 0 & 0 & 0 & 0 & 0 & \cdots\\
0 & 1 & 0 & 0 & 0 & 0 & 0 & 0\\
0 & 0 & 0 & 0 & 0 & 0 & 1 & 0\\
0 & 0 & 0 & 1 & 0 & 0 & 0 & 0\\
0 & 0 & 0 & 0 & 1 & 0 & 0 & 0\\
0 & 0 & 0 & 0 & 0 & 1 & 0 & 0\\
-1 & 0 & 0 & 0 & 0 & 0 & 0 & 0\\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 1\\
\vdots & & & & & & & & \ddots
\end{bmatrix}$$

Only $T_6$ and $R_2$ are non-zero; corresponds to far field $\vec{F}_6 = \vec{F}^{(3)}_{211}$ only (RHCP at $\theta = 0$, LHCP at $\theta = \pi$).

### Composition rule

Receiving/transmitting elements of an antenna built from multiple dipoles = sum of the individual elements weighted by excitation. Scattering elements **cannot** be simply summed — must be derived from incident-field considerations.

### Receiving-formula special case (z-electric dipole, matched load, $\Gamma_l = 0$)

$$
w = a_4 = \frac{\sqrt{6\pi}}{2}\,\frac{\sqrt{\eta}}{k}\,E_z(0,\theta,\phi) \tag{2.152}
$$

Scattered field for that dipole: $b_i = a_i$ for $i \ne 4$, $b_4 = 0$ (Eq. 2.153).

---

## 8. Alternative ("Source") Scattering Matrix — Yaghjian (§2.3.5)

Same waveguide modes, but use $\vec{F}^{(1)}_j$ (standing) instead of $\vec{F}^{(4)}_j$ for the second basis:

$$
\vec{E} = \frac{k}{\sqrt{\eta}}\sum_j \bigl\{a'_j\,\vec{F}^{(1)}_j + b'_j\,\vec{F}^{(3)}_j\bigr\} \tag{2.164}
$$

Using $\vec{F}^{(1)}_j = \tfrac{1}{2}\{\vec{F}^{(3)}_j + \vec{F}^{(4)}_j\}$ (2.165):

$$
\begin{bmatrix}\Gamma' & \mathbf{R}'\\\mathbf{T}' & \mathbf{S}'\end{bmatrix} = \begin{bmatrix}\Gamma & \tfrac{1}{2}\mathbf{R}\\\mathbf{T} & \tfrac{1}{2}(\mathbf{S} - \mathbf{I})\end{bmatrix} \tag{2.166}
$$

Implications:
- $\Gamma' = \Gamma$, $\mathbf{T}' = \mathbf{T}$.
- $\mathbf{R}' = \tfrac{1}{2}\mathbf{R}$.
- Empty space: $\mathbf{S}' = \mathbf{0}$ (vs. $\mathbf{S} = \mathbf{I}$ in the classical form).

For pure scatterers (no local port), $\mathbf{S}'$ reduces to Waterman's T-matrix [25]. The book retains the classical formulation throughout.

---

## 9. Far-Field Patterns (§2.4.1)

### Definition

$$
\vec{K}_{smn}(\theta,\phi) = \lim_{kr\to\infty}\Bigl[\sqrt{4\pi}\,\frac{kr}{e^{ikr}}\,\vec{F}^{(3)}_{smn}(r,\theta,\phi)\Bigr] \tag{2.175}
$$

### Explicit pattern functions

$$
\boxed{\;\vec{K}_{1mn}(\theta,\phi) = \sqrt{\frac{2}{n(n+1)}}\Bigl(-\frac{m}{|m|}\Bigr)^{m} e^{im\phi}(-i)^{n+1}\left\{\frac{im\,\bar{P}_n^{|m|}(\cos\theta)}{\sin\theta}\,\hat{\theta} - \frac{d\bar{P}_n^{|m|}(\cos\theta)}{d\theta}\,\hat{\phi}\right\}\;} \tag{2.176}
$$

$$
\boxed{\;\vec{K}_{2mn}(\theta,\phi) = \sqrt{\frac{2}{n(n+1)}}\Bigl(-\frac{m}{|m|}\Bigr)^{m} e^{im\phi}(-i)^{n}\left\{\frac{d\bar{P}_n^{|m|}(\cos\theta)}{d\theta}\,\hat{\theta} + \frac{im\,\bar{P}_n^{|m|}(\cos\theta)}{\sin\theta}\,\hat{\phi}\right\}\;} \tag{2.177}
$$

### Useful identity

$$
\vec{K}_{smn} = i\,\hat{r}\times\vec{K}_{3-s,m,n} \tag{2.178}
$$

### Far-field expressions

$$
\vec{E}(r,\theta,\phi) \to \frac{k}{\sqrt{\eta}}\,\frac{1}{\sqrt{4\pi}}\,\frac{e^{ikr}}{kr}\sum_{smn} Q^{(3)}_{smn}\,\vec{K}_{smn}(\theta,\phi) = \frac{k}{\sqrt{\eta}}\,\frac{1}{\sqrt{4\pi}}\,\frac{e^{ikr}}{kr}\,v\sum_{smn} T_{smn}\,\vec{K}_{smn}(\theta,\phi) = \frac{k}{\sqrt{\eta}}\,\frac{1}{\sqrt{4\pi}}\,\frac{e^{ikr}}{kr}\,v\,\vec{K}(\theta,\phi) \tag{2.179,\,2.180}
$$

$$
\vec{H} \to \eta\,\hat{r}\times\vec{E} = k\sqrt{\eta}\,\frac{1}{\sqrt{4\pi}}\,\frac{e^{ikr}}{kr}\,v\,\hat{r}\times\vec{K}(\theta,\phi) \tag{2.181}
$$

### Absolute far-field pattern

$$
\vec{K}(\theta,\phi) = \sum_{smn} T_{smn}\,\vec{K}_{smn}(\theta,\phi) \tag{2.182}
$$

Dimensionless. $C\vec{K}(\theta,\phi)$ for arbitrary $C$ is a **relative far-field pattern**.

---

## 10. Polarization (§2.4.2)

### Orthogonal polarization unit vector

$$
\hat{i}_{\text{cross}}(\theta,\phi) = \hat{r}\times\hat{i}_{co}^{*}(\theta,\phi) \tag{2.183}
$$

with $\hat{i}_{co}\cdot\hat{i}_{\text{cross}}^{*} = 0$ (2.184).

### Decomposition

$$
K_{co}(\theta,\phi) = \vec{K}\cdot\hat{i}_{co}^{*} = \sum_{smn} T_{smn}\,\vec{K}_{smn}\cdot\hat{i}_{co}^{*} \tag{2.185, 2.186}
$$

$$
K_{\text{cross}}(\theta,\phi) = \vec{K}\cdot\hat{i}_{\text{cross}}^{*} \tag{2.187, 2.188}
$$

$$
\vec{K} = K_{co}\,\hat{i}_{co} + K_{\text{cross}}\,\hat{i}_{\text{cross}} \tag{2.189}
$$

### Ludwig's "Definition 3" — Linear (boresight along $+z$, reference angle $\phi_o$)

$$
\hat{i}_{co,3L}(\theta,\phi) = \hat{\theta}\cos(\phi-\phi_o) - \hat{\phi}\sin(\phi-\phi_o),\quad 0\le\theta<\pi \tag{2.190}
$$

$$
\hat{i}_{\text{cross},3L}(\theta,\phi) = \hat{\theta}\sin(\phi-\phi_o) + \hat{\phi}\cos(\phi-\phi_o) \tag{2.191}
$$

With $\phi_o = 0$: identical to a $\hat{x}$-electric + $\hat{y}$-magnetic Huygens source.

### Circular polarization unit vectors ($\phi_o = 0$)

**Right-hand (RCP):**

$$
\hat{i}_{co,RC}(\theta,\phi) = \tfrac{1}{\sqrt 2}\bigl[\hat{i}_{co,3L} + i\,\hat{i}_{\text{cross},3L}\bigr]_{\phi_o=0} = \tfrac{1}{\sqrt 2}\,e^{i\phi}(\hat{\theta} + i\hat{\phi}) \tag{2.192}
$$

$$
\hat{i}_{\text{cross},RC}(\theta,\phi) = \tfrac{i}{\sqrt 2}\,e^{-i\phi}(\hat{\theta} - i\hat{\phi}) \tag{2.193}
$$

**Left-hand (LCP):**

$$
\hat{i}_{co,LC}(\theta,\phi) = \tfrac{1}{\sqrt 2}\,e^{-i\phi}(\hat{\theta} - i\hat{\phi}) \tag{2.194}
$$

$$
\hat{i}_{\text{cross},LC}(\theta,\phi) = \tfrac{-i}{\sqrt 2}\,e^{i\phi}(\hat{\theta} + i\hat{\phi}) \tag{2.195}
$$

The factors $i$ may be dropped in practice but the $e^{\pm i\phi}$ factors must remain for continuity at $\theta = 0$. All unit-vector distributions are discontinuous at $\theta = \pi$.

### Polarization-ellipse parameters

Decomposition (without the cross-polar $i$):

$$
\vec{K}(\theta,\phi) = K_R\,\hat{i}_{co,RC} + K_L\,\hat{i}_{co,LC} \tag{2.196}
$$

$$
K_R = |K_R|\,e^{i\psi_R},\quad K_L = |K_L|\,e^{i\psi_L} \tag{2.197, 2.198}
$$

$$
Q = K_R / K_L \tag{2.199}
$$

**Axial ratio** $r = \tan\alpha$:

$$
\tan\alpha = \frac{|K_R| - |K_L|}{|K_R| + |K_L|} = \frac{|Q| - 1}{|Q| + 1},\quad -\tfrac{\pi}{4}\le\alpha\le\tfrac{\pi}{4} \tag{2.200, 2.201}
$$

$\alpha = 0$ ↔ linear; $\alpha > 0$ ↔ RH-elliptical; $\alpha < 0$ ↔ LH-elliptical.

**Tilt angle** $\beta$ (relative to $\hat{i}_{co,3L}|_{\phi_o=0}$):

$$
\beta = \frac{\psi_L - \psi_R}{2} = -\tfrac{1}{2}\arg(Q) \tag{2.202, 2.203}
$$

---

## 11. Directivity and Gain (§2.4.3)

### Power per unit solid angle (far field)

$$
r^2\,\tfrac{1}{2}\,\mathrm{Re}\{\vec{E}\times\vec{H}^*\}\cdot\hat{r} = r^2\,\tfrac{1}{2}\eta|\vec{E}|^2 = \frac{1}{2}\,\frac{1}{4\pi}\Bigl|\sum_{smn}Q^{(3)}_{smn}\vec{K}_{smn}(\theta,\phi)\Bigr|^2 \tag{2.204, 2.205}
$$

### Isotropic reference

$$
\frac{P}{4\pi} = \frac{1}{4\pi}\,\frac{1}{2}\sum_{smn}|Q^{(3)}_{smn}|^2 \tag{2.206}
$$

### Directivity

$$
\boxed{\;D(\theta,\phi) = \frac{\Bigl|\sum_{smn} Q^{(3)}_{smn}\,\vec{K}_{smn}(\theta,\phi)\Bigr|^2}{\sum_{smn}|Q^{(3)}_{smn}|^2}\;} \tag{2.207}
$$

Equivalently (since $Q^{(3)}_{smn} = b_{smn} = v\,T_{smn}$):

$$
D(\theta,\phi) = \frac{\Bigl|\sum_{smn} T_{smn}\,\vec{K}_{smn}(\theta,\phi)\Bigr|^2}{\sum_{smn}|T_{smn}|^2} = \frac{|\vec{K}(\theta,\phi)|^2}{\sum_{smn}|T_{smn}|^2} = |\vec{K}(\theta,\phi)|^2 \tag{2.208}
$$

The last equality requires **matched and lossless** ($\sum |T_{smn}|^2 = 1$).

### Gain

Input power $P_{\text{in}} = \tfrac{1}{2}|v|^2(1 - |\Gamma|^2)$ (2.209).

$$
\boxed{\;G(\theta,\phi) = \frac{\Bigl|\sum_{smn} Q^{(3)}_{smn}\,\vec{K}_{smn}(\theta,\phi)\Bigr|^2}{|v|^2(1-|\Gamma|^2)} = \frac{|\vec{K}(\theta,\phi)|^2}{1-|\Gamma|^2} = \frac{|\vec{K}(\theta,\phi)|^2}{P_{\text{loss}}/P_{\text{inc}} + \sum|T_{smn}|^2}\;} \tag{2.210, 2.211, 2.212}
$$

### Polarization additivity

$D$ (or $G$) in any direction = sum of $D$ (or $G$) for any two orthogonal polarizations.

---

## 12. Maximum Directivity (§2.4.4)

### Component directivities

$$
D(\theta,\phi) = D_{co}(\theta,\phi) + D_{\text{cross}}(\theta,\phi) \tag{2.213}
$$

$$
D_{co}(\theta,\phi) = \frac{|\sum_{smn} Q_{smn}\,\vec{K}_{smn}(\theta,\phi)\cdot\hat{i}_{co}^{*}|^2}{\sum_{smn}|Q_{smn}|^2} = |\vec{K}(\theta,\phi)\cdot\hat{i}_{co}^{*}|^2 \tag{2.214}
$$

$$
D_{\text{cross}}(\theta,\phi) = |\vec{K}(\theta,\phi)\cdot\hat{i}_{\text{cross}}^{*}|^2 \tag{2.215}
$$

### Cauchy–Schwartz bound

$$
D_{co}(\theta',\phi') \le \sum_{smn}|\vec{K}_{smn}(\theta',\phi')\cdot\hat{i}_{co}^{*}|^2 = D_{co,\max}(\theta',\phi') \tag{2.216, 2.217}
$$

Equality (i.e. the maximum) is achieved by:

$$
\boxed{\;Q_{smn} = c\,\bigl(\vec{K}_{smn}(\theta',\phi')\cdot\hat{i}_{co}^{*}\bigr)^*\;} \tag{2.218}
$$

with arbitrary constant $c$.

### Maximum directivity value (independent of direction and polarization)

For $\hat{i}_{co} = \alpha\hat{\theta} + \beta\hat{\phi}$, $|\alpha|^2 + |\beta|^2 = 1$ (2.219, 2.220):

$$
\boxed{\;D_{co,\max}(\theta',\phi') = N^2 + 2N\;} \tag{2.225, 2.226}
$$

where $N$ is the truncation in $n$.

### Coefficients for $\hat{x}$-polarized peak at $(\theta',\phi')=(0,0)$

$$
Q_{1,1,n} = Q_{1,-1,n} = Q_{2,1,n} = -Q_{2,-1,n} = c\,(-i^{n})\,\tfrac{1}{2}\sqrt{2n+1} \tag{2.227}
$$

All other coefficients zero. Cross-polarization is zero in every direction.

> These are proportional to the coefficients for an $\hat{x}$-polarized plane wave travelling along $+z$, but evaluated with $\vec{K}_{smn}$ rather than $\vec{F}^{(1)}_{smn}$.

### Reference values

| $N$ | $D_{\max} = N^2+2N$ | $D_{\max}$ (dB) |
|---|---|---|
| 1 | 3   | 4.77 |
| 2 | 8   | 9.03 |
| 3 | 15  | 11.76 |
| 4 | 24  | 13.80 |
| 10 | 120 | 20.79 |
| 20 | 440 | 26.43 |

> The Huygens source ($N=1$, $D=3$) is the simplest maximum-directivity antenna.

---

## 13. Implementation Checklist for Python Code

Treat each of these as a unit test or assertion:

1. **Time convention.** Code must use $e^{-i\omega t}$. Outgoing modes carry $h_n^{(1)} \sim e^{+ikr}/(kr)$. *Validate with `(2.13)`–`(2.16)`.*
2. **Hankel-function identity.** `j_n(z) == 0.5 * (h_n^(1)(z) + h_n^(2)(z))` — Eq. (2.30).
3. **Wronskian.** `R1n^(1) * R2n^(2) - R2n^(1) * R1n^(2) == 1/(kr)^2` — Eq. (2.52).
4. **$\bar P_n^{|m|}$ normalization.** Use Belousov's convention (not Schmidt, not unnormalized). Spot-check $\int_{-1}^{1}\bar P_n^m(\mu)^2 d\mu = 2/(2n+1)$ × normalization factor consistent with Eq. (2.18).
5. **Phase factor.** `(-m/abs(m))**m == 1` when `m == 0` — Eq. (2.19). For $m \ne 0$: equals $(-1)^m$ if $m > 0$, $(-1)^{|m|}\cdot(\text{sign})$ if $m < 0$. Reduces to $1$ for $m=0$ and to $\pm 1$ otherwise per the Edmonds convention.
6. **Index swap.** $\vec{H}$ uses $\vec{F}^{(c)}_{3-s,m,n}$ — Eq. (2.23). A common bug is mirroring $\vec{E}$ in $\vec{H}$.
7. **Single index.** `j = 2*(n*(n+1) + m - 1) + s` — Eq. (2.27). Spot-check the $n=1$ table above.
8. **Power normalization.** Pure $c=3$ field with one mode of unit $|Q|$ radiates $0.5$ W — Eq. (2.24).
9. **Truncation.** `N = floor(k*r0) + n1`, default `n1 = 10` — Eq. (2.31).
10. **Reciprocity.** For reciprocal antennas: `R[s,m,n] == (-1)**m * T[s,-m,n]` — Eq. (2.107).
11. **Unitarity (lossless).** `S_hat^H @ S_hat == I` and `|Γ|^2 + sum |T|^2 == 1` — Eq. (2.61).
12. **Empty space.** `S == I` (classical) but `S' == 0` (source-form) — §2.3.1 / Eq. (2.166).
13. **z-electric dipole.** Only $T_4 = 1$, $R_4 = 1$ — Eq. (2.148). Sanity check via `Q_201 = -k*d_e/(sqrt(6*pi)*sqrt(eta))` — Eq. (2.117).
14. **z-magnetic dipole.** $T_3 = -i$, $R_3 = -i$ — Eq. (2.156). $Q_{101} = -i\,Q_{201}$ when $d_m = -\zeta d_e$ — Eq. (2.142).
15. **Far-field K identity.** `K_smn = i * r_hat × K_{3-s,m,n}` — Eq. (2.178).
16. **Directivity self-consistency.** When the antenna is matched and lossless: `D(θ,φ) == |K(θ,φ)|^2` — Eq. (2.208).
17. **Max directivity.** A spherical-wave field truncated to $N$ has $D_{\max} = N^2 + 2N$ — Eq. (2.226). Test by maximizing the Cauchy–Schwartz inequality with coefficients (2.227).
18. **Adjoint vs reciprocal.** For non-reciprocal antennas, the $(-1)^m T_{s,-m,n}$ relation gives the **adjoint** $R'$, not $R$ — Eqs. (2.104) vs. (2.107).
19. **Scattered field.** Use $\mathbf{b}' = (\mathbf{S}-\mathbf{I})\mathbf{a}$ under matched load — Eq. (2.78). Don't forget to subtract $\mathbf{I}$.
20. **Axial ratio sign.** `α > 0` for right-handed, `α < 0` for left-handed, under the $e^{-i\omega t}$ time convention — §2.4.2.

---

*End of Chapter 2 reference.*
