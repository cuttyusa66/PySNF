# Chapter 2 — Scattering Matrix Description of an Antenna

**Source:** Hansen, *Spherical Near-Field Antenna Measurements*, Chapter 2 (pp. 8–60).
**Cross-reference:** Symbols defined in [A5_List_of_principal_symbols_and_uses.md](A5_List_of_principal_symbols_and_uses.md).

> **Purpose of this document:** structured reference for validating Python implementations of spherical-wave / scattering-matrix antenna code. Every equation is tagged with its book number `(2.N)` so that test cases and code comments can cite them directly.

---

## 0. Global Conventions (must be honoured by any implementation)

| Convention | Value | Comment |
|---|---|---|
| Time factor | $\\exp(-i\\omega t)$ | Suppressed throughout. **Sign flip = wrong sign of $i$ in every radial-function expression.** |
| Imaginary unit | $i$ | Not $j$. Outgoing wave $\\propto e^{+ikr}$ under this convention. |
| Spherical coords | $r, \\theta, \\phi$ | $0 \\le r < \\infty$, $0 \\le \\theta \\le \\pi$, $0 \\le \\phi < 2\\pi$ |
| Wavenumber | $k = \\omega\\sqrt{\\mu\\varepsilon} = 2\\pi/\\lambda$ | Real for loss-free; complex allowed (except power & directivity formulas). |
| Specific admittance | $\\eta = \\sqrt{\\varepsilon/\\mu}$ | Used in field-coefficient prefactors. |
| Specific impedance | $\\zeta = \\sqrt{\\mu/\\varepsilon}$ | Equals $1/\\eta$. |
| Index $s$ | $s \\in \\{1, 2\\}$ | $s=1$ ↔ TE in coefficient $Q$; $s=2$ ↔ TM in coefficient $Q$. **In the wave-function $\\vec{F}^{(c)}\_{smn}$ the meaning is reversed**: $\\vec{F}^{(c)}\_{1mn}$ is transverse (no radial component); $\\vec{F}^{(c)}\_{2mn}$ carries the radial component. |
| Index $n$ | $n = 1, 2, 3, \\dots$ | No $n=0$. |
| Index $m$ | $m = -n, -n+1, \\dots, n-1, n$ | $\\|m\\| \\le n$. |
| Index $c$ | $c \\in \\{1, 2, 3, 4\\}$ | $c=1$: $j\_n$ (finite at origin). $c=2$: $n\_n$ (singular at origin). $c=3$: $h\_n^{(1)}$ (outgoing). $c=4$: $h\_n^{(2)}$ (incoming). |
| Power normalization | A single $c=3$ mode with $\\|Q^{(3)}\_{smn}\\| = 1$ radiates $\\tfrac{1}{2}$ W. | The factor $k/\\sqrt{\\eta}$ in $\\vec{E}$ is what makes this hold. |
| Wave-coefficient units | $[Q^{(c)}\_{smn}] = \\text{W}^{1/2}$ | The functions $\\vec{F}^{(c)}\_{smn}$ are dimensionless. |
| Waveguide-port units | $[v] = [w] = \\text{W}^{1/2}$ | Incident power at local port is $\\tfrac{1}{2}\\|v\\|^2$. |
| Convention $(-m/\\|m\\|)^m$ | $= 1$ when $m = 0$ (Eq. 2.19). | Ensures the Edmonds phase. |

---

## 1. Maxwell's Equations and the Vector Wave Equation (§2.2.1)

In a linear, isotropic, homogeneous medium with assumed sources $\\vec{J}$ (electric), $\\vec{M}$ (magnetic):

|  Eq.  |   |
| :---: | :-- |
| (2.1) | $\\nabla \\times \\vec{H} = -i\\omega\\varepsilon \\vec{E} + \\vec{J}$ |
| (2.2) | $\\nabla \\times \\vec{E} = i\\omega\\mu \\vec{H} - \\vec{M}$ |

In a source-free region both $\\vec{E}$ and $\\vec{H}$ satisfy:

|  Eq.  |   |
| :---: | :-- |
| (2.3) | $\\nabla \\times (\\nabla \\times \\vec{C}) - k^2 \\vec{C} = 0$ |

### Generating function and Hansen vector functions

The scalar Helmholtz equation $(\\nabla^2 + k^2) f = 0$ (Eq. 2.4) generates:

|  Eq.  |   |
| :---: | :-- |
| (2.5) | $\\vec{m} = \\nabla f \\times \\vec{r}$ |
| (2.6) | $\\vec{n} = k^{-1} \\nabla \\times \\vec{m}$ |

Reciprocal curl relations:

|  Eq.  |   |
| :---: | :-- |
| (2.7) | $\\vec{m} = k^{-2} \\nabla \\times (\\nabla \\times \\vec{m})$ |
| (2.8) | $\\vec{m} = k^{-1} \\nabla \\times \\vec{n}$ |

> The third (irrotational) Hansen function $\\vec{l} = \\nabla f$ is not needed here.

### Stratton generating function

**(2.9)**

```math
f^{(c)}_{\sigma mn}(r,\theta,\phi) = z^{(c)}_n(kr)\, P_n^m(\cos\theta)\, \begin{cases}\cos m\phi\\\sin m\phi\end{cases}
```

(With $\\sigma \\in \\{e, o\\}$ for even/odd trig.) $P\_n^m$ is the unnormalized associated Legendre function.

### Radial functions $z\_n^{(c)}(kr)$

| $c$ | $z\_n^{(c)}$ | Wave type |
|---|---|---|
| 1 | $j\_n(kr)$ — spherical Bessel | Standing, finite at origin (2.10a) |
| 2 | $n\_n(kr)$ — spherical Neumann | Standing, singular at origin (2.10b) |
| 3 | $h\_n^{(1)}(kr) = j\_n(kr) + i\\, n\_n(kr)$ | **Outward travelling** (2.10c) |
| 4 | $h\_n^{(2)}(kr) = j\_n(kr) - i\\, n\_n(kr)$ | **Inward travelling** (2.10d) |

**Bessel identity (used for matching at origin):**

|  Eq.  |   |
| :---: | :-- |
| (2.30) | $j\_n(kr) = \\tfrac{1}{2}\\bigl(h\_n^{(1)}(kr) + h\_n^{(2)}(kr)\\bigr)$ |

### Large-argument (far-field) asymptotics for $kr \\to \\infty$, $kr \\gg n$

|  Eq.  |   |
| :---: | :-- |
| (2.13) | $z\_n^{(3)}(kr) \\to (-i)^{n+1} \\frac{e^{ikr}}{kr}$ |
| (2.14) | $z\_n^{(4)}(kr) \\to i^{n+1} \\frac{e^{-ikr}}{kr}$ |
| (2.15) | $\\frac{1}{kr}\\frac{d}{d(kr)}\\bigl\\{kr\\, z\_n^{(3)}(kr)\\bigr\\} \\to (-i)^{n}\\frac{e^{ikr}}{kr}$ |
| (2.16) | $\\frac{1}{kr}\\frac{d}{d(kr)}\\bigl\\{kr\\, z\_n^{(4)}(kr)\\bigr\\} \\to i^{n}\\frac{e^{-ikr}}{kr}$ |

---

## 2. Power-Normalized Spherical Wave Functions (§2.2.2)

### Modified scalar generating function (Jensen [13])

|  Eq.  |   |
| :---: | :-- |
| (2.18) | $F^{(c)}\_{mn}(r,\\theta,\\phi) = \\frac{1}{\\sqrt{2\\pi}}\\,\\frac{1}{\\sqrt{n(n+1)}}\\Bigl(-\\frac{m}{\|m\|}\\Bigr)^{m} z\_n^{(c)}(kr)\\, \\bar{P}\_n^{\|m\|}(\\cos\\theta)\\, e^{im\\phi}$ |

with $\\bar{P}\_n^{|m|}$ the **normalized** associated Legendre function (Belousov [14]; see Appendix A1).

### Vector wave functions

**$s=1$ (transverse — no radial component):**

|  Eq.  |   |
| :---: | :-- |
| (2.20) | $\\vec{F}^{(c)}\_{1mn}(r,\\theta,\\phi) = \\nabla F^{(c)}\_{mn} \\times \\vec{r} = \\frac{1}{\\sqrt{2\\pi}}\\frac{1}{\\sqrt{n(n+1)}}\\Bigl(-\\frac{m}{\|m\|}\\Bigr)^{m}\\left\\{ \\,z\_n^{(c)}(kr)\\,\\frac{im\\,\\bar{P}\_n^{\|m\|}(\\cos\\theta)}{\\sin\\theta}\\,e^{im\\phi}\\,\\hat{\\theta} -z\_n^{(c)}(kr)\\,\\frac{d\\bar{P}\_n^{\|m\|}(\\cos\\theta)}{d\\theta}\\,e^{im\\phi}\\,\\hat{\\phi} \\right\\}$ |

**$s=2$ (carries radial component):**

|  Eq.  |   |
| :---: | :-- |
| (2.21) | $\\vec{F}^{(c)}\_{2mn}(r,\\theta,\\phi) = k^{-1}\\nabla \\times \\vec{F}^{(c)}\_{1mn} = \\frac{1}{\\sqrt{2\\pi}}\\frac{1}{\\sqrt{n(n+1)}}\\Bigl(-\\frac{m}{\|m\|}\\Bigr)^{m}\\Biggl\\{ \\frac{n(n+1)}{kr}\\,z\_n^{(c)}(kr)\\,\\bar{P}\_n^{\|m\|}(\\cos\\theta)\\,e^{im\\phi}\\,\\hat{r} +\\,\\frac{1}{kr}\\frac{d}{d(kr)}\\bigl\\{kr z\_n^{(c)}(kr)\\bigr\\}\\,\\frac{d\\bar{P}\_n^{\|m\|}(\\cos\\theta)}{d\\theta}\\,e^{im\\phi}\\,\\hat{\\theta} +\\,\\frac{1}{kr}\\frac{d}{d(kr)}\\bigl\\{kr z\_n^{(c)}(kr)\\bigr\\}\\,\\frac{im\\,\\bar{P}\_n^{\|m\|}(\\cos\\theta)}{\\sin\\theta}\\,e^{im\\phi}\\,\\hat{\\phi} \\Biggr\\}$ |

> Both $\\vec{F}^{(c)}\_{1mn}$ and $\\vec{F}^{(c)}\_{2mn}$ are **dimensionless**.

### Field expansion (in a source-free region such as region 2 of Fig. 2.2)

|  Eq.  |   |
| :---: | :-- |
| (2.22) | $\\vec{E}(r,\\theta,\\phi) = \\frac{k}{\\sqrt{\\eta}}\\sum\_{c\\,s\\,m\\,n} Q^{(c)}\_{smn}\\,\\vec{F}^{(c)}\_{smn}(r,\\theta,\\phi)$ |
| (2.23) | $\\vec{H}(r,\\theta,\\phi) = (i\\omega\\mu)^{-1}\\nabla\\times\\vec{E} = -ik\\sqrt{\\eta}\\sum\_{c\\,s\\,m\\,n} Q^{(c)}\_{smn}\\,\\vec{F}^{(c)}\_{3-s,m,n}(r,\\theta,\\phi)$ |

**Note the index swap $s \\to 3-s$ in the magnetic field.** A correct implementation must mirror this.

### Summation conventions

|  Eq.  |   |
| :---: | :-- |
| (2.25) | $\\sum\_{smn} \\;=\\; \\sum\_{s=1}^{2}\\sum\_{n=1}^{\\infty}\\sum\_{m=-n}^{n}$ |
| (2.26) | $\\sum\_{csmn} \\;=\\; \\sum\_{c=3}^{4}\\sum\_{s=1}^{2}\\sum\_{n=1}^{\\infty}\\sum\_{m=-n}^{n}$ |

### Single-index transformation (when $M = N$)

|  Eq.  |   |
| :---: | :-- |
| (2.27) | $j = 2\\{n(n+1) + m - 1\\} + s$ |
| (2.28) | $\\sum\_{s=1}^{2}\\sum\_{n=1}^{N}\\sum\_{m=-n}^{n} = \\sum\_{j=1}^{J}$ |
| (2.29) | $J = 2N(N+2)$ |

Substitutions: $Q^{(c)}\_{smn} = Q^{(c)}\_j$, $\\vec{F}^{(c)}\_{smn} = \\vec{F}^{(c)}\_j$.

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

|  Eq.  |   |
| :---: | :-- |
| (2.24, 2.55) | $P = \\frac{1}{2}\\sum\_{smn}\|Q^{(3)}\_{smn}\|^2 \\quad\\text{[watts]}$ |

### Mode-type identification

- For a TM wave: $\\vec{H} \\propto \\vec{F}^{(c)}\_{1mn}$ (no radial), $\\vec{E} \\propto \\vec{F}^{(c)}\_{2mn}$. Coefficient $Q^{(c)}\_{2mn}$ ($s=2$).
- For a TE wave: $\\vec{E} \\propto \\vec{F}^{(c)}\_{1mn}$, $\\vec{H} \\propto \\vec{F}^{(c)}\_{2mn}$. Coefficient $Q^{(c)}\_{1mn}$ ($s=1$).
- Sum of wave-function $s$ indices for $\\vec{E}$ and $\\vec{H}$ of the same mode = 3.

### Truncation rule

Convergence is obtained for:

|  Eq.  |   |
| :---: | :-- |
| (2.31) | $N = \\lfloor kr\_0 \\rfloor + n\_1$ |

with $r\_0$ = minimum-sphere radius, $n\_1$ depending on accuracy/geometry. **Empirical default: $n\_1 = 10$** when the field point is more than a few wavelengths from the minimum sphere and 4 correct digits suffice. $n\_1$ grows roughly as $(kr\_0)^{1/3}$ for fixed accuracy.

For an interior (point closer than the source-free inner radius $r\_i$):

|  Eq.  |   |
| :---: | :-- |
| (2.32) | $N = \\lfloor kr \\rfloor + n\_1$ |

### Region definitions for a multi-mode antenna ($N$ = truncation)

| Region | Range |
|---|---|
| Evanescent | $r\_0 \\lesssim r \\lesssim N/k$ |
| Fresnel (near-field) | $N/k \\lesssim r \\lesssim 4N^2/(\\pi k)$ |
| Fraunhofer (far-field) | $4N^2/(\\pi k) \\lesssim r < \\infty$ |

### Rayleigh distance

|  Eq.  |   |
| :---: | :-- |
| (2.33), (2.34) | $R = \\frac{2D^2}{\\lambda} = \\frac{2}{\\lambda}(2 r\_0)^2 \\approx \\frac{2}{\\lambda}\\Bigl(2\\frac{N}{k}\\Bigr)^2 = \\frac{4}{\\pi}\\frac{N^2}{k}$ |

### On-axis behaviour (singular at poles, used in dipole derivations)

|  Eq.  |   |
| :---: | :-- |
| (2.35) | $[\\vec{F}^{(c)}\_{smn}(r,0,\\phi)]\_r = [\\vec{F}^{(c)}\_{smn}(r,\\pi,\\phi)]\_r = 0 \\quad\\text{for } m \\ne 0$ |
| (2.36) | $[\\vec{F}^{(c)}\_{smn}(r,0,\\phi)]\_{\\theta,\\phi} = [\\vec{F}^{(c)}\_{smn}(r,\\pi,\\phi)]\_{\\theta,\\phi} = 0 \\quad\\text{for } m \\ne \\pm 1$ |

---

## 3. Power Flow (§2.2.4)

### Complex Poynting vector

|  Eq.  |   |
| :---: | :-- |
| (2.37) | $\\vec{S} = \\tfrac{1}{2}\\,\\vec{E}\\times\\vec{H}^\*$ |

### Power radiated through a sphere of radius $r$

|  Eq.  |   |
| :---: | :-- |
| (2.38) | $P = \\int\_0^{2\\pi}\\!\\!\\int\_0^{\\pi} \\mathrm{Re}(\\hat{r}\\cdot\\vec{S})\\, r^2\\sin\\theta\\, d\\theta\\, d\\phi$ |

### Reactive-energy relation

|  Eq.  |   |
| :---: | :-- |
| (2.39) | $2\\omega(W\_e - W\_m) = \\int\_0^{2\\pi}\\!\\!\\int\_0^{\\pi} \\mathrm{Im}(\\hat{r}\\cdot\\vec{S})\\, r^2\\sin\\theta\\, d\\theta\\, d\\phi$ |

For a TE mode ($s=1$): $W\_m > W\_e$. For a TM mode ($s=2$): $W\_e > W\_m$.

### Far-field plane-wave relation

|  Eq.  |   |
| :---: | :-- |
| (2.40) | $\\vec{H} = \\eta\\,\\hat{r}\\times\\vec{E}$ |
| (2.41) | $\\vec{S} = \\tfrac{1}{2}\\eta\|\\vec{E}\|^2\\,\\hat{r} \\quad\\text{[W/m²]}$ |
| (2.42) | $P\_1(\\theta,\\phi) = \\tfrac{1}{2}\\eta\|\\vec{E}\|^2 r^2 \\quad\\text{[W/sr]}$ |

### Complex-conjugate identity for $\\vec{F}^{(c)}\_{smn}$

|  Eq.  |   |
| :---: | :-- |
| (2.45) | $\\vec{F}^{(3)\*}\_{smn}(r,\\theta,\\phi) = (-1)^m\\,\\vec{F}^{(4)}\_{s,-m,n}(r,\\theta,\\phi)$ |

### Orthogonality integral (used for power and reciprocity)

|  Eq.  |   |
| :---: | :-- |
| (2.46) | $\\int\_0^{2\\pi}\\!\\!\\int\_0^{\\pi} \\bigl\\{\\vec{F}^{(c)}\_{smn}\\times\\vec{F}^{(\\gamma)}\_{\\sigma\\mu\\nu}\\bigr\\}\\cdot\\hat{r}\\,\\sin\\theta\\,d\\theta\\,d\\phi = \\delta\_{3-s,\\sigma}\\,\\delta\_{m,-\\mu}\\,\\delta\_{n\\nu}\\,(-1)^{3-s}\\,(-1)^m\\,R^{(c)}\_{sn}(kr)\\,R^{(\\gamma)}\_{3-s,n}(kr)$ |

with the radial-function abbreviation

**(2.47)**

```math
R^{(c)}_{sn}(kr) = \begin{cases} z_n^{(c)}(kr), & s = 1 \\ \dfrac{1}{kr}\dfrac{d}{d(kr)}\{kr\,z_n^{(c)}(kr)\}, & s = 2 \end{cases}
```

### Kronecker delta and Wronskian

**(2.48)**

```math
\delta_{ij} = \begin{cases}0 & i \ne j \\ 1 & i = j\end{cases}
```

|  Eq.  |   |
| :---: | :-- |
| (2.50) | $R^{(3)}\_{sn}(kr) = R^{(1)}\_{sn}(kr) + i\\,R^{(2)}\_{sn}(kr)$ |
| (2.51) | $R^{(4)}\_{sn}(kr) = R^{(1)}\_{sn}(kr) - i\\,R^{(2)}\_{sn}(kr)$ |
| (2.52) | $R^{(1)}\_{1n}(kr)R^{(2)}\_{2n}(kr) - R^{(1)}\_{2n}(kr)R^{(2)}\_{1n}(kr) = (kr)^{-2}$ |

### Complex-power-flux result (outgoing $c=3$ field)

|  Eq.  |   |
| :---: | :-- |
| (2.53) | $\\tfrac{1}{2}\\int(\\vec{E}\\times\\vec{H}^\*)\\cdot\\hat{r}\\,r^2\\sin\\theta\\,d\\theta\\,d\\phi = \\sum\_{smn}\\Bigl\\{\\tfrac{1}{2} + \\tfrac{1}{2}i(-1)^{3-s}(kr)^2 V\_n(kr)\\Bigr\\}\\bigl\|Q^{(3)}\_{smn}\\bigr\|^2$ |

with cross-product

|  Eq.  |   |
| :---: | :-- |
| (2.54) | $V\_n(kr) = R^{(1)}\_{1n}(kr)R^{(2)}\_{1n}(kr) + R^{(1)}\_{2n}(kr)R^{(2)}\_{2n}(kr) = \\Bigl(\\frac{1}{kr} + \\frac{1}{2}\\frac{d}{d(kr)}\\Bigr)\|h\_n^{(1)}(kr)\|^2$ |

$V\_n(kr)$ is always negative (Abramowitz & Stegun 10.1.27).

**Real power:** independent of $r$, as in (2.55) above.

---

## 4. The Antenna Scattering Matrix (§2.3.1)

### Field outside the minimum sphere

|  Eq.  |   |
| :---: | :-- |
| (2.56) | $\\vec{E}(r,\\theta,\\phi) = \\frac{k}{\\sqrt{\\eta}}\\sum\_{j=1}^{J}\\bigl\\{a\_j\\,\\vec{F}^{(4)}\_j(r,\\theta,\\phi) + b\_j\\,\\vec{F}^{(3)}\_j(r,\\theta,\\phi)\\bigr\\},\\quad r > r\_0$ |

Index relation: $Q^{(4)}\_{smn} = a\_j = a\_{smn}$ (incoming), $Q^{(3)}\_{smn} = b\_j = b\_{smn}$ (outgoing), with $j$ from (2.27).

### Antenna scattering equation (order $J+1$)

**(2.57)**

```math
\begin{bmatrix}\Gamma & \mathbf{R}\\ \mathbf{T} & \mathbf{S}\end{bmatrix}\begin{bmatrix}v\\ \mathbf{a}\end{bmatrix} = \begin{bmatrix}w\\ \mathbf{b}\end{bmatrix}
```

Expanded:

|  Eq.  |   |
| :---: | :-- |
| (2.58) | $\\Gamma v + \\sum\_{j=1}^{J} R\_j a\_j = w$ |
| (2.59) | $T\_i v + \\sum\_{j=1}^{J} S\_{ij} a\_j = b\_i,\\quad i = 1, 2, \\dots, J$ |

Condensed: $\\hat{\\mathbf{S}}\\hat{\\mathbf{a}} = \\hat{\\mathbf{b}}$ (2.60), where $\\hat{\\mathbf{S}}$ is the **total** scattering matrix.

| Symbol | Shape | Meaning |
|---|---|---|
| $\\Gamma$ | scalar | Antenna reflection coefficient |
| $\\mathbf{R}$ | $1 \\times J$ row | Receiving coefficients $R\_j$ |
| $\\mathbf{T}$ | $J \\times 1$ column | Transmitting coefficients $T\_i$ |
| $\\mathbf{S}$ | $J \\times J$ | Scattering coefficients $S\_{ij}$ |

All elements of $\\hat{\\mathbf{S}}$ are dimensionless.

### Lossless antenna — unitarity

$\\hat{\\mathbf{S}}^{+}\\hat{\\mathbf{S}} = \\hat{\\mathbf{I}}$ (unit matrix of order $J+1$).

In particular, for the first column:

|  Eq.  |   |
| :---: | :-- |
| (2.61) | $\|\\Gamma\|^2 + \|\\mathbf{T}\|^2 = 1$ |

### Lossy antenna

|  Eq.  |   |
| :---: | :-- |
| (2.62) | $\\tfrac{1}{2}\\sum\_{i=1}^{J}\|b\_i\|^2 = \\tfrac{1}{2}\|v\|^2 - \\bigl(\\tfrac{1}{2}\|w\|^2 + P\_{\\text{loss}}\\bigr)$ |
| (2.63) | $\|\\Gamma\|^2 + \|\\mathbf{T}\|^2 = 1 - \\frac{P\_{\\text{loss}}}{P\_{\\text{inc}}}$ |

In the lossy case, every row and column of $\\hat{\\mathbf{S}}$ has norm $\\le 1$.

### Empty space

$\\hat{\\mathbf{S}}$ undefined, but $\\mathbf{S} = \\mathbf{I}$ (unit matrix of infinite order). Outgoing equals incoming: $\\mathbf{b} = \\mathbf{a}$.

### Generator coupling — radiated waves

Generator: $v = v\_g + \\Gamma\_g w$ (2.64). With $\\mathbf{a} = 0$:

|  Eq.  |   |
| :---: | :-- |
| (2.65) | $\\Gamma v = w$ |
| (2.66) | $\\mathbf{T} v = \\mathbf{b}$ |
| (2.67) | $\\boxed{\\;\\mathbf{b} = \\frac{v\_g}{1 - \\Gamma\_g \\Gamma}\\,\\mathbf{T}\\;}$ |
| (2.68) | $\\vec{E}(r,\\theta,\\phi) = \\frac{k}{\\sqrt{\\eta}}\\sum\_{i=1}^{J} b\_i\\,\\vec{F}^{(3)}\_i(r,\\theta,\\phi)$ |

### Load coupling — received wave & scattering

Load: $v = \\Gamma\_l w$ (2.69). With $\\Gamma v + \\mathbf{R}\\mathbf{a} = w$ (2.70):

|  Eq.  |   |
| :---: | :-- |
| (2.71) | $\\boxed{\\;w = \\frac{1}{1-\\Gamma\_l\\Gamma}\\,\\mathbf{R}\\mathbf{a}\\;}$ |

Power accepted by the load:

|  Eq.  |   |
| :---: | :-- |
| (2.72) | $P = \\tfrac{1}{2}(1-\|\\Gamma\_l\|^2)\|w\|^2 = \\tfrac{1}{2}(1-\|\\Gamma\_l\|^2)\\,\\frac{\|\\mathbf{R}\\mathbf{a}\|^2}{\|1-\\Gamma\_l\\Gamma\|^2}$ |

**Matched load** ($\\Gamma\_l = 0$): $P = P' = \\tfrac{1}{2}|\\mathbf{R}\\mathbf{a}|^2$ (2.73).

**Conjugate-matched load** ($\\Gamma\_l = \\Gamma^\*$): maximum, "available" power:

|  Eq.  |   |
| :---: | :-- |
| (2.74) | $P\_a = \\frac{1}{2}\\,\\frac{\|\\mathbf{R}\\mathbf{a}\|^2}{1-\|\\Gamma\|^2}$ |

### Scattered field

From $\\mathbf{T} v + \\mathbf{S}\\mathbf{a} = \\mathbf{b}$ (2.75) and (2.69), (2.71):

|  Eq.  |   |
| :---: | :-- |
| (2.76) | $\\mathbf{b} = \\{\\mathbf{T}\\Gamma\_l(1-\\Gamma\\Gamma\_l)^{-1}\\mathbf{R} + \\mathbf{S}\\}\\mathbf{a}$ |

**Subtract empty-space scattering** to get the field actually scattered by the antenna:

|  Eq.  |   |
| :---: | :-- |
| (2.77) | $\\mathbf{b}' = \\{\\mathbf{T}\\Gamma\_l(1-\\Gamma\\Gamma\_l)^{-1}\\mathbf{R} + (\\mathbf{S} - \\mathbf{I})\\}\\mathbf{a}$ |

Matched load ($\\Gamma\_l = 0$): $\\mathbf{b}' = (\\mathbf{S} - \\mathbf{I})\\mathbf{a}$ (2.78).

|  Eq.  |   |
| :---: | :-- |
| (2.79) | $\\vec{E}'(r,\\theta,\\phi) = \\frac{k}{\\sqrt{\\eta}}\\sum\_{i=1}^{J} b'\_i\\,\\vec{F}^{(3)}\_i(r,\\theta,\\phi)$ |

---

## 5. Reciprocity (§2.3.2)

### Generalized Lorentz form (unprimed and primed = adjoint device)

|  Eq.  |   |
| :---: | :-- |
| (2.93) | $vw' - v'w = \\sum\_{smn}(-1)^m\\bigl(b\_{smn}\\,a'\_{s,-m,n} - a\_{smn}\\,b'\_{s,-m,n}\\bigr)$ |

### Adjoint-antenna relations (with $\\mathbf{R}'$, $\\mathbf{T}$, etc.)

|  Eq.  |   |
| :---: | :-- |
| (2.104) | $R'\_{smn} = (-1)^m\\,T\_{s,-m,n}$ |
| (2.105) | $S'^{\\sigma\\mu\\nu}\_{smn} = (-1)^{m+\\mu}\\,S^{s,-m,n}\_{\\sigma,-\\mu,\\nu}$ |
| (2.106) | $\\Gamma' = \\Gamma$ |

### Reciprocal antenna (special case — antenna is its own adjoint)

|  Eq.  |   |
| :---: | :-- |
| (2.107) | $\\boxed{\\;R\_{smn} = (-1)^m\\,T\_{s,-m,n}\\;}$ |
| (2.108) | $S^{\\sigma\\mu\\nu}\_{smn} = (-1)^{m+\\mu}\\,S^{s,-m,n}\_{\\sigma,-\\mu,\\nu}$ |

---

## 6. Fields of Electric and Magnetic Dipoles (§2.3.3)

### z-directed electric dipole at origin

Dipole moment $d\_e = I\\ell$.

|  Eq.  |   |
| :---: | :-- |
| (2.109) | $\\vec{E}^{z}\_e = -\\frac{\\zeta k^2}{2\\pi}d\_e\\,\\frac{h\_1^{(1)}(kr)}{kr}\\cos\\theta\\,\\hat{r} + \\frac{\\zeta k^2}{4\\pi}d\_e\\,\\frac{1}{kr}\\frac{d}{d(kr)}\\{kr\\,h\_1^{(1)}(kr)\\}\\sin\\theta\\,\\hat{\\theta}$ |
| (2.110) | $\\vec{H}^{z}\_e = \\frac{ik^2}{4\\pi}d\_e\\,h\_1^{(1)}(kr)\\sin\\theta\\,\\hat{\\phi}$ |

with

|  Eq.  |   |
| :---: | :-- |
| (2.111) | $h\_1^{(1)}(kr) = -\\frac{e^{ikr}}{kr}\\Bigl(1 + \\frac{i}{kr}\\Bigr)$ |
| (2.112) | $\\frac{1}{kr}\\frac{d}{d(kr)}\\{kr\\,h\_1^{(1)}(kr)\\} = \\frac{e^{ikr}}{kr}\\Bigl\\{-i + \\frac{1}{kr} + \\frac{i}{(kr)^2}\\Bigr\\}$ |

### Spherical-wave form

|  Eq.  |   |
| :---: | :-- |
| (2.115) | $\\vec{E}^{z}\_e = \\frac{k}{\\sqrt{\\eta}}\\,Q\_{201}\\,\\vec{F}^{(3)}\_{201}(r,\\theta,\\phi)$ |
| (2.116) | $\\vec{H}^{z}\_e = -ik\\sqrt{\\eta}\\,Q\_{201}\\,\\vec{F}^{(3)}\_{101}(r,\\theta,\\phi)$ |
| (2.117) | $\\boxed{\\;Q\_{201} = -\\frac{1}{\\sqrt{6\\pi}}\\,\\frac{k}{\\sqrt{\\eta}}\\,d\_e\\;}$ |

### z-directed magnetic dipole

Moment $d\_m = I\_m \\ell = -i\\omega\\mu S I'$ (2.130, 2.131).

|  Eq.  |   |
| :---: | :-- |
| (2.132) | $\\vec{E}^{z}\_m = \\frac{k}{\\sqrt{\\eta}}\\,Q\_{101}\\,\\vec{F}^{(3)}\_{101}(r,\\theta,\\phi)$ |
| (2.133) | $\\vec{H}^{z}\_m = -ik\\sqrt{\\eta}\\,Q\_{101}\\,\\vec{F}^{(3)}\_{201}(r,\\theta,\\phi)$ |
| (2.138) | $\\boxed{\\;Q\_{101} = -\\frac{i}{\\sqrt{6\\pi}}\\,k\\sqrt{\\eta}\\,d\_m\\;}$ |

### Duality

If $d\_m = -\\zeta d\_e$ (2.139):

|  Eq.  |   |
| :---: | :-- |
| (2.140, 2.141) | $\\vec{E}\_m = \\zeta\\vec{H}\_e,\\quad \\vec{H}\_m = -\\eta\\vec{E}\_e$ |
| (2.142) | $Q\_{101} = -i\\,Q\_{201}$ |

### Rotation result for x-directed dipole (Euler angles $(0,-\\pi/2,0)$)

Rotation coefficients used:

|  Eq.  |   |
| :---: | :-- |
| (2.120-2.122) | $d^{1}\_{-1,0}(-\\pi/2) = \\tfrac{\\sqrt{2}}{2},\\quad d^{1}\_{0,0}(-\\pi/2) = 0,\\quad d^{1}\_{1,0}(-\\pi/2) = -\\tfrac{\\sqrt{2}}{2}$ |

Yielding:

|  Eq.  |   |
| :---: | :-- |
| (2.123) | $\\vec{F}^{(3)}\_{201}(r,\\theta,\\phi) = \\tfrac{\\sqrt{2}}{2}\\vec{F}^{(3)}\_{2,-1,1}(r',\\theta',\\phi') - \\tfrac{\\sqrt{2}}{2}\\vec{F}^{(3)}\_{211}(r',\\theta',\\phi')$ |

### Fields of x- and y-directed dipoles

**x-directed electric dipole:**

|  Eq.  |   |
| :---: | :-- |
| (2.124) | $\\vec{E}^{x}\_e = \\frac{k}{\\sqrt{\\eta}}Q\_{201}\\,\\tfrac{\\sqrt{2}}{2}\\{\\vec{F}^{(3)}\_{2,-1,1} - \\vec{F}^{(3)}\_{211}\\}$ |
| (2.125) | $\\vec{H}^{x}\_e = -ik\\sqrt{\\eta}\\,Q\_{201}\\,\\tfrac{\\sqrt{2}}{2}\\{\\vec{F}^{(3)}\_{1,-1,1} - \\vec{F}^{(3)}\_{111}\\}$ |

**y-directed electric dipole:**

|  Eq.  |   |
| :---: | :-- |
| (2.128) | $\\vec{E}^{y}\_e = \\frac{k}{\\sqrt{\\eta}}Q\_{201}\\,\\tfrac{i\\sqrt{2}}{2}\\{\\vec{F}^{(3)}\_{2,-1,1} + \\vec{F}^{(3)}\_{211}\\}$ |
| (2.129) | $\\vec{H}^{y}\_e = -ik\\sqrt{\\eta}\\,Q\_{201}\\,\\tfrac{i\\sqrt{2}}{2}\\{\\vec{F}^{(3)}\_{1,-1,1} + \\vec{F}^{(3)}\_{111}\\}$ |

**x-directed magnetic dipole:**

|  Eq.  |   |
| :---: | :-- |
| (2.134) | $\\vec{E}^{x}\_m = \\frac{k}{\\sqrt{\\eta}}Q\_{101}\\,\\tfrac{\\sqrt{2}}{2}\\{\\vec{F}^{(3)}\_{1,-1,1} - \\vec{F}^{(3)}\_{111}\\}$ |
| (2.135) | $\\vec{H}^{x}\_m = -ik\\sqrt{\\eta}\\,Q\_{101}\\,\\tfrac{\\sqrt{2}}{2}\\{\\vec{F}^{(3)}\_{2,-1,1} - \\vec{F}^{(3)}\_{211}\\}$ |

**y-directed magnetic dipole:**

|  Eq.  |   |
| :---: | :-- |
| (2.136) | $\\vec{E}^{y}\_m = \\frac{k}{\\sqrt{\\eta}}Q\_{101}\\,\\tfrac{i\\sqrt{2}}{2}\\{\\vec{F}^{(3)}\_{1,-1,1} + \\vec{F}^{(3)}\_{111}\\}$ |
| (2.137) | $\\vec{H}^{y}\_m = -ik\\sqrt{\\eta}\\,Q\_{101}\\,\\tfrac{i\\sqrt{2}}{2}\\{\\vec{F}^{(3)}\_{2,-1,1} + \\vec{F}^{(3)}\_{211}\\}$ |

---

## 7. Scattering Matrices for Electric and Magnetic Dipoles (§2.3.4)

> All dipoles assumed lossless and matched ($\\Gamma = 0$, $|\\mathbf{T}|^2 = 1$).
> Matrices indexed with the single-index convention (2.27). Layout: full $(J+1)\\times(J+1)$ block of $\\hat{\\mathbf{S}}$, columns/rows $0\\dots 7$ (where row/column 0 corresponds to the local-port reflection / receiving / transmitting; rows/cols 1..7 are the first seven spherical-mode ports = all $n=1$ modes plus the first $n=2$ mode).
>
> **Index reminder for $n=1$:** $j=1: (1,-1,1)$, $j=2: (2,-1,1)$, $j=3: (1,0,1)$, $j=4: (2,0,1)$, $j=5: (1,1,1)$, $j=6: (2,1,1)$, $j=7: (1,-2,2)$.

### z-directed electric dipole — Eq. (2.148)

Only mode $j=4$ ($s=2, m=0, n=1$) couples to the port.
- $T\_4 = 1$, all other $T\_i = 0$.
- $R\_4 = 1$, all other $R\_j = 0$ (consistent with reciprocity (2.107)).
- $S\_{ij} = \\delta\_{ij}$ for $i,j \\ne 4$; $S\_{4j} = S\_{i4} = 0$.

```math
\hat{\mathbf{S}}^{z}_e = \begin{bmatrix} 0 & 0 & 0 & 0 & 1 & 0 & 0 & \cdots\\ 0 & 1 & 0 & 0 & 0 & 0 & 0\\ 0 & 0 & 1 & 0 & 0 & 0 & 0\\ 0 & 0 & 0 & 1 & 0 & 0 & 0\\ 1 & 0 & 0 & 0 & 0 & 0 & 0\\ 0 & 0 & 0 & 0 & 0 & 1 & 0\\ 0 & 0 & 0 & 0 & 0 & 0 & 1\\ \vdots & & & & & & & \ddots \end{bmatrix}
```

### x-directed electric dipole — Eq. (2.154)

Non-zero $\\mathbf{T}$: $T\_2 = \\sqrt{2}/2$, $T\_6 = -\\sqrt{2}/2$.
Non-zero $\\mathbf{R}$: $R\_2 = \\sqrt{2}/2$, $R\_6 = -\\sqrt{2}/2$.
Non-zero $\\mathbf{S}$ entries (within $n=1$ block): $S\_{22} = S\_{26} = S\_{62} = S\_{66} = 1/2$.
All other entries: $S\_{ii} = 1$ outside the $\\{2,6\\}$ subspace; rest zero.

```math
\hat{\mathbf{S}}^{x}_e = \begin{bmatrix} 0 & 0 & \tfrac{\sqrt 2}{2} & 0 & 0 & 0 & -\tfrac{\sqrt 2}{2} & 0 & \cdots\\ 0 & 1 & 0 & 0 & 0 & 0 & 0 & 0\\ \tfrac{\sqrt 2}{2} & 0 & \tfrac{1}{2} & 0 & 0 & 0 & \tfrac{1}{2} & 0\\ 0 & 0 & 0 & 1 & 0 & 0 & 0 & 0\\ 0 & 0 & 0 & 0 & 1 & 0 & 0 & 0\\ 0 & 0 & 0 & 0 & 0 & 1 & 0 & 0\\ -\tfrac{\sqrt 2}{2} & 0 & \tfrac{1}{2} & 0 & 0 & 0 & \tfrac{1}{2} & 0\\ 0 & 0 & 0 & 0 & 0 & 0 & 0 & 1\\ \vdots & & & & & & & & \ddots \end{bmatrix}
```

### y-directed electric dipole — Eq. (2.155)

Non-zero $\\mathbf{T}$: $T\_2 = i\\sqrt{2}/2$, $T\_6 = i\\sqrt{2}/2$.
Non-zero $\\mathbf{R}$: $R\_2 = -i\\sqrt{2}/2$, $R\_6 = -i\\sqrt{2}/2$.
Non-zero $\\mathbf{S}$ in $\\{2,6\\}$ block: $S\_{22} = 1/2$, $S\_{26} = -1/2$, $S\_{62} = -1/2$, $S\_{66} = 1/2$.

```math
\hat{\mathbf{S}}^{y}_e = \begin{bmatrix} 0 & 0 & -\tfrac{i\sqrt 2}{2} & 0 & 0 & 0 & -\tfrac{i\sqrt 2}{2} & 0 & \cdots\\ 0 & 1 & 0 & 0 & 0 & 0 & 0 & 0\\ \tfrac{i\sqrt 2}{2} & 0 & \tfrac{1}{2} & 0 & 0 & 0 & -\tfrac{1}{2} & 0\\ 0 & 0 & 0 & 1 & 0 & 0 & 0 & 0\\ 0 & 0 & 0 & 0 & 1 & 0 & 0 & 0\\ 0 & 0 & 0 & 0 & 0 & 1 & 0 & 0\\ \tfrac{i\sqrt 2}{2} & 0 & -\tfrac{1}{2} & 0 & 0 & 0 & \tfrac{1}{2} & 0\\ 0 & 0 & 0 & 0 & 0 & 0 & 0 & 1\\ \vdots & & & & & & & & \ddots \end{bmatrix}
```

### z-directed magnetic dipole — Eq. (2.156)

Only mode $j=3$ ($s=1, m=0, n=1$) couples. $T\_3 = -i$, $R\_3 = -i$, $S\_{33} = 0$.

```math
\hat{\mathbf{S}}^{z}_m = \begin{bmatrix} 0 & 0 & 0 & -i & 0 & 0 & \cdots\\ 0 & 1 & 0 & 0 & 0 & 0\\ 0 & 0 & 1 & 0 & 0 & 0\\ -i & 0 & 0 & 0 & 0 & 0\\ 0 & 0 & 0 & 0 & 1 & 0\\ 0 & 0 & 0 & 0 & 0 & 1\\ \vdots & & & & & & \ddots \end{bmatrix}
```

### x-directed magnetic dipole — Eq. (2.157)

Non-zero $\\mathbf{T}$: $T\_1 = -i\\sqrt{2}/2$, $T\_5 = i\\sqrt{2}/2$.
$\\{1,5\\}$ block: $S\_{11} = S\_{55} = 1/2$, $S\_{15} = S\_{51} = 1/2$.

```math
\hat{\mathbf{S}}^{x}_m = \begin{bmatrix} 0 & -\tfrac{i\sqrt 2}{2} & 0 & 0 & 0 & \tfrac{i\sqrt 2}{2} & 0 & \cdots\\ -\tfrac{i\sqrt 2}{2} & \tfrac{1}{2} & 0 & 0 & 0 & \tfrac{1}{2} & 0\\ 0 & 0 & 1 & 0 & 0 & 0 & 0\\ 0 & 0 & 0 & 1 & 0 & 0 & 0\\ 0 & 0 & 0 & 0 & 1 & 0 & 0\\ \tfrac{i\sqrt 2}{2} & \tfrac{1}{2} & 0 & 0 & 0 & \tfrac{1}{2} & 0\\ 0 & 0 & 0 & 0 & 0 & 0 & 1\\ \vdots & & & & & & & \ddots \end{bmatrix}
```

### y-directed magnetic dipole — Eq. (2.158)

Non-zero $\\mathbf{T}$: $T\_1 = \\sqrt{2}/2$, $T\_5 = \\sqrt{2}/2$.
$\\{1,5\\}$ block: $S\_{11} = 1/2$, $S\_{15} = -1/2$, $S\_{51} = -1/2$, $S\_{55} = 1/2$.

```math
\hat{\mathbf{S}}^{y}_m = \begin{bmatrix} 0 & -\tfrac{\sqrt 2}{2} & 0 & 0 & 0 & -\tfrac{\sqrt 2}{2} & 0 & \cdots\\ \tfrac{\sqrt 2}{2} & \tfrac{1}{2} & 0 & 0 & 0 & -\tfrac{1}{2} & 0\\ 0 & 0 & 1 & 0 & 0 & 0 & 0\\ 0 & 0 & 0 & 1 & 0 & 0 & 0\\ 0 & 0 & 0 & 0 & 1 & 0 & 0\\ \tfrac{\sqrt 2}{2} & -\tfrac{1}{2} & 0 & 0 & 0 & \tfrac{1}{2} & 0\\ 0 & 0 & 0 & 0 & 0 & 0 & 1\\ \vdots & & & & & & & \ddots \end{bmatrix}
```

### Combined elements (Huygens, turnstile)

**z-directed Huygens source** = $\\hat{x}$-electric + $\\hat{y}$-magnetic-dual dipoles (the magnetic dipole is the dual source of the electric one). Total matrix (Eq. 2.159):

```math
\hat{\mathbf{S}}^{z}_H = \begin{bmatrix} 0 & -\tfrac{1}{2} & \tfrac{1}{2} & 0 & 0 & -\tfrac{1}{2} & -\tfrac{1}{2} & 0 & \cdots\\ \tfrac{1}{2} & \tfrac{3}{4} & \tfrac{1}{4} & 0 & 0 & -\tfrac{1}{4} & -\tfrac{1}{4} & 0\\ \tfrac{1}{2} & -\tfrac{1}{4} & \tfrac{1}{4} & 0 & 0 & -\tfrac{1}{4} & \tfrac{3}{4} & 0\\ 0 & 0 & 0 & 1 & 0 & 0 & 0 & 0\\ 0 & 0 & 0 & 0 & 1 & 0 & 0 & 0\\ \tfrac{1}{2} & -\tfrac{1}{4} & \tfrac{1}{4} & 0 & 0 & \tfrac{3}{4} & -\tfrac{1}{4} & 0\\ -\tfrac{1}{2} & \tfrac{1}{4} & \tfrac{3}{4} & 0 & 0 & \tfrac{1}{4} & \tfrac{1}{4} & 0\\ 0 & 0 & 0 & 0 & 0 & 0 & 0 & 1\\ \vdots & & & & & & & & \ddots \end{bmatrix}
```

**z-directed turnstile** = $\\hat{x}$-electric + $\\hat{y}$-electric dipoles in phase quadrature. Total matrix (Eq. 2.160):

```math
\hat{\mathbf{S}}^{z}_T = \begin{bmatrix} 0 & 0 & 1 & 0 & 0 & 0 & 0 & 0 & \cdots\\ 0 & 1 & 0 & 0 & 0 & 0 & 0 & 0\\ 0 & 0 & 0 & 0 & 0 & 0 & 1 & 0\\ 0 & 0 & 0 & 1 & 0 & 0 & 0 & 0\\ 0 & 0 & 0 & 0 & 1 & 0 & 0 & 0\\ 0 & 0 & 0 & 0 & 0 & 1 & 0 & 0\\ -1 & 0 & 0 & 0 & 0 & 0 & 0 & 0\\ 0 & 0 & 0 & 0 & 0 & 0 & 0 & 1\\ \vdots & & & & & & & & \ddots \end{bmatrix}
```

Only $T\_6$ and $R\_2$ are non-zero; corresponds to far field $\\vec{F}\_6 = \\vec{F}^{(3)}\_{211}$ only (RHCP at $\\theta = 0$, LHCP at $\\theta = \\pi$).

### Composition rule

Receiving/transmitting elements of an antenna built from multiple dipoles = sum of the individual elements weighted by excitation. Scattering elements **cannot** be simply summed — must be derived from incident-field considerations.

### Receiving-formula special case (z-electric dipole, matched load, $\\Gamma\_l = 0$)

|  Eq.  |   |
| :---: | :-- |
| (2.152) | $w = a\_4 = \\frac{\\sqrt{6\\pi}}{2}\\,\\frac{\\sqrt{\\eta}}{k}\\,E\_z(0,\\theta,\\phi)$ |

Scattered field for that dipole: $b\_i = a\_i$ for $i \\ne 4$, $b\_4 = 0$ (Eq. 2.153).

---

## 8. Alternative ("Source") Scattering Matrix — Yaghjian (§2.3.5)

Same waveguide modes, but use $\\vec{F}^{(1)}\_j$ (standing) instead of $\\vec{F}^{(4)}\_j$ for the second basis:

|  Eq.  |   |
| :---: | :-- |
| (2.164) | $\\vec{E} = \\frac{k}{\\sqrt{\\eta}}\\sum\_j \\bigl\\{a'\_j\\,\\vec{F}^{(1)}\_j + b'\_j\\,\\vec{F}^{(3)}\_j\\bigr\\}$ |

Using $\\vec{F}^{(1)}\_j = \\tfrac{1}{2}\\{\\vec{F}^{(3)}\_j + \\vec{F}^{(4)}\_j\\}$ (2.165):

**(2.166)**

```math
\begin{bmatrix}\Gamma' & \mathbf{R}'\\\mathbf{T}' & \mathbf{S}'\end{bmatrix} = \begin{bmatrix}\Gamma & \tfrac{1}{2}\mathbf{R}\\\mathbf{T} & \tfrac{1}{2}(\mathbf{S} - \mathbf{I})\end{bmatrix}
```

Implications:
- $\\Gamma' = \\Gamma$, $\\mathbf{T}' = \\mathbf{T}$.
- $\\mathbf{R}' = \\tfrac{1}{2}\\mathbf{R}$.
- Empty space: $\\mathbf{S}' = \\mathbf{0}$ (vs. $\\mathbf{S} = \\mathbf{I}$ in the classical form).

For pure scatterers (no local port), $\\mathbf{S}'$ reduces to Waterman's T-matrix [25]. The book retains the classical formulation throughout.

---

## 9. Far-Field Patterns (§2.4.1)

### Definition

|  Eq.  |   |
| :---: | :-- |
| (2.175) | $\\vec{K}\_{smn}(\\theta,\\phi) = \\lim\_{kr\\to\\infty}\\Bigl[\\sqrt{4\\pi}\\,\\frac{kr}{e^{ikr}}\\,\\vec{F}^{(3)}\_{smn}(r,\\theta,\\phi)\\Bigr]$ |

### Explicit pattern functions

|  Eq.  |   |
| :---: | :-- |
| (2.176) | $\\boxed{\\;\\vec{K}\_{1mn}(\\theta,\\phi) = \\sqrt{\\frac{2}{n(n+1)}}\\Bigl(-\\frac{m}{\|m\|}\\Bigr)^{m} e^{im\\phi}(-i)^{n+1}\\left\\{\\frac{im\\,\\bar{P}\_n^{\|m\|}(\\cos\\theta)}{\\sin\\theta}\\,\\hat{\\theta} - \\frac{d\\bar{P}\_n^{\|m\|}(\\cos\\theta)}{d\\theta}\\,\\hat{\\phi}\\right\\}\\;}$ |
| (2.177) | $\\boxed{\\;\\vec{K}\_{2mn}(\\theta,\\phi) = \\sqrt{\\frac{2}{n(n+1)}}\\Bigl(-\\frac{m}{\|m\|}\\Bigr)^{m} e^{im\\phi}(-i)^{n}\\left\\{\\frac{d\\bar{P}\_n^{\|m\|}(\\cos\\theta)}{d\\theta}\\,\\hat{\\theta} + \\frac{im\\,\\bar{P}\_n^{\|m\|}(\\cos\\theta)}{\\sin\\theta}\\,\\hat{\\phi}\\right\\}\\;}$ |

### Useful identity

|  Eq.  |   |
| :---: | :-- |
| (2.178) | $\\vec{K}\_{smn} = i\\,\\hat{r}\\times\\vec{K}\_{3-s,m,n}$ |

### Far-field expressions

|  Eq.  |   |
| :---: | :-- |
| (2.179,\,2.180) | $\\vec{E}(r,\\theta,\\phi) \\to \\frac{k}{\\sqrt{\\eta}}\\,\\frac{1}{\\sqrt{4\\pi}}\\,\\frac{e^{ikr}}{kr}\\sum\_{smn} Q^{(3)}\_{smn}\\,\\vec{K}\_{smn}(\\theta,\\phi) = \\frac{k}{\\sqrt{\\eta}}\\,\\frac{1}{\\sqrt{4\\pi}}\\,\\frac{e^{ikr}}{kr}\\,v\\sum\_{smn} T\_{smn}\\,\\vec{K}\_{smn}(\\theta,\\phi) = \\frac{k}{\\sqrt{\\eta}}\\,\\frac{1}{\\sqrt{4\\pi}}\\,\\frac{e^{ikr}}{kr}\\,v\\,\\vec{K}(\\theta,\\phi)$ |
| (2.181) | $\\vec{H} \\to \\eta\\,\\hat{r}\\times\\vec{E} = k\\sqrt{\\eta}\\,\\frac{1}{\\sqrt{4\\pi}}\\,\\frac{e^{ikr}}{kr}\\,v\\,\\hat{r}\\times\\vec{K}(\\theta,\\phi)$ |

### Absolute far-field pattern

|  Eq.  |   |
| :---: | :-- |
| (2.182) | $\\vec{K}(\\theta,\\phi) = \\sum\_{smn} T\_{smn}\\,\\vec{K}\_{smn}(\\theta,\\phi)$ |

Dimensionless. $C\\vec{K}(\\theta,\\phi)$ for arbitrary $C$ is a **relative far-field pattern**.

---

## 10. Polarization (§2.4.2)

### Orthogonal polarization unit vector

|  Eq.  |   |
| :---: | :-- |
| (2.183) | $\\hat{i}\_{\\text{cross}}(\\theta,\\phi) = \\hat{r}\\times\\hat{i}\_{co}^{\*}(\\theta,\\phi)$ |

with $\\hat{i}\_{co}\\cdot\\hat{i}\_{\\text{cross}}^{\*} = 0$ (2.184).

### Decomposition

|  Eq.  |   |
| :---: | :-- |
| (2.185, 2.186) | $K\_{co}(\\theta,\\phi) = \\vec{K}\\cdot\\hat{i}\_{co}^{\*} = \\sum\_{smn} T\_{smn}\\,\\vec{K}\_{smn}\\cdot\\hat{i}\_{co}^{\*}$ |
| (2.187, 2.188) | $K\_{\\text{cross}}(\\theta,\\phi) = \\vec{K}\\cdot\\hat{i}\_{\\text{cross}}^{\*}$ |
| (2.189) | $\\vec{K} = K\_{co}\\,\\hat{i}\_{co} + K\_{\\text{cross}}\\,\\hat{i}\_{\\text{cross}}$ |

### Ludwig's "Definition 3" — Linear (boresight along $+z$, reference angle $\\phi\_o$)

|  Eq.  |   |
| :---: | :-- |
| (2.190) | $\\hat{i}\_{co,3L}(\\theta,\\phi) = \\hat{\\theta}\\cos(\\phi-\\phi\_o) - \\hat{\\phi}\\sin(\\phi-\\phi\_o),\\quad 0\\le\\theta<\\pi$ |
| (2.191) | $\\hat{i}\_{\\text{cross},3L}(\\theta,\\phi) = \\hat{\\theta}\\sin(\\phi-\\phi\_o) + \\hat{\\phi}\\cos(\\phi-\\phi\_o)$ |

With $\\phi\_o = 0$: identical to a $\\hat{x}$-electric + $\\hat{y}$-magnetic Huygens source.

### Circular polarization unit vectors ($\\phi\_o = 0$)

**Right-hand (RCP):**

|  Eq.  |   |
| :---: | :-- |
| (2.192) | $\\hat{i}\_{co,RC}(\\theta,\\phi) = \\tfrac{1}{\\sqrt 2}\\bigl[\\hat{i}\_{co,3L} + i\\,\\hat{i}\_{\\text{cross},3L}\\bigr]\_{\\phi\_o=0} = \\tfrac{1}{\\sqrt 2}\\,e^{i\\phi}(\\hat{\\theta} + i\\hat{\\phi})$ |
| (2.193) | $\\hat{i}\_{\\text{cross},RC}(\\theta,\\phi) = \\tfrac{i}{\\sqrt 2}\\,e^{-i\\phi}(\\hat{\\theta} - i\\hat{\\phi})$ |

**Left-hand (LCP):**

|  Eq.  |   |
| :---: | :-- |
| (2.194) | $\\hat{i}\_{co,LC}(\\theta,\\phi) = \\tfrac{1}{\\sqrt 2}\\,e^{-i\\phi}(\\hat{\\theta} - i\\hat{\\phi})$ |
| (2.195) | $\\hat{i}\_{\\text{cross},LC}(\\theta,\\phi) = \\tfrac{-i}{\\sqrt 2}\\,e^{i\\phi}(\\hat{\\theta} + i\\hat{\\phi})$ |

The factors $i$ may be dropped in practice but the $e^{\\pm i\\phi}$ factors must remain for continuity at $\\theta = 0$. All unit-vector distributions are discontinuous at $\\theta = \\pi$.

### Polarization-ellipse parameters

Decomposition (without the cross-polar $i$):

|  Eq.  |   |
| :---: | :-- |
| (2.196) | $\\vec{K}(\\theta,\\phi) = K\_R\\,\\hat{i}\_{co,RC} + K\_L\\,\\hat{i}\_{co,LC}$ |
| (2.197, 2.198) | $K\_R = \|K\_R\|\\,e^{i\\psi\_R},\\quad K\_L = \|K\_L\|\\,e^{i\\psi\_L}$ |
| (2.199) | $Q = K\_R / K\_L$ |

**Axial ratio** $r = \\tan\\alpha$:

|  Eq.  |   |
| :---: | :-- |
| (2.200, 2.201) | $\\tan\\alpha = \\frac{\|K\_R\| - \|K\_L\|}{\|K\_R\| + \|K\_L\|} = \\frac{\|Q\| - 1}{\|Q\| + 1},\\quad -\\tfrac{\\pi}{4}\\le\\alpha\\le\\tfrac{\\pi}{4}$ |

$\\alpha = 0$ ↔ linear; $\\alpha > 0$ ↔ RH-elliptical; $\\alpha < 0$ ↔ LH-elliptical.

**Tilt angle** $\\beta$ (relative to $\\hat{i}\_{co,3L}|\_{\\phi\_o=0}$):

|  Eq.  |   |
| :---: | :-- |
| (2.202, 2.203) | $\\beta = \\frac{\\psi\_L - \\psi\_R}{2} = -\\tfrac{1}{2}\\arg(Q)$ |

---

## 11. Directivity and Gain (§2.4.3)

### Power per unit solid angle (far field)

|  Eq.  |   |
| :---: | :-- |
| (2.204, 2.205) | $r^2\\,\\tfrac{1}{2}\\,\\mathrm{Re}\\{\\vec{E}\\times\\vec{H}^\*\\}\\cdot\\hat{r} = r^2\\,\\tfrac{1}{2}\\eta\|\\vec{E}\|^2 = \\frac{1}{2}\\,\\frac{1}{4\\pi}\\Bigl\|\\sum\_{smn}Q^{(3)}\_{smn}\\vec{K}\_{smn}(\\theta,\\phi)\\Bigr\|^2$ |

### Isotropic reference

|  Eq.  |   |
| :---: | :-- |
| (2.206) | $\\frac{P}{4\\pi} = \\frac{1}{4\\pi}\\,\\frac{1}{2}\\sum\_{smn}\|Q^{(3)}\_{smn}\|^2$ |

### Directivity

|  Eq.  |   |
| :---: | :-- |
| (2.207) | $\\boxed{\\;D(\\theta,\\phi) = \\frac{\\Bigl\|\\sum\_{smn} Q^{(3)}\_{smn}\\,\\vec{K}\_{smn}(\\theta,\\phi)\\Bigr\|^2}{\\sum\_{smn}\|Q^{(3)}\_{smn}\|^2}\\;}$ |

Equivalently (since $Q^{(3)}\_{smn} = b\_{smn} = v\\,T\_{smn}$):

|  Eq.  |   |
| :---: | :-- |
| (2.208) | $D(\\theta,\\phi) = \\frac{\\Bigl\|\\sum\_{smn} T\_{smn}\\,\\vec{K}\_{smn}(\\theta,\\phi)\\Bigr\|^2}{\\sum\_{smn}\|T\_{smn}\|^2} = \\frac{\|\\vec{K}(\\theta,\\phi)\|^2}{\\sum\_{smn}\|T\_{smn}\|^2} = \|\\vec{K}(\\theta,\\phi)\|^2$ |

The last equality requires **matched and lossless** ($\\sum |T\_{smn}|^2 = 1$).

### Gain

Input power $P\_{\\text{in}} = \\tfrac{1}{2}|v|^2(1 - |\\Gamma|^2)$ (2.209).

|  Eq.  |   |
| :---: | :-- |
| (2.210, 2.211, 2.212) | $\\boxed{\\;G(\\theta,\\phi) = \\frac{\\Bigl\|\\sum\_{smn} Q^{(3)}\_{smn}\\,\\vec{K}\_{smn}(\\theta,\\phi)\\Bigr\|^2}{\|v\|^2(1-\|\\Gamma\|^2)} = \\frac{\|\\vec{K}(\\theta,\\phi)\|^2}{1-\|\\Gamma\|^2} = \\frac{\|\\vec{K}(\\theta,\\phi)\|^2}{P\_{\\text{loss}}/P\_{\\text{inc}} + \\sum\|T\_{smn}\|^2}\\;}$ |

### Polarization additivity

$D$ (or $G$) in any direction = sum of $D$ (or $G$) for any two orthogonal polarizations.

---

## 12. Maximum Directivity (§2.4.4)

### Component directivities

|  Eq.  |   |
| :---: | :-- |
| (2.213) | $D(\\theta,\\phi) = D\_{co}(\\theta,\\phi) + D\_{\\text{cross}}(\\theta,\\phi)$ |
| (2.214) | $D\_{co}(\\theta,\\phi) = \\frac{\|\\sum\_{smn} Q\_{smn}\\,\\vec{K}\_{smn}(\\theta,\\phi)\\cdot\\hat{i}\_{co}^{\*}\|^2}{\\sum\_{smn}\|Q\_{smn}\|^2} = \|\\vec{K}(\\theta,\\phi)\\cdot\\hat{i}\_{co}^{\*}\|^2$ |
| (2.215) | $D\_{\\text{cross}}(\\theta,\\phi) = \|\\vec{K}(\\theta,\\phi)\\cdot\\hat{i}\_{\\text{cross}}^{\*}\|^2$ |

### Cauchy–Schwartz bound

|  Eq.  |   |
| :---: | :-- |
| (2.216, 2.217) | $D\_{co}(\\theta',\\phi') \\le \\sum\_{smn}\|\\vec{K}\_{smn}(\\theta',\\phi')\\cdot\\hat{i}\_{co}^{\*}\|^2 = D\_{co,\\max}(\\theta',\\phi')$ |

Equality (i.e. the maximum) is achieved by:

|  Eq.  |   |
| :---: | :-- |
| (2.218) | $\\boxed{\\;Q\_{smn} = c\\,\\bigl(\\vec{K}\_{smn}(\\theta',\\phi')\\cdot\\hat{i}\_{co}^{\*}\\bigr)^\*\\;}$ |

with arbitrary constant $c$.

### Maximum directivity value (independent of direction and polarization)

For $\\hat{i}\_{co} = \\alpha\\hat{\\theta} + \\beta\\hat{\\phi}$, $|\\alpha|^2 + |\\beta|^2 = 1$ (2.219, 2.220):

|  Eq.  |   |
| :---: | :-- |
| (2.225, 2.226) | $\\boxed{\\;D\_{co,\\max}(\\theta',\\phi') = N^2 + 2N\\;}$ |

where $N$ is the truncation in $n$.

### Coefficients for $\\hat{x}$-polarized peak at $(\\theta',\\phi')=(0,0)$

|  Eq.  |   |
| :---: | :-- |
| (2.227) | $Q\_{1,1,n} = Q\_{1,-1,n} = Q\_{2,1,n} = -Q\_{2,-1,n} = c\\,(-i^{n})\\,\\tfrac{1}{2}\\sqrt{2n+1}$ |

All other coefficients zero. Cross-polarization is zero in every direction.

> These are proportional to the coefficients for an $\\hat{x}$-polarized plane wave travelling along $+z$, but evaluated with $\\vec{K}\_{smn}$ rather than $\\vec{F}^{(1)}\_{smn}$.

### Reference values

| $N$ | $D\_{\\max} = N^2+2N$ | $D\_{\\max}$ (dB) |
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

1. **Time convention.** Code must use $e^{-i\\omega t}$. Outgoing modes carry $h\_n^{(1)} \\sim e^{+ikr}/(kr)$. *Validate with `(2.13)`–`(2.16)`.*
2. **Hankel-function identity.** `j_n(z) == 0.5 * (h_n^(1)(z) + h_n^(2)(z))` — Eq. (2.30).
3. **Wronskian.** `R1n^(1) * R2n^(2) - R2n^(1) * R1n^(2) == 1/(kr)^2` — Eq. (2.52).
4. **$\\bar P\_n^{|m|}$ normalization.** Use Belousov's convention (not Schmidt, not unnormalized). Spot-check $\\int\_{-1}^{1}\\bar P\_n^m(\\mu)^2 d\\mu = 2/(2n+1)$ × normalization factor consistent with Eq. (2.18).
5. **Phase factor.** `(-m/abs(m))**m == 1` when `m == 0` — Eq. (2.19). For $m \\ne 0$: equals $(-1)^m$ if $m > 0$, $(-1)^{|m|}\\cdot(\\text{sign})$ if $m < 0$. Reduces to $1$ for $m=0$ and to $\\pm 1$ otherwise per the Edmonds convention.
6. **Index swap.** $\\vec{H}$ uses $\\vec{F}^{(c)}\_{3-s,m,n}$ — Eq. (2.23). A common bug is mirroring $\\vec{E}$ in $\\vec{H}$.
7. **Single index.** `j = 2*(n*(n+1) + m - 1) + s` — Eq. (2.27). Spot-check the $n=1$ table above.
8. **Power normalization.** Pure $c=3$ field with one mode of unit $|Q|$ radiates $0.5$ W — Eq. (2.24).
9. **Truncation.** `N = floor(k*r0) + n1`, default `n1 = 10` — Eq. (2.31).
10. **Reciprocity.** For reciprocal antennas: `R[s,m,n] == (-1)**m * T[s,-m,n]` — Eq. (2.107).
11. **Unitarity (lossless).** `S_hat^H @ S_hat == I` and `|Γ|^2 + sum |T|^2 == 1` — Eq. (2.61).
12. **Empty space.** `S == I` (classical) but `S' == 0` (source-form) — §2.3.1 / Eq. (2.166).
13. **z-electric dipole.** Only $T\_4 = 1$, $R\_4 = 1$ — Eq. (2.148). Sanity check via `Q_201 = -k*d_e/(sqrt(6*pi)*sqrt(eta))` — Eq. (2.117).
14. **z-magnetic dipole.** $T\_3 = -i$, $R\_3 = -i$ — Eq. (2.156). $Q\_{101} = -i\\,Q\_{201}$ when $d\_m = -\\zeta d\_e$ — Eq. (2.142).
15. **Far-field K identity.** `K_smn = i * r_hat × K_{3-s,m,n}` — Eq. (2.178).
16. **Directivity self-consistency.** When the antenna is matched and lossless: `D(θ,φ) == |K(θ,φ)|^2` — Eq. (2.208).
17. **Max directivity.** A spherical-wave field truncated to $N$ has $D\_{\\max} = N^2 + 2N$ — Eq. (2.226). Test by maximizing the Cauchy–Schwartz inequality with coefficients (2.227).
18. **Adjoint vs reciprocal.** For non-reciprocal antennas, the $(-1)^m T\_{s,-m,n}$ relation gives the **adjoint** $R'$, not $R$ — Eqs. (2.104) vs. (2.107).
19. **Scattered field.** Use $\\mathbf{b}' = (\\mathbf{S}-\\mathbf{I})\\mathbf{a}$ under matched load — Eq. (2.78). Don't forget to subtract $\\mathbf{I}$.
20. **Axial ratio sign.** `α > 0` for right-handed, `α < 0` for left-handed, under the $e^{-i\\omega t}$ time convention — §2.4.2.

---

*End of Chapter 2 reference.*
