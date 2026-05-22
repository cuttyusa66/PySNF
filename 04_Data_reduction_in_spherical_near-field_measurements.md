# Chapter 4 — Data Reduction in Spherical Near-Field Measurements

**Source:** Hansen, *Spherical Near-Field Antenna Measurements*, Chapter 4 (pp. 89–146).

**Cross-references:**
- Symbols: [A5_List_of_principal_symbols_and_uses.md](A5_List_of_principal_symbols_and_uses.md)
- Chapter 2 (Scattering matrix, $\\vec{F}^{(c)}\_{smn}$): [02_Scattering_matrix_description_of_an_antenna.md](02_Scattering_matrix_description_of_an_antenna.md)
- Chapter 3 (Transmission formula): [03_Scattering_matrix_description_of_antenna_coupling.md](03_Scattering_matrix_description_of_antenna_coupling.md)
- Appendix A1 (wave functions, orthogonality integrals): [A1_Spherical_wave_functions_notation_and_properties.md](A1_Spherical_wave_functions_notation_and_properties.md)
- Appendix A2 (rotation, $d^n\_{\\mu m}$, $\\Delta^n\_{m'm}$): [A2_Rotation_of_spherical_waves.md](A2_Rotation_of_spherical_waves.md)
- Appendix A3 (translation $C^{sn(c)}\_{\\sigma\\mu\\nu}$): [A3_Translation_of_spherical_waves.md](A3_Translation_of_spherical_waves.md)
- Appendix A4 (DFT, sampling, reconstruction): [A4_Data_processing_in_antenna_measurements.md](A4_Data_processing_in_antenna_measurements.md)

> **Purpose.** Validation reference for **the inversion algorithm** (the SNIFT-style transformation): given measured probe signals $w(A, \\chi, \\theta, \\phi)$, recover the test-antenna transmitting coefficients $T\_{smn}$, then evaluate fields at any sphere via an output probe. Every formula tagged with its book number. The §4.4.5 worked example provides explicit intermediate values for unit tests.

---

## 1. Problem Statement and Big Picture (§4.1)

An antenna's total scattering matrix $\\hat{\\mathbf{S}}$ contains four blocks:

| Block | Count | Status |
|---|---|---|
| Reflection $\\Gamma$ | 1 | Direct measurement (trivial) |
| Scattering $S\_{ij}$ | $J^2$ | Hard; not pursued here |
| Receiving $R\_j$ | $J$ | Recovered if reciprocity used |
| Transmitting $T\_j$ | $J$ | **Main goal** of this chapter |

where $J = 2N(N+2)$ and $N = \\lceil kr\_0 \\rceil + \\sim 10$ (Eq. 2.31).

### Two approaches

| Approach | Probe model | When valid |
|---|---|---|
| **With probe correction** (§4.3) | Realistic horn; full $R^p\_{\\sigma\\mu\\nu}$ | All cases; the standard method |
| **Without probe correction** (§4.2) | Ideal short electric/magnetic dipoles | Low-directivity probes at moderate distance |

The book's pragmatic recommendation: **always use probe correction**. §4.2 is included for academic completeness and historical perspective.

---

## 2. Measurement Without Probe Correction (§4.2)

These methods assume the probe is an ideal dipole — signal is **proportional** to a field component. Then the unknown $T\_{smn}$ are extracted by orthogonality integrals from the **field** $\\vec{E}$ or $\\vec{H}$ sampled on the measurement sphere.

### 2.1 Setup

|  Eq.  |   |
| :---: | :-- |
| (4.2) | $\\vec{E}(A, \\theta, \\phi) = \\frac{k}{\\sqrt\\eta}\\sum\_{s=1}^2\\sum\_{n=1}^N\\sum\_{m=-n}^n v\\,T\_{smn}\\,\\vec{F}^{(3)}\_{smn}(A, \\theta, \\phi),\\quad A > r\_0$ |
| (4.3) | $\\vec{H}(A, \\theta, \\phi) = -ik\\sqrt\\eta\\sum v\\,T\_{smn}\\,\\vec{F}^{(3)}\_{3-s,m,n}(A, \\theta, \\phi)$ |

### 2.2 What you can recover from full-sphere $\\vec{E}$ or $\\vec{H}$ measurements (§4.2.2, Table 4.1)

Using orthogonality (A1.70) applied to the appropriate region/wave-type:

| Region (sources) | Measured | Recoverable |
|---|---|---|
| Both inner and outer have sources; both $Q^{(3)}, Q^{(4)}$ in expansion | $\\vec{E}(A)$ | **All TM** ($Q^{(3)}\_{2mn}$ + $Q^{(4)}\_{2mn}$) |
| Same | $\\vec{H}(A)$ | **All TE** ($Q^{(3)}\_{1mn}$ + $Q^{(4)}\_{1mn}$) |
| Inner is source-free | $\\vec{E}(A)$ | All $Q^{(1)}\_{smn}$ except possibly a single TE (if $j\_n(kA) = 0$) |
| Same | $\\vec{H}(A)$ | All $Q^{(1)}\_{smn}$ except possibly a single TM |
| Outer is source-free (typical antenna) | $\\vec{E}(A)$ | **All $Q^{(3)}\_{smn}$** |
| Same | $\\vec{H}(A)$ | **All $Q^{(3)}\_{smn}$** |

### 2.3 Radial-component measurement (§4.2.3) — uses orthogonality (A1.68)

Radial electric dipole probe:

|  Eq.  |   |
| :---: | :-- |
| (4.24) | $\\boxed{\\;v\\,T\_{2mn} = \\frac{2}{\\sqrt{6\\pi}}\\,(-1)^m\\,\\frac{1}{n(n+1)}\\,\\Bigl(\\frac{kA}{h\_n^{(1)}(kA)}\\Bigr)^2\\int\_0^{2\\pi}\\!\\!\\int\_0^\\pi w\_r^e(A, 0, \\theta, \\phi)\\,\\vec{F}^{(3)}\_{2,-m,n}(A,\\theta,\\phi)\\cdot\\hat r\\,\\sin\\theta\\,d\\theta\\,d\\phi\\;}$ |

Radial magnetic dipole probe (recovers TE):

|  Eq.  |   |
| :---: | :-- |
| (4.26) | $v\\,T\_{1mn} = \\frac{2i}{\\sqrt{6\\pi}}\\,(-1)^m\\,\\frac{1}{n(n+1)}\\,\\Bigl(\\frac{kA}{h\_n^{(1)}(kA)}\\Bigr)^2\\int w\_r^m\\,\\vec{F}^{(3)}\_{2,-m,n}\\cdot\\hat r\\,\\sin\\theta\\,d\\theta\\,d\\phi$ |

> **Drawback.** Radial field components decay as $r^{-2}$ while tangential decay as $r^{-1}$. At large measurement distances, the radial measurement is dominated by tangential leakage — **probe-alignment-sensitive**.

### 2.4 Two-tangential-component measurement (§4.2.4) — uses (A1.69)

$\\hat\\theta$ and $\\hat\\phi$ dipole probes:

|  Eq.  |   |
| :---: | :-- |
| (4.30) | $\\boxed{\\;v\\,T\_{smn} = \\frac{2}{\\sqrt{6\\pi}}\\,(-1)^m\\,\\{R^{(3)}\_{sn}(kA)\\}^{-2}\\int\_0^{2\\pi}\\!\\!\\int\_0^\\pi \\{w^e\_\\theta\\,\\hat\\theta + w^e\_\\phi\\,\\hat\\phi\\}\\cdot\\vec{F}^{(3)}\_{s,-m,n}(A,\\theta,\\phi)\\,\\sin\\theta\\,d\\theta\\,d\\phi\\;}$ |

Magnetic-dipole version:

|  Eq.  |   |
| :---: | :-- |
| (4.31) | $v\\,T\_{smn} = \\frac{2i}{\\sqrt{6\\pi}}\\,(-1)^m\\,\\{R^{(3)}\_{3-s,n}(kA)\\}^{-2}\\int \\{w^m\_\\theta\\,\\hat\\theta + w^m\_\\phi\\,\\hat\\phi\\}\\cdot\\vec{F}^{(3)}\_{3-s,-m,n}\\,\\sin\\theta\\,d\\theta\\,d\\phi$ |

### 2.5 Wood's method (§4.2.5) — eccentric measurement sphere

Based on reciprocity integral (A1.74). Allows test-antenna minimum sphere to be offset from measurement-sphere center:

|  Eq.  |   |
| :---: | :-- |
| (4.36, 4.37) | $v\\,T\_{smn} = \\frac{ik^2}{\\sqrt{6\\pi}}\\,(-1)^m\\int\_S\\Bigl\\{[w\_\\theta^e \\hat\\phi - w\_\\phi^e \\hat\\theta]\\cdot\\vec{F}^{(4)}\_{3-s,-m,n}(r,\\theta,\\phi) + i[w\_\\theta^m \\hat\\phi - w\_\\phi^m \\hat\\theta]\\cdot\\vec{F}^{(4)}\_{s,-m,n}(r,\\theta,\\phi)\\Bigr\\}\\,A^2\\sin\\theta'\\,d\\theta'\\,d\\phi'$ |

> **Requires** simultaneous measurement of two electric and two magnetic tangential components — practically replaced by Huygens-source approximations.

---

## 3. Measurement With Probe Correction — Analytical Solution (§4.3.2)

### 3.1 Starting point — transmission formula

|  Eq.  |   |
| :---: | :-- |
| (4.1) | $\\boxed{\\;w(A, \\chi, \\theta, \\phi) = \\frac{v}{2}\\sum\_{smn,\\sigma\\mu\\nu} T\_{smn}\\,e^{im\\phi}\\,d^n\_{\\mu m}(\\theta)\\,e^{i\\mu\\chi}\\,C^{sn(3)}\_{\\sigma\\mu\\nu}(kA)\\,R^p\_{\\sigma\\mu\\nu}\\;}$ |

### 3.2 Compactification via probe response constants

|  Eq.  |   |
| :---: | :-- |
| (4.39) | $\\boxed{\\;P\_{s\\mu n}(kA) = \\frac{1}{2}\\sum\_{\\sigma\\nu} C^{sn(3)}\_{\\sigma\\mu\\nu}(kA)\\,R^p\_{\\sigma\\mu\\nu}\\;}$ |

Transmission formula becomes:

|  Eq.  |   |
| :---: | :-- |
| (4.40) | $\\boxed{\\;w(A, \\chi, \\theta, \\phi) = v\\sum\_{smn,\\mu} T\_{smn}\\,e^{im\\phi}\\,d^n\_{\\mu m}(\\theta)\\,e^{i\\mu\\chi}\\,P\_{s\\mu n}(kA)\\;}$ |

> **Implementation.** Precompute $P\_{s\\mu n}(kA)$ once per probe/distance; it doesn't depend on $(\\chi, \\theta, \\phi)$. **Major speedup.**

### 3.3 The three-step Fourier inversion

Three orthogonalities are used in sequence:
- (4.41) $\\int e^{i(m-m')\\psi}\\,d\\psi = 2\\pi\\delta\_{mm'}$ — twice
- (A2.10) $\\int d^n\_{\\mu m}(\\theta)\\,d^{n'}\_{\\mu m}(\\theta)\\,\\sin\\theta\\,d\\theta = \\frac{2}{2n+1}\\delta\_{nn'}$ — once

**Step 1: $\\chi$ integral.** Write the transmission formula as a Fourier series in $\\chi$:

|  Eq.  |   |
| :---: | :-- |
| (4.42) | $w(A, \\chi, \\theta, \\phi) = \\sum\_{\\mu = -\\nu\_{\\max}}^{\\nu\_{\\max}} w\_\\mu(A, \\theta, \\phi)\\,e^{i\\mu\\chi}$ |
| (4.43) | $w\_\\mu(A, \\theta, \\phi) = v\\sum\_{s=1}^2 \\sum\_{n=1}^N \\sum\_{m=-n}^n T\_{smn}\\,e^{im\\phi}\\,d^n\_{\\mu m}(\\theta)\\,P\_{s\\mu n}(kA)$ |

Inversion:

|  Eq.  |   |
| :---: | :-- |
| (4.45) | $\\boxed{\\;w\_\\mu(A, \\theta, \\phi) = \\frac{1}{2\\pi}\\int\_0^{2\\pi} w(A, \\chi, \\theta, \\phi)\\,e^{-i\\mu\\chi}\\,d\\chi\\;}$ |

**Step 2: $\\phi$ integral.** Now $w\_\\mu(A, \\theta, \\phi) = \\sum\_{m=-N}^N w\_{\\mu m}(A, \\theta)\\,e^{im\\phi}$ (4.46), with:

|  Eq.  |   |
| :---: | :-- |
| (4.47) | $w\_{\\mu m}(A, \\theta) = v\\sum\_{s=1}^2 \\sum\_{\\substack{n=\|m\|\\\\ n \\ne 0}}^N T\_{smn}\\,d^n\_{\\mu m}(\\theta)\\,P\_{s\\mu n}(kA)$ |

Inversion:

|  Eq.  |   |
| :---: | :-- |
| (4.48) | $\\boxed{\\;w\_{\\mu m}(A, \\theta) = \\frac{1}{2\\pi}\\int\_0^{2\\pi} w\_\\mu(A, \\theta, \\phi)\\,e^{-im\\phi}\\,d\\phi\\;}$ |

**Step 3: $\\theta$ integral.** Now $w\_{\\mu m}(A, \\theta) = \\sum\_{n = |m|, n \\ne 0}^N w^n\_{\\mu m}(A)\\,d^n\_{\\mu m}(\\theta)$ (4.49), with:

|  Eq.  |   |
| :---: | :-- |
| (4.50) | $w^n\_{\\mu m}(A) = v\\sum\_{s=1}^2 T\_{smn}\\,P\_{s\\mu n}(kA)$ |

Inversion via (A2.10):

|  Eq.  |   |
| :---: | :-- |
| (4.51) | $\\boxed{\\;w^n\_{\\mu m}(A) = \\frac{2n+1}{2}\\int\_0^\\pi w\_{\\mu m}(A, \\theta)\\,d^n\_{\\mu m}(\\theta)\\,\\sin\\theta\\,d\\theta\\;}$ |

### 3.4 Final solve for $T\_{smn}$

For each $(m, n)$, write (4.50) explicitly:

|  Eq.  |   |
| :---: | :-- |
| (4.52) | $v\\,T\_{1mn}\\,P\_{1\\mu n}(kA) + v\\,T\_{2mn}\\,P\_{2\\mu n}(kA) = w^n\_{\\mu m}(A)$ |

For a **$\\mu = \\pm 1$ probe** (the standard case), this gives a 2×2 system:

|  Eq.  |   |
| :---: | :-- |
| (4.53) | $\\boxed{\\;v\\,T\_{1mn}\\,P\_{11n}(kA) + v\\,T\_{2mn}\\,P\_{21n}(kA) = w^n\_{1m}(A)\\;}$ |
| (4.54) | $\\boxed{\\;v\\,T\_{1mn}\\,P\_{1,-1,n}(kA) + v\\,T\_{2mn}\\,P\_{2,-1,n}(kA) = w^n\_{-1, m}(A)\\;}$ |

Solved by Cramer's rule (or direct 2×2 inverse) for each $(m, n)$ pair.

---

## 4. Discrete Solution — Practical Implementation (§4.3.3)

### 4.1 Chi integration (§4.3.3.2)

For a $\\mu = \\pm 1$ probe, only $w\_1$ and $w\_{-1}$ are non-zero. Two samples $\\chi = 0$ and $\\chi = \\pi/2$ suffice:

|  Eq.  |   |
| :---: | :-- |
| (4.65, 4.66) | $\\boxed{\\;w\_{\\pm 1}(A, \\theta, \\phi) = \\frac{1}{2}\\Bigl\\{w(A, 0, \\theta, \\phi) \\mp i\\,w\\Bigl(A, \\frac{\\pi}{2}, \\theta, \\phi\\Bigr)\\Bigr\\}\\;}$ |

(Verify: $\\chi = 0$ value gives $w\_1 + w\_{-1}$; $\\chi = \\pi/2$ gives $i(w\_1 - w\_{-1})$; solve.)

> **Sign check.** $w\_{+1}$ corresponds to **$-i$** in the bracket; $w\_{-1}$ to **$+i$**. Mnemonic: $e^{i\\chi}$ at $\\chi = \\pi/2$ gives $+i$, so to extract $w\_{+1}$ subtract; to extract $w\_{-1}$ add.

Alternative four-sample DFT scheme (Eq. 4.62) — slight oversampling, but exploits FFT.

### 4.2 Phi integration (§4.3.3.3)

$J\_\\phi$ equispaced samples in $0 \\le \\phi < 2\\pi$, $\\Delta\\phi = 2\\pi/J\_\\phi$, $J\_\\phi \\ge 2N+1$:

|  Eq.  |   |
| :---: | :-- |
| (4.68) | $\\boxed{\\;\\{w\_{\\mu m}(A, \\theta) \\mid m = 0, 1, \\dots, N, -N, \\dots, -1\\} = \\text{IDFT}\\{w\_\\mu(A, \\theta, j\\Delta\\phi) \\mid j = 0, 1, \\dots, J\_\\phi - 1\\}\\;}$ |

**Layout convention** (per Appendix A4): IDFT output places positive $m$ in bins $0, 1, \\dots, N$ and negative $m$ in bins $J\_\\phi - N, \\dots, J\_\\phi - 1$. **The IDFT here is Hansen's IDFT (numpy's FFT / $J$).**

### 4.3 Theta integration (§4.3.3.4) — the hard part

Theta integrand is **not periodic** on $[0, \\pi]$. Trick: extend $w\_{\\mu m}(A, \\theta)$ to $[0, 2\\pi]$ with matching parity:

**(4.72)**

```math
\tilde w_{\mu m}(A, \theta) = \begin{cases} w_{\mu m}(A, \theta), & 0 \le \theta \le \pi \\ w_{\mu m}(A, 2\pi - \theta), & \pi < \theta < 2\pi,\ (\mu - m)\text{ even} \\ -w_{\mu m}(A, 2\pi - \theta), & \pi < \theta < 2\pi,\ (\mu - m)\text{ odd}\end{cases}
```

> **Parity rule.** Each $w\_{\\mu m}$ sequence has parity $(\\mu - m)$ under $\\theta \\to 2\\pi - \\theta$. This follows from the parity of $d^n\_{\\mu m}$ about $\\pi$ (Appendix A2 — same parity as $\\mu + m$, which matches $\\mu - m$ mod 2).

Now expand into a Fourier series:

|  Eq.  |   |
| :---: | :-- |
| (4.73) | $\\tilde w\_{\\mu m}(A, \\theta) = \\sum\_{l=-N}^N b\_l^{\\mu m}\\,e^{il\\theta},\\quad 0 \\le \\theta < 2\\pi$ |

Compute coefficients via IDFT with $J\_\\theta \\ge 2N + 1$ samples:

|  Eq.  |   |
| :---: | :-- |
| (4.77) | $\\boxed{\\;\\{b\_l^{\\mu m} \\mid l = 0, 1, \\dots, N, -N, \\dots, -1\\} = \\text{IDFT}\\{\\tilde w\_{\\mu m}(A, j\\Delta\\theta) \\mid j = 0, 1, \\dots, J\_\\theta - 1\\}\\;}$ |

### 4.4 Closed-form theta integral via $\\Pi(l - m')$ (§4.3.3.5)

Substituting (4.73) and the Fourier expansion (4.70) of $d^n\_{\\mu m}$:

|  Eq.  |   |
| :---: | :-- |
| (4.70) | $d^n\_{\\mu m}(\\theta) = i^{\\mu - m}\\sum\_{m'=-n}^n \\Delta^n\_{m'\\mu}\\,\\Delta^n\_{m'm}\\,e^{-im'\\theta}$ |

into (4.51), the theta integral becomes algebraic:

|  Eq.  |   |
| :---: | :-- |
| (4.75) | $w^n\_{\\mu m}(A) = \\frac{2n+1}{2}\\,i^{\\mu-m}\\sum\_{l=-N}^N b\_l^{\\mu m}\\sum\_{m'=-n}^n \\Delta^n\_{m'\\mu}\\,\\Delta^n\_{m'm}\\,G(l - m')$ |

where:

**(4.76)**

```math
G(l - m') = \int_0^\pi e^{i(l-m')\theta}\sin\theta\,d\theta = \begin{cases} \pm i\pi/2, & l - m' = \pm 1 \\ 0, & |l - m'| = 3, 5, 7, \dots \\ 2/[1 - (l-m')^2], & |l - m'| = 0, 2, 4, \dots\end{cases}
```

### 4.5 Parity reduction and final form (§4.3.3.5)

Using $b\_l^{\\mu m} = (-1)^{\\mu + m}\\,b\_{-l}^{\\mu m}$ (Eq. 4.79) and $\\Delta^n\_{m'\\mu}\\Delta^n\_{m'm} = (-1)^{\\mu+m}\\Delta^n\_{-m',\\mu}\\Delta^n\_{-m',m}$ (Eq. 4.80, 4.81), the $l = \\pm 1$ terms in $G$ **cancel**, giving:

|  Eq.  |   |
| :---: | :-- |
| (4.83) | $\\boxed{\\;w^n\_{\\mu m}(A) = \\frac{2n+1}{2}\\,i^{\\mu-m}\\sum\_{m'=-n}^n \\Delta^n\_{m'\\mu}\\,\\Delta^n\_{m'm}\\,\\sum\_{l=-N}^N \\Pi(l - m')\\,b\_l^{\\mu m}\\;}$ |

with:

**(4.84)**

```math
\Pi(l - m') = \begin{cases} 0, & (l - m')\text{ odd} \\ 2 / [1 - (l-m')^2], & (l - m')\text{ even}\end{cases}
```

> $\\Pi(l - m') = \\Pi(m' - l)$ — even function of its argument.

### 4.6 The $K(m')$ convolution (§4.3.3.5)

Define the inner $l$-sum:

|  Eq.  |   |
| :---: | :-- |
| (4.85) | $K(m') = \\sum\_{l=-N}^N \\Pi(l - m')\\,b\_l^{\\mu m},\\quad -N \\le m' \\le N$ |

This is a **convolution**. To use FFT, extend both sequences to period $4N$:

|  Eq.  |   |
| :---: | :-- |
| (4.86) | $\\tilde\\Pi(j) = \\Pi(j),\\quad -2N < j \\le 2N,\\quad \\tilde\\Pi(j) = \\tilde\\Pi(j + c\\cdot 4N)$ |

**(4.87)**

```math
\tilde b_l^{\mu m} = \begin{cases} b_l^{\mu m}, & -N \le l \le N \\ 0, & -2N < l < -N\text{ and }N < l \le 2N\end{cases},\quad \tilde b_l^{\mu m} = \tilde b_{l + c\cdot 4N}^{\mu m}
```

Then:

|  Eq.  |   |
| :---: | :-- |
| (4.88) | $K(m') = \\sum\_{l=0}^{4N-1} \\tilde\\Pi(l - m')\\,\\tilde b\_l^{\\mu m}$ |

Computed via DFT:

|  Eq.  |   |
| :---: | :-- |
| (4.89) | $\\boxed{\\;K(m') = \\text{IDFT}\\{\\text{DFT}\\{\\tilde\\Pi(i)\\}\\cdot \\text{DFT}\\{\\tilde b\_j^{\\mu m}\\}\\}\\;}$ |

(Both sequences of length $4N$. $\\text{DFT}\\{\\tilde\\Pi\\}$ is precomputed **once** per $N$.)

**Parity exploitation:** $\\tilde b\_l^{\\mu m} = (-1)^{\\mu + m}\\,\\tilde b\_{-l}^{\\mu m}$ (Eq. 4.90) ⟹ $K(m') = (-1)^{\\mu+m}\\,K(-m')$ (Eq. 4.91). Halves storage and work.

### 4.7 Final $m'$ summation

|  Eq.  |   |
| :---: | :-- |
| (4.92) | $w^n\_{\\mu m}(A) = \\frac{2n+1}{2}\\,i^{\\mu - m}\\sum\_{m'=-n}^n \\Delta^n\_{m'\\mu}\\,\\Delta^n\_{m'm}\\,K(m'),\\quad \\mu = \\pm 1$ |

Using $\\Delta^n\_{-m',\\mu}\\Delta^n\_{-m',m}\\,K(-m') = \\Delta^n\_{m'\\mu}\\Delta^n\_{m'm}\\,K(m')$ (Eq. 4.93), reduce the sum to $0 \\le m' \\le n$ (multiply by 2 for $m' > 0$).

### 4.8 Delta recurrence

|  Eq.  |   |
| :---: | :-- |
| (4.94) | $\\sqrt{(n+m'+1)(n-m')}\\,\\Delta^n\_{m'+1,m} + \\sqrt{(n+m')(n-m'+1)}\\,\\Delta^n\_{m'-1, m} + 2m\\,\\Delta^n\_{m'm} = 0$ |

**Run backwards from $m' = n$ to $m' = 0$** (stable; see Appendix A2 §5). Seed value:

|  Eq.  |   |
| :---: | :-- |
| (4.95) | $\\Delta^n\_{nm} = 2^{-n}\\sqrt{\\frac{2n(2n-1)\\cdots(n-m+1)}{(n-m)!}}$ |

Compute $\\Delta^n\_{m'\\mu}$ and $\\Delta^n\_{m'm}$ **simultaneously** (same square roots) — saves ~half the work.

---

## 5. Field Transformations and Antenna Parameters (§4.3.4)

### 5.1 Input/output probe concept (§4.3.4.1)

Once $T\_{smn}$ is recovered, the transmission formula (4.40) re-applied with a new probe and distance $A'$ computes the field on any other sphere:

|  Eq.  |   |
| :---: | :-- |
| (4.97) | $w'(A', \\chi, \\theta, \\phi) = \\sum\_{smn,\\mu=\\pm 1} v\\,T\_{smn}\\,e^{im\\phi}\\,d^n\_{\\mu m}(\\theta)\\,e^{i\\mu\\chi}\\,P'\_{s\\mu n}(kA')$ |

| Probe role | Symbol | Purpose |
|---|---|---|
| Input probe | $P\_{s\\mu n}(kA)$ | Used in inversion (recovers $T\_{smn}$) |
| Output probe | $P'\_{s\\mu n}(kA')$ | Used in synthesis (evaluates field at $A'$) |

### 5.2 Far-field limits (§4.3.4.2)

For $kA \\to \\infty$, the translation coefficient asymptotics (A3.22–A3.24) give:

|  Eq.  |   |
| :---: | :-- |
| (4.100) | $P\_{s,\\pm 1, n}(kA) \\to \\frac{e^{ikA}}{kA}\\,P^\\infty\_{s,\\pm 1, n},\\quad kA \\to \\infty$ |

**General linearly-polarized $\\hat x'$ probe:**

|  Eq.  |   |
| :---: | :-- |
| (4.101) | $P^\\infty\_{s1n} = \\tfrac{1}{4}\\sqrt{2n+1}\\,i^{-n-1}\\sum\_{\\nu=1}^{\\nu\_{\\max}}\\sqrt{2\\nu+1}\\,i^\\nu\\,\\{R^p\_{11\\nu} + R^p\_{21\\nu}\\}$ |
| (4.102) | $P^\\infty\_{s,-1,n} = (-1)^{s+1}\\,P^\\infty\_{s1n}$ |

**Closed forms for electric and magnetic Hertzian-dipole probes:**

| Probe | $P^\\infty\_{s1n}$ | $P^\\infty\_{s,-1,n}$ |
|---|---|---|
| $\\hat x'$ electric dipole | $-\\tfrac{\\sqrt 6}{8}\\sqrt{2n+1}\\,i^{-n}$ (4.103) | $(-1)^s\\,\\tfrac{\\sqrt 6}{8}\\sqrt{2n+1}\\,i^{-n}$ (4.104) |
| $\\hat x'$ magnetic dipole (note: $\\hat y'$-polarized!) | $\\tfrac{\\sqrt 6}{8}\\sqrt{2n+1}\\,i^{-n+1}$ (4.105) | $(-1)^s\\,\\tfrac{\\sqrt 6}{8}\\sqrt{2n+1}\\,i^{-n+1}$ (4.106) |

**Normalized far-field signal:**

|  Eq.  |   |
| :---: | :-- |
| (4.107) | $W(\\chi, \\theta, \\phi) = \\lim\_{kA\\to\\infty}\\Bigl[w(A, \\chi, \\theta, \\phi)\\,\\frac{kA}{e^{ikA}}\\Bigr]$ |

Far-field transmission formula:

|  Eq.  |   |
| :---: | :-- |
| (4.108) | $W(\\chi, \\theta, \\phi) = \\sum\_{smn, \\mu=\\pm 1} v\\,T\_{smn}\\,e^{im\\phi}\\,d^n\_{\\mu m}(\\theta)\\,e^{i\\mu\\chi}\\,P^\\infty\_{s\\mu n}$ |
| (4.109) | $W'(\\chi, \\theta, \\phi) = \\sum\_{smn, \\mu=\\pm 1} v\\,T\_{smn}\\,e^{im\\phi}\\,d^n\_{\\mu m}(\\theta)\\,e^{i\\mu\\chi}\\,P'^\\infty\_{s\\mu n}$ |

### 5.3 Relative measurements (§4.3.4.3) — propagation of unknown constants

| Input | Probe constants | Algorithm produces |
|---|---|---|
| Absolute $w$ | Absolute $P$ | $v\\,T\_{smn}$, $w'$ |
| Relative ($c\_1 w$) | Absolute ($c\_2 = 1$) | $c\_1\\,v T\_{smn}$, $c\_1\\,w'$ |
| Absolute | Relative ($c\_2$) | $v T\_{smn}/c\_2$, $w'/c\_2$ |
| Both relative | | $vT\_{smn}\\,c\_1/c\_2$, $w'\\,c\_1/c\_2 \\cdot c\_3$ |

(With $c\_3$ from unknown output-probe scaling.)

### 5.4 Antenna parameter extraction (§4.3.4.4)

**Directivity (relative measurement, output = far-field $\\hat x'$ electric dipole):**

|  Eq.  |   |
| :---: | :-- |
| (4.112) | $\\boxed{\\;D\_t(\\theta, \\phi) = \\frac{8}{3}\\cdot\\frac{\|W'\\,c\_1/c\_2\|^2}{\\sum\_{smn}\|v T\_{smn}\\,c\_1/c\_2\|^2}\\;}$ |

The unknown $c\_1/c\_2$ **cancels** — directivity is recoverable from purely relative measurements.

**Gain (requires measurement of $c\_1 v$ separately):**

|  Eq.  |   |
| :---: | :-- |
| (4.115) | $\\boxed{\\;G\_t = \\frac{8}{3}\\cdot\\frac{\|c\_1 W'\|^2}{\|c\_1 v\|^2}\\;}$ |

**EIRP:**

|  Eq.  |   |
| :---: | :-- |
| (4.116) | $\\text{EIRP}(\\theta, \\phi) = \\tfrac{1}{2}\|v\|^2\\,G\_t(\\theta, \\phi)$ |

**Absolute far-field pattern:**

|  Eq.  |   |
| :---: | :-- |
| (4.117) | $\\boxed{\\;\\vec K(\\theta, \\phi) = \\frac{2\\sqrt 6}{3v}\\bigl\\{W'(0, \\theta, \\phi)\\,\\hat\\theta + W'(\\pi/2, \\theta, \\phi)\\,\\hat\\phi\\bigr\\}\\;}$ |

With a general directive output probe of gain $G\_p$:

|  Eq.  |   |
| :---: | :-- |
| (4.118) | $\\vec K(\\theta, \\phi) = \\frac{2\\sqrt 6}{3v}\\bigl\\{W'(0)\\,\\hat\\theta + W'(\\pi/2)\\,\\hat\\phi\\bigr\\}\\sqrt{G\_e/G\_p},\\quad G\_e = 3/2$ |

**Relative far-field pattern** (no $v$ knowledge needed):

|  Eq.  |   |
| :---: | :-- |
| (4.119) | $\\vec K\_\\text{rel}(\\theta, \\phi) = c\\,\\{W'(0, \\theta, \\phi)\\,\\hat\\theta + W'(\\pi/2, \\theta, \\phi)\\,\\hat\\phi\\}$ |

**Tangential field at radius $A'$:**

|  Eq.  |   |
| :---: | :-- |
| (4.120) | $\\vec E\_\\text{tang}(A', \\theta, \\phi) = \\frac{2k}{\\sqrt{6\\pi\\eta}}\\bigl\\{w'^e(A', 0)\\,\\hat\\theta + w'^e(A', \\pi/2)\\,\\hat\\phi\\bigr\\}$ |
| (4.121) | $\\vec H\_\\text{tang}(A', \\theta, \\phi) = \\frac{2k\\sqrt\\eta}{\\sqrt{6\\pi}}\\bigl\\{w'^m(A', 0)\\,\\hat\\theta + w'^m(A', \\pi/2)\\,\\hat\\phi\\bigr\\}$ |

---

## 6. SNIFT Algorithm — Block-by-Block (§4.4)

### 6.1 Four modes of operation

| Mode | Input | Output | Use case |
|---|---|---|---|
| (1) $w \\to w'$ | Near-field samples | Near-field at $A'$ | Field at different radius |
| (2) $w \\to W'$ | Near-field samples | Far-field pattern | **Standard NF-to-FF** |
| (3) $W \\to w'$ | Far-field samples | Near-field at $A'$ | FF-to-NF |
| (4) $W \\to W'$ | Far-field samples | Far-field pattern | Coordinate/polarization filtering |

### 6.2 First part — recovery of $T\_{smn}$ (Steps 1–7)

Pipeline:

| Step | Operation | Equation |
|---|---|---|
| 1 | $\\chi$ integration: $w\_{\\pm 1}(A, \\theta, \\phi) = \\tfrac{1}{2}\\{w(A, 0, \\theta, \\phi) \\mp i\\,w(A, \\pi/2, \\theta, \\phi)\\}$ | (4.126) |
| 2 | $\\phi$ IDFT: $\\{w\_{\\pm 1, m}(A, \\theta)\\} = \\text{IDFT}\\{w\_{\\pm 1}(A, \\theta, j\\Delta\\phi)\\}$ | (4.127) |
| 3 | Extend data $\\tilde w\_{\\mu m}$ per (4.72); $\\theta$ IDFT: $\\{b\_l^{\\pm 1, m}\\} = \\text{IDFT}\\{\\tilde w\_{\\pm 1, m}(A, j\\Delta\\theta)\\}$ | (4.128) |
| 4 | Convolve via FFT: $K(m') = \\text{IDFT}\\{\\text{DFT}\\{\\tilde\\Pi\\}\\cdot\\text{DFT}\\{\\tilde b^{\\pm 1, m}\_j\\}\\}$ | (4.131) |
| 5 | For each $n$: compute $w^n\_{\\mu m}(A) = \\tfrac{2n+1}{2}\\,i^{\\mu-m}\\sum\_{m'=-n}^n \\Delta^n\_{m'\\mu}\\Delta^n\_{m'm}\\,K(m')$ for $\\mu = \\pm 1$ | (4.132) |
| 6 | Recurrence-compute $\\Delta^n\_{m'\\mu}$, $\\Delta^n\_{m'm}$ for $0 \\le m' \\le n$ (using 4.94, 4.95) | — |
| 7 | Solve 2×2 system: $T\_{1mn}, T\_{2mn}$ from (4.53, 4.54) | (4.133, 4.134) |

### 6.3 Second part — evaluation of $w'$ (Steps 8–11)

Reverse pipeline using already-computed delta products:

|  Eq.  |   |
| :---: | :-- |
| (4.135) | $w'(A', \\chi, \\theta, \\phi) = \\sum\_{\\mu = \\pm 1} e^{i\\mu\\chi}\\sum\_{m=-N}^N e^{im\\phi}\\sum\_{m'=-N}^N e^{im'\\theta}\\sum\_{n=\\max(\|m'\|, \|m\|, 1)}^N \\Delta^n\_{m'\\mu}\\Delta^n\_{m'm}\\,i^{m-\\mu}\\sum\_{s=1}^2 v\\,T\_{smn}\\,P'\_{s\\mu n}(kA')$ |

| Step | Operation |
|---|---|
| 8 | For each $n$: accumulate $\\Delta^n\_{m'\\mu}\\Delta^n\_{m'm}\\,i^{m-\\mu}\\sum\_s v T\_{smn} P'\_{s\\mu n}$ into the $n$-sum array |
| 9 | $m'$ summation via DFT (yields $\\theta$-grid output) |
| 10 | $m$ summation via DFT (yields $\\phi$-grid output) |
| 11 | $\\mu$ summation (yields final $\\chi$-dependence: $\\chi = 0, \\pi/2, \\dots$) |

### 6.4 Input parameters (Table 4.4)

| Parameter | Meaning |
|---|---|
| $J\_\\theta$ | # samples in $0 \\le \\theta < 2\\pi$ of input |
| $J\_\\phi$ | # samples in $0 \\le \\phi < 2\\pi$ of input |
| $A/\\lambda$ | Input sphere radius (skipped if input is $W$) |
| $\\nu\_{\\max}$ | Max $\\nu$ for input probe |
| $N$ | Max $n$ for test antenna ($1 \\le N \\le (J\_\\theta-1)/2$) |
| $M$ | Max $\\|m\\|$ (default $N$; $0 \\le M \\le \\min((J\_\\phi-1)/2, N)$) |
| $J'\_\\theta, J'\_\\phi, A'/\\lambda, \\nu'\_{\\max}$ | Same for output |
| Mode | (1), (2), (3), or (4) |

### 6.5 Sampling counts (§4.4.3)

| Sampling strategy | # samples $N\_s$ on full sphere |
|---|---|
| Uniform $(2N+1)$ × $(2N+1)$, 2 χ values | $N\_1 = 2(2N+1)(N+1)$ — Eq. (4.137) |
| Thinned (sin θ scaling near poles) | $N\_2 = (4/\\pi)(2N+1)(N+1)$ — Eq. (4.138) |
| Cylindrical-symmetric antenna ($M = kr\_c + 10$) | $N\_3 = 2(2M+1)(N+1)$ — Eq. (4.139) |

Total modes $N\_0 = 2N(N+2)$ (Eq. 4.136). Uniform sampling is ~50% efficient; thinning yields ~$N\_0$.

---

## 7. Worked Example (§4.4.5) — Critical Regression Test

**Input:** $W(\\chi, \\theta, \\phi) = (20\\cos 2\\theta + 32\\cos\\theta + 12)\\cos(\\chi + \\phi)$, band-limited with $N = 2, M = 1$.

**Setup:**
- $\\theta$ samples at $0°, 60°, 120°, 180°$ (4 samples in natural range; 2 generated to fill $0°$–$360°$)
- $\\phi$ samples at $0°, 90°, 180°, 270°$
- $\\chi$ samples at $0°, 90°$
- Input/output probes: $\\hat x'$-polarized maximum-directivity antenna with $N = 2$
- All radii infinite — only $P^\\infty$ and $P'^\\infty$ used (no Hankel functions)

**Probe response constants** (Eq. 4.151):

|  Eq.  |   |
| :---: | :-- |
|  | $P^\\infty\_{111} = P^\\infty\_{211} = P^\\infty\_{1,-1,1} = -\\frac{\\sqrt 6}{2}i,\\quad P^\\infty\_{2,-1,1} = \\frac{\\sqrt 6}{2}i$ |
|  | $P^\\infty\_{112} = P^\\infty\_{212} = P^\\infty\_{1,-1,2} = -\\frac{\\sqrt{10}}{2},\\quad P^\\infty\_{2,-1,2} = \\frac{\\sqrt{10}}{2}$ |

**Delta products** (Eq. 4.152):

|  Eq.  |   |
| :---: | :-- |
|  | $\\Delta^1\_{01}\\Delta^1\_{01} = \\tfrac{1}{2},\\ \\Delta^2\_{01}\\Delta^2\_{01} = 0$ |
|  | $\\Delta^1\_{11}\\Delta^1\_{11} = \\tfrac{1}{4},\\ \\Delta^2\_{11}\\Delta^2\_{11} = \\tfrac{1}{4}$ |
|  | $\\Delta^2\_{21}\\Delta^2\_{21} = \\tfrac{1}{4}$ |

**Step 1 — χ integration** (Table 4.5b — partial):

| | $\\phi=0°, \\mu=+1$ | $\\phi=0°, \\mu=-1$ |
|---|---|---|
| $\\theta = 0°$ | 32 | 32 |
| $\\theta = 60°$ | 9 | 9 |
| $\\theta = 120°$ | $-7$ | $-7$ |
| $\\theta = 180°$ | 0 | 0 |

**Step 2 — φ integration** (Table 4.5c) yields non-zero only at $m = \\pm 1$:

| $\\theta$ | $W\_{11}(\\theta)$ | $W\_{-1,-1}(\\theta)$ |
|---|---|---|
| $0°$ | 32 | 32 |
| $60°$ | 9 | 9 |
| $120°$ | $-7$ | $-7$ |
| $180°$ | 0 | 0 |

**Step 3 — extend and IDFT** (with parity reasoning):

|  Eq.  |   |
| :---: | :-- |
| (4.143) | $\\tilde W\_{11} = \\{32, 9, -7, 0, -7, 9\\}\\quad\\Rightarrow\\quad b\_l^{11} = \\{6, 8, 5, 0, 5, 8\\}$ |

Zero-padded to length 12:

|  Eq.  |   |
| :---: | :-- |
| (4.144) | $\\tilde b\_l^{11} = \\{6, 8, 5, 0, 0, 0, 0, 0, 0, 0, 5, 8\\}$ |

**Step 12 — Π sequence (length 12)**:

|  Eq.  |   |
| :---: | :-- |
| (4.145) | $\\tilde\\Pi = \\{2, 0, -\\tfrac{2}{3}, 0, -\\tfrac{2}{15}, 0, -\\tfrac{2}{35}, 0, -\\tfrac{2}{15}, 0, -\\tfrac{2}{3}, 0\\}$ |

**Step 4 — Convolution via FFT:**

|  Eq.  |   |
| :---: | :-- |
| (4.150) | $K(m') \\to \\{K(0), K(1), K(2)\\} = \\{\\tfrac{16}{3}, \\tfrac{64}{3}, \\tfrac{32}{3}\\}\\quad\\text{(after doubling for } m' > 0\\text{)}$ |

**Step 5–7 — solve for $T\_{smn}$** at each $n$:

| $n$ | Equation system | Solution |
|---|---|---|
| 1 | $\\{T\_{111}, T\_{211}\\} = \\{12, 0\\}/(- \\sqrt 6/2\\cdot i)$-matrix | $v T\_{111} = v T\_{211} = v T\_{1,-1,1} = -v T\_{2,-1,1} = (\\sqrt 6/2)i$ |
| 2 | $\\{T\_{112}, T\_{212}\\} = \\{20, 0\\}/\\dots$ | $v T\_{112} = v T\_{212} = v T\_{1,-1,2} = -v T\_{2,-1,2} = -2\\sqrt{10}$ |

**Step 8–11 — back-substitution to compute $W'$**:

The final output table (Table 4.7c) reproduces the input $W$ exactly:

| | $\\phi=0°, \\chi=0$ | $\\phi=0°, \\chi=\\pi/2$ |
|---|---|---|
| $\\theta = 0°$ | 64 | 0 |
| $\\theta = 60°$ | 18 | 0 |
| $\\theta = 120°$ | $-14$ | 0 |
| $\\theta = 180°$ | 0 | 0 |

> **Validation power.** Reproducing every table value (4.5a-c, 4.6a-c, 4.7a-c) of the worked example is the **gold-standard regression test** for a Python implementation. Every intermediate quantity is integer or simple radical — exact arithmetic is possible.

---

## 8. Reference Test Case Errors (§4.4.6)

12-element dipole array antenna (50λ diameter), 2-element probe, $A = 50\\lambda$:

| $J$ (samples/great circle) | 100% error | 90% error | 50% error | 10% error |
|---|---|---|---|---|
| 260 | ~0 dB | $-15$ dB | $-23$ dB | $-25$ dB |
| 300 | $-20$ | $-30$ | $-35$ | $-40$ |
| 350 | $-40$ | $-65$ | $-75$ | $-85$ |
| 360 | $-50$ | $-90$ | $-105$ | $-115$ |
| 400 | $-50$ | $-90$ | $-105$ | $-115$ |

> **Conclusion.** With $J \\ge \\lceil kr\_t \\rceil + 10$ samples per great circle, sampling error vanishes; remaining errors are floating-point rounding (in single precision, ~$-50$ dB; in double, ~$-100$ dB+).

---

## 9. Algorithm Variants — Review (§4.5)

| Method | Reference | Notes |
|---|---|---|
| Simpson rule on θ, φ | James & Longdon | Only practical for $N \\le 3$ |
| Trapezoidal θ + Ludwig iteration | Ludwig | Iterative refinement; corrects rounding |
| Simpson + Wood approximations | Wood | Eccentric measurements via Huygens probes |
| Gaussian quadrature | — | Exact, but non-equidistant samples |
| Ricardi–Burrows recursion | — | Fourier extension + recursion |
| Wacker | [1, 2] | Fourier expand $d^n\_{\\mu m}\\sin\\theta$ jointly; key insight |
| Wacker–Larsen | [3, 4] | Keep $\\sin\\theta$ in integral; reuse $\\Delta$'s for synthesis; convolution form for $K(m')$ |
| Yaghjian–Wittmann | [20, 21, 22] | Tangential-component formulation; halves $\\theta$-integral terms |
| Lewis hemispherical | [12, 23] | Halves work when back hemisphere is zero |

The **Wacker–Larsen** method (§4.3) is what SNIFT implements and what this document treats as canonical.

---

## 10. Implementation Checklist for Python Code

### Setup
1. **Truncation $N$**: choose $N = \\lfloor kr\_0 \\rfloor + 10$ where $r\_0$ is the test-antenna minimum-sphere radius (Eq. 2.31).
2. **Sample counts**: $J\_\\theta \\ge 2N + 1$ for theta in $[0, 2\\pi)$ extended, $J\_\\phi \\ge 2N + 1$ for phi, and 2 chi values ($\\chi = 0, \\pi/2$).
3. **Equispaced grids**: $\\theta\_j = j\\cdot 2\\pi/J\_\\theta$, $\\phi\_j = j\\cdot 2\\pi/J\_\\phi$. **Measurements are taken only at $\\theta\_j \\in [0, \\pi]$** — the second half is computer-generated via (4.72).

### Probe response constants
4. **Precompute $P\_{s\\mu n}(kA)$** once via (4.39); use sympy or recurrence-based 3-j computation per Appendix A3.
5. **Symmetry check**: $P\_{s,-1,n} = (-1)^{s+1}\\,P\_{s,1,n}$ (Eq. 3.27).
6. **Far-field limit**: $|P\_{s,\\pm 1, n}(kA)\\cdot kA/e^{ikA}|$ stable as $kA \\to \\infty$, matching $P^\\infty\_{s,\\pm 1, n}$ from (4.101)–(4.106).

### Step 1: χ integration
7. **Verify (4.65)/(4.66)**: with synthetic $w(\\chi) = c\_+ e^{i\\chi} + c\_- e^{-i\\chi}$, recovery should give $w\_1 = c\_+$, $w\_{-1} = c\_-$.

### Step 2: φ IDFT
8. **Layout convention**: Hansen-IDFT output places $m = 0, 1, \\dots, N$ in bins 0..N and $m = -N, \\dots, -1$ in bins $J\_\\phi - N, \\dots, J\_\\phi - 1$. (Middle bins should be zero for band-limited input.)
9. **Hansen-vs-numpy mapping**: Implement `hansen_idft(x) = np.fft.fft(x) / len(x)` per Appendix A4 §2.

### Step 3: θ extension and IDFT
10. **Parity rule (4.72)**: extension parity is $(\\mu - m)$. Test for each pair: $(\\mu = +1, m = 0) \\to$ odd; $(\\mu = +1, m = 1) \\to$ even; etc.
11. **Verify extension preserves period**: $\\tilde w\_{\\mu m}(0) = \\tilde w\_{\\mu m}(2\\pi)$ (after periodic wrap).
12. **IDFT layout**: same as step 8 but for $l$ bins.

### Step 4: K(m') convolution
13. **$\\tilde\\Pi$ table**: precompute for the chosen $N$. Verify $\\tilde\\Pi(0) = 2$, $\\tilde\\Pi(\\pm 1) = 0$, $\\tilde\\Pi(\\pm 2) = -2/3$, $\\tilde\\Pi(\\pm 4) = -2/15$, etc.
14. **Parity check after convolution**: $K(-m') = (-1)^{\\mu+m}\\,K(m')$ (Eq. 4.91).
15. **Length-$4N$ zero-padding**: assert that bins $\\{-2N..-N-1\\} \\cup \\{N+1..2N\\}$ of $\\tilde b$ are zero.

### Steps 5–6: Delta recurrence
16. **Seed $\\Delta^n\_{nm}$** with (4.95). Verify $\\Delta^n\_{nn} = 2^{-n}$ (Appendix A2 Eq. A2.42).
17. **Backwards recurrence (4.94)**: stable from $m' = n$ down to $m' = 0$. Compare to Appendix A2 §6 explicit tables for $n \\le 5$.
18. **Simultaneous computation**: $\\Delta^n\_{m'\\mu}$ and $\\Delta^n\_{m'm}$ recur with the same square-root coefficients — compute both in a single loop.

### Step 7: 2×2 solve
19. **2×2 system**: condition number should be $O(1)$ for typical probes; flag if $|P\_{11n}P\_{2,-1,n} - P\_{21n}P\_{1,-1,n}|$ is small.

### Synthesis (Steps 8–11)
20. **Delta-product reuse**: same $\\Delta^n\_{m'\\mu}\\Delta^n\_{m'm}$ are used in both analysis and synthesis — don't recompute. Cache during the n-loop.
21. **DFT vs IDFT consistency**: analysis uses IDFT (samples → coefficients); synthesis uses DFT (coefficients → samples). Wrap both around `hansen_idft`/`hansen_dft` per A4 §2.

### End-to-end validation
22. **Worked-example regression test (§7)**: implement the explicit $W = (20\\cos 2\\theta + 32\\cos\\theta + 12)\\cos(\\chi + \\phi)$ input. Verify every table in §4.4.5 reproduces exactly (integer/radical arithmetic permits zero tolerance).
23. **Round-trip identity**: for a synthetic $T\_{smn}$ set, generate $w$ via (4.40), invert to recover $T'\_{smn}$, assert $\\|T' - T\\|/\\|T\\| < 10^{-12}$ (double precision).
24. **Friis far-field check**: at $kA = 10000$ with an electric-dipole probe, the recovered $T\_{smn}$ should produce a far field matching the Friis formula $|w|^2/|v|^2 = G\_p G\_t/(4(kA)^2)$ to 3+ significant figures.
25. **Hertzian-dipole regression**: with a $\\hat z$-electric dipole as test antenna ($T\_4 = 1$ only), all other recovered $T\_{smn}$ must be $0$ to floating-point precision; the recovered $T\_4$ must equal 1.

### Sampling and convergence
26. **Truncation N**: scan $J\_\\theta = J\_\\phi$ from $2N - 10$ up to $4N$ and plot error. Should plateau at $J \\approx 2N + 1$ (Eq. 4.137) and become floating-point-limited beyond.
27. **Antenna parameters**: directivity, gain, EIRP, far-field patterns recovered per (4.112)–(4.118) — closed-form sanity check for known antennas (dipole: $D = 1.5$; Huygens: $D = 3$).

### Performance
28. **FFT length**: $4N$ for $K(m')$ convolution (Eq. 4.89). Choose $4N$ as a product of small primes; pad to $4N' \\ge 4N$ if needed.
29. **Half-storage**: exploit $K(-m') = \\pm K(m')$ (Eq. 4.91) to halve memory.
30. **Parity packing**: pair $m$ even/odd sequences through Fourier ops (Lewis trick, [12]) — cuts FFT count in half.

---

*End of Chapter 4 reference.*
