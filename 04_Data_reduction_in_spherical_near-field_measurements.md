# Chapter 4 — Data Reduction in Spherical Near-Field Measurements

**Source:** Hansen, *Spherical Near-Field Antenna Measurements*, Chapter 4 (pp. 89–146).

**Cross-references:**
- Symbols: [A5_List_of_principal_symbols_and_uses.md](A5_List_of_principal_symbols_and_uses.md)
- Chapter 2 (Scattering matrix, $\vec{F}^{(c)}_{smn}$): [02_Scattering_matrix_description_of_an_antenna.md](02_Scattering_matrix_description_of_an_antenna.md)
- Chapter 3 (Transmission formula): [03_Scattering_matrix_description_of_antenna_coupling.md](03_Scattering_matrix_description_of_antenna_coupling.md)
- Appendix A1 (wave functions, orthogonality integrals): [A1_Spherical_wave_functions_notation_and_properties.md](A1_Spherical_wave_functions_notation_and_properties.md)
- Appendix A2 (rotation, $d^n_{\mu m}$, $\Delta^n_{m'm}$): [A2_Rotation_of_spherical_waves.md](A2_Rotation_of_spherical_waves.md)
- Appendix A3 (translation $C^{sn(c)}_{\sigma\mu\nu}$): [A3_Translation_of_spherical_waves.md](A3_Translation_of_spherical_waves.md)
- Appendix A4 (DFT, sampling, reconstruction): [A4_Data_processing_in_antenna_measurements.md](A4_Data_processing_in_antenna_measurements.md)

> **Purpose.** Validation reference for **the inversion algorithm** (the SNIFT-style transformation): given measured probe signals $w(A, \chi, \theta, \phi)$, recover the test-antenna transmitting coefficients $T_{smn}$, then evaluate fields at any sphere via an output probe. Every formula tagged with its book number. The §4.4.5 worked example provides explicit intermediate values for unit tests.

---

## 1. Problem Statement and Big Picture (§4.1)

An antenna's total scattering matrix $\hat{\mathbf{S}}$ contains four blocks:

| Block | Count | Status |
|---|---|---|
| Reflection $\Gamma$ | 1 | Direct measurement (trivial) |
| Scattering $S_{ij}$ | $J^2$ | Hard; not pursued here |
| Receiving $R_j$ | $J$ | Recovered if reciprocity used |
| Transmitting $T_j$ | $J$ | **Main goal** of this chapter |

where $J = 2N(N+2)$ and $N = \lceil kr_0 \rceil + \sim 10$ (Eq. 2.31).

### Two approaches

| Approach | Probe model | When valid |
|---|---|---|
| **With probe correction** (§4.3) | Realistic horn; full $R^p_{\sigma\mu\nu}$ | All cases; the standard method |
| **Without probe correction** (§4.2) | Ideal short electric/magnetic dipoles | Low-directivity probes at moderate distance |

The book's pragmatic recommendation: **always use probe correction**. §4.2 is included for academic completeness and historical perspective.

---

## 2. Measurement Without Probe Correction (§4.2)

These methods assume the probe is an ideal dipole — signal is **proportional** to a field component. Then the unknown $T_{smn}$ are extracted by orthogonality integrals from the **field** $\vec{E}$ or $\vec{H}$ sampled on the measurement sphere.

### 2.1 Setup

$$\vec{E}(A, \theta, \phi) = \frac{k}{\sqrt\eta}\sum_{s=1}^2\sum_{n=1}^N\sum_{m=-n}^n v\,T_{smn}\,\vec{F}^{(3)}_{smn}(A, \theta, \phi),\quad A > r_0 \tag{4.2}$$

$$\vec{H}(A, \theta, \phi) = -ik\sqrt\eta\sum v\,T_{smn}\,\vec{F}^{(3)}_{3-s,m,n}(A, \theta, \phi) \tag{4.3}$$

### 2.2 What you can recover from full-sphere $\vec{E}$ or $\vec{H}$ measurements (§4.2.2, Table 4.1)

Using orthogonality (A1.70) applied to the appropriate region/wave-type:

| Region (sources) | Measured | Recoverable |
|---|---|---|
| Both inner and outer have sources; both $Q^{(3)}, Q^{(4)}$ in expansion | $\vec{E}(A)$ | **All TM** ($Q^{(3)}_{2mn}$ + $Q^{(4)}_{2mn}$) |
| Same | $\vec{H}(A)$ | **All TE** ($Q^{(3)}_{1mn}$ + $Q^{(4)}_{1mn}$) |
| Inner is source-free | $\vec{E}(A)$ | All $Q^{(1)}_{smn}$ except possibly a single TE (if $j_n(kA) = 0$) |
| Same | $\vec{H}(A)$ | All $Q^{(1)}_{smn}$ except possibly a single TM |
| Outer is source-free (typical antenna) | $\vec{E}(A)$ | **All $Q^{(3)}_{smn}$** |
| Same | $\vec{H}(A)$ | **All $Q^{(3)}_{smn}$** |

### 2.3 Radial-component measurement (§4.2.3) — uses orthogonality (A1.68)

Radial electric dipole probe:

$$\boxed{\;v\,T_{2mn} = \frac{2}{\sqrt{6\pi}}\,(-1)^m\,\frac{1}{n(n+1)}\,\Bigl(\frac{kA}{h_n^{(1)}(kA)}\Bigr)^2\int_0^{2\pi}\!\!\int_0^\pi w_r^e(A, 0, \theta, \phi)\,\vec{F}^{(3)}_{2,-m,n}(A,\theta,\phi)\cdot\hat r\,\sin\theta\,d\theta\,d\phi\;} \tag{4.24}$$

Radial magnetic dipole probe (recovers TE):

$$v\,T_{1mn} = \frac{2i}{\sqrt{6\pi}}\,(-1)^m\,\frac{1}{n(n+1)}\,\Bigl(\frac{kA}{h_n^{(1)}(kA)}\Bigr)^2\int w_r^m\,\vec{F}^{(3)}_{2,-m,n}\cdot\hat r\,\sin\theta\,d\theta\,d\phi \tag{4.26}$$

> **Drawback.** Radial field components decay as $r^{-2}$ while tangential decay as $r^{-1}$. At large measurement distances, the radial measurement is dominated by tangential leakage — **probe-alignment-sensitive**.

### 2.4 Two-tangential-component measurement (§4.2.4) — uses (A1.69)

$\hat\theta$ and $\hat\phi$ dipole probes:

$$\boxed{\;v\,T_{smn} = \frac{2}{\sqrt{6\pi}}\,(-1)^m\,\{R^{(3)}_{sn}(kA)\}^{-2}\int_0^{2\pi}\!\!\int_0^\pi \{w^e_\theta\,\hat\theta + w^e_\phi\,\hat\phi\}\cdot\vec{F}^{(3)}_{s,-m,n}(A,\theta,\phi)\,\sin\theta\,d\theta\,d\phi\;} \tag{4.30}$$

Magnetic-dipole version:

$$v\,T_{smn} = \frac{2i}{\sqrt{6\pi}}\,(-1)^m\,\{R^{(3)}_{3-s,n}(kA)\}^{-2}\int \{w^m_\theta\,\hat\theta + w^m_\phi\,\hat\phi\}\cdot\vec{F}^{(3)}_{3-s,-m,n}\,\sin\theta\,d\theta\,d\phi \tag{4.31}$$

### 2.5 Wood's method (§4.2.5) — eccentric measurement sphere

Based on reciprocity integral (A1.74). Allows test-antenna minimum sphere to be offset from measurement-sphere center:

$$v\,T_{smn} = \frac{ik^2}{\sqrt{6\pi}}\,(-1)^m\int_S\Bigl\{[w_\theta^e \hat\phi - w_\phi^e \hat\theta]\cdot\vec{F}^{(4)}_{3-s,-m,n}(r,\theta,\phi) + i[w_\theta^m \hat\phi - w_\phi^m \hat\theta]\cdot\vec{F}^{(4)}_{s,-m,n}(r,\theta,\phi)\Bigr\}\,A^2\sin\theta'\,d\theta'\,d\phi' \tag{4.36, 4.37}$$

> **Requires** simultaneous measurement of two electric and two magnetic tangential components — practically replaced by Huygens-source approximations.

---

## 3. Measurement With Probe Correction — Analytical Solution (§4.3.2)

### 3.1 Starting point — transmission formula

$$\boxed{\;w(A, \chi, \theta, \phi) = \frac{v}{2}\sum_{smn,\sigma\mu\nu} T_{smn}\,e^{im\phi}\,d^n_{\mu m}(\theta)\,e^{i\mu\chi}\,C^{sn(3)}_{\sigma\mu\nu}(kA)\,R^p_{\sigma\mu\nu}\;} \tag{4.1}$$

### 3.2 Compactification via probe response constants

$$\boxed{\;P_{s\mu n}(kA) = \frac{1}{2}\sum_{\sigma\nu} C^{sn(3)}_{\sigma\mu\nu}(kA)\,R^p_{\sigma\mu\nu}\;} \tag{4.39}$$

Transmission formula becomes:

$$\boxed{\;w(A, \chi, \theta, \phi) = v\sum_{smn,\mu} T_{smn}\,e^{im\phi}\,d^n_{\mu m}(\theta)\,e^{i\mu\chi}\,P_{s\mu n}(kA)\;} \tag{4.40}$$

> **Implementation.** Precompute $P_{s\mu n}(kA)$ once per probe/distance; it doesn't depend on $(\chi, \theta, \phi)$. **Major speedup.**

### 3.3 The three-step Fourier inversion

Three orthogonalities are used in sequence:
- (4.41) $\int e^{i(m-m')\psi}\,d\psi = 2\pi\delta_{mm'}$ — twice
- (A2.10) $\int d^n_{\mu m}(\theta)\,d^{n'}_{\mu m}(\theta)\,\sin\theta\,d\theta = \frac{2}{2n+1}\delta_{nn'}$ — once

**Step 1: $\chi$ integral.** Write the transmission formula as a Fourier series in $\chi$:

$$w(A, \chi, \theta, \phi) = \sum_{\mu = -\nu_{\max}}^{\nu_{\max}} w_\mu(A, \theta, \phi)\,e^{i\mu\chi} \tag{4.42}$$

$$w_\mu(A, \theta, \phi) = v\sum_{s=1}^2 \sum_{n=1}^N \sum_{m=-n}^n T_{smn}\,e^{im\phi}\,d^n_{\mu m}(\theta)\,P_{s\mu n}(kA) \tag{4.43}$$

Inversion:

$$\boxed{\;w_\mu(A, \theta, \phi) = \frac{1}{2\pi}\int_0^{2\pi} w(A, \chi, \theta, \phi)\,e^{-i\mu\chi}\,d\chi\;} \tag{4.45}$$

**Step 2: $\phi$ integral.** Now $w_\mu(A, \theta, \phi) = \sum_{m=-N}^N w_{\mu m}(A, \theta)\,e^{im\phi}$ (4.46), with:

$$w_{\mu m}(A, \theta) = v\sum_{s=1}^2 \sum_{\substack{n=|m|\\ n \ne 0}}^N T_{smn}\,d^n_{\mu m}(\theta)\,P_{s\mu n}(kA) \tag{4.47}$$

Inversion:

$$\boxed{\;w_{\mu m}(A, \theta) = \frac{1}{2\pi}\int_0^{2\pi} w_\mu(A, \theta, \phi)\,e^{-im\phi}\,d\phi\;} \tag{4.48}$$

**Step 3: $\theta$ integral.** Now $w_{\mu m}(A, \theta) = \sum_{n = |m|, n \ne 0}^N w^n_{\mu m}(A)\,d^n_{\mu m}(\theta)$ (4.49), with:

$$w^n_{\mu m}(A) = v\sum_{s=1}^2 T_{smn}\,P_{s\mu n}(kA) \tag{4.50}$$

Inversion via (A2.10):

$$\boxed{\;w^n_{\mu m}(A) = \frac{2n+1}{2}\int_0^\pi w_{\mu m}(A, \theta)\,d^n_{\mu m}(\theta)\,\sin\theta\,d\theta\;} \tag{4.51}$$

### 3.4 Final solve for $T_{smn}$

For each $(m, n)$, write (4.50) explicitly:

$$v\,T_{1mn}\,P_{1\mu n}(kA) + v\,T_{2mn}\,P_{2\mu n}(kA) = w^n_{\mu m}(A) \tag{4.52}$$

For a **$\mu = \pm 1$ probe** (the standard case), this gives a 2×2 system:

$$\boxed{\;v\,T_{1mn}\,P_{11n}(kA) + v\,T_{2mn}\,P_{21n}(kA) = w^n_{1m}(A)\;} \tag{4.53}$$

$$\boxed{\;v\,T_{1mn}\,P_{1,-1,n}(kA) + v\,T_{2mn}\,P_{2,-1,n}(kA) = w^n_{-1, m}(A)\;} \tag{4.54}$$

Solved by Cramer's rule (or direct 2×2 inverse) for each $(m, n)$ pair.

---

## 4. Discrete Solution — Practical Implementation (§4.3.3)

### 4.1 Chi integration (§4.3.3.2)

For a $\mu = \pm 1$ probe, only $w_1$ and $w_{-1}$ are non-zero. Two samples $\chi = 0$ and $\chi = \pi/2$ suffice:

$$\boxed{\;w_{\pm 1}(A, \theta, \phi) = \frac{1}{2}\Bigl\{w(A, 0, \theta, \phi) \mp i\,w\Bigl(A, \frac{\pi}{2}, \theta, \phi\Bigr)\Bigr\}\;} \tag{4.65, 4.66}$$

(Verify: $\chi = 0$ value gives $w_1 + w_{-1}$; $\chi = \pi/2$ gives $i(w_1 - w_{-1})$; solve.)

> **Sign check.** $w_{+1}$ corresponds to **$-i$** in the bracket; $w_{-1}$ to **$+i$**. Mnemonic: $e^{i\chi}$ at $\chi = \pi/2$ gives $+i$, so to extract $w_{+1}$ subtract; to extract $w_{-1}$ add.

Alternative four-sample DFT scheme (Eq. 4.62) — slight oversampling, but exploits FFT.

### 4.2 Phi integration (§4.3.3.3)

$J_\phi$ equispaced samples in $0 \le \phi < 2\pi$, $\Delta\phi = 2\pi/J_\phi$, $J_\phi \ge 2N+1$:

$$\boxed{\;\{w_{\mu m}(A, \theta) \mid m = 0, 1, \dots, N, -N, \dots, -1\} = \text{IDFT}\{w_\mu(A, \theta, j\Delta\phi) \mid j = 0, 1, \dots, J_\phi - 1\}\;} \tag{4.68}$$

**Layout convention** (per Appendix A4): IDFT output places positive $m$ in bins $0, 1, \dots, N$ and negative $m$ in bins $J_\phi - N, \dots, J_\phi - 1$. **The IDFT here is Hansen's IDFT (numpy's FFT / $J$).**

### 4.3 Theta integration (§4.3.3.4) — the hard part

Theta integrand is **not periodic** on $[0, \pi]$. Trick: extend $w_{\mu m}(A, \theta)$ to $[0, 2\pi]$ with matching parity:

$$\tilde w_{\mu m}(A, \theta) = \begin{cases} w_{\mu m}(A, \theta), & 0 \le \theta \le \pi \\ w_{\mu m}(A, 2\pi - \theta), & \pi < \theta < 2\pi,\ (\mu - m)\text{ even} \\ -w_{\mu m}(A, 2\pi - \theta), & \pi < \theta < 2\pi,\ (\mu - m)\text{ odd}\end{cases} \tag{4.72}$$

> **Parity rule.** Each $w_{\mu m}$ sequence has parity $(\mu - m)$ under $\theta \to 2\pi - \theta$. This follows from the parity of $d^n_{\mu m}$ about $\pi$ (Appendix A2 — same parity as $\mu + m$, which matches $\mu - m$ mod 2).

Now expand into a Fourier series:

$$\tilde w_{\mu m}(A, \theta) = \sum_{l=-N}^N b_l^{\mu m}\,e^{il\theta},\quad 0 \le \theta < 2\pi \tag{4.73}$$

Compute coefficients via IDFT with $J_\theta \ge 2N + 1$ samples:

$$\boxed{\;\{b_l^{\mu m} \mid l = 0, 1, \dots, N, -N, \dots, -1\} = \text{IDFT}\{\tilde w_{\mu m}(A, j\Delta\theta) \mid j = 0, 1, \dots, J_\theta - 1\}\;} \tag{4.77}$$

### 4.4 Closed-form theta integral via $\Pi(l - m')$ (§4.3.3.5)

Substituting (4.73) and the Fourier expansion (4.70) of $d^n_{\mu m}$:

$$d^n_{\mu m}(\theta) = i^{\mu - m}\sum_{m'=-n}^n \Delta^n_{m'\mu}\,\Delta^n_{m'm}\,e^{-im'\theta} \tag{4.70}$$

into (4.51), the theta integral becomes algebraic:

$$w^n_{\mu m}(A) = \frac{2n+1}{2}\,i^{\mu-m}\sum_{l=-N}^N b_l^{\mu m}\sum_{m'=-n}^n \Delta^n_{m'\mu}\,\Delta^n_{m'm}\,G(l - m') \tag{4.75}$$

where:

$$G(l - m') = \int_0^\pi e^{i(l-m')\theta}\sin\theta\,d\theta = \begin{cases} \pm i\pi/2, & l - m' = \pm 1 \\ 0, & |l - m'| = 3, 5, 7, \dots \\ 2/[1 - (l-m')^2], & |l - m'| = 0, 2, 4, \dots\end{cases} \tag{4.76}$$

### 4.5 Parity reduction and final form (§4.3.3.5)

Using $b_l^{\mu m} = (-1)^{\mu + m}\,b_{-l}^{\mu m}$ (Eq. 4.79) and $\Delta^n_{m'\mu}\Delta^n_{m'm} = (-1)^{\mu+m}\Delta^n_{-m',\mu}\Delta^n_{-m',m}$ (Eq. 4.80, 4.81), the $l = \pm 1$ terms in $G$ **cancel**, giving:

$$\boxed{\;w^n_{\mu m}(A) = \frac{2n+1}{2}\,i^{\mu-m}\sum_{m'=-n}^n \Delta^n_{m'\mu}\,\Delta^n_{m'm}\,\sum_{l=-N}^N \Pi(l - m')\,b_l^{\mu m}\;} \tag{4.83}$$

with:

$$\Pi(l - m') = \begin{cases} 0, & (l - m')\text{ odd} \\ 2 / [1 - (l-m')^2], & (l - m')\text{ even}\end{cases} \tag{4.84}$$

> $\Pi(l - m') = \Pi(m' - l)$ — even function of its argument.

### 4.6 The $K(m')$ convolution (§4.3.3.5)

Define the inner $l$-sum:

$$K(m') = \sum_{l=-N}^N \Pi(l - m')\,b_l^{\mu m},\quad -N \le m' \le N \tag{4.85}$$

This is a **convolution**. To use FFT, extend both sequences to period $4N$:

$$\tilde\Pi(j) = \Pi(j),\quad -2N < j \le 2N,\quad \tilde\Pi(j) = \tilde\Pi(j + c\cdot 4N) \tag{4.86}$$

$$\tilde b_l^{\mu m} = \begin{cases} b_l^{\mu m}, & -N \le l \le N \\ 0, & -2N < l < -N\text{ and }N < l \le 2N\end{cases},\quad \tilde b_l^{\mu m} = \tilde b_{l + c\cdot 4N}^{\mu m} \tag{4.87}$$

Then:

$$K(m') = \sum_{l=0}^{4N-1} \tilde\Pi(l - m')\,\tilde b_l^{\mu m} \tag{4.88}$$

Computed via DFT:

$$\boxed{\;K(m') = \text{IDFT}\{\text{DFT}\{\tilde\Pi(i)\}\cdot \text{DFT}\{\tilde b_j^{\mu m}\}\}\;} \tag{4.89}$$

(Both sequences of length $4N$. $\text{DFT}\{\tilde\Pi\}$ is precomputed **once** per $N$.)

**Parity exploitation:** $\tilde b_l^{\mu m} = (-1)^{\mu + m}\,\tilde b_{-l}^{\mu m}$ (Eq. 4.90) ⟹ $K(m') = (-1)^{\mu+m}\,K(-m')$ (Eq. 4.91). Halves storage and work.

### 4.7 Final $m'$ summation

$$w^n_{\mu m}(A) = \frac{2n+1}{2}\,i^{\mu - m}\sum_{m'=-n}^n \Delta^n_{m'\mu}\,\Delta^n_{m'm}\,K(m'),\quad \mu = \pm 1 \tag{4.92}$$

Using $\Delta^n_{-m',\mu}\Delta^n_{-m',m}\,K(-m') = \Delta^n_{m'\mu}\Delta^n_{m'm}\,K(m')$ (Eq. 4.93), reduce the sum to $0 \le m' \le n$ (multiply by 2 for $m' > 0$).

### 4.8 Delta recurrence

$$\sqrt{(n+m'+1)(n-m')}\,\Delta^n_{m'+1,m} + \sqrt{(n+m')(n-m'+1)}\,\Delta^n_{m'-1, m} + 2m\,\Delta^n_{m'm} = 0 \tag{4.94}$$

**Run backwards from $m' = n$ to $m' = 0$** (stable; see Appendix A2 §5). Seed value:

$$\Delta^n_{nm} = 2^{-n}\sqrt{\frac{2n(2n-1)\cdots(n-m+1)}{(n-m)!}} \tag{4.95}$$

Compute $\Delta^n_{m'\mu}$ and $\Delta^n_{m'm}$ **simultaneously** (same square roots) — saves ~half the work.

---

## 5. Field Transformations and Antenna Parameters (§4.3.4)

### 5.1 Input/output probe concept (§4.3.4.1)

Once $T_{smn}$ is recovered, the transmission formula (4.40) re-applied with a new probe and distance $A'$ computes the field on any other sphere:

$$w'(A', \chi, \theta, \phi) = \sum_{smn,\mu=\pm 1} v\,T_{smn}\,e^{im\phi}\,d^n_{\mu m}(\theta)\,e^{i\mu\chi}\,P'_{s\mu n}(kA') \tag{4.97}$$

| Probe role | Symbol | Purpose |
|---|---|---|
| Input probe | $P_{s\mu n}(kA)$ | Used in inversion (recovers $T_{smn}$) |
| Output probe | $P'_{s\mu n}(kA')$ | Used in synthesis (evaluates field at $A'$) |

### 5.2 Far-field limits (§4.3.4.2)

For $kA \to \infty$, the translation coefficient asymptotics (A3.22–A3.24) give:

$$P_{s,\pm 1, n}(kA) \to \frac{e^{ikA}}{kA}\,P^\infty_{s,\pm 1, n},\quad kA \to \infty \tag{4.100}$$

**General linearly-polarized $\hat x'$ probe:**

$$P^\infty_{s1n} = \tfrac{1}{4}\sqrt{2n+1}\,i^{-n-1}\sum_{\nu=1}^{\nu_{\max}}\sqrt{2\nu+1}\,i^\nu\,\{R^p_{11\nu} + R^p_{21\nu}\} \tag{4.101}$$

$$P^\infty_{s,-1,n} = (-1)^{s+1}\,P^\infty_{s1n} \tag{4.102}$$

**Closed forms for electric and magnetic Hertzian-dipole probes:**

| Probe | $P^\infty_{s1n}$ | $P^\infty_{s,-1,n}$ |
|---|---|---|
| $\hat x'$ electric dipole | $-\tfrac{\sqrt 6}{8}\sqrt{2n+1}\,i^{-n}$ (4.103) | $(-1)^s\,\tfrac{\sqrt 6}{8}\sqrt{2n+1}\,i^{-n}$ (4.104) |
| $\hat x'$ magnetic dipole (note: $\hat y'$-polarized!) | $\tfrac{\sqrt 6}{8}\sqrt{2n+1}\,i^{-n+1}$ (4.105) | $(-1)^s\,\tfrac{\sqrt 6}{8}\sqrt{2n+1}\,i^{-n+1}$ (4.106) |

**Normalized far-field signal:**

$$W(\chi, \theta, \phi) = \lim_{kA\to\infty}\Bigl[w(A, \chi, \theta, \phi)\,\frac{kA}{e^{ikA}}\Bigr] \tag{4.107}$$

Far-field transmission formula:

$$W(\chi, \theta, \phi) = \sum_{smn, \mu=\pm 1} v\,T_{smn}\,e^{im\phi}\,d^n_{\mu m}(\theta)\,e^{i\mu\chi}\,P^\infty_{s\mu n} \tag{4.108}$$

$$W'(\chi, \theta, \phi) = \sum_{smn, \mu=\pm 1} v\,T_{smn}\,e^{im\phi}\,d^n_{\mu m}(\theta)\,e^{i\mu\chi}\,P'^\infty_{s\mu n} \tag{4.109}$$

### 5.3 Relative measurements (§4.3.4.3) — propagation of unknown constants

| Input | Probe constants | Algorithm produces |
|---|---|---|
| Absolute $w$ | Absolute $P$ | $v\,T_{smn}$, $w'$ |
| Relative ($c_1 w$) | Absolute ($c_2 = 1$) | $c_1\,v T_{smn}$, $c_1\,w'$ |
| Absolute | Relative ($c_2$) | $v T_{smn}/c_2$, $w'/c_2$ |
| Both relative | | $vT_{smn}\,c_1/c_2$, $w'\,c_1/c_2 \cdot c_3$ |

(With $c_3$ from unknown output-probe scaling.)

### 5.4 Antenna parameter extraction (§4.3.4.4)

**Directivity (relative measurement, output = far-field $\hat x'$ electric dipole):**

$$\boxed{\;D_t(\theta, \phi) = \frac{8}{3}\cdot\frac{|W'\,c_1/c_2|^2}{\sum_{smn}|v T_{smn}\,c_1/c_2|^2}\;} \tag{4.112}$$

The unknown $c_1/c_2$ **cancels** — directivity is recoverable from purely relative measurements.

**Gain (requires measurement of $c_1 v$ separately):**

$$\boxed{\;G_t = \frac{8}{3}\cdot\frac{|c_1 W'|^2}{|c_1 v|^2}\;} \tag{4.115}$$

**EIRP:**

$$\text{EIRP}(\theta, \phi) = \tfrac{1}{2}|v|^2\,G_t(\theta, \phi) \tag{4.116}$$

**Absolute far-field pattern:**

$$\boxed{\;\vec K(\theta, \phi) = \frac{2\sqrt 6}{3v}\bigl\{W'(0, \theta, \phi)\,\hat\theta + W'(\pi/2, \theta, \phi)\,\hat\phi\bigr\}\;} \tag{4.117}$$

With a general directive output probe of gain $G_p$:

$$\vec K(\theta, \phi) = \frac{2\sqrt 6}{3v}\bigl\{W'(0)\,\hat\theta + W'(\pi/2)\,\hat\phi\bigr\}\sqrt{G_e/G_p},\quad G_e = 3/2 \tag{4.118}$$

**Relative far-field pattern** (no $v$ knowledge needed):

$$\vec K_\text{rel}(\theta, \phi) = c\,\{W'(0, \theta, \phi)\,\hat\theta + W'(\pi/2, \theta, \phi)\,\hat\phi\} \tag{4.119}$$

**Tangential field at radius $A'$:**

$$\vec E_\text{tang}(A', \theta, \phi) = \frac{2k}{\sqrt{6\pi\eta}}\bigl\{w'^e(A', 0)\,\hat\theta + w'^e(A', \pi/2)\,\hat\phi\bigr\} \tag{4.120}$$

$$\vec H_\text{tang}(A', \theta, \phi) = \frac{2k\sqrt\eta}{\sqrt{6\pi}}\bigl\{w'^m(A', 0)\,\hat\theta + w'^m(A', \pi/2)\,\hat\phi\bigr\} \tag{4.121}$$

---

## 6. SNIFT Algorithm — Block-by-Block (§4.4)

### 6.1 Four modes of operation

| Mode | Input | Output | Use case |
|---|---|---|---|
| (1) $w \to w'$ | Near-field samples | Near-field at $A'$ | Field at different radius |
| (2) $w \to W'$ | Near-field samples | Far-field pattern | **Standard NF-to-FF** |
| (3) $W \to w'$ | Far-field samples | Near-field at $A'$ | FF-to-NF |
| (4) $W \to W'$ | Far-field samples | Far-field pattern | Coordinate/polarization filtering |

### 6.2 First part — recovery of $T_{smn}$ (Steps 1–7)

Pipeline:

| Step | Operation | Equation |
|---|---|---|
| 1 | $\chi$ integration: $w_{\pm 1}(A, \theta, \phi) = \tfrac{1}{2}\{w(A, 0, \theta, \phi) \mp i\,w(A, \pi/2, \theta, \phi)\}$ | (4.126) |
| 2 | $\phi$ IDFT: $\{w_{\pm 1, m}(A, \theta)\} = \text{IDFT}\{w_{\pm 1}(A, \theta, j\Delta\phi)\}$ | (4.127) |
| 3 | Extend data $\tilde w_{\mu m}$ per (4.72); $\theta$ IDFT: $\{b_l^{\pm 1, m}\} = \text{IDFT}\{\tilde w_{\pm 1, m}(A, j\Delta\theta)\}$ | (4.128) |
| 4 | Convolve via FFT: $K(m') = \text{IDFT}\{\text{DFT}\{\tilde\Pi\}\cdot\text{DFT}\{\tilde b^{\pm 1, m}_j\}\}$ | (4.131) |
| 5 | For each $n$: compute $w^n_{\mu m}(A) = \tfrac{2n+1}{2}\,i^{\mu-m}\sum_{m'=-n}^n \Delta^n_{m'\mu}\Delta^n_{m'm}\,K(m')$ for $\mu = \pm 1$ | (4.132) |
| 6 | Recurrence-compute $\Delta^n_{m'\mu}$, $\Delta^n_{m'm}$ for $0 \le m' \le n$ (using 4.94, 4.95) | — |
| 7 | Solve 2×2 system: $T_{1mn}, T_{2mn}$ from (4.53, 4.54) | (4.133, 4.134) |

### 6.3 Second part — evaluation of $w'$ (Steps 8–11)

Reverse pipeline using already-computed delta products:

$$w'(A', \chi, \theta, \phi) = \sum_{\mu = \pm 1} e^{i\mu\chi}\sum_{m=-N}^N e^{im\phi}\sum_{m'=-N}^N e^{im'\theta}\sum_{n=\max(|m'|, |m|, 1)}^N \Delta^n_{m'\mu}\Delta^n_{m'm}\,i^{m-\mu}\sum_{s=1}^2 v\,T_{smn}\,P'_{s\mu n}(kA') \tag{4.135}$$

| Step | Operation |
|---|---|
| 8 | For each $n$: accumulate $\Delta^n_{m'\mu}\Delta^n_{m'm}\,i^{m-\mu}\sum_s v T_{smn} P'_{s\mu n}$ into the $n$-sum array |
| 9 | $m'$ summation via DFT (yields $\theta$-grid output) |
| 10 | $m$ summation via DFT (yields $\phi$-grid output) |
| 11 | $\mu$ summation (yields final $\chi$-dependence: $\chi = 0, \pi/2, \dots$) |

### 6.4 Input parameters (Table 4.4)

| Parameter | Meaning |
|---|---|
| $J_\theta$ | # samples in $0 \le \theta < 2\pi$ of input |
| $J_\phi$ | # samples in $0 \le \phi < 2\pi$ of input |
| $A/\lambda$ | Input sphere radius (skipped if input is $W$) |
| $\nu_{\max}$ | Max $\nu$ for input probe |
| $N$ | Max $n$ for test antenna ($1 \le N \le (J_\theta-1)/2$) |
| $M$ | Max $\|m\|$ (default $N$; $0 \le M \le \min((J_\phi-1)/2, N)$) |
| $J'_\theta, J'_\phi, A'/\lambda, \nu'_{\max}$ | Same for output |
| Mode | (1), (2), (3), or (4) |

### 6.5 Sampling counts (§4.4.3)

| Sampling strategy | # samples $N_s$ on full sphere |
|---|---|
| Uniform $(2N+1)$ × $(2N+1)$, 2 χ values | $N_1 = 2(2N+1)(N+1)$ — Eq. (4.137) |
| Thinned (sin θ scaling near poles) | $N_2 = (4/\pi)(2N+1)(N+1)$ — Eq. (4.138) |
| Cylindrical-symmetric antenna ($M = kr_c + 10$) | $N_3 = 2(2M+1)(N+1)$ — Eq. (4.139) |

Total modes $N_0 = 2N(N+2)$ (Eq. 4.136). Uniform sampling is ~50% efficient; thinning yields ~$N_0$.

---

## 7. Worked Example (§4.4.5) — Critical Regression Test

**Input:** $W(\chi, \theta, \phi) = (20\cos 2\theta + 32\cos\theta + 12)\cos(\chi + \phi)$, band-limited with $N = 2, M = 1$.

**Setup:**
- $\theta$ samples at $0°, 60°, 120°, 180°$ (4 samples in natural range; 2 generated to fill $0°$–$360°$)
- $\phi$ samples at $0°, 90°, 180°, 270°$
- $\chi$ samples at $0°, 90°$
- Input/output probes: $\hat x'$-polarized maximum-directivity antenna with $N = 2$
- All radii infinite — only $P^\infty$ and $P'^\infty$ used (no Hankel functions)

**Probe response constants** (Eq. 4.151):

$$P^\infty_{111} = P^\infty_{211} = P^\infty_{1,-1,1} = -\frac{\sqrt 6}{2}i,\quad P^\infty_{2,-1,1} = \frac{\sqrt 6}{2}i$$

$$P^\infty_{112} = P^\infty_{212} = P^\infty_{1,-1,2} = -\frac{\sqrt{10}}{2},\quad P^\infty_{2,-1,2} = \frac{\sqrt{10}}{2}$$

**Delta products** (Eq. 4.152):

$$\Delta^1_{01}\Delta^1_{01} = \tfrac{1}{2},\ \Delta^2_{01}\Delta^2_{01} = 0$$

$$\Delta^1_{11}\Delta^1_{11} = \tfrac{1}{4},\ \Delta^2_{11}\Delta^2_{11} = \tfrac{1}{4}$$

$$\Delta^2_{21}\Delta^2_{21} = \tfrac{1}{4}$$

**Step 1 — χ integration** (Table 4.5b — partial):

| | $\phi=0°, \mu=+1$ | $\phi=0°, \mu=-1$ |
|---|---|---|
| $\theta = 0°$ | 32 | 32 |
| $\theta = 60°$ | 9 | 9 |
| $\theta = 120°$ | $-7$ | $-7$ |
| $\theta = 180°$ | 0 | 0 |

**Step 2 — φ integration** (Table 4.5c) yields non-zero only at $m = \pm 1$:

| $\theta$ | $W_{11}(\theta)$ | $W_{-1,-1}(\theta)$ |
|---|---|---|
| $0°$ | 32 | 32 |
| $60°$ | 9 | 9 |
| $120°$ | $-7$ | $-7$ |
| $180°$ | 0 | 0 |

**Step 3 — extend and IDFT** (with parity reasoning):

$$\tilde W_{11} = \{32, 9, -7, 0, -7, 9\}\quad\Rightarrow\quad b_l^{11} = \{6, 8, 5, 0, 5, 8\} \tag{4.143}$$

Zero-padded to length 12:

$$\tilde b_l^{11} = \{6, 8, 5, 0, 0, 0, 0, 0, 0, 0, 5, 8\} \tag{4.144}$$

**Step 12 — Π sequence (length 12)**:

$$\tilde\Pi = \{2, 0, -\tfrac{2}{3}, 0, -\tfrac{2}{15}, 0, -\tfrac{2}{35}, 0, -\tfrac{2}{15}, 0, -\tfrac{2}{3}, 0\} \tag{4.145}$$

**Step 4 — Convolution via FFT:**

$$K(m') \to \{K(0), K(1), K(2)\} = \{\tfrac{16}{3}, \tfrac{64}{3}, \tfrac{32}{3}\}\quad\text{(after doubling for } m' > 0\text{)} \tag{4.150}$$

**Step 5–7 — solve for $T_{smn}$** at each $n$:

| $n$ | Equation system | Solution |
|---|---|---|
| 1 | $\{T_{111}, T_{211}\} = \{12, 0\}/(- \sqrt 6/2\cdot i)$-matrix | $v T_{111} = v T_{211} = v T_{1,-1,1} = -v T_{2,-1,1} = (\sqrt 6/2)i$ |
| 2 | $\{T_{112}, T_{212}\} = \{20, 0\}/\dots$ | $v T_{112} = v T_{212} = v T_{1,-1,2} = -v T_{2,-1,2} = -2\sqrt{10}$ |

**Step 8–11 — back-substitution to compute $W'$**:

The final output table (Table 4.7c) reproduces the input $W$ exactly:

| | $\phi=0°, \chi=0$ | $\phi=0°, \chi=\pi/2$ |
|---|---|---|
| $\theta = 0°$ | 64 | 0 |
| $\theta = 60°$ | 18 | 0 |
| $\theta = 120°$ | $-14$ | 0 |
| $\theta = 180°$ | 0 | 0 |

> **Validation power.** Reproducing every table value (4.5a-c, 4.6a-c, 4.7a-c) of the worked example is the **gold-standard regression test** for a Python implementation. Every intermediate quantity is integer or simple radical — exact arithmetic is possible.

---

## 8. Reference Test Case Errors (§4.4.6)

12-element dipole array antenna (50λ diameter), 2-element probe, $A = 50\lambda$:

| $J$ (samples/great circle) | 100% error | 90% error | 50% error | 10% error |
|---|---|---|---|---|
| 260 | ~0 dB | $-15$ dB | $-23$ dB | $-25$ dB |
| 300 | $-20$ | $-30$ | $-35$ | $-40$ |
| 350 | $-40$ | $-65$ | $-75$ | $-85$ |
| 360 | $-50$ | $-90$ | $-105$ | $-115$ |
| 400 | $-50$ | $-90$ | $-105$ | $-115$ |

> **Conclusion.** With $J \ge \lceil kr_t \rceil + 10$ samples per great circle, sampling error vanishes; remaining errors are floating-point rounding (in single precision, ~$-50$ dB; in double, ~$-100$ dB+).

---

## 9. Algorithm Variants — Review (§4.5)

| Method | Reference | Notes |
|---|---|---|
| Simpson rule on θ, φ | James & Longdon | Only practical for $N \le 3$ |
| Trapezoidal θ + Ludwig iteration | Ludwig | Iterative refinement; corrects rounding |
| Simpson + Wood approximations | Wood | Eccentric measurements via Huygens probes |
| Gaussian quadrature | — | Exact, but non-equidistant samples |
| Ricardi–Burrows recursion | — | Fourier extension + recursion |
| Wacker | [1, 2] | Fourier expand $d^n_{\mu m}\sin\theta$ jointly; key insight |
| Wacker–Larsen | [3, 4] | Keep $\sin\theta$ in integral; reuse $\Delta$'s for synthesis; convolution form for $K(m')$ |
| Yaghjian–Wittmann | [20, 21, 22] | Tangential-component formulation; halves $\theta$-integral terms |
| Lewis hemispherical | [12, 23] | Halves work when back hemisphere is zero |

The **Wacker–Larsen** method (§4.3) is what SNIFT implements and what this document treats as canonical.

---

## 10. Implementation Checklist for Python Code

### Setup
1. **Truncation $N$**: choose $N = \lfloor kr_0 \rfloor + 10$ where $r_0$ is the test-antenna minimum-sphere radius (Eq. 2.31).
2. **Sample counts**: $J_\theta \ge 2N + 1$ for theta in $[0, 2\pi)$ extended, $J_\phi \ge 2N + 1$ for phi, and 2 chi values ($\chi = 0, \pi/2$).
3. **Equispaced grids**: $\theta_j = j\cdot 2\pi/J_\theta$, $\phi_j = j\cdot 2\pi/J_\phi$. **Measurements are taken only at $\theta_j \in [0, \pi]$** — the second half is computer-generated via (4.72).

### Probe response constants
4. **Precompute $P_{s\mu n}(kA)$** once via (4.39); use sympy or recurrence-based 3-j computation per Appendix A3.
5. **Symmetry check**: $P_{s,-1,n} = (-1)^{s+1}\,P_{s,1,n}$ (Eq. 3.27).
6. **Far-field limit**: $|P_{s,\pm 1, n}(kA)\cdot kA/e^{ikA}|$ stable as $kA \to \infty$, matching $P^\infty_{s,\pm 1, n}$ from (4.101)–(4.106).

### Step 1: χ integration
7. **Verify (4.65)/(4.66)**: with synthetic $w(\chi) = c_+ e^{i\chi} + c_- e^{-i\chi}$, recovery should give $w_1 = c_+$, $w_{-1} = c_-$.

### Step 2: φ IDFT
8. **Layout convention**: Hansen-IDFT output places $m = 0, 1, \dots, N$ in bins 0..N and $m = -N, \dots, -1$ in bins $J_\phi - N, \dots, J_\phi - 1$. (Middle bins should be zero for band-limited input.)
9. **Hansen-vs-numpy mapping**: Implement `hansen_idft(x) = np.fft.fft(x) / len(x)` per Appendix A4 §2.

### Step 3: θ extension and IDFT
10. **Parity rule (4.72)**: extension parity is $(\mu - m)$. Test for each pair: $(\mu = +1, m = 0) \to$ odd; $(\mu = +1, m = 1) \to$ even; etc.
11. **Verify extension preserves period**: $\tilde w_{\mu m}(0) = \tilde w_{\mu m}(2\pi)$ (after periodic wrap).
12. **IDFT layout**: same as step 8 but for $l$ bins.

### Step 4: K(m') convolution
13. **$\tilde\Pi$ table**: precompute for the chosen $N$. Verify $\tilde\Pi(0) = 2$, $\tilde\Pi(\pm 1) = 0$, $\tilde\Pi(\pm 2) = -2/3$, $\tilde\Pi(\pm 4) = -2/15$, etc.
14. **Parity check after convolution**: $K(-m') = (-1)^{\mu+m}\,K(m')$ (Eq. 4.91).
15. **Length-$4N$ zero-padding**: assert that bins $\{-2N..-N-1\} \cup \{N+1..2N\}$ of $\tilde b$ are zero.

### Steps 5–6: Delta recurrence
16. **Seed $\Delta^n_{nm}$** with (4.95). Verify $\Delta^n_{nn} = 2^{-n}$ (Appendix A2 Eq. A2.42).
17. **Backwards recurrence (4.94)**: stable from $m' = n$ down to $m' = 0$. Compare to Appendix A2 §6 explicit tables for $n \le 5$.
18. **Simultaneous computation**: $\Delta^n_{m'\mu}$ and $\Delta^n_{m'm}$ recur with the same square-root coefficients — compute both in a single loop.

### Step 7: 2×2 solve
19. **2×2 system**: condition number should be $O(1)$ for typical probes; flag if $|P_{11n}P_{2,-1,n} - P_{21n}P_{1,-1,n}|$ is small.

### Synthesis (Steps 8–11)
20. **Delta-product reuse**: same $\Delta^n_{m'\mu}\Delta^n_{m'm}$ are used in both analysis and synthesis — don't recompute. Cache during the n-loop.
21. **DFT vs IDFT consistency**: analysis uses IDFT (samples → coefficients); synthesis uses DFT (coefficients → samples). Wrap both around `hansen_idft`/`hansen_dft` per A4 §2.

### End-to-end validation
22. **Worked-example regression test (§7)**: implement the explicit $W = (20\cos 2\theta + 32\cos\theta + 12)\cos(\chi + \phi)$ input. Verify every table in §4.4.5 reproduces exactly (integer/radical arithmetic permits zero tolerance).
23. **Round-trip identity**: for a synthetic $T_{smn}$ set, generate $w$ via (4.40), invert to recover $T'_{smn}$, assert $\|T' - T\|/\|T\| < 10^{-12}$ (double precision).
24. **Friis far-field check**: at $kA = 10000$ with an electric-dipole probe, the recovered $T_{smn}$ should produce a far field matching the Friis formula $|w|^2/|v|^2 = G_p G_t/(4(kA)^2)$ to 3+ significant figures.
25. **Hertzian-dipole regression**: with a $\hat z$-electric dipole as test antenna ($T_4 = 1$ only), all other recovered $T_{smn}$ must be $0$ to floating-point precision; the recovered $T_4$ must equal 1.

### Sampling and convergence
26. **Truncation N**: scan $J_\theta = J_\phi$ from $2N - 10$ up to $4N$ and plot error. Should plateau at $J \approx 2N + 1$ (Eq. 4.137) and become floating-point-limited beyond.
27. **Antenna parameters**: directivity, gain, EIRP, far-field patterns recovered per (4.112)–(4.118) — closed-form sanity check for known antennas (dipole: $D = 1.5$; Huygens: $D = 3$).

### Performance
28. **FFT length**: $4N$ for $K(m')$ convolution (Eq. 4.89). Choose $4N$ as a product of small primes; pad to $4N' \ge 4N$ if needed.
29. **Half-storage**: exploit $K(-m') = \pm K(m')$ (Eq. 4.91) to halve memory.
30. **Parity packing**: pair $m$ even/odd sequences through Fourier ops (Lewis trick, [12]) — cuts FFT count in half.

---

*End of Chapter 4 reference.*
