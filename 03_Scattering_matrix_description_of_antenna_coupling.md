# Chapter 3 — Scattering Matrix Description of Antenna Coupling

**Source:** Hansen, *Spherical Near-Field Antenna Measurements*, Chapter 3 (pp. 61–88).

**Cross-references:**
- Symbols: [A5_List_of_principal_symbols_and_uses.md](A5_List_of_principal_symbols_and_uses.md)
- Chapter 2 (Scattering matrix theory): [02_Scattering_matrix_description_of_an_antenna.md](02_Scattering_matrix_description_of_an_antenna.md)
- Appendix A1 (Wave functions, $\vec{K}_{smn}$): [A1_Spherical_wave_functions_notation_and_properties.md](A1_Spherical_wave_functions_notation_and_properties.md)
- Appendix A2 (Rotation, $d^n_{\mu m}$, Euler angles): [A2_Rotation_of_spherical_waves.md](A2_Rotation_of_spherical_waves.md)
- Appendix A3 (Translation, $C^{sn(c)}_{\sigma\mu\nu}$): [A3_Translation_of_spherical_waves.md](A3_Translation_of_spherical_waves.md)

> **Purpose.** Computational reference for the **spherical transmission formula** — the central equation of probe-corrected spherical near-field measurements. Every formula is tagged by book number `(3.N)`. The transmission formula (3.10) and its companions (3.17), (3.20) are the boxed equations that any Python implementation must reproduce exactly.

---

## 1. Geometry (§3.2.1)

| Frame | Origin | Coordinates | Holds |
|---|---|---|---|
| Unprimed $(x, y, z)$ | Test antenna | $(r, \theta, \phi)$ | Test antenna fixed; minimum sphere radius $r_0$ |
| Primed $(x', y', z')$ | Probe | $(r', \theta', \phi')$ | Probe pointing **at the unprimed origin**; minimum sphere radius $r'_0$ |
| Probe position in unprimed | — | $(A, \theta, \phi)$ | Probe origin at distance $A$ from test-antenna origin, in direction $(\theta, \phi)$ |
| Probe roll | — | $\chi$ | Rotation of probe about its own $z'$ axis |

**Probe minimum sphere does not intersect test-antenna minimum sphere:** $A > r_0 + r'_0$.

Typical scan values: $\chi \in \{0, \pi/2\}$ (two roll positions sufficient for full polarization recovery).

---

## 2. Test Antenna Transmitting, Probe Receiving (§3.2.2)

### Starting field (test antenna only)

$$
\vec{E}_t(r,\theta,\phi) = \frac{k}{\sqrt{\eta}}\sum_{smn} Q^{(3)}_{smn}\,\vec{F}^{(3)}_{smn}(r,\theta,\phi) = \frac{k}{\sqrt{\eta}}\sum_{smn} v\,T_{smn}\,\vec{F}^{(3)}_{smn}(r,\theta,\phi),\quad r > r_0 \tag{3.1}
$$

(Using $Q^{(3)}_{smn} = b_{smn} = v\,T_{smn}$ from Eq. 2.66.)

### The four-step coordinate cascade (test-antenna frame → probe frame)

Each spherical mode $\vec{F}^{(3)}_{smn}$ is re-expressed in the probe frame by:

| Step | Transform | Mathematical action |
|---|---|---|
| 1 | Rotate frame 1 by $\phi_o$ about $\hat z$ | Phase $e^{im\phi_o}$ — Eq. (3.2) |
| 2 | Rotate frame 2 by $\theta_o$ about $\hat y_1$ | Mix over $\mu$ via rotation coefficient $d^n_{\mu m}(\theta_o)$ — Eq. (3.3) |
| 3 | Rotate frame 3 by $\chi_o$ about $\hat z_2$ | Phase $e^{i\mu\chi_o}$ — Eq. (3.4) |
| 4 | Translate by $A$ along $\hat z_3$ | Mix over $(\sigma, \nu)$ via $C^{sn(3)}_{\sigma\mu\nu}(kA)$, with $\vec{F}^{(3)} \to \tfrac{1}{2}(\vec{F}^{(3)} + \vec{F}^{(4)})$ — Eq. (3.5) |

> **Why the $\tfrac{1}{2}(\vec{F}^{(3)} + \vec{F}^{(4)})$ split?** After axial translation, the original outgoing wave $\vec{F}^{(3)}$ from the test-antenna origin appears in the probe frame as a standing wave $\vec{F}^{(1)} = \tfrac{1}{2}(\vec{F}^{(3)} + \vec{F}^{(4)})$ inside the sphere of radius $A - r_0$ around the probe origin (where the probe is located). The probe sees the $\vec{F}^{(4)}$ component as the **incoming wave**.

### Composite transform formula

$$
\boxed{\;\vec{F}^{(3)}_{smn}(r,\theta,\phi) = \sum_{\sigma\mu\nu} e^{im\phi_o}\,d^n_{\mu m}(\theta_o)\,e^{i\mu\chi_o}\,C^{sn(3)}_{\sigma\mu\nu}(kA)\,\tfrac{1}{2}\bigl\{\vec{F}^{(3)}_{\sigma\mu\nu}(r',\theta',\phi') + \vec{F}^{(4)}_{\sigma\mu\nu}(r',\theta',\phi')\bigr\}\;} \tag{3.6}
$$

### Field expansion in the probe frame

$$
\vec{E}_t = \frac{k}{\sqrt{\eta}}\sum_{\substack{smn\\\sigma\mu\nu}} v\,T_{smn}\,e^{im\phi_o}\,d^n_{\mu m}(\theta_o)\,e^{i\mu\chi_o}\,C^{sn(3)}_{\sigma\mu\nu}(kA)\,\tfrac{1}{2}\{\vec{F}^{(3)}_{\sigma\mu\nu} + \vec{F}^{(4)}_{\sigma\mu\nu}\} \tag{3.7}
$$

### Incoming-wave coefficients at the probe

The probe sees $\vec{F}^{(4)}_{\sigma\mu\nu}$ as incoming with amplitude:

$$
a_{\sigma\mu\nu} = \frac{v}{2}\sum_{smn} T_{smn}\,e^{im\phi_o}\,d^n_{\mu m}(\theta_o)\,e^{i\mu\chi_o}\,C^{sn(3)}_{\sigma\mu\nu}(kA) \tag{3.9}
$$

### Probe-received signal (matched load: $\Gamma_l = 0$)

Using $w = \mathbf{R}^p\,\mathbf{a}$ (Eq. 2.71 with $\Gamma_l = 0$) and dropping the $o$-subscript on probe angles:

$$
\boxed{\;w(A, \chi, \theta, \phi) = \frac{v}{2}\sum_{\substack{smn\\\sigma\mu\nu}} T_{smn}\,e^{im\phi}\,d^n_{\mu m}(\theta)\,e^{i\mu\chi}\,C^{sn(3)}_{\sigma\mu\nu}(kA)\,R^p_{\sigma\mu\nu}\;} \tag{3.10}
$$

> **The Transmission Formula.** This is the central equation of the entire book. It expresses the complex probe-received signal as a function of:
> - **Inputs (knowns):** probe coordinates $(A, \chi, \theta, \phi)$, probe receiving coefficients $R^p_{\sigma\mu\nu}$, source amplitude $v$.
> - **Unknowns:** test-antenna transmitting coefficients $T_{smn}$.
>
> The role of Chapter 4 is to **invert (3.10)** to recover $T_{smn}$ from sampled $w$-data. **No reciprocity assumed for either antenna.**

### Assumptions used in deriving (3.10)

1. **No multiple reflections** between probe and test antenna (probe-scattered field does not perturb the test antenna). Relaxed in §9 below.
2. **Matched load** at the probe ($\Gamma_l = 0$). Relaxed by the gain-measurement discussion in Chapter 5.

---

## 3. Test Antenna Receiving, Probe Transmitting (§3.2.3)

### Probe-radiated field in the probe frame

$$
\vec{E}_p(r',\theta',\phi') = \frac{k}{\sqrt{\eta}}\sum_{\sigma\mu\nu} v_p\,T^p_{\sigma\mu\nu}\,\vec{F}^{(3)}_{\sigma\mu\nu}(r',\theta',\phi'),\quad r' > r'_0 \tag{3.11}
$$

### Inverse coordinate cascade

Apply the same four-step cascade in **reverse order with negated angles** $(-A, -\chi_o, -\theta_o, -\phi_o)$. Result expressing $\vec{E}_p$ in the test-antenna frame:

$$
\vec{E}_p = \frac{k}{\sqrt{\eta}}\sum_{\substack{\sigma\mu\nu\\smn}} v_p\,T^p_{\sigma\mu\nu}\,C^{\sigma\nu(3)}_{s\mu n}(-kA)\,e^{-i\mu\chi_o}\,d^n_{m\mu}(-\theta_o)\,e^{-im\phi_o}\,\tfrac{1}{2}\{\vec{F}^{(3)}_{smn} + \vec{F}^{(4)}_{smn}\} \tag{3.12}
$$

### Test antenna received signal (matched load: $\Gamma = 0$)

Initial form:

$$
a_{smn} = \frac{v_p}{2}\sum_{\sigma\mu\nu} T^p_{\sigma\mu\nu}\,C^{\sigma\nu(3)}_{s\mu n}(-kA)\,e^{-i\mu\chi_o}\,d^n_{m\mu}(-\theta_o)\,e^{-im\phi_o} \tag{3.14}
$$

$$
w_t = \frac{v_p}{2}\sum_{\substack{\sigma\mu\nu\\smn}} T^p_{\sigma\mu\nu}\,C^{\sigma\nu(3)}_{s\mu n}(-kA)\,e^{-i\mu\chi}\,d^n_{m\mu}(-\theta)\,e^{-im\phi}\,R_{smn} \tag{3.15}
$$

### Simplification via symmetry relations (A3.13) and (A2.7)

$$
w_t = \frac{v_p}{2}\sum_{\substack{\sigma\mu\nu\\smn}} T^p_{\sigma\mu\nu}\,C^{sn(3)}_{\sigma,-\mu,\nu}(kA)\,e^{-i\mu\chi}\,d^n_{\mu m}(\theta)\,e^{-im\phi}\,R_{smn} \tag{3.16}
$$

### Final form via (A2.9) and relabeling

$$
\boxed{\;w_t = \frac{v_p}{2}\sum_{\substack{\sigma\mu\nu\\smn}}(-1)^\mu\,T^p_{\sigma,-\mu,\nu}\,C^{sn(3)}_{\sigma\mu\nu}(kA)\,e^{i\mu\chi}\,d^n_{\mu m}(\theta)\,e^{im\phi}\,(-1)^m\,R_{s,-m,n}\;} \tag{3.17}
$$

---

## 4. Adjoint Substitution — Recasting (3.17) (§3.2.3)

Using the adjoint relations (Eqs. 2.103, 2.104) between original and adjoint antenna coefficients:

$$
T^p_{\sigma,-\mu,\nu} = (-1)^\mu\,R^{p'}_{\sigma\mu\nu} \tag{3.18}
$$

$$
R_{s,-m,n} = (-1)^m\,T'_{smn} \tag{3.19}
$$

The receiving formula (3.17) becomes:

$$
\boxed{\;w_t = \frac{v_p}{2}\sum_{\substack{\sigma\mu\nu\\smn}} R^{p'}_{\sigma\mu\nu}\,C^{sn(3)}_{\sigma\mu\nu}(kA)\,e^{i\mu\chi}\,d^n_{\mu m}(\theta)\,e^{im\phi}\,T'_{smn}\;} \tag{3.20}
$$

> **Key observation.** Equation (3.20) has the **same structural form** as the transmitting-case formula (3.10) — only with primes (adjoint quantities) substituted. The same Python algorithm that inverts (3.10) will also invert (3.20), provided the user supplies adjoint-probe coefficients $R^{p'}_{\sigma\mu\nu}$ instead of $R^p_{\sigma\mu\nu}$.

### Reciprocal antennas

If both antennas are reciprocal, primes drop and the two formulas (3.10), (3.20) become identical via:

$$
T^p_{\sigma,-\mu,\nu} = (-1)^\mu\,R^p_{\sigma\mu\nu} \tag{3.21}
$$

$$
R_{s,-m,n} = (-1)^m\,T_{smn} \tag{3.22}
$$

---

## 5. Direction of Transmission — Practical Application (§3.2.4)

### Transmitting case (3.23) — same as (3.10):

$$
w(A, \chi, \theta, \phi) = \frac{v}{2}\sum T_{smn}\,e^{im\phi}\,d^n_{\mu m}(\theta)\,e^{i\mu\chi}\,C^{sn(3)}_{\sigma\mu\nu}(kA)\,R^p_{\sigma\mu\nu}
$$

**Inputs needed by inversion algorithm:** measured $w(A, \chi, \theta, \phi)$, probe $R^p_{\sigma\mu\nu}$, and excitation $v$.
**Output:** test-antenna $T_{smn}$.

### Receiving case (3.24) — same as (3.20):

$$
w_t(A, \chi, \theta, \phi) = \frac{v_p}{2}\sum T'_{smn}\,e^{im\phi}\,d^n_{\mu m}(\theta)\,e^{i\mu\chi}\,C^{sn(3)}_{\sigma\mu\nu}(kA)\,R^{p'}_{\sigma\mu\nu}
$$

**Inputs needed:** measured $w_t$, adjoint-probe $R^{p'}_{\sigma\mu\nu}$, and excitation $v_p$.
**Output:** adjoint test-antenna $T'_{smn}$. Apply (3.19) for the physical receiving coefficients $R_{smn}$.

### Single-algorithm philosophy

| Case | Probe coefficient input | Algorithm output | Physical answer |
|---|---|---|---|
| Test antenna transmits | $R^p_{\sigma\mu\nu}$ | $T_{smn}$ | $T_{smn}$ directly; $R_{smn}$ via reciprocity (3.22) if reciprocal |
| Test antenna receives | $R^{p'}_{\sigma\mu\nu}$ (adjoint) | $T'_{smn}$ | $R_{smn}$ via (3.19) |

**Same algorithm in both cases.**

---

## 6. Iterative Probe Calibration (§3.2.5)

A two-scan iterative scheme calibrates the probe **without** assuming knowledge of the auxiliary probe.

### Procedure

1. **Scan A:** Probe-under-calibration as "test antenna", auxiliary probe as "probe". Record probe-side data set $A$.
2. **Scan B:** Swap roles + reverse direction. Auxiliary probe as "test antenna", probe-under-calibration as "probe". Record data set $B$.

### Initial seed

$c_0$ = theoretical coefficients of a **simple source** (Hertzian dipole or Huygens source) with the correct on-axis polarization for the auxiliary probe.

### Iteration table

| Run | Measured data | Probe coefficients (input) | Output coefficients |
|---|---|---|---|
| 1 | $A$ | $c_0$ | $c_1$ |
| 2 | $B$ | $c_1$ | $c_2$ |
| 3 | $A$ | $c_2$ | $c_3$ |
| 4 | $B$ | $c_3$ | $c_4$ |
| … | … | … | … |

Process converges in **~2 steps** if $c_0$ has correct on-axis polarization. The output simultaneously calibrates **both** the probe and the auxiliary probe.

### Coefficient-type choice depends on probe nature

| Probe type | Calibrated by | Output | Conversion to $R^p$ |
|---|---|---|---|
| Non-reciprocal receiving | Receiving experiment | $T^{p'}_{\sigma\mu\nu}$ (adjoint transmitting) | $R^p_{\sigma\mu\nu}$ via (3.18) |
| Non-reciprocal transmitting | Transmitting experiment | $T^p_{\sigma\mu\nu}$ | $R^{p'}_{\sigma\mu\nu}$ via (3.18) |
| Reciprocal | Either | $T^p_{\sigma\mu\nu}$ | $R^p_{\sigma\mu\nu}$ via (3.21) |

---

## 7. Special Case: Linearly Polarized $\mu = \pm 1$ Probe (§3.3.1)

The most common probe class — conical horn fed by circular waveguide with only $\text{TE}_{11}$. Probe axis along $\hat z'$, linear polarization in the $x'z'$-plane.

### Probe coefficient structure

Only $\mu = \pm 1$ entries are non-zero, with:

$$
R^p_{1,-1,\nu} = R^p_{1,1,\nu},\quad R^p_{2,-1,\nu} = -R^p_{2,1,\nu} \tag{Eqs. 5.8, 5.9}
$$

### Probe response constants $P_{s\mu n}(kA)$

$$
\boxed{\;P_{s\mu n}(kA) = \frac{1}{2}\sum_{\sigma\nu} C^{sn(3)}_{\sigma\mu\nu}(kA)\,R^p_{\sigma\mu\nu}\;} \tag{3.26}
$$

> **Implementation note.** $P_{s\mu n}(kA)$ is precomputed once per probe and per measurement distance, then reused across all scan points $(\chi, \theta, \phi)$. Massive speedup.

Symmetry from (A3.12):

$$
P_{s,-1,n}(kA) = (-1)^{s+1}\,P_{s1n}(kA) \tag{3.27}
$$

Vanishes for $\mu \ne \pm 1$ (since the probe coefficients do).

### Reduced transmission formula

$$
w = v\sum_{smn} T_{smn}\,e^{im\phi}\bigl\{d^n_{1m}(\theta)\,e^{i\chi} + (-1)^{s+1}\,d^n_{-1,m}(\theta)\,e^{-i\chi}\bigr\}\,P_{s1n}(kA) \tag{3.28}
$$

### Far-field-pattern formulation via (A2.20)–(A2.23)

$$
\boxed{\;w(A, \chi, \theta, \phi) = v\sum_{smn} T_{smn}\,P_{s1n}(kA)\,\frac{-2i^n}{\sqrt{2n+1}}\,\vec{K}_{smn}(\theta,\phi)\cdot\hat x'\;} \tag{3.32}
$$

with $\hat x' = \hat\theta$ when $\chi = 0$ and $\hat x' = \hat\phi$ when $\chi = \pi/2$.

### Vector form (two-$\chi$ measurement)

Let $w_\theta = w(A, 0, \theta, \phi)$, $w_\phi = w(A, \pi/2, \theta, \phi)$:

$$
\boxed{\;w_\theta\,\hat\theta + w_\phi\,\hat\phi = v\sum_{smn} T_{smn}\,\Bigl\{P_{s1n}(kA)\,\frac{-2i^n}{\sqrt{2n+1}}\Bigr\}\,\vec{K}_{smn}(\theta,\phi)\;} \tag{3.33}
$$

### Direct far-field probe condition

If the probe satisfies

$$
P_{s1n}(kA) = -\tfrac{1}{2}\,i^{-n}\,\sqrt{2n+1} \tag{3.34}
$$

the probe signal **is** the test antenna far field. (These match the plane-wave expansion coefficients per Eq. A1.108 — far-field measurement is the antenna's response to a plane wave.)

### Alternative tangential-component form

Using the spherical-wave-function tangential decomposition (Eq. 3.35):

$$
[\vec{F}^{(3)}_{smn}(A, \theta, \phi)]_\text{tang} = \frac{-i^{n-s}}{\sqrt{4\pi}}\,R^{(3)}_{sn}(kA)\,\vec{K}_{smn}(\theta, \phi)
$$

and introducing the electric-dipole-probe response constant $P^e_{s1n}(kA)$ (Eq. 3.36) as a normalization:

$$
P^e_{s1n}(kA) = \tfrac{\sqrt{6}}{8}\,i^{-s}\,\sqrt{2n+1}\,R^{(3)}_{sn}(kA) \tag{3.36}
$$

we get the **tangential transmission formula**:

$$
\boxed{\;w(A, \chi, \theta, \phi) = \frac{\sqrt{6\pi}}{2}\,v\sum_{smn} \frac{P_{s1n}(kA)}{P^e_{s1n}(kA)}\,T_{smn}\,[\vec{F}^{(3)}_{smn}(A, \theta, \phi)]\cdot\hat x'\;} \tag{3.37}
$$

Vector form:

$$
w_\theta\,\hat\theta + w_\phi\,\hat\phi = \frac{\sqrt{6\pi}}{2}\,v\sum_{smn} \frac{P_{s1n}(kA)}{P^e_{s1n}(kA)}\,T_{smn}\,[\vec{F}^{(3)}_{smn}(A, \theta, \phi)]_\text{tang} \tag{3.38}
$$

> **Interpretation.** The probe acts as a "filter" applying weight $P_{s1n}/P^e_{s1n}$ to each spherical wave function in the tangential expansion of the test-antenna field at the probe location.

---

## 8. Hertzian Dipole Probes (§3.3.2)

### $\hat x'$-directed electric dipole

Non-zero receiving coefficients: $R^p_{211} = -R^p_{2,-1,1} = -\sqrt{2}/2$ (Eq. 2.154).

Response constants:

$$
P^e_{s1n}(kA) = \tfrac{\sqrt{6}}{8}\,i^{-s}\,\sqrt{2n+1}\,R^{(3)}_{sn}(kA) \tag{3.39}
$$

$$
P^e_{s,-1,n}(kA) = -\tfrac{\sqrt{6}}{8}\,i^s\,\sqrt{2n+1}\,R^{(3)}_{sn}(kA) \tag{3.40}
$$

### $\hat x'$-directed magnetic dipole (**radiates $\hat y'$-polarized**)

Non-zero receiving coefficients: $R^p_{111} = -R^p_{1,-1,1} = -i\sqrt{2}/2$ (Eq. 2.157).

$$
P^m_{s1n}(kA) = \tfrac{\sqrt{6}}{8}\,i^s\,\sqrt{2n+1}\,R^{(3)}_{3-s,n}(kA) \tag{3.41}
$$

$$
P^m_{s,-1,n}(kA) = \tfrac{\sqrt{6}}{8}\,i^{-s}\,\sqrt{2n+1}\,R^{(3)}_{3-s,n}(kA) \tag{3.42}
$$

### Received signal for electric dipole

$$
w^e(A, \chi, \theta, \phi) = \frac{\sqrt{6\pi\eta}}{2k}\Bigl\{\frac{k}{\sqrt\eta}\,v\sum_{smn} T_{smn}\,\vec{F}^{(3)}_{smn}(A, \theta, \phi)\Bigr\}\cdot\hat x' \tag{3.43}
$$

The curly bracket equals $\vec{E}(A, \theta, \phi)$ of the test antenna field, so:

$$
\boxed{\;w^e(A, \chi, \theta, \phi) = \frac{\sqrt{6\pi\eta}}{2k}\,\vec{E}(A, \theta, \phi)\cdot\hat x'\;}
$$

### Magnetic dipole

$$
w^m(A, \chi, \theta, \phi) = \frac{\sqrt{6\pi}}{2k\sqrt\eta}\,\vec{H}(A, \theta, \phi)\cdot\hat x' \tag{3.44}
$$

### On-axis special cases ($\chi = 0$ for $\hat\theta$, $\chi = \pi/2$ for $\hat\phi$)

| Equation | Component |
|---|---|
| $w^e(A, 0, \theta, \phi) = \dfrac{\sqrt{6\pi\eta}}{2k}\,E_\theta(A, \theta, \phi)$ | (3.45) |
| $w^e(A, \pi/2, \theta, \phi) = \dfrac{\sqrt{6\pi\eta}}{2k}\,E_\phi(A, \theta, \phi)$ | (3.46) |
| $w^m(A, 0, \theta, \phi) = \dfrac{\sqrt{6\pi}}{2k\sqrt\eta}\,H_\theta(A, \theta, \phi)$ | (3.47) |
| $w^m(A, \pi/2, \theta, \phi) = \dfrac{\sqrt{6\pi}}{2k\sqrt\eta}\,H_\phi(A, \theta, \phi)$ | (3.48) |

### Normalized far-field probe signal

$$
W^e(\chi, \theta, \phi) = \lim_{kA \to \infty}\Bigl[w^e(A, \chi, \theta, \phi)\,\frac{kA}{e^{ikA}}\Bigr] \tag{3.50}
$$

For electric-dipole probe:

$$
W^e(\chi, \theta, \phi) = \tfrac{\sqrt{6}}{4}\,v\,\vec{K}(\theta, \phi)\cdot\hat x' \tag{3.51}
$$

$$
W^e(0, \theta, \phi) = \tfrac{\sqrt{6}}{4}\,v\,K_\theta(\theta, \phi),\quad W^e(\pi/2, \theta, \phi) = \tfrac{\sqrt{6}}{4}\,v\,K_\phi(\theta, \phi) \tag{3.52, 3.53}
$$

### Test antenna far-field pattern from dipole-probe data

$$
\boxed{\;\vec{K}(\theta, \phi) = \frac{2\sqrt{6}}{3v}\Bigl\{W^e(0, \theta, \phi)\,\hat\theta + W^e(\pi/2, \theta, \phi)\,\hat\phi\Bigr\}\;} \tag{3.54}
$$

### Gain (matched, lossless test antenna)

$$
G(\theta, \phi) = |\vec{K}(\theta, \phi)|^2 = \frac{8}{3|v|^2}\bigl\{|W^e(0, \theta, \phi)|^2 + |W^e(\pi/2, \theta, \phi)|^2\bigr\} \tag{3.55}
$$

---

## 9. Friis' Transmission Formula (§3.3.3)

### Far-field limit of the transmission formula

For $kA \to \infty$, only $\mu = \pm 1$ terms in (3.10) survive (by A3.22–A3.24). With probe linearly $\hat x'$-polarized ($R^p_{1,-1,\nu} = R^p_{11\nu}$, $R^p_{2,-1,\nu} = -R^p_{2,1,\nu}$, $\chi = 0$) and rotation $d^n_{\mu m}(0) = \delta_{\mu m}$ when probe on $\hat z$ axis ($\theta = \phi = 0$):

$$
w \to \frac{v\,e^{ikA}}{2kA}\,\{i\vec{K}^p(\pi, \phi)\cdot\hat x\}\,\{\vec{K}(0, \phi)\cdot\hat x\}\quad\text{as } kA \to \infty \tag{3.64}
$$

### Friis' formula

Define:

$$
G_p = |\vec{K}^p(\pi, \phi)\cdot\hat x|^2,\quad G_t = |\vec{K}(0, \phi)\cdot\hat x|^2 \tag{3.65, 3.66}
$$

Taking $|\cdot|^2$ of (3.64):

$$
\boxed{\;\frac{\tfrac{1}{2}|w|^2}{\tfrac{1}{2}|v|^2} = \frac{G_p\,G_t}{4(kA)^2} = \Bigl(\frac{\lambda}{4\pi A}\Bigr)^2 G_p\,G_t\;} \tag{3.67, 3.57}
$$

Using normalized far-field signal $W = \lim_{kA\to\infty}[w\,kA/e^{ikA}]$ (Eq. 3.68):

$$
\frac{\tfrac{1}{2}|W|^2}{\tfrac{1}{2}|v|^2} = \frac{G_p\,G_t}{4} \tag{3.69}
$$

### Directivity version (using radiated power instead of accepted power)

$$
\frac{\tfrac{1}{2}|w|^2}{\tfrac{1}{2}\sum|v T_{smn}|^2} = \frac{G_p\,D_t}{4(kA)^2} \tag{3.70}
$$

$$
\frac{\tfrac{1}{2}|W|^2}{\tfrac{1}{2}\sum|v T_{smn}|^2} = \frac{G_p\,D_t}{4} \tag{3.71}
$$

> **Footnote:** The factor $i$ in (3.64) gives a 90° phase shift that also appears in Brown's [15] generalized reciprocity theorem.

---

## 10. Transmission Formula with Multiple Reflections (§3.4)

### Setup

Both antennas described by their full scattering matrices (Eqs. 3.72, 3.73):

$$
\begin{bmatrix}\Gamma_p & \mathbf{R}^p\\\mathbf{T}^p & \mathbf{S}^p\end{bmatrix}\begin{bmatrix}v_p\\\mathbf{a}^p\end{bmatrix} = \begin{bmatrix}w\\\mathbf{b}^p\end{bmatrix},\qquad \begin{bmatrix}\Gamma & \mathbf{R}\\\mathbf{T} & \mathbf{S}\end{bmatrix}\begin{bmatrix}v\\\mathbf{a}\end{bmatrix} = \begin{bmatrix}w_t\\\mathbf{b}\end{bmatrix}
$$

Truncation: $J = 2N(N+2)$ for test antenna, $J_p = 2N_p(N_p+2)$ for probe.

### Inter-antenna coupling matrices

$\mathbf{G}^+$ takes test-antenna outgoing modes to probe-incoming modes:

$$
G^+_{\beta i} = \tfrac{1}{2}\,e^{im\phi}\,d^n_{\mu m}(\theta)\,e^{i\mu\chi}\,C^{sn(3)}_{\sigma\mu\nu}(kA) \tag{3.80}
$$

(index $i$ ↔ $(s,m,n)$ for test antenna; index $\beta$ ↔ $(\sigma,\mu,\nu)$ for probe)

$\mathbf{G}^-$ takes probe-outgoing modes to test-antenna-incoming modes:

$$
G^-_{j\alpha} = \tfrac{1}{2}\,C^{\sigma\nu(3)}_{s\mu n}(-kA)\,e^{-im\phi}\,d^n_{m\mu}(-\theta)\,e^{-i\mu\chi} \tag{3.84}
$$

Couplings:

$$
\mathbf{a}^p = \mathbf{G}^+(\mathbf{b} - \mathbf{a}) \tag{3.81}
$$

$$
\mathbf{a} = \mathbf{G}^-(\mathbf{b}^p - \mathbf{a}^p) \tag{3.83}
$$

### Reduced 2-port system

After eliminating $\mathbf{a}^p, \mathbf{b}^p, \mathbf{a}, \mathbf{b}$:

$$
\begin{bmatrix}M_{11} & M_{12}\\M_{21} & M_{22}\end{bmatrix}\begin{bmatrix}v\\v_p\end{bmatrix} = \begin{bmatrix}w_t\\w\end{bmatrix} \tag{3.85}
$$

with

$$
M_{11} = \Gamma + \mathbf{R}\{\mathbf{I} - \mathbf{G}^-(\mathbf{S}^p - \mathbf{I}^p)\mathbf{G}^+(\mathbf{S} - \mathbf{I})\}^{-1}\mathbf{G}^-(\mathbf{S}^p - \mathbf{I}^p)\mathbf{G}^+\mathbf{T} \tag{3.86}
$$

$$
M_{12} = \mathbf{R}\{\mathbf{I} - \mathbf{G}^-(\mathbf{S}^p - \mathbf{I}^p)\mathbf{G}^+(\mathbf{S} - \mathbf{I})\}^{-1}\mathbf{G}^-\mathbf{T}^p \tag{3.87}
$$

$$
M_{21} = \mathbf{R}^p\{\mathbf{I}^p - \mathbf{G}^+(\mathbf{S} - \mathbf{I})\mathbf{G}^-(\mathbf{S}^p - \mathbf{I}^p)\}^{-1}\mathbf{G}^+\mathbf{T} \tag{3.88}
$$

$$
M_{22} = \Gamma_p + \mathbf{R}^p\{\mathbf{I}^p - \mathbf{G}^+(\mathbf{S} - \mathbf{I})\mathbf{G}^-(\mathbf{S}^p - \mathbf{I}^p)\}^{-1}\mathbf{G}^+(\mathbf{S} - \mathbf{I})\mathbf{G}^-\mathbf{T}^p \tag{3.89}
$$

### Full transmission formula with multi-reflection (including generator and load)

Adding generator ($v = v_g + \Gamma_g w_t$) and load ($v_p = \Gamma_l w$):

$$
\boxed{\;w = \frac{M_{21}}{(1 - M_{11}\Gamma_g)(1 - M_{22}\Gamma_l) - M_{21}M_{12}\Gamma_g\Gamma_l}\,v_g\;} \tag{3.90}
$$

### Neumann series expansion

$$
\{\mathbf{I} - \mathbf{G}^-(\mathbf{S}^p - \mathbf{I}^p)\mathbf{G}^+(\mathbf{S} - \mathbf{I})\}^{-1} = \mathbf{I} + \{\mathbf{G}^-(\mathbf{S}^p - \mathbf{I}^p)\mathbf{G}^+(\mathbf{S} - \mathbf{I})\} + \{\mathbf{G}^-(\mathbf{S}^p - \mathbf{I}^p)\mathbf{G}^+(\mathbf{S} - \mathbf{I})\}^2 + \dots \tag{3.91}
$$

- First term ($\mathbf{I}$): direct path (recovers the no-multi-reflection (3.10))
- Second term: triply-scattered field (probe → test antenna → probe)
- Higher orders: multiple-bounce paths
- Series converges for $A$ larger than some finite distance; → 0 as $kA \to \infty$

### Empty space case

If no probe present, $\mathbf{S}^p = \mathbf{I}^p$ → $M_{11} = \Gamma$ (only test-antenna reflection survives, as expected).

### Mitigation: microwave isolator

Inserting a matched isolator between probe and load makes the back-scattering ~10 dB lower than for short-circuited probe. The flow-graph node $v_p$ vanishes, simplifying the system (Fig. 3.9).

---

## 11. Implementation Checklist for Python Code

### Geometry and rotation cascade
1. **Coordinate convention.** Test antenna at origin (unprimed); probe at $(A, \theta, \phi)$ in unprimed, with its own primed frame; probe points to origin. Probe roll $\chi$ about probe $\hat z'$.
2. **Euler-angle pipeline (3.6).** Implement as four sequential operations matching Appendix A2: phase $e^{im\phi}$ → mix via $d^n_{\mu m}(\theta)$ → phase $e^{i\mu\chi}$ → mix via $C^{sn(3)}_{\sigma\mu\nu}(kA)$. **Critical:** the printed Euler tuple is $(\chi, \theta, \phi)$ but the operations apply in order $\phi \to \theta \to \chi$.

### Core transmission formula (3.10)
3. **Forward synthesis test.** Given test antenna $T_{smn}$ (e.g., from Section 2.3.4 dipole matrices), probe $R^p_{\sigma\mu\nu}$ (e.g., from Eq. 2.154), and scan parameters $(A, \chi, \theta, \phi)$, compute $w$ via (3.10). Compare to:
 - Direct evaluation of (3.43) for a Hertzian-dipole probe (electric or magnetic)
 - Direct evaluation of (3.45)–(3.48) on-axis
 - Far-field Friis (3.67) at very large $kA$
4. **Order-of-summation indifference.** $(s, m, n)$ and $(\sigma, \mu, \nu)$ can be summed in any order; verify with summation-axis permutation tests.
5. **Index range.** $s, \sigma \in \{1, 2\}$; $n, \nu \ge 1$; $m \in [-n, n]$; $\mu \in [-\nu, \nu]$. Test that values outside these ranges are explicitly zero (catches off-by-one bugs).

### Probe response constants (3.26)
6. **Precompute $P_{s\mu n}(kA)$ once.** Build a lookup table for all $(s, \mu, n)$ at the measurement distance. Verify symmetry (3.27): `P[s, -1, n] == (-1)**(s+1) * P[s, 1, n]`.
7. **Vanishing for $|\mu| \ne 1$ probes.** A linearly polarized $\mu = \pm 1$ probe must satisfy $P_{s\mu n}(kA) = 0$ for $\mu \notin \{-1, +1\}$.
8. **Reduced formula equivalence.** For a $\mu = \pm 1$ probe, verify that (3.10) gives the same result as the reduced form (3.28) — and that (3.28) gives the same as (3.32) via the $\vec{K}_{smn}$ identity.

### Hertzian dipole probe regression
9. **Electric dipole response (3.39).** Closed-form: $P^e_{s1n}(kA) = (\sqrt{6}/8)\,i^{-s}\sqrt{2n+1}\,R^{(3)}_{sn}(kA)$. Test for $n = 1\dots 5$ at $kA = 5, 10, 20$.
10. **Magnetic dipole response (3.41).** Note the $R^{(3)}_{3-s,n}$ (index-swapped) and the $i^s$ (not $i^{-s}$).
11. **On-axis special cases.** Equations (3.45)–(3.48): probe-received signal equals $E_\theta, E_\phi, H_\theta, H_\phi$ scaled by $\sqrt{6\pi\eta}/(2k)$ or $\sqrt{6\pi}/(2k\sqrt\eta)$. Direct integration check.
12. **Far-field gain reconstruction (3.55).** $G = (8/(3|v|^2))(|W^e_0|^2 + |W^e_{\pi/2}|^2)$.

### Direction-of-transmission equivalence
13. **Same algorithm both ways.** Verify that the inversion algorithm for (3.10) (test transmits) and (3.20) (test receives, via adjoint) returns the **same** physical pattern (when reciprocal). Critical regression: build a reciprocal test antenna, simulate both forward and reverse measurements, recover $T_{smn}$ in each — both should match the input.
14. **Adjoint substitution (3.18, 3.19).** `T^p[σ, -μ, ν] = (-1)**μ * R^{p'}[σ, μ, ν]` and `R[s, -m, n] = (-1)**m * T'[s, m, n]`.
15. **Reciprocal collapse.** When both antennas are reciprocal, (3.10) and (3.20) must yield bit-identical results given the same scan geometry and excitation.

### Multiple reflections (§3.4)
16. **No-reflection limit.** Setting $\mathbf{S}^p = \mathbf{I}^p$ and $\mathbf{S} = \mathbf{I}$ in (3.86)–(3.89) collapses $M_{21}$ to the single-reflection $\mathbf{R}^p\mathbf{G}^+\mathbf{T}$, recovering (3.10) (after the appropriate generator/load handling).
17. **Neumann series convergence.** For a chosen $A$, verify that successive terms in (3.91) decay; report the largest term as a "multi-reflection figure of merit."
18. **Empty-space limit.** $M_{11} = \Gamma$ when $\mathbf{S}^p \to \mathbf{I}^p$.

### Numerical sanity
19. **Far-field limit (Friis).** Take $kA = 1000$ with a Hertzian-dipole probe. The full transmission formula (3.10) should match Friis (3.67) to ~3 significant figures (limited by the $1/(kA)$ corrections).
20. **Probe-scaling invariance.** Multiplying $R^p_{\sigma\mu\nu}$ by a scalar $\alpha$ should multiply the predicted $w$ by $\alpha$. Useful for catching unit-confusion bugs.
21. **Tangential consistency.** For a $\mu = \pm 1$ probe with two measurements at $\chi = 0, \pi/2$, the inferred field should satisfy $\vec E \cdot \hat r = 0$ on the measurement sphere (no radial component in the tangential frame).
22. **Sign of $A$.** $A > 0$ in (3.10). Use Eq. (A3.13) if you need translation coefficients with $kA < 0$ — never recompute (3.6) for $A < 0$.

---

*End of Chapter 3 reference.*
