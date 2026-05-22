# Chapter 3 — Scattering Matrix Description of Antenna Coupling

**Source:** Hansen, *Spherical Near-Field Antenna Measurements*, Chapter 3 (pp. 61–88).

**Cross-references:**
- Symbols: [A5_List_of_principal_symbols_and_uses.md](A5_List_of_principal_symbols_and_uses.md)
- Chapter 2 (Scattering matrix theory): [02_Scattering_matrix_description_of_an_antenna.md](02_Scattering_matrix_description_of_an_antenna.md)
- Appendix A1 (Wave functions, $\\vec{K}\_{smn}$): [A1_Spherical_wave_functions_notation_and_properties.md](A1_Spherical_wave_functions_notation_and_properties.md)
- Appendix A2 (Rotation, $d^n\_{\\mu m}$, Euler angles): [A2_Rotation_of_spherical_waves.md](A2_Rotation_of_spherical_waves.md)
- Appendix A3 (Translation, $C^{sn(c)}\_{\\sigma\\mu\\nu}$): [A3_Translation_of_spherical_waves.md](A3_Translation_of_spherical_waves.md)

> **Purpose.** Computational reference for the **spherical transmission formula** — the central equation of probe-corrected spherical near-field measurements. Every formula is tagged by book number `(3.N)`. The transmission formula (3.10) and its companions (3.17), (3.20) are the boxed equations that any Python implementation must reproduce exactly.

---

## 1. Geometry (§3.2.1)

| Frame | Origin | Coordinates | Holds |
|---|---|---|---|
| Unprimed $(x, y, z)$ | Test antenna | $(r, \\theta, \\phi)$ | Test antenna fixed; minimum sphere radius $r\_0$ |
| Primed $(x', y', z')$ | Probe | $(r', \\theta', \\phi')$ | Probe pointing **at the unprimed origin**; minimum sphere radius $r'\_0$ |
| Probe position in unprimed | — | $(A, \\theta, \\phi)$ | Probe origin at distance $A$ from test-antenna origin, in direction $(\\theta, \\phi)$ |
| Probe roll | — | $\\chi$ | Rotation of probe about its own $z'$ axis |

**Probe minimum sphere does not intersect test-antenna minimum sphere:** $A > r\_0 + r'\_0$.

Typical scan values: $\\chi \\in \\{0, \\pi/2\\}$ (two roll positions sufficient for full polarization recovery).

---

## 2. Test Antenna Transmitting, Probe Receiving (§3.2.2)

### Starting field (test antenna only)

|  Eq.  |   |
| :---: | :-- |
| (3.1) | $\\vec{E}\_t(r,\\theta,\\phi) = \\frac{k}{\\sqrt{\\eta}}\\sum\_{smn} Q^{(3)}\_{smn}\\,\\vec{F}^{(3)}\_{smn}(r,\\theta,\\phi) = \\frac{k}{\\sqrt{\\eta}}\\sum\_{smn} v\\,T\_{smn}\\,\\vec{F}^{(3)}\_{smn}(r,\\theta,\\phi),\\quad r > r\_0$ |

(Using $Q^{(3)}\_{smn} = b\_{smn} = v\\,T\_{smn}$ from Eq. 2.66.)

### The four-step coordinate cascade (test-antenna frame → probe frame)

Each spherical mode $\\vec{F}^{(3)}\_{smn}$ is re-expressed in the probe frame by:

| Step | Transform | Mathematical action |
|---|---|---|
| 1 | Rotate frame 1 by $\\phi\_o$ about $\\hat z$ | Phase $e^{im\\phi\_o}$ — Eq. (3.2) |
| 2 | Rotate frame 2 by $\\theta\_o$ about $\\hat y\_1$ | Mix over $\\mu$ via rotation coefficient $d^n\_{\\mu m}(\\theta\_o)$ — Eq. (3.3) |
| 3 | Rotate frame 3 by $\\chi\_o$ about $\\hat z\_2$ | Phase $e^{i\\mu\\chi\_o}$ — Eq. (3.4) |
| 4 | Translate by $A$ along $\\hat z\_3$ | Mix over $(\\sigma, \\nu)$ via $C^{sn(3)}\_{\\sigma\\mu\\nu}(kA)$, with $\\vec{F}^{(3)} \\to \\tfrac{1}{2}(\\vec{F}^{(3)} + \\vec{F}^{(4)})$ — Eq. (3.5) |

> **Why the $\\tfrac{1}{2}(\\vec{F}^{(3)} + \\vec{F}^{(4)})$ split?** After axial translation, the original outgoing wave $\\vec{F}^{(3)}$ from the test-antenna origin appears in the probe frame as a standing wave $\\vec{F}^{(1)} = \\tfrac{1}{2}(\\vec{F}^{(3)} + \\vec{F}^{(4)})$ inside the sphere of radius $A - r\_0$ around the probe origin (where the probe is located). The probe sees the $\\vec{F}^{(4)}$ component as the **incoming wave**.

### Composite transform formula

|  Eq.  |   |
| :---: | :-- |
| (3.6) | $\\boxed{\\;\\vec{F}^{(3)}\_{smn}(r,\\theta,\\phi) = \\sum\_{\\sigma\\mu\\nu} e^{im\\phi\_o}\\,d^n\_{\\mu m}(\\theta\_o)\\,e^{i\\mu\\chi\_o}\\,C^{sn(3)}\_{\\sigma\\mu\\nu}(kA)\\,\\tfrac{1}{2}\\bigl\\{\\vec{F}^{(3)}\_{\\sigma\\mu\\nu}(r',\\theta',\\phi') + \\vec{F}^{(4)}\_{\\sigma\\mu\\nu}(r',\\theta',\\phi')\\bigr\\}\\;}$ |

### Field expansion in the probe frame

|  Eq.  |   |
| :---: | :-- |
| (3.7) | $\\vec{E}\_t = \\frac{k}{\\sqrt{\\eta}}\\sum\_{\\substack{smn\\\\\\sigma\\mu\\nu}} v\\,T\_{smn}\\,e^{im\\phi\_o}\\,d^n\_{\\mu m}(\\theta\_o)\\,e^{i\\mu\\chi\_o}\\,C^{sn(3)}\_{\\sigma\\mu\\nu}(kA)\\,\\tfrac{1}{2}\\{\\vec{F}^{(3)}\_{\\sigma\\mu\\nu} + \\vec{F}^{(4)}\_{\\sigma\\mu\\nu}\\}$ |

### Incoming-wave coefficients at the probe

The probe sees $\\vec{F}^{(4)}\_{\\sigma\\mu\\nu}$ as incoming with amplitude:

|  Eq.  |   |
| :---: | :-- |
| (3.9) | $a\_{\\sigma\\mu\\nu} = \\frac{v}{2}\\sum\_{smn} T\_{smn}\\,e^{im\\phi\_o}\\,d^n\_{\\mu m}(\\theta\_o)\\,e^{i\\mu\\chi\_o}\\,C^{sn(3)}\_{\\sigma\\mu\\nu}(kA)$ |

### Probe-received signal (matched load: $\\Gamma\_l = 0$)

Using $w = \\mathbf{R}^p\\,\\mathbf{a}$ (Eq. 2.71 with $\\Gamma\_l = 0$) and dropping the $o$-subscript on probe angles:

|  Eq.  |   |
| :---: | :-- |
| (3.10) | $\\boxed{\\;w(A, \\chi, \\theta, \\phi) = \\frac{v}{2}\\sum\_{\\substack{smn\\\\\\sigma\\mu\\nu}} T\_{smn}\\,e^{im\\phi}\\,d^n\_{\\mu m}(\\theta)\\,e^{i\\mu\\chi}\\,C^{sn(3)}\_{\\sigma\\mu\\nu}(kA)\\,R^p\_{\\sigma\\mu\\nu}\\;}$ |

> **The Transmission Formula.** This is the central equation of the entire book. It expresses the complex probe-received signal as a function of:
> - **Inputs (knowns):** probe coordinates $(A, \\chi, \\theta, \\phi)$, probe receiving coefficients $R^p\_{\\sigma\\mu\\nu}$, source amplitude $v$.
> - **Unknowns:** test-antenna transmitting coefficients $T\_{smn}$.
>
> The role of Chapter 4 is to **invert (3.10)** to recover $T\_{smn}$ from sampled $w$-data. **No reciprocity assumed for either antenna.**

### Assumptions used in deriving (3.10)

1. **No multiple reflections** between probe and test antenna (probe-scattered field does not perturb the test antenna). Relaxed in §9 below.
2. **Matched load** at the probe ($\\Gamma\_l = 0$). Relaxed by the gain-measurement discussion in Chapter 5.

---

## 3. Test Antenna Receiving, Probe Transmitting (§3.2.3)

### Probe-radiated field in the probe frame

|  Eq.  |   |
| :---: | :-- |
| (3.11) | $\\vec{E}\_p(r',\\theta',\\phi') = \\frac{k}{\\sqrt{\\eta}}\\sum\_{\\sigma\\mu\\nu} v\_p\\,T^p\_{\\sigma\\mu\\nu}\\,\\vec{F}^{(3)}\_{\\sigma\\mu\\nu}(r',\\theta',\\phi'),\\quad r' > r'\_0$ |

### Inverse coordinate cascade

Apply the same four-step cascade in **reverse order with negated angles** $(-A, -\\chi\_o, -\\theta\_o, -\\phi\_o)$. Result expressing $\\vec{E}\_p$ in the test-antenna frame:

|  Eq.  |   |
| :---: | :-- |
| (3.12) | $\\vec{E}\_p = \\frac{k}{\\sqrt{\\eta}}\\sum\_{\\substack{\\sigma\\mu\\nu\\\\smn}} v\_p\\,T^p\_{\\sigma\\mu\\nu}\\,C^{\\sigma\\nu(3)}\_{s\\mu n}(-kA)\\,e^{-i\\mu\\chi\_o}\\,d^n\_{m\\mu}(-\\theta\_o)\\,e^{-im\\phi\_o}\\,\\tfrac{1}{2}\\{\\vec{F}^{(3)}\_{smn} + \\vec{F}^{(4)}\_{smn}\\}$ |

### Test antenna received signal (matched load: $\\Gamma = 0$)

Initial form:

|  Eq.  |   |
| :---: | :-- |
| (3.14) | $a\_{smn} = \\frac{v\_p}{2}\\sum\_{\\sigma\\mu\\nu} T^p\_{\\sigma\\mu\\nu}\\,C^{\\sigma\\nu(3)}\_{s\\mu n}(-kA)\\,e^{-i\\mu\\chi\_o}\\,d^n\_{m\\mu}(-\\theta\_o)\\,e^{-im\\phi\_o}$ |
| (3.15) | $w\_t = \\frac{v\_p}{2}\\sum\_{\\substack{\\sigma\\mu\\nu\\\\smn}} T^p\_{\\sigma\\mu\\nu}\\,C^{\\sigma\\nu(3)}\_{s\\mu n}(-kA)\\,e^{-i\\mu\\chi}\\,d^n\_{m\\mu}(-\\theta)\\,e^{-im\\phi}\\,R\_{smn}$ |

### Simplification via symmetry relations (A3.13) and (A2.7)

|  Eq.  |   |
| :---: | :-- |
| (3.16) | $w\_t = \\frac{v\_p}{2}\\sum\_{\\substack{\\sigma\\mu\\nu\\\\smn}} T^p\_{\\sigma\\mu\\nu}\\,C^{sn(3)}\_{\\sigma,-\\mu,\\nu}(kA)\\,e^{-i\\mu\\chi}\\,d^n\_{\\mu m}(\\theta)\\,e^{-im\\phi}\\,R\_{smn}$ |

### Final form via (A2.9) and relabeling

|  Eq.  |   |
| :---: | :-- |
| (3.17) | $\\boxed{\\;w\_t = \\frac{v\_p}{2}\\sum\_{\\substack{\\sigma\\mu\\nu\\\\smn}}(-1)^\\mu\\,T^p\_{\\sigma,-\\mu,\\nu}\\,C^{sn(3)}\_{\\sigma\\mu\\nu}(kA)\\,e^{i\\mu\\chi}\\,d^n\_{\\mu m}(\\theta)\\,e^{im\\phi}\\,(-1)^m\\,R\_{s,-m,n}\\;}$ |

---

## 4. Adjoint Substitution — Recasting (3.17) (§3.2.3)

Using the adjoint relations (Eqs. 2.103, 2.104) between original and adjoint antenna coefficients:

|  Eq.  |   |
| :---: | :-- |
| (3.18) | $T^p\_{\\sigma,-\\mu,\\nu} = (-1)^\\mu\\,R^{p'}\_{\\sigma\\mu\\nu}$ |
| (3.19) | $R\_{s,-m,n} = (-1)^m\\,T'\_{smn}$ |

The receiving formula (3.17) becomes:

|  Eq.  |   |
| :---: | :-- |
| (3.20) | $\\boxed{\\;w\_t = \\frac{v\_p}{2}\\sum\_{\\substack{\\sigma\\mu\\nu\\\\smn}} R^{p'}\_{\\sigma\\mu\\nu}\\,C^{sn(3)}\_{\\sigma\\mu\\nu}(kA)\\,e^{i\\mu\\chi}\\,d^n\_{\\mu m}(\\theta)\\,e^{im\\phi}\\,T'\_{smn}\\;}$ |

> **Key observation.** Equation (3.20) has the **same structural form** as the transmitting-case formula (3.10) — only with primes (adjoint quantities) substituted. The same Python algorithm that inverts (3.10) will also invert (3.20), provided the user supplies adjoint-probe coefficients $R^{p'}\_{\\sigma\\mu\\nu}$ instead of $R^p\_{\\sigma\\mu\\nu}$.

### Reciprocal antennas

If both antennas are reciprocal, primes drop and the two formulas (3.10), (3.20) become identical via:

|  Eq.  |   |
| :---: | :-- |
| (3.21) | $T^p\_{\\sigma,-\\mu,\\nu} = (-1)^\\mu\\,R^p\_{\\sigma\\mu\\nu}$ |
| (3.22) | $R\_{s,-m,n} = (-1)^m\\,T\_{smn}$ |

---

## 5. Direction of Transmission — Practical Application (§3.2.4)

### Transmitting case (3.23) — same as (3.10):

|  Eq.  |   |
| :---: | :-- |
|  | $w(A, \\chi, \\theta, \\phi) = \\frac{v}{2}\\sum T\_{smn}\\,e^{im\\phi}\\,d^n\_{\\mu m}(\\theta)\\,e^{i\\mu\\chi}\\,C^{sn(3)}\_{\\sigma\\mu\\nu}(kA)\\,R^p\_{\\sigma\\mu\\nu}$ |

**Inputs needed by inversion algorithm:** measured $w(A, \\chi, \\theta, \\phi)$, probe $R^p\_{\\sigma\\mu\\nu}$, and excitation $v$.
**Output:** test-antenna $T\_{smn}$.

### Receiving case (3.24) — same as (3.20):

|  Eq.  |   |
| :---: | :-- |
|  | $w\_t(A, \\chi, \\theta, \\phi) = \\frac{v\_p}{2}\\sum T'\_{smn}\\,e^{im\\phi}\\,d^n\_{\\mu m}(\\theta)\\,e^{i\\mu\\chi}\\,C^{sn(3)}\_{\\sigma\\mu\\nu}(kA)\\,R^{p'}\_{\\sigma\\mu\\nu}$ |

**Inputs needed:** measured $w\_t$, adjoint-probe $R^{p'}\_{\\sigma\\mu\\nu}$, and excitation $v\_p$.
**Output:** adjoint test-antenna $T'\_{smn}$. Apply (3.19) for the physical receiving coefficients $R\_{smn}$.

### Single-algorithm philosophy

| Case | Probe coefficient input | Algorithm output | Physical answer |
|---|---|---|---|
| Test antenna transmits | $R^p\_{\\sigma\\mu\\nu}$ | $T\_{smn}$ | $T\_{smn}$ directly; $R\_{smn}$ via reciprocity (3.22) if reciprocal |
| Test antenna receives | $R^{p'}\_{\\sigma\\mu\\nu}$ (adjoint) | $T'\_{smn}$ | $R\_{smn}$ via (3.19) |

**Same algorithm in both cases.**

---

## 6. Iterative Probe Calibration (§3.2.5)

A two-scan iterative scheme calibrates the probe **without** assuming knowledge of the auxiliary probe.

### Procedure

1. **Scan A:** Probe-under-calibration as "test antenna", auxiliary probe as "probe". Record probe-side data set $A$.
2. **Scan B:** Swap roles + reverse direction. Auxiliary probe as "test antenna", probe-under-calibration as "probe". Record data set $B$.

### Initial seed

$c\_0$ = theoretical coefficients of a **simple source** (Hertzian dipole or Huygens source) with the correct on-axis polarization for the auxiliary probe.

### Iteration table

| Run | Measured data | Probe coefficients (input) | Output coefficients |
|---|---|---|---|
| 1 | $A$ | $c\_0$ | $c\_1$ |
| 2 | $B$ | $c\_1$ | $c\_2$ |
| 3 | $A$ | $c\_2$ | $c\_3$ |
| 4 | $B$ | $c\_3$ | $c\_4$ |
| … | … | … | … |

Process converges in **~2 steps** if $c\_0$ has correct on-axis polarization. The output simultaneously calibrates **both** the probe and the auxiliary probe.

### Coefficient-type choice depends on probe nature

| Probe type | Calibrated by | Output | Conversion to $R^p$ |
|---|---|---|---|
| Non-reciprocal receiving | Receiving experiment | $T^{p'}\_{\\sigma\\mu\\nu}$ (adjoint transmitting) | $R^p\_{\\sigma\\mu\\nu}$ via (3.18) |
| Non-reciprocal transmitting | Transmitting experiment | $T^p\_{\\sigma\\mu\\nu}$ | $R^{p'}\_{\\sigma\\mu\\nu}$ via (3.18) |
| Reciprocal | Either | $T^p\_{\\sigma\\mu\\nu}$ | $R^p\_{\\sigma\\mu\\nu}$ via (3.21) |

---

## 7. Special Case: Linearly Polarized $\\mu = \\pm 1$ Probe (§3.3.1)

The most common probe class — conical horn fed by circular waveguide with only $\\text{TE}\_{11}$. Probe axis along $\\hat z'$, linear polarization in the $x'z'$-plane.

### Probe coefficient structure

Only $\\mu = \\pm 1$ entries are non-zero, with:

|  Eq.  |   |
| :---: | :-- |
| (Eqs. 5.8, 5.9) | $R^p\_{1,-1,\\nu} = R^p\_{1,1,\\nu},\\quad R^p\_{2,-1,\\nu} = -R^p\_{2,1,\\nu}$ |

### Probe response constants $P\_{s\\mu n}(kA)$

|  Eq.  |   |
| :---: | :-- |
| (3.26) | $\\boxed{\\;P\_{s\\mu n}(kA) = \\frac{1}{2}\\sum\_{\\sigma\\nu} C^{sn(3)}\_{\\sigma\\mu\\nu}(kA)\\,R^p\_{\\sigma\\mu\\nu}\\;}$ |

> **Implementation note.** $P\_{s\\mu n}(kA)$ is precomputed once per probe and per measurement distance, then reused across all scan points $(\\chi, \\theta, \\phi)$. Massive speedup.

Symmetry from (A3.12):

|  Eq.  |   |
| :---: | :-- |
| (3.27) | $P\_{s,-1,n}(kA) = (-1)^{s+1}\\,P\_{s1n}(kA)$ |

Vanishes for $\\mu \\ne \\pm 1$ (since the probe coefficients do).

### Reduced transmission formula

|  Eq.  |   |
| :---: | :-- |
| (3.28) | $w = v\\sum\_{smn} T\_{smn}\\,e^{im\\phi}\\bigl\\{d^n\_{1m}(\\theta)\\,e^{i\\chi} + (-1)^{s+1}\\,d^n\_{-1,m}(\\theta)\\,e^{-i\\chi}\\bigr\\}\\,P\_{s1n}(kA)$ |

### Far-field-pattern formulation via (A2.20)–(A2.23)

|  Eq.  |   |
| :---: | :-- |
| (3.32) | $\\boxed{\\;w(A, \\chi, \\theta, \\phi) = v\\sum\_{smn} T\_{smn}\\,P\_{s1n}(kA)\\,\\frac{-2i^n}{\\sqrt{2n+1}}\\,\\vec{K}\_{smn}(\\theta,\\phi)\\cdot\\hat x'\\;}$ |

with $\\hat x' = \\hat\\theta$ when $\\chi = 0$ and $\\hat x' = \\hat\\phi$ when $\\chi = \\pi/2$.

### Vector form (two-$\\chi$ measurement)

Let $w\_\\theta = w(A, 0, \\theta, \\phi)$, $w\_\\phi = w(A, \\pi/2, \\theta, \\phi)$:

|  Eq.  |   |
| :---: | :-- |
| (3.33) | $\\boxed{\\;w\_\\theta\\,\\hat\\theta + w\_\\phi\\,\\hat\\phi = v\\sum\_{smn} T\_{smn}\\,\\Bigl\\{P\_{s1n}(kA)\\,\\frac{-2i^n}{\\sqrt{2n+1}}\\Bigr\\}\\,\\vec{K}\_{smn}(\\theta,\\phi)\\;}$ |

### Direct far-field probe condition

If the probe satisfies

|  Eq.  |   |
| :---: | :-- |
| (3.34) | $P\_{s1n}(kA) = -\\tfrac{1}{2}\\,i^{-n}\\,\\sqrt{2n+1}$ |

the probe signal **is** the test antenna far field. (These match the plane-wave expansion coefficients per Eq. A1.108 — far-field measurement is the antenna's response to a plane wave.)

### Alternative tangential-component form

Using the spherical-wave-function tangential decomposition (Eq. 3.35):

|  Eq.  |   |
| :---: | :-- |
|  | $[\\vec{F}^{(3)}\_{smn}(A, \\theta, \\phi)]\_\\text{tang} = \\frac{-i^{n-s}}{\\sqrt{4\\pi}}\\,R^{(3)}\_{sn}(kA)\\,\\vec{K}\_{smn}(\\theta, \\phi)$ |

and introducing the electric-dipole-probe response constant $P^e\_{s1n}(kA)$ (Eq. 3.36) as a normalization:

|  Eq.  |   |
| :---: | :-- |
| (3.36) | $P^e\_{s1n}(kA) = \\tfrac{\\sqrt{6}}{8}\\,i^{-s}\\,\\sqrt{2n+1}\\,R^{(3)}\_{sn}(kA)$ |

we get the **tangential transmission formula**:

|  Eq.  |   |
| :---: | :-- |
| (3.37) | $\\boxed{\\;w(A, \\chi, \\theta, \\phi) = \\frac{\\sqrt{6\\pi}}{2}\\,v\\sum\_{smn} \\frac{P\_{s1n}(kA)}{P^e\_{s1n}(kA)}\\,T\_{smn}\\,[\\vec{F}^{(3)}\_{smn}(A, \\theta, \\phi)]\\cdot\\hat x'\\;}$ |

Vector form:

|  Eq.  |   |
| :---: | :-- |
| (3.38) | $w\_\\theta\\,\\hat\\theta + w\_\\phi\\,\\hat\\phi = \\frac{\\sqrt{6\\pi}}{2}\\,v\\sum\_{smn} \\frac{P\_{s1n}(kA)}{P^e\_{s1n}(kA)}\\,T\_{smn}\\,[\\vec{F}^{(3)}\_{smn}(A, \\theta, \\phi)]\_\\text{tang}$ |

> **Interpretation.** The probe acts as a "filter" applying weight $P\_{s1n}/P^e\_{s1n}$ to each spherical wave function in the tangential expansion of the test-antenna field at the probe location.

---

## 8. Hertzian Dipole Probes (§3.3.2)

### $\\hat x'$-directed electric dipole

Non-zero receiving coefficients: $R^p\_{211} = -R^p\_{2,-1,1} = -\\sqrt{2}/2$ (Eq. 2.154).

Response constants:

|  Eq.  |   |
| :---: | :-- |
| (3.39) | $P^e\_{s1n}(kA) = \\tfrac{\\sqrt{6}}{8}\\,i^{-s}\\,\\sqrt{2n+1}\\,R^{(3)}\_{sn}(kA)$ |
| (3.40) | $P^e\_{s,-1,n}(kA) = -\\tfrac{\\sqrt{6}}{8}\\,i^s\\,\\sqrt{2n+1}\\,R^{(3)}\_{sn}(kA)$ |

### $\\hat x'$-directed magnetic dipole (**radiates $\\hat y'$-polarized**)

Non-zero receiving coefficients: $R^p\_{111} = -R^p\_{1,-1,1} = -i\\sqrt{2}/2$ (Eq. 2.157).

|  Eq.  |   |
| :---: | :-- |
| (3.41) | $P^m\_{s1n}(kA) = \\tfrac{\\sqrt{6}}{8}\\,i^s\\,\\sqrt{2n+1}\\,R^{(3)}\_{3-s,n}(kA)$ |
| (3.42) | $P^m\_{s,-1,n}(kA) = \\tfrac{\\sqrt{6}}{8}\\,i^{-s}\\,\\sqrt{2n+1}\\,R^{(3)}\_{3-s,n}(kA)$ |

### Received signal for electric dipole

|  Eq.  |   |
| :---: | :-- |
| (3.43) | $w^e(A, \\chi, \\theta, \\phi) = \\frac{\\sqrt{6\\pi\\eta}}{2k}\\Bigl\\{\\frac{k}{\\sqrt\\eta}\\,v\\sum\_{smn} T\_{smn}\\,\\vec{F}^{(3)}\_{smn}(A, \\theta, \\phi)\\Bigr\\}\\cdot\\hat x'$ |

The curly bracket equals $\\vec{E}(A, \\theta, \\phi)$ of the test antenna field, so:

|  Eq.  |   |
| :---: | :-- |
|  | $\\boxed{\\;w^e(A, \\chi, \\theta, \\phi) = \\frac{\\sqrt{6\\pi\\eta}}{2k}\\,\\vec{E}(A, \\theta, \\phi)\\cdot\\hat x'\\;}$ |

### Magnetic dipole

|  Eq.  |   |
| :---: | :-- |
| (3.44) | $w^m(A, \\chi, \\theta, \\phi) = \\frac{\\sqrt{6\\pi}}{2k\\sqrt\\eta}\\,\\vec{H}(A, \\theta, \\phi)\\cdot\\hat x'$ |

### On-axis special cases ($\\chi = 0$ for $\\hat\\theta$, $\\chi = \\pi/2$ for $\\hat\\phi$)

| Equation | Component |
|---|---|
| $w^e(A, 0, \\theta, \\phi) = \\dfrac{\\sqrt{6\\pi\\eta}}{2k}\\,E\_\\theta(A, \\theta, \\phi)$ | (3.45) |
| $w^e(A, \\pi/2, \\theta, \\phi) = \\dfrac{\\sqrt{6\\pi\\eta}}{2k}\\,E\_\\phi(A, \\theta, \\phi)$ | (3.46) |
| $w^m(A, 0, \\theta, \\phi) = \\dfrac{\\sqrt{6\\pi}}{2k\\sqrt\\eta}\\,H\_\\theta(A, \\theta, \\phi)$ | (3.47) |
| $w^m(A, \\pi/2, \\theta, \\phi) = \\dfrac{\\sqrt{6\\pi}}{2k\\sqrt\\eta}\\,H\_\\phi(A, \\theta, \\phi)$ | (3.48) |

### Normalized far-field probe signal

|  Eq.  |   |
| :---: | :-- |
| (3.50) | $W^e(\\chi, \\theta, \\phi) = \\lim\_{kA \\to \\infty}\\Bigl[w^e(A, \\chi, \\theta, \\phi)\\,\\frac{kA}{e^{ikA}}\\Bigr]$ |

For electric-dipole probe:

|  Eq.  |   |
| :---: | :-- |
| (3.51) | $W^e(\\chi, \\theta, \\phi) = \\tfrac{\\sqrt{6}}{4}\\,v\\,\\vec{K}(\\theta, \\phi)\\cdot\\hat x'$ |
| (3.52, 3.53) | $W^e(0, \\theta, \\phi) = \\tfrac{\\sqrt{6}}{4}\\,v\\,K\_\\theta(\\theta, \\phi),\\quad W^e(\\pi/2, \\theta, \\phi) = \\tfrac{\\sqrt{6}}{4}\\,v\\,K\_\\phi(\\theta, \\phi)$ |

### Test antenna far-field pattern from dipole-probe data

|  Eq.  |   |
| :---: | :-- |
| (3.54) | $\\boxed{\\;\\vec{K}(\\theta, \\phi) = \\frac{2\\sqrt{6}}{3v}\\Bigl\\{W^e(0, \\theta, \\phi)\\,\\hat\\theta + W^e(\\pi/2, \\theta, \\phi)\\,\\hat\\phi\\Bigr\\}\\;}$ |

### Gain (matched, lossless test antenna)

|  Eq.  |   |
| :---: | :-- |
| (3.55) | $G(\\theta, \\phi) = \|\\vec{K}(\\theta, \\phi)\|^2 = \\frac{8}{3\|v\|^2}\\bigl\\{\|W^e(0, \\theta, \\phi)\|^2 + \|W^e(\\pi/2, \\theta, \\phi)\|^2\\bigr\\}$ |

---

## 9. Friis' Transmission Formula (§3.3.3)

### Far-field limit of the transmission formula

For $kA \\to \\infty$, only $\\mu = \\pm 1$ terms in (3.10) survive (by A3.22–A3.24). With probe linearly $\\hat x'$-polarized ($R^p\_{1,-1,\\nu} = R^p\_{11\\nu}$, $R^p\_{2,-1,\\nu} = -R^p\_{2,1,\\nu}$, $\\chi = 0$) and rotation $d^n\_{\\mu m}(0) = \\delta\_{\\mu m}$ when probe on $\\hat z$ axis ($\\theta = \\phi = 0$):

|  Eq.  |   |
| :---: | :-- |
| (3.64) | $w \\to \\frac{v\\,e^{ikA}}{2kA}\\,\\{i\\vec{K}^p(\\pi, \\phi)\\cdot\\hat x\\}\\,\\{\\vec{K}(0, \\phi)\\cdot\\hat x\\}\\quad\\text{as } kA \\to \\infty$ |

### Friis' formula

Define:

|  Eq.  |   |
| :---: | :-- |
| (3.65, 3.66) | $G\_p = \|\\vec{K}^p(\\pi, \\phi)\\cdot\\hat x\|^2,\\quad G\_t = \|\\vec{K}(0, \\phi)\\cdot\\hat x\|^2$ |

Taking $|\\cdot|^2$ of (3.64):

|  Eq.  |   |
| :---: | :-- |
| (3.67, 3.57) | $\\boxed{\\;\\frac{\\tfrac{1}{2}\|w\|^2}{\\tfrac{1}{2}\|v\|^2} = \\frac{G\_p\\,G\_t}{4(kA)^2} = \\Bigl(\\frac{\\lambda}{4\\pi A}\\Bigr)^2 G\_p\\,G\_t\\;}$ |

Using normalized far-field signal $W = \\lim\_{kA\\to\\infty}[w\\,kA/e^{ikA}]$ (Eq. 3.68):

|  Eq.  |   |
| :---: | :-- |
| (3.69) | $\\frac{\\tfrac{1}{2}\|W\|^2}{\\tfrac{1}{2}\|v\|^2} = \\frac{G\_p\\,G\_t}{4}$ |

### Directivity version (using radiated power instead of accepted power)

|  Eq.  |   |
| :---: | :-- |
| (3.70) | $\\frac{\\tfrac{1}{2}\|w\|^2}{\\tfrac{1}{2}\\sum\|v T\_{smn}\|^2} = \\frac{G\_p\\,D\_t}{4(kA)^2}$ |
| (3.71) | $\\frac{\\tfrac{1}{2}\|W\|^2}{\\tfrac{1}{2}\\sum\|v T\_{smn}\|^2} = \\frac{G\_p\\,D\_t}{4}$ |

> **Footnote:** The factor $i$ in (3.64) gives a 90° phase shift that also appears in Brown's [15] generalized reciprocity theorem.

---

## 10. Transmission Formula with Multiple Reflections (§3.4)

### Setup

Both antennas described by their full scattering matrices (Eqs. 3.72, 3.73):

|  Eq.  |   |
| :---: | :-- |
|  | $\\begin{bmatrix}\\Gamma\_p & \\mathbf{R}^p\\\\\\mathbf{T}^p & \\mathbf{S}^p\\end{bmatrix}\\begin{bmatrix}v\_p\\\\\\mathbf{a}^p\\end{bmatrix} = \\begin{bmatrix}w\\\\\\mathbf{b}^p\\end{bmatrix},\\qquad \\begin{bmatrix}\\Gamma & \\mathbf{R}\\\\\\mathbf{T} & \\mathbf{S}\\end{bmatrix}\\begin{bmatrix}v\\\\\\mathbf{a}\\end{bmatrix} = \\begin{bmatrix}w\_t\\\\\\mathbf{b}\\end{bmatrix}$ |

Truncation: $J = 2N(N+2)$ for test antenna, $J\_p = 2N\_p(N\_p+2)$ for probe.

### Inter-antenna coupling matrices

$\\mathbf{G}^+$ takes test-antenna outgoing modes to probe-incoming modes:

|  Eq.  |   |
| :---: | :-- |
| (3.80) | $G^+\_{\\beta i} = \\tfrac{1}{2}\\,e^{im\\phi}\\,d^n\_{\\mu m}(\\theta)\\,e^{i\\mu\\chi}\\,C^{sn(3)}\_{\\sigma\\mu\\nu}(kA)$ |

(index $i$ ↔ $(s,m,n)$ for test antenna; index $\\beta$ ↔ $(\\sigma,\\mu,\\nu)$ for probe)

$\\mathbf{G}^-$ takes probe-outgoing modes to test-antenna-incoming modes:

|  Eq.  |   |
| :---: | :-- |
| (3.84) | $G^-\_{j\\alpha} = \\tfrac{1}{2}\\,C^{\\sigma\\nu(3)}\_{s\\mu n}(-kA)\\,e^{-im\\phi}\\,d^n\_{m\\mu}(-\\theta)\\,e^{-i\\mu\\chi}$ |

Couplings:

|  Eq.  |   |
| :---: | :-- |
| (3.81) | $\\mathbf{a}^p = \\mathbf{G}^+(\\mathbf{b} - \\mathbf{a})$ |
| (3.83) | $\\mathbf{a} = \\mathbf{G}^-(\\mathbf{b}^p - \\mathbf{a}^p)$ |

### Reduced 2-port system

After eliminating $\\mathbf{a}^p, \\mathbf{b}^p, \\mathbf{a}, \\mathbf{b}$:

|  Eq.  |   |
| :---: | :-- |
| (3.85) | $\\begin{bmatrix}M\_{11} & M\_{12}\\\\M\_{21} & M\_{22}\\end{bmatrix}\\begin{bmatrix}v\\\\v\_p\\end{bmatrix} = \\begin{bmatrix}w\_t\\\\w\\end{bmatrix}$ |

with

|  Eq.  |   |
| :---: | :-- |
| (3.86) | $M\_{11} = \\Gamma + \\mathbf{R}\\{\\mathbf{I} - \\mathbf{G}^-(\\mathbf{S}^p - \\mathbf{I}^p)\\mathbf{G}^+(\\mathbf{S} - \\mathbf{I})\\}^{-1}\\mathbf{G}^-(\\mathbf{S}^p - \\mathbf{I}^p)\\mathbf{G}^+\\mathbf{T}$ |
| (3.87) | $M\_{12} = \\mathbf{R}\\{\\mathbf{I} - \\mathbf{G}^-(\\mathbf{S}^p - \\mathbf{I}^p)\\mathbf{G}^+(\\mathbf{S} - \\mathbf{I})\\}^{-1}\\mathbf{G}^-\\mathbf{T}^p$ |
| (3.88) | $M\_{21} = \\mathbf{R}^p\\{\\mathbf{I}^p - \\mathbf{G}^+(\\mathbf{S} - \\mathbf{I})\\mathbf{G}^-(\\mathbf{S}^p - \\mathbf{I}^p)\\}^{-1}\\mathbf{G}^+\\mathbf{T}$ |
| (3.89) | $M\_{22} = \\Gamma\_p + \\mathbf{R}^p\\{\\mathbf{I}^p - \\mathbf{G}^+(\\mathbf{S} - \\mathbf{I})\\mathbf{G}^-(\\mathbf{S}^p - \\mathbf{I}^p)\\}^{-1}\\mathbf{G}^+(\\mathbf{S} - \\mathbf{I})\\mathbf{G}^-\\mathbf{T}^p$ |

### Full transmission formula with multi-reflection (including generator and load)

Adding generator ($v = v\_g + \\Gamma\_g w\_t$) and load ($v\_p = \\Gamma\_l w$):

|  Eq.  |   |
| :---: | :-- |
| (3.90) | $\\boxed{\\;w = \\frac{M\_{21}}{(1 - M\_{11}\\Gamma\_g)(1 - M\_{22}\\Gamma\_l) - M\_{21}M\_{12}\\Gamma\_g\\Gamma\_l}\\,v\_g\\;}$ |

### Neumann series expansion

|  Eq.  |   |
| :---: | :-- |
| (3.91) | $\\{\\mathbf{I} - \\mathbf{G}^-(\\mathbf{S}^p - \\mathbf{I}^p)\\mathbf{G}^+(\\mathbf{S} - \\mathbf{I})\\}^{-1} = \\mathbf{I} + \\{\\mathbf{G}^-(\\mathbf{S}^p - \\mathbf{I}^p)\\mathbf{G}^+(\\mathbf{S} - \\mathbf{I})\\} + \\{\\mathbf{G}^-(\\mathbf{S}^p - \\mathbf{I}^p)\\mathbf{G}^+(\\mathbf{S} - \\mathbf{I})\\}^2 + \\dots$ |

- First term ($\\mathbf{I}$): direct path (recovers the no-multi-reflection (3.10))
- Second term: triply-scattered field (probe → test antenna → probe)
- Higher orders: multiple-bounce paths
- Series converges for $A$ larger than some finite distance; → 0 as $kA \\to \\infty$

### Empty space case

If no probe present, $\\mathbf{S}^p = \\mathbf{I}^p$ → $M\_{11} = \\Gamma$ (only test-antenna reflection survives, as expected).

### Mitigation: microwave isolator

Inserting a matched isolator between probe and load makes the back-scattering ~10 dB lower than for short-circuited probe. The flow-graph node $v\_p$ vanishes, simplifying the system (Fig. 3.9).

---

## 11. Implementation Checklist for Python Code

### Geometry and rotation cascade
1. **Coordinate convention.** Test antenna at origin (unprimed); probe at $(A, \\theta, \\phi)$ in unprimed, with its own primed frame; probe points to origin. Probe roll $\\chi$ about probe $\\hat z'$.
2. **Euler-angle pipeline (3.6).** Implement as four sequential operations matching Appendix A2: phase $e^{im\\phi}$ → mix via $d^n\_{\\mu m}(\\theta)$ → phase $e^{i\\mu\\chi}$ → mix via $C^{sn(3)}\_{\\sigma\\mu\\nu}(kA)$. **Critical:** the printed Euler tuple is $(\\chi, \\theta, \\phi)$ but the operations apply in order $\\phi \\to \\theta \\to \\chi$.

### Core transmission formula (3.10)
3. **Forward synthesis test.** Given test antenna $T\_{smn}$ (e.g., from Section 2.3.4 dipole matrices), probe $R^p\_{\\sigma\\mu\\nu}$ (e.g., from Eq. 2.154), and scan parameters $(A, \\chi, \\theta, \\phi)$, compute $w$ via (3.10). Compare to:
 - Direct evaluation of (3.43) for a Hertzian-dipole probe (electric or magnetic)
 - Direct evaluation of (3.45)–(3.48) on-axis
 - Far-field Friis (3.67) at very large $kA$
4. **Order-of-summation indifference.** $(s, m, n)$ and $(\\sigma, \\mu, \\nu)$ can be summed in any order; verify with summation-axis permutation tests.
5. **Index range.** $s, \\sigma \\in \\{1, 2\\}$; $n, \\nu \\ge 1$; $m \\in [-n, n]$; $\\mu \\in [-\\nu, \\nu]$. Test that values outside these ranges are explicitly zero (catches off-by-one bugs).

### Probe response constants (3.26)
6. **Precompute $P\_{s\\mu n}(kA)$ once.** Build a lookup table for all $(s, \\mu, n)$ at the measurement distance. Verify symmetry (3.27): `P[s, -1, n] == (-1)**(s+1) * P[s, 1, n]`.
7. **Vanishing for $|\\mu| \\ne 1$ probes.** A linearly polarized $\\mu = \\pm 1$ probe must satisfy $P\_{s\\mu n}(kA) = 0$ for $\\mu \\notin \\{-1, +1\\}$.
8. **Reduced formula equivalence.** For a $\\mu = \\pm 1$ probe, verify that (3.10) gives the same result as the reduced form (3.28) — and that (3.28) gives the same as (3.32) via the $\\vec{K}\_{smn}$ identity.

### Hertzian dipole probe regression
9. **Electric dipole response (3.39).** Closed-form: $P^e\_{s1n}(kA) = (\\sqrt{6}/8)\\,i^{-s}\\sqrt{2n+1}\\,R^{(3)}\_{sn}(kA)$. Test for $n = 1\\dots 5$ at $kA = 5, 10, 20$.
10. **Magnetic dipole response (3.41).** Note the $R^{(3)}\_{3-s,n}$ (index-swapped) and the $i^s$ (not $i^{-s}$).
11. **On-axis special cases.** Equations (3.45)–(3.48): probe-received signal equals $E\_\\theta, E\_\\phi, H\_\\theta, H\_\\phi$ scaled by $\\sqrt{6\\pi\\eta}/(2k)$ or $\\sqrt{6\\pi}/(2k\\sqrt\\eta)$. Direct integration check.
12. **Far-field gain reconstruction (3.55).** $G = (8/(3|v|^2))(|W^e\_0|^2 + |W^e\_{\\pi/2}|^2)$.

### Direction-of-transmission equivalence
13. **Same algorithm both ways.** Verify that the inversion algorithm for (3.10) (test transmits) and (3.20) (test receives, via adjoint) returns the **same** physical pattern (when reciprocal). Critical regression: build a reciprocal test antenna, simulate both forward and reverse measurements, recover $T\_{smn}$ in each — both should match the input.
14. **Adjoint substitution (3.18, 3.19).** `T^p[σ, -μ, ν] = (-1)**μ * R^{p'}[σ, μ, ν]` and `R[s, -m, n] = (-1)**m * T'[s, m, n]`.
15. **Reciprocal collapse.** When both antennas are reciprocal, (3.10) and (3.20) must yield bit-identical results given the same scan geometry and excitation.

### Multiple reflections (§3.4)
16. **No-reflection limit.** Setting $\\mathbf{S}^p = \\mathbf{I}^p$ and $\\mathbf{S} = \\mathbf{I}$ in (3.86)–(3.89) collapses $M\_{21}$ to the single-reflection $\\mathbf{R}^p\\mathbf{G}^+\\mathbf{T}$, recovering (3.10) (after the appropriate generator/load handling).
17. **Neumann series convergence.** For a chosen $A$, verify that successive terms in (3.91) decay; report the largest term as a "multi-reflection figure of merit."
18. **Empty-space limit.** $M\_{11} = \\Gamma$ when $\\mathbf{S}^p \\to \\mathbf{I}^p$.

### Numerical sanity
19. **Far-field limit (Friis).** Take $kA = 1000$ with a Hertzian-dipole probe. The full transmission formula (3.10) should match Friis (3.67) to ~3 significant figures (limited by the $1/(kA)$ corrections).
20. **Probe-scaling invariance.** Multiplying $R^p\_{\\sigma\\mu\\nu}$ by a scalar $\\alpha$ should multiply the predicted $w$ by $\\alpha$. Useful for catching unit-confusion bugs.
21. **Tangential consistency.** For a $\\mu = \\pm 1$ probe with two measurements at $\\chi = 0, \\pi/2$, the inferred field should satisfy $\\vec E \\cdot \\hat r = 0$ on the measurement sphere (no radial component in the tangential frame).
22. **Sign of $A$.** $A > 0$ in (3.10). Use Eq. (A3.13) if you need translation coefficients with $kA < 0$ — never recompute (3.6) for $A < 0$.

---

*End of Chapter 3 reference.*
