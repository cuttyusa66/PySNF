# Appendix A4 — Data Processing in Antenna Measurements

**Source:** Hansen, *Spherical Near-Field Antenna Measurements*, Appendix A4 (pp. 361–373).
**Cross-reference:** Symbols: [A5_List_of_principal_symbols_and_uses.md](A5_List_of_principal_symbols_and_uses.md)

> **Purpose.** Computational reference for the DFT-based sampling, reconstruction, and integration schemes used throughout the book. Hansen's DFT/IDFT convention is **opposite-signed** to numpy's — flagged with explicit mappings.

---

## 1. Setup and Sampling Bound (§A4.1)

For an antenna whose spherical-wave expansion is truncated at polar index $N$, the field sampled along **any** great circle on the measurement sphere is a periodic band-limited function $f(\\psi)$ of period $2\\pi$ with highest harmonics $e^{\\pm i N\\psi}$. (The bound $N$ is invariant under rotation, so the choice of great circle doesn't matter — Hansen takes the equator $\\psi = \\phi$ without loss of generality.)

### Sampling lower bounds

| Signal | Real or complex | Number of samples in $[0, 2\\pi)$ |
|---|---|---|
| Complex pattern $f(\\phi)$ | complex | $J \\ge 2N + 1$ |
| Power pattern $g(\\phi) = \\|f(\\phi)\\|^2$ | real | $J \\ge 4N + 1$ |

> The bound $J \\ge 2N + 1$ corresponds to a **sampling interval of order one half-power beamwidth** — coarser than typically used for plotting.

---

## 2. The Discrete Fourier Transform (§A4.2.1)

### Hansen's convention

Let $g(m)$, $m = 0, 1, \\dots, J-1$, be a complex sequence of length $J$ and define $\\varepsilon\_J = e^{2\\pi i / J}$ (A4.2).

**Discrete Fourier Transform (DFT):**

|  Eq.  |   |
| :---: | :-- |
| (A4.1) | $\\boxed{\\;G(l) = \\text{DFT}\\{g(m)\\}\_l = \\sum\_{m=0}^{J-1} g(m)\\,\\varepsilon\_J^{ml},\\quad l = 0, 1, \\dots, J-1\\;}$ |

**Inverse Discrete Fourier Transform (IDFT):**

|  Eq.  |   |
| :---: | :-- |
| (A4.4) | $\\boxed{\\;g(m) = \\text{IDFT}\\{G(l)\\}\_m = \\frac{1}{J}\\sum\_{l=0}^{J-1} G(l)\\,\\varepsilon\_J^{-lm},\\quad m = 0, 1, \\dots, J-1\\;}$ |

### ⚠️ Sign convention — critical for Python interop

| Step | Sign in exponent | Normalization | numpy equivalent |
|---|---|---|---|
| Hansen DFT (A4.1) | **$+2\\pi i\\,ml/J$** | none | `J * np.fft.ifft(g)` |
| Hansen IDFT (A4.4) | $-2\\pi i\\,lm/J$ | $1/J$ | `np.fft.fft(G) / J` |

> **Hansen's DFT goes coefficients → samples (synthesis); his IDFT goes samples → coefficients (analysis).** This is the **opposite** of numpy's naming. Build a wrapper:
>
> ```python
> def hansen_dft(g):
>     return len(g) * np.fft.ifft(g)
> def hansen_idft(G):
>     return np.fft.fft(G) / len(G)
> ```
>
> Test by verifying $f(l\\Delta\\phi) = \\sum\_k c'\_k e^{ikl\\Delta\\phi}$ with $c'\_k$ from `hansen_idft`.

### Orthogonality (used to derive A4.4)

**(A4.3)**

```math
\sum_{l=0}^{J-1}\varepsilon_J^{l(s-m)} = \begin{cases}J, & s \equiv m \pmod J \\ 0, & s \not\equiv m \pmod J\end{cases}
```

### Periodic extension

Both $G(l)$ and $g(m)$ are implicitly $J$-periodic:

|  Eq.  |   |
| :---: | :-- |
| (A4.5, A4.6) | $G(l) = G(kJ + l),\\quad g(m) = g(kJ + m),\\quad k \\in \\mathbb{Z}$ |

### Even / odd sequences

A sequence is **even** if $g(m) = g(J-m)$ for $m = 1,\\dots,J-1$ (A4.7).
A sequence is **odd** if $g(0) = 0$ and $g(m) = -g(J-m)$ for $m = 1,\\dots,J-1$ (A4.8). When $J$ is even, an odd sequence must also have $g(J/2) = 0$.

Decomposition $g = g^e + g^o$ (A4.9, A4.10):

```math
g^e(m) = \begin{cases} g(0), & m = 0 \\ \tfrac{1}{2}(g(m) + g(J-m)), & m \ge 1\end{cases}
```

```math
g^o(m) = \begin{cases} 0, & m = 0 \\ \tfrac{1}{2}(g(m) - g(J-m)), & m \ge 1\end{cases}
```

---

## 3. DFT Pair Table (§A4.2.2)

| Sequence | DFT |
|---|---|
| $g(m)$ | $G(l)$ |
| $a\\,g\_1(m) + b\\,g\_2(m)$ | $a\\,G\_1(l) + b\\,G\_2(l)$ |
| $g(-m)$ | $G(-l)$ |
| $g^\*(m)$ | $G^\*(-l)$ |
| $g^\*(-m)$ | $G^\*(l)$ |
| $g^e(m)$ (even part) | $G^e(l)$ |
| $g^o(m)$ (odd part) | $G^o(l)$ |
| $g(m - n)$ (shift in $m$) | $\\varepsilon\_J^{ln}\\,G(l)$ |
| $\\varepsilon\_J^{-km}\\,g(m)$ (modulation) | $G(l - k)$ |
| $\\sum\_{j=0}^{J-1} g\_1(m - j)\\,g\_2(j)$ (convolution) | $G\_1(l)\\,G\_2(l)$ |
| $g\_1(m)\\,g\_2(m)$ | $\\tfrac{1}{J}\\sum\_{j=0}^{J-1} G\_1(l - j)\\,G\_2(j)$ (convolution) |

> **Sanity test.** Convolution-vs-product: pick two short sequences, compute Hansen-convolution by hand, compute Hansen DFT of each, multiply, then Hansen IDFT — should match.

---

## 4. Reconstruction of Periodic Functions (§A4.3)

### Fourier series

|  Eq.  |   |
| :---: | :-- |
| (A4.11, A4.12) | $f(\\phi) = \\sum\_{k=-\\infty}^{\\infty} c\_k\\,e^{ik\\phi},\\quad c\_k = \\frac{1}{2\\pi}\\int\_0^{2\\pi} f(\\phi)\\,e^{-ik\\phi}\\,d\\phi$ |

Sample at $\\phi = l\\Delta\\phi$, $\\Delta\\phi = 2\\pi/J$.

### Aliased coefficients

|  Eq.  |   |
| :---: | :-- |
| (A4.20, A4.27) | $\\boxed{\\;\\bar{c}\_n = \\sum\_{r=-\\infty}^{\\infty} c\_{n+rJ},\\quad n = 0, 1, \\dots, J-1\\;}$ |

Samples expressed via aliased coefficients:

|  Eq.  |   |
| :---: | :-- |
| (A4.21) | $f(l\\Delta\\phi) = \\sum\_{n=0}^{J-1}\\bar{c}\_n\\,\\varepsilon\_J^{nl} = \\text{DFT}\\{\\bar{c}\_n\\}\_l$ |

Solving for $\\bar{c}\_n$:

|  Eq.  |   |
| :---: | :-- |
| (A4.22) | $\\boxed{\\;\\bar{c}\_n = \\frac{1}{J}\\sum\_{l=0}^{J-1} f(l\\Delta\\phi)\\,\\varepsilon\_J^{-ln} = \\text{IDFT}\\{f(l\\Delta\\phi)\\}\_n\\;}$ |

### 4.1. Quasi-band-limited case (§A4.3.2)

If $c\_k$ decays for $|k| > N$ but isn't strictly zero:

**(A4.23)**

```math
c_k \approx \begin{cases} \bar{c}_k, & 0 \le k \le N \\ \bar{c}_{k+J}, & -N \le k < 0 \\ 0, & |k| > N\end{cases} \quad (J \ge 2N + 1)
```

Approximate reconstruction:

|  Eq.  |   |
| :---: | :-- |
| (A4.24) | $f(\\phi) \\approx \\sum\_{k=-N}^{N} c'\_k\\,e^{ik\\phi}$ |

with $c'\_k$ given by the **zero-stuffed** sequence:

|  Eq.  |   |
| :---: | :-- |
| (A4.25) | $\\boxed{\\;\\{c'\_0, c'\_1, \\dots, c'\_N, 0, \\dots, 0, c'\_{-N}, \\dots, c'\_{-1}\\} = \\text{IDFT}\\{f(l\\Delta\\phi) \\mid l = 0, \\dots, J-1\\}\\;}$ |

(Both sequences are length $J$.)

> **Layout.** The IDFT output occupies $J$ bins; the **first $N+1$ bins are $c\_0, c\_1, \\dots, c\_N$**, the **last $N$ bins are $c\_{-N}, c\_{-N+1}, \\dots, c\_{-1}$**, and the middle $J - (2N+1)$ bins are nominally zero. The "aliasing error" lives in those middle bins.

### 4.2. Band-limited case, sufficient sampling (§A4.3.3.1)

If $c\_k = 0$ strictly for $|k| > N$ **and** $J \\ge 2N+1$, (A4.23) becomes exact equality, and:

|  Eq.  |   |
| :---: | :-- |
| (A4.29) | $f(\\phi) = \\sum\_{k=-N}^{N} c\_k\\,e^{ik\\phi}$ |

with $c\_k$ from (A4.30):

|  Eq.  |   |
| :---: | :-- |
|  | $\\{c\_0, c\_1, \\dots, c\_N, 0, \\dots, 0, c\_{-N}, \\dots, c\_{-1}\\} = \\text{IDFT}\\{f(l\\Delta\\phi)\\}$ |

### 4.3. Undersampling (§A4.3.3.2)

If $J \\le 2N$, aliasing fills the middle bins. From (A4.31):

**Case $N < J \\le 2N$:**

```math
c_k = \begin{cases} \bar{c}_k, & k = 0, 1, \dots, J - N - 1 \\ \bar{c}_{k+J}, & k = -(J - N - 1), \dots, -2, -1 \\ \text{undetermined}, & |k| = J - N, J - N + 1, \dots, N\end{cases}
```

**Case $J \\le N$:** all $c\_k$ are undetermined.

### 4.4. Special case $J = 2N$ with $c_N = c_{-N}$ — Eq. (A4.32)

If you know a priori that $c\_N = c\_{-N}$ (e.g. from a symmetry), $J = 2N$ samples suffice:

```math
c_k = \begin{cases} \bar{c}_k, & k = 0, 1, \dots, N - 1 \\ \tfrac{1}{2}\bar{c}_k, & k = N \\ \tfrac{1}{2}\bar{c}_{k+J}, & k = -N \\ \bar{c}_{k+J}, & k = -(N-1), \dots, -2, -1 \end{cases}
```

---

## 5. Special Cases of Reconstruction (§A4.4)

> All operations below assume $J \\ge 2N + 1$ — the input is fully recoverable.

### 5.1. Shifting sample points (§A4.4.2)

To obtain shifted samples $f(l\\Delta\\phi + \\phi\_o)$ from $\\{f(l\\Delta\\phi)\\}$:

1. Compute the aliased coefficients $\\bar{c}\_k$ via IDFT (A4.22).
2. Form the phase sequence $\\bar{e}\_k$ of period $J$ with elements $e^{ik\\phi\_o}$ for $k = -N, \\dots, N$ and zeros in the middle (analogous to A4.37):

   |  Eq.  |   |
   | :---: | :-- |
   |  | $\\{\\bar{e}\_k\\} = \\{e^{i\\cdot 0\\cdot \\phi\_o}, e^{i\\cdot 1\\cdot \\phi\_o}, \\dots, e^{iN\\phi\_o}, 0, \\dots, 0, e^{-iN\\phi\_o}, \\dots, e^{-i\\phi\_o}\\}$ |

3. Term-by-term multiply: $\\bar{e}\_k \\bar{c}\_k$.
4. Forward DFT to get the shifted samples:

   |  Eq.  |   |
   | :---: | :-- |
   | (A4.38) | $\\boxed{\\;\\{f(l\\Delta\\phi + \\phi\_o)\\} = \\text{DFT}\\{\\bar{e}\_k \\bar{c}\_k\\}\\;}$ |

> **Tip.** If $|\\phi\_o| \\gg \\Delta\\phi$, first apply a cyclic shift of the samples (an integer number of $\\Delta\\phi$ steps) to reduce $|\\phi\_o|$ — this is numerically cleaner.

### 5.2. Increased resolution / Whittaker interpolation (§A4.4.3)

To obtain samples at finer spacing $\\Delta'\\phi = 2\\pi/J'$ with $J' > J \\ge 2N + 1$:

1. Compute $\\{\\bar{c}\_k\\}$ of period $J$ via IDFT (A4.22).
2. **Zero-pad in the middle** to length $J'$ (A4.39):

   |  Eq.  |   |
   | :---: | :-- |
   |  | $\\{\\bar{b}\_k\\} = \\{\\bar{c}\_0, \\bar{c}\_1, \\dots, \\bar{c}\_N, \\underbrace{0, 0, \\dots, 0}\_{J' - J \\text{ zeros}}, \\bar{c}\_{-N}, \\dots, \\bar{c}\_{-1}\\}$ |

   (Length $J'$.)
3. DFT to get fine samples:

   |  Eq.  |   |
   | :---: | :-- |
   | (A4.40) | $\\boxed{\\;\\{f(l\\Delta'\\phi)\\} = \\text{DFT}\\{\\bar{b}\_k \\mid k = 0, \\dots, J'-1\\}\\;}$ |

> **Edge case (footnote on p. 370):** If $J$ is even and $N$ is uncertain, the original middle element $\\bar{c}\_{J/2}$ should be **split** into two elements of $\\tfrac{1}{2}\\bar{c}\_{J/2}$ each, with $(J' - J - 1)$ zeros inserted between them — preserves spectral symmetry for non-strict bandlimits.

### 5.3. Noise suppression by zero-stuffing (§A4.4.4)

Used when $J \\gg 2N + 1$ and the samples carry additive noise:

1. IDFT the noisy samples → coefficient sequence $\\bar{c}'\_k$ where all $J$ bins are non-zero (noise has leaked into the middle).
2. **Zero out the middle $J - (2N + 1)$ bins** to form $\\bar{c}''\_k$.
3. DFT back to get denoised samples:

   |  Eq.  |   |
   | :---: | :-- |
   | (A4.42) | $\\{f''(l\\Delta\\phi)\\} = \\text{DFT}\\{\\bar{c}''\_k\\}$ |

> Effective noise reduction ratio scales with $J / (2N+1)$ for white noise. Works because the signal occupies only $2N + 1$ bins by hypothesis; the remaining bins are "noise only."

### 5.4. Decreased resolution (§A4.4.5)

To obtain samples at coarser spacing $\\Delta'\\phi = 2\\pi/J'$ with $J' < J$:

1. IDFT to get $\\{\\bar{c}\_k\\}$ of period $J$, extract the $2N+1$ "real" coefficients $\\{c\_k \\mid k = -N, \\dots, N\\}$ (A4.43).
2. **Re-alias** to a new sequence $\\{\\bar{c}'\_k\\}$ of period $J'$ via Eq. (A4.27) — overlap is expected if $J' < 2N+1$.
3. DFT back:

   |  Eq.  |   |
   | :---: | :-- |
   | (A4.44) | $\\{f(l\\Delta'\\phi)\\} = \\text{DFT}\\{\\bar{c}'\_k\\mid k = 0, \\dots, J'-1\\}$ |

Valid for any $J' \\in [1, J]$, including the trivial $J' = 1$ (the mean value).

---

## 6. Numerical Integration of Periodic Band-Limited Functions (§A4.5)

### 6.1. The trapezoidal rule is exact (§A4.5.1)

For a band-limited $f$ with $2N + 1$ Fourier coefficients:

|  Eq.  |   |
| :---: | :-- |
| (A4.45) | $I\_0 = \\frac{1}{2\\pi}\\int\_0^{2\\pi} f(\\phi)\\,d\\phi = c\_0$ |

**With $J \\ge N + 1$** samples in $[0, 2\\pi)$:

|  Eq.  |   |
| :---: | :-- |
| (A4.46) | $\\boxed{\\;I\_0 = c\_0 = \\frac{1}{J}\\sum\_{l=0}^{J-1} f(l\\Delta\\phi)\\;}$ |

This is the **trapezoidal rule** with equal weights — it is **exact** for periodic band-limited functions. Higher-order rules (Simpson, Romberg) do **not** improve accuracy for this class.

### 6.2. Generalised integrals (§A4.5.2)

For the $2N+1$ moment integrals:

|  Eq.  |   |
| :---: | :-- |
| (A4.47) | $I\_m = \\frac{1}{2\\pi}\\int\_0^{2\\pi} f(\\phi)\\,e^{-im\\phi}\\,d\\phi,\\quad m = -N, \\dots, 0, \\dots, N$ |

the integrand has highest harmonics $e^{\\pm i 2N\\phi}$, so **$J \\ge 2N + 1$** samples suffice (and the trapezoidal rule remains exact):

|  Eq.  |   |
| :---: | :-- |
| (A4.48) | $I\_m = \\frac{1}{J}\\sum\_{l=0}^{J-1} f(l\\Delta\\phi)\\,e^{-iml\\Delta\\phi}$ |

Recognising the sums as an IDFT:

|  Eq.  |   |
| :---: | :-- |
| (A4.49) | $\\boxed{\\;\\{I\_0, I\_1, \\dots, I\_N, I\_{-N}, \\dots, I\_{-1}\\} = \\text{IDFT}\\{f(l\\Delta\\phi) \\mid l = 0, \\dots, J-1\\}\\;}$ |

> **Implementation power-tool.** A single IDFT computes all $2N+1$ moment integrals at once. Excellent for batch-extracting Fourier coefficients of measured probe-output signals.

---

## 7. Implementation Checklist for Python Code

### DFT convention
1. **Sign convention.** Hansen's $\\text{DFT}$ has the **$+i$** exponent and no normalization; numpy's `np.fft.fft` has the $-i$ exponent. Implement explicit wrappers `hansen_dft`/`hansen_idft` mapping to `len(g) * np.fft.ifft(g)` and `np.fft.fft(G) / len(G)` respectively.
2. **Round-trip identity.** `hansen_idft(hansen_dft(g)) == g` to machine precision. Test with random complex `g` of length 8, 16, 17 (odd).
3. **Periodicity.** Verify $G(l) = G(l + J)$ and $g(m) = g(m + J)$ on the output of each transform.
4. **Orthogonality (A4.3).** Numerically verify $\\sum\_l \\varepsilon\_J^{l(s-m)} = J\\,\\delta\_{s \\bmod J,\\,m \\bmod J}$ for a few $(s, m, J)$.
5. **Even/odd decomposition (A4.9, A4.10).** Build a sequence with known even+odd parts; verify `g = g_even + g_odd` and that `g_even`/`g_odd` satisfy (A4.7)/(A4.8).
6. **DFT pair table (§3).** Spot-check at least three pairs: shift theorem, convolution-vs-product, and conjugation $G^\*(-l)$.

### Reconstruction
7. **Aliased coefficient identity (A4.20, A4.22).** For a known $c\_k$ band-limited to $|k| \\le N$, sample at $J \\ge 2N+1$ and verify that `hansen_idft(samples)` recovers $\\{c\_0, c\_1, \\dots, c\_N, 0, \\dots, 0, c\_{-N}, \\dots, c\_{-1}\\}$.
8. **Layout convention.** The IDFT-of-samples output places non-zero coefficients at indices $\\{0, 1, \\dots, N\\}$ and $\\{J-N, \\dots, J-1\\}$. Indices $\\{N+1, \\dots, J-N-1\\}$ should be zero up to numerical noise. Test with $N=4$, $J = 16$.
9. **Reconstruction (A4.24).** Given samples, reconstruct $f(\\phi)$ at arbitrary $\\phi$ (not just sample points) by summing $\\sum\_{k=-N}^{N} c'\_k e^{ik\\phi}$. Compare to ground-truth analytic $f(\\phi)$ at $\\phi = \\pi/7$, etc.
10. **Undersampling regression (A4.31).** With $N = 5$, $J = 8$: only coefficients $c\_0, c\_1, c\_2, c\_{-1}, c\_{-2}$ should be recoverable; $c\_3, c\_4, c\_5, c\_{-3}, c\_{-4}, c\_{-5}$ are aliased into the recoverable bins. Verify the aliasing pattern matches (A4.31).
11. **$J = 2N$ symmetric special case (A4.32).** If $c\_N = c\_{-N}$ is enforced, $J = 2N$ samples should suffice. Construct an explicit example with $N = 4, J = 8, c\_4 = c\_{-4}$, verify recovery via (A4.32). Note the $\\tfrac{1}{2}\\bar{c}$ factors at $k = \\pm N$.

### Sample operations
12. **Shift (A4.38).** Generate samples at $\\{l\\Delta\\phi\\}$, apply the shift recipe, compare to direct evaluation at $\\{l\\Delta\\phi + \\phi\_o\\}$. Try $\\phi\_o = \\Delta\\phi / 3$ (sub-sample), $\\phi\_o = 1.5\\Delta\\phi$ (super-sample), and $\\phi\_o = 100\\Delta\\phi$ (test the cyclic-shift optimization).
13. **Whittaker interpolation (A4.40).** Generate $J = 17$ samples of a known $N = 8$ band-limited function. Zero-pad to $J' = 64$, IDFT, compare to direct evaluation at the fine grid. Should match to ~$10^{-13}$.
14. **Whittaker edge case.** When $J$ is even, verify the $\\tfrac{1}{2}\\bar{c}\_{J/2}$ split prescription gives better continuity than naive zero-padding at one end (compare reconstruction error in two strategies).
15. **Noise suppression (A4.42).** Add white Gaussian noise to samples. Apply zero-stuffing with $J = 8(2N+1)$. Verify residual noise variance reduces by factor $\\approx J/(2N+1) = 8$.
16. **Decreased resolution (A4.44).** Start from $J = 32$ samples of an $N = 5$ pattern. Decimate to $J' = 11 (\\ge 2N + 1)$ — should be exact. Decimate to $J' = 8 (< 2N + 1)$ — should produce aliased samples consistent with the re-aliased $\\bar{c}'\_k$.

### Numerical integration
17. **Trapezoidal exactness (A4.46).** For a band-limited $f$ with $N = 5$, evaluate $\\int\_0^{2\\pi} f / (2\\pi)$ via direct quadrature (scipy `quad`) and via trapezoidal with $J = N + 1 = 6$ samples. Should agree to machine precision.
18. **Higher-order rules don't help.** Repeat (17) with Simpson's rule at $J = N + 1$. **Should be no more accurate** (and may be slightly worse due to a non-band-limited weighting). Counter-intuitive but key.
19. **Batch moments (A4.49).** For an $N = 5$ band-limited $f$, compute all 11 moments $I\_m$ for $m = -5, \\dots, 5$ via a single IDFT of $J = 11$ samples. Cross-check each against a direct quadrature.

### Convention gotchas
20. **Sample indexing.** Hansen samples at $l\\Delta\\phi$ for $l = 0, 1, \\dots, J-1$ — i.e., the first sample is at $\\phi = 0$, **not** $\\phi = \\Delta\\phi/2$. Centered-bin samplers must be shifted before applying these formulas.
21. **Periodic wrap.** Samples at $\\phi = 0$ and $\\phi = 2\\pi$ are the same point; the sample array has $J$ entries (not $J+1$). Double-counting the endpoint is a classic bug.

---

*End of Appendix A4 reference.*
