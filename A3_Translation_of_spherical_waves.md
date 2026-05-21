# Appendix A3 — Translation of Spherical Waves

**Source:** Hansen, *Spherical Near-Field Antenna Measurements*, Appendix A3 (pp. 355–360).

**Cross-references:**
- Symbols: [A5_List_of_principal_symbols_and_uses.md](A5_List_of_principal_symbols_and_uses.md)
- Chapter 2 (Translation coefficient first introduced): [02_Scattering_matrix_description_of_an_antenna.md](02_Scattering_matrix_description_of_an_antenna.md)
- Appendix A1 ($\vec{F}^{(c)}\_{smn}$, $R^{(c)}\_{sn}$, $h\_n^{(1)}$, $P\_n^m$): [A1_Spherical_wave_functions_notation_and_properties.md](A1_Spherical_wave_functions_notation_and_properties.md)
- Appendix A2 (For combined rotation + axial translation): [A2_Rotation_of_spherical_waves.md](A2_Rotation_of_spherical_waves.md)

> **Purpose.** Computational reference for axial-translation of spherical waves. The translation coefficient $C^{sn(c)}\_{\sigma\mu\nu}(kA)$ is the most numerically delicate object in the book — overflow is the default unless the $B(J) \to B'(J)$ trick is applied. Every formula is tagged by book number `(A3.N)`.

---

## 1. Setup (§A3.1)

Two right-handed Cartesian frames, initially coincident. The **primed** frame's origin $O'$ is translated a distance $A$ along the **positive $z$-axis** of the unprimed frame.

|  Eq.  |   |
| :---: | :-- |
|  | $A = \\|\\overrightarrow{OO'}\\|$ |

with **$A > 0$** when $O'$ lies on the positive-$z$ side of $O$, and **$A < 0$** otherwise.

Only **axial translation** (along $\hat{z}$) is covered. A general translation in an arbitrary direction is decomposed as:

|  Eq.  |   |
| :---: | :-- |
|  | $\\boxed{\\text{general translation} = (\\text{rotation}) \\to (\\text{axial translation}) \\to (\\text{inverse rotation})}$ |

> See Stein [1] and Cruzan [2] for direct (non-axial) formulations; Hansen prefers the decomposition for numerical efficiency.

---

## 2. The Translation Equations

A wave function defined in the *unprimed* frame is re-expressed in the *primed* frame. The form of the expansion depends on whether the field point is inside or outside the sphere of radius $|A|$ centered at $O'$:

### Near-region $r' < |A|$ (use $c = 1$ basis in primed frame)

|  Eq.  |   |
| :---: | :-- |
| (A3.1) | $\\boxed{\\;\\vec{F}^{(c)}\_{s\\mu n}(r,\\theta,\\phi) = \\sum\_{\\sigma=1}^{2}\\sum\_{\\substack{\\nu=\\|\\mu\\|\\\\ \\nu \\ne 0}}^{\\infty} C^{sn(c)}\_{\\sigma\\mu\\nu}(kA)\\,\\vec{F}^{(1)}\_{\\sigma\\mu\\nu}(r',\\theta',\\phi'),\\quad r' < \\|A\\|\\;}$ |

### Far-region $r' > |A|$ (preserve $c$ in primed frame)

|  Eq.  |   |
| :---: | :-- |
| (A3.2) | $\\boxed{\\;\\vec{F}^{(c)}\_{s\\mu n}(r,\\theta,\\phi) = \\sum\_{\\sigma=1}^{2}\\sum\_{\\substack{\\nu=\\|\\mu\\|\\\\ \\nu \\ne 0}}^{\\infty} C^{sn(1)}\_{\\sigma\\mu\\nu}(kA)\\,\\vec{F}^{(c)}\_{\\sigma\\mu\\nu}(r',\\theta',\\phi'),\\quad r' > \\|A\\|\\;}$ |

> **Why two regions.** Wave functions of type $c \in \{2,3,4\}$ are singular at their origin. In (A3.1), the primed-frame expansion sees the original source as "external" (it's outside the sphere $r' < |A|$ that contains the new origin), so we use the regular $c=1$ Bessel basis. In (A3.2), the primed-frame expansion sees the original source as "internal" — so we keep the original character $c$.

### Index conventions

| | |
|---|---|
| $s$, $n$ | Indices of original wave function (LHS) |
| $\sigma$, $\nu$ | Indices of new (primed-frame) wave function (summed over) |
| $\mu$ | **Shared** azimuthal index — preserved by axial translation (the $\phi$-dependence $e^{i\mu\phi}$ doesn't change) |
| Upper index $(c)$ on $C$ | Tells you which radial function $z\_p^{(c)}(kA)$ is used to **compute** $C$ itself — see (A3.3) below. In (A3.1) it matches the LHS; in (A3.2) it is $1$ regardless. |

> **Sums.** $\sigma \in \{1, 2\}$ (TE/TM both contribute), $\nu \ge |\mu|$ with $\nu \ge 1$ enforced. **Translation mixes $s, n$ but preserves $\mu$.**

---

## 3. General Formula for $C^{sn(c)}\_{\sigma\mu\nu}(kA)$ (§A3.2)

Bruning & Lo [3]; Larsen [4]. For $kA > 0$:

|  Eq.  |   |
| :---: | :-- |
|  | $\\boxed{\\;C^{sn(c)}\_{\\sigma\\mu\\nu}(kA) = \\sqrt{\\frac{(2n+1)(2\\nu+1)}{n(n+1)\\,\\nu(\\nu+1)}}\\,\\sqrt{\\frac{(\\nu+\\mu)!\\,(n-\\mu)!}{(\\nu-\\mu)!\\,(n+\\mu)!}}\\,(-1)^{\\mu}\\,\\tfrac{1}{2}\\,i^{n-\\nu}\\;}$ |
| (A3.3) | $\\boxed{\\;\\times \\sum\_{p=\\|n-\\nu\\|}^{n+\\nu} i^{-p}\\Bigl[\\delta\_{s\\sigma}\\{n(n+1) + \\nu(\\nu+1) - p(p+1)\\} + \\delta\_{3-s,\\sigma}\\{2i\\mu kA\\}\\Bigr]\\,a(\\mu,\\, n,\\, -\\mu,\\, \\nu,\\, p)\\,z\_p^{(c)}(kA)\\;}$ |

**Parity rule (A3.7):** the linearization coefficient $a(\mu, n, -\mu, \nu, p)$ vanishes when $(n + \nu + p)$ is odd. Therefore the summation in (A3.3) effectively runs only over

|  Eq.  |   |
| :---: | :-- |
|  | $p = \\|n-\\nu\\|,\\, \\|n-\\nu\\|+2,\\, \\|n-\\nu\\|+4,\\, \\dots,\\, n+\\nu-2,\\, n+\\nu$ |

> $C^{sn(c)}\_{\sigma\mu\nu}(kA)$ for $kA < 0$ is obtained from the symmetry relation (A3.13) — never recompute from scratch.

---

## 4. Linearization Coefficients $a(m, n, \mu, \nu, p)$

Defined by the expansion of a product of two unnormalized associated Legendre functions:

|  Eq.  |   |
| :---: | :-- |
| (A3.4) | $P\_n^m(x)\\,P\_\\nu^\\mu(x) = \\sum\_{p=\\|n-\\nu\\|}^{n+\\nu} a(m, n, \\mu, \\nu, p)\\,P\_p^{m+\\mu}(x)$ |

The form needed in (A3.3) is the special case $m \to \mu, \mu \to -\mu$:

|  Eq.  |   |
| :---: | :-- |
| (A3.5) | $P\_n^\\mu(x)\\,P\_\\nu^{-\\mu}(x) = \\sum\_{p=\\|n-\\nu\\|}^{n+\\nu} a(\\mu, n, -\\mu, \\nu, p)\\,P\_p(x)$ |

### Wigner-3j form

|  Eq.  |   |
| :---: | :-- |
| (A3.6) | $a(\\mu, n, -\\mu, \\nu, p) = (2p + 1)\\sqrt{\\frac{(n+\\mu)!\\,(\\nu-\\mu)!}{(n-\\mu)!\\,(\\nu+\\mu)!}}\\,\\begin{pmatrix}n & \\nu & p\\\\ 0 & 0 & 0\\end{pmatrix}\\begin{pmatrix}n & \\nu & p\\\\ \\mu & -\\mu & 0\\end{pmatrix}$ |

### Parity selection rule

|  Eq.  |   |
| :---: | :-- |
| (A3.7) | $\\begin{pmatrix}n & \\nu & p\\\\ 0 & 0 & 0\\end{pmatrix} = 0 \\quad\\text{for } (n+\\nu+p)\\text{ odd}$ |

So $p$ steps in 2's, halving the work.

> **Implementation tip.** Bruning & Lo [6] give recurrences in $p$ for efficient batched computation of $a(\mu, n, -\mu, \nu, p)$ across all valid $p$. For Python: `sympy.physics.wigner.wigner_3j` is exact but slow; a custom recurrence implementation in Cython/numba is the typical production choice. For SciPy users: `scipy.special` does **not** ship 3-j symbols; consider `py3nj` or `sympy`.

---

## 5. Symmetries (§A3.2)

### Elementary block symmetries

|  Eq.  |   |
| :---: | :-- |
| (A3.8) | $C^{1n(c)}\_{1\\mu\\nu}(kA) = C^{2n(c)}\_{2\\mu\\nu}(kA)$ |
| (A3.9) | $C^{1n(c)}\_{2\\mu\\nu}(kA) = C^{2n(c)}\_{1\\mu\\nu}(kA)$ |

> So the $2\times 2$ "s-σ" block has only **two** independent entries: "same" ($s=\sigma$) and "swap" ($s \ne \sigma$).

### Index-permutation symmetries (power-normalized form)

|  Eq.  |   |
| :---: | :-- |
| (A3.10), (A3.11), (A3.12) | $C^{sn(c)}\_{\\sigma\\mu\\nu}(kA) = (-1)^{n+\\nu}\\,C^{s\\nu(c)}\_{\\sigma\\mu n}(kA) = (-1)^{s+\\sigma}\\,(-1)^{n+\\nu}\\,C^{\\sigma\\nu(c)}\_{s,-\\mu,n}(kA) = (-1)^{s+\\sigma}\\,C^{sn(c)}\_{\\sigma,-\\mu,\\nu}(kA)$ |

### Sign of $A$

|  Eq.  |   |
| :---: | :-- |
| (A3.13), (A3.14), (A3.15) | $C^{sn(c)}\_{\\sigma\\mu\\nu}(-kA) = (-1)^{s+\\sigma}\\,(-1)^{n+\\nu}\\,C^{sn(c)}\_{\\sigma\\mu\\nu}(kA) = (-1)^{s+\\sigma}\\,C^{s\\nu(c)}\_{\\sigma\\mu n}(kA) = C^{\\sigma\\nu(c)}\_{s,-\\mu,n}(kA)$ |

### Complex conjugation

|  Eq.  |   |
| :---: | :-- |
| (A3.16) | $\\bigl[C^{sn(3)}\_{\\sigma\\mu\\nu}(kA)\\bigr]^{\*} = (-1)^{s+\\sigma}\\,C^{sn(4)}\_{\\sigma\\mu\\nu}(kA)$ |

> Use to translate **incoming** waves (c=4) without rewriting; in particular, only $c = 1, 3$ ever need explicit evaluation.

---

## 6. Special Cases (§A3.2)

> All special cases below use the convenient binomial-coefficient combination
>
> $\\frac{\\dbinom{-n+\\nu+p}{(-n+\\nu+p)/2}\\,\\dbinom{n-\\nu+p}{(n-\\nu+p)/2}\\,\\dbinom{n+\\nu-p}{(n+\\nu-p)/2}}{\\dbinom{n+\\nu+p}{(n+\\nu+p)/2}}$
>
> arising from explicit 3-j evaluation. The four binomials are all $\binom{J}{J/2}$ with even $J$ — replace each by $B'(J) = \binom{J}{J/2}\,2^{-J}$ to avoid overflow (see §8).

### 1. $\mu = 0$, $\nu \ge 1$  — Eq. (A3.17)

|  Eq.  |   |
| :---: | :-- |
|  | $C^{sn(3)}\_{\\sigma 0 \\nu}(kA) = \\delta\_{s\\sigma}\\,\\tfrac{1}{2}\\,i^{n}\\sqrt{\\frac{2n+1}{n(n+1)}}\\,i^{-\\nu}\\sqrt{\\frac{2\\nu+1}{\\nu(\\nu+1)}}$ |
|  | $\\times \\sum\_{p=\\|n-\\nu\\|,\\,\\text{step }2}^{n+\\nu}\\Biggl[\\frac{n(n+1)+\\nu(\\nu+1)-p(p+1)}{n+\\nu+p+1}\\cdot \\frac{\\dbinom{-n+\\nu+p}{(-n+\\nu+p)/2}\\,\\dbinom{n-\\nu+p}{(n-\\nu+p)/2}\\,\\dbinom{n+\\nu-p}{(n+\\nu-p)/2}}{\\dbinom{n+\\nu+p}{(n+\\nu+p)/2}}\\,i^{-p}\\,(2p+1)\\,h\_p^{(1)}(kA)\\Biggr]$ |

Note: the $\delta\_{s\sigma}$ prefactor means **only diagonal blocks ($s = \sigma$) survive** for $\mu = 0$.

### 2. $\mu = 0$, $\nu = 1$  — Eq. (A3.18)

|  Eq.  |   |
| :---: | :-- |
| (A3.18) | $\\boxed{\\;C^{sn(3)}\_{\\sigma 0 1}(kA) = \\delta\_{s\\sigma}\\,\\sqrt{\\tfrac{3}{2}}\\,\\sqrt{n(n+1)(2n+1)}\\,\\frac{h\_n^{(1)}(kA)}{kA}\\;}$ |

### 3. $\mu = 1$, $\nu \ge 1$  — Eq. (A3.19)

|  Eq.  |   |
| :---: | :-- |
|  | $C^{sn(3)}\_{\\sigma 1 \\nu}(kA) = \\tfrac{1}{4}\\,i^{n}\\,\\frac{\\sqrt{2n+1}}{n(n+1)}\\,i^{-\\nu}\\,\\frac{\\sqrt{2\\nu+1}}{\\nu(\\nu+1)}$ |
|  | $\\times \\sum\_{p=\\|n-\\nu\\|,\\,\\text{step }2}^{n+\\nu}\\Biggl[\\Bigl\\{\\delta\_{s\\sigma}\\,\\frac{(n(n+1)+\\nu(\\nu+1)-p(p+1))^2}{n+\\nu+p+1} + \\delta\_{3-s,\\sigma}\\,2ikA\\,\\frac{n(n+1)+\\nu(\\nu+1)-p(p+1)}{n+\\nu+p+1}\\Bigr\\}$ |
| (A3.19) | $\\times\\,\\frac{\\dbinom{-n+\\nu+p}{(-n+\\nu+p)/2}\\,\\dbinom{n-\\nu+p}{(n-\\nu+p)/2}\\,\\dbinom{n+\\nu-p}{(n+\\nu-p)/2}}{\\dbinom{n+\\nu+p}{(n+\\nu+p)/2}}\\,i^{-p}\\,(2p+1)\\,h\_p^{(1)}(kA)\\Biggr]$ |

Both diagonal ($\delta\_{s\sigma}$) and **off-diagonal** ($\delta\_{3-s,\sigma}$) blocks contribute when $\mu \ne 0$.

### 4. $\mu = 1$, $\nu = 1$  — Eq. (A3.20)

|  Eq.  |   |
| :---: | :-- |
| (A3.20) | $\\boxed{\\;C^{sn(3)}\_{\\sigma 1 1}(kA) = \\tfrac{\\sqrt{3}}{2}\\,\\sqrt{2n+1}\\,\\bigl\\{\\delta\_{s\\sigma}\\,R^{(3)}\_{2n}(kA) + \\delta\_{3-s,\\sigma}\\,i\\,R^{(3)}\_{1n}(kA)\\bigr\\}\\;}$ |

(Using $R^{(c)}\_{sn}$ from Eq. A1.6: $R^{(3)}\_{1n} = h\_n^{(1)}$, $R^{(3)}\_{2n} = (kr)^{-1}\,d/d(kr)\{kr\,h\_n^{(1)}\}$.)

### 5. $\mu = -1$, $\nu \ge 1$  — Eq. (A3.21)

|  Eq.  |   |
| :---: | :-- |
| (A3.21) | $\\boxed{\\;C^{sn(3)}\_{\\sigma,-1,\\nu}(kA) = (-1)^{s+\\sigma}\\,C^{sn(3)}\_{\\sigma 1 \\nu}(kA)\\;}$ |

> So $\mu = -1$ doesn't need a separate computation — apply (A3.21) after evaluating $\mu = +1$. By (A3.12) more generally, any $-\mu$ entry follows from $+\mu$.

---

## 7. Asymptotic Behaviour (Bruning & Lo [6])

For $kA \to \infty$:

|  Eq.  |   |
| :---: | :-- |
| (A3.22) | $C^{sn(3)}\_{\\sigma\\mu\\nu}(kA) = o\\!\\Bigl(\\frac{1}{kA}\\Bigr)\\quad\\text{for } \\mu \\ne \\pm 1$ |
| (A3.23) | $C^{sn(3)}\_{\\sigma 1 \\nu}(kA) = \\frac{\\sqrt{(2n+1)(2\\nu+1)}}{2}\\,i^{\\nu - n - 1}\\,\\frac{e^{ikA}}{kA} + o\\!\\Bigl(\\frac{1}{kA}\\Bigr)$ |
| (A3.24) | $C^{sn(3)}\_{\\sigma,-1,\\nu}(kA) = \\frac{\\sqrt{(2n+1)(2\\nu+1)}}{2}\\,i^{\\nu - n - 1}\\,(-1)^{s+\\sigma}\\,\\frac{e^{ikA}}{kA} + o\\!\\Bigl(\\frac{1}{kA}\\Bigr)$ |

> **Practical implication.** Far-field probe correction depends primarily on the $\mu = \pm 1$ translation coefficients. Other $\mu$ contributions vanish faster than $1/(kA)$ at large probe distance — useful for sanity-checking truncations.

---

## 8. Numerical Computation Note (§A3.2, end)

The combination of four central binomials in (A3.17) and (A3.19) involves $B(J) = \binom{J}{J/2}$ for even $J$. Direct computation overflows for $J \gtrsim 60$ in IEEE double.

**Replace each $B(J)$ with**

|  Eq.  |   |
| :---: | :-- |
| (A3.26) | $B'(J) = \\binom{J}{J/2}\\,2^{-J}$ |

The four-binomial **ratio is invariant** under this scaling (all four factors gain $2^{-J\_i}$ and the four $J\_i$ cancel pairwise — verify: in the (A3.17) combination, the exponents are $(-n+\nu+p) + (n-\nu+p) + (n+\nu-p) - (n+\nu+p) = 0$).

**Range:** $0.025 < B'(J) < 0.5$ for $J = 2, 4, \dots, 1000$.

**Recurrence** (initial $B'(0) = 1$):

|  Eq.  |   |
| :---: | :-- |
| (A3.27) | $\\boxed{\\;B'(J+2) = \\frac{J+1}{J+2}\\,B'(J)\\;}$ |

> So $B'(2) = 1/2$, $B'(4) = 3/8$, $B'(6) = 5/16$, $B'(8) = 35/128$, ... — verify against $\binom{2}{1}/4 = 1/2$, $\binom{4}{2}/16 = 6/16 = 3/8$, $\binom{6}{3}/64 = 20/64 = 5/16$, $\binom{8}{4}/256 = 70/256 = 35/128$. ✓

---

## 9. Implementation Checklist for Python Code

### Setup
1. **Sign of $A$.** Always pass signed $A$ (or $kA$). Use (A3.13) internally if you need $kA < 0$ — never recompute (A3.3) for negative arguments.
2. **Region split.** Implement two separate methods: `translate_near(r' < |A|)` returning $C^{sn(c)}\_{\sigma\mu\nu}$ with the LHS $c$, and `translate_far(r' > |A|)` returning $C^{sn(1)}\_{\sigma\mu\nu}$.

### Linearization coefficients $a(\mu, n, -\mu, \nu, p)$
3. **Parity selection.** $a = 0$ for $(n + \nu + p)$ odd — Eq. (A3.7). Skip those $p$ values.
4. **Wigner-3j check.** Compare your $a$ implementation to the explicit 3-j form (A3.6) using `sympy.physics.wigner.wigner_3j(n, ν, p, 0, 0, 0)` and `wigner_3j(n, ν, p, μ, -μ, 0)`. Match for $n, \nu \le 5$ and all valid $p$.
5. **Recurrence in $p$.** For production code, implement Bruning & Lo's $p$-recurrence (their 1969 report [6]); test it against the 3-j path.
6. **Product identity.** $P\_n^\mu(x)\,P\_\nu^{-\mu}(x) = \sum\_p a(\mu, n, -\mu, \nu, p)\,P\_p(x)$ — Eq. (A3.5). Numerically verify at $x = 0.3, 0.7$ for $n, \nu \le 4$ using `scipy.special.lpmv`.

### Translation coefficient $C^{sn(c)}\_{\sigma\mu\nu}(kA)$
7. **Block structure.** Only two independent $(s,\sigma)$ blocks exist (A3.8, A3.9). Compute only "same" ($s=\sigma$) and "swap" ($s \ne \sigma$); fill the other entries by symmetry.
8. **$\mu = -\mu$ reduction.** Use (A3.12) to get $\mu < 0$ from $\mu > 0$. Equivalent to (A3.21) for the special case $\mu = \pm 1$.
9. **Symmetry $(A3.10)$:** $C^{sn(c)}\_{\sigma\mu\nu}(kA) = (-1)^{n+\nu}\,C^{s\nu(c)}\_{\sigma\mu n}(kA)$. Halves the (n, ν) table.
10. **Conjugation $(A3.16)$:** $C^{sn(4)} = (-1)^{s+\sigma}\,[C^{sn(3)}]^{\*}$. Don't compute $c = 4$ from scratch.
11. **Special-case regression.**
    - $C^{sn(3)}\_{\sigma 0 1}(kA) = \delta\_{s\sigma}\sqrt{3/2}\sqrt{n(n+1)(2n+1)}\,h\_n^{(1)}(kA)/(kA)$ — Eq. (A3.18). Test for $n = 1\dots 5$, $kA = 1, 5, 10$.
    - $C^{sn(3)}\_{\sigma 1 1}(kA) = (\sqrt{3}/2)\sqrt{2n+1}\{\delta\_{s\sigma}R^{(3)}\_{2n}(kA) + \delta\_{3-s,\sigma}\,i\,R^{(3)}\_{1n}(kA)\}$ — Eq. (A3.20).
    - $C^{sn(3)}\_{\sigma,-1,\nu}(kA) = (-1)^{s+\sigma}C^{sn(3)}\_{\sigma 1 \nu}(kA)$ — Eq. (A3.21).
12. **Asymptotic regression.** At $kA = 1000$, $\mu = 1$ entries match $\sqrt{(2n+1)(2\nu+1)}/2 \cdot i^{\nu-n-1}\,e^{ikA}/(kA)$ to 3 significant digits — Eq. (A3.23). Test for $n, \nu \le 5$.
13. **Vanishing-at-large-$kA$ for $|\mu| > 1$.** $|C^{sn(3)}\_{\sigma\mu\nu}(kA)| \cdot kA \to 0$ as $kA \to \infty$ — Eq. (A3.22). Plot vs $kA$ to confirm faster-than-$1/(kA)$ decay.

### Numerical stability
14. **$B'(J)$ recurrence.** Eq. (A3.27): $B'(0) = 1$, $B'(J+2) = ((J+1)/(J+2))\,B'(J)$. Verify $B'(2) = 1/2$, $B'(8) = 35/128 \approx 0.2734$.
15. **Range check.** Assert $0.025 < B'(J) < 0.5$ for $J \in [2, 1000]$.
16. **Combination invariance.** The four-binomial ratio in (A3.17) is unchanged by the $B(J) \to B'(J)$ substitution. Verify symbolically for $(n, \nu, p) = (3, 4, 5)$.
17. **Overflow regression.** Compute $C^{sn(3)}\_{\sigma 0 \nu}$ at $n = \nu = 30$, $kA = 10$ using both naive $\binom{J}{J/2}$ and the $B'(J)$ formulation. The naive form should overflow (or at least lose precision); the $B'$ form should be finite and reasonable.

### Round-trip / consistency
18. **Translation by zero.** $C^{sn(c)}\_{\sigma\mu\nu}(0)$ should reduce to $\delta\_{s\sigma}\,\delta\_{n\nu}$ (identity) — confirms normalization. Caveat: $z\_p^{(3)}(0)$ is singular for $p \ge 1$, so use the analytic $kA \to 0$ limit, not direct evaluation of (A3.3) at $kA = 0$.
19. **Translation then inverse.** Translating by $+A$ and then $-A$ should compose to identity:

    |  Eq.  |   |
    | :---: | :-- |
    |  | $\\sum\_{\\sigma',\\nu'} C^{sn(1)}\_{\\sigma'\\mu\\nu'}(+kA)\\,C^{\\sigma'\\nu'(1)}\_{\\sigma\\mu\\nu}(-kA) = \\delta\_{s\\sigma}\\,\\delta\_{n\\nu}$ |

    Truncate $\nu'$ generously (rule of thumb: $\nu' = N + \lceil kA \rceil$) and check residual.
20. **Wave-function reconstruction.** Synthesize $\vec{F}^{(3)}\_{2,0,1}(r,\theta,\phi)$ centered at $O$ via direct evaluation, then via (A3.2) applied at a $z$-translated frame with chosen $A = \lambda/4$ and $\nu\_{\max} = 20$. Compare values at $r = 5\lambda$, $\theta = \pi/3$, $\phi = \pi/4$ — should match to 4+ digits.

### Combined rotation + translation
21. **Arbitrary translation.** For a general displacement vector $\vec{d}$ (not aligned with $\hat{z}$): (a) rotate frame so $\hat{z}'$ aligns with $\vec{d}$ using $D^n\_{\mu m}(\chi\_o, \theta\_o, \phi\_o)$ from Appendix A2; (b) axial translate by $|\vec{d}|$; (c) inverse rotate. Coefficients are matrix products. Verify against direct displacement using a single test wave function and a small $|\vec{d}|$.
22. **Reciprocity check.** Translation coefficients satisfy reciprocity through the symmetries (A3.10–A3.15); apply both sides to a single mode to confirm self-consistency.

---

*End of Appendix A3 reference.*
