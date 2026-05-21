# Appendix A1 — Spherical Wave Functions, Notation and Properties

**Source:** Hansen, *Spherical Near-Field Antenna Measurements*, Appendix A1 (pp. 312–342).
**Cross-references:**
- Symbols: [A5_List_of_principal_symbols_and_uses.md](A5_List_of_principal_symbols_and_uses.md)
- Chapter 2: [02_Scattering_matrix_description_of_an_antenna.md](02_Scattering_matrix_description_of_an_antenna.md)

> **Purpose.** Computational reference for Python code that builds spherical-wave expansions. Every formula is tagged with its book number `(A1.N)` so unit tests can cite it. Explicit low-order tables are included as concrete validation targets.

> **OCR caveat.** The published tables render negative exponents ambiguously (e.g. `x⁻²` printed as `x²` in places). Where this affects a formula, the corrected form is noted; trust the algebraic identities in §1 over any inconsistent literal.

---

## 1. Global Conventions (recap)

| Convention | Value |
|---|---|
| Time factor | $e^{-i\omega t}$ |
| $k$ | $\omega\sqrt{\mu\varepsilon} = 2\pi/\lambda$ |
| $\eta$ | $\sqrt{\varepsilon/\mu}$ (specific admittance) |
| $\zeta$ | $\sqrt{\mu/\varepsilon}$ (specific impedance) |
| $s$ | $s \in \{1, 2\}$. $\vec{F}^{(c)}_{1mn}$ purely transverse; $\vec{F}^{(c)}_{2mn}$ has a radial component. In coefficient $Q^{(c)}_{smn}$: $s=1 \to$ TE, $s=2 \to$ TM. |
| $n$ | $1, 2, 3, \dots, \infty$ (truncate at $N$). |
| $m$ | $-n, \dots, n$. |
| $c$ | 1: $j_n$; 2: $n_n$; 3: $h_n^{(1)}$ outgoing; 4: $h_n^{(2)}$ incoming. |

---

## 2. The General Spherical Wave Expansion (§A1.1.1)

In a source-free region bounded by spheres centred at the origin:

|  Eq.  |   |
| :---: | :-- |
| (A1.1) | $\vec{E}(r,\theta,\phi) = \frac{k}{\sqrt{\eta}}\sum_{csmn} Q^{(c)}_{smn}\,\vec{F}^{(c)}_{smn}(r,\theta,\phi)$ |
| (A1.2) | $\vec{H}(r,\theta,\phi) = -ik\sqrt{\eta}\sum_{csmn} Q^{(c)}_{smn}\,\vec{F}^{(c)}_{3-s,m,n}(r,\theta,\phi)$ |

- $\vec{F}^{(c)}_{smn}$ are **dimensionless**.
- $[Q^{(c)}_{smn}] = \text{W}^{1/2}$.
- $Q^{(c)}_{1mn}$ → TE part; $Q^{(c)}_{2mn}$ → TM part.

Outgoing-wave radiated power:

|  Eq.  |   |
| :---: | :-- |
| (A1.3) | $P = \tfrac{1}{2}\sum_{smn}\|Q^{(3)}_{smn}\|^2 \quad\text{[W]}$ |

---

## 3. Single-Index Convention (§A1.1.2)

### Forward: $(s,m,n) \to j$

|  Eq.  |   |
| :---: | :-- |
|  | $\boxed{\;j = 2\{n(n+1) + m - 1\} + s\;}$ |

with $j = 1, 2, \dots, J$ and $J = 2N(N+2)$.

Two equivalent triple summation orders:

|  Eq.  |   |
| :---: | :-- |
| (A1.4) | $\sum_{smn} = \sum_{s=1}^{2}\sum_{m=-N}^{N}\sum_{\substack{n=\|m\|\\ n \ne 0}}^{N} = \sum_{s=1}^{2}\sum_{n=1}^{N}\sum_{m=-n}^{n}$ |

### Inverse: $j \to (s,m,n)$ — three-step algorithm

1. $s = \begin{cases} 1, & j \text{ odd} \\ 2, & j \text{ even}\end{cases}$
2. $n = \left\lfloor \sqrt{(j-s)/2 + 1\,} \right\rfloor$
3. $m = (j-s)/2 + 1 - n(n+1)$

### Validation table (small $j$)

| $j$ | $s$ | $m$ | $n$ |
|---|---|---|---|
| 1 | 1 | -1 | 1 |
| 2 | 2 | -1 | 1 |
| 3 | 1 |  0 | 1 |
| 4 | 2 |  0 | 1 |
| 5 | 1 | +1 | 1 |
| 6 | 2 | +1 | 1 |
| 7 | 1 | -2 | 2 |
| 8 | 2 | -2 | 2 |
| 9 | 1 | -1 | 2 |
| 10 | 2 | -1 | 2 |
| 11 | 1 |  0 | 2 |
| 12 | 2 |  0 | 2 |
| 24 | 2 | +2 | 3 |

---

## 4. Far-Field Pattern Functions (§A1.1.3)

|  Eq.  |   |
| :---: | :-- |
| (A1.5) | $\vec{K}_{smn}(\theta,\phi) = \lim_{kr\to\infty}\left\{\sqrt{4\pi}\,\frac{kr}{e^{ikr}}\,\vec{F}^{(3)}_{smn}(r,\theta,\phi)\right\}$ |

The $\sqrt{4\pi}$ factor simplifies the gain/directivity formulas.

---

## 5. Radial Functions $R^{(c)}_{sn}(kr)$ (§A1.2.1)

### Definition

|  Eq.  |   |
| :---: | :-- |
| (A1.6) | $R^{(c)}_{sn}(kr) = \begin{cases} z_n^{(c)}(kr), & s = 1 \\[6pt] \dfrac{1}{kr}\dfrac{d}{d(kr)}\{kr\,z_n^{(c)}(kr)\}, & s = 2 \end{cases}$ |

with

| $c$ | $z_n^{(c)}$ | Description |
|---|---|---|
| 1 | $j_n(kr)$ | Spherical Bessel |
| 2 | $n_n(kr)$ | Spherical Neumann |
| 3 | $h_n^{(1)}(kr) = j_n + i\,n_n$ | Outgoing Hankel |
| 4 | $h_n^{(2)}(kr) = j_n - i\,n_n$ | Incoming Hankel |

### Recurrence relations (let $x = kr$)

|  Eq.  |   |
| :---: | :-- |
| (A1.7) | $\frac{z_n^{(c)}}{x} = \frac{1}{2n+1}\bigl\{z_{n-1}^{(c)} + z_{n+1}^{(c)}\bigr\}$ |
| (A1.8), (A1.9), (A1.10) | $\frac{1}{x}\frac{d}{dx}\bigl\{x\,z_n^{(c)}\bigr\} = z_{n-1}^{(c)} - n\,\frac{z_n^{(c)}}{x} = (n+1)\,\frac{z_n^{(c)}}{x} - z_{n+1}^{(c)} = \frac{1}{2n+1}\bigl\{(n+1)\,z_{n-1}^{(c)} - n\,z_{n+1}^{(c)}\bigr\}$ |

### Wronskian

|  Eq.  |   |
| :---: | :-- |
| (A1.11) | $R^{(c)}_{sn}(kr)\,R^{(\gamma)}_{3-s,n}(kr) - R^{(\gamma)}_{sn}(kr)\,R^{(c)}_{3-s,n}(kr) = -(-1)^s\,\frac{A^{(c,\gamma)}}{(kr)^2}$ |

The $A^{(c,\gamma)}$ table:

| $A^{(c,\gamma)}$ | $\gamma=1$ | $\gamma=2$ | $\gamma=3$ | $\gamma=4$ |
|---|---|---|---|---|
| $c=1$ | 0  | 1  | $i$  | $-i$ |
| $c=2$ | $-1$ | 0 | $-1$ | $-1$ |
| $c=3$ | $-i$ | 1 | 0 | $-2i$ |
| $c=4$ | $i$ | 1 | $2i$ | 0 |

**Spot check:** $s=1, c=1, \gamma=2$ gives $R^{(1)}_{1n}R^{(2)}_{2n} - R^{(2)}_{1n}R^{(1)}_{2n} = +1/(kr)^2$ (matches Eq. (2.52) in Chapter 2).

### Cross-product identity

|  Eq.  |   |
| :---: | :-- |
| (A1.13) | $R^{(1)}_{1n}(kr)\,R^{(1)}_{2n}(kr) + R^{(2)}_{1n}(kr)\,R^{(2)}_{2n}(kr) = \Bigl\{\frac{1}{kr} + \frac{1}{2}\frac{d}{d(kr)}\Bigr\}\|h_n^{(1)}(kr)\|^2$ |

### Asymptotic ($kr \to \infty$, $kr \gg n$)

|  Eq.  |   |
| :---: | :-- |
| (A1.14) | $R^{(3)}_{1n}(kr) = z_n^{(3)}(kr) \to (-i)^{n+1}\,\frac{e^{ikr}}{kr}$ |
| (A1.15) | $R^{(4)}_{1n}(kr) = z_n^{(4)}(kr) \to i^{n+1}\,\frac{e^{-ikr}}{kr}$ |
| (A1.16) | $R^{(3)}_{2n}(kr) = \frac{1}{kr}\frac{d}{d(kr)}\{kr\,z_n^{(3)}\} \to (-i)^{n}\,\frac{e^{ikr}}{kr}$ |
| (A1.17) | $R^{(4)}_{2n}(kr) = \frac{1}{kr}\frac{d}{d(kr)}\{kr\,z_n^{(4)}\} \to i^{n}\,\frac{e^{-ikr}}{kr}$ |

### Special values at $x = 0$

|  Eq.  |   |
| :---: | :-- |
| (A1.18) | $j_n(0) = 0,\quad n \ge 1$ |
| (A1.19) | $\lim_{x\to 0}\frac{j_n(x)}{x} = \begin{cases} 1/3, & n = 1 \\ 0, & n > 1\end{cases}$ |
| (A1.20) | $\lim_{x\to 0}\frac{1}{x}\frac{d}{dx}\{x\,j_n(x)\} = \begin{cases} 2/3, & n = 1 \\ 0, & n > 1\end{cases}$ |

### Low-order $j_n(x)$ ($s = 1$) — Eq. (A1.21)

| $n$ | $j_n(x)$ |
|---|---|
| 0 | $x^{-1}\sin x$ |
| 1 | $x^{-1}\{-\cos x + x^{-1}\sin x\}$ |
| 2 | $x^{-1}\{-3x^{-1}\cos x + (-1 + 3x^{-2})\sin x\}$ |
| 3 | $x^{-1}\{(1 - 15x^{-2})\cos x + (-6x^{-1} + 15x^{-3})\sin x\}$ |
| 4 | $x^{-1}\{(10x^{-1} - 105x^{-3})\cos x + (1 - 45x^{-2} + 105x^{-4})\sin x\}$ |
| 5 | $x^{-1}\{(-1 + 105x^{-2} - 945x^{-4})\cos x + (15x^{-1} - 420x^{-3} + 945x^{-5})\sin x\}$ |

### Low-order $h_n^{(1)}(x)$ ($s = 1$) — Eq. (A1.22)

| $n$ | $h_n^{(1)}(x)$ |
|---|---|
| 0 | $x^{-1}e^{ix}\{-i\}$ |
| 1 | $x^{-1}e^{ix}\{-1 - ix^{-1}\}$ |
| 2 | $x^{-1}e^{ix}\{-3x^{-1} + i(1 - 3x^{-2})\}$ |
| 3 | $x^{-1}e^{ix}\{(1 - 15x^{-2}) + i(6x^{-1} - 15x^{-3})\}$ |
| 4 | $x^{-1}e^{ix}\{(10x^{-1} - 105x^{-3}) + i(-1 + 45x^{-2} - 105x^{-4})\}$ |
| 5 | $x^{-1}e^{ix}\{(-1 + 105x^{-2} - 945x^{-4}) + i(-15x^{-1} + 420x^{-3} - 945x^{-5})\}$ |

### Low-order $R^{(1)}_{2n}(x) = \dfrac{1}{x}\dfrac{d}{dx}\{x\,j_n(x)\}$ — Eq. (A1.23)

| $n$ | $R^{(1)}_{2n}(x)$ |
|---|---|
| 0 | $x^{-1}\cos x$ |
| 1 | $x^{-1}\{x^{-1}\cos x + (1 - x^{-2})\sin x\}$ |
| 2 | $x^{-1}\{(-1 + 6x^{-2})\cos x + (3x^{-1} - 6x^{-3})\sin x\}$ |
| 3 | $x^{-1}\{(-6x^{-1} + 45x^{-3})\cos x + (-1 + 21x^{-2} - 45x^{-4})\sin x\}$ |
| 4 | $x^{-1}\{(1 - 55x^{-2} + 420x^{-4})\cos x + (-10x^{-1} + 195x^{-3} - 420x^{-5})\sin x\}$ |
| 5 | $x^{-1}\{(15x^{-1} - 630x^{-3} + 4725x^{-5})\cos x + (1 - 120x^{-2} + 2205x^{-4} - 4725x^{-6})\sin x\}$ |

### Low-order $R^{(3)}_{2n}(x) = \dfrac{1}{x}\dfrac{d}{dx}\{x\,h_n^{(1)}(x)\}$ — Eq. (A1.24)

| $n$ | $R^{(3)}_{2n}(x)$ |
|---|---|
| 0 | $x^{-1}e^{ix}$ |
| 1 | $x^{-1}e^{ix}\{x^{-1} + i(-1 + x^{-2})\}$ |
| 2 | $x^{-1}e^{ix}\{(-1 + 6x^{-2}) + i(-3x^{-1} + 6x^{-3})\}$ |
| 3 | $x^{-1}e^{ix}\{(-6x^{-1} + 45x^{-3}) + i(1 - 21x^{-2} + 45x^{-4})\}$ |
| 4 | $x^{-1}e^{ix}\{(1 - 55x^{-2} + 420x^{-4}) + i(10x^{-1} - 195x^{-3} + 420x^{-5})\}$ |
| 5 | $x^{-1}e^{ix}\{(15x^{-1} - 630x^{-3} + 4725x^{-5}) + i(-1 + 120x^{-2} - 2205x^{-4} + 4725x^{-6})\}$ |

---

## 6. Angular Functions (§A1.2.2)

### Definitions

Normalized associated Legendre function (Belousov [2] convention):

|  Eq.  |   |
| :---: | :-- |
| (A1.25) | $\bar{P}_n^m(\cos\theta) = \sqrt{\frac{2n+1}{2}\cdot\frac{(n-m)!}{(n+m)!}}\,P_n^m(\cos\theta)$ |

Unnormalized Stratton/Hansen convention:

|  Eq.  |   |
| :---: | :-- |
|  | $P_n^m(\cos\theta) = (\sin\theta)^m\,\frac{d^m P_n(\cos\theta)}{d(\cos\theta)^m}$ |
|  | $P_n(\cos\theta) = \frac{1}{2^n n!}\,\frac{d^n}{d(\cos\theta)^n}(\cos^2\theta - 1)^n$ |

> **Sign convention warning.** Abramowitz & Stegun [3] and Harrington [4] include an extra $(-1)^m$ — Hansen and Stratton [1] do **not**. **scipy's `scipy.special.lpmv` follows the Hansen/Stratton convention; `scipy.special.sph_harm` does NOT.** Always check before use.

Hansen uses only $|m| \ge 0$ in $\bar{P}_n^{|m|}$.

### Orthogonality

|  Eq.  |   |
| :---: | :-- |
| (A1.26) | $\int_{-1}^{1} P_k^m(\mu)\,P_n^m(\mu)\,d\mu = \frac{2}{2n+1}\,\frac{(n+m)!}{(n-m)!}\,\delta_{nk}$ |
| (A1.27) | $\int_{-1}^{1} \bar{P}_k^m(\mu)\,\bar{P}_n^m(\mu)\,d\mu = \delta_{nk}$ |
| (A1.28) | $\int_{-1}^{1} P_n^m(\mu)\,P_n^k(\mu)\,\sin^{-2}\theta\,d\mu = \frac{1}{m}\,\frac{(n+m)!}{(n-m)!}\,\delta_{mk},\quad (m,k) \ne (0,0)$ |
| (A1.29) | $\int_{-1}^{1}\Bigl\{\frac{dP_n^m}{d\theta}\frac{dP_k^m}{d\theta} + \frac{m^2}{\sin^2\theta}P_n^m P_k^m\Bigr\}d\mu = \frac{2}{2n+1}\frac{(n+m)!}{(n-m)!}\,n(n+1)\,\delta_{nk}$ |
| (A1.30) | $\int_{-1}^{1}\Bigl\{\frac{P_n^m}{\sin\theta}\frac{dP_k^m}{d\theta} + \frac{P_k^m}{\sin\theta}\frac{dP_n^m}{d\theta}\Bigr\}\sin\theta\,d\mu = 0$ |

### Recurrence relations

|  Eq.  |   |
| :---: | :-- |
| (A1.31) | $(n-m+1)P_{n+1}^m - (2n+1)\cos\theta\,P_n^m + (n+m)P_{n-1}^m = 0$ |
| (A1.32) | $\sin\theta\,P_n^{m+1} - 2m\cos\theta\,P_n^m + (n+m)(n-m+1)\sin\theta\,P_n^{m-1} = 0$ |
| (A1.33) | $P_n^n(\cos\theta) - (2n-1)\sin\theta\,P_{n-1}^{n-1}(\cos\theta) = 0$ |
| (A1.34a) | $\frac{m\,P_n^m(\cos\theta)}{\sin\theta} = \begin{cases} 0, & m = 0 \\ \tfrac{1}{2}\cos\theta\{(n-m+1)(n+m)P_n^{m-1} + P_n^{m+1}\} + m\sin\theta\,P_n^m, & m > 0\end{cases}$ |
| (A1.34b) | $\frac{dP_n^m(\cos\theta)}{d\theta} = \begin{cases} -P_n^1(\cos\theta), & m = 0 \\ \tfrac{1}{2}\{(n-m+1)(n+m)P_n^{m-1} - P_n^{m+1}\}, & m > 0\end{cases}$ |

> **OCR note.** Eq. (A1.34a) as printed includes a curious $+m\sin\theta\,P_n^m$ tail; check carefully against `scipy.special.lpmv` before relying on it for code.

### Fourier expansion in $\theta$

|  Eq.  |   |
| :---: | :-- |
| (A1.35a) | $P_n^m(\cos\theta) = \sum_{m'=-n}^{n} c_{m'}\,e^{im'\theta}$ |

with $c_{m'} = 0$ when $(m'+n)$ is odd. Recurrence (see Eq. A2.37):

|  Eq.  |   |
| :---: | :-- |
| (A1.35b) | $(n+m'+2)(n-m'-1)c_{m'+2} - 2(n^2 - m'^2 + n - 2m^2)\,c_{m'} + (n+m'-1)(n-m'+2)\,c_{m'-2} = 0$ |

Initial values:

|  Eq.  |   |
| :---: | :-- |
|  | $c_n = \Bigl(-i\,\frac{m}{\|m\|}\Bigr)^m\,\frac{1}{2^{2n}}\,\frac{(2n)!}{(n-m)!\,n!},\quad c_{n-2} = \frac{n - 2m^2}{2n - 1}\,c_n$ |

For $m = 0$ the recurrence simplifies to:

|  Eq.  |   |
| :---: | :-- |
| (A1.35c) | $(n+m')(n-m'+1)c_{m'} - (n+m'-1)(n-m'+2)c_{m'-2} = 0$ |

### Special values

Define $n!! = n(n-2)(n-4)\cdots$ (terminating at 1 or 2).

|  Eq.  |   |
| :---: | :-- |
| (A1.36) | $P_n^{\|m\|}(\cos\theta)\bigm\|_{\theta=0} = \begin{cases} 1, & \|m\| = 0 \\ 0, & \|m\| > 0\end{cases}$ |
| (A1.37) | $P_n^{\|m\|}(\cos\theta)\bigm\|_{\theta=\pi/2} = \begin{cases} (-1)^{(n-\|m\|)/2}\,\dfrac{(n+\|m\|-1)!!}{(n-\|m\|)!!}, & (n+\|m\|)\text{ even} \\ 0, & (n+\|m\|)\text{ odd}\end{cases}$ |
| (A1.38) | $P_n^{\|m\|}(\cos\theta)\bigm\|_{\theta=\pi} = \begin{cases} (-1)^n, & \|m\| = 0 \\ 0, & \|m\| > 0\end{cases}$ |
| (A1.39) | $\frac{m\,P_n^{\|m\|}(\cos\theta)}{\sin\theta}\biggm\|_{\theta=0} = \begin{cases} 0, & m \ne \pm 1 \\ \pm\,\dfrac{n(n+1)}{2}, & m = \pm 1\end{cases}$ |
| (A1.40) | $\frac{m\,P_n^{\|m\|}(\cos\theta)}{\sin\theta}\biggm\|_{\theta=\pi/2} = \begin{cases} m\,(-1)^{(n-\|m\|)/2}\,\dfrac{(n+\|m\|-1)!!}{(n-\|m\|)!!}, & (n+\|m\|)\text{ even} \\ 0, & (n+\|m\|)\text{ odd}\end{cases}$ |
| (A1.41) | $\frac{m\,P_n^{\|m\|}(\cos\theta)}{\sin\theta}\biggm\|_{\theta=\pi} = \begin{cases} 0, & m \ne \pm 1 \\ \pm(-1)^{n+1}\,\dfrac{n(n+1)}{2}, & m = \pm 1\end{cases}$ |
| (A1.42) | $\frac{dP_n^{\|m\|}(\cos\theta)}{d\theta}\biggm\|_{\theta=0} = \begin{cases} 0, & \|m\| \ne 1 \\ \dfrac{n(n+1)}{2}, & \|m\| = 1\end{cases}$ |
| (A1.43) | $\frac{dP_n^{\|m\|}(\cos\theta)}{d\theta}\biggm\|_{\theta=\pi/2} = \begin{cases} 0, & (n+\|m\|)\text{ even} \\ (-1)^{(n-\|m\|+1)/2}\,\dfrac{(n+\|m\|)!!}{(n-\|m\|-1)!!}, & (n+\|m\|)\text{ odd}\end{cases}$ |
| (A1.44) | $\frac{dP_n^{\|m\|}(\cos\theta)}{d\theta}\biggm\|_{\theta=\pi} = \begin{cases} 0, & \|m\| \ne 1 \\ (-1)^n\,\dfrac{n(n+1)}{2}, & \|m\| = 1\end{cases}$ |

### Low-order $\bar{P}_n^{|m|}(\cos\theta)$ table

|  | $n=1$ | $n=2$ | $n=3$ | $n=4$ | $n=5$ |
|---|---|---|---|---|---|
| $|m|=0$ | $\tfrac{\sqrt{6}}{2}\cos\theta$ | $\tfrac{\sqrt{10}}{8}(3\cos 2\theta + 1)$ | $\tfrac{\sqrt{14}}{16}(5\cos 3\theta + 3\cos\theta)$ | $\tfrac{3\sqrt{2}}{128}(35\cos 4\theta + 20\cos 2\theta + 9)$ | $\tfrac{\sqrt{22}}{256}(63\cos 5\theta + 35\cos 3\theta + 30\cos\theta)$ |
| $|m|=1$ | $\tfrac{\sqrt{3}}{2}\sin\theta$ | $\tfrac{\sqrt{15}}{4}\sin 2\theta$ | $\tfrac{\sqrt{42}}{32}(5\sin 3\theta + \sin\theta)$ | $\tfrac{3\sqrt{10}}{64}(7\sin 4\theta + 2\sin 2\theta)$ | $\tfrac{\sqrt{165}}{256}(21\sin 5\theta + 7\sin 3\theta + 2\sin\theta)$ |
| $|m|=2$ |  | $-\tfrac{\sqrt{15}}{8}(\cos 2\theta - 1)$ | $\tfrac{\sqrt{105}}{16}(\cos 3\theta - \cos\theta)$ | $-\tfrac{3\sqrt{5}}{64}(7\cos 4\theta - 4\cos 2\theta - 3)$ | $-\tfrac{\sqrt{1155}}{128}(3\cos 5\theta - \cos 3\theta - 2\cos\theta)$ |
| $|m|=3$ |  |  | $-\tfrac{\sqrt{70}}{32}(\sin 3\theta - 3\sin\theta)$ | $-\tfrac{3\sqrt{70}}{64}(\sin 4\theta - 2\sin 2\theta)$ | $-\tfrac{\sqrt{770}}{512}(9\sin 5\theta - 13\sin 3\theta - 6\sin\theta)$ |
| $|m|=4$ |  |  |  | $\tfrac{3\sqrt{35}}{128}(\cos 4\theta - 4\cos 2\theta + 3)$ | $\tfrac{3\sqrt{385}}{256}(\cos 5\theta - 3\cos 3\theta + 2\cos\theta)$ |
| $|m|=5$ |  |  |  |  | $\tfrac{3\sqrt{154}}{512}(\sin 5\theta - 5\sin 3\theta + 10\sin\theta)$ |

### Low-order $\dfrac{m}{\sin\theta}\bar{P}_n^{|m|}(\cos\theta)$ — used in the $\hat{\theta}$ component of $\vec{F}^{(c)}_{1mn}$

|  | $n=1$ | $n=2$ | $n=3$ | $n=4$ | $n=5$ |
|---|---|---|---|---|---|
| $m=0$ | 0 | 0 | 0 | 0 | 0 |
| $m=1$ | $\tfrac{\sqrt{3}}{2}$ | $\tfrac{\sqrt{15}}{2}\cos\theta$ | $\tfrac{\sqrt{42}}{16}(5\cos 2\theta + 3)$ | $\tfrac{3\sqrt{10}}{32}(7\cos 3\theta + 9\cos\theta)$ | $\tfrac{\sqrt{165}}{128}(21\cos 4\theta + 28\cos 2\theta + 15)$ |
| $m=2$ |  | $\tfrac{\sqrt{15}}{2}\sin\theta$ | $\tfrac{\sqrt{42}}{4}\sin 2\theta$ | $\tfrac{3\sqrt{5}}{16}(7\sin 3\theta + 3\sin\theta)$ | $\tfrac{\sqrt{1155}}{32}(3\sin 4\theta + 2\sin 2\theta)$ |
| $m=3$ |  |  | $-\tfrac{3\sqrt{70}}{16}(\cos 2\theta - 1)$ | $\tfrac{3\sqrt{5}}{16}(7\sin 3\theta + 3\sin\theta)$ | $-\tfrac{3\sqrt{770}}{256}(9\cos 4\theta - 4\cos 2\theta - 5)$ |
| $m=4$ |  |  |  | $-\tfrac{3\sqrt{35}}{16}(\sin 3\theta - 3\sin\theta)$ | $-\tfrac{3\sqrt{385}}{32}(\sin 4\theta - 2\sin 2\theta)$ |
| $m=5$ |  |  |  |  | $\tfrac{15\sqrt{154}}{256}(\cos 4\theta - 4\cos 2\theta + 3)$ |

### Low-order $\dfrac{d}{d\theta}\bar{P}_n^{|m|}(\cos\theta)$

|  | $n=1$ | $n=2$ | $n=3$ | $n=4$ | $n=5$ |
|---|---|---|---|---|---|
| $|m|=0$ | $-\tfrac{\sqrt{6}}{2}\sin\theta$ | $-\tfrac{3\sqrt{10}}{4}\sin 2\theta$ | $-\tfrac{3\sqrt{14}}{16}(5\sin 3\theta + \sin\theta)$ | $-\tfrac{15\sqrt{2}}{32}(7\sin 4\theta + 2\sin 2\theta)$ | $-\tfrac{15\sqrt{22}}{256}(21\sin 5\theta + 7\sin 3\theta + 2\sin\theta)$ |
| $|m|=1$ | $\tfrac{\sqrt{3}}{2}\cos\theta$ | $\tfrac{\sqrt{15}}{2}\cos 2\theta$ | $\tfrac{3\sqrt{42}}{32}(5\cos 3\theta + \cos\theta)$ | $\tfrac{3\sqrt{10}}{16}(7\cos 4\theta + \cos 2\theta)$ | $\tfrac{\sqrt{165}}{256}(105\cos 5\theta + 21\cos 3\theta + 2\cos\theta)$ |
| $|m|=2$ |  | $\tfrac{\sqrt{15}}{4}\sin 2\theta$ | $\tfrac{\sqrt{105}}{16}(3\sin 3\theta - \sin\theta)$ | $\tfrac{3\sqrt{5}}{16}(7\sin 4\theta - 2\sin 2\theta)$ | $\tfrac{\sqrt{1155}}{128}(15\sin 5\theta - 3\sin 3\theta - 2\sin\theta)$ |
| $|m|=3$ |  |  | $-\tfrac{3\sqrt{70}}{32}(\cos 3\theta - \cos\theta)$ | $\tfrac{3\sqrt{70}}{16}(7\sin 4\theta - 2\sin 2\theta)$ | $-\tfrac{3\sqrt{770}}{512}(15\cos 5\theta - 13\cos 3\theta - 2\cos\theta)$ |
| $|m|=4$ |  |  |  | $-\tfrac{3\sqrt{35}}{32}(\sin 4\theta - 2\sin 2\theta)$ | $-\tfrac{3\sqrt{385}}{256}(5\sin 5\theta - 9\sin 3\theta + 2\sin\theta)$ |
| $|m|=5$ |  |  |  |  | $\tfrac{15\sqrt{154}}{512}(\cos 5\theta - 3\cos 3\theta + 2\cos\theta)$ |

---

## 7. Spherical Wave Functions $\vec{F}^{(c)}_{smn}$ (§A1.3.1)

### General expressions

|  Eq.  |   |
| :---: | :-- |
| (A1.45) | $\boxed{\;\vec{F}^{(c)}_{1mn}(r,\theta,\phi) = \frac{1}{\sqrt{2\pi}}\,\frac{1}{\sqrt{n(n+1)}}\,\Bigl(-\frac{m}{\|m\|}\Bigr)^m \times \left[ z_n^{(c)}(kr)\,\frac{im\,\bar{P}_n^{\|m\|}(\cos\theta)}{\sin\theta}\,e^{im\phi}\,\hat{\theta} - z_n^{(c)}(kr)\,\frac{d\bar{P}_n^{\|m\|}(\cos\theta)}{d\theta}\,e^{im\phi}\,\hat{\phi}\right]\;}$ |
| (A1.46) | $\boxed{\;\vec{F}^{(c)}_{2mn}(r,\theta,\phi) = \frac{1}{\sqrt{2\pi}}\,\frac{1}{\sqrt{n(n+1)}}\,\Bigl(-\frac{m}{\|m\|}\Bigr)^m \times \left[\frac{n(n+1)}{kr}\,z_n^{(c)}(kr)\,\bar{P}_n^{\|m\|}(\cos\theta)\,e^{im\phi}\,\hat{r} + \frac{1}{kr}\frac{d}{d(kr)}\{kr\,z_n^{(c)}(kr)\}\,\frac{d\bar{P}_n^{\|m\|}(\cos\theta)}{d\theta}\,e^{im\phi}\,\hat{\theta} + \frac{1}{kr}\frac{d}{d(kr)}\{kr\,z_n^{(c)}(kr)\}\,\frac{im\,\bar{P}_n^{\|m\|}(\cos\theta)}{\sin\theta}\,e^{im\phi}\,\hat{\phi}\right]\;}$ |

Convention: $(-m/|m|)^m = 1$ when $m = 0$.

### Special values at $\theta = 0$ (A1.47, A1.48)

For $\vec{F}^{(c)}_{1mn}$:

|  Eq.  |   |
| :---: | :-- |
| (A1.47) | $\vec{F}^{(c)}_{1mn}(r,0,\phi) = \begin{cases} 0, & m > 1 \\ -\tfrac{1}{4}\sqrt{\tfrac{2n+1}{\pi}}\,z_n^{(c)}(kr)\,i\,e^{i\phi}(\hat{\theta} + i\hat{\phi}), & m = 1 \\ 0, & m = 0 \\ -\tfrac{1}{4}\sqrt{\tfrac{2n+1}{\pi}}\,z_n^{(c)}(kr)\,i\,e^{-i\phi}(\hat{\theta} - i\hat{\phi}), & m = -1 \\ 0, & m < -1 \end{cases}$ |

For $\vec{F}^{(c)}_{2mn}$:

|  Eq.  |   |
| :---: | :-- |
| (A1.48) | $\vec{F}^{(c)}_{2mn}(r,0,\phi) = \begin{cases} 0, & m > 1 \\ -\tfrac{1}{4}\sqrt{\tfrac{2n+1}{\pi}}\,\tfrac{1}{kr}\tfrac{d}{d(kr)}\{kr\,z_n^{(c)}\}\,e^{i\phi}(\hat{\theta} + i\hat{\phi}), & m = 1 \\ \sqrt{\tfrac{n(n+1)(2n+1)}{4\pi}}\,\dfrac{z_n^{(c)}(kr)}{kr}\,\hat{r}, & m = 0 \\ \tfrac{1}{4}\sqrt{\tfrac{2n+1}{\pi}}\,\tfrac{1}{kr}\tfrac{d}{d(kr)}\{kr\,z_n^{(c)}\}\,e^{-i\phi}(\hat{\theta} - i\hat{\phi}), & m = -1 \\ 0, & m < -1 \end{cases}$ |

### Special values at $\theta = \pi$ (A1.49, A1.50)

For $\vec{F}^{(c)}_{1mn}$ (same structure as A1.47 but with $(-1)^n$ prefactor and sign flip in second component):

|  Eq.  |   |
| :---: | :-- |
| (A1.49) | $\vec{F}^{(c)}_{1mn}(r,\pi,\phi) = \begin{cases} 0, & m > 1 \\ (-1)^n\tfrac{1}{4}\sqrt{\tfrac{2n+1}{\pi}}\,z_n^{(c)}(kr)\,i\,e^{i\phi}(\hat{\theta} - i\hat{\phi}), & m = 1 \\ 0, & m = 0 \\ (-1)^n\tfrac{1}{4}\sqrt{\tfrac{2n+1}{\pi}}\,z_n^{(c)}(kr)\,i\,e^{-i\phi}(\hat{\theta} + i\hat{\phi}), & m = -1 \\ 0, & m < -1 \end{cases}$ |

For $\vec{F}^{(c)}_{2mn}$:

|  Eq.  |   |
| :---: | :-- |
| (A1.50) | $\vec{F}^{(c)}_{2mn}(r,\pi,\phi) = \begin{cases} 0, & m > 1 \\ -(-1)^n\tfrac{1}{4}\sqrt{\tfrac{2n+1}{\pi}}\,\tfrac{1}{kr}\tfrac{d}{d(kr)}\{kr\,z_n^{(c)}\}\,e^{i\phi}(\hat{\theta} - i\hat{\phi}), & m = 1 \\ (-1)^n\sqrt{\tfrac{n(n+1)(2n+1)}{4\pi}}\,\dfrac{z_n^{(c)}(kr)}{kr}\,\hat{r}, & m = 0 \\ (-1)^n\tfrac{1}{4}\sqrt{\tfrac{2n+1}{\pi}}\,\tfrac{1}{kr}\tfrac{d}{d(kr)}\{kr\,z_n^{(c)}\}\,e^{-i\phi}(\hat{\theta} + i\hat{\phi}), & m = -1 \\ 0, & m < -1 \end{cases}$ |

### Value at the origin $r = 0$ — only $n = 1$, $s = 2$ survives

|  Eq.  |   |
| :---: | :-- |
| (A1.51) | $\vec{F}^{(1)}_{1,1,1}(0,\theta,\phi) = 0,\quad \vec{F}^{(1)}_{1,0,1}(0,\theta,\phi) = 0,\quad \vec{F}^{(1)}_{1,-1,1}(0,\theta,\phi) = 0$ |
| (A1.52a) | $\vec{F}^{(1)}_{2,1,1}(0,\theta,\phi) = -\tfrac{\sqrt{3}}{6\sqrt{\pi}}\,e^{i\phi}(\sin\theta\,\hat{r} + \cos\theta\,\hat{\theta} + i\hat{\phi}) = -\tfrac{\sqrt{3}}{6\sqrt{\pi}}(\hat{x} + i\hat{y})$ |
| (A1.52b) | $\vec{F}^{(1)}_{2,0,1}(0,\theta,\phi) = \tfrac{\sqrt{6}}{6\sqrt{\pi}}(\cos\theta\,\hat{r} - \sin\theta\,\hat{\theta}) = \tfrac{\sqrt{6}}{6\sqrt{\pi}}\,\hat{z}$ |
| (A1.52c) | $\vec{F}^{(1)}_{2,-1,1}(0,\theta,\phi) = \tfrac{\sqrt{3}}{6\sqrt{\pi}}\,e^{-i\phi}(\sin\theta\,\hat{r} + \cos\theta\,\hat{\theta} - i\hat{\phi}) = \tfrac{\sqrt{3}}{6\sqrt{\pi}}(\hat{x} - i\hat{y})$ |
| (A1.53) | $\vec{F}^{(1)}_{smn}(0,\theta,\phi) = 0,\quad n \ne 1$ |

### Conjugation identity

|  Eq.  |   |
| :---: | :-- |
| (A1.54) | $\boxed{\;\vec{F}^{(3)*}_{smn}(r,\theta,\phi) = (-1)^m\,\vec{F}^{(4)}_{s,-m,n}(r,\theta,\phi)\;}$ |

### Explicit $n = 1$ wave functions

**$\vec{F}^{(1)}_{1,m,1}$** (with $j_1$ form $-\cos(kr) + \sin(kr)/(kr)$):

|  Eq.  |   |
| :---: | :-- |
| (A1.55a) | $\vec{F}^{(1)}_{1,1,1} = \frac{-\sqrt{3}}{4\sqrt{\pi}}\,e^{i\phi}\,\frac{1}{kr}\Bigl\{-\cos(kr) + \frac{\sin(kr)}{kr}\Bigr\}(i\hat{\theta} - \cos\theta\,\hat{\phi})$ |
| (A1.55b) | $\vec{F}^{(1)}_{1,0,1} = \frac{\sqrt{6}}{4\sqrt{\pi}}\,\frac{1}{kr}\Bigl\{-\cos(kr) + \frac{\sin(kr)}{kr}\Bigr\}\sin\theta\,\hat{\phi}$ |
| (A1.55c) | $\vec{F}^{(1)}_{1,-1,1} = \frac{-\sqrt{3}}{4\sqrt{\pi}}\,e^{-i\phi}\,\frac{1}{kr}\Bigl\{-\cos(kr) + \frac{\sin(kr)}{kr}\Bigr\}(i\hat{\theta} + \cos\theta\,\hat{\phi})$ |

**$\vec{F}^{(1)}_{2,m,1}$:**

|  Eq.  |   |
| :---: | :-- |
| (A1.56a) | $\vec{F}^{(1)}_{2,1,1} = \frac{-\sqrt{3}}{2\sqrt{\pi}}\,e^{i\phi}\,\frac{1}{(kr)^2}\Bigl\{-\cos(kr) + \frac{\sin(kr)}{kr}\Bigr\}\sin\theta\,\hat{r} - \frac{\sqrt{3}}{4\sqrt{\pi}}\,e^{i\phi}\,\frac{1}{kr}\Bigl\{\sin(kr) + \frac{\cos(kr)}{kr} - \frac{\sin(kr)}{(kr)^2}\Bigr\}(\cos\theta\,\hat{\theta} + i\hat{\phi})$ |
| (A1.56b) | $\vec{F}^{(1)}_{2,0,1} = \frac{\sqrt{6}}{2\sqrt{\pi}}\,\frac{1}{(kr)^2}\Bigl\{-\cos(kr) + \frac{\sin(kr)}{kr}\Bigr\}\cos\theta\,\hat{r} - \frac{\sqrt{6}}{4\sqrt{\pi}}\,\frac{1}{kr}\Bigl\{\sin(kr) + \frac{\cos(kr)}{kr} - \frac{\sin(kr)}{(kr)^2}\Bigr\}\sin\theta\,\hat{\theta}$ |
| (A1.56c) | $\vec{F}^{(1)}_{2,-1,1} = \frac{\sqrt{3}}{2\sqrt{\pi}}\,e^{-i\phi}\,\frac{1}{(kr)^2}\Bigl\{-\cos(kr) + \frac{\sin(kr)}{kr}\Bigr\}\sin\theta\,\hat{r} + \frac{\sqrt{3}}{4\sqrt{\pi}}\,e^{-i\phi}\,\frac{1}{kr}\Bigl\{\sin(kr) + \frac{\cos(kr)}{kr} - \frac{\sin(kr)}{(kr)^2}\Bigr\}(\cos\theta\,\hat{\theta} - i\hat{\phi})$ |

**$\vec{F}^{(3)}_{1,m,1}$** (outgoing Hankel):

|  Eq.  |   |
| :---: | :-- |
| (A1.57a) | $\vec{F}^{(3)}_{1,1,1} = \frac{\sqrt{3}}{4\sqrt{\pi}}\,e^{i\phi}\,\frac{e^{ikr}}{kr}\Bigl(1 + \frac{i}{kr}\Bigr)(i\hat{\theta} - \cos\theta\,\hat{\phi})$ |
| (A1.57b) | $\vec{F}^{(3)}_{1,0,1} = -\frac{\sqrt{6}}{4\sqrt{\pi}}\,\frac{e^{ikr}}{kr}\Bigl(1 + \frac{i}{kr}\Bigr)\sin\theta\,\hat{\phi}$ |
| (A1.57c) | $\vec{F}^{(3)}_{1,-1,1} = \frac{\sqrt{3}}{4\sqrt{\pi}}\,e^{-i\phi}\,\frac{e^{ikr}}{kr}\Bigl(1 + \frac{i}{kr}\Bigr)(i\hat{\theta} + \cos\theta\,\hat{\phi})$ |

**$\vec{F}^{(3)}_{2,m,1}$** (outgoing Hankel):

|  Eq.  |   |
| :---: | :-- |
| (A1.58a) | $\vec{F}^{(3)}_{2,1,1} = \frac{\sqrt{3}}{2\sqrt{\pi}}\,e^{i\phi}\,\frac{e^{ikr}}{(kr)^2}\Bigl(1 + \frac{i}{kr}\Bigr)\sin\theta\,\hat{r} - \frac{\sqrt{3}}{4\sqrt{\pi}}\,e^{i\phi}\,\frac{e^{ikr}}{kr}\Bigl(-i + \frac{1}{kr} + \frac{i}{(kr)^2}\Bigr)(\cos\theta\,\hat{\theta} + i\hat{\phi})$ |
| (A1.58b) | $\vec{F}^{(3)}_{2,0,1} = -\frac{\sqrt{6}}{2\sqrt{\pi}}\,\frac{e^{ikr}}{(kr)^2}\Bigl(1 + \frac{i}{kr}\Bigr)\cos\theta\,\hat{r} - \frac{\sqrt{6}}{4\sqrt{\pi}}\,\frac{e^{ikr}}{kr}\Bigl(-i + \frac{1}{kr} + \frac{i}{(kr)^2}\Bigr)\sin\theta\,\hat{\theta}$ |
| (A1.58c) | $\vec{F}^{(3)}_{2,-1,1} = -\frac{\sqrt{3}}{2\sqrt{\pi}}\,e^{-i\phi}\,\frac{e^{ikr}}{(kr)^2}\Bigl(1 + \frac{i}{kr}\Bigr)\sin\theta\,\hat{r} + \frac{\sqrt{3}}{4\sqrt{\pi}}\,e^{-i\phi}\,\frac{e^{ikr}}{kr}\Bigl(-i + \frac{1}{kr} + \frac{i}{(kr)^2}\Bigr)(\cos\theta\,\hat{\theta} - i\hat{\phi})$ |

---

## 8. Far-Field Pattern Functions $\vec{K}_{smn}$ (§A1.3.2)

### General expressions

|  Eq.  |   |
| :---: | :-- |
| (A1.59) | $\boxed{\;\vec{K}_{1mn}(\theta,\phi) = \sqrt{\frac{2}{n(n+1)}}\Bigl(-\frac{m}{\|m\|}\Bigr)^m e^{im\phi}(-i)^{n+1}\left[\frac{im\bar{P}_n^{\|m\|}(\cos\theta)}{\sin\theta}\hat{\theta} - \frac{d\bar{P}_n^{\|m\|}(\cos\theta)}{d\theta}\hat{\phi}\right]\;}$ |
| (A1.60) | $\boxed{\;\vec{K}_{2mn}(\theta,\phi) = \sqrt{\frac{2}{n(n+1)}}\Bigl(-\frac{m}{\|m\|}\Bigr)^m e^{im\phi}(-i)^{n}\left[\frac{d\bar{P}_n^{\|m\|}(\cos\theta)}{d\theta}\hat{\theta} + \frac{im\bar{P}_n^{\|m\|}(\cos\theta)}{\sin\theta}\hat{\phi}\right]\;}$ |

### Special values

**At $\theta = 0$:**

|  Eq.  |   |
| :---: | :-- |
| (A1.61, A1.62) | $\vec{K}_{smn}(0,\phi) = \begin{cases} 0, & \|m\| > 1 \\ -(-i)^n\tfrac{1}{2}\sqrt{2n+1}\,e^{i\phi}(\hat{\theta} + i\hat{\phi}), & s=1,\,m=+1 \\ 0, & s=1,\,m=0 \\ -(-i)^n\tfrac{1}{2}\sqrt{2n+1}\,e^{-i\phi}(\hat{\theta} - i\hat{\phi}), & s=1,\,m=-1 \\ -(-i)^n\tfrac{1}{2}\sqrt{2n+1}\,e^{i\phi}(\hat{\theta} + i\hat{\phi}), & s=2,\,m=+1 \\ 0, & s=2,\,m=0 \\ (-i)^n\tfrac{1}{2}\sqrt{2n+1}\,e^{-i\phi}(\hat{\theta} - i\hat{\phi}), & s=2,\,m=-1 \end{cases}$ |

**At $\theta = \pi$:**

|  Eq.  |   |
| :---: | :-- |
| (A1.63, A1.64) | $\vec{K}_{smn}(\pi,\phi) = \begin{cases} 0, & \|m\| > 1 \\ i^n\tfrac{1}{2}\sqrt{2n+1}\,e^{i\phi}(\hat{\theta} - i\hat{\phi}), & s=1,\,m=+1 \\ 0, & s=1,\,m=0 \\ i^n\tfrac{1}{2}\sqrt{2n+1}\,e^{-i\phi}(\hat{\theta} + i\hat{\phi}), & s=1,\,m=-1 \\ -i^n\tfrac{1}{2}\sqrt{2n+1}\,e^{i\phi}(\hat{\theta} - i\hat{\phi}), & s=2,\,m=+1 \\ 0, & s=2,\,m=0 \\ i^n\tfrac{1}{2}\sqrt{2n+1}\,e^{-i\phi}(\hat{\theta} + i\hat{\phi}), & s=2,\,m=-1 \end{cases}$ |

### Squared magnitude on z-axis

|  Eq.  |   |
| :---: | :-- |
| (A1.65) | $\|\vec{K}_{smn}(0,\phi)\|^2 = \begin{cases} \dfrac{2n+1}{2}, & m = \pm 1 \\ 0, & m \ne \pm 1\end{cases}$ |

### Explicit $n = 1$ patterns

|  Eq.  |   |
| :---: | :-- |
| (A1.66a) | $\vec{K}_{1,1,1} = \tfrac{\sqrt{3}}{2}\,e^{i\phi}(i\hat{\theta} - \cos\theta\,\hat{\phi})$ |
| (A1.66b) | $\vec{K}_{1,0,1} = -\tfrac{\sqrt{6}}{2}\sin\theta\,\hat{\phi}$ |
| (A1.66c) | $\vec{K}_{1,-1,1} = \tfrac{\sqrt{3}}{2}\,e^{-i\phi}(i\hat{\theta} + \cos\theta\,\hat{\phi})$ |
| (A1.67a) | $\vec{K}_{2,1,1} = i\tfrac{\sqrt{3}}{2}\,e^{i\phi}(\cos\theta\,\hat{\theta} + i\hat{\phi})$ |
| (A1.67b) | $\vec{K}_{2,0,1} = i\tfrac{\sqrt{6}}{2}\sin\theta\,\hat{\theta}$ |
| (A1.67c) | $\vec{K}_{2,-1,1} = i\tfrac{\sqrt{3}}{2}\,e^{-i\phi}(-\cos\theta\,\hat{\theta} + i\hat{\phi})$ |

---

## 9. Orthogonality Integrals (§A1.4)

Notation: $\delta_{ij} = 1$ if $i=j$, else $0$. Integration is over the unit sphere ($\theta \in [0,\pi]$, $\phi \in [0, 2\pi]$).

### Radial-component product (A1.68)

|  Eq.  |   |
| :---: | :-- |
|  | $\int_0^{2\pi}\!\!\int_0^{\pi}\{\vec{F}^{(c)}_{smn}\cdot\hat{r}\}\{\vec{F}^{(\gamma)}_{\sigma\mu\nu}\cdot\hat{r}\}\sin\theta\,d\theta\,d\phi = \delta_{s\sigma}\,\delta_{s2}\,\delta_{m,-\mu}\,\delta_{n\nu}\,(-1)^m\,n(n+1)\,\frac{z_n^{(c)}(kr)}{kr}\,\frac{z_n^{(\gamma)}(kr)}{kr}$ |

### Tangential-component products (A1.69)

|  Eq.  |   |
| :---: | :-- |
|  | $\int_0^{2\pi}\!\!\int_0^{\pi}\bigl[\{\vec{F}^{(c)}_{smn}\cdot\hat{\theta}\}\{\vec{F}^{(\gamma)}_{\sigma\mu\nu}\cdot\hat{\theta}\} + \{\vec{F}^{(c)}_{smn}\cdot\hat{\phi}\}\{\vec{F}^{(\gamma)}_{\sigma\mu\nu}\cdot\hat{\phi}\}\bigr]\sin\theta\,d\theta\,d\phi = \delta_{s\sigma}\,\delta_{m,-\mu}\,\delta_{n\nu}\,(-1)^m\,R^{(c)}_{sn}(kr)\,R^{(\gamma)}_{sn}(kr)$ |

### Scalar product (A1.70)

|  Eq.  |   |
| :---: | :-- |
|  | $\int_0^{2\pi}\!\!\int_0^{\pi}\vec{F}^{(c)}_{smn}\cdot\vec{F}^{(\gamma)}_{\sigma\mu\nu}\,\sin\theta\,d\theta\,d\phi = \delta_{s\sigma}\,\delta_{m,-\mu}\,\delta_{n\nu}\,(-1)^m\left\{R^{(c)}_{sn}(kr)R^{(\gamma)}_{sn}(kr) + \delta_{s2}\,n(n+1)\,\frac{z_n^{(c)}(kr)}{kr}\,\frac{z_n^{(\gamma)}(kr)}{kr}\right\}$ |

### Vector (cross) product — used for power flux (A1.71)

|  Eq.  |   |
| :---: | :-- |
|  | $\int_0^{2\pi}\!\!\int_0^{\pi}\{\vec{F}^{(c)}_{smn}\times\vec{F}^{(\gamma)}_{\sigma\mu\nu}\}\cdot\hat{r}\,\sin\theta\,d\theta\,d\phi = -\delta_{s,3-\sigma}\,\delta_{m,-\mu}\,\delta_{n\nu}\,(-1)^{m+s}\,R^{(c)}_{sn}(kr)\,R^{(\gamma)}_{3-s,n}(kr)$ |

> **Sign discrepancy note.** Chapter 2 (Eq. 2.46) wrote this with $(-1)^{3-s}(-1)^m$ — equivalent to $(-1)^{m+s+1}$, which differs from the $(-1)^{m+s}$ here by an overall sign. Both forms appear in the book; resolve by checking against power-flow (must give $+\tfrac{1}{2}|Q^{(3)}|^2$ for outgoing modes). **Recommended unit test:** integrate the cross-product for a single outgoing mode and require $P > 0$.

### Reciprocity integral (A1.74)

For two single-mode fields each carrying $1\text{ W}^{1/2}$:

|  Eq.  |   |
| :---: | :-- |
| (A1.72) | $(\vec{E}^{(c)}_{smn},\vec{H}^{(c)}_{smn}) = \Bigl(\tfrac{k}{\sqrt{\eta}}\,\vec{F}^{(c)}_{smn},\,-ik\sqrt{\eta}\,\vec{F}^{(c)}_{3-s,m,n}\Bigr)$ |
| (A1.74) | $\int_S\{\vec{E}^{(c)}_{smn}\times\vec{H}^{(\gamma)}_{\sigma\mu\nu} - \vec{E}^{(\gamma)}_{\sigma\mu\nu}\times\vec{H}^{(c)}_{smn}\}\cdot d\vec{S} = \delta_{s\sigma}\,\delta_{m,-\mu}\,\delta_{n\nu}\,(-1)^m\,(-i)\,A^{(c,\gamma)}$ |

$A^{(c,\gamma)}$ from the table in §5. The surface $S$ can be deformed across any source-free region.

---

## 10. Spherical-Wave Coefficients from Current Distributions (§A1.5)

### General reciprocity theorem (A1.75)

For source pairs $(\vec{J}_1,\vec{M}_1)$ and $(\vec{J}_2,\vec{M}_2)$ producing $(\vec{E}_1,\vec{H}_1)$ and $(\vec{E}_2,\vec{H}_2)$:

|  Eq.  |   |
| :---: | :-- |
|  | $\int_S(\vec{E}_1\times\vec{H}_2 - \vec{E}_2\times\vec{H}_1)\cdot d\vec{S} = \int_V(\vec{E}_2\cdot\vec{J}_1 - \vec{E}_1\cdot\vec{J}_2 - \vec{H}_2\cdot\vec{M}_1 + \vec{H}_1\cdot\vec{M}_2)\,dV$ |

### Exterior field (sources inside $S$, field outside)

|  Eq.  |   |
| :---: | :-- |
| (A1.76) | $\vec{E}(r,\theta,\phi) = \frac{k}{\sqrt{\eta}}\sum_{smn} Q^{(3)}_{smn}\,\vec{F}^{(3)}_{smn}(r,\theta,\phi)$ |
| (A1.78) | $\boxed{\;Q^{(3)}_{smn} = (-1)^{m+1}\int_V\left(\frac{k}{\sqrt{\eta}}\,\vec{F}^{(1)}_{s,-m,n}\cdot\vec{J} + ik\sqrt{\eta}\,\vec{F}^{(1)}_{3-s,-m,n}\cdot\vec{M}\right)dV\;}$ |

### Interior field (sources outside, field inside)

|  Eq.  |   |
| :---: | :-- |
| (A1.79) | $\vec{E}(r,\theta,\phi) = \frac{k}{\sqrt{\eta}}\sum_{smn} Q^{(1)}_{smn}\,\vec{F}^{(1)}_{smn}(r,\theta,\phi)$ |
| (A1.81) | $\boxed{\;Q^{(1)}_{smn} = (-1)^{m+1}\int_V\left(\frac{k}{\sqrt{\eta}}\,\vec{F}^{(3)}_{s,-m,n}\cdot\vec{J} + ik\sqrt{\eta}\,\vec{F}^{(3)}_{3-s,-m,n}\cdot\vec{M}\right)dV\;}$ |

### Source symmetries (§A1.5.2) — electric currents only

**Rotational symmetry of order $m$:**

|  Eq.  |   |
| :---: | :-- |
| (A1.83) | $\vec{J}(r,\theta,\phi) = e^{im\phi}\{J_r(r,\theta)\hat{r} + J_\theta(r,\theta)\hat{\theta} + J_\phi(r,\theta)\hat{\phi}\} \;\Longrightarrow\; Q^{(c)}_{sm'n} = 0 \text{ for } m' \ne m$ |

**Planar image** in plane $\phi = w$, $\phi = w + \pi$:

|  Eq.  |   |
| :---: | :-- |
| (A1.87) | $\vec{J}^w(r,\theta,\phi) = J_r(r,\theta,2w-\phi)\hat{r} + J_\theta(r,\theta,2w-\phi)\hat{\theta} - J_\phi(r,\theta,2w-\phi)\hat{\phi}$ |
| (A1.88) | $Q^{(c)w}_{smn} = (-1)^{m+s}\,e^{-im2w}\,Q^{(c)}_{s,-m,n}$ |

If $\phi = u$ is a symmetry plane:

|  Eq.  |   |
| :---: | :-- |
| (A1.89) | $Q^{(c)}_{s,-m,n} = (-1)^{m+s}\,e^{im2u}\,Q^{(c)}_{smn}$ |

If $\phi = v$ is an anti-symmetry plane:

|  Eq.  |   |
| :---: | :-- |
| (A1.90) | $Q^{(c)}_{s,-m,n} = (-1)^{m+s+1}\,e^{im2v}\,Q^{(c)}_{smn}$ |

**Both symmetries simultaneously** (requires $v = u + \pi/2 + p\pi$ for integer $p$):

|  Eq.  |   |
| :---: | :-- |
| (A1.91a) | $Q^{(c)}_{smn} = 0 \text{ for } m \text{ even}$ |
| (A1.91b) | $Q^{(c)}_{s,-m,n} = (-1)^{s+1}\,e^{im2u}\,Q^{(c)}_{smn} \text{ for } m \text{ odd}$ |

This is the symmetry of a linearly polarized plane wave on the $z$-axis.

**Rotation by $\phi_o$ about $z$:** $\vec{J}_{\phi_o}(r,\theta,\phi) = \vec{J}(r,\theta,\phi - \phi_o)$ gives:

|  Eq.  |   |
| :---: | :-- |
| (A1.94) | $\boxed{\;Q^{(c)}_{smn,\phi_o} = e^{-im\phi_o}\,Q^{(c)}_{smn}\;}$ |

### Continuous x-polarized planar current ring (§A1.5.3)

Position: ring of radius $r_o$ in plane $z = r_o\cos\theta_o$, $\theta = \theta_o$.

|  Eq.  |   |
| :---: | :-- |
| (A1.95) | $\vec{J}(r,\theta,\phi) = J_o(r,\theta)\,\hat{x},\quad J_o(r,\theta) = \frac{\delta(r-r_o)\,\delta(\theta - \theta_o)}{2\pi r_o^2 \sin\theta_o}\,d_e$ |

Normalized so that $\int_V \vec{J}\,dV = d_e\,\hat{x}$. As $\theta_o \to 0$ or $\pi$, it becomes a short dipole on the $z$-axis.

**Symmetries:** rotational orders $m = \pm 1$ only; double planar symmetry with $u = 0$.

**Spherical mode coefficients** ($c = 1$ interior, $c = 3$ exterior):

|  Eq.  |   |
| :---: | :-- |
| (A1.97a) | $Q^{(c)}_{smn} = 0 \text{ for } m \ne \pm 1$ |
| (A1.97b) | $Q^{(c)}_{1mn} = \frac{-ik\,d_e}{\sqrt{8\pi n(n+1)\eta}}\,R^{(4-c)}_{1n}(kr_o)\left\{\cos\theta_o\,\frac{\bar{P}_n^1(\cos\theta_o)}{\sin\theta_o} + \frac{d\bar{P}_n^1(\cos\theta)}{d\theta}\biggm\|_{\theta=\theta_o}\right\},\quad m = \pm 1$ |
| (A1.97c) | $Q^{(c)}_{2mn} = \frac{m\,k\,d_e}{\sqrt{8\pi n(n+1)\eta}}\Bigl[\frac{n(n+1)}{kr_o}R^{(4-c)}_{1n}(kr_o)\,\bar{P}_n^1(\cos\theta_o)\,\sin\theta_o + R^{(4-c)}_{2n}(kr_o)\Bigl\{\cos\theta_o\,\frac{d\bar{P}_n^1(\cos\theta)}{d\theta}\biggm\|_{\theta=\theta_o} + \frac{\bar{P}_n^1(\cos\theta_o)}{\sin\theta_o}\Bigr\}\Bigr],\quad m = \pm 1$ |

### Sampled x-polarized planar current ring (§A1.5.4)

$L$ equally-spaced samples (uniform ring array) at $\phi_l = \phi_o + (l-1)\,2\pi/L$, $l = 1, \dots, L$:

|  Eq.  |   |
| :---: | :-- |
| (A1.98) | $\vec{J}^L(r,\theta,\phi) = \frac{2\pi}{L}\sum_{l=1}^{L}\delta(\phi - \phi_l)\,\vec{J}(r,\theta,\phi)$ |

Each sample is a short dipole of moment $d_e/L$ — total moment $\int_V\vec{J}^L\,dV = d_e\,\hat{x}$.

**Symmetries:** for $\phi_o = p\pi/L$ ($p$ integer): symmetry about $\phi = 0$; if also $L$ is even, anti-symmetry about $\phi = \pi/2$.

**Mode coefficients for $\phi_o = 0$:**

|  Eq.  |   |
| :---: | :-- |
| (A1.100a) | $Q^{(c)}_{1mn} = \frac{ik\,d_e}{\sqrt{8\pi n(n+1)\eta}}\Bigl(-\frac{m}{\|m\|}\Bigr)^m\,R^{(4-c)}_{1n}(kr_o)\Bigl\{\frac{m\bar{P}_n^{\|m\|}(\cos\theta_o)}{\sin\theta_o}\cos\theta_o\,(\delta^{L}_{m1} + \delta^{L}_{m,-1}) + \frac{d\bar{P}_n^{\|m\|}(\cos\theta)}{d\theta}\biggm\|_{\theta=\theta_o}(\delta^{L}_{m1} - \delta^{L}_{m,-1})\Bigr\}$ |
| (A1.100b) | $Q^{(c)}_{2mn} = \frac{-k\,d_e}{\sqrt{8\pi n(n+1)\eta}}\Bigl(-\frac{m}{\|m\|}\Bigr)^m\Bigl[\frac{n(n+1)}{kr_o}R^{(4-c)}_{1n}(kr_o)\bar{P}_n^{\|m\|}(\cos\theta_o)\sin\theta_o\,(\delta^{L}_{m1} + \delta^{L}_{m,-1}) + R^{(4-c)}_{2n}(kr_o)\Bigl\{\frac{d\bar{P}_n^{\|m\|}(\cos\theta)}{d\theta}\biggm\|_{\theta=\theta_o}\cos\theta_o\,(\delta^{L}_{m1} + \delta^{L}_{m,-1}) + \frac{m\bar{P}_n^{\|m\|}(\cos\theta_o)}{\sin\theta_o}(\delta^{L}_{m1} - \delta^{L}_{m,-1})\Bigr\}\Bigr]$ |

with the **periodic Kronecker delta**:

|  Eq.  |   |
| :---: | :-- |
| (A1.101) | $\delta^{L}_{m\mu} = \begin{cases} 1, & m \equiv \mu \pmod L \\ 0, & \text{otherwise}\end{cases}$ |

> **Sampling sidebands:** the sampled ring excites $m = \pm 1 + jL$ for $j = \pm 1, \pm 2, \dots$. If $L \ge 3$, sidebands do not interfere with the $m = \pm 1$ "main band," and the main-band coefficients are identical to the continuous-ring values.

For $\phi_o \ne 0$, replace $\delta^{L}_{m\mu}$ with $\delta^{L}_{m\mu}\,e^{i(m-\mu)\phi_o}$.

### Ring of currents tangential to a sphere (§A1.5.5)

Approximate $x$-polarization for $\theta_o$ near $\pi$:

|  Eq.  |   |
| :---: | :-- |
| (A1.102) | $\vec{J}_t(r,\theta,\phi) = J_o(r,\theta)\,\hat{p}_t,\quad \hat{p}_t = -\cos\phi\,\hat{\theta} - \sin\phi\,\hat{\phi}$ |

(Approaches $\hat{x}$-polarized short dipole on the negative $z$-axis as $\theta \to \pi$.)

**Mode coefficients for $\phi_o = 0$, $L \ge 3$:**

|  Eq.  |   |
| :---: | :-- |
| (A1.103a) | $Q^{(c)}_{1mn} = \frac{ik\,d_e}{\sqrt{8\pi n(n+1)\eta}}\Bigl(-\frac{m}{\|m\|}\Bigr)^m R^{(4-c)}_{1n}(kr_o)\Bigl\{\frac{-m\bar{P}_n^{\|m\|}(\cos\theta_o)}{\sin\theta_o}(\delta^{L}_{m1} + \delta^{L}_{m,-1}) + \frac{d\bar{P}_n^{\|m\|}(\cos\theta)}{d\theta}\biggm\|_{\theta=\theta_o}(\delta^{L}_{m1} - \delta^{L}_{m,-1})\Bigr\}$ |
| (A1.103b) | $Q^{(c)}_{2mn} = \frac{-k\,d_e}{\sqrt{8\pi n(n+1)\eta}}\Bigl(-\frac{m}{\|m\|}\Bigr)^m R^{(4-c)}_{2n}(kr_o)\Bigl\{-\frac{d\bar{P}_n^{\|m\|}(\cos\theta)}{d\theta}\biggm\|_{\theta=\theta_o}(\delta^{L}_{m1} + \delta^{L}_{m,-1}) + \frac{m\bar{P}_n^{\|m\|}(\cos\theta_o)}{\sin\theta_o}(\delta^{L}_{m1} - \delta^{L}_{m,-1})\Bigr\}$ |

---

## 11. The Plane Wave (§A1.6)

A plane wave $\vec{E}_o\,e^{i\vec{k}_o\cdot\vec{r}}$ arriving from direction $(\theta_o,\phi_o)$:

|  Eq.  |   |
| :---: | :-- |
| (A1.104) | $\vec{k}_o = -k\sin\theta_o\cos\phi_o\,\hat{x} - k\sin\theta_o\sin\phi_o\,\hat{y} - k\cos\theta_o\,\hat{z}$ |

**Standing-wave expansion** (valid for $r < R \approx N/k$):

|  Eq.  |   |
| :---: | :-- |
| (A1.105) | $\vec{E}_o\,e^{i\vec{k}_o\cdot\vec{r}} = \frac{k}{\sqrt{\eta}}\sum_{s=1}^{2}\sum_{n=1}^{N}\sum_{m=-n}^{n}Q_{smn}\,\vec{F}^{(1)}_{smn}(r,\theta,\phi)$ |

**Coefficient formula:**

|  Eq.  |   |
| :---: | :-- |
| (A1.106) | $\boxed{\;Q_{smn} = \frac{\sqrt{\eta}}{k}\,(-1)^m\sqrt{4\pi}\,i\,\vec{E}_o\cdot\vec{K}_{s,-m,n}(\theta_o,\phi_o)\;}$ |

### Validation example: $\hat{x}$-polarized plane wave traveling along $+\hat{z}$

|  Eq.  |   |
| :---: | :-- |
| (A1.107) | $\hat{x}\,E_o\,e^{ikz} = \frac{k}{\sqrt{\eta}}\sum_{s,n,m}Q_{smn}\,\vec{F}^{(1)}_{smn}(r,\theta,\phi)$ |

All coefficients vanish except:

|  Eq.  |   |
| :---: | :-- |
| (A1.108) | $\boxed{\;Q_{1,1,n} = Q_{1,-1,n} = Q_{2,1,n} = -Q_{2,-1,n} = \frac{\sqrt{\eta}}{k}\,E_o\,\sqrt{4\pi}\,i^{n+1}\,\tfrac{1}{2}\sqrt{2n+1}\;}$ |

Use this as a **golden test** in code: build an $\hat{x}$-polarized plane wave by other means (Cartesian expansion of $e^{ikz}$) and verify the spherical-mode expansion matches.

---

## 12. Coordinate / Unit-Vector Relations (§A1.7)

### Coordinate transformation

| Spherical → Rectangular | Rectangular → Spherical |
|---|---|
| $x = r\sin\theta\cos\phi$ | $r = \sqrt{x^2 + y^2 + z^2}$ |
| $y = r\sin\theta\sin\phi$ | $\theta = \texttt{atan2}(\sqrt{x^2+y^2},\,z)$ |
| $z = r\cos\theta$ | $\phi = \texttt{atan2}(y,\,x)$ |

### Unit-vector scalar products

| | $\hat{x}$ | $\hat{y}$ | $\hat{z}$ |
|---|---|---|---|
| $\hat{r}$ | $\sin\theta\cos\phi$ | $\sin\theta\sin\phi$ | $\cos\theta$ |
| $\hat{\theta}$ | $\cos\theta\cos\phi$ | $\cos\theta\sin\phi$ | $-\sin\theta$ |
| $\hat{\phi}$ | $-\sin\phi$ | $\cos\phi$ | 0 |

### Behaviour of vector distributions at $\theta = 0$ and $\theta = \pi$

| Distribution on sphere | $\theta = 0$ | $\theta = \pi$ |
|---|---|---|
| $\hat{r}$ | $\hat{z}$ | $-\hat{z}$ |
| $\hat{\theta}$ | discontinuous | discontinuous |
| $\hat{\phi}$ | discontinuous | discontinuous |
| $\hat{\theta}\cos\phi - \hat{\phi}\sin\phi$ | $\hat{x}$ | discontinuous |
| $\hat{\theta}\sin\phi + \hat{\phi}\cos\phi$ | $\hat{y}$ | discontinuous |
| $\hat{\theta}\cos\phi + \hat{\phi}\sin\phi$ | discontinuous | $-\hat{x}$ |
| $\hat{\theta}\sin\phi - \hat{\phi}\cos\phi$ | discontinuous | $-\hat{y}$ |
| $e^{i\phi}(\hat{\theta} + i\hat{\phi})$ | $\hat{x} + i\hat{y}$ | discontinuous |
| $e^{-i\phi}(\hat{\theta} - i\hat{\phi})$ | $\hat{x} - i\hat{y}$ | discontinuous |
| $e^{i\phi}(\hat{\theta} - i\hat{\phi})$ | discontinuous | $-\hat{x} - i\hat{y}$ |
| $e^{-i\phi}(\hat{\theta} + i\hat{\phi})$ | discontinuous | $-\hat{x} + i\hat{y}$ |

---

## 13. Implementation Checklist for Python Code

Each item below should map to a unit test or assertion.

### Indexing
1. **Single-index forward/inverse.** Round-trip $(s,m,n) \to j \to (s,m,n)$ for all $j \in [1, 2N(N+2)]$ with $N = 10$. Validate against the table in §3.
2. **Total mode count.** $J = 2N(N+2)$.

### Radial functions
3. **Bessel/Neumann definitions.** Match `scipy.special.spherical_jn`, `scipy.special.spherical_yn`. Note scipy's $y_n = n_n$.
4. **Hankel via combination.** $h_n^{(1)} = j_n + i\,n_n$, $h_n^{(2)} = j_n - i\,n_n$.
5. **R-function distinction.** $R^{(c)}_{1n} = z_n^{(c)}$, $R^{(c)}_{2n} = (1/x)\,d/dx\{x z_n^{(c)}\}$. Use **stable derivative** form (e.g. via recurrence A1.8 or A1.9), not numerical differentiation.
6. **Recurrence identity check.** $z_n^{(c)}/x = (z_{n-1}^{(c)} + z_{n+1}^{(c)})/(2n+1)$ — Eq. (A1.7).
7. **Cross-product identity.** Eq. (A1.13): $R^{(1)}_{1n}R^{(1)}_{2n} + R^{(2)}_{1n}R^{(2)}_{2n}$ must match the RHS form involving $|h_n^{(1)}|^2$.
8. **Wronskian.** Eq. (A1.11) with the $A^{(c,\gamma)}$ table from §5; pick $c=1,\gamma=2,s=1$ and verify $R^{(1)}_{1n}R^{(2)}_{2n} - R^{(2)}_{1n}R^{(1)}_{2n} = (kr)^{-2}$ for varied $n, kr$.
9. **Asymptotic forms.** For $kr \gg n$, the four limits in Eqs. (A1.14)–(A1.17). Test at $kr = 100\,n$.
10. **Origin limits.** $\lim_{x\to 0} j_n(x)/x = 1/3$ for $n=1$, $0$ for $n > 1$ — Eq. (A1.19). Same with derivative — Eq. (A1.20).
11. **Low-order explicit forms.** Compare against the closed-form tables in §5 for $n = 0\dots5$ at, e.g., $x = 1, 5, 10$.

### Angular functions
12. **Legendre sign convention.** Hansen does **not** include the $(-1)^m$ Condon–Shortley phase. Use `scipy.special.lpmv` which matches.
13. **Normalization (A1.27).** $\int_{-1}^1 \bar{P}_n^m(\mu)^2 d\mu = 1$. Verify by Gauss–Legendre quadrature.
14. **Special values.** $P_n^{|m|}(\cos 0) = \delta_{m0}$ (A1.36); $P_n^{|m|}(\cos\pi) = (-1)^n \delta_{m0}$ (A1.38); on-axis $m/\sin\theta$ values from (A1.39); derivative on-axis from (A1.42).
15. **Low-order tables.** Compare $\bar{P}_n^{|m|}$ against the multiple-angle expressions in §6 for $n,|m| \le 5$.

### Wave functions
16. **$\vec{F}^{(c)}_{1mn}$ has no radial component.** $\vec{F}^{(c)}_{1mn}\cdot\hat{r} \equiv 0$.
17. **Phase factor.** $(-m/|m|)^m = 1$ when $m=0$ (A1.46 convention).
18. **Origin behaviour.** $\vec{F}^{(1)}_{smn}(0,\theta,\phi) = 0$ unless $s=2, n=1$ — Eq. (A1.53). The non-zero values are concrete vectors in Cartesian form: $\vec{F}^{(1)}_{201}(0) = (\sqrt{6}/(6\sqrt{\pi}))\hat{z}$, $\vec{F}^{(1)}_{2,\pm 1,1}(0) = \mp(\sqrt{3}/(6\sqrt{\pi}))(\hat{x}\pm i\hat{y})$ — Eq. (A1.52).
19. **On-axis values.** For $m \notin \{-1,0,+1\}$, $\vec{F}^{(c)}_{smn}(r,0,\phi) = 0$. For $m=0$, only the $\hat{r}$ component of $\vec{F}^{(c)}_{2,0,n}$ survives. See Eqs. (A1.47–A1.50).
20. **Conjugation identity.** $\vec{F}^{(3)*}_{smn} = (-1)^m\,\vec{F}^{(4)}_{s,-m,n}$ — Eq. (A1.54). Easy point-wise test.
21. **Explicit $n=1$ forms.** Compare against Eqs. (A1.55)–(A1.58) at several $(r,\theta,\phi)$ points.

### Far-field patterns
22. **Definition.** $\vec{K}_{smn}(\theta,\phi) = \lim_{kr\to\infty}\{\sqrt{4\pi}\,kr/e^{ikr}\,\vec{F}^{(3)}_{smn}\}$ — Eq. (A1.5). Test by evaluating $\vec{F}^{(3)}_{smn}$ at large $kr$ and stripping the asymptotic factor.
23. **Helicity identity.** $\vec{K}_{smn} = i\hat{r}\times\vec{K}_{3-s,m,n}$ (Eq. 2.178 in Chapter 2).
24. **Z-axis magnitude.** $|\vec{K}_{smn}(0,\phi)|^2 = (2n+1)/2$ when $m = \pm 1$, zero otherwise — Eq. (A1.65). One-line unit test for normalization.
25. **Explicit $n=1$ patterns.** Closed forms in Eqs. (A1.66), (A1.67).

### Orthogonality
26. **Tangential-product orthogonality.** Numerically integrate Eq. (A1.69) for two outgoing modes ($c = \gamma = 3$) on a unit sphere: result should be zero unless $s=\sigma,m=-\mu,n=\nu$.
27. **Vector-product (power) integral.** Eq. (A1.71). For $(s,m,n)$ matched with $(3-s,-m,n)$, the integral gives a non-zero radial-flux density. Test consistency with the power formula $P = \tfrac{1}{2}\sum|Q^{(3)}|^2$ for a single mode (must give exactly $\tfrac{1}{2}$ when $|Q|=1$).

### Plane wave
28. **x-polarized plane wave along +z.** Use Eq. (A1.108) for coefficients, then synthesize $\vec{E}$ via Eq. (A1.107) and compare to $\hat{x}\,E_o\,e^{ikz}$. Convergence should be very fast for points near the origin and slower as $r$ approaches $N/k$.
29. **General plane wave.** Eq. (A1.106). Verify by reconstructing $\vec{E}_o\,e^{i\vec{k}_o\cdot\vec{r}}$ from coefficients for a chosen $(\theta_o, \phi_o)$.

### Current sources
30. **Exterior coefficient.** Use Eq. (A1.78) to compute the $Q^{(3)}_{smn}$ of a known source (e.g. short $\hat{z}$ dipole), then verify the radiated field matches the analytic dipole field — cross-check with Eqs. (2.109), (2.110), (2.115)–(2.117) of Chapter 2.
31. **Interior coefficient.** Same with Eq. (A1.81); good test using a current ring whose interior field is known.
32. **Rotation rule.** Verify Eq. (A1.94): rotating a source by $\phi_o$ multiplies $Q^{(c)}_{smn}$ by $e^{-im\phi_o}$.
33. **Symmetric source.** For a source with mirror symmetry in the $x$-$z$ plane and anti-symmetry in $y$-$z$, verify that all even-$m$ coefficients vanish (Eq. A1.91a).
34. **Current ring (continuous).** Compute Eqs. (A1.97) and compare radiated field to direct dipole-sum integration.
35. **Sampling sidebands.** For the $L$-sample ring at $\phi_o = 0$, verify Eq. (A1.100): main band at $m = \pm 1$ identical to continuous ring (when $L \ge 3$), additional bands at $m = \pm 1 + jL$.

### Frame transformations
36. **Unit-vector matrix.** Round-trip a Cartesian vector $(\hat{x},\hat{y},\hat{z})$ through the spherical-basis transformation in §12 and back. Should be identity to machine precision.
37. **Pole behaviour.** Use the §12 table to handle the $\theta = 0, \pi$ singularities of $\hat{\theta}, \hat{\phi}$ when building Cartesian outputs from spherical-basis fields.

---

*End of Appendix A1 reference.*
