# PySNF

**Py**thon **S**pherical **N**ear-**F**ield: a single-file library for
spherical near-field to far-field transformation of antenna measurement
data, implementing the SNIFT-style algorithm from Hansen [1].

Version 1.2.0 — adds opt-in gain normalization of the probe pattern and
of the far-field output (see
[Gain-normalized (ideal) transform](#gain-normalized-ideal-transform)).
Version 1.1.0 added latitude (sin θ) and Maximum Fourier Content (MFC)
sampling grids, and fixed a loss of accuracy above N ≈ 86 (see
[What's new in 1.1.0](#whats-new-in-110)).

## Installation

Install from [PyPI](https://pypi.org/project/pysnf/):

```
pip install pysnf
```

This installs `pysnf` and the optional G/T helpers, `pysnf_gt`.  For
the plotting helper, use `pip install "pysnf[plot]"` to pull in
`matplotlib` as well.

PySNF is also still a single module (`pysnf.py`): you can instead drop
it into your project, or add the containing directory to `PYTHONPATH`.
`pysnf_gt.py` imports `pysnf`; copy it alongside if you need it.  The
example data file used below, `wr284_horn_Tsmn.csv`, is in the
[GitHub repository](https://github.com/cuttyusa66/PySNF) and is not
installed by `pip`.

Dependencies:

- Python 3.10+
- `numpy`
- `scipy` (uses `scipy.special.assoc_legendre_p_all`, requires
  SciPy >= 1.15)
- `matplotlib` (optional — only needed for
  `plot_spherical_wave_coefficients_mag_db`; lazy-imported)

## Quick start

Full probe-corrected round trip on a realistic antenna: load the
spherical-wave coefficients of a WR-284 standard-gain pyramidal horn
at 3.95 GHz (shipped as [wr284_horn_Tsmn.csv](wr284_horn_Tsmn.csv),
obtained via full-wave electromagnetic simulation), simulate the
probe signal at a finite measurement radius via the SNF transmission
formula, then recover the antenna coefficients with the
probe-corrected inverse transform:

```python
import numpy as np
import pysnf as sf

# ---- Antenna T_smn: load the WR-284 horn from CSV --------------------
# Columns in the CSV: n, m, s, Re_T, Im_T (triangle-valid modes only).
raw = np.loadtxt('wr284_horn_Tsmn.csv', delimiter=',', skiprows=5)
n_idx, m_idx, s_idx = raw[:, 0].astype(int), raw[:, 1].astype(int), raw[:, 2].astype(int)
n_max = int(n_idx.max())
m_max = n_max                                       # triangle storage: |m| <= n
T_ant = np.zeros((n_max, 2*m_max + 1, 2), dtype='complex128')
T_ant[n_idx - 1, m_idx + m_max, s_idx - 1] = raw[:, 3] + 1j*raw[:, 4]

# ---- Probe T_smn: short +x electric Hertzian dipole ------------------
T_probe = np.zeros((1, 3, 2), dtype='complex128')
T_probe[0, 1 + (-1), 1] = +np.sqrt(2.0)/2.0
T_probe[0, 1 + (+1), 1] = -np.sqrt(2.0)/2.0

# Convert the probe's transmitting T to receiving R_p (Hansen (5.67),
# (2.104)); transmission_formula wants the (sigma, mu, nu) layout.
R_probe = sf.reciprocity(sf.rotate_wavecoeffs_about_axis(T_probe, 'x'))
r_p = np.swapaxes(R_probe, 0, 2)

# ---- Forward SNF transmission at a realistic measurement radius -----
kA = 100.0                       # probe-to-origin distance * wavenumber
dT = dP = 2.0                    # scan grid, degrees
thetas = np.linspace(0.0, np.pi, int(round(180.0/dT)) + 1)
phis = np.linspace(0.0, 2.0*np.pi, int(round(360.0/dP)), endpoint=False)
T_grid, P_grid = np.meshgrid(thetas, phis, indexing='ij')

# Probe signal at chi = 0 (theta-pol arm) and chi = pi/2 (phi-pol arm).
E_theta = sf.transmission_formula(T_ant, 0.0,     T_grid, P_grid, r_p=r_p, ka=kA)
E_phi   = sf.transmission_formula(T_ant, np.pi/2, T_grid, P_grid, r_p=r_p, ka=kA)
scan = {'E_theta': E_theta, 'E_phi': E_phi, 'dT': dT, 'dP': dP}

# ---- Probe-corrected recovery of the antenna coefficients ------------
# field2wavecoeffs takes the probe as its own measured far-field pattern.
pE_theta, pE_phi = sf.wavecoeffs2farfield_uniform(T_probe, dT, dP)
probe_scan = {'E_theta': pE_theta, 'E_phi': pE_phi, 'dT': dT, 'dP': dP}

T_ant_recovered = sf.field2wavecoeffs(scan, probe=probe_scan, probe_pol='x',
                                      ka=kA, n_max=n_max, m_max=m_max)

print('relative max error =',
      np.max(np.abs(T_ant_recovered - T_ant)) / np.max(np.abs(T_ant)))
```

Expected output: `relative max error ~ 1e-15` — the round trip
recovers the horn's modes to numerical precision.

For the end-to-end **measurement** path — go from a measured scan on
a sphere of radius A directly to the test antenna's far field — use
the one-shot helper:

```python
farfield = sf.nearfield2farfield(nearfield, dtheta_out, dphi_out,
                                 probe=probe_scan, probe_pol='x', ka=k*A)
```

where `nearfield` is the measured probe signal and `probe_scan` is the
probe's own measured far-field pattern, both as dicts with keys
`'E_theta'`, `'E_phi'`, `'dT'`, `'dP'` (probe signal at chi = 0 and
chi = pi/2, plus theta/phi sampling intervals in degrees).

For the ideal dipole-probe limit (no custom probe, `ka = inf`), the
minimum synthesize + recover round trip in the self-test at the bottom
of `pysnf.py` remains the shortest example.

## Latitude and MFC sampling

A uniform scan takes the same number of phi samples on every theta
ring, which oversamples the rings near the poles.  A *latitude* grid
keeps the same rings but thins the phi samples ring by ring.  PySNF
supports two rules for the per-ring count:

- `'latitude'` — the classical sin θ weighting (Hansen [1],
  pp. 128-130): the equator count `2M + 1` scaled by sin θ.
- `'mfc'` — the Maximum Fourier Content rule of Cutshall & Lawrence [2],
  Eq. (15), which keeps aliasing below -100 dB for any antenna of the
  same `N`.

`sampling_grid` generates the sample locations, and `field2wavecoeffs`
(or `nearfield2farfield`) accepts the scan through its `sampling`
argument.  Continuing from the Quick start (`T_ant`, `n_max`, `m_max`):

```python
# ---- Where to sample ------------------------------------------------
grid = sf.sampling_grid(n_max, m_max, sampling='mfc')   # or 'latitude'
print(grid['n_points'], 'points;', grid['ring_counts'][:4], '... per ring')

# ---- Probe signal at those points -----------------------------------
# Here: simulated with an ideal far-field dipole probe.  In a measurement
# these are the chi = 0 and chi = pi/2 probe signals at each point.
k_theta, k_phi = sf.wavecoeffs2farfield(T_ant, np.deg2rad(grid['theta']),
                                        np.deg2rad(grid['phi']))
scan = {'theta': grid['theta'], 'phi': grid['phi'],     # degrees
        'E_theta': np.sqrt(6.0)/4.0 * k_theta,          # dipole-probe scaling
        'E_phi':   np.sqrt(6.0)/4.0 * k_phi}

# ---- Recover the coefficients ---------------------------------------
T_mfc = sf.field2wavecoeffs(scan, n_max=n_max, m_max=m_max, sampling='mfc')
print('relative max error =',
      np.max(np.abs(T_mfc - T_ant)) / np.max(np.abs(T_ant)))
```

Expected output: 2842 points (the uniform grid needs 3444) and a
relative max error of order 1e-15.

A latitude / MFC scan is four flat 1-D arrays of equal length, one
entry per sample point: `'theta'` and `'phi'` (degrees), `'E_theta'` and
`'E_phi'`.  The samples may be in any order; rings are found by grouping
on theta.  What is required (and checked, with an error naming the
offending ring):

- ring thetas equally spaced from 0 to 180 degrees, poles included;
- on each ring, phi equally spaced over 360 degrees starting at phi = 0;
- at least 3 samples per ring.

Each ring is Fourier-interpolated in phi onto a common uniform grid and
the standard transform is applied, so `'latitude'` and `'mfc'` scans are
processed identically; the two names differ only in `sampling_grid`.

Choosing between them: the sin θ grid has fewer points but is an inexact
scheme.  On the shipped horn it leaves about -72 dB of aliasing, and for
a directive antenna pointed at the horizon it fails outright (about
-25 dB at N = 40).  The MFC grid recovers both to better than -200 dB.
Both figures are reproduced by the self-test.

## Gain-normalized (ideal) transform

By default every far-field array PySNF reads or writes is the signal of
an infinitely remote Hertzian-dipole probe, `W = (√6/4)·K`, not the
absolute far-field pattern `K` itself (`|K|²` is gain for
power-normalized coefficients).  That is 4.26 dB below gain, and it
applies in two places:

- the **output** of `wavecoeffs2farfield_uniform` / `nearfield2farfield`;
- the **probe pattern** passed as `probe=` to `field2wavecoeffs` /
  `nearfield2farfield`.

Two opt-in arguments switch each to `K`:

| Argument | On | `'dipole'` (default) | `'gain'` |
|---|---|---|---|
| `normalization` | `wavecoeffs2farfield_uniform`, `nearfield2farfield` | output is `(√6/4)·K` | output is `K` |
| `probe_normalization` | `field2wavecoeffs`, `nearfield2farfield` | probe pattern is `(√6/4)·K_p` | probe pattern is `K_p` (gain-normalized) |

With both set to `'gain'` the transform adds no scaling of its own to
the measured data — an *ideal* SNF-to-FF transform in the sense of [3],
which is what absolute quantities such as G/T need.  The scan itself is
never rescaled.

The Appendix of [3] gives a closed-form check: a z-directed Hertzian
dipole measured with a Hertzian-dipole probe must come out at exactly
1.76 dBi, at any measurement distance.

```python
import numpy as np
import pysnf as sf

kA = 10.0
dT = dP = 5.0
thetas = np.linspace(0.0, np.pi, int(round(180.0/dT)) + 1)
phis = np.linspace(0.0, 2.0*np.pi, int(round(360.0/dP)), endpoint=False)
T, P = np.meshgrid(thetas, phis, indexing='ij')

# Probe signal, [3] Eq. (26): chi = 0 (theta-pol) and chi = pi/2 (phi-pol).
w = 0.75*np.exp(1j*kA)*(-1j/kA + 1/kA**2 + 1j/kA**3)*np.sin(T)
nearfield = {'E_theta': w, 'E_phi': np.zeros_like(w), 'dT': dT, 'dP': dP}

# Gain-normalized pattern of the x'-polarized dipole probe, [3] Eq. (27).
probe = {'E_theta': np.sqrt(6)/2*np.cos(T)*np.cos(P),
         'E_phi':  -np.sqrt(6)/2*np.sin(P), 'dT': dT, 'dP': dP}

ff = sf.nearfield2farfield(nearfield, dT, dP, probe=probe, probe_pol='x', ka=kA,
                           n_max=3, m_max=3, nu_max=1,
                           normalization='gain', probe_normalization='gain')
gain = np.abs(ff['E_theta'])**2 + np.abs(ff['E_phi'])**2
print('peak gain = {:.4f} dBi'.format(10*np.log10(gain.max())))
```

Expected output: `peak gain = 1.7609 dBi`.  With the defaults the same
call reports −6.76 dB: one factor of √6/4 from the probe pattern and one
from the output.

Two things to keep in mind:

- **Probe roll.** `'E_theta'` is the signal at χ = 0, where the probe's
  x′ axis lies along θ̂, and `'E_phi'` the signal at χ = π/2 (x′ along φ̂,
  y′ along −θ̂).  The probe pattern is always given in the probe's own
  frame.  For a y′-polarized probe ([3] Eq. (28)) the scan above
  therefore becomes `'E_theta': 0, 'E_phi': -w`; leaving `w` under
  `'E_theta'` gives a wrong answer without any error.
- **Phase.** [3] Eqs. (27)-(28) equal `i·K` of the Hertzian-dipole
  coefficients PySNF takes from Hansen [1].  A constant phase on the probe
  pattern appears, conjugated, on the recovered coefficients; magnitudes
  and gain are unaffected.

## G/T measurements (`pysnf_gt.py`)

The companion module [pysnf_gt.py](pysnf_gt.py) computes Gain over Noise
Temperature of an active unit under test by the CW-Ambient technique of
[3].  It needs four measurements: the power into the probe, the ambient
noise power out of the UUT (P1), the power out of the UUT with the probe
radiating a CW tone from one pose P (P2), and a full dual-polarized scan
that contains P.

```python
import numpy as np
import pysnf as sf
import pysnf_gt as gt

# --- A synthetic measurement: z-dipole UUT (gain 1.5), T_sys = 437 K ---
k_B, T_sys, B, P_in = gt.BOLTZMANN, 437.0, 2.0e6, 3.0e-4
kA = 10.0
dT = dP = 5.0
thetas = np.linspace(0.0, np.pi, int(round(180.0/dT)) + 1)
phis = np.linspace(0.0, 2.0*np.pi, int(round(360.0/dP)), endpoint=False)
T, P = np.meshgrid(thetas, phis, indexing='ij')
w = 0.75*np.exp(1j*kA)*(-1j/kA + 1/kA**2 + 1j/kA**3)*np.sin(T)   # [3] Eq. (26)
G_elex, a_sf = 1.0e3, 0.37*np.exp(0.7j)        # unknown to the method
ref = (18, 0)                                  # CW pose P: theta = 90 deg, phi = 0
nearfield = {'E_theta': np.sqrt(2*P_in*G_elex)*a_sf*w,           # receiver values
             'E_phi': np.zeros_like(w), 'dT': dT, 'dP': dP}
p_ambient = k_B*T_sys*B*G_elex                                   # P1
p_cw = P_in*abs(w[ref])**2*G_elex + p_ambient                    # P2

# --- Gain-normalized probe pattern, [3] Eq. (27) ---
probe = {'E_theta': np.sqrt(6)/2*np.cos(T)*np.cos(P),
         'E_phi':  -np.sqrt(6)/2*np.sin(P), 'dT': dT, 'dP': dP}

res = gt.g_over_t(nearfield, probe, kA, P_in, p_ambient, p_cw, B,
                  ref_index=ref, ref_pol='theta', n_max=3, m_max=3, nu_max=1)
print('peak G/T = {:.4f} dB/K'.format(10*np.log10(res['G_T'].max())))
```

Expected output: `peak G/T = -24.6439 dB/K`, i.e. `1.5 / 437`.

- **Probe gain.** The probe pattern must carry its gain
  (`|E_theta|² + |E_phi|²` = gain).  For a pattern normalized to 0 dB,
  pass `probe_gain=` (linear) and it is rescaled.  No separate division
  by the probe gain is needed, unlike [3] Eq. (19).
- **Receiver.** The scan is complex receiver values with any constant
  scale factor; it cancels.  Because the reference is taken from the scan
  itself, there is no RMS / peak-power choice to make ([3] Eqs. (19) and
  (23) coincide).  `p_in` must be true power in watts.
- **Output.** `res['G_T_theta']`, `res['G_T_phi']`, `res['G_T']` (total)
  in 1/K; pass `copol=(c_theta, c_phi)` for a specific polarization.
- `gt.near_zone_correction` returns the factor N_SNF of [3] Eq. (24).

`python pysnf_gt.py` runs the module's own four checks.  They are
synthetic: they verify the scale factors, not the assumptions of the
CW-Ambient technique (UUT noise dominant, matched lines, stable ambient),
which are listed in the module docstring.

## Conventions

A condensed summary; see the module docstring of `pysnf.py` for full
detail.

- **Time convention.** `exp(-i*omega*t)` (Hansen / physics).  Outgoing
  waves carry `h_n^{(1)}(kr) ~ e^{+ikr}/(kr)` at large `kr`.
- **Angles.** Public functions accept `dtheta`, `dphi` in **degrees**;
  internal spherical-wave-function evaluation uses **radians**.
- **DFT.** Hansen uses the positive-exponent / no-normalization
  forward transform — opposite of `numpy.fft`.  Use the internal
  `_fft` / `_ifft` wrappers when porting formulas straight from
  Hansen's text.
- **Probe.** Default is an infinitely remote +x-polarized electric
  Hertzian dipole.  Custom probes carry mu = +/-1 only (linearly
  polarized convention); mu = 0 (z-oriented dipole) is not yet
  supported.
- **Normalization.** Far-field arrays are dipole-probe signals,
  `(√6/4)·K`, unless `normalization='gain'` /
  `probe_normalization='gain'` is passed (see
  [Gain-normalized (ideal) transform](#gain-normalized-ideal-transform)).

## Array layout cheat sheet

PySNF stores spherical-wave-domain arrays in dense layouts indexed by
`(n, m, s)` triplets (or their `(n, mu, s)` / `(s, mu, nu)` variants):

| Variable in code      | Hansen symbol                        | Shape                            | Indexing                                   |
|-----------------------|--------------------------------------|----------------------------------|--------------------------------------------|
| `q_n_m_s` / `t_n_m_s` | `T_{smn}` (transmitting)             | `(n_max, 2*m_max+1, 2)`          | `[n-1, m + m_max, s-1]`                    |
| `p_n_mu_s`            | `P_{s mu n}(kA)` (probe response)    | `(n_max, 2, 2)`                  | `[n-1, mu_idx, s-1]`, mu in {-1, +1}       |
| `r_p`                 | `R^p_{sigma mu nu}` (probe receive)  | `(2, 2*mu_max+1, nu_max)`        | `[sigma-1, mu+mu_max, nu-1]`               |
| `c_s_n_sig_nu_mu`     | `C^{sn(c)}_{sigma mu nu}(kA)` (translation) | `(2, n_max, 2, nu_max, 2)` | `[s-1, n-1, sigma-1, nu-1, mu_idx]` |
| `deltas`              | `Delta^n_{m' m}`                     | `(n_max+1, 2*n_max+1, 2*n_max+1)`| `[n, m'+n_max, m+n_max]`                   |
| `d` (rotation coeff)  | `d^n_{mu m}(theta)`                  | `(n_max, 2*mu_max+1, 2*m_max+1, T)` | `[n-1, mu+mu_max, m+m_max, t]`         |
| `w_n_m_mu`            | `w_{nm}^{(mu)}` (4.92)               | `(n_max, 2*m_max+1, 2)`          | `[n-1, m + m_max, mu_idx]`                 |

`mu_idx`: `0` for `mu = -1`, `1` for `mu = +1`.  Spherical-wave coefficients are only
defined for `|m| <= n`; out-of-triangle entries of `q_n_m_s` must be zero
(see `_check_m_le_n`).

## Public API

**Transform pipeline**
- `nearfield2farfield` — measured scan -> far-field pattern (one-shot);
  takes both `normalization` and `probe_normalization`.
- `field2wavecoeffs` — scan -> spherical-wave coefficients `T_{smn}`;
  `sampling='uniform' | 'latitude' | 'mfc'` selects the scan grid,
  `probe_normalization='dipole' | 'gain'` the probe-pattern scaling.
- `wavecoeffs2farfield_uniform` — `T_{smn}` -> far-field on a uniform
  (dtheta, dphi) grid (dipole output probe);
  `normalization='dipole' | 'gain'` selects `(√6/4)·K` or `K`.
- `wavecoeffs2farfield` — `T_{smn}` -> far-field at arbitrary paired
  (theta, phi) directions (no probe normalization).

**G/T (`pysnf_gt.py`)**
- `g_over_t` — G/T of an active unit under test by the CW-Ambient
  technique of [3].
- `near_zone_correction` — the near-zone factor N_SNF of [3] Eq. (24).

**Sampling grids**
- `sampling_grid` — sample locations of an MFC, latitude (sin θ) or
  uniform grid.
- `mfc_phi_counts` — phi samples per ring, MFC rule ([2] Eq. (15)).
- `sin_theta_phi_counts` — phi samples per ring, sin θ rule.

**Building blocks**
- `transmission_formula` — Hansen Eq. (3.10) probe-signal evaluation.
- `farfield_pattern_functions` /
  `incomplete_farfield_pattern_functions` — `K_{smn}(theta, phi)`.
- `rotation_coefficients` — `d^n_{mu m}(theta)`.
- `translation_coefficients` — `C^{sn(3)}_{sigma mu nu}(kA)`.
- `delta_pyramid` — `Delta^n_{m' m}`.
- `lpmn_norm` — normalized associated Legendre functions
  (Hansen A1.25).
- `sph_hankel_first_kind` — `h_n^{(1)}(x)` for n = 0 .. N.

**Probes**
- `dipole_probe_response_constants` — ideal Hertzian dipole.
- `probe_response_constants` — generic probe from `R^p` coefficients.
- `reciprocity` — `T` <-> `R` via Hansen Eq. (2.104).
- `rotate_wavecoeffs_about_axis` — rotate `T_{smn}` (currently x-axis
  only).

**Utilities**
- `singlesphere2doublesphere` — symmetry-extend a single-sphere scan
  to a double sphere for FFT processing.
- `interpft` — FFT-based phi-axis resample.
- `db` — `20 log10 |a|`.
- `plot_spherical_wave_coefficients_mag_db` — quick-look plot of
  `|T_{smn}|` magnitudes.

## References

[1] J. E. Hansen (ed.), *Spherical Near-Field Antenna Measurements*,
IEE Electromagnetic Waves Series, vol. 26, Peter Peregrinus Ltd., 1988.

[2] R. T. Cutshall and J. D. Lawrence, "Experimental Measurements Using
the Uniform, Latitude, and Equally-Spaced Spherical Near-Field
Measurement Grids," *AMTA 2015 Proceedings*.

[3] R. T. Cutshall, J. Dobbins and M. N. Barr, "Measuring G/T with a
Spherical Near-Field Antenna Measurement System via the CW-Ambient
Technique," *2020 Antenna Measurement Techniques Association Symposium
(AMTA)*, Newport, RI, USA, 2020, pp. 1-6.

Cross-references in the source cite Hansen equation numbers (e.g.
`[1], (4.135)`).  See [HANSEN_INDEX.md](HANSEN_INDEX.md) for a
Hansen-equation → PySNF-function lookup table.

## What's new in 1.2.0

- **Gain normalization.** `wavecoeffs2farfield_uniform` and
  `nearfield2farfield` take `normalization='gain'` to return the absolute
  far-field pattern `K` instead of the dipole-probe signal `(√6/4)·K`.
  `field2wavecoeffs` and `nearfield2farfield` take
  `probe_normalization='gain'` to accept a gain-normalized probe pattern.
  Both default to `'dipole'`, so existing results are unchanged.
- **Self-test.** Six checks added, including the closed-form
  Hertzian-dipole case of [3].
- **G/T helpers.** New companion module `pysnf_gt.py` with `g_over_t`
  and `near_zone_correction`, implementing the CW-Ambient technique of
  [3] (see [G/T measurements](#gt-measurements-pysnf_gtpy)).

## What's new in 1.1.0

- **Latitude and MFC sampling.** `field2wavecoeffs` and
  `nearfield2farfield` take a `sampling` argument; `sampling_grid`,
  `mfc_phi_counts` and `sin_theta_phi_counts` are new.  The default
  (`'uniform'`) is unchanged.
- **Accuracy above N ≈ 86 fixed.** `lpmn_norm` scaled unnormalized
  Legendre functions by a factorial ratio that underflows double
  precision once `n + m` exceeds about 170.  Results were silently wrong
  from N = 87 and meaningless from N = 89.  It now uses SciPy's
  normalized recurrence; the uniform round trip is accurate to about
  2e-13 at N = 180.  Results below N = 86 are unchanged to rounding.
- **`lpmn_norm` at the poles.** Exactly at theta = 0 and theta = pi the
  function now returns Hansen's closed-form values.  Two outputs differ
  from 1.0.1: the `m = 0` value at theta = pi for odd `n` (the sign was
  wrong), and the `|m| = 1` theta-derivative at both poles (was NaN).
  Neither was used by the transforms, so no transform result changes.
- **Output grid size.** `wavecoeffs2farfield_uniform` sized its grid by
  truncation, so at spacings such as `180/169` degrees (N = 168) it
  returned one theta ring and one phi sample too few.  It now rounds.
- **Spurious warning removed.** `translation_coefficients` no longer
  raises "invalid value encountered in divide"; its output is unchanged.

## Self-test

```
python pysnf.py
```

runs 27 checks and ends with `=== 27 of 27 checks passed ===` (exit
code 1 otherwise).  They cover the original round-trip and single-mode
tests, the normalized Legendre functions against Hansen's formulas, a
uniform round trip at N = 180, the latitude / MFC grids — including
the worst-case antenna of [2] at N = 40 and N = 180 — and the
gain-normalized Hertzian-dipole case of [3].  It takes a few
seconds and needs about 1.6 GB of memory for the N = 180 cases.

The G/T module has its own, separate self-test: `python pysnf_gt.py`
runs 4 checks and ends with `=== 4 of 4 checks passed ===`.

## License

PySNF is released under the MIT License.  See [LICENSE](LICENSE) for
the full text.  `SPDX-License-Identifier: MIT` applies to all source
files in this repository.
