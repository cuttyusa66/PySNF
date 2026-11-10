# PySNF

**Py**thon **S**pherical **N**ear-**F**ield: a single-file library for
spherical near-field to far-field transformation of antenna measurement
data, implementing the SNIFT-style algorithm from Hansen [1].

Version 1.1.0 — adds latitude (sin θ) and Maximum Fourier Content (MFC)
sampling grids, and fixes a loss of accuracy above N ≈ 86 (see
[What's new in 1.1.0](#whats-new-in-110)).

## Installation

PySNF is a single module (`pysnf.py`).  Drop it into your project, or
add the containing directory to `PYTHONPATH`.

Dependencies:

- Python 3.9+
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
- `nearfield2farfield` — measured scan -> far-field pattern (one-shot).
- `field2wavecoeffs` — scan -> spherical-wave coefficients `T_{smn}`;
  `sampling='uniform' | 'latitude' | 'mfc'` selects the scan grid.
- `wavecoeffs2farfield_uniform` — `T_{smn}` -> far-field on a uniform
  (dtheta, dphi) grid (dipole output probe).
- `wavecoeffs2farfield` — `T_{smn}` -> far-field at arbitrary paired
  (theta, phi) directions (no probe normalization).

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

Cross-references in the source cite Hansen equation numbers (e.g.
`[1], (4.135)`).  See [HANSEN_INDEX.md](HANSEN_INDEX.md) for a
Hansen-equation → PySNF-function lookup table.

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

runs 21 checks and ends with `=== 21 of 21 checks passed ===` (exit
code 1 otherwise).  They cover the original round-trip and single-mode
tests, the normalized Legendre functions against Hansen's formulas, a
uniform round trip at N = 180, and the latitude / MFC grids — including
the worst-case antenna of [2] at N = 40 and N = 180.  It takes a few
seconds and needs about 1.6 GB of memory for the N = 180 cases.

## License

PySNF is released under the MIT License.  See [LICENSE](LICENSE) for
the full text.  `SPDX-License-Identifier: MIT` applies to all source
files in this repository.
