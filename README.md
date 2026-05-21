# PySNF

**Py**thon **S**pherical **N**ear-**F**ield: a single-file library for
spherical near-field to far-field transformation of antenna measurement
data, implementing the SNIFT-style algorithm from Hansen [1].

Version 0.1 — early release; the public API is stable for the headline
entry points but several helper functions are still lightly documented.

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

The minimum round-trip — synthesize a far-field from known spherical-wave
coefficients, then recover the coefficients back from that far-field
(this mirrors the self-test at the bottom of `pysnf.py`):

```python
import numpy as np
import pysnf as sf

# Build T_{smn} for a small antenna (x-directed electric dipole +
# z-directed magnetic dipole, n_max = 2, m_max = 1).
n_max, m_max = 2, 1
T_in = np.zeros((n_max, 2*m_max+1, 2), dtype='complex128')
T_in[0, m_max + (-1), 1] = +np.sqrt(2.0)/2.0   # T_{s=2, m=-1, n=1}
T_in[0, m_max + (+1), 1] = -np.sqrt(2.0)/2.0   # T_{s=2, m=+1, n=1}
T_in[0, m_max + 0,    0] = 1.0                 # T_{s=1, m= 0, n=1}

# Synthesize a far-field pattern on a 5-degree grid (dipole output probe).
dT = dP = 5.0   # degrees
E_theta, E_phi = sf.wavecoeffs2farfield_uniform(T_in, dT, dP)

# Recover T_{smn} from the synthesized pattern (dipole input probe, ka = inf).
ff = {'E_theta': E_theta, 'E_phi': E_phi, 'dT': dT, 'dP': dP}
T_out = sf.field2wavecoeffs(ff, probe=None, probe_pol='x', ka=np.inf,
                            n_max=n_max, m_max=m_max)

print('max |T_out - T_in| =', np.max(np.abs(T_out - T_in)))
```

For a measured (or simulated) scan with a non-ideal probe, swap the
two-call sequence above for the one-shot helper:

```python
farfield = sf.nearfield2farfield(nearfield, dtheta_out, dphi_out,
                                 probe=probe, probe_pol='x', ka=k*A)
```

where `nearfield` and `probe` are dicts with keys
`'E_theta'`, `'E_phi'`, `'dT'`, `'dP'` (probe signal at chi = 0 and
chi = pi/2, plus theta/phi sampling intervals in degrees).

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
- `field2wavecoeffs` — scan -> spherical-wave coefficients `T_{smn}`.
- `wavecoeffs2farfield_uniform` — `T_{smn}` -> far-field on a uniform
  (dtheta, dphi) grid (dipole output probe).
- `wavecoeffs2farfield` — `T_{smn}` -> far-field at arbitrary paired
  (theta, phi) directions (no probe normalization).

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

A set of companion markdown summaries of Hansen Chapters 2-4 and
Appendices A1-A5 is maintained alongside this code under
`Hansen_SphericalNearFieldAntennaMeasurements/` (separate repository).
Cross-references in the source cite Hansen equation numbers (e.g.
`[1], (4.135)`) that resolve to those summaries.

## Self-test

```
python pysnf.py
```

runs a round-trip regression test and a single-mode purity test; both
should report `PASS` at a 1e-10 tolerance.
