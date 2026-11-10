# Hansen → PySNF cross-reference index

This is a cross-reference index only: it maps equation numbers from

> [1] J. E. Hansen (ed.), *Spherical Near-Field Antenna Measurements*,
> IEE Electromagnetic Waves Series, vol. 26, Peter Peregrinus Ltd., 1988

to the PySNF function (or source location) that uses them.  No
formulas, tables, or figures from Hansen's book are reproduced here.
Readers who want the equations themselves should consult a copy of
Hansen directly.

## Chapter 2 — Scattering-matrix description of a single antenna

| Hansen eq. | PySNF function | Purpose |
|------------|----------------|---------|
| (2.19)     | `incomplete_farfield_pattern_functions`, `rotation_coefficients` | `(-m/m)**m` sign convention at `m = 0` |
| (2.104)    | `reciprocity`  | Transmitting ↔ receiving coefficient conversion |
| (2.124)-(2.125) | `pysnf.py` module docstring / self-test | Ideal Hertzian dipole at the origin |
| (2.154)-(2.155) | `dipole_probe_response_constants` | `R^p` for `+x` / `+y` electric dipoles |
| (2.157)-(2.158) | `dipole_probe_response_constants` | `R^p` for `+x` / `+y` magnetic dipoles |

## Chapter 3 — Antenna coupling

| Hansen eq. | PySNF function | Purpose |
|------------|----------------|---------|
| (3.10)     | `transmission_formula`, `_transmission_formula_kbased` | Spherical near-field transmission formula |
| (3.19)     | `translation_coefficients` | Series form of the translation coefficients |
| (3.26)     | `probe_response_constants` | Build `P_{s mu n}(kA)` from `R^p` coefficients |
| (3.27)-(3.28), (3.32), (3.36) | `transmission_formula` | x-polarized-probe reduction of (3.10) |
| (3.52)-(3.54) | `wavecoeffs2farfield_uniform` | Dipole output-probe signal at `chi = 0, pi/2` |

## Chapter 4 — Data reduction in spherical near-field measurements

| Hansen eq. | PySNF function | Purpose |
|------------|----------------|---------|
| (4.84), (4.86) | `pi_wiggle` | `pi_l` coefficients for the `K(m')` convolution |
| (4.87)     | `b_wiggle`     | Periodic embedding of `b_{l, m, mu}` |
| (4.89)     | `field2w_n_m_mu` | FFT-based `K(m')` convolution (typo-corrected) |
| (4.92)     | `field2w_n_m_mu` | Assembly of `w_{nm}^{(mu)}` |
| (4.100)-(4.102) | `probe_response_constants` | Far-field limit (`ka = inf`) of `P_{s mu n}` |
| (4.103)-(4.104) | `dipole_probe_response_constants` | Far-field dipole `P_{s mu n}` (reference values) |
| (4.126)-(4.128) | `field2w_n_m_mu` | Intermediate quantities in the data-reduction chain |
| (4.133)-(4.135) | `field2wavecoeffs`, `wavecoeffs2farfield_uniform` | Solve for `T_{smn}`; far-field M summation |
| pg. 192 symmetry | `singlesphere2doublesphere` | Single- → double-sphere scan extension |

## Chapter 5 — Rotation of spherical waves (rigid-body rotations of `T_{smn}`)

| Hansen eq. | PySNF function | Purpose |
|------------|----------------|---------|
| (5.67)     | `rotate_wavecoeffs_about_axis` | Rotate `T_{smn}` about a coordinate axis |

## Appendix A1 — Spherical wave functions

| Hansen eq. | PySNF function | Purpose |
|------------|----------------|---------|
| (A1.25)    | `lpmn_norm`    | Normalized associated Legendre functions |
| (A1.36), (A1.38) | `lpmn_norm` | Closed-form pole values of the Legendre functions |
| (A1.42), (A1.44) | `lpmn_norm` | Closed-form pole values of their theta-derivatives |
| (A1.59)-(A1.60) | `incomplete_farfield_pattern_functions`, `farfield_pattern_functions`, `wavecoeffs2farfield` | Far-field pattern functions `K_{smn}` |
| (A1.61)-(A1.64) | `incomplete_farfield_pattern_functions` | Closed-form pole values of `K_{smn}` |

## Appendix A2 — Rotation of spherical waves

| Hansen eq. | PySNF function | Purpose |
|------------|----------------|---------|
| (A2.11)    | `rotation_coefficients` | General Fourier expansion of `d^n_{mu m}` |
| (A2.17)-(A2.19) | `rotation_coefficients` | Closed-form `d^n_{mu m}` for `mu in {-1, 0, +1}` |
| (A2.26), (A2.28), (A2.30), (A2.32) | `delta_pyramid` | Delta-pyramid symmetry relations |
| (A2.35)    | `delta_pyramid` | Delta-pyramid recurrence |
| (A2.41)    | `delta_pyramid` | Delta-pyramid Region-I bottom edge |

## Appendix A3 — Translation of spherical waves

| Hansen eq. | PySNF function | Purpose |
|------------|----------------|---------|
| (A3.7)     | `translation_coefficients` | Parity / triangle validity rule |
| (A3.8), (A3.9) | `translation_coefficients` | Elementary block symmetries |
| (A3.21)    | `translation_coefficients` | `mu = -1` column from `mu = +1` |
| (A3.23)-(A3.24) | `_translation_coefficients_far_field` | Asymptotic (`ka = inf`) translation coefficients |
| (A3.27)    | `translation_coefficients` | `B'(J)` recurrence prefactor |

## Appendix A4 — Data processing in antenna measurements

Used for the DFT conventions in `_fft` / `_ifft` and the `interpft`
band-limited resampler.  No individual equations are cross-referenced;
see the module docstring of `pysnf.py` under "Conventions / DFT" for
the sign convention, and `interpft` for the interpolation basis.

The same interpolation (Section A4.4.3) is what brings latitude / MFC
scans back to the uniform grid in `field2wavecoeffs(..., sampling=...)`
(internal helper `_rings2uniform`).  Thinning the phi samples near the
poles is discussed on pp. 128-130; the per-ring count rules are
`sin_theta_phi_counts` and `mfc_phi_counts`, and `sampling_grid`
generates the sample locations.  The MFC rule itself is not from Hansen;
see reference [2] in the `pysnf.py` module docstring.

## Appendix A5 — List of principal symbols

Hansen's symbol glossary is not reproduced here.  The
"Symbol cross-reference (Hansen Appendix A5 → PySNF)" section of the
`pysnf.py` module docstring maps each Hansen symbol used in this
library onto its PySNF array name and shape.
