# -*- coding: utf-8 -*-
# SPDX-License-Identifier: MIT
"""pysnf_gt.py -- G/T from spherical near-field data (CW-Ambient technique).

Helpers built on `pysnf` for the method of

    R. T. Cutshall, J. Dobbins and M. N. Barr, "Measuring G/T with a
    Spherical Near-Field Antenna Measurement System via the CW-Ambient
    Technique," AMTA Proceedings.                       (cited as [GT])

Functions
---------
g_over_t
    G/T of an active unit under test over the far-field sphere.
near_zone_correction
    The near-zone correction factor N_SNF of [GT] Eq. (24), which replaces
    (lambda*d)^2 in the far-field G/T and Friis expressions.

The measurements ([GT] Sec. III.B)
----------------------------------
1. ``p_in``      power into the probe, at the probe's gain reference plane.
2. ``p_ambient`` noise power out of the UUT, probe not radiating (P1).
3. ``p_cw``      power out of the UUT with the probe radiating a CW tone
                 from location P = (A, theta_o, phi_o, chi_o) (P2).
4. ``nearfield`` a full dual-polarized SNF scan of the UUT, as complex
                 receiver values, that contains the pose P.

The expression used
-------------------
Written in terms of the complex receiver value at P, ``w(P)``::

    G/T (theta, phi) = k B (X - 1) / (P_in |w(P)|^2) * |F(theta, phi)|^2

with ``X = p_cw / p_ambient``, ``k`` Boltzmann's constant, ``B`` the noise
bandwidth and ``F`` the gain-normalized far field returned by
``pysnf.nearfield2farfield(..., normalization='gain',
probe_normalization='gain')`` from a probe pattern that carries the probe's
gain (``|K_p|^2`` = gain).  Three points where this differs in form, not in
substance, from [GT] Eqs. (19) and (23):

* **No separate G_p.**  [GT] assumes a probe pattern normalized to 0 dB and
  divides by the probe gain G_p.  Here the probe gain is carried inside the
  probe pattern, so the transform has already divided it out.  Pass
  ``probe_gain`` when your probe pattern is not gain-normalized; the pattern
  is then rescaled so that its peak equals that gain.
* **No RMS / peak switch.**  [GT] (19) uses ``P_RCVR = |w|^2 / 2`` and (23)
  uses ``|w|^2``; both reduce to the expression above once written in terms
  of the complex value ``w(P)`` taken from the scan itself.  Any constant
  complex receiver scale factor cancels between ``|F|^2`` and ``|w(P)|^2``.
  ``p_in`` must be true (RMS) power, as a power meter reports it.
* **No transform scale factor.**  With both normalizations set to ``'gain'``
  PySNF adds no scaling of its own (see the PySNF README).

Assumptions inherited from the CW-Ambient technique ([GT] Sec. III.C-D and
its reference [2]); none can be checked here: the UUT noise dominates the
measurement device's own noise (k T_M B negligible), lines and devices are
matched, cable / attenuator losses are accounted for in the supplied
powers, the ambient brightness temperature is stable, and ``p_ambient`` and
``p_cw`` are measured with the same device, bandwidth and averaging.

Conventions are PySNF's: exp(-i*omega*t), angles in degrees at the API,
``'E_theta'`` = probe signal at chi = 0, ``'E_phi'`` = signal at chi = pi/2.
"""

import numpy as np

import pysnf as sf

__all__ = ['g_over_t', 'near_zone_correction', 'BOLTZMANN']

BOLTZMANN = 1.380649e-23      # J/K (exact, SI 2019)

_POL_KEY = {'theta': 'E_theta', 'phi': 'E_phi'}


def _reference_value(nearfield, ref_index, ref_pol, w_ref):
    """Complex receiver value w(P) at the CW-measurement pose."""
    if w_ref is not None:
        if ref_index is not None:
            raise ValueError("give either ref_index or w_ref, not both")
        return complex(w_ref)
    if ref_index is None:
        raise ValueError("ref_index (the scan sample taken at the CW pose P) "
                         "or w_ref is required")
    if ref_pol not in _POL_KEY:
        raise ValueError("ref_pol must be 'theta' (chi = 0) or 'phi' "
                         "(chi = pi/2); got {!r}".format(ref_pol))
    return complex(np.asarray(nearfield[_POL_KEY[ref_pol]])[ref_index])


def _gain_carrying_probe(probe, probe_gain):
    """Return (probe dict with |K_p|^2 = gain, peak probe gain G_p)."""
    if probe is None:
        # PySNF's built-in probe: an ideal Hertzian dipole, gain 1.5.
        if probe_gain is not None:
            raise ValueError("probe_gain cannot be set for the built-in dipole "
                             "probe (probe=None); its gain is 1.5")
        return None, 1.5
    e_th = np.asarray(probe['E_theta'])
    e_ph = np.asarray(probe['E_phi'])
    peak = float(np.max(np.abs(e_th)**2 + np.abs(e_ph)**2))
    if not peak > 0.0:
        raise ValueError("probe pattern is identically zero")
    if probe_gain is None:
        return probe, peak
    if not probe_gain > 0.0:
        raise ValueError("probe_gain must be positive (linear, not dB)")
    scale = np.sqrt(probe_gain/peak)
    scaled = dict(probe)
    scaled['E_theta'] = scale*e_th
    scaled['E_phi'] = scale*e_ph
    return scaled, float(probe_gain)


def _transform(nearfield, probe, ka, dtheta_out, dphi_out, probe_gain,
               ref_index, ref_pol, w_ref, transform_kwargs):
    """Gain-normalized far field, w(P), and the peak probe gain."""
    for reserved in ('normalization', 'probe_normalization'):
        if reserved in transform_kwargs:
            raise TypeError("{} is fixed to 'gain' by this function".format(reserved))
    if dtheta_out is None or dphi_out is None:
        if 'dT' not in nearfield or 'dP' not in nearfield:
            raise ValueError("dtheta_out and dphi_out are required for a "
                             "latitude / MFC scan")
        dtheta_out = nearfield['dT'] if dtheta_out is None else dtheta_out
        dphi_out = nearfield['dP'] if dphi_out is None else dphi_out

    w_p = _reference_value(nearfield, ref_index, ref_pol, w_ref)
    if not abs(w_p) > 0.0:
        raise ValueError("the reference value w(P) is zero; choose a pose and "
                         "probe roll with signal (the co-polarized one)")
    probe, g_p = _gain_carrying_probe(probe, probe_gain)
    ff = sf.nearfield2farfield(nearfield, dtheta_out, dphi_out, probe=probe, ka=ka,
                               normalization='gain', probe_normalization='gain',
                               **transform_kwargs)
    return ff, w_p, g_p


def _copol_power(e_theta, e_phi, copol):
    """|F . conj(i_co)|^2 for a constant co-polarization vector."""
    c_th, c_ph = complex(copol[0]), complex(copol[1])
    norm = np.sqrt(abs(c_th)**2 + abs(c_ph)**2)
    if not norm > 0.0:
        raise ValueError("copol must be a non-zero (theta, phi) pair")
    return np.abs(e_theta*np.conj(c_th/norm) + e_phi*np.conj(c_ph/norm))**2


def g_over_t(nearfield, probe, ka, p_in, p_ambient, p_cw, bandwidth,
             ref_index=None, ref_pol='theta', w_ref=None, probe_gain=None,
             dtheta_out=None, dphi_out=None, copol=None, **transform_kwargs):
    """G/T of an active UUT from SNF data, by the CW-Ambient technique.

    Parameters
    ----------
    nearfield : dict
        Full dual-polarized SNF scan of the UUT as complex receiver values,
        in the layout `pysnf.nearfield2farfield` takes (``'E_theta'``,
        ``'E_phi'``, ``'dT'``, ``'dP'``; or the flat latitude / MFC layout
        with ``sampling=`` passed through).  Absolute calibration of the
        receiver is not needed.
    probe : dict or None
        Probe far-field pattern, same dict layout as a uniform scan.  It
        must be gain-normalized (``|E_theta|^2 + |E_phi|^2`` = probe gain)
        unless ``probe_gain`` is given.  ``None`` selects PySNF's ideal
        Hertzian-dipole probe (gain 1.5), available at ``ka = inf`` only.
    ka : float
        Probe distance times wavenumber, ``k*A``.
    p_in : float
        Power into the probe during the CW measurement, in **watts**, at
        the reference plane to which the probe gain is defined.
    p_ambient, p_cw : float
        The ambient (P1) and CW (P2) powers out of the UUT.  Only their
        ratio is used, so any common linear unit will do.
    bandwidth : float
        Noise bandwidth of the P1 / P2 measurement, in hertz.
    ref_index : int or tuple of int
        Index into the scan arrays of the sample taken at the CW pose P --
        ``(i_theta, i_phi)`` for a uniform scan, a flat index for a
        latitude / MFC scan.
    ref_pol : {'theta', 'phi'}
        Which scan array holds the probe roll used for the CW measurement:
        ``'theta'`` for chi = 0, ``'phi'`` for chi = pi/2.
    w_ref : complex, optional
        The receiver value at P, given directly instead of ``ref_index``
        (for a pose that is not a scan sample).  It must be measured with
        the same receiver configuration as the scan.
    probe_gain : float, optional
        Peak gain of the probe, **linear**.  When given, the probe pattern
        is rescaled so that its peak ``|E|^2`` equals this value; use it
        for a pattern normalized to 0 dB or otherwise arbitrarily scaled.
    dtheta_out, dphi_out : float, optional
        Output grid spacing in degrees.  Default: the scan's own spacing
        (required for latitude / MFC scans).
    copol : pair of complex, optional
        Constant co-polarization vector ``(c_theta, c_phi)`` in the far
        field, e.g. ``(1, 0)`` for theta, ``(1, 1j)`` or ``(1, -1j)`` for
        the two circular senses under exp(-i*omega*t).  Normalized
        internally.  Adds ``'G_T_copol'`` to the result.
    **transform_kwargs
        Passed to `pysnf.nearfield2farfield`: ``probe_pol``, ``n_max``,
        ``m_max``, ``nu_max``, ``dtype``, ``sampling``.

    Returns
    -------
    result : dict
        ``'G_T_theta'``, ``'G_T_phi'`` : partial G/T for the two far-field
        components, linear, in 1/K, shape (numthetas, numphis).
        ``'G_T'`` : their sum (total G/T).
        ``'G_T_copol'`` : G/T for the ``copol`` vector, when requested.
        ``'E_theta'``, ``'E_phi'`` : the complex gain-normalized far field
        ``F``, so that ``scale*|F . conj(i_co)|^2`` gives G/T for any other
        polarization.
        ``'scale'`` : ``k B (X - 1) / (P_in |w(P)|^2)``.
        ``'y_factor'`` : ``X = p_cw / p_ambient``.
        ``'dT'``, ``'dP'`` : output spacing, degrees.

        Use ``10*log10`` of a G/T entry for dB/K.
    """
    if not (p_in > 0.0 and p_ambient > 0.0 and bandwidth > 0.0):
        raise ValueError("p_in, p_ambient and bandwidth must be positive "
                         "(linear units, not dB)")
    y_factor = float(p_cw)/float(p_ambient)
    if not y_factor > 1.0:
        raise ValueError("p_cw must exceed p_ambient (X = {:.6g}); the CW "
                         "signal is not above the noise".format(y_factor))

    ff, w_p, _ = _transform(nearfield, probe, ka, dtheta_out, dphi_out, probe_gain,
                            ref_index, ref_pol, w_ref, transform_kwargs)

    scale = BOLTZMANN*bandwidth*(y_factor - 1.0)/(p_in*abs(w_p)**2)
    gt_theta = scale*np.abs(ff['E_theta'])**2
    gt_phi = scale*np.abs(ff['E_phi'])**2
    result = {'G_T_theta': gt_theta, 'G_T_phi': gt_phi, 'G_T': gt_theta + gt_phi,
              'E_theta': ff['E_theta'], 'E_phi': ff['E_phi'],
              'scale': scale, 'y_factor': y_factor,
              'dT': ff['dT'], 'dP': ff['dP']}
    if copol is not None:
        result['G_T_copol'] = scale*_copol_power(ff['E_theta'], ff['E_phi'], copol)
    return result


def near_zone_correction(nearfield, probe, ka, wavelength,
                         ref_index=None, ref_pol='theta', w_ref=None,
                         probe_gain=None, dtheta_out=None, dphi_out=None,
                         copol=None, **transform_kwargs):
    """Near-zone correction factor N_SNF(P) of [GT] Eq. (24).

    ``N`` replaces ``(lambda*d)^2`` in the far-field G/T expression and in
    the Friis equation, for the probe at the scan pose P and the UUT gain
    taken in the direction ``(theta, phi)``::

        P_out(P) / P_in = G_p * G(theta, phi) * lambda^4 / ((4 pi)^2 N(theta, phi))

    In PySNF terms ``N = (lambda^2 / 4 pi)^2 * G_p * |F|^2 / |w(P)|^2``.
    As ``ka -> inf`` with (theta, phi) the direction of P, ``N`` tends to
    ``(lambda*A)^2``.

    Parameters are as for `g_over_t`, plus:

    wavelength : float
        Free-space wavelength; ``N`` is returned in that unit to the fourth
        power.

    Returns
    -------
    result : dict
        ``'N_theta'``, ``'N_phi'`` (for the two far-field components),
        ``'N'`` (for the total gain), ``'N_copol'`` when ``copol`` is given,
        ``'probe_gain'`` (the peak G_p used), ``'dT'``, ``'dP'``.
    """
    if not wavelength > 0.0:
        raise ValueError("wavelength must be positive")
    ff, w_p, g_p = _transform(nearfield, probe, ka, dtheta_out, dphi_out, probe_gain,
                              ref_index, ref_pol, w_ref, transform_kwargs)
    scale = (wavelength**2/(4.0*np.pi))**2*g_p/abs(w_p)**2
    n_theta = scale*np.abs(ff['E_theta'])**2
    n_phi = scale*np.abs(ff['E_phi'])**2
    result = {'N_theta': n_theta, 'N_phi': n_phi, 'N': n_theta + n_phi,
              'probe_gain': g_p, 'dT': ff['dT'], 'dP': ff['dP']}
    if copol is not None:
        result['N_copol'] = scale*_copol_power(ff['E_theta'], ff['E_phi'], copol)
    return result


# ---------------------------------------------------------------------------
# Self-test:  python pysnf_gt.py
# ---------------------------------------------------------------------------
if __name__ == '__main__':
    import os
    import sys

    results = []

    def _check(name, passed, detail):
        results.append((name, bool(passed)))
        print('  {}  {}'.format('PASS' if passed else 'FAIL', detail))

    # A synthetic CW-Ambient measurement.  Everything the method is supposed
    # to be blind to is given an arbitrary value: electronics gain, receiver
    # scale factor (magnitude and phase), measurement-device gain.
    T_SYS = 437.0                    # K   -- the answer is G / T_SYS
    BW = 2.0e6                       # Hz
    P_IN = 3.0e-4                    # W into the probe
    G_ELEX = 10.0**(31.0/10.0)       # UUT electronics gain
    G_M = 10.0**(-4.0/10.0)          # power meter / analyzer gain
    A_SF = 0.37*np.exp(0.7j)         # uncalibrated receiver scale factor
    V = np.sqrt(2.0*P_IN)            # probe excitation, P_in = |v|^2 / 2

    def _measure(w_unit, ref_index, ref_pol):
        """Receiver scan and the two powers, from the v = 1 probe signal."""
        scan = {k: V*np.sqrt(G_ELEX)*A_SF*w_unit[k] for k in ('E_theta', 'E_phi')}
        g_snf = abs(w_unit[_POL_KEY[ref_pol]][ref_index])**2     # [GT] Eq. (8)
        p1 = BOLTZMANN*T_SYS*BW*G_ELEX*G_M                       # [GT] Eq. (12)
        p2 = P_IN*g_snf*G_ELEX*G_M + p1                          # [GT] Eq. (13)
        return scan, p1, p2

    # ------------------------------------------------------------------
    # Test 1: Hertzian dipoles, closed form ([GT] Appendix).
    # AUT: z-directed dipole, gain 1.5.  Probe: x'-polarized dipole, gain
    # 1.5.  The answer is known exactly: peak G/T = 1.5 / T_sys at
    # theta = 90 deg, at any distance, for any receiver scale factor, and
    # whether the probe gain is carried by the pattern or passed separately.
    # ------------------------------------------------------------------
    print('--- Test 1: Hertzian-dipole UUT, closed-form G/T = 1.5 / T_sys ---')
    dT = dP = 5.0
    th = np.linspace(0.0, np.pi, int(round(180.0/dT)) + 1)
    ph = np.linspace(0.0, 2.0*np.pi, int(round(360.0/dP)), endpoint=False)
    TH, PH = np.meshgrid(th, ph, indexing='ij')
    probe27 = {'E_theta': np.sqrt(6.0)/2.0*np.cos(TH)*np.cos(PH),     # [GT] Eq. (27)
               'E_phi': -np.sqrt(6.0)/2.0*np.sin(PH), 'dT': dT, 'dP': dP}
    probe_0db = dict(probe27, E_theta=probe27['E_theta']/np.sqrt(1.5),
                     E_phi=probe27['E_phi']/np.sqrt(1.5))             # peak = 0 dB
    i90 = int(np.argmin(np.abs(th - np.pi/2.0)))
    kw = dict(n_max=3, m_max=3, nu_max=1)
    expected = 1.5/T_SYS
    worst = 0.0
    for ka in (3.0, 10.0, 100.0):
        w26 = 0.75*np.exp(1j*ka)*(-1j/ka + 1.0/ka**2 + 1j/ka**3)*np.sin(TH)  # Eq. (26)
        unit = {'E_theta': w26, 'E_phi': np.zeros_like(w26)}
        for ref in ((i90, 0), (i90 - 6, 17)):          # on and off the peak
            scan, p1, p2 = _measure(unit, ref, 'theta')
            scan.update(dT=dT, dP=dP)
            for prb, g in ((probe27, None), (probe_0db, 1.5)):
                gt = g_over_t(scan, prb, ka, P_IN, p1, p2, BW, ref_index=ref,
                              probe_gain=g, **kw)
                row = int(np.unravel_index(np.argmax(gt['G_T']), gt['G_T'].shape)[0])
                err = abs(gt['G_T'].max()/expected - 1.0)
                worst = max(worst, err if row == i90 else 1.0)
    _check('Dipole G/T', worst < 1e-10,
           'peak G/T = {:.4f} dB/K expected; worst relative error {:.2e} over '
           'kA = 3, 10, 100, two CW poses, two probe-gain conventions '
           '(tolerance 1e-10)'.format(10.0*np.log10(expected), worst))

    # Far-field limit of the near-zone correction: N -> (lambda*A)^2.
    ka_far, lam = 1.0e4, 0.03
    w26 = 0.75*np.exp(1j*ka_far)*(-1j/ka_far + 1.0/ka_far**2 + 1j/ka_far**3)*np.sin(TH)
    scan = {'E_theta': w26, 'E_phi': np.zeros_like(w26), 'dT': dT, 'dP': dP}
    nz = near_zone_correction(scan, probe27, ka_far, lam, ref_index=(i90, 0), **kw)
    ratio = nz['N'][i90, 0]/(lam*ka_far*lam/(2.0*np.pi))**2
    _check('Near-zone factor, far-field limit', abs(ratio - 1.0) < 1e-6,
           'N / (lambda*A)^2 = {:.9f} at kA = 1e4  (tolerance 1e-6)'.format(ratio))

    # ------------------------------------------------------------------
    # Test 2: a real antenna in the near field.
    # The shipped WR-284 horn (N = 40) scanned at kA = 100 with a dipole
    # probe via pysnf.transmission_formula.  The whole G/T pattern must
    # equal |K_horn(theta, phi)|^2 / T_sys, with K from the horn's own
    # coefficients, and must not depend on which scan sample is used as
    # the CW pose P.
    # ------------------------------------------------------------------
    print('--- Test 2: WR-284 horn at kA = 100, full G/T pattern ---')
    csv_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'wr284_horn_Tsmn.csv')
    if os.path.exists(csv_path):
        raw = np.loadtxt(csv_path, delimiter=',', skiprows=5)
        n_i, m_i, s_i = (raw[:, k].astype(int) for k in range(3))
        n_h = int(n_i.max())
        T_horn = np.zeros((n_h, 2*n_h + 1, 2), dtype='complex128')
        T_horn[n_i - 1, m_i + n_h, s_i - 1] = raw[:, 3] + 1j*raw[:, 4]

        T_probe = np.zeros((1, 3, 2), dtype='complex128')         # +x dipole
        T_probe[0, 0, 1] = +np.sqrt(2.0)/2.0
        T_probe[0, 2, 1] = -np.sqrt(2.0)/2.0
        r_p = np.swapaxes(sf.reciprocity(sf.rotate_wavecoeffs_about_axis(T_probe, 'x')), 0, 2)

        ka, dT = 100.0, 2.0
        th = np.linspace(0.0, np.pi, int(round(180.0/dT)) + 1)
        ph = np.linspace(0.0, 2.0*np.pi, int(round(360.0/dT)), endpoint=False)
        TH, PH = np.meshgrid(th, ph, indexing='ij')
        unit = {'E_theta': sf.transmission_formula(T_horn, 0.0, TH, PH, r_p=r_p, ka=ka),
                'E_phi': sf.transmission_formula(T_horn, np.pi/2.0, TH, PH, r_p=r_p, ka=ka)}
        pe_th, pe_ph = sf.wavecoeffs2farfield_uniform(T_probe, dT, dT, normalization='gain')
        probe = {'E_theta': pe_th, 'E_phi': pe_ph, 'dT': dT, 'dP': dT}

        k_th, k_ph = sf.wavecoeffs2farfield(T_horn, TH.ravel(), PH.ravel())
        ref_th = (np.abs(k_th)**2).reshape(TH.shape)/T_SYS
        ref_ph = (np.abs(k_ph)**2).reshape(TH.shape)/T_SYS
        peak = float((ref_th + ref_ph).max())

        # CW pose: the strongest scan sample, then a weak one 30 deg off.
        both = np.stack([np.abs(unit['E_theta']), np.abs(unit['E_phi'])])
        ip, it, iph = np.unravel_index(np.argmax(both), both.shape)
        poses = [((int(it), int(iph)), ('theta', 'phi')[ip]),
                 ((int(it) + 15, (int(iph) + 20) % ph.size), ('theta', 'phi')[ip])]
        worst, peaks = 0.0, []
        for ref, pol in poses:
            scan, p1, p2 = _measure(unit, ref, pol)
            scan.update(dT=dT, dP=dT)
            gt = g_over_t(scan, probe, ka, P_IN, p1, p2, BW, ref_index=ref, ref_pol=pol,
                          n_max=n_h, m_max=n_h, nu_max=1, copol=(1.0, 0.0))
            err = max(np.max(np.abs(gt['G_T_theta'] - ref_th)),
                      np.max(np.abs(gt['G_T_phi'] - ref_ph)),
                      np.max(np.abs(gt['G_T_copol'] - ref_th)))/peak
            worst = max(worst, float(err))
            peaks.append(10.0*np.log10(gt['G_T'].max()))
        _check('Horn G/T pattern', worst < 1e-9,
               'peak G/T = {:.4f} dB/K (expected {:.4f}); worst pattern error / peak = '
               '{:.2e} over two CW poses, both components and copol  '
               '(tolerance 1e-9)'.format(peaks[0], 10.0*np.log10(peak), worst))
    else:
        print('  SKIPPED (wr284_horn_Tsmn.csv not found next to pysnf_gt.py)')

    # ------------------------------------------------------------------
    # Input validation
    # ------------------------------------------------------------------
    print('--- Input validation ---')
    scan0 = {'E_theta': np.ones((37, 72), complex), 'E_phi': np.zeros((37, 72), complex),
             'dT': 5.0, 'dP': 5.0}
    bad = {'CW not above ambient': dict(p_ambient=2.0, p_cw=1.0, ref_index=(18, 0)),
           'no CW pose given': dict(p_ambient=1.0, p_cw=2.0),
           'zero reference value': dict(p_ambient=1.0, p_cw=2.0, ref_index=(18, 0), ref_pol='phi'),
           'bad ref_pol': dict(p_ambient=1.0, p_cw=2.0, ref_index=(18, 0), ref_pol='x')}
    raised = 0
    for label, args in bad.items():
        try:
            g_over_t(scan0, probe27, 10.0, 1.0e-3, bandwidth=1.0e6, **args)
        except ValueError:
            raised += 1
    _check('Validation', raised == len(bad),
           '{} of {} bad inputs rejected with ValueError'.format(raised, len(bad)))

    n_fail = sum(1 for _, ok in results if not ok)
    print('=== {} of {} checks passed ==='.format(len(results) - n_fail, len(results)))
    sys.exit(1 if n_fail else 0)
