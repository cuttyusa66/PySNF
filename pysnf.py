# -*- coding: utf-8 -*-
""" PySNF Version 0.1

A spherical near-field antenna measurement data-reduction library
implementing the SNIFT-style transformation algorithm from Hansen [1],
Chapter 4 / Section 4.4.

References
----------
[1] J.E. Hansen, Spherical Near-Field Antenna Measurements, 1988
    (IEE Electromagnetic Waves Series, vol. 26)

Conventions
-----------

Time convention
    The library uses the physics/engineering convention exp(-i*omega*t)
    throughout, matching Hansen [1].  Outgoing waves therefore carry
    h_n^{(1)}(kr) ~ e^{+ikr}/(kr) at large kr.  This shows up explicitly
    in factors like ``1j**(-n_array)`` in the dipole and far-field probe
    response constants.  If you port code from a convention that uses
    exp(+i*omega*t), invert the sign of every imaginary exponent.

Wave-coefficient array layout
    Test-antenna transmitting coefficients T_{s, m, n} are stored in
    the array ``q_n_m_s`` with shape ``(n_max, 2*m_max+1, 2)`` and
    indexing convention::

        q_n_m_s[n-1, m + m_max, s-1]  ==  T_{s, m, n}

    where ``s in {1, 2}`` (TE / TM coefficient block per Hansen
    Section 2.2.2), ``n in {1, ..., n_max}``, and
    ``m in {-m_max, ..., m_max}``.

DFT / IDFT convention
    Hansen [1] uses the convention (Appendix A4)::

        DFT  : G(l) = sum_m  g(m) * eps_J^{ m * l}      (+i exponent, no 1/J)
        IDFT : g(m) = (1/J) sum_l G(l) * eps_J^{-l * m} (-i exponent, with 1/J)

    where ``eps_J = exp(+2*pi*i/J)``.  This is the *opposite* sign of
    np.fft.  The module-internal helpers ``_fft`` and ``_ifft`` adapt
    numpy's routines to Hansen's convention::

        _fft(g)  == J * np.fft.ifft(g)
        _ifft(G) == np.fft.fft(G) / J

    All internal pipeline code uses these wrappers (rather than np.fft
    directly) to remain consistent with the formulas in [1].  The leading
    underscore marks them as internal helpers; external callers should
    interact with the library via `nearfield2farfield`, `field2wavecoeffs`,
    and `wavecoeffs2farfield_uniform`.

Angle convention
    Public-facing functions accept theta- and phi-spacing parameters
    (``dtheta``, ``dphi``, etc.) in **degrees**.  Internally,
    spherical-wave-function evaluation uses **radians**.

Probe convention
    A "+x" linearly-polarized probe with mu = +/- 1 is assumed unless
    a custom probe is supplied.  See ``dipole_probe_response_constants``
    and ``probe_response_constants`` for the precise definitions.
"""

# ---------------------
# IMPORTS
# ---------------------
import math
import numpy as np
import scipy.special as sps
import fractions
# Note: matplotlib is intentionally NOT imported at module load.  It is
# lazy-imported inside ``plot_spherical_wave_coefficients_mag_db`` so that
# users on headless systems (or those who never plot) don't pay the
# matplotlib import cost or need it installed.

# ---------------------
# GLOBAL VARIABLES
# ---------------------
INF = np.inf
PI = np.pi

# ---------------------
# FUNCTIONS AND METHODS
# ---------------------
# ------------------------------------------------------------------------


def nearfield2farfield(nearfield, dtheta_out, dphi_out, probe=None, probe_pol='x', ka=INF,
                       n_max=None, m_max=None, nu_max=None, dtype='complex128'):
    """Compute the far-field of the UUT from a sampled near-field scan.

    Parameters
    ----------
    nearfield : dict
        Near-field measurement data with keys:

        * ``'E_theta'`` : complex 2-D array of shape (numthetas, numphis)
            The theta-polarization probe signal (probe roll chi = 0).
            This is the polarization component, *not* the theta-axis grid.
        * ``'E_phi'`` : complex 2-D array of shape (numthetas, numphis)
            The phi-polarization probe signal (probe roll chi = pi/2).
        * ``'dT'``, ``'dP'`` : float
            Theta and phi sampling intervals, in degrees.
    dtheta_out, dphi_out : float
        Output sampling intervals in degrees.
    probe : dict or None
        Custom probe data with the same dictionary layout as `nearfield`
        ('E_theta', 'E_phi', 'dT', 'dP'); or None for an ideal Hertzian
        dipole.
    probe_pol : {'x', 'y'}
        Probe polarization (only 'x' supported for ka=INF dipole today).
    ka : float
        Measurement distance times wavenumber; use np.inf for far-field.
    n_max, m_max, nu_max : int or None
        Truncation indices; auto-detected from input sizes if None.
    dtype : {'complex128', 'complex64'}
        Precision used for the input and output arrays.  The internal
        computation always runs in complex128 for numerical headroom;
        this parameter controls only the storage precision at the API
        boundary, which is where large arrays live.  Use 'complex64' to
        halve memory for input scans and output patterns.

    Returns
    -------
    farfield : dict
        Far-field pattern with keys 'E_theta', 'E_phi', 'dT', 'dP' in
        the same layout as `nearfield`.

    Notes
    -----
    The 'E_theta' / 'E_phi' keys are *polarization components* of the
    measured probe signal, NOT spherical-coordinate grids.  The
    theta-axis grid is determined implicitly by ``dT`` (and the number
    of rows of E_theta); similarly for the phi-axis grid and ``dP``.
    """

    np_dtype = np.dtype(dtype)
    if np_dtype not in (np.dtype('complex64'), np.dtype('complex128')):
        raise ValueError(
            "dtype must be 'complex64' or 'complex128'; got {!r}".format(dtype))

    # Cast the input to the requested dtype at entry so the whole pipeline
    # propagates the user's precision.  All major allocations downstream
    # honor `dtype`; arithmetic with literal complex constants (e.g. 1j**n)
    # may briefly upcast to complex128, but the result is downcast back to
    # `dtype` at storage.
    nearfield_internal = nearfield
    if (nearfield['E_theta'].dtype != np_dtype
            or nearfield['E_phi'].dtype != np_dtype):
        nearfield_internal = nearfield.copy()
        nearfield_internal['E_theta'] = nearfield['E_theta'].astype(np_dtype)
        nearfield_internal['E_phi']   = nearfield['E_phi'].astype(np_dtype)

    # Calculate the spherical wave coefficients of the UUT
    q_n_m_s = field2wavecoeffs(nearfield_internal, probe, probe_pol, ka,
                               n_max, m_max, nu_max, dtype=dtype)

    # From the spherical wave coefficients, calculate the far-field
    E_theta, E_phi = wavecoeffs2farfield_uniform(q_n_m_s, dtheta_out, dphi_out,
                                                 dtype=dtype)

    # Organize into a dictionary
    farfield = {'E_theta': E_theta, 'E_phi': E_phi,
                'dT': dtheta_out, 'dP': dphi_out}

    return farfield


# ------------------------------------------------------------------------


def field2wavecoeffs(uut, probe=None, probe_pol='x', ka=INF,
                     n_max=None, m_max=None, nu_max=None, dtype='complex128'):

    np_dtype = np.dtype(dtype)

    # Get the values for w_n_m_mu
    w_n_m_mu = field2w_n_m_mu(uut, n_max, m_max, probe_pol, dtype=dtype)

    # If n_max was passed into this function as None, determine n_max
    if n_max is None:
        n_max = np.size(w_n_m_mu, 0)

    # If m_max was passed into this function as None, determine m_max
    if m_max is None:
        m_max = (np.size(w_n_m_mu, 1) - 1)//2

    # Get the probe response constants
    if (probe is None) and (probe_pol == 'x') and (ka == INF):
        p_n_mu_s = dipole_probe_response_constants(n_max)
    elif (probe is None) and (probe_pol == 'y') and (ka == INF):
        raise NotImplementedError("field2wavecoeffs: the infinitely remote "
                                  "y-polarized electric dipole option is not "
                                  "yet implemented.")
    elif (probe is None) and (probe_pol == 'x') and (ka < INF):
        raise NotImplementedError("field2wavecoeffs: the finitely remote "
                                  "x-polarized electric dipole option is not "
                                  "yet implemented.")
    elif (probe is None) and (probe_pol == 'y') and (ka < INF):
        raise NotImplementedError("field2wavecoeffs: the finitely remote "
                                  "y-polarized electric dipole option is not "
                                  "yet implemented.")
    else:
        # If nu_max was passed into this function as None, determine nu_max
        if nu_max is None:
            nu_max = np.size(probe['E_theta'], axis=0)//2
        t_p = field2wavecoeffs(probe, n_max=nu_max, m_max=1)  # t_p is nu, mu, sig
        t_p_rot = rotate_wavecoeffs_about_axis(t_p,'x')  # t_p_rot is nu, mu, sig
        r_p = reciprocity(t_p_rot)  # r_p is nu, mu, sig
        r_p = np.swapaxes(r_p, 0, 2)  # r_p is now sig, mu, nu
        p_n_mu_s = probe_response_constants(r_p, n_max, ka)

    # Initialize the wave coefficient matrix in the requested precision.
    q_n_m_s = np.zeros((n_max, 2*m_max+1, 2), dtype=np_dtype)

    # Pull out and reshape the necessary values of the probe response constants
    p_n_neg1_1 = np.reshape(p_n_mu_s[:, 0, 0], (n_max,1))  # mu = -1 , s = 1
    p_n_pos1_1 = np.reshape(p_n_mu_s[:, 1, 0], (n_max,1))  # mu = +1 , s = 1
    p_n_neg1_2 = np.reshape(p_n_mu_s[:, 0, 1], (n_max,1))  # mu = -1 , s = 2
    p_n_pos1_2 = np.reshape(p_n_mu_s[:, 1, 1], (n_max,1))  # mu = +1 , s = 2

    # Solve for q_n_m_s as described in [1],(4.133) and [1],(4.134), AND A FUTURE BLOG POST
    # Note: Give Jeffrey Hyman developer credit here.
    determinant = p_n_pos1_1*p_n_neg1_2 - p_n_neg1_1*p_n_pos1_2
    q_n_m_s[:, :, 0] = (p_n_neg1_2*w_n_m_mu[:, :, 1] - p_n_pos1_2*w_n_m_mu[:, :, 0])/determinant
    q_n_m_s[:, :, 1] = (p_n_pos1_1*w_n_m_mu[:, :, 0] - p_n_neg1_1*w_n_m_mu[:, :, 1])/determinant

    return q_n_m_s

# ------------------------------------------------------------------------


def wavecoeffs2farfield_uniform(q_n_m_s, dtheta, dphi, dtype='complex128'):
    """Synthesize the test-antenna far field from spherical wave coefficients.

    Limitation
    ----------
    The output probe is hardcoded to an infinitely remote, +x-polarized,
    electric Hertzian dipole.  Two outputs are returned: `theta_pol` (the
    normalized far-field signal in a theta-oriented short electric dipole,
    chi = 0) and `phi_pol` (the same for a phi-oriented short electric
    dipole, chi = pi/2).  Per Hansen, Spherical Near-Field Antenna
    Measurements, Eqs. (3.52)-(3.53) and Eq. (3.54), this is equivalent to
    sampling the test antenna's absolute far-field pattern K(theta, phi).

    If you need to predict the signal that would be measured by a different
    output probe (e.g. a calibrated horn), use the lower-level transmission
    formula machinery directly with custom probe response constants.

    Parameters
    ----------
    q_n_m_s : ndarray, shape (n_max, 2*m_max+1, 2), complex
        Test-antenna transmitting coefficients v * T_{s, m, n} indexed as
        q[n-1, m + m_max, s-1].
    dtheta, dphi : float
        Output sampling intervals in *degrees*.

    Returns
    -------
    theta_pol, phi_pol : ndarray of complex
        Normalized far-field probe signals on a uniform (theta, phi) grid
        with numthetas = int(180/dtheta + 1) theta samples in [0, pi]
        inclusive and numphis = int(360/dphi) phi samples in [0, 2*pi).
    """

    np_dtype = np.dtype(dtype)

    # Determine n_max and m_max from q_n_m_s
    n_max, m_max, s_max = np.shape(q_n_m_s)
    m_max = (m_max - 1)//2

    # Determine the number of theta and phis points required based on dT and dP
    numthetas = int(180.0/dtheta + 1)
    numphis = int(360.0/dphi)

    # Get the dipole probe response constants
    # (N x 2 x 2, where dim 0 = n, dim 1 = mu, dim 2 = s)
    p_n_mu_s = dipole_probe_response_constants(n_max)

    # Copy out the mu==1, s==1 probe response constants
    p_n = np.reshape(p_n_mu_s[:, 0, 0], (n_max,1,1))  # N x M x Theta

    # Calculate the rotation coefficients
    # (N x 3 x M x numThetas)
    d_n_mu_m = rotation_coefficients(n_max, m_max, 1, np.linspace(0, PI, numthetas))

    # Calculate the addition and subtraction of the rotation coefficients
    # (N x M x numThetas)
    dp1_plus_dm1 = d_n_mu_m[:, 2, :, :] + d_n_mu_m[:, 0, :, :]  # N x M x Theta
    dp1_minus_dm1 = d_n_mu_m[:, 2, :, :] - d_n_mu_m[:, 0, :, :]

    # Reshape the spherical wave coefficients to get ready for multiplication
    q_n_m_s = np.reshape(q_n_m_s, (n_max, 2*m_max+1, s_max, 1))  # N x M x S x Theta

    # ------------------------------------------------------------------
    # N summation per [1], (4.135), batched for chi=0 (theta-pol) and
    # chi=pi/2 (phi-pol).  The two cases share the same delta-product
    # structure but pair the q_s=1, q_s=2 modes with the opposite
    # rotation-coefficient sum/difference.
    # ------------------------------------------------------------------
    chi_0_inner  = q_n_m_s[:, :, 0, :]*dp1_plus_dm1  + q_n_m_s[:, :, 1, :]*dp1_minus_dm1
    chi_90_inner = q_n_m_s[:, :, 0, :]*dp1_minus_dm1 + q_n_m_s[:, :, 1, :]*dp1_plus_dm1
    stacked = np.stack([p_n * chi_0_inner,
                           1j * p_n * chi_90_inner], axis=0)  # (2, N, M, Theta)
    n_sums = np.sum(stacked, axis=1)                       # (2, M, Theta)
    n_sums = np.swapaxes(n_sums, 1, 2)                     # (2, Theta, M)

    # ------------------------------------------------------------------
    # M summation per [1], (4.135), via a single batched FFT.
    # Reorder the m axis with ifftshift so that bin 0 corresponds to m=0,
    # zero-pad to numphis, and FFT over phi.
    # ------------------------------------------------------------------
    n_sums = np.fft.ifftshift(n_sums, axes=2)
    # Allocate the zero-padded FFT input in the requested precision.
    # Assignments from `n_sums` (which may be complex128 due to literal
    # arithmetic) will be implicitly downcast here when np_dtype is
    # complex64, giving us the FFT speedup.
    temp = np.zeros((2, numthetas, numphis), dtype=np_dtype)
    temp[:, :, 0:m_max+1] = n_sums[:, :, 0:m_max+1]
    temp[:, :, numphis-m_max:] = n_sums[:, :, m_max+1:]
    fft_out = _fft(temp, axis=2)   # (2, numthetas, numphis)

    theta_pol = fft_out[0]
    phi_pol   = fft_out[1]

    return theta_pol, phi_pol

# ------------------------------------------------------------------------


def field2w_n_m_mu(uut, n_max=None, m_max=None, probe_pol='x', dtype='complex128'):

    np_dtype = np.dtype(dtype)

    # Make a copy of the UUT data so that we don't change the original data.
    # This is needed since dictionaries are mutable objects.
    aut = uut.copy()

    # AUT is a single sphere of data, i.e., 0 <= theta <= 180 and 0 <= phi < 360.
    # Use symmetry relationship in [1],Page 192 to create a double sphere of data,
    # i.e., 0 <= theta < 360 and 0 <= phi < 360.
    # Note: 'E_theta' / 'E_phi' are polarization components of the measured
    # probe signal, *not* spherical-coordinate grids.

    aut['E_theta'] = singlesphere2doublesphere(aut['E_theta'], dtype=np_dtype)
    aut['E_phi']   = singlesphere2doublesphere(aut['E_phi'],   dtype=np_dtype)

    num_th = np.size(aut['E_theta'], axis=0)
    num_ph = np.size(aut['E_theta'], axis=1)

    if n_max is None:
        n_max = num_th//2  # This works because num_th should always be even after
                           # the singlesphere2doublesphere function
    elif n_max > num_th//2:
        raise ValueError("n_max must be less than or equal to num_th/2")

    if m_max is None:
        m_max = num_ph//2  # This works because num_ph should always be even after
                           # the singlesphere2doublesphere function
    elif m_max > num_ph//2:
        raise ValueError("m_max must be less than or equal to num_ph/2")

    # Perform (4.126).  E_theta is the chi=0 probe signal, E_phi is chi=pi/2.
    w_th_ph_mu = np.zeros((num_th, num_ph, 2), dtype=np_dtype)
    if probe_pol == 'x':
        w_th_ph_mu[:, :, 1] = (1./2.)*(aut['E_theta'] - 1j*aut['E_phi'])  # mu = +1
        w_th_ph_mu[:, :, 0] = (1./2.)*(aut['E_theta'] + 1j*aut['E_phi'])  # mu = -1
    elif probe_pol == 'y':
        w_th_ph_mu[:, :, 1] = (1./2.)*(aut['E_phi'] - 1j*-aut['E_theta'])  # mu = +1
        w_th_ph_mu[:, :, 0] = (1./2.)*(aut['E_phi'] + 1j*-aut['E_theta'])  # mu = -1
    else:
        raise ValueError("probePol must either be 'x' or 'y'")

    # Perform (4.127)
    temp = _ifft(w_th_ph_mu, axis=1)

    # Reorganize such that m runs from -m_max to m_max
    w_th_m_mu = np.zeros((num_th, 2*m_max+1, 2), dtype=np_dtype)
    temp = np.fft.fftshift(temp, axes=(1,))
    if num_ph/2.0 == m_max:
        w_th_m_mu[:, :-1, :] = temp
        w_th_m_mu[:, -1, :] = temp[:, 0, :]
    elif num_ph/2.0 > m_max:
        w_th_m_mu[:, :, :] = temp[:, num_ph//2-m_max:num_ph//2+m_max+1, :]
    else:
        raise ValueError("Bad value for num_ph.")

    # Perform (4.128)
    temp = _ifft(w_th_m_mu, axis=0)

    # Reorganize such that n runs from -n_max to n_max
    b_l_m_mu = np.zeros((2*n_max+1, 2*m_max+1, 2), dtype=np_dtype)
    temp = np.fft.fftshift(temp, axes=(0,))
    if num_th/2.0 == n_max:
        b_l_m_mu[:-1, :, :] = temp
        b_l_m_mu[-1, :, :] = temp[0, :, :]
    elif num_th/2.0 > n_max:
        b_l_m_mu[:, :, :] = temp[num_th//2-n_max:num_th//2+n_max+1, :, :]
    else:
        raise ValueError("Bad value for num_th.")

    # Calculate the pi_wiggle array with [1],(4.84) and [1],(4.86)
    pi_wig = pi_wiggle(n_max)

    # Calculate the b_l_m_mu_wiggle array from b_l_m_mu using [1],(4.87)
    b_l_m_mu_wiggle = b_wiggle(b_l_m_mu)

    # Calculate k_mp with fast convolution via FFT methods as explained in [1],(4.89).
    # However, note that [1],(4.89) has a typo. If correct, [1],(4.89) should read:
    #
    #   K(m') = IDFT{ DFT{ PI_wiggle(i) | i = 0,1,...,4N-1 } *
    #                 DFT{ b_j_m_mu_wiggle | j = 0,1,...,4N-1 } }
    #
    k_mp = _ifft(_fft(pi_wig, axis=0) * _fft(b_l_m_mu_wiggle, axis=0), axis=0)

    # Keep only the values of k_mp where -n_max <= m' <= n_max.
    # This is required prior to the evaluation of [1],(4.92)
    temp = np.zeros((2*n_max+1, 2*m_max+1, 2), dtype=np_dtype)
    temp[0:n_max, :, :] = k_mp[3*n_max:, :, :]
    temp[n_max:, :, :] = k_mp[0:n_max+1, :, :]
    k_mp = temp

    # Pull out k_mp for mu == -1 and mu == +1 (2-D, shape (2*n_max+1, 2*m_max+1))
    k_mp_m1 = k_mp[:, :, 0]
    k_mp_p1 = k_mp[:, :, 1]

    # Initialize the n and m arrays
    n_array = np.reshape(np.linspace(1, n_max, n_max), (n_max, 1))
    m_array = np.reshape(np.linspace(-m_max, m_max, 2*m_max+1), (1, 2*m_max+1))

    # Initialize a delta pyramid helper function
    def get_mi(m_): return m_ + n_max  # This returns the m index of the deltas

    # Get the deltas for 1 <= n <= n_max and -m_max <= m <= m_max
    deltas = delta_pyramid(n_max)
    deltas = deltas[1:, :, get_mi(-m_max):get_mi(m_max)+1]

    # Update the delta pyramid helper function
    def get_mi(m_): return m_ + m_max  # This returns the m index of the revised deltas

    # Get the deltas for only mu == -1 or 1 (2-D, shape (n_max, 2*n_max+1))
    deltas_mu_m1 = deltas[:, :, get_mi(-1)]
    deltas_mu_p1 = deltas[:, :, get_mi(+1)]

    # Compute the m'-sum (axis "i" in the einsum) per [1],(4.92):
    #   result[n, m] = sum_{m'} Delta^n_{m', mu} * Delta^n_{m', m} * K(m')
    # Using einsum avoids materializing the (n, m', m) intermediate tensor.
    sum_m1 = np.einsum('ni,nim,im->nm', deltas_mu_m1, deltas, k_mp_m1)
    sum_p1 = np.einsum('ni,nim,im->nm', deltas_mu_p1, deltas, k_mp_p1)

    # Assemble w_n_m_mu per [1],(4.92).  Phase factor is i^(mu - m).
    w_n_m_mu = np.zeros((n_max, 2*m_max+1, 2), dtype=np_dtype)
    w_n_m_mu[:, :, 0] = ((2.0*n_array + 1.0)/2.0 * (1j**(-1 - m_array)) * sum_m1)
    w_n_m_mu[:, :, 1] = ((2.0*n_array + 1.0)/2.0 * (1j**(+1 - m_array)) * sum_p1)

    return w_n_m_mu


# ------------------------------------------------------------------------


def incomplete_farfield_pattern_functions(n_max, m_max, thetas, delta=1.0e-6):

    # Make a local copy so we don't mutate the caller's array when we
    # nudge values away from the poles.
    thetas = np.asarray(thetas, dtype=float).copy()

    # Check to make sure that thetas are between 0 and pi
    if np.any(thetas < 0) or np.any(thetas > PI):
        raise ValueError("thetas must all be greater than or equal to 0, and less than or equal to pi.")

    # If any theta value is less than the delta/10.0, add the delta
    zero_inds = thetas < delta/10.0
    thetas[zero_inds] += delta

    # If any theta value is within delta/10.0 of pi, subtract the delta
    pi_inds = thetas > (PI - delta/10.0)
    thetas[pi_inds] -= delta

    # Calculate the normalized associated Legendre function values,
    # and the values of the derivative with respect to cos(theta)
    lpmn, dlpmn = lpmn_norm(n_max, m_max, thetas)  # n_max+1 by m_max+1 by len(thetas)

    # Extend the Lpmn_norm array to negative values of m,
    # such that Lpmn_norm of -m is equal to Lpmn_norm of m.
    # Also remove the n == 0 part of the array.
    lpmn_norm_extended = np.zeros((n_max, 2*m_max+1, len(thetas)))
    lpmn_norm_extended[:, 0:m_max, :] = np.fliplr(lpmn[1:, 1:, :])
    lpmn_norm_extended[:, m_max:, :] = lpmn[1:, :, :]  # n_max by 2*m_max+1 by len(thetas)
    lpmn_norm_extended = lpmn_norm_extended.reshape((n_max, 2*m_max+1, len(thetas)))

    # Extend the dLpmn_norm array to negative values of m,
    # such that dLpmn_norm of -m is equal to dLpmn_norm of m.
    # Also remove the n == 0 part of the array.
    dlpmn_norm_extended = np.zeros((n_max, 2*m_max+1, len(thetas)))
    dlpmn_norm_extended[:, 0:m_max, :] = np.fliplr(dlpmn[1:, 1:, :])
    dlpmn_norm_extended[:, m_max:, :] = dlpmn[1:, :, :]  # n_max by 2*m_max+1 by len(thetas)
    dlpmn_norm_extended = dlpmn_norm_extended.reshape((n_max, 2*m_max+1, len(thetas)))

    # Create arrays that hold the n and m indices. Also promote
    # the thetas and phis to 4-dimensional arrays. This will come in handy
    # when later multiplying the matrices (due to numpy's "broadcasting" feature).
    n_vec = np.linspace(1, n_max, n_max)
    n_array = np.reshape(n_vec, (n_max, 1, 1))
    m_vec = np.linspace(-m_max, m_max, 2*m_max+1)
    m_array = np.reshape(m_vec, (1, 2*m_max+1, 1))
    thetas = np.reshape(thetas, (1, 1, len(thetas)))

    # Calculate the (-m/m)**m factor, using [1],(2.19) to handle
    # the case when m==0
    nonzeroinds = m_array != 0
    mm = np.ones(np.shape(m_array))
    mm[nonzeroinds] = (-m_array[nonzeroinds] /
                       abs(m_array[nonzeroinds]))**m_array[nonzeroinds]

    # Initialize the incomplete farfield pattern function arrays,
    # one for theta polarization and one for phi polarization
    k_n_m_s_theta = np.zeros((n_max, 2*m_max+1, np.size(thetas), 2), dtype='complex128')
    k_n_m_s_phi = np.zeros((n_max, 2*m_max+1, np.size(thetas), 2), dtype='complex128')

    # Calculate some common terms
    term1 = np.sqrt(2.0/(n_array*(n_array+1.0))) * mm
    term2a = (-1j)**(n_array+1.0)
    term2b = (-1j)**n_array
    term3a = 1j*m_array*lpmn_norm_extended/np.sin(thetas)
    term3b = dlpmn_norm_extended
    term_1_2a = term1*term2a
    term_1_2b = term1*term2b

    # Calculate the incomplete farfield pattern function for theta-pol when s = 1, according to [1],(A1.59)
    # and a FUTURE BLOG POST
    k_n_m_s_theta[:, :, :, 0] = term_1_2a * term3a

    # Calculate the incomplete farfield pattern function for theta-pol when s = 2, according to [1],(A1.60)
    # and a FUTURE BLOG POST
    k_n_m_s_theta[:, :, :, 1] = term_1_2b * term3b

    # Calculate the incomplete farfield pattern function for phi-pol when s = 1, according to [1],(A1.59)
    # and a FUTURE BLOG POST
    k_n_m_s_phi[:, :, :, 0] = term_1_2a * -term3b

    # Calculate the incomplete farfield pattern function for phi-pol when s = 2, according to [1],(A1.60)
    # and a FUTURE BLOG POST
    k_n_m_s_phi[:, :, :, 1] = term_1_2b * term3a

    # Swap axes
    k_n_m_s_theta = np.swapaxes(k_n_m_s_theta, 2, 3)
    k_n_m_s_phi = np.swapaxes(k_n_m_s_phi, 2, 3)

    return k_n_m_s_theta, k_n_m_s_phi

# ------------------------------------------------------------------------


def delta_pyramid(n_max):
    # Function used to calculate the delta pyramid.
    # The delta pyramid is defined in [1], Section A2.4

    # Region Diagram
    #
    #      |<--(-m)---m---(+m)-->|
    #  -    ---------------------
    #  ˄   |          |          |
    #  |   |          |          |
    # -m'  |    IV    |   III    |
    #  |   |          |          |
    #  |   |          |          |
    #  m'  |---------------------|
    #  |   |          |          |
    #  |   |          |          |
    # +m'  |    II    |    I     |
    #  |   |          |          |
    #  ˅   |          |          |
    #  -    ---------------------

    # Set mp_max and m_max equal to n_max in order to simplify calculations
    mp_max = m_max = n_max

    # Initialize the matrix to store the entire delta pyramid
    deltas = np.zeros((n_max+1, 2*mp_max+1, 2*m_max+1))

    # Initialize two "helper" functions
    def get_mp_index(mp_): return mp_max + mp_

    def get_m_index(m_): return m_max + m_

    # Initialize the delta values for n==0 and n==1, [1],Section A2.6
    deltas[0, get_mp_index(0), get_m_index(0)] = 1.0
    deltas[1, get_mp_index(0), get_m_index(1)] = -np.sqrt(2.0)/2.0
    deltas[1, get_mp_index(1), get_m_index(0)] = np.sqrt(2.0)/2.0
    deltas[1, get_mp_index(1), get_m_index(1)] = 1.0/2.0

    # Initialize positive values of mp and m
    mp_vec = np.linspace(0, mp_max, mp_max+1)
    m_vec = np.linspace(0, m_max, m_max+1)

    # Promote mp_vec and m_vec to 2 dimensions
    # (this will facilitate broadcasting later)
    mp_array = np.reshape(mp_vec, (mp_max+1, 1))
    m_array = np.reshape(m_vec, (1, m_max+1))

    # Copy two common indices to their own variables
    ind_mp_0 = get_mp_index(0)
    ind_m_0 = get_m_index(0)

    # Traverse through each n level from n==2 to n==N
    for n in range(2, n_max+1):

        # Copy two more common indices to their own variables
        ind_mp_n = get_mp_index(n)
        ind_m_n = get_m_index(n)

        # Calculate this particular n level using the previous two n levels.
        # [1],(A2.35)
        term1 = -1.0/(np.sqrt((n+mp_array[0:n, :]) *
                                 (n-mp_array[0:n, :]) *
                                 (n+m_array[:, 0:n]) *
                                 (n-m_array[:, 0:n]))*(n-1))
        term2 = np.sqrt((n+mp_array[0:n, :]-1) *
                           (n-mp_array[0:n, :]-1) *
                           (n+m_array[:, 0:n]-1) *
                           (n-m_array[:, 0:n]-1))*n
        term3 = (2*n-1)*mp_array[0:n, :]*m_array[:, 0:n]

        deltas[n, ind_mp_0:ind_mp_n, ind_m_0:ind_mp_n] = term1*(
            term2*deltas[n-2, ind_mp_0:ind_mp_n, ind_m_0:ind_mp_n] +
            term3*deltas[n-1, ind_mp_0:ind_mp_n, ind_m_0:ind_mp_n])

        # Initialize the vector full of this iteration's n value
        n_vec = n*np.ones(n+1)
        # Calculate the bottom edge of Region I, [1],(A2.41)
        temp = 1/(2.0**n)*np.sqrt(sps.binom(2*n_vec, n_vec - m_vec[0:n+1]))
        deltas[n, ind_mp_n, ind_m_0:ind_m_n+1] = temp
        # Copy the bottom edge of Region I to the right edge of Region I
        # using the symmetry relationship [1],(A2.26) with m set equal to n.
        deltas[n, ind_mp_0:ind_mp_n, ind_m_n] = (-1)**(mp_vec[0:n] + n)*temp[0:n]

    # Create new n_vec, mp_vec, and m_vec vectors
    n_vec = np.linspace(0, n_max, n_max+1)
    mp_vec = np.linspace(-mp_max, mp_max, 2*mp_max+1)
    m_vec = np.linspace(-m_max, m_max, 2*m_max+1)

    # Promote n_vec, mp_vec, and m_vec to 3 dimensions
    # (this will facilitate broadcasting later)
    n_array = np.reshape(n_vec, (n_max+1, 1, 1))
    mp_array = np.reshape(mp_vec, (1, 2*mp_max+1, 1))
    m_array = np.reshape(m_vec, (1, 1, 2*m_max+1))

    # Compute Region II using symmetry relation [1],(A2.32)
    temp1 = ((-1)**(n_array + mp_array[:, get_mp_index(0):, :]) *
             deltas[:, get_mp_index(0):, get_m_index(1):])
    temp1 = temp1[:, :, ::-1]
    deltas[:, get_mp_index(0):, 0:get_m_index(0)] = temp1

    # Compute Region III using symmetry relation [1],(A2.28)
    temp2 = ((-1)**(n_array + m_array[:, :, get_m_index(0):]) *
             deltas[:, get_mp_index(1):, get_m_index(0):])
    temp2 = temp2[:, ::-1, :]
    deltas[:, 0:get_mp_index(0), get_m_index(0):] = temp2

    # Compute Region IV using symmetry relation [1],(A2.30)
    temp3 = deltas[:, get_mp_index(1):, get_m_index(1):]
    temp3 = temp3[:, ::-1, ::-1].transpose((0, 2, 1))
    deltas[:, 0:get_mp_index(0), 0:get_m_index(0)] = temp3

    # Return the deltas from the function
    return deltas


# ------------------------------------------------------------------------


def lpmn_norm(n_max, m_max, thetas):
    """Normalized associated Legendre functions per Hansen Eq. (A1.25).

    Returns
    -------
    legendre_m_n_thetas : ndarray of shape (n_max+1, m_max+1, len(thetas))
        Values of the normalized associated Legendre function bar{P}_n^|m|(cos theta).
    dlegendre_m_n_thetas : ndarray of same shape
        Derivative with respect to theta (not cos theta).

    Notes
    -----
    Uses scipy.special.assoc_legendre_p_all, the modern replacement for
    scipy.special.lpmn (deprecated in SciPy 1.15).  The assoc_legendre_p_all
    routine returns derivatives with respect to z = cos(theta); we apply
    the chain rule to convert to d/d(theta).

    Both lpmn and assoc_legendre_p_all use the Condon-Shortley convention
    (i.e., include a (-1)^m phase factor relative to Hansen).  We undo this
    by including the same (-1)^m factor in the normalization, matching the
    convention in Hansen [1], pg. 322.

    At theta = 0 or pi, dP/dz is singular for some (m, n).  Multiplying by
    -sin(theta) = 0 produces NaN entries in the derivative.  Callers of this
    function that depend on derivative values at the poles must handle the
    pole values separately (`rotation_coefficients` does so).
    """

    # Make sure that the thetas variable is a numpy array
    if isinstance(thetas, float):
        thetas = np.array([thetas])
    else:
        thetas = np.array(thetas)

    # Make sure that M <= N
    if m_max > n_max:
        raise ValueError('M must be less than or equal to N.')

    # Make sure that all theta values are between 0 and pi
    if np.any(thetas < 0) or np.any(thetas > PI):
        raise ValueError('theta values must be between 0 and pi, inclusive.')

    # Modern replacement for scipy.special.lpmn (deprecated in SciPy 1.15).
    # assoc_legendre_p_all is vectorized over the z input, so we can compute
    # all theta values in a single call (no Python theta loop needed).
    #
    # Output shape with diff_n=1: (2, n_max+1, 2*m_max+1, T) where
    #   axis 0 : derivative order (0 = function, 1 = first derivative w.r.t. z)
    #   axis 1 : degree n in [0, n_max]
    #   axis 2 : order m in [0, 1, ..., m_max, -m_max, -m_max+1, ..., -1] (FFT layout)
    #   axis 3 : theta sample index
    #
    # The order m runs over both signs (FFT layout); Hansen Eq. (A1.25) uses
    # only non-negative orders, so we slice axis 2 to indices [0, m_max].
    z = np.cos(thetas)
    p_and_deriv = sps.assoc_legendre_p_all(n_max, m_max, z, diff_n=1)
    p_values = p_and_deriv[0][:, :m_max+1, :]   # (n_max+1, m_max+1, T)
    p_dz     = p_and_deriv[1][:, :m_max+1, :]   # (n_max+1, m_max+1, T)

    # Convert derivative from d/dz to d/dtheta via chain rule:
    # d/dtheta = (dz/dtheta) * d/dz = -sin(theta) * d/dz.
    # At theta = 0 or pi, d/dz may be Inf and -sin(theta) = 0, producing NaN.
    # This is acceptable because callers handle the pole values separately.
    with np.errstate(invalid='ignore'):
        p_dtheta = p_dz * (-np.sin(thetas))   # broadcast (n, m, T) * (T,)

    # The (n_max+1, m_max+1, T) layout matches what downstream callers expect.
    legendre_m_n_thetas = p_values
    dlegendre_m_n_thetas = p_dtheta

    # Normalization per Hansen Eq. (A1.25):
    #   bar{P}_n^|m|(cos theta) = (-1)^m * sqrt((2n+1)/2 * (n-m)!/(n+m)!) * P_n^|m|(cos theta)
    # The (-1)^m factor undoes the Condon-Shortley phase included by scipy.
    normfactor = np.zeros((n_max+1, m_max+1, 1))
    for m in range(m_max+1):
        for n in range(n_max+1):
            if m > n:
                continue
            temp1 = math.factorial(n - m)
            temp2 = math.factorial(n + m)
            temp3 = fractions.Fraction(temp1, temp2)
            normfactor[n, m, 0] = (-1)**m * math.sqrt((2*n + 1)/2.0 * float(temp3))

    legendre_m_n_thetas *= normfactor
    dlegendre_m_n_thetas *= normfactor

    return legendre_m_n_thetas, dlegendre_m_n_thetas


# ------------------------------------------------------------------------


def rotation_coefficients(n_max, m_max, mu_max, thetas):
    # Function used to calculate the rotation coefficients.
    # The rotation coefficients are defined in [1],Section A2.3

    # Make sure that the thetas variable is a numpy array
    if isinstance(thetas, float):
        thetas = np.array([thetas])
    else:
        thetas = np.array(thetas)

    # Make sure that MU <= M <= N
    if not(mu_max <= m_max <= n_max):
        raise ValueError('MU must be less than or equal to M, '
                         'which in turn must be less than or '
                         'equal to N.')

    # Also make sure that MU >= 1
    if mu_max < 1:
        raise ValueError('MU must be greater than or equal to 1.')

    # Make sure that all theta values are between 0 and pi
    if np.any(thetas < 0) or np.any(thetas > PI):
        raise ValueError('theta values must be between 0 and pi, inclusive.')

    # For mu == -1, 0, or 1, calculate the rotation coefficients using [1],(A2.17),
    # (A2.18), and (A2.19). This is done to improve computation speed.
    ## print 'Calculating rotation coefficients for mu = -1, 0, and 1'

    # Calculate the normalized associated Legendre function values,
    # and the values of the derivative with respect to cos(theta)
    legendre_norm, dlegendre_norm = lpmn_norm(n_max, m_max, thetas)

    # Extend the legendre_norm array to negative values
    # of m (to ease future calculations), but set the arrays
    # such that legendre_norm of -m is equal to legendre_norm of m.
    # Also remove the n == 0 part of the array.
    lpmn_norm_extended = np.zeros((n_max, 2*m_max+1, len(thetas)))
    lpmn_norm_extended[:, 0:m_max, :] = np.fliplr(legendre_norm[1:, 1:, :])
    lpmn_norm_extended[:, m_max:, :] = legendre_norm[1:, :, :]

    # Extend the dlegendre_norm array to negative values
    # of m (to ease future calculations), but set the arrays
    # such that legendre_norm of -m is equal to legendre_norm of m.
    # Also remove the n == 0 part of the array.
    dlpmn_norm_extended = np.zeros((n_max, 2*m_max+1, len(thetas)))
    dlpmn_norm_extended[:, 0:m_max, :] = np.fliplr(dlegendre_norm[1:, 1:, :])
    dlpmn_norm_extended[:, m_max:, :] = dlegendre_norm[1:, :, :]

    # Initialize the matrix that will hold the rotation coefficients
    d = np.zeros((n_max, 3, 2*m_max+1, len(thetas)))

    # Create arrays that hold the n, m, and mu indices. Also promote
    # the thetas to a 4-dimensional array. This will come in handy
    # when later multiplying the matrices. No need to perform
    # a matrix replication because of the way that Numpy "broadcasts"
    # the arrays during multiplication.
    n_vec = np.linspace(1, n_max, n_max)
    n_array = np.reshape(n_vec, (n_max, 1, 1, 1))
    mu_vec = np.linspace(-1, 1, 3)
    mu_array = np.reshape(mu_vec, (1, 3, 1, 1))
    m_vec = np.linspace(-m_max, m_max, 2*m_max+1)
    m_array = np.reshape(m_vec, (1, 1, 2*m_max+1, 1))

    # Calculate the (-m/m)**m factor, using [1],(2.19) to handle
    # the case when m==0
    nonzeroinds = m_array != 0
    mm = np.ones(np.shape(m_array))
    mm[nonzeroinds] = (-m_array[nonzeroinds] /
                       abs(m_array[nonzeroinds]))**m_array[nonzeroinds]

    # Calculate leading terms that are used in [1],(A2.17),(A2.18),
    # and (A2.19)
    c1 = mm*np.sqrt(2.0/(2.0*n_array+1))
    c2 = -2.0/np.sqrt(n_array*(n_array+1.0))

    # Calculate the mu==0 rotation coefficients according to (A2.17)
    d_0 = c1[:, 0, :, :]*lpmn_norm_extended

    # Initialize the d_plus1 and d_minus1 arrays
    d_plus1 = np.zeros([n_max, 2*m_max+1, len(thetas)])
    d_minus1 = np.zeros([n_max, 2*m_max+1, len(thetas)])

    # Calculate the mu==-1 and mu==1 rotation coefficients by
    # solving for d_-1 and d_+1 using equations (A2.18) and (A2.19).
    # Only calculate for values of theta where 1e-6 <= theta <= pi-1e-6.
    inds = np.logical_and((thetas >= 1e-6), (thetas <= PI-1e-6))
    d_plus1__plus__d_minus1 = (c2[:, 0, :, :]*m_array[:, 0, :, :])*d_0[:, :, inds]/np.sin(thetas[inds])
    d_plus1__minus__d_minus1 = (c1[:, 0, :, :]*c2[:, 0, :, :])*dlpmn_norm_extended[:, :, inds]
    d_minus1[:, :, inds] = (d_plus1__plus__d_minus1 - d_plus1__minus__d_minus1)/2.0
    d_plus1[:, :, inds] = (d_plus1__plus__d_minus1 + d_plus1__minus__d_minus1)/2.0

    # Assign d_minus1, d_0, and d_plus1 back into the d array
    d[:, 0, :, :] = d_minus1
    d[:, 1, :, :] = d_0
    d[:, 2, :, :] = d_plus1

    # Calculate the special cases when theta < 1e-6 or theta > pi-1e-6
    inds = thetas < 1e-6
    d[:, :, :, inds] = mu_array == m_array
    inds = thetas > PI-1e-6
    d[:, :, :, inds] = (-1)**(n_array+m_array)*(mu_array == -m_array)

    # If MU > 1, calculate the mu != -1, 0, or 1 rotation coefficients using [1],(A2.11)
    if mu_max > 1:

        # Store d for mu == -1,0,1 to a temp variable
        d_temp = d

        # Output shaped like this: (n,mu,m,theta)
        d = np.zeros((n_max, 2*mu_max+1, 2*m_max+1, len(thetas)))

        # Define a helper function
        def get_mu_index(mu_): return mu_max + mu_

        # Place d_temp back inside of d
        d[:, get_mu_index(-1):get_mu_index(1)+1, :, :] = d_temp

        # Calculate the delta pyramid
        deltas = delta_pyramid(n_max)

        # Create the mp and m vectors
        mp_vec = np.linspace(-n_max, n_max, 2*n_max+1)
        m_vec = np.linspace(-m_max, m_max, 2*m_max+1)

        # Remove the extra m values from the delta pyramid, if necessary
        extra = (len(mp_vec)-len(m_vec))//2
        if extra > 0:
            deltas = deltas[:, :, extra:-extra]

        # promote mp_vec to 4 dimensions
        mp_array = np.reshape(mp_vec, (1, 2*n_max+1, 1, 1))

        # Assign len(thetas) to a variable before we promote it to 4 dimensions
        numthetas = len(thetas)

        # promote thetas to 4 dimensions
        thetas = np.reshape(thetas, (1, 1, 1, len(thetas)))

        # Delete n==0 in deltas
        deltas = np.delete(deltas, 0, 0)

        # promote the deltas to 4 dimensions
        deltas4d = np.reshape(deltas, (n_max, 2*n_max+1, 2*m_max+1, 1))

        # Calculate the exponent matrix
        expon = np.exp(-1j*mp_array*thetas)

        # promote the m_vec to 2 dimensions
        m_array = np.reshape(m_vec, (1, 2*m_max+1))

        # define a helper function
        def get_m_index(m_): return m_max + m_

        # Vectorized summation for mu != -1, 0, 1.
        # Implements Hansen Eq. (A2.11):
        #   d^n_{mu, m}(theta) = Re[ i^(mu - m) *
        #                            sum_{m'} Delta^n_{m', mu} * Delta^n_{m', m} *
        #                            exp(-i * m' * theta) ]
        mu_vec = np.linspace(-mu_max, mu_max, 2*mu_max+1)
        mu_vec = mu_vec[np.logical_or(mu_vec < -1, mu_vec > 1)]
        if mu_vec.size > 0:
            # Indices in the deltas array's m-axis for each mu in mu_vec.
            mu_indices_in_deltas = (m_max + mu_vec).astype(int)        # (K_mu,)
            # Indices in the output array's mu-axis for each mu in mu_vec.
            mu_indices_in_d = (mu_max + mu_vec).astype(int)             # (K_mu,)

            # Delta^n_{m', mu} for each mu in mu_vec.  Shape (n_max, 2*n_max+1, K_mu).
            delta_mu = deltas[:, :, mu_indices_in_deltas]

            # Squeeze the singleton axes from the previously broadcasted
            # exponent so we can pass it to einsum as (m', theta).
            expon_2d = expon[0, :, 0, :]                                 # (2*n_max+1, num_thetas)

            # Sum over m' via einsum:
            #   result[n, K, m, t] = sum_i delta_mu[n, i, K] * deltas[n, i, m] * expon_2d[i, t]
            result = np.einsum('niK,nim,it->nKmt',
                                  delta_mu, deltas, expon_2d)

            # Multiply by i^(mu - m); shape (K_mu, 2*m_max+1).
            mu_phase = (1j ** (mu_vec[:, None] - m_vec[None, :]))
            result = result * mu_phase[None, :, :, None]

            # Place the real part into d at the correct mu indices.
            d[:, mu_indices_in_d, :, :] = result.real

    return d

# ------------------------------------------------------------------------


def translation_coefficients(n_max, nu_max, ka):
    # Vectorized computation of the translation coefficients
    # C^{sn(3)}_{sigma, mu, nu}(kA) per Hansen, Spherical Near-Field
    # Antenna Measurements, Eq. (3.19) / Appendix A3.
    #
    # Returns an array of shape (s, n, sigma, nu, mu) with
    #   axis 0 (s)     : s in {1, 2}   (s_idx = s - 1)
    #   axis 1 (n)     : n in {1, ..., n_max}
    #   axis 2 (sigma) : sigma in {1, 2}
    #   axis 3 (nu)    : nu in {1, ..., nu_max}
    #   axis 4 (mu)    : mu in {-1, +1} (mu_idx = 0 for -1, 1 for +1)
    #
    # Notes on the (s, sigma) block structure
    # ---------------------------------------
    # Only the sigma = 1 columns are computed directly via Eq. (3.19).
    # The sigma = 2 columns are then filled in via the elementary block
    # symmetries (Hansen Eqs. A3.8, A3.9):
    #     C^{2n(3)}_{2, mu, nu} = C^{1n(3)}_{1, mu, nu}     (Eq. A3.8)
    #     C^{1n(3)}_{2, mu, nu} = C^{2n(3)}_{1, mu, nu}     (Eq. A3.9)
    # Similarly, only mu = +1 is computed directly; mu = -1 is obtained
    # by multiplying by the (-1)^(s + sigma) factor from Eq. (A3.21).

    p_total_max = n_max + nu_max

    # B'(J) array per Hansen Eq. (A3.27): B'(0) = 1, B'(J+2) = (J+1)/(J+2) * B'(J)
    b_prime = np.zeros(2 * p_total_max + 1)
    b_prime[0] = 1.0
    for J in range(0, 2 * p_total_max, 2):
        b_prime[J + 2] = (J + 1.0) / (J + 2.0) * b_prime[J]

    # Spherical Hankel function of the first kind for orders 0..p_total_max
    hn = sph_hankel_first_kind(p_total_max, ka)

    # n and nu grids
    n_int = np.arange(1, n_max + 1)[:, None]            # (n_max, 1)
    nu_int = np.arange(1, nu_max + 1)[None, :]          # (1, nu_max)
    n_g = n_int.astype(float)
    nu_g = nu_int.astype(float)

    # Front terms (n_max, nu_max) per Hansen Eq. (3.19) prefactor
    front_terms = (
        0.25
        * (1j ** n_g) * (np.sqrt(2.0*n_g + 1.0) / (n_g * (n_g + 1.0)))
        * (1j ** (-nu_g)) * (np.sqrt(2.0*nu_g + 1.0) / (nu_g * (nu_g + 1.0)))
    )

    # Accumulators for the p-sum.
    # diag_sum picks up the s == sigma contribution; off_diag_sum picks up
    # the s != sigma contribution.
    diag_sum = np.zeros((n_max, nu_max), dtype='complex128')
    off_diag_sum = np.zeros((n_max, nu_max), dtype='complex128')

    # Loop only over p; everything else is vectorized over (n, nu).
    for p in range(p_total_max + 1):
        # Validity mask per Hansen Eq. (A3.7) parity rule and triangle inequality.
        valid = (np.abs(n_int - nu_int) <= p) & (p <= n_int + nu_int) \
                & (((n_int + nu_int + p) % 2) == 0)
        if not valid.any():
            continue

        d1 = n_g * (n_g + 1.0) + nu_g * (nu_g + 1.0) - p * (p + 1.0)
        d2 = n_g + nu_g + p + 1.0

        # B'(J) indices; clip negatives to 0 (those entries will be masked out).
        idx_a = np.clip(-n_int + nu_int + p, 0, None)
        idx_b = np.clip( n_int - nu_int + p, 0, None)
        idx_c = np.clip( n_int + nu_int - p, 0, None)
        idx_d = n_int + nu_int + p   # always non-negative

        term3 = b_prime[idx_a] * b_prime[idx_b] * b_prime[idx_c] / b_prime[idx_d]
        term4 = (1j ** (-p)) * (2.0 * p + 1.0)

        contrib = term3 * term4 * hn[p]                # (n_max, nu_max)

        diag_sum += np.where(valid, (d1 * d1 / d2) * contrib, 0.0)
        off_diag_sum += np.where(valid, (2.0 * 1j * ka * d1 / d2) * contrib, 0.0)

    # Assemble C[s, n, sigma, nu, mu], shape (2, n_max, 2, nu_max, 2).
    c_s_n_sig_nu_mu = np.zeros((2, n_max, 2, nu_max, 2), dtype='complex128')

    # sigma = 1, mu = +1 (coeff = 1)
    c_s_n_sig_nu_mu[0, :, 0, :, 1] = front_terms * diag_sum       # s=1, sigma=1 (diag)
    c_s_n_sig_nu_mu[1, :, 0, :, 1] = front_terms * off_diag_sum   # s=2, sigma=1 (off-diag)

    # sigma = 1, mu = -1 (coeff = (-1)^(s + sigma))
    c_s_n_sig_nu_mu[0, :, 0, :, 0] = front_terms * diag_sum * ((-1.0) ** (1 + 1))   # +1
    c_s_n_sig_nu_mu[1, :, 0, :, 0] = front_terms * off_diag_sum * ((-1.0) ** (2 + 1))  # -1

    # Apply elementary block symmetries (Hansen Eqs. A3.8, A3.9):
    #   C^{2n(3)}_{2,mu,nu} = C^{1n(3)}_{1,mu,nu}
    #   C^{1n(3)}_{2,mu,nu} = C^{2n(3)}_{1,mu,nu}
    c_s_n_sig_nu_mu[1, :, 1, :, :] = c_s_n_sig_nu_mu[0, :, 0, :, :]
    c_s_n_sig_nu_mu[0, :, 1, :, :] = c_s_n_sig_nu_mu[1, :, 0, :, :]

    return c_s_n_sig_nu_mu

# ------------------------------------------------------------------------


def sph_hankel_first_kind(n, x):

    # Returns h_n^{(1)}(x) = j_n(x) + i y_n(x) for orders 0, 1, ..., n
    # as a 1D array of length n+1, matching the legacy sph_jn/sph_yn API
    # used previously. (scipy.special.sph_jn / sph_yn were removed in
    # scipy 1.0; we use the vectorized spherical_jn / spherical_yn instead.)

    orders = np.arange(n + 1)

    jn = sps.spherical_jn(orders, x)
    yn = sps.spherical_yn(orders, x)

    return jn + 1j*yn

# ------------------------------------------------------------------------


def _fft(c, n=None, axis=-1):
    """Hansen-convention DFT (positive exponent, no normalization).

    Equivalent to ``n * np.fft.ifft(c, n, axis)``.  See the module
    docstring (DFT / IDFT convention) for details.  Marked with a leading
    underscore to make clear that this is an internal helper, *not* the
    standard `np.fft.fft`.
    """
    if n is None:
        n = np.size(c, axis)

    c_k = n*np.fft.ifft(c, n, axis)

    return c_k


# ------------------------------------------------------------------------


def _ifft(c, n=None, axis=-1):
    """Hansen-convention IDFT (negative exponent, with 1/n normalization).

    Equivalent to ``np.fft.fft(c, n, axis) / n``.  See the module
    docstring (DFT / IDFT convention) for details.  Marked with a
    leading underscore as an internal helper.
    """
    if n is None:
        n = np.size(c, axis)

    c_m = (1.0/n)*np.fft.fft(c, n, axis)

    return c_m


# ------------------------------------------------------------------------


# def PI(l, mp):
#     if (l-mp) % 2 == 1:
#         return 0.0
#     elif (l-mp) % 2 == 0:
#         return 2.0/(1.0 - (l - mp)**2.0)
#     else:
#         raise Exception("error in pysnf.PI")


# ------------------------------------------------------------------------


def pi_wiggle(n_max):

    jj = np.linspace(-2*n_max+1, 2*n_max, 4*n_max)

    even_indices = jj % 2 == 0

    pi_wig = np.zeros(4*n_max)

    pi_wig[even_indices] = 2.0/(1.0 - jj[even_indices]**2)

    pi_wig = np.roll(pi_wig, -2*n_max+1)

    pi_wig = np.reshape(pi_wig,(4*n_max,1,1))

    return pi_wig

# ------------------------------------------------------------------------


def b_wiggle(b_l_m_mu):

    numrows, numcolumns, numpages = np.shape(b_l_m_mu)
    n_max = (numrows-1)//2
    m_max = (numcolumns-1)//2

    # Inherit dtype from the input array so dtype-propagation works end-to-end.
    b_l_m_mu_wiggle = np.zeros((4*n_max, 2*m_max+1, 2), dtype=b_l_m_mu.dtype)

    def get_l_index(l_): return l_ + 2*n_max - 1

    b_l_m_mu_wiggle[get_l_index(-n_max):get_l_index(n_max)+1, :, :] = b_l_m_mu

    b_l_m_mu_wiggle = np.roll(b_l_m_mu_wiggle, shift=-2*n_max+1, axis=0)

    return b_l_m_mu_wiggle

# ------------------------------------------------------------------------


def db(a):
    return 20*np.log10(np.absolute(a))

# ------------------------------------------------------------------------


def plot_spherical_wave_coefficients_mag_db(q_in):

        # Lazy import - matplotlib is only needed for plotting.
        import matplotlib.pyplot as plt

        q = q_in.copy()

        qmax = np.max((np.absolute(q[:,:,0])**2 +
                          np.absolute(q[:,:,1])**2).ravel())
        q = q/qmax

        fig1, (ax1, ax2) = plt.subplots(nrows=2)

        m_max = (np.size(q, 1)-1)//2
        n_max = np.size(q, 0)

        ax1.imshow(db(q[:, :, 0]).T,
                   extent=[1, n_max, -m_max, m_max], aspect='auto', interpolation='none', vmin=-120, vmax=0)
        ax2.imshow(db(q[:, :, 1]).T,
                   extent=[1, n_max, -m_max, m_max], aspect='auto', interpolation='none', vmin=-120, vmax=0)

        fig2, ax3 = plt.subplots(nrows=1)
        ax3.imshow((10*np.log10(np.absolute(q[:,:,0])**2 + np.absolute(q[:,:,1])**2)).T,
                   extent=[1, n_max, -m_max, m_max], aspect='auto', interpolation='none', vmin=-120, vmax=0)

        plt.show()

# ------------------------------------------------------------------------


def dipole_probe_response_constants(n_max, ka=INF, direction='+x', dipole_type='electric'):
    """Return the probe response constants for an ideal Hertzian dipole.

    Currently only the (``direction='+x'``, ``dipole_type='electric'``,
    ``ka=INF``) configuration is implemented; any other combination
    raises ``NotImplementedError``.
    """

    # Validate inputs.  Only one (direction, dipole_type, ka) combination
    # is implemented today; raise immediately for everything else so the
    # caller doesn't get a silent zero array.  The detailed error tells
    # the caller exactly which parameter(s) are out of range.
    unsupported = []
    if direction != '+x':
        unsupported.append("direction={!r} (only '+x' is implemented)".format(direction))
    if dipole_type != 'electric':
        unsupported.append("dipole_type={!r} (only 'electric' is implemented)".format(dipole_type))
    if ka != INF:
        unsupported.append("ka={!r} (only INF is implemented; use INF for the "
                           "far-field limit)".format(ka))
    if unsupported:
        raise NotImplementedError(
            "dipole_probe_response_constants: unsupported parameter(s): "
            + "; ".join(unsupported))

    # Far-field, +x-directed, electric Hertzian dipole probe.
    # See Hansen, Spherical Near-Field Antenna Measurements, Eqs. (4.103)-(4.104).
    p_n_mu_s = np.zeros((n_max, 2, 2), dtype='complex128')
    n_array = np.linspace(1, n_max, n_max)
    temp = -np.sqrt(6.0)/8.0 * np.sqrt(2.0*n_array + 1.0) * 1j**(-n_array)
    p_n_mu_s[:, 1, 0] = temp     # mu=+1, s=1
    p_n_mu_s[:, 1, 1] = temp     # mu=+1, s=2
    p_n_mu_s[:, 0, 0] = temp     # mu=-1, s=1
    p_n_mu_s[:, 0, 1] = -temp    # mu=-1, s=2

    return p_n_mu_s


# ------------------------------------------------------------------------


def reciprocity(trans_coeffs):

    # [1],(2.104)

    # T is nu, mu, sig

    n_max, m_max, s_max = np.shape(trans_coeffs)

    if m_max > 2:
        m_max = (m_max-1)//2

        receive_coeffs = np.zeros((n_max, 2*m_max+1, 2), dtype='complex128')

        for m in range(-m_max, m_max+1):

            mi = m + m_max
            neg_mi = m_max - m

            receive_coeffs[:, mi, :] = (-1)**m * trans_coeffs[:, neg_mi, :]

    else:

        receive_coeffs = np.zeros((n_max, 2, 2), dtype='complex128')

        receive_coeffs[:, 1, :] = -1*trans_coeffs[:, 0, :]
        receive_coeffs[:, 0, :] = -1*trans_coeffs[:, 1, :]

    return receive_coeffs

# ------------------------------------------------------------------------


def rotate_wavecoeffs_about_axis(q, ax):

    # [1],(5.67)

    # Q is n, m, s

    n_max, m_max, s_max = np.shape(q)

    n_vec = np.linspace(1.0, n_max, n_max)
    n_vec = np.reshape(n_vec, (n_max, 1))

    if ax == 'x':

        if m_max > 2:

            m_max = (m_max-1)//2

            q_rot = np.zeros((n_max, 2*m_max+1, 2), dtype='complex128')

            for m in range(-m_max, m_max+1):

                mi = m + m_max
                neg_mi = m_max - m

                q_rot[:, mi, :] = (-1)**n_vec * q[:, neg_mi, :]

        else:

            q_rot = np.zeros((n_max, 2, 2), dtype='complex128')

            q_rot[:, 1, :] = (-1)**n_vec * q[:, 0, :]
            q_rot[:, 0, :] = (-1)**n_vec * q[:, 1, :]

    else:

        raise NotImplementedError("Code not yet implemented for rotations about y or z axes.")

    return q_rot

# ------------------------------------------------------------------------


def _translation_coefficients_far_field(n_max, nu_max):
    # Asymptotic translation coefficients per Hansen Eqs. (A3.23) / (A3.24),
    # with the common e^{ikA} / (kA) factor *stripped*.
    #
    # Returns an array of shape (s, n, sigma, nu, mu) matching
    # `translation_coefficients`, where:
    #   axis 0 (s)     : s in {1, 2}   (s_idx = s - 1)
    #   axis 1 (n)     : n in {1, ..., n_max}
    #   axis 2 (sigma) : sigma in {1, 2}
    #   axis 3 (nu)    : nu in {1, ..., nu_max}
    #   axis 4 (mu)    : mu in {-1, +1} (mu_idx = 0 for -1, 1 for +1)
    n_arr = np.arange(1, n_max + 1, dtype=float)
    nu_arr = np.arange(1, nu_max + 1, dtype=float)

    # Base value depends only on (n, nu): shape (n_max, nu_max)
    base = (np.sqrt(np.outer(2.0*n_arr + 1.0, 2.0*nu_arr + 1.0)) / 2.0 *
            (1j ** (nu_arr[None, :] - n_arr[:, None] - 1.0)))

    c_far = np.zeros((2, n_max, 2, nu_max, 2), dtype='complex128')

    # mu = +1: all (s, sigma) entries are equal to `base`
    c_far[:, :, :, :, 1] = base[None, :, None, :]

    # mu = -1: multiply by (-1)^(s + sigma)
    s_sig_signs = np.array([[+1.0, -1.0],
                               [-1.0, +1.0]])  # [s_idx, sig_idx]
    c_far[:, :, :, :, 0] = base[None, :, None, :] * s_sig_signs[:, None, :, None]

    return c_far


def probe_response_constants(r_p, n_max, ka, c=None):

    # Assumes that R_p is entered as sig, mu, nu, with the mu axis having an
    # *odd* number of entries centered on mu=0 (i.e. entries are mu = -mu_max,
    # ..., 0, ..., +mu_max for a total length 2*mu_max + 1).
    #
    # When ka == INF, the returned probe response constants are the normalized
    # far-field values P_{s, mu, n}^infinity (per Hansen Eqs. (4.100)-(4.102)),
    # i.e. with the common e^{ikA} / (kA) factor stripped.

    if (ka != INF) and (ka > 1e4):
        raise NotImplementedError(
            "probe_response_constants: finite kA > 1e4 is not yet supported. "
            "Use ka = INF (np.inf) for the far-field limit, or call with "
            "ka <= 1e4 for a finite-range probe.")

    sig_count, mu_count, nu_max = np.shape(r_p)

    if sig_count != 2:
        raise ValueError("First dimension (sigma) of receiving coefficients must have a size of exactly 2.")

    if mu_count < 2:
        raise ValueError("Second dimension (mu) of receiving coefficients must have a size of at least 2.")
    elif mu_count > 2:
        # mu axis must be odd-length and centered on mu=0
        if mu_count % 2 == 0:
            raise ValueError("mu axis must have odd length (centered on mu=0).")
        mu_center_idx = (mu_count - 1) // 2     # index of mu = 0
        mu_neg1_idx = mu_center_idx - 1          # index of mu = -1
        mu_pos1_idx = mu_center_idx + 1          # index of mu = +1
        temp = r_p.copy()
        r_p = np.zeros((2, 2, nu_max), dtype='complex128')
        r_p[:, 0, :] = temp[:, mu_neg1_idx, :]
        r_p[:, 1, :] = temp[:, mu_pos1_idx, :]

    if c is None:
        if ka == INF:
            c = _translation_coefficients_far_field(n_max, nu_max)  # s, n, sig, nu, mu
        else:
            c = translation_coefficients(n_max, nu_max, ka)         # s, n, sig, nu, mu

    r_p = np.swapaxes(r_p, 1, 2)  # R_p is now sig, nu, mu

    r_p = np.reshape(r_p, (1, 1, 2, nu_max, 2))  # R_p is now s, n, sig, nu, mu

    temp = np.sum(c*r_p, axis=3)  # Sum over nu. temp is s, n, sig, mu

    p = 0.5*np.sum(temp, axis=2)  # Sum over sigma. P is s, n, mu

    p = np.swapaxes(p, 1, 2)  # P is now s, mu, n

    p = np.swapaxes(p, 0, 2)  # P is now n, mu, s

    return p

# ------------------------------------------------------------------------


def singlesphere2doublesphere(singlesphere, dtype='complex128'):

    np_dtype = np.dtype(dtype)

    numthetas = np.size(singlesphere, axis=0)
    numphis = np.size(singlesphere, axis=1)

    if numphis % 2 == 1:
        ss = interpft(singlesphere, numphis-1)
        numphis -= 1
    else:
        ss = singlesphere.copy()

    doublesphere = np.zeros((2*(numthetas-1), numphis), dtype=np_dtype)

    doublesphere[:numthetas, :] = ss[:, :]
    doublesphere[numthetas:, :] = -np.roll(ss[-2:0:-1, :], numphis//2, axis=1)

    return doublesphere


# ------------------------------------------------------------------------


def interpft(a, ny):

    # Operates on the last axis of a 2D array
    axis = -1

    # Get initial length and width of the 2D array
    n, m = np.shape(a)

    # Ensure that ny is an integer
    ny = int(np.floor(ny))

    # If necessary, increase ny by an integer multiple to make ny > size(a,axis)
    if ny <= 0:
        raise ValueError("n must be an integer greater than 0.")
    elif np.size(a, axis) > m:
        incr = 1
    else:
        incr = int(m // ny) + 1
        ny *= incr

    b = _fft(a, axis=axis)

    nyqst = int(np.ceil((m + 1.0)/2.0))

    c = np.zeros((n, ny), dtype='complex128')
    c[:, 0:nyqst] = b[:, 0:nyqst]
    c[:, nyqst+(ny-m):] = b[:, nyqst:]

    if np.remainder(m, 2) == 0:
        c[:, nyqst-1] = c[:, nyqst-1]/2.0
        c[:, nyqst+ny-m-1] = c[:, nyqst-1]

    d = _ifft(c, axis=axis)

    d *= float(ny)/float(m)

    d = d[:, 0::incr]  # Skip over extra points when original ny <= m

    return d


# ---------
# SELF TEST
# ---------
if __name__ == '__main__':

    # -------------------------------------------------------------------
    # Round-trip regression test
    # -------------------------------------------------------------------
    # Construct a known T_{smn}, synthesize the far field via
    # wavecoeffs2farfield_uniform (using a dipole output probe), then
    # recover T_{smn} via field2wavecoeffs (using a dipole input probe).
    # The recovered coefficients should match the input to within
    # floating-point precision.
    #
    # This is the closest practical analog to the worked example in
    # Hansen, Spherical Near-Field Antenna Measurements, Section 4.4.5,
    # adapted to the dipole-probe assumption built into the current code.
    # -------------------------------------------------------------------

    n_max_test = 2
    m_max_test = 1

    # In the q_n_m_s array, indexing is q[n-1, m + m_max, s-1].
    # The four modes below correspond to an x-directed electric Hertzian
    # dipole at the origin (per Hansen Eqs. 2.124-2.125):
    #     T_{s=2, m=-1, n=1} = +sqrt(2)/2
    #     T_{s=2, m=+1, n=1} = -sqrt(2)/2
    # plus a z-directed magnetic dipole (Eqs. 2.132-2.133):
    #     T_{s=1, m= 0, n=1} = +1  (arbitrary non-zero choice)
    T_in = np.zeros((n_max_test, 2*m_max_test+1, 2), dtype='complex128')
    T_in[0, m_max_test + (-1), 2 - 1] = +np.sqrt(2.0)/2.0
    T_in[0, m_max_test + (+1), 2 - 1] = -np.sqrt(2.0)/2.0
    T_in[0, m_max_test + 0,    1 - 1] = 1.0

    # Synthesize far-field on a 5-degree grid using a dipole output probe
    dT_test = 5.0   # degrees
    dP_test = 5.0   # degrees
    theta_pol, phi_pol = wavecoeffs2farfield_uniform(T_in, dT_test, dP_test)

    # Package as input for analysis
    ff_in = {'E_theta': theta_pol, 'E_phi': phi_pol,
             'dT': dT_test, 'dP': dP_test}

    # Recover T_smn using a dipole input probe at infinity
    T_out = field2wavecoeffs(ff_in, probe=None, probe_pol='x', ka=INF,
                             n_max=n_max_test, m_max=m_max_test)

    # Compare
    max_abs_err = float(np.max(np.abs(T_out - T_in)))
    print('--- Round-trip regression test ---')
    print('  Input  T (non-zero entries):')
    for n_idx in range(n_max_test):
        for m_idx in range(2*m_max_test+1):
            for s_idx in range(2):
                v_in = T_in[n_idx, m_idx, s_idx]
                if abs(v_in) > 0:
                    print('    T[n={}, m={:+d}, s={}] = {}'.format(
                        n_idx+1, m_idx - m_max_test, s_idx+1, v_in))
    print('  Output T (non-zero entries, |T| > 1e-10):')
    for n_idx in range(n_max_test):
        for m_idx in range(2*m_max_test+1):
            for s_idx in range(2):
                v_out = T_out[n_idx, m_idx, s_idx]
                if abs(v_out) > 1e-10:
                    print('    T[n={}, m={:+d}, s={}] = {}'.format(
                        n_idx+1, m_idx - m_max_test, s_idx+1, v_out))
    print('  Max |T_out - T_in| = {:.3e}'.format(max_abs_err))
    if max_abs_err < 1e-10:
        print('  PASS')
    else:
        print('  FAIL (tolerance 1e-10)')

    # -------------------------------------------------------------------
    # Single-mode purity test
    # -------------------------------------------------------------------
    # Set only T_{s=2, m=0, n=1} = 1 (a pure z-directed electric dipole).
    # Synthesize, recover, and confirm that *every other* T_{smn} entry
    # is zero to within ~1e-10.  Mode-leakage bugs that a multi-mode test
    # might mask should show up here as a non-zero entry at a position
    # that started at zero.
    # -------------------------------------------------------------------

    T_pure = np.zeros((n_max_test, 2*m_max_test+1, 2), dtype='complex128')
    # Index for T_{s=2, m=0, n=1}: n_idx=0, m_idx=m_max_test, s_idx=1
    T_pure[0, m_max_test, 1] = 1.0

    theta_pol_p, phi_pol_p = wavecoeffs2farfield_uniform(T_pure, dT_test, dP_test)
    ff_pure = {'E_theta': theta_pol_p, 'E_phi': phi_pol_p,
               'dT': dT_test, 'dP': dP_test}
    T_pure_out = field2wavecoeffs(ff_pure, probe=None, probe_pol='x', ka=INF,
                                  n_max=n_max_test, m_max=m_max_test)

    # The recovered T should equal T_pure exactly (one non-zero entry of value 1)
    mode_err = float(np.max(np.abs(T_pure_out - T_pure)))
    # Mask out the one expected non-zero entry and check the leakage
    leakage = np.abs(T_pure_out).copy()
    leakage[0, m_max_test, 1] = 0.0
    max_leakage = float(np.max(leakage))
    recovered_value = T_pure_out[0, m_max_test, 1]
    print('--- Single-mode purity test (T_{s=2, m=0, n=1} = 1) ---')
    print('  Recovered T[s=2, m=0, n=1] = {}  (expected 1+0j)'.format(recovered_value))
    print('  Max leakage into other modes = {:.3e}'.format(max_leakage))
    print('  Max |T_out - T_in|            = {:.3e}'.format(mode_err))
    if (max_leakage < 1e-10) and (mode_err < 1e-10):
        print('  PASS')
    else:
        print('  FAIL (tolerance 1e-10)')
