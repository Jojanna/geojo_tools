from scipy.interpolate import interp1d
import numpy as np
import pandas as pd
import math
import sys

def ray_param(v, theta):
    '''
    Calculates the ray parameter p

    Usage:
    ------
        p = ray_param(v, theta)

    Inputs:
    -------
            v = interval velocity
        theta = incidence angle of ray (degrees)

    Output:
    -------
        p = ray parameter (i.e. sin(theta)/v )
    '''

    import math

    # Cast inputs to floats
    theta = float(theta)
    v = float(v)

    p = math.sin(math.radians(theta)) / v  # ray parameter calculation

    return p



def rc_zoep_2layer(vp1, vs1, rho1, vp2, vs2, rho2, theta1):
    '''
    Reflection & Transmission coefficients calculated using full Zoeppritz
    equations.

    Usage:
    ------
    R = rc_zoep(vp1, vs1, rho1, vp2, vs2, rho2, theta1)

    Reference:
    ----------
    The Rock Physics Handbook, Dvorkin et al.
    '''

    import math

    # Cast inputs to floats
    vp1 = float(vp1)
    vp2 = float(vp2)
    vs1 = float(vs1)
    vs2 = float(vs2)
    rho1 = float(rho1)
    rho2 = float(rho2)
    theta1 = float(theta1)

    # conversion to radians occurs within ray param function
    theta_deg = theta1
    theta_rad = math.radians(theta_deg)
    p = ray_param(vp1, theta_deg)
    # Calculate reflection & transmission angles
    # Convert theta1 to radians
    # Ray parameter
    if ((p * vp2) > 1.0) | ((p * vp2) < -1.0) | ((p * vs1) > 1.0) | ((p * vs1) < -1.0) |((p * vs2) > 1.0) | ((p * vs2) < -1.0):

        R = 0.0  # np.nan
    else:
        theta2 = math.asin(p * vp2);  # Transmission angle of P-wave
        phi1 = math.asin(p * vs1);  # Reflection angle of converted S-wave
        phi2 = math.asin(p * vs2);  # Transmission angle of converted S-wave

        # Matrix form of Zoeppritz Equations... M & N are two of the matricies
        M = np.array([ \
            [-math.sin(theta_rad), -math.cos(phi1), math.sin(theta2), math.cos(phi2)], \
            [math.cos(theta_rad), -math.sin(phi1), math.cos(theta2), -math.sin(phi2)], \
            [2 * rho1 * vs1 * math.sin(phi1) * math.cos(theta_rad), rho1 * vs1 * (1 - 2 * math.sin(phi1) ** 2), \
             2 * rho2 * vs2 * math.sin(phi2) * math.cos(theta2), rho2 * vs2 * (1 - 2 * math.sin(phi2) ** 2)], \
            [-rho1 * vp1 * (1 - 2 * math.sin(phi1) ** 2), rho1 * vs1 * math.sin(2 * phi1), \
             rho2 * vp2 * (1 - 2 * math.sin(phi2) ** 2), -rho2 * vs2 * math.sin(2 * phi2)]
        ], dtype='float')

        N = np.array([ \
            [math.sin(theta_rad), math.cos(phi1), -math.sin(theta2), -math.cos(phi2)], \
            [math.cos(theta_rad), -math.sin(phi1), math.cos(theta2), -math.sin(phi2)], \
            [2 * rho1 * vs1 * math.sin(phi1) * math.cos(theta_rad), rho1 * vs1 * (1 - 2 * math.sin(phi1) ** 2), \
             2 * rho2 * vs2 * math.sin(phi2) * math.cos(theta2), rho2 * vs2 * (1 - 2 * math.sin(phi2) ** 2)], \
            [rho1 * vp1 * (1 - 2 * math.sin(phi1) ** 2), -rho1 * vs1 * math.sin(2 * phi1), \
             -rho2 * vp2 * (1 - 2 * math.sin(phi2) ** 2), rho2 * vs2 * math.sin(2 * phi2)] \
            ], dtype='float')

        # This is the important step, calculating coefficients for all modes and rays
        if np.linalg.cond(M) < 1 / sys.float_info.epsilon:
            R = np.dot(np.linalg.inv(M), N)
        else:
            R = 0.0  # np.nan

    return R




def rc_zoep_trace(vp_log, vs_log, rhob_log, angles):
    test = False

    '''
    Reflection & Transmission coefficients calculated using full Zoeppritz
    equations.

    Usage:
    ------
    R = rc_zoep(vp1, vs1, rho1, vp2, vs2, rho2, theta1)

    Reference:
    ----------
    The Rock Physics Handbook, Dvorkin et al.
    '''

    import math

    # Cast inputs to floats

    vp = [float(x) for x in vp_log]
    vs = [float(x) for x in vs_log]
    rhob = [float(x) for x in rhob_log]
    angles = [float(x) for x in angles]

    angle_ix = list(range(0, len(angles)))
    samp = list(range(0, len(vp) - 1))

    refs = np.full(shape=(len(angles), len(vp)), fill_value=np.nan)

    for i, theta in zip(angle_ix, angles):

        theta_deg = theta
        theta_rad = math.radians(theta)
        # Convert theta1 to radians

        # modifications to try and handle critical rays
        r_log = [0.0]

        for j in samp:
            vp1 = vp[j]
            vp2 = vp[j + 1]
            vs1 = vs[j]
            vs2 = vs[j + 1]
            rho1 = rhob[j]
            rho2 = rhob[j + 1]

            ref_out = rc_zoep_2layer(vp1, vs1, rho1, vp2, vs2, rho2, theta_deg)
            if type(ref_out) == float:
                R = (0.0)
            else:
                R = ref_out[0, 0]
            r_log.append(R)

        # print (i, theta)
        refs[i, :] = r_log

    return refs


def convolve_volume(vol, nmodel, nangles, wvlt_amp, dt):
    # print (nmodel, nangle)
    sgy = np.ndarray(shape=(1, len(vol[0, :, 0, 0]), len(vol[0, 0, :, 0]), len(vol[0, 0, 0, :])))
    # print (np.shape(sgy))

    for i in range(0, nmodel):
        # print ("Model no. = %i" % i)

        for j in range(0, nangles):
            # print ("Angle no. = %i" % (j))
            rc = vol[0, i, j, :]
            syn = np.convolve(rc, wvlt_amp, mode='same')
            sgy[0, i, j, :] = syn

    return sgy