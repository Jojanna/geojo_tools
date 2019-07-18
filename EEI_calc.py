import numpy as np
import math
from sys import exit
import lasio
import pandas as pd


def eei_calc(md, vp, vs, rhob, zoi_upper, zoi_lower, eei_angle):

    null = -999.25
    a = math.radians(eei_angle)

    # estimate K from the mean value of(Vs / Vp) ^ 2 over ZOI
    vs_vp = list(s / p for s, p in zip(vs, vp))
    k_a = list(pow(n, 2) for n in vs_vp)

    nsamp_k = len(k_a)
    i = 0
    k_sum = 0
    k_samp = 0

    while (i < nsamp_k):
        z = md[i]
        if (z > zoi_upper and z < zoi_lower):
            k_sum = k_sum + k_a[i]
            k_samp = k_samp + 1
        else:
            pass
        i = i + 1
    if k_samp == 0:
        print ("check bounds for ZOI")
        exit(0)

    else:
        k = k_sum / k_samp
    #print(k)

    # calculate the normalised values for each log

    nsamp_vp = len(vp)
    i = 0
    vp_sum = 0
    vp_samp = 0
    while (i < nsamp_vp):

        z = md[i]
        if (z > zoi_upper and z < zoi_lower):
            vp_sum = vp_sum + vp[i]
            vp_samp = vp_samp + 1
        else:
            pass

        i = i + 1

    a0 = vp_sum / vp_samp
    #print(a0)

    nsamp_vs = len(vs)
    i = 0
    vs_sum = 0
    vs_samp = 0

    while (i < nsamp_vs):

        z = md[i]
        if (z > zoi_upper and z < zoi_lower):
            vs_sum = vs_sum + vs[i]
            vs_samp = vs_samp + 1
        else:
            pass
        i = i + 1

    b0 = vs_sum / vs_samp
    #print(b0)

    nsamp_rhob = len(rhob)
    rhob_sum = 0
    rhob_samp = 0

    i = 0
    while (i < nsamp_rhob):
        z = md[i]
        if (z > zoi_upper and z < zoi_lower):
            rhob_sum = rhob_sum + rhob[i]
            rhob_samp = rhob_samp + 1
        else:
            pass
        i = i + 1

    r0 = rhob_sum / rhob_samp
    #print(r0)
    # a0=mean(Vp)
    # b0=mean(Vs)
    # r0=mean(DEN)

    # define exponential functions     # calculate the EEI


    p = np.cos(a) + np.sin(a)
    q = (-8) * k * np.sin(a)
    r = np.cos(a) - 4 * k * np.sin(a)
    #print (p, q, r)
    eei = []
    for vp_n, vs_n, rho_n in zip(vp, vs, rhob):
        if rho_n < 0 or vp_n < 0 or vs_n < 0:
            eei.append(null)
        else:
            eei.append(a0 * r0 * pow(((vp_n) / a0), p) * pow((vs_n / b0), q) * pow((rho_n / r0), r))
    return eei