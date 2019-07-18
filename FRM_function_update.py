import numpy as np
import pandas as pd
from rppy.fluid import batzle_wang


"""
Contents:
* k_wet
* fluids_calc
* reuss fluids
*rho_dry_calc


"""


def k_wet (k_dry, k_fluid, phiT, k_HS): # gassmanns
    #print ("k_fluid")
    #print (k_fluid)
    rhs = k_dry / (k_HS - k_dry) + k_fluid / (phiT * (k_HS - k_fluid))
    k_sat = rhs * k_HS / (1 + rhs)
    return k_sat


def fluids_calc(pressure, temp, fl, s, g, api_oil, ratio):
    from geojo.moduli import k_mu, calc_vp_vs
    w = batzle_wang(pressure, temp, 'brine', S=s, G=g, api=api_oil, Rg=ratio)
    # print (w)
    rho_brine = w["rho"]
    vp_brine = w["Vp"]
    k_brine, mu_brine = k_mu(vp_brine, 0, rho_brine)
    # print ("Kw: %f" % k_brine)
    if fl == 'oil':
        h = batzle_wang(pressure, temp, 'oil', S=s, G=g, api=api_oil, Rg=ratio)
        # print (h)
        rho_h = h["rho"]
        vp_h = h["Vp"]
        k_h, mu_h = k_mu(vp_h, 0, rho_h)
        # print ("KH: %f" % k_h)
    elif fl == 'gas':
        h = batzle_wang(pressure, temp, 'gas', S=s, G=g, api=api_oil, Rg=ratio)
        # print (h)
        rho_h = h["rho"]
        k_h = h["K"] / np.array(10 ** 3) # output from BW is in MPA --> supply in GPa
        # print ("KH: %f" % k_h)

    else:
        print('check fluids!')
        rho_h, vp_h = 0, 0
        exit(0)

    # k_brine, k_h = [k_brine, k_h] * np.array(10**9)

    return rho_brine, k_brine, rho_h, k_h


def reuss_fluids(rho_brine, k_brine, rho_h, k_h, sw_data):
    k_f = []
    rho_f = []
    for s in sw_data:
        # s_we = 2 * np.exp(-11 * p)
        # print (s_we)
        k_f_s = 1 / (s / k_brine + (1 - s) / k_h)  # Reuss average
        k_f.append(k_f_s)
        rho_f_p = s * rho_brine + (1 - s) * rho_h
        rho_f.append(rho_f_p)

    return k_f, rho_f

def rho_dry_calc (rho_data, rho_f_ini, phie_data):
    rho_dry = []
    for rho_rock, rho_f, p in zip(rho_data, rho_f_ini, phie_data):
        rho_dry.append(rho_rock - p * rho_f)

    return rho_dry

def dry_rock(phie_data,k_sat,k_f, k_ma):
    k_pore_data = []
    k_d_data = []

    pore_params = zip(phie_data,k_sat,k_f, k_ma)

    for phie_n, k_sat_n, k_f_n, k_ma_n in pore_params:
        k_pore_n = phie_n /(1/k_sat_n - 1/k_ma_n) - 1 / (1/(k_f_n) - 1/k_ma_n)
        k_pore_data.append(k_pore_n)

    dry_params = zip(phie_data, k_pore_data, k_ma)
    for phie_n, k_pore_n, k_ma_n in dry_params:
        k_d_n = 1/(phie_n/k_pore_n + 1/k_ma_n)
        k_d_data.append(k_d_n)

    return k_pore_data, k_d_data

def fluid_sub_k(k_ma, phie_data, k_pore_data, k_f):

    sub_params = zip(k_ma, phie_data, k_pore_data, k_f)
    k_out = []
    #phie_exc_0 = []
    for k_ma_n, phie_n, k_pore_n, k_f_n in sub_params:
        if phie_n > 0:
            k_out_n = 1 / (1/k_ma_n + phie_n/(k_pore_n + (k_ma_n * k_f_n)/(k_ma_n - k_f_n)))
            #k_out.append(k_out_n)
            #phie_exc_0.append(phie_n)
        else:
            k_out_n = 1 / (1/k_ma_n)
        k_out.append(k_out_n)

    return k_out


def multiple_FRM(phie_data, sw_out, k_ma, kphi_set, rho_dry, pressure_out, temp, fl_out, s, g, api_oil, ratio):
    rho_brine, k_brine, rho_h_out, k_h_out = fluids_calc(pressure_out, temp, fl_out, s, g, api_oil, ratio)

    k_f_out = list(map(lambda sw_out_n: 1 / (sw_out_n / k_brine + (1 - sw_out_n) / k_h_out), sw_out))
    rho_f_out = list(map(lambda sw_out_n: sw_out_n * rho_brine + (1 - sw_out_n) * rho_h_out, sw_out))

    k_out = []
    for k_ma_n, phie_n, k_pore_n, k_f_n in zip(k_ma, phie_data, kphi_set, k_f_out):
        if phie_n > 0:
            k_out_n = 1 / (1 / k_ma_n + phie_n / (k_pore_n + (k_ma_n * k_f_n) / (k_ma_n - k_f_n)))
            # k_out.append(k_out_n)
            # phie_exc_0.append(phie_n)
        else:
            k_out_n = 1 / (1 / k_ma_n)
        k_out.append(k_out_n)

    rho_out = []
    for rhod, phie, rhof in zip(rho_dry, phie_data, rho_f_out):
        rho_out.append(rhod + phie * rhof)

    k_out_norm = [(k / k_ma_n) for k, k_ma_n in zip(k_out, k_ma)]

    return k_f_out, rho_f_out, k_out, rho_out, k_out_norm