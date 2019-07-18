import numpy as np
import pandas as pd
#from rppy.moduli import poissons

from rppy.media import hashin_shtrikman # soft_sand

def hertz_mindlin(u, v,  P, phi, n=None, f = 1):
    """
    Elastic moduli of an elastic sphere pack subject to confining pressure
    given by the Hertz-Mindlin [Mindlin, 1949] theory.

    If the coordination number n is not given, the function uses the empirical
    dependency of n on porosity shown by Murphy [1982]
    """

    if not n:
        n = 20 - 34*phi + 14*phi**2

    Khm = ((n**2*(1-phi)**2*u**2*P) / (18*np.pi**2*(1-v)**2))**(1/3)

    if f == 1:
        uhm = ((5-4*v)/(10-5*v)) * ((3*n**2*(1-phi)**2*u**2*P) /
                                (2*np.pi**2*(1-v)**2))**(1/3)
    else:
        uhm = ((2 + 3 * f - v * (1 + 3*f))/(10-5*v)) * ((3*n**2*(1-phi)**2*u**2*P) /
                                (2*np.pi**2*(1-v)**2))**(1/3)

    return(Khm, uhm)

def soft_sand(Kg, ug, phi, phi_0=0.45, n = None, P=0.3, f = 1):
    """

    :param Kg:
    :param ug:
    :param phi:
    :param phi_0:
    :param n:
    :param P:
    :return:
    """

    """
    The soft-sand (or uncemented-sand, or friable-sand) model calculates the
    bulk and shear moduli of dry sand in which cement is deposited away from
    grain contacts. The framework is a dense random pack of identical spherical
    grains, with critical porosity ~0.36 and average coordination number
    (contacts per grain) C ~ 5-9.

    Hertz-Mindlin theory gives the effective bulk and shear moduli of the pack
    at this porosity, and a heuristic modified Hashin-Shtrikman lower bound is
    used to interpolate at lower porosities.
    """
    from rppy.moduli import poissons
    if n == None:
        n = 20 - 34 * phi_0 + 14 * np.power(phi_0, 2)

    # Moduli at critical porosity
    vg = poissons(K=Kg, u=ug)
    Khm, uhm = hertz_mindlin(ug, vg, P, phi_0, n = n, f = f)

    # Moduli at sub-critical porosity.
    A = (phi/phi_0) / (Khm + 4/3*uhm)
    B = (1 - phi/phi_0) / (Kg + 4/3*uhm)
    C = 4/3*uhm

    D = (phi/phi_0) / (uhm + (uhm/6) * ((9*Khm + 8*uhm) / (Khm + 2*uhm)))
    E = (1 - phi/phi_0) / (ug + (uhm/6) * ((9*Khm + 8*uhm) / (Khm + 2*uhm)))
    F = (uhm/6) * ((9*Khm + 8*uhm) / (Khm + 2*uhm))

    Keff = (A + B)**-1 - C
    ueff = (D + E)**-1 - F

    return(Keff, ueff)


def stiff_sand(Kg, ug, phi, phi_0=0.45, n=None, P=0.3, f = 1):
    """

    :param Kg:
    :param ug:
    :param phi:
    :param phi_0:
    :param n:
    :param P:
    :return:
    """
    """
    The soft-sand (or uncemented-sand, or friable-sand) model calculates the
    bulk and shear moduli of dry sand in which cement is deposited away from
    grain contacts. The framework is a dense random pack of identical spherical
    grains, with critical porosity ~0.36 and average coordination number
    (contacts per grain) C ~ 5-9.

    Hertz-Mindlin theory gives the effective bulk and shear moduli of the pack
    at this porosity, and a heuristic modified Hashin-Shtrikman lower bound is
    used to interpolate at lower porosities.
    """
    from rppy.moduli import poissons

    if n == None:
        n = 20 - 34 * phi_0 + 14 * np.power(phi_0, 2)

    # Moduli at critical porosity
    vg = poissons(K=Kg, u=ug)
    Khm, uhm = hertz_mindlin(ug, vg, P, phi_0, n=n, f=f)

    zhm = (uhm / 6) * ((9 * Khm + 8 * uhm) / (Khm + 2 * uhm))
    z = (ug / 6) * ((9 * Kg + 8 * ug) / (Kg + 2 * ug))

    # Moduli at sub-critical porosity.
    A = (phi / phi_0) / (Khm + 4 / 3 * ug)
    B = (1 - phi / phi_0) / (Kg + 4 / 3 * ug)
    C = 4 / 3 * ug

    D = (phi / phi_0) / (uhm + z)
    E = (1 - phi / phi_0) / (ug + z)
    F = z

    Keff = (A + B) ** -1 - C
    ueff = (D + E) ** -1 - F

    return (Keff, ueff)

def cemented_sand(u, v, p, uc, vc, pc, phi, phi0=0.36, n=None, style='constant'):
    """
    Bulk and shear moduli of dry sand in which dement is deposity at grain
    contacts. The cement is elastic, and may differ from that of the spherical
    pack.

    This may also referred to as the Dvorkin-Nur cement model
    [Dvorkin and Nur, 1996]

    style = 'contact' - All cement at grain contacts (default)
    style = 'constant' - Cement deposited evenly on grain surface

    If the critical porosity os not given, the functions uses 0.36

    If the coordination number n is not given, the function uses the empirical
    dependency of n on porosity shown by Murphy [1982]
    """
    import rppy.moduli as rpmod

    if n == None:
        n = 20 - 34*phi + 14*phi**2

    if style == 'contact':
        a = 2*((phi0 - phi) / (3*n*(1 - phi0)))**0.25
    elif style == 'constant':
        a = ((2*(phi0 - phi)) / (3*(1 - phi0)))**0.5
    else:
        raise ValueError('You must specify either contact or constant cement.')

    Vpc = rpmod.Vp(pc, u=uc, v=vc)
    Vsc = rpmod.Vs(pc, u=uc, v=vc)

    Ln = (2*uc*(1 - v)*(1 - vc)) / (np.pi*u*(1 - 2*vc))
    Lt = uc / (np.pi*u)

    An = -0.024153*Ln**(-1.3646)
    Bn = 0.20405*Ln**(-0.89008)
    Cn = 0.00024649*Ln**(-1.9864)
    Sn = An*a**2 + Bn*a + Cn

    At = -10**(-2)*(2.26*v**2 + 2.07*v + 2.3)*Lt**(0.079*v**2 + 0.1754*v - 1.342)
    Bt = (0.0573*v**2 + 0.0937*v + 0.202)*Lt**(0.0274*v**2 + 0.0529*v - 0.8765)
    Ct = 10**(-4)*(9.654*v**2 + 4.945*v + 3.1)*Lt**(0.01867*v**2 + 0.4011*v - 1.8186)
    St = At*a**2 + Bt*a + Ct

    Mc = pc*Vpc**2
    uc = pc*Vsc**2

    Keff = (1/6)*n*(1 - phi0)*Mc*Sn
    ueff = (3/5)*Keff + (3/20)*n*(1 - phi0)*uc*St

    return (Keff, ueff)


#kw = 2.974050 #2.25 #* pow(10, 9)
#rhow = 1.070539#* pow(10, 3)
#uw = 0

def media_models_df(phi_0 = 0.45, P = 0.3, kg = 36, ug = 35, rhog = 2.65, kc = 36, uc = 35, rhoc = 2.65, kw = 2.25, rhow = 1.0, uw = 0, n = None, cement = None, f = 1):
    """

    :param phi_crit:
    :param P:
    :param kg:
    :param ug:
    :param rhog:
    :param kc:
    :param uc:
    :param rhoc:
    :param kw:
    :param rhow:
    :param uw:
    :param C:
    :param cement:
    :return:
    """

    from .moduli import calc_vp_vs
    from rppy.moduli import poissons
    from.FRM_function_update import k_wet

    if n == None:
        n = 20 - 34 * phi_0 + 14 * np.power(phi_0, 2)

    f1 = np.arange(1 - phi_0, 1, 0.001)
    f2 = 1 - f1

    media_models = pd.DataFrame()
    media_models["PHIE"] = f2

    hs_k_upper = []
    hs_k_lower = []
    hs_mu_upper = []
    hs_mu_lower = []

    # calculate upper/lower Hashin-Shtrickman bounds

    for f1a, f2a in zip(f1, f2):
        K_hi, K_lo, u_hi, u_lo = hashin_shtrikman(np.array([kg, kw]), np.array([ug, uw]), np.array([f1a, f2a]))

        hs_k_upper.append(K_hi)  # / np.power(10, 9))
        hs_k_lower.append(K_lo)  # / np.power(10, 9))
        hs_mu_upper.append(u_hi)  # / np.power(10, 9))
        hs_mu_lower.append(u_lo)  # / np.power(10, 9))

    media_models["RHOB"] = media_models["PHIE"] * rhow + (1 - media_models["PHIE"]) * rhog
    media_models["K_HS_upper"] = hs_k_upper
    media_models["K_HS_lower"] = hs_k_lower
    media_models["MU_HS_upper"] = hs_mu_upper
    media_models["MU_HS_lower"] = hs_mu_lower

    media_models["VP_HS_upper"], media_models["VS_HS_upper"] = calc_vp_vs(media_models["K_HS_upper"],
                                                                          media_models["MU_HS_upper"],
                                                                          media_models["RHOB"])
    media_models["VP_HS_lower"], media_models["VS_HS_lower"] = calc_vp_vs(media_models["K_HS_lower"],
                                                                          media_models["MU_HS_lower"],
                                                                          media_models["RHOB"])

    media_models["AI_HS_upper"] = media_models["VP_HS_upper"] * media_models["RHOB"]
    media_models["AI_HS_lower"] = media_models["VP_HS_lower"] * media_models["RHOB"]

    vg = poissons(K=kg, u=ug)  # poissons ratio for grains/matrix (assumed 100% qtz)
    vc = poissons(K=kc, u=uc)  # poissions ratio for cement (assumed 100% qtz)
    Khm, uhm = hertz_mindlin(ug, vg, P, phi_0, n=n, f = f)

    # soft sand model

    media_models["K_SS_dry"], media_models["MU_SS"] = soft_sand(kg, ug, media_models["PHIE"], phi_0=phi_0, n=n,
                                                                P=P)
    media_models["K_SS"] = k_wet(media_models["K_SS_dry"], kw, media_models["PHIE"], kg)
    media_models["VP_SS"], media_models["VS_SS"] = calc_vp_vs(media_models["K_SS"], media_models["MU_SS"],
                                                              media_models["RHOB"])
    media_models["AI_SS"] = media_models["VP_SS"] * media_models["RHOB"]

    # stiff sand model

    media_models["K_stiff_dry"], media_models["MU_stiff"] = stiff_sand(kg, ug, media_models["PHIE"], phi_0=phi_0,
                                                                       n=n, P=P)
    media_models["K_stiff"] = k_wet(media_models["K_stiff_dry"], kw, media_models["PHIE"], kg)
    media_models["VP_stiff"], media_models["VS_stiff"] = calc_vp_vs(media_models["K_stiff"], media_models["MU_stiff"],
                                                                    media_models["RHOB"])
    media_models["AI_stiff"] = media_models["VP_stiff"] * media_models["RHOB"]

    C = 20 - 34 * phi_0 + 14 * np.power(phi_0, 2)  # Murphy's law?)

    media_models["K_contactc_dry"], media_models["MU_contactc"] = cemented_sand(ug, vg, rhog, uc, vc, rhoc, media_models["PHIE"], phi0=phi_0, n=6, style='constant')
    media_models["K_contactc"] = k_wet(media_models["K_contactc_dry"], kw, media_models["PHIE"], kg)
    media_models["VP_contactc"], media_models["VS_contactc"] = calc_vp_vs(media_models["K_contactc"], media_models["MU_contactc"], media_models["RHOB"])
    media_models["AI_contactc"] = media_models["VP_contactc"] * media_models["RHOB"]

    print ("I'm here!")
    # constant cement model
    if cement == None:

        cement = np.arange(0.01, 0.1, 0.01)

    for c in cement:
        # find high porosity end member
        #cemented_sand(u, v, p, uc, vc, pc, phi, phi0=0.36, n=None, style='constant')
        phi = phi_0 - c
        #print ("ug= %f, vg = %f, rhog = %f, uc = %f, vc = %f, rhoc= %f, phi = %f, phi_0 = %f, n = %f" % (ug, vg, rhog, uc, vc, rhoc, phi, phi_0, n))
        k_contactc_dry, mu_contactc = cemented_sand(ug, vg, rhog, uc, vc, rhoc, phi, phi0=phi_0, n = n, style='constant')
        k_contactc_wet = k_wet(k_contactc_dry, kw, phi - c, kg)

        hs_k_upper = []
        hs_k_lower = []
        hs_mu_upper = []
        hs_mu_lower = []

        for phi_n in media_models["PHIE"].tolist():
            # hi_fraction =
            K_hi, K_lo, u_hi, u_lo = hashin_shtrikman(np.array([kg, k_contactc_wet]), np.array([ug, mu_contactc]),
                                                      np.array([1 - phi_n / (phi_0 - c), phi_n / (phi_0 - c)]))
            hs_k_upper.append(K_hi)
            hs_k_lower.append(K_lo)
            hs_mu_upper.append(u_hi)
            hs_mu_lower.append(u_lo)

        media_models["RHOB_constant_" + str(c)] = c * rhoc + media_models["PHIE"] * rhow + (
                1 - c - media_models["PHIE"]) * rhog

        media_models["K_constant_" + str(c)] = hs_k_lower
        media_models["MU_constant_" + str(c)] = hs_mu_lower
        media_models["VP_constant_" + str(c)], media_models["VS_constant_" + str(c)] = calc_vp_vs(
            media_models["K_constant_" + str(c)], media_models["MU_constant_" + str(c)],
            media_models["RHOB_constant_" + str(c)])
        media_models["AI_constant_" + str(c)] = media_models["VP_constant_" + str(c)] * media_models[
            "RHOB_constant_" + str(c)]

    return media_models


def plot_media_models(media_models):

    fig, ((ax1, ax2, ax3), (ax4, ax5, ax6)) = plt.subplots(2, 3, figsize=(30, 20))
    fig.set_size_inches(18, 12)

    ax3.plot(media_models["PHIE"], media_models["RHOB"], linewidth=0.5)
    for c in cement:
        p = c * 100
        ax1.plot(media_models["PHIE"], media_models["K_constant_" + str(c)], linewidth=0.5,
                 label="Constant Cement %.0f%%" % p)
        ax2.plot(media_models["PHIE"], media_models["MU_constant_" + str(c)], linewidth=0.5,
                 label="Constant Cement %.0f%%" % p)

        ax4.plot(media_models["PHIE"], media_models["VP_constant_" + str(c)], linewidth=0.5,
                 label="Constant Cement %.0f%%" % p)
        ax5.plot(media_models["PHIE"], media_models["VS_constant_" + str(c)], linewidth=0.5,
                 label="Constant Cement %.0f%%" % p)
        ax6.plot(media_models["PHIE"], media_models["AI_constant_" + str(c)], linewidth=0.5,
                 label="Constant Cement %.0f%%" % p)

    for model, label in zip(["contactc", "HS_upper", "HS_lower", "SS", "stiff"],
                            ["Contact Cement, wet", "HS Upper", "HS Lower", "Soft Sand, wet", "Stiff Sand, wet"]):
        ax1.plot(media_models["PHIE"], media_models["K_" + model], label=label, linewidth=0.5)
        ax2.plot(media_models["PHIE"], media_models["MU_" + model], label=label, linewidth=0.5)

        ax4.plot(media_models["PHIE"], media_models["VP_" + model], label=label, linewidth=0.5)
        ax5.plot(media_models["PHIE"], media_models["VS_" + model], label=label, linewidth=0.5)
        ax6.plot(media_models["PHIE"], media_models["AI_" + model], label=label, linewidth=0.5)

    return fig

