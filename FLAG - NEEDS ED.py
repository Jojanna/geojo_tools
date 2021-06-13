import numpy as np




def api (rho0):
    a = 141.5 / rho0 - 131.5
    return a

def rho0 (api):
    r = 141.5 /(api + 131.5)
    return r



def batzle_wang(P, T, fluid, S=None, G=None, api=None, Rg=None):
    """
    Calculate the elastic properties of reservoir fluids using the
    Batzle & Wang [1992] equations.
    :param P: Pressure (MPa)
    :param T: Temperature {deg C)
    :param fluid: Fluid type to calculate: brine, gas, or oil
    :param S: Salinity (brine only, in ppm)
    :param G: Gas gravity (gas mode only, ratio of gas density to air density
              at 15.6C and atmospheric pressure)
    :param api: American Petroleum Insitute (API) oil gravity
    :param Rg: Gas-oil ratio
            #GOR ratio (bbl/bbl = GOR(scf/bbl) * 17810760667903525
    """

    if fluid == 'brine':
        S = S / pow(10, 6)  # ppm to fraction of one
        w = np.array([
            [1402.85, 1.524, 3.437e-3, -1.197e-5],
            [4.871, -0.0111, 1.739e-4, -1.628e-6],
            [-0.04783, 2.747e-4, -2.135e-6, 1.237e-8],
            [1.487e-4, -6.503e-7, -1.455e-8, 1.327e-10],
            [-2.197e-7, 7.987e-10, 5.230e-11, -4.614e-13],
        ])

        rhow = (1 + (10 ** -6) * (-80 * T - 3.3 * pow(T, 2) + 0.00175 * pow(T, 3) +
                                  489 * P - 2 * T * P + 0.016 * pow(T, 2) * P - (1.3e-5) * pow(T, 3) * P -
                                  0.333 * pow(P, 2) - 0.002 * T * pow(P, 2)))

        rhob = rhow + S * (0.668 + 0.44 * S + pow(10, -6) * (300 * P - 2400 * P * S +
                                                             T * (80 + 3 * T - 3300 * S - 13 * P + 47 * P * S)))

        Vw = 0
        for i in range(4):
            for j in range(3):
                Vw = Vw + w[i][j] * pow(T, i) * pow(P, j)

        Vb = (Vw + S * (1170 - 9.8 * T + 0.055 * pow(T, 2) - 8.5e-5 * pow(T, 3) + 2.6 * P -
                        0.0029 * T * P - 0.0476 * pow(P, 2) + pow(S, (3 / 2)) * (780 - 10 * P + 0.16 * pow(P, 2)) -
                        1820 * pow(S, 2)))

        kb = (pow(Vb, 2) * rhob * 1000) / (4 / 3)

        out = (rhob, kb)

    elif fluid == 'oil':
        rho0 = 141.5 / (api + 131.5) # <-- oil reference density, derived from api?
        # print(rho0)
        # print(G)
        # print(Rg)

        V = A - B * T + C * T+ D * T * P
        #V = velocity of dead oil, i.e. no dissolved gas. "gas free", GOR = 0
        A = 2090 * pow(rho0/ (2.6 - rho0), 0.5)
        B = 3.7
        C = 4.64
        D = 0.0115

        # Rg = GOR

        #rho_pv = velocity pseudo density

        rho_a = 0.61731 * pow(10, -0.00326 * api) + 1.5177 - 0.54349 *
        vg = Rg #??
        e = 0.113 # gas parameter
        rho_pv = rho0 * (1 - vg) + e * rho_a * vg


        B0 = 0.972 + 0.00038 * pow((2.4 * Rg * pow((G / rho0), 0.5) + T + 17.8), (1.175))

        # input to calculation of velocicty
        rho_r = (rho0 / B0) * (1 + 0.001 * Rg) ** -1  # pseudo-density of oil

        # input to calculation of density
        rhog = (rho0 + 0.0012 * G * Rg) / B0  # density of oil with gas
        rhop = (rhog + (0.00277 * P -  # correct for pressure
                        1.71e-7 * P ** 3) * (rhog - 1.15) ** 2 + 3.49e-4 * P)

        rho = rhop / (0.972 + 3.81e-4 * (T + 17.78) ** 1.175)  # correct for temp
        Vp = 2096 * (rho_r / (2.6 - rho_r)) ** 0.5 - 3.7 * T + 4.64 * P
        ko = (pow(Vp, 2) * rho * 1000) / (4 / 3)
        # print (Vp)
        # print (ko)

        out = (rho, ko)

    elif fluid == 'gas':
        Ta = T + 273.15  # absolute temperature
        Pr = P / (4.892 - 0.4048 * G)  # pseudo-pressure
        Tr = Ta / (94.72 + 170.75 * G)  # pseudo-temperature

        R = 8.31441
        d = np.exp(-(0.45 + 8 * (0.56 - 1 / Tr) ** 2) * Pr ** 1.2 / Tr)
        c = 0.109 * (3.85 - Tr) ** 2
        b = 0.642 * Tr - 0.007 * Tr ** 4 - 0.52
        a = 0.03 + 0.00527 * (3.5 - Tr) ** 3
        m = 1.2 * (-(0.45 + 8 * (0.56 - 1 / Tr) ** 2) * Pr ** 0.2 / Tr)
        y = (0.85 + 5.6 / (Pr + 2) + 27.1 / (Pr + 3.5) ** 2 -
             8.7 * np.exp(-0.65 * (Pr + 1)))
        f = c * d * m + a
        E = c * d
        Z = a * Pr + b + E

        rhog = (28.8 * G * P) / (Z * R * Ta)
        Kg = P * y / (1 - Pr * f / Z)

        out = (rhog, Kg)
    else:
        out = None

    return (out)
