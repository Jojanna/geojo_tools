import numpy as np

def k_mu(vp, vs, rhob):
    """

    :param vp:
    :param vs:
    :param rhob:
    :return:
    """
    mu = np.power(vs, 2) * rhob * 1000
    k = np.power(vp, 2) * rhob * 1000 - 4 / 3 * mu
    mu = mu / np.power(10, 9)
    k = k / np.power(10, 9)

    return k, mu


def poisson(k, mu):
    """

    :param k:
    :param mu:
    :return:
    """
    pr = (3 * k - 2 * mu) / (2 * (3 * k + mu))
    return pr


def calc_vp_vs(k, mu, rho):
    """

    :param k:
    :param mu:
    :param rho:
    :return:
    """
    k = k * 10 ** 9
    mu = mu * 10 ** 9
    vp = []
    vs = []
    for k_n, mu_n, rho_n in zip(k, mu, rho):
        vs.append((pow(mu_n / (rho_n * 1000), 0.5)))
        vp.append((pow((mu_n * 4 / 3 + k_n) / (rho_n * 1000), 0.5)))
    return vp, vs