import pandas as pd
import numpy as np

def slow_to_velocity(dt):
    vp = list(map(lambda d: ((1/d) * 0.3048 * 10 ** 6), dt))

    return vp

def ft_to_m (feet):
    metres = feet * 0.3048
    return metres


def psi_to_pa(psi):
    pa = psi * 6894.76
    return psi

def pa_to_psi(pa):
    psi = pa * 0.000145038

