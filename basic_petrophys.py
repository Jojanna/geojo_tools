import numpy as np
import pandas as pd

"""
https://www.cgg.com/data/1/rec_docs/883_shale_volume_calculation.pdf
The gamma ray log has several nonlinear empirical responses as well a linear
responses. The non linear responses are based on geographic area or formation age.
All non linear relationships are more optimistic that is they produce a shale volume
value lower than that from the linear equation.
"""

# linear estimation of Vsh from GR

def vsh_from_gr_linear(gr_data, gr_shale, gr_sand):

    vsh_linear = list(map(lambda gr: (gr - gr_sand) / (gr_shale - gr_sand), gr_data))
    vsh_linear = list(map(lambda vsh: np.where(vsh > float(1), float(1), vsh), vsh_linear))
    print (np.max(vsh_linear), np.min(vsh_linear))
    vsh_linear = list(map(lambda vsh: np.where(vsh < 0, 0, vsh), vsh_linear))

    return vsh_linear

def vsh_from_gr_steiber (gr_data, gr_shale, gr_sand):

    vsh_linear = vsh_from_gr_linear(gr_data, gr_shale, gr_sand)

    vsh_steiber = list(map(lambda vsh: vsh / (3 - 2 * vsh), vsh_linear))

    return vsh_steiber

#from Larinov 1969 for Tertiary rocks
def vsh_from_gr_larinov (gr_data, gr_shale, gr_sand):

    vsh_linear = vsh_from_gr_linear(gr_data, gr_shale, gr_sand)

    vsh_larinov = 0.083 * (pow(2, 2.71 * vsh_linear))

    return vsh_larinov
