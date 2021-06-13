import numpy as np

#  castagana
# gardner
#... what else...



def rhob_gardner_relation(vp_data, litho = "shale"):
    """

    :param litho: "shale", "sst", "lst", "avg"
    :return:
    """
    if litho == "shale":
        rhob_gardner = rhob_from_gardner_shale(vp_data)
    elif litho == "sst":
        rhob_gardner = rhob_from_gardner_sst(vp_data)
    elif litho == "lst":
        rhob_gardner = rhob_from_gardner_lst(vp_data)
    elif litho == "avg":
        rhob_gardner = rhob_from_gardner_avg(vp_data)
    else:
        print ("Specify lithology!")

    return rhob_gardner


def rhob_from_gardner_avg(vp_data, a = 1.741,  b = 0.25):

    rhob_gardner = a * np.power(vp_data / 1000, b)

    return rhob_gardner

def rhob_from_gardner_sst(vp_data, a = 1.66,  b = 0.261):

    rhob_gardner = a * np.power(vp_data / 1000, b)

    return rhob_gardner

def rhob_from_gardner_shale(vp_data, a = 1.75,  b = 0.265):

    rhob_gardner = a * np.power(vp_data / 1000, b)

    return rhob_gardner

def rhob_from_gardner_lst(vp_data, a = 1.36,  b = 0.386):

    rhob_gardner = a * np.power(vp_data / 1000, b)

    return rhob_gardner

def rhob_from_gardner_dolomite(vp_data, a = 1.74,  b = 0.252):
    rhob_gardner = a * np.power(vp_data / 1000, b)
    return rhob_gardner


def rhob_gardner_avg(vp_data, a = 1.741,  b = 0.25):

    rhob_gardner = a * np.power(vp_data / 1000, b)

    return rhob_gardner

def rhob_gardner_sst(vp_data, a = 1.66,  b = 0.261):

    rhob_gardner = a * np.power(vp_data / 1000, b)

    return rhob_gardner

def rhob_gardner_shale(vp_data, a = 1.75,  b = 0.265):

    rhob_gardner = a * np.power(vp_data / 1000, b)

    return rhob_gardner

def rhob_gardner_lst(vp_data, a = 1.36,  b = 0.386):

    rhob_gardner = a * np.power(vp_data / 1000, b)

    return rhob_gardner

def rhob_gardner_lst2(vp_data, a = 1.5,  b = 0.225):

    rhob_gardner = a * np.power(vp_data / 1000, b)

    return rhob_gardner

def rhob_gardner_dolomite(vp_data, a = 1.74,  b = 0.252):
    rhob_gardner = a * np.power(vp_data / 1000, b)
    return rhob_gardner



def vp_gardner_relation(rhob_data, litho = "shale"):
    """

    :param litho: "shale", "sst", "lst", "avg"
    :return:
    """
    if litho == "shale":
        vp_gardner = vp_from_gardner_shale(rhob_data)
    elif litho == "sst":
        vp_gardner = vp_from_gardner_sst(rhob_data)
    elif litho == "lst":
        vp_gardner = vp_from_gardner_lst(rhob_data)
    elif litho == "avg":
        vp_gardner = vp_from_gardner_avg(rhob_data)
    else:
        print ("Specify lithology!")

    return vp_gardner


def vp_from_gardner_avg(rhob_data, a=1.741, b=0.25):
    vp_gardner = np.power((rhob_data / a), 1/b) * 1000

    return vp_gardner

def vp_from_gardner_sst(rhob_data, a=1.66, b=0.261):
    vp_gardner = np.power((rhob_data / a), 1 / b) * 1000

    return vp_gardner

def vp_from_gardner_shale(rhob_data, a=1.75, b=0.265):
    vp_gardner = np.power((rhob_data / a), 1 / b) * 1000

    return vp_gardner

def vp_from_gardner_lst(rhob_data, a=1.36, b=0.386):
    vp_gardner = np.power((rhob_data / a), 1 / b) * 1000

    return vp_gardner



def vs_castagna_relation(vp_data, litho = "shale"):
    """

    :param litho: "shale", "sst", "lst", "dolomite", "han_sst"
    :return:
    """
    if litho == "shale":
        vs_castagna = vs_from_castagna_shale(vp_data)
    elif litho == "sst":
        vs_castagna = vs_from_castagna_sst(vp_data)
    elif litho == "lst":
        vs_castagna = vs_from_castagna_lst(vp_data)
    elif litho == "dolomite":
        vs_castagna = vs_from_castagna_dolomite(vp_data)
    elif litho == "han_sst":
        vs_castagna = vs_from_han_sst(vp_data)
    else:
        print ("Specify lithology!")

    return vs_castagna


# VP-VS relations. VS from VP

def vs_from_castagna_shale(vp_data, m = 0.862, c = -1.172):

    vs_castagna_shale = ((vp_data / 1000) * m + c) * 1000

    return vs_castagna_shale

def vs_from_castagna_sst(vp_data, m = 0.804, c = -0.856):

    vs_castagna_sst = ((vp_data / 1000) * m + c) * 1000

    return vs_castagna_sst

def vs_from_castagna_dolomite(vp_data, m = 0.583, c = -0.078):

    vs_castagna_dolo = ((vp_data / 1000) * m + c) * 1000

    return vs_castagna_dolo

def vs_from_castagna_lst(vp_data, a = -0.055, b = 1.017, c = -1.031):

    vs_castagna_lst = ((pow((vp_data / 1000), 2) * a + b * (vp_data / 1000) + c)) * 1000

    return vs_castagna_lst

def vs_from_han_sst(vp_data, m = 0.794, c = -0.787):

    vs_from_han_sst = ((vp_data / 1000) * m + c) * 1000

    return vs_from_han_sst



def quadratic_equ (y, a, b, c):

    c = c - y
    x = -b + (np.power(b, 2) - 4 * a * c) / 2 * a
    return x

def vp_castaga_relation(vs_data, litho = "shale"):
    """

    :param litho: "shale", "sst", "lst", "dolomite", "han_sst"
    :return:
    """
    if litho == "shale":
        vp_castagna = vp_from_castagna_shale(vs_data)
    elif litho == "sst":
        vp_castagna = vp_from_castagna_sst(vs_data)
    elif litho == "lst":
        vp_castagna = vp_from_castagna_lst(vs_data)
    elif litho == "dolomite":
        vp_castagna = vp_from_castagna_dolomite(vs_data)
    elif litho == "han_sst":
        vp_castagna = vp_from_han_sst(vs_data)
    else:
        print ("Specify lithology!")

    return vp_castagna


def vp_from_castagna_shale(vs_data, m = 0.862, c = -1.172):

    vp_castagna_shale = (vs_data / 1000 - c) / m * 1000

    return vp_castagna_shale

def vp_from_castagna_sst(vs_data, m = 0.804, c = -0.856):

    vp_castagna_sst = (vs_data / 1000 - c) / m * 1000

    return vp_castagna_sst

def vp_from_castagna_dolomite(vs_data, m = 0.583, c = -0.078):

    vp_castagna_dolo = (vs_data / 1000 - c) / m * 1000

    return vp_castagna_dolo

def vp_from_castagna_lst(vs_data, a = -0.055, b = 1.017, c = -1.031):
    vs_data = vs_data / 1000
    x = quadratic_equ(vs_data, a, b, c)
    vp_castagna_lst = x * 1000

    return vp_castagna_lst

def vp_from_han_sst(vs_data, m = 0.794, c = -0.787):

    vp_from_han_sst= (vs_data / 1000 - c) / m * 1000

    return vp_from_han_sst


# calculation of co ordination number from porosity
def murphys_law(phi):
    n = 20 - 34 * phi + 14 * phi ** 2
    return n


