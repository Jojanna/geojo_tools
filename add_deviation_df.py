from scipy.interpolate import interp1d
import numpy as np
import pandas as pd

from geojo.las_handling import load_las

def import_deviation(df, wd, kb, deviation = "deviated", dev_format = "txt", dev_path = None, dev_file = None):

    if deviation == "vertical":
        df["TVDKB"] = df["MDKB"]
        df["TVDSS"] = df["TVDKB"] - kb
        df["TVDML"] = df["TVDKB"] - wd


    elif deviation == "deviated":

        if dev_format == "las":
            out, df["TVDKB"] = dev_las(dev_path, dev_file, df)

        elif dev_format == "txt":
            out, df["TVDKB"] = dev_txt(dev_path, dev_file, df)
        else:
            print("build me!")

        df["TVDSS"] = df["TVDKB"] - kb
        df["TVDml"] = df["TVDSS"] - wd

    else:
        print("no deviation specified")

    return df


## adding deviation to dataframe of well data
def dev_las(dev_path, dev_file, data_df, mdkb = "DEPTH", tvdkb = "TVD"):

    dev, u = load_las(dev_path, dev_file, -999.25)
    dev.dropna(how = "any", inplace = True)

    f = interp1d(dev[mdkb], dev[tvdkb], kind = "cubic", fill_value = "extrapolate", assume_sorted=True)
    out = pd.DataFrame()
    out["MDKB"] = data_df["MDKB"]
    out["TVDKB"] = f(out["MDKB"])
    # print (dev)

    return out, out["TVDKB"]

def dev_txt(dev_path, dev_file, data_df, mdkb = "MD", tvdkb = "TVD", wd = 0, kb = 0):
    dev = pd.read_csv(dev_path + "\\" + dev_file, delim_whitespace = True, skiprows = 0, header = 0 ,encoding = "utf-16")
    # print (dev)
    dev = dev.loc[(dev[mdkb] >= wd + kb)]
    f = interp1d(dev[mdkb], dev[tvdkb], kind = "cubic", fill_value = "extrapolate", assume_sorted=True)

    out = pd.DataFrame()
    out["MDKB"] = data_df["MDKB"]
    out["TVDKB"] = f(out["MDKB"])

    return out, out["TVDKB"]

## adding TWT conversion to dataframe of well data from cs/well tie
def geoview_cs_las(cs_path, cs_file, data_df, mdkb = "DEPTH", vint = "DPTM"):
    cs, u = load_las(cs_path, cs_file, -999.25)

    f = interp1d(cs[mdkb], cs[vint], kind="cubic", fill_value="extrapolate", assume_sorted=True)
    out = pd.DataFrame()
    out["MDKB"] = data_df["MDKB"]
    out["TVDSS"] = data_df["TVDSS"]
    out["TWT"] = f(out["MDKB"])
    out["VInt"] = np.nan
    for i in range(1, len(out["MDKB"]), 1):
        out["VInt"].iloc[i] = abs(
            (out["TVDSS"][i] - out["TVDSS"][i - 1]) / ((out["TWT"][i] - out["TWT"][i - 1]) / 2000))
    out["VInt"].iloc[0] = out["VInt"].iloc[1]

    return out, out["TWT"], out["VInt"]
