from scipy.interpolate import interp1d
import numpy as np
import pandas as pd

from geojo.las_handling import load_las



def import_tdr(df, time_depth = "checkshot", table_format = "ascii_3col", start_time = 0, cs_path = None, cs_file = None):

    if time_depth == "VP":
        df["VInt"] = df["VP"]

        df["TWT"].iloc[0] = start_time
        # data["TWT"].iloc[1] = wd_twt #(data["TVDSS"].iloc[0] / data["VInt"].iloc[0]) * 2000
        samp = range(1, len(df.index))
        df["TWT"] = np.nan
        df["dz"] = np.nan
        df["dt"] = np.nan

        for d in samp:
            dz = df["TVDSS"].iloc[d] - df["TVDSS"].iloc[d - 1]
            df["dz"].iloc[d] = dz
            df["dt"].iloc[d] = dz / df["VInt"].iloc[d] * 2000
            df["TWT"].iloc[d] = df["TWT"].iloc[d - 1] + df["dt"].iloc[d]

    elif time_depth == "checkshot":
        if table_format == "geoview_txt":
            cs, df["TWT"], df["VInt"] = geoview_cs_txt(cs_path, cs_file, df, start_time)
        elif table_format == "las":
            cs, df["TWT"], df["VInt"] = geoview_cs_las(cs_path, cs_file, df)
        elif table_format == "dug_vint_txt":
            cs, df["TWT"], df["VInt"] = dug_vint_txt(cs_path, cs_file, df)
        elif table_format == "ascii_3col":
            cs, df["TWT"], df["VInt"] = ascii_3col(cs_path, cs_file, df)
        elif table_format == "ascii_3col_trace":
            cs, df["TWT"], df["VInt"] = ascii_3col_trace(cs_path, cs_file, df)

        else:
            print("build me!")



    else:
        print("no time-depth specified")

    return df




def geoview_cs_txt(cs_path, cs_file, data, start_time = 0):
    skiprows = 16
    header = None
    cs = pd.read_csv(filepath_or_buffer=cs_path + "\\" + cs_file, delim_whitespace=True, header=header,
                     skiprows=skiprows)
    cs.columns = ["MD (from Surface)", "MDKB", "TVDSS", "TVDKB", "TVD (from SRD)", "TWT", "x (m)", "y (m)",
                  "Absolute x (m)", "Absolute y (m)"]
    cs["VInt"] = np.nan
    for i in range(1, len(cs["VInt"]), 1):
        cs["VInt"].iloc[i] = abs(
            (cs["TVDSS"].iloc[i] - cs["TVDSS"].iloc[i - 1]) / ((cs["TWT"].iloc[i] - cs["TWT"].iloc[i - 1]) / 2000))
    cs["VInt"].iloc[0] = cs["VInt"].iloc[1]

    upper_fill, lower_fill = cs["VInt"].iloc[0], cs["VInt"].iloc[-1]
    f = interp1d(cs["TVDSS"], cs["VInt"], kind="cubic", bounds_error=False, fill_value=(upper_fill, lower_fill))

    out = pd.DataFrame()
    out["MDKB"] = data["MDKB"]
    out["TVDSS"] = data["TVDSS"]
    out["VInt"] = abs(f(out["TVDSS"]))
    out["TWT"] = np.nan

    out["TWT"].iloc[0] = start_time
    for i in range(1, len(out["VInt"]), 1):
        out["TWT"].iloc[i] = (((out["TVDSS"].iloc[i] - out["TVDSS"].iloc[i - 1]) / out["VInt"].iloc[i]) * 2000) + \
                             out["TWT"].iloc[i - 1]

    return out, out["TWT"], out["VInt"]


def dug_vint_txt(cs_path, cs_file, data):
    cs = pd.read_csv(filepath_or_buffer=cs_path + "\\" + cs_file, delim_whitespace=True, header=None, skiprows=1)
    cs.columns = ["TVDSS", "VInt"]
    cs["TWT"] = np.nan

    cs["TWT"].iloc[0] = twt_zero_md
    for i in range(1, len(cs["VInt"]), 1):
        cs["TWT"].iloc[i] = (((cs["TVDSS"].iloc[i] - cs["TVDSS"].iloc[i - 1]) / cs["VInt"].iloc[i]) * 2000) + \
                            cs["TWT"].iloc[i - 1]

    out = pd.DataFrame()
    out["MDKB"] = data["MDKB"]
    out["TVDSS"] = data["TVDSS"]
    out["VInt"] = np.nan
    out["TWT"] = np.nan

    # upper_fill, lower_fill = cs["TWT"].iloc[0], cs["TWT"].iloc[-1]
    f = interp1d(cs["TVDSS"], cs["TWT"], kind="cubic", bounds_error=False, fill_value="extrapolate", assume_sorted=True)
    out["TWT"] = f(out["TVDSS"])
    for i in range(1, len(out["VInt"]), 1):
        out["VInt"].iloc[i] = (out["TVDSS"].iloc[i] - out["TVDSS"].iloc[i - 1]) / (
                (out["TWT"].iloc[i] - out["TWT"].iloc[i - 1]) / 2000)
    out["VInt"].iloc[0] = out["VInt"].iloc[1]

    return out, out["TWT"], out["VInt"]


def acsii_3col_dug(cs_path, cs_file, data):
    skiprows = 1
    header = None
    cs = pd.read_csv(filepath_or_buffer=cs_path + "\\" + cs_file, delim_whitespace=True, header=header,
                     skiprows=skiprows)
    cs.columns = ["MDSS", "TVDSS", "TWT"]

    cs["VInt"] = np.nan
    for i in range(1, len(cs["VInt"]), 1):
        cs["VInt"].iloc[i] = abs(
            (cs["TVDSS"].iloc[i] - cs["TVDSS"].iloc[i - 1]) / ((cs["TWT"].iloc[i] - cs["TWT"].iloc[i - 1]) / 2000))
    cs["VInt"].iloc[0] = cs["VInt"].iloc[1]

    upper_fill, lower_fill = cs["VInt"].iloc[0], cs["VInt"].iloc[-1]
    f = interp1d(cs["TVDSS"], cs["VInt"], kind="linear", bounds_error=False, fill_value=(upper_fill, lower_fill))

    out = pd.DataFrame()
    out["MDKB"] = data["MDKB"]
    out["TVDSS"] = data["TVDSS"]
    out["VInt"] = abs(f(out["TVDSS"]))
    out["TWT"] = np.nan

    out["TWT"].iloc[0] = twt_zero_md
    for i in range(1, len(out["VInt"]), 1):
        out["TWT"].iloc[i] = (((out["TVDSS"].iloc[i] - out["TVDSS"].iloc[i - 1]) / out["VInt"].iloc[i]) * 2000) + \
                             out["TWT"].iloc[i - 1]

    return out, out["TWT"], out["VInt"]


def ascii_3col(cs_path, cs_file, data, header=None, skiprows=1):
    cs = pd.read_csv(filepath_or_buffer=cs_path + "\\" + cs_file, delim_whitespace=True, header=header,
                     skiprows=skiprows)
    cs.columns = ["MDKB", "TVDSS", "TWT"]
    cs["VInt"] = np.nan

    for i in range(1, len(cs["VInt"]), 1):
        cs["VInt"].iloc[i] = abs(
            (cs["TVDSS"].iloc[i] - cs["TVDSS"].iloc[i - 1]) / ((cs["TWT"].iloc[i] - cs["TWT"].iloc[i - 1]) / 2000))
    cs["VInt"].iloc[0] = cs["VInt"].iloc[1]

    f_twt = interp1d(cs["TVDSS"], cs["TWT"], kind="linear", bounds_error=False, fill_value="extrapolate")
    f_vint = interp1d(cs["TVDSS"], cs["VInt"], kind="linear", bounds_error=False, fill_value="extrapolate")

    out = pd.DataFrame()
    out["MDKB"] = data["MDKB"]
    out["TVDSS"] = data["TVDSS"]
    out["VInt"] = abs(f_vint(out["TVDSS"]))
    out["TWT"] = abs(f_twt(out["TVDSS"]))

    # for i in range(1, len(out["VInt"]), 1):
    # out["VInt"].iloc[i] = (((out["TVDSS"].iloc[i] - out["TVDSS"].iloc[i-1]) /  ((out["TWT"].iloc[i] - out["TWT"].iloc[i-1]) / 2000)))

    return out, out["TWT"], out["VInt"]


def ascii_3col_trace(cs_path, cs_file, data):
    import pandas as pd
    skiprows = 1
    header = None
    cs = pd.read_csv(filepath_or_buffer=cs_path + "\\" + cs_file, delim_whitespace=True, header=header,
                     skiprows=skiprows)
    cs.columns = ["TVDSS", "VInt", "TWT"]

    f_twt = interp1d(cs["TVDSS"], cs["TWT"], kind="linear", bounds_error=False, fill_value="extrapolate")
    f_vint = interp1d(cs["TVDSS"], cs["VInt"], kind="linear", bounds_error=False, fill_value="extrapolate")

    out = pd.DataFrame()
    out["MDKB"] = data["MDKB"].values.tolist()
    out["TVDSS"] = data["TVDSS"].values.tolist()
    out["VInt"] = abs(f_vint(out["TVDSS"]))
    out["TWT"] = abs(f_twt(out["TVDSS"]))

    # for i in range(1, len(out["VInt"]), 1):
    # out["VInt"].iloc[i] = (((out["TVDSS"].iloc[i] - out["TVDSS"].iloc[i-1]) /  ((out["TWT"].iloc[i] - out["TWT"].iloc[i-1]) / 2000)))

    return out, out["TWT"], out["VInt"]