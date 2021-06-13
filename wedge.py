from scipy.interpolate import interp1d
import pandas as pd
import numpy as np




dz_min = 0
dz_max = 100
dz_step = 10
def tvdss_to_twt(tvdss_in, df, tvdss = "TVDSS", twt = "TWT"):
    f = interp1d(df[tvdss], df[twt], kind = "cubic", fill_value = "extrapolate", assume_sorted=True)
    twt_out = f(tvdss_in)
    return twt_out

def twt_to_tvdss(twt_in, df, twt = "TWT", tvdss = "TVDSS"):
    g = interp1d(df[twt], df[tvdss], kind = "cubic", fill_value = "extrapolate", assume_sorted=True)
    tvdss_out = g(twt_in)
    return tvdss_out

def md_to_tvdss(md_in, md_log, tvdss_log):
    f = interp1d(md_log, tvdss_log, kind = "cubic", fill_value = "extrapolate", assume_sorted=True)
    out = f(md_in)
    return out

def md_to_twt(md_in, md_log, twt_log):
    f = interp1d(md_log, twt_log, kind = "cubic", fill_value = "extrapolate", assume_sorted=True)
    out = f(md_in)
    return out


def wedge(wedge_logs, tvdss = "TVDSS", mdkb = "MDKB", twt = "TWT", vp_log = "VP", vs_log = "VS", rhob_log = "RHOB", model_top = None, model_base = None, dz_min = 0, dz_max = 100, dz_step = 10):
    import sys

    thick0 = model_base - model_top
    nmodel = int((dz_max - dz_min) / dz_step + 1)
    thickness = np.linspace(dz_min, dz_max, nmodel)

    zmin0 = model_top
    zmax0 = model_base


    tmin0 = tvdss_to_twt(model_top, wedge_logs, tvdss = tvdss, twt = twt)
    tmax0 = tvdss_to_twt(model_base, wedge_logs, tvdss = tvdss, twt = twt)
    print ("TWT min/max: ")
    print (tmin0, tmax0)



    # create time-depth conversion for wedge. Uses mean velocity within wedge section

    vint_wedge = (zmax0 - zmin0) / (tmax0 - tmin0) * 2000

    wedge_logs["TWT_Wedge"] = np.nan
    wedge_logs["Vint"] = np.nan
    wedge_logs["Vint_Wedge"] = np.nan

    ## This Vint wedge bit is dubious and should maybe be avoided...

    wedge_logs["TWT_Wedge"].iloc[0] = wedge_logs["TWT"].iloc[0]
    print ("no. samples: %.0f" % len(wedge_logs[twt]))
    for i in range(1, len(wedge_logs[twt])):
        wedge_logs["Vint"].iloc[i] = (wedge_logs[tvdss].iloc[i] - wedge_logs[tvdss].iloc[i - 1]) / (
                wedge_logs[twt].iloc[i] - wedge_logs[twt].iloc[i - 1]) * 2000

        if wedge_logs[twt].iloc[i] <= tmin0:
            wedge_logs["TWT_Wedge"].iloc[i] = wedge_logs[twt].iloc[i]
        else:
            wedge_logs["TWT_Wedge"].iloc[i] = wedge_logs["TWT_Wedge"].iloc[i - 1] + (
                    wedge_logs[tvdss].iloc[i] - wedge_logs[tvdss].iloc[i - 1]) / vint_wedge * 2000
        #print("check3")
        wedge_logs["Vint_Wedge"].iloc[i] = (wedge_logs[tvdss].iloc[i] - wedge_logs[tvdss].iloc[i - 1]) / (
                wedge_logs["TWT_Wedge"].iloc[i] - wedge_logs["TWT_Wedge"].iloc[i - 1]) * 2000
    wedge_logs["Vint"].iloc[0] = wedge_logs["Vint"].iloc[1]
    wedge_logs["Vint_Wedge"].iloc[0] = wedge_logs["Vint_Wedge"].iloc[1]

    # create wedge logs

    # logs lists - retaining a list of column labels
    vp_list = []
    vs_list = []
    rhob_list = []

    i = 0
    for thick in thickness:

        if thick == 0:
            thick = thick + 0.1
        elif thick < 0:
            print("Warning! Negative Thickness!")
            sys.exit()

        #z0 = wedge_logs[tvdss]
        z_factor = thick / thick0
        #print (z_factor)
        # print (wedge_logs["TVDSS"].max())
        z = ((wedge_logs[tvdss] - zmin0) * z_factor + zmin0).tolist()
        dz = thick - thick0
        # print (thick)
        # print (max(z))

        vp = wedge_logs[vp_log].tolist()

        f = interp1d(z, vp, kind="nearest", fill_value="extrapolate", assume_sorted=True)
        wedge_logs["VP_%.0fm" % thick] = f(wedge_logs["TVDSS"].tolist())
        wedge_logs.loc[(wedge_logs[tvdss] < zmin0), "VP_%.0fm" % thick] = wedge_logs.loc[(wedge_logs[tvdss] < zmin0), vp_log]
        wedge_logs.loc[(wedge_logs[tvdss] >= (zmin0 + thick)), "VP_%.0fm" % thick] = np.nan

        wedge_list = wedge_logs["VP_%.0fm" % thick].dropna(how="any").tolist()
        wedge_list.extend(wedge_logs.loc[(wedge_logs[tvdss] >= (zmin0 + thick0)), vp_log].tolist())

        if len(wedge_logs.index) > len(wedge_list):
            tail = [np.nan] * (len(wedge_logs.index) - len(wedge_list))
            wedge_list.extend(tail)

        x = len(wedge_logs.index)
        wedge_logs["VP_%.0fm" % thick] = wedge_list[0:x]

        vp_list.append("VP_%.0fm" % thick)

        vs = wedge_logs[vs_log].tolist()
        f = interp1d(z, vs, kind="nearest", fill_value="extrapolate", assume_sorted=True)
        wedge_logs["VS_%.0fm" % thick] = f(wedge_logs["TVDSS"])
        wedge_logs.loc[(wedge_logs[tvdss] < zmin0), "VS_%.0fm" % thick] = wedge_logs.loc[
            wedge_logs[tvdss] < zmin0, vs_log]

        wedge_logs.loc[wedge_logs[tvdss] >= (zmin0 + thick), "VS_%.0fm" % thick] = np.nan

        wedge_list = wedge_logs["VS_%.0fm" % thick].dropna(how="any").tolist()
        wedge_list.extend(wedge_logs.loc[wedge_logs[tvdss] >= (zmin0 + thick0), vs_log].tolist())

        if len(wedge_logs.index) > len(wedge_list):
            tail = [np.nan] * (len(wedge_logs.index) - len(wedge_list))
            wedge_list.extend(tail)
        wedge_logs["VS_%.0fm" % thick] = wedge_list[0:x]

        vs_list.append("VS_%.0fm" % thick)

        rho = wedge_logs[rhob_log].tolist()

        f = interp1d(z, rho, kind="nearest", fill_value="extrapolate", assume_sorted=True)
        wedge_logs["RHOB_%.0fm" % thick] = f(wedge_logs["TVDSS"])
        wedge_logs.loc[wedge_logs[tvdss] < zmin0, "RHOB_%.0fm" % thick] = wedge_logs.loc[
            wedge_logs[tvdss] < zmin0, rhob_log]
        wedge_logs.loc[wedge_logs[tvdss] >=(zmin0 + thick), "RHOB_%.0fm" % thick] = np.nan
        wedge_list = wedge_logs["RHOB_%.0fm" % thick].dropna(how="any").tolist()
        # print (min(wedge_list))
        wedge_list.extend(wedge_logs.loc[wedge_logs[tvdss] >= (zmin0 + thick0), rhob_log].tolist())

        if len(wedge_logs.index) > len(wedge_list):
            tail = [np.nan] * (len(wedge_logs.index) - len(wedge_list))
            wedge_list.extend(tail)
        wedge_logs["RHOB_%.0fm" % thick] = wedge_list[0:x]

        rhob_list.append("RHOB_%.0fm" % thick)
        i = i + 1

    return wedge_logs


def wedge(wedge_logs, tvdss = "TVDSS", mdkb = "MDKB", twt = "TWT", vp_log = "VP", vs_log = "VS", rhob_log = "RHOB", model_top = None, model_base = None, dz_min = 0, dz_max = 100, dz_step = 10):
    import sys

    thick0 = model_base - model_top
    nmodel = int((dz_max - dz_min) / dz_step + 1)
    thickness = np.linspace(dz_min, dz_max, nmodel)

    zmin0 = model_top
    zmax0 = model_base


    tmin0 = tvdss_to_twt(model_top, wedge_logs, tvdss = tvdss, twt = twt)
    tmax0 = tvdss_to_twt(model_base, wedge_logs, tvdss = tvdss, twt = twt)
    print ("TWT min/max: ")
    print (tmin0, tmax0)



    # create time-depth conversion for wedge. Uses mean velocity within wedge section

    vint_wedge = (zmax0 - zmin0) / (tmax0 - tmin0) * 2000

    wedge_logs["TWT_Wedge"] = np.nan
    wedge_logs["Vint"] = np.nan
    wedge_logs["Vint_Wedge"] = np.nan

    ## This Vint wedge bit is dubious and should maybe be avoided...

    wedge_logs["TWT_Wedge"].iloc[0] = wedge_logs["TWT"].iloc[0]
    print ("no. samples: %.0f" % len(wedge_logs[twt]))
    for i in range(1, len(wedge_logs[twt])):
        wedge_logs["Vint"].iloc[i] = (wedge_logs[tvdss].iloc[i] - wedge_logs[tvdss].iloc[i - 1]) / (
                wedge_logs[twt].iloc[i] - wedge_logs[twt].iloc[i - 1]) * 2000

        if wedge_logs[twt].iloc[i] <= tmin0:
            wedge_logs["TWT_Wedge"].iloc[i] = wedge_logs[twt].iloc[i]
        else:
            wedge_logs["TWT_Wedge"].iloc[i] = wedge_logs["TWT_Wedge"].iloc[i - 1] + (
                    wedge_logs[tvdss].iloc[i] - wedge_logs[tvdss].iloc[i - 1]) / vint_wedge * 2000
        #print("check3")
        wedge_logs["Vint_Wedge"].iloc[i] = (wedge_logs[tvdss].iloc[i] - wedge_logs[tvdss].iloc[i - 1]) / (
                wedge_logs["TWT_Wedge"].iloc[i] - wedge_logs["TWT_Wedge"].iloc[i - 1]) * 2000
    wedge_logs["Vint"].iloc[0] = wedge_logs["Vint"].iloc[1]
    wedge_logs["Vint_Wedge"].iloc[0] = wedge_logs["Vint_Wedge"].iloc[1]

    # create wedge logs

    # logs lists - retaining a list of column labels
    vp_list = []
    vs_list = []
    rhob_list = []

    i = 0
    for thick in thickness:

        if thick == 0:
            thick = thick + 0.1
        elif thick < 0:
            print("Warning! Negative Thickness!")
            sys.exit()

        #z0 = wedge_logs[tvdss]
        z_factor = thick / thick0
        #print (z_factor)
        # print (wedge_logs["TVDSS"].max())
        z = ((wedge_logs[tvdss] - zmin0) * z_factor + zmin0).tolist()
        dz = thick - thick0
        # print (thick)
        # print (max(z))

        vp = wedge_logs[vp_log].tolist()

        f = interp1d(z, vp, kind="nearest", fill_value="extrapolate", assume_sorted=True)
        wedge_logs["VP_%.0fm" % thick] = f(wedge_logs["TVDSS"].tolist())
        wedge_logs.loc[(wedge_logs[tvdss] < zmin0), "VP_%.0fm" % thick] = wedge_logs.loc[(wedge_logs[tvdss] < zmin0), vp_log]
        wedge_logs.loc[(wedge_logs[tvdss] >= (zmin0 + thick)), "VP_%.0fm" % thick] = np.nan

        wedge_list = wedge_logs["VP_%.0fm" % thick].dropna(how="any").tolist()
        wedge_list.extend(wedge_logs.loc[(wedge_logs[tvdss] >= zmax0), vp_log].tolist()) # (zmin0 + thick0)

        if len(wedge_logs.index) > len(wedge_list):
            tail = [np.nan] * (len(wedge_logs.index) - len(wedge_list))
            wedge_list.extend(tail)

        x = len(wedge_logs.index)
        wedge_logs["VP_%.0fm" % thick] = wedge_list[0:x]

        vp_list.append("VP_%.0fm" % thick)

        vs = wedge_logs[vs_log].tolist()
        f = interp1d(z, vs, kind="nearest", fill_value="extrapolate", assume_sorted=True)
        wedge_logs["VS_%.0fm" % thick] = f(wedge_logs["TVDSS"])
        wedge_logs.loc[(wedge_logs[tvdss] < zmin0), "VS_%.0fm" % thick] = wedge_logs.loc[
            wedge_logs[tvdss] < zmin0, vs_log]

        wedge_logs.loc[wedge_logs[tvdss] >= (zmin0 + thick), "VS_%.0fm" % thick] = np.nan

        wedge_list = wedge_logs["VS_%.0fm" % thick].dropna(how="any").tolist()
        wedge_list.extend(wedge_logs.loc[wedge_logs[tvdss] >= zmax0, vs_log].tolist()) #(zmin0 + thick0)

        if len(wedge_logs.index) > len(wedge_list):
            tail = [np.nan] * (len(wedge_logs.index) - len(wedge_list))
            wedge_list.extend(tail)
        wedge_logs["VS_%.0fm" % thick] = wedge_list[0:x]

        vs_list.append("VS_%.0fm" % thick)

        rho = wedge_logs[rhob_log].tolist()

        f = interp1d(z, rho, kind="nearest", fill_value="extrapolate", assume_sorted=True)
        wedge_logs["RHOB_%.0fm" % thick] = f(wedge_logs["TVDSS"])
        wedge_logs.loc[wedge_logs[tvdss] < zmin0, "RHOB_%.0fm" % thick] = wedge_logs.loc[
            wedge_logs[tvdss] < zmin0, rhob_log]
        wedge_logs.loc[wedge_logs[tvdss] >=(zmin0 + thick), "RHOB_%.0fm" % thick] = np.nan
        wedge_list = wedge_logs["RHOB_%.0fm" % thick].dropna(how="any").tolist()
        # print (min(wedge_list))
        wedge_list.extend(wedge_logs.loc[wedge_logs[tvdss] >= zmax0, rhob_log].tolist()) #(zmin0 + thick0)

        if len(wedge_logs.index) > len(wedge_list):
            tail = [np.nan] * (len(wedge_logs.index) - len(wedge_list))
            wedge_list.extend(tail)
        wedge_logs["RHOB_%.0fm" % thick] = wedge_list[0:x]

        rhob_list.append("RHOB_%.0fm" % thick)
        i = i + 1

    return wedge_logs

def wedge_build_in_twt(wedge_logs, tvdss = "TVDSS", mdkb = "MDKB", twt = "TWT", vp_log = "VP", vs_log = "VS", rhob_log = "RHOB", model_top = None, model_base = None, dz_min = 0, dz_max = 100, dz_step = 10):

    # This isn't working
    import sys

    thick0 = model_base - model_top
    nmodel = int((dz_max - dz_min) / dz_step + 1)
    thickness = np.linspace(dz_min, dz_max, nmodel)

    zmin0 = model_top
    zmax0 = model_base


    tmin0 = tvdss_to_twt(model_top, wedge_logs, tvdss = tvdss, twt = twt)
    tmax0 = tvdss_to_twt(model_base, wedge_logs, tvdss = tvdss, twt = twt)
    print ("TWT min/max: ")
    print (tmin0, tmax0)



    # create time-depth conversion for wedge. Uses mean velocity within wedge section

    vint_wedge = (zmax0 - zmin0) / (tmax0 - tmin0) * 2000

    wedge_logs["TWT_Wedge"] = np.nan
    wedge_logs["Vint"] = np.nan
    wedge_logs["Vint_Wedge"] = np.nan

    ## This Vint wedge bit is dubious and should maybe be avoided...

    wedge_logs["TWT_Wedge"].iloc[0] = wedge_logs["TWT"].iloc[0]
    print ("no. samples: %.0f" % len(wedge_logs[twt]))
    for i in range(1, len(wedge_logs[twt])):
        wedge_logs["Vint"].iloc[i] = (wedge_logs[tvdss].iloc[i] - wedge_logs[tvdss].iloc[i - 1]) / (
                wedge_logs[twt].iloc[i] - wedge_logs[twt].iloc[i - 1]) * 2000

    # create wedge logs

    # logs lists - retaining a list of column labels
    vp_list = []
    vs_list = []
    rhob_list = []

    i = 0
    for thick in thickness:

        if thick == 0:
            thick = thick + 0.1
        elif thick < 0:
            print("Warning! Negative Thickness!")
            sys.exit()

        #z0 = wedge_logs[tvdss]
        z_factor = thick / thick0
        #print (z_factor)
        # print (wedge_logs["TVDSS"].max())
        z = ((wedge_logs[tvdss] - zmin0) * z_factor + zmin0).tolist()
        dz = thick - thick0
        # print (thick)
        # print (max(z))
        thick_twt = (thick/vint_wedge)*2000
        #print (thick_twt)
        #t = ((wedge_logs[twt] - tmin0) * z_factor + tmin0).tolist()
        t = z/vint_wedge * 2000

        vint_log = wedge_logs["Vint_Wedge"].tolist()

        for log, log_str, log_list in zip([vp_log, vs_log, rhob_log],["VP", "VS", "RHOB"], [vp_list, vs_list, rhob_list]):

            log_ini = wedge_logs[log_str].tolist()

            f = interp1d(t, log_ini, kind="nearest", fill_value="extrapolate", assume_sorted=True)
            wedge_logs[log_str + "_%.0fm" % thick] = f(wedge_logs["TWT"].tolist())
            wedge_logs.loc[(wedge_logs[twt] < tmin0), log_str + "_%.0fm" % thick] = wedge_logs.loc[(wedge_logs[twt] < tmin0), log_str]
            wedge_logs.loc[(wedge_logs[twt] >= (tmin0 + thick_twt)), log_str + "_%.0fm" % thick] = np.nan

            wedge_list = wedge_logs[log_str + "_%.0fm" % thick].dropna(how="any").tolist()
            wedge_list.extend(wedge_logs.loc[(wedge_logs[twt] >= (tmax0)), log_str].tolist()) #??? tmax0 vs (zmin0 + thick0)??

            if len(wedge_logs.index) > len(wedge_list):
                tail = [np.nan] * (len(wedge_logs.index) - len(wedge_list))
                wedge_list.extend(tail)

            x = len(wedge_logs.index)
            wedge_logs[log_str + "_%.0fm" % thick] = wedge_list[0:x]

            log_list.append(log_str + "_%.0fm" % thick)

        i = i + 1

    return wedge_logs


def wedge_vint(wedge_logs, tvdss="TVDSS", mdkb="MDKB", twt="TWT", vp_log="VP", vs_log="VS", rhob_log="RHOB",
                       model_top=None, model_base=None, dz_min=0, dz_max=100, dz_step=10):
    #make sure input logs are regularly sampled in TVDSS

    #test - checking for rgular sampling
    wedge_logs["dz"] = np.nan
    for i in range(1, len(wedge_logs[tvdss])):
        z1 = wedge_logs[tvdss].iloc[i - 1]
        z2 = wedge_logs[tvdss].iloc[i]
        wedge_logs["dz"].iloc[i] = (z2 - z1)

    wedge_logs["dz"].iloc[0] = wedge_logs["dz"].iloc[1]

    for i in range(1, len(wedge_logs[tvdss])):
        #print(wedge_logs["dz"].iloc[i].round(3), np.round(wedge_logs["dz"].mean(), 3))
        if wedge_logs["dz"].iloc[i].round(3) != np.round(wedge_logs["dz"].mean(), 3):
            import sys
            print ("irregular TVDSS sampling!")
            print (wedge_logs["dz"].iloc[i], wedge_logs["dz"].mean())
            sys.exit()


    import sys

    thick0 = model_base - model_top
    nmodel = int((dz_max - dz_min) / dz_step + 1)
    thickness = np.linspace(dz_min, dz_max, nmodel)

    zmin0 = model_top
    zmax0 = model_base

    tmin0 = tvdss_to_twt(model_top, wedge_logs, tvdss=tvdss, twt=twt)
    tmax0 = tvdss_to_twt(model_base, wedge_logs, tvdss=tvdss, twt=twt)
    print("TWT min/max: ")
    print(tmin0, tmax0)


    wedge_logs["Vint"] = np.nan

    # calculate vint log, stretch/squeeze with vp/vs/rhob logs, then calculate twt at each thickness with this
    print("no. samples: %.0f" % len(wedge_logs[twt]))
    for i in range(1, len(wedge_logs[tvdss])):
        z1 = wedge_logs[tvdss].iloc[i - 1]
        z2 = wedge_logs[tvdss].iloc[i]
        t2 = wedge_logs[twt].iloc[i]
        t1 = wedge_logs[twt].iloc[i - 1]
        wedge_logs["Vint"].iloc[i] = (z2 - z1) / ((t2 - t1)/2000)

    wedge_logs["Vint"].iloc[0] = wedge_logs["Vint"].iloc[1]

    # create wedge logs
    #vint_log = wedge_logs["Vint"].tolist()
    vint_log = "Vint"
    # logs lists - retaining a list of column labels
    vp_list = []
    vs_list = []
    rhob_list = []
    vint_list = []


    i = 0
    for thick in thickness:

        if thick == 0:
            thick = thick + 0.1
        elif thick < 0:
            print("Warning! Negative Thickness!")
            sys.exit()

        z_factor = thick / thick0
        # this log represents the stretched log
        z = ((wedge_logs[tvdss] - zmin0) * z_factor + zmin0).tolist()
        #df_stretch = wedge_logs[[tvdss, vint_log, vp_log, vs_log, rhob_log]].copy()
        #df_stretch["TVDSS_ss"] = ((df_stretch[tvdss] - zmin0) * z_factor + zmin0).tolist()
        #df_stretch = df_stretch.loc[(df_stretch[tvdss]> float(df_stretch["TVDSS_ss"].min())) & (df_stretch[tvdss] < float(df_stretch["TVDSS_ss"].max()))].dropna()

        #print(wedge_logs[tvdss].min(), wedge_logs[tvdss].max())
        #print(df_stretch["TVDSS_ss"].min(), df_stretch["TVDSS_ss"].max())
        #print(df_stretch[tvdss].min(), df_stretch[tvdss].max())

        #z = df_stretch["TVDSS_ss"].tolist()

        for log, log_str, log_list in zip([vint_log, vp_log, vs_log, rhob_log], ["Vint", "VP", "VS", "RHOB"],
                                          [vint_list, vp_list, vs_list, rhob_list]):

            log_ini = wedge_logs[log_str].tolist()

            if log_str == "Vint":
                kind = "linear"
            else:
                kind = "nearest"

            fill_value = (log_ini[0], log_ini[-1])
            #print ("fill value")
            #print (fill_value)
            f = interp1d(z, log_ini, kind=kind, fill_value="extrapolate", assume_sorted=True)
            wedge_logs[log_str + "_%.0fm" % thick] = f(wedge_logs["TVDSS"].tolist())
            wedge_logs.loc[(wedge_logs[tvdss] < zmin0), log_str + "_%.0fm" % thick] = wedge_logs.loc[
                (wedge_logs[tvdss] < zmin0), log_str]
            wedge_logs.loc[(wedge_logs[tvdss] >= (zmin0 + thick)), log_str + "_%.0fm" % thick] = np.nan

            wedge_list = wedge_logs[log_str + "_%.0fm" % thick].dropna(how="any").tolist()
            wedge_list.extend(wedge_logs.loc[(wedge_logs[tvdss] >= zmax0), log_str].tolist())

            if len(wedge_logs.index) > len(wedge_list):
                tail = [log_ini[-1]] * (len(wedge_logs.index) - len(wedge_list))
                wedge_list.extend(tail)

            x = len(wedge_logs.index)
            wedge_logs[log_str + "_%.0fm" % thick] = wedge_list[0:x]

            log_list.append(log_str + "_%.0fm" % thick)

        wedge_logs["TWT" + "_%.0fm" % thick] = np.nan

        wedge_logs.loc[wedge_logs["TVDSS"] < zmin0, "TWT" + "_%.0fm" % thick] = wedge_logs.loc[
            wedge_logs["TVDSS"] < zmin0, "TWT"]
        for i in range(1, len(wedge_logs.index)):
            t1 = wedge_logs["TWT" + "_%.0fm" % thick].iloc[i-1]
            z1 = wedge_logs["TVDSS"].iloc[i-1]
            z2 = wedge_logs["TVDSS"].iloc[i]
            vint =  wedge_logs["Vint" + "_%.0fm" % thick].iloc[i]

            wedge_logs["TWT" + "_%.0fm" % thick].iloc[i] = (z2-z1)/vint * 2000 + t1

        i = i + 1

    return wedge_logs


def zone_tvdss(md_log, tvdss_log, zone):
    f = interp1d(md_log, tvdss_log, kind = "cubic", fill_value = "extrapolate", assume_sorted=True)
    out = [f(zone[0]), f(zone[1]), zone[2], zone[3]]
    return out

def zone_twt(md_log, twt_log, zone):
    f = interp1d(md_log, twt_log, kind = "cubic", fill_value = "extrapolate", assume_sorted=True)
    out = [f(zone[0]), f(zone[1]), zone[2], zone[3]]
    return out

def plot_wedge_z(wedge_logs, tvdss = "TVDSS", mdkb = "MDKB", vp = "VP", vs = "VS", rhob = "RHOB", model_top = None, model_base = None, dz_min = 0, dz_max = 100, dz_step = 10, tops = None, fill = True):
    # where tops are defined zone = [upper, lower, name (str), colour]
    nmodel = int((dz_max - dz_min) / dz_step + 1)
    thickness = np.linspace(dz_min, dz_max, nmodel)
    thick0 = model_base - model_top

    zmin0 = model_top
    zmax0 = model_base

    fig, ax = plt.subplots(1, nmodel, sharey=True, figsize=(16, 10))
    i = 0
    for thick in thickness:
        # ax[i].set_facecolor('tan')
        # ax[0][i].plot(wedge_logs["VP_Blocky"], wedge_logs["TVDSS"], label = "")
        ax[i].plot(wedge_logs["VP_%.0fm" % thick], wedge_logs[tvdss], label="VP, m/s")
        ax[i].plot(wedge_logs["VS_%.0fm" % thick], wedge_logs[tvdss], label="VS, m/s")
        ax[i].plot(wedge_logs["RHOB_%.0fm" % thick] * 1000, wedge_logs[tvdss], label="RHOB, kg/cm3")
        ax[i].set_ylim(zmin0 - 100, zmax0 + 200)
        # ax[i].set_ylim(1900, 2000)
        ax[i].grid()
        ax[i].axhline(zmin0, c="k")
        ax[i].axhline(zmin0 + thick, c="k")

        ax[i].set_xlim(0, 4000)
        ax[i].set_title("Thickness = %.0fm" % thick, fontsize="small")
        ax[i].tick_params(labelsize="small")

        if tops != None:
            z_factor = thick / thick0
            for zone in tops:
                z = zone_tvdss(wedge_logs[mdkb], wedge_logs[tvdss], zone)
                z1 = (z[0] - zmin0) * z_factor + zmin0
                z2 = (z[1] - zmin0) * z_factor + zmin0
                ax[i].axhline(z1, c="k")
                ax[i].axhline(z2, c="k")
                if fill == True:
                    ax[i].fill_between([0, 4000], z1, z2, color=z[3])

                if i == nmodel - 1:
                    ax[nmodel - 1].text(4000, np.mean([z1, z2]), z[2], color="r")
                # print (np.mean((z1, z2), z[2])

        i = i + 1
    #fig.suptitle("Tuning Wedge Logs (Anchois-ST1, B Sand)")
    # ax1.text(0, np.mean([zone[0], zone[1]]), str(zone[2]), color = "r")
    ax[0].invert_yaxis()
    ax[0].set_ylabel("TVDSS, m")
    ax[nmodel - 1].legend()

    fig.subplots_adjust(wspace=0)
    plt.tight_layout()

    return fig

def plot_wedge_properties_mesh(wedge_logs, z = "TVDSS", dz_min = 0, dz_max = 100, dz_step = 10, yrange = [1900, 2100]):

    nmodel = int((dz_max - dz_min) / dz_step + 1)
    thickness = np.linspace(dz_min, dz_max, nmodel)

    vp_list = []
    vs_list = []
    rhob_list = []

    for thick in thickness:
        vp_list.append("VP_%.0fm" % thick)
        vs_list.append("VS_%.0fm" % thick)
        rhob_list.append("RHOB_%.0fm" % thick)

    vp_vol_matrix = wedge_logs[vp_list].to_numpy()
    vs_vol_matrix = wedge_logs[vs_list].to_numpy()
    rhob_vol_matrix = wedge_logs[rhob_list].to_numpy()

    ai_vol_matrix = vp_vol_matrix * rhob_vol_matrix
    vpvs_vol_matrix = vp_vol_matrix / vs_vol_matrix

    fig, ((ax1, ax2, ax3), (ax4, ax5, ax6)) = plt.subplots(2, 3, figsize=(30, 20))
    fig.set_size_inches(24, 12)

    cb = ax1.pcolormesh(thickness, wedge_logs[z].tolist(), vp_vol_matrix, cmap="cool")
    ax1.set_title("Vp, m/s\n", )

    cb2 = ax2.pcolormesh(thickness, wedge_logs[z].tolist(), vs_vol_matrix, cmap="cool")
    ax2.set_xticks(np.arange(dz_min, dz_max, dz_step))
    ax2.set_title("Vs, m/s\n")
    cb3 = ax3.pcolormesh(thickness, wedge_logs[z].tolist(), rhob_vol_matrix, cmap="cool")
    ax3.set_xticks(np.arange(dz_min, dz_max, dz_step))
    ax3.set_title("RhoB, g/cc\n")
    cb4 = ax4.pcolormesh(thickness, wedge_logs[z].tolist(), ai_vol_matrix, cmap="cool")
    ax4.set_title("P-Impedance, m/s * g/cc\n")
    ax4.set_xticks(np.arange(dz_min, dz_max, dz_step))
    cb5 = ax5.pcolormesh(thickness, wedge_logs[z].tolist(), vpvs_vol_matrix, cmap="cool")
    ax5.set_xticks(np.arange(dz_min, dz_max, dz_step))
    ax5.set_title("Vp/Vs, ratio\n")
    ax6.remove()

    # text_box = TextBox(ax6, 'Evaluate', initial="initial_text")

    # a=ax1.get_yticks().tolist()
    # print (a)
    # a = [b * dt for b in a]

    plt.colorbar(cb, ax=ax1, label="Vp, m/s")
    plt.colorbar(cb2, ax=ax2, label="Vs, m/s")
    plt.colorbar(cb3, ax=ax3, label="RhoB, g/cc")
    plt.colorbar(cb4, ax=ax4, label="AI, m/s * g/cc")
    plt.colorbar(cb5, ax=ax5, label="Vp/Vs, ratio")

    for ax in (ax1, ax2, ax3, ax4, ax5):
        ax.set_ylim(yrange[0], yrange[1])
        ax.set_xlabel("Thickness, m")
        ax.set_ylabel(str(z))
        ax.set_xticks(np.arange(dz_min, dz_max, dz_step))
        ax.xaxis.set_tick_params(labeltop=True, labelbottom=False, bottom=False, top=True)

        ax.invert_yaxis()

    plt.tight_layout()
    return fig

def convert_resamp_wedge_twt(wedge_logs, twt = "TWT", tmin = 0, tmax = 1000, dt = 0.5, logs_list = None):

    if logs_list == None:
        logs_list = wedge_logs.columns

    wedge_logs_twt = pd.DataFrame()
    wedge_logs_twt[twt] = np.linspace(tmin, tmax, (tmax - tmin + dt) / dt)
    #print((tmax - tmin + dt) / dt)

    for col in logs_list:
        f = interp1d(wedge_logs["TWT_Wedge"], wedge_logs[col], kind="nearest", fill_value="extrapolate",
                     bounds_error=False)
        wedge_logs_twt[col] = f(wedge_logs_twt["TWT"])

    return wedge_logs_twt

def zoep_wedge(wedge_logs_twt, angles, dz_min = 0, dz_max = 100, dz_step = 10):

    from .zoep_synth import rc_zoep_trace
    nmodel = int((dz_max - dz_min) / dz_step + 1)
    thickness = np.linspace(dz_min, dz_max, nmodel)

    vp_list = []
    vs_list = []
    rhob_list = []

    for thick in thickness:
        vp_list.append("VP_%.0fm" % thick)
        vs_list.append("VS_%.0fm" % thick)
        rhob_list.append("RHOB_%.0fm" % thick)

    vol_refs = np.ndarray(shape=(1, nmodel, len(angles), len(wedge_logs_twt.index)))

    for vp_log, vs_log, rhob_log, n in zip(vp_list, vs_list, rhob_list, range(0, nmodel)):
        print(vp_log, vs_log, rhob_log)
        vp = wedge_logs_twt[vp_log].tolist()
        vs = wedge_logs_twt[vs_log].tolist()
        rhob = wedge_logs_twt[rhob_log].tolist()
        refs = rc_zoep_trace(vp, vs, rhob, angles)

        vol_refs[0, n, :, :] = refs


def convolve_wedge_volume(vol, nmodel, nangles, wvlt_amp, dt):
    # print (nmodel, nangle)
    sgy = np.ndarray(shape=(1, len(vol[0, :, 0, 0]), len(vol[0, 0, :, 0]), len(vol[0, 0, 0, :])))
    # print (np.shape(sgy))

    for i in range(0, nmodel):
        # print ("Model no. = %i" % i)

        for j in range(0, nangles):
            # print ("Angle no. = %i" % (j))
            rc = vol[0, i, j, :]
            syn = np.convolve(rc, wvlt_amp, mode='same')
            sgy[0, i, j, :] = syn

    return sgy


def resample_wedge_volume(sgy, nmodel, nangle, dt_in, dt_out):
    from scipy.signal import resample
    nsamp_in = len(sgy[0, 0, 0, :])
    print(nsamp_in)
    nsamp_out = int(nsamp_in / (dt_out / dt_in))
    print(nsamp_out)
    print(dt_out / dt_in)
    sgy_out = np.ndarray(shape=(1, len(sgy[0, :, 0, 0]), len(sgy[0, 0, :, 0]), nsamp_out))
    for i in range(nmodel):
        for j in range(nangle):
            syn_resamp = signal.resample(sgy[0, i, j, :], num=nsamp_out)
            sgy_out[0, i, j, :] = syn_resamp

    return sgy_out

def plot_traces_wedge(sgy, wedge_logs_twt, dz_min = 0, dz_max = 100, dz_step = 10, mdkb = "MDKB", twt = "TWT", model_top = None, model_base = None, markers = None, plot_markers = True):

    nmodel = int((dz_max - dz_min) / dz_step + 1)
    thickness = np.linspace(dz_min, dz_max, nmodel)
    thick0 = model_base - model_top
    tmin0 = tvdss_to_twt(model_top, wedge_logs_twt, tvdss=tvdss, twt=twt)

    fig, ax = plt.subplots(1, nmodel, figsize=(16, 10), sharey=True)


    for n, thick in zip(range(0, nmodel), thickness):

        ax[n].plot(sgy[0, n, 0, :], wedge_logs_twt[twt], c="k", linewidth=0.5)

        cmap = plt.get_cmap(seismic_petrel)
        norm = plt.Normalize(vmin=-10.0, vmax=10.0)
        # print (n)

        for j in range(len(wedge_logs_twt["TWT"]) - 1):
            a = float(wedge_logs_twt["TWT"].tolist()[j])
            b = float(wedge_logs_twt["TWT"].tolist()[j + 1])
            c = float(sgy[0, n, 0, j])
            d = float(sgy[0, n, 0, j + 1])
            ax[n].fill_betweenx([a, b], [c, d], color=cmap(norm(c)))  #

            z_factor = thick / thick0
        if markers == True:
            for zone in (markers):
                z = zone_twt(wedge_logs_twt[mdkb], wedge_logs_twt[twt], zone)
                z1 = (z[0] - tmin0) * z_factor + tmin0
                z2 = (z[1] - tmin0) * z_factor + tmin0
                ax[n].axhline(z1, c="r")
                ax[n].axhline(z2, c="r")

                if n == nmodel - 1:
                    ax[nmodel - 1].text(0.1, np.mean([z1, z2]), z[2], color="r")

        ax[n].grid()
        # ax[n].set_xlim(-0.2, 0.2)
        ax[n].set_title("Thickness = %.0fm\n" % thick, fontsize="small")
        ax[n].tick_params(labelsize="small")
        ax[n].set_xlim(-10, 10)
    ax[0].set_ylim(1900, 2100)
    ax[0].set_ylabel("TWT, ms")
    ax[0].invert_yaxis()
    fig.subplots_adjust(wspace=0)

    return fig