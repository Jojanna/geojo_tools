import lasio
import pandas as pd
import numpy as np
from EEI_calc import eei_calc

# root for las file
root = r'F:\Namibia\Data\las'
wells = ["Murombe-1"]#["Shark-1"] #
file_suffix = "_merge_interp"
depth_column = "DEPTH"
null = -999.25

output_suffix = "_calc_logs"

vp_log = "VP"
vs_log = "VS"
rhob_log = "RHOB"

# for EEI calc
angles = [0, 90]
zoi_upper = None
zoi_lower = None



def calc_logs(wells, AI = True, K = True, MU = True, VPVS = True, LMR = True, PR = True, EEI = True, zoi_upper = None, zoi_lower = None):
    for well in wells:

        las_file = root + "\\" + well + "\\" + well + file_suffix + ".las"
        las = lasio.read(las_file)
        data_df = las.df()
        data_df.reset_index(inplace = True)

        logs_dict = {}
        for curve in las.curves:
            logs_dict[curve.mnemonic] = curve.unit

        if AI == True:
            data_df["AI"] = data_df[vp_log] * data_df[rhob_log]
            logs_dict["AI"] = "m/s*g/cc"

        if MU == True:
            data_df["MU"] = ((data_df[vs_log] ** 2) * (data_df[rhob_log] * 1000)) * pow(10, -9)
            logs_dict["MU"] = "GPa"

        if K == True:
            if "MU" not in data_df:
                print ("MU is necessary to calculate K from Vp!")
            else:
                data_df["K"] =  ((data_df[vp_log] ** 2) * (data_df[rhob_log] * 1000)) * pow(10, -9) - 4/3 * data_df["MU"]
                logs_dict["K"] = "GPa"

        if VPVS == True:
            data_df["VPVS"] = data_df[vp_log] / data_df[vs_log]
            logs_dict["VPVS"] = "ratio"

        if LMR == True:
            data_df["LAMBDA"] = data_df["K"] - 2/3 * data_df["MU"]
            data_df["LAMBDA-RHO"] = data_df["LAMBDA"] * data_df[rhob_log]
            data_df["MU-RHO"] = data_df["MU"] * data_df[rhob_log]
            logs_dict["LAMBDA"] = "GPa"
            logs_dict["LAMBDA-RHO"] = "GPa*g/cc"
            logs_dict["MU-RHO"] = "GPa*g/cc"

        if PR == True:
            data_df["POISSONS RATIO"] = (3 * data_df["K"] - 2 * data_df["MU"]) / (2 * (3 * data_df["K"] + data_df["MU"]))
            logs_dict["POISSONS RATIO"] = "ratio"

        if EEI == True:
            for eei_angle in angles:

                filt = data_df[["DEPTH", vp_log, vs_log, rhob_log]].dropna(how = "any")
                filt.reset_index(inplace = True)
                #print (filt)
                md = filt["DEPTH"].tolist()
                vp = filt[vp_log].tolist()
                vs = filt[vs_log].tolist()
                rhob = filt[rhob_log].tolist()
                rhob = list(map(lambda x: x * 1000, rhob))

                if zoi_upper == None:
                    zoi_upper = filt["DEPTH"].iloc[0]
                if zoi_lower == None:
                    zoi_lower = filt["DEPTH"].iloc[-1]

                eei = eei_calc(md, vp, vs, rhob, zoi_upper, zoi_lower, eei_angle)

                eei_df = pd.DataFrame()

                eei_df["EEI_" + str(eei_angle)] = eei
                eei_df["DEPTH"] = md

                logs_dict["EEI_" + str(eei_angle)] = "kg * m2/s"
                data_df = pd.merge(data_df, eei_df, on = "DEPTH", how = "outer")

        # write las
        output_file = root + "\\" + well + "\\" + well + output_suffix + ".las"

        with open(output_file, mode = 'w') as lasfile:
            las = lasio.LASFile()
            las.depth = ["DEPTH"]
            las.well["WELL"].value = str(well)
            las.well["NULL"].value = null
            for log in list(data_df.columns.values):
                las.add_curve(log, data_df[log], unit = logs_dict.get(log))
            las.write(lasfile, version=2, fmt = "%10.10g")



    return


calc_logs(wells, AI = True, K = True, MU = True, VPVS = True, LMR = True, PR = True, EEI = True, zoi_upper = None, zoi_lower = None)