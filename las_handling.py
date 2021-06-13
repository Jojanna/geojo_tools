import lasio
import pandas as pd
import numpy as np

def load_las(root, file, null):
    path = root + "\\" + file + ".las"

    data = lasio.read(path)

    units_dict = {}
    logs = []
    units = []

    for curve in data.curves:
        logs.append(curve.mnemonic)
        units.append(curve.unit)

        units_dict[curve.mnemonic] = curve.unit

    logs_df = pd.DataFrame()

    for log in logs:
        logs_df[str(log)] = data[str(log)]
        logs_df[str(log)] = np.where(logs_df[str(log)] == null, np.nan, logs_df[str(log)])
        #logs_df[str(log)] = logs_df[str(log)].interpolate(method="nearest")
    logs_df.set_index(logs[0], inplace = True)

    return logs_df, units_dict


def write_las(las_path, file, df, logs_list, depth_log = "DEPTH", units_dict = None):

    df = df[logs_list].dropna(how="all")

    output = las_path + "\\" + file + ".las"
    with open(output, mode='w') as lasfile:
        las = lasio.LASFile()
        las.well["NULL"].value = -999.25
        las.depth = depth_log
        las.well["WELL"].value = str(file)
        las.add_curve(depth_log, df.index)

        for log in logs_list:
            if units_dict == None:
                unit = None
            else:
                if log in units_dict:
                    unit = units_dict[log]
                else:
                    unit = None
            # print (unit)
            las.add_curve(str(log), df[log].tolist(), unit=unit)

        las.write(lasfile, version=2, fmt='%10.9g')

    return output