import lasio
import pandas as pd
import numpy as np

# root for main las file
root = r'F:\Namibia\Data\las'
wells = ["Wingat-1"]#["Shark-1"] #
file_name = "Wingat-1_KINGDOM"
depth_column = "DEPTH" #this will be the column that the files are merged on. These files should be sampled at the same intervals.

# root for las file to be merged with main las
merge_root = root
merge_name = "Wingat-1_comp"
merge_depth_col = "DEPT" #this will be the column that the files are merged on. These files should be sampled at the same intervals.

#root for output las file
out_root = root
out_suffix = "_merge"

# Null value
null = -999.25

rename = {} #{"Rho_FINAL": "RHOB", "DT_FINAL": "DTCO", "DTS_FINAL": "DTSM",  }



def merge_las(root, wells, file_name, depth_column, merge_root, merge_name, merge_depth_col, out_root, null):
    remove = False
    remove_logs = []#["TVD", "TVDSS", "X_Coord", "Y_Coord", "TWT", "HTEN", "PEF"]

    select = True
    select_logs_main = None #[]# [depth_column, "AI", "CALI", "DTCO", "GR", "NPHI", "PR", "RES_D", "RES_S", "RHOB", "SW", "SO", "VClay", "VQuartz", "VCalcite", "VP", "VS"]
    select_logs_merge = ["HDRA"] #["AT90", "AT10", "BADHOLE", "N_P"]

    filepaths = []
    output_filenames = []

    for well in wells:

        # load merge file

        logs_dict = {}
        merge_path = merge_root + "\\" + well + "\\" + merge_name + ".las"
        print (merge_path)

        data = lasio.read(merge_path)
        logs_merge = []
        units_merge = []
        for curve in data.curves:
            logs_merge.append(curve.mnemonic)
            units_merge.append(curve.unit)
            if curve.unit != "UNKNOWN":
                logs_dict[curve.mnemonic] = curve.unit

        if select == True:
            if select_logs_merge != None:
                logs_merge = select_logs_merge
                logs_merge.append(merge_depth_col)

        merge_df = pd.DataFrame()


        for log in logs_merge:
            merge_df[str(log)] = data[str(log)]
            merge_df[str(log)] = np.where(merge_df[str(log)] == null, np.nan, merge_df[str(log)])

        merge_df.rename(columns = {merge_depth_col: depth_column}, inplace=True)

        # load main las file

        path = root + "\\" + well + "\\" + file_name + ".las"
        print (path)
        filepaths.append(path)

        data_main = lasio.read(path)
        logs_main = []
        units_main = []
        for curve in data_main.curves:
            logs_main.append(curve.mnemonic)
            units_main.append(curve.unit)
            if curve.unit != "UNKNOWN":
                logs_dict[curve.mnemonic] = curve.unit

            # to select logs
            if select == True:
                if select_logs_main != None:
                    logs_main = select_logs_main
                    logs_main.append(depth_column)

            main_df = pd.DataFrame()

            for log in logs_main:
                main_df[str(log)] = data_main[str(log)]
                main_df[str(log)] = np.where(main_df[str(log)] == null, np.nan, main_df[str(log)])

        # merge and print
        output = pd.merge(main_df, merge_df, on = depth_column, how = "outer")

        if remove == True:
            for rem in remove_logs:
                print (list(output.columns.values))
                print (output.head())
                if rem in list(output.columns.values):
                    output.drop(rem, axis = 1, inplace = True)

        output.dropna(inplace = True, how = "all")

        output_file = out_root + "\\" + well + "\\" + well + out_suffix + ".las"
        output_filenames.append(output_file)

        output.rename(columns={depth_column: "DEPTH"}, inplace=True)
        output.set_index("DEPTH", inplace = True)
        output.dropna(how = "all", axis = 0, inplace = True)#, subset = ss)
        output.sort_index(axis = 1, ascending = True, inplace = True)


        with open(output_file, mode = 'w') as lasfile:
            las = lasio.LASFile()
            las.depth = ["DEPTH"]
            las.well["WELL"].value = str(well)
            las.well["NULL"].value = null
            las.add_curve("DEPTH", output.index.values, unit = logs_dict.get(depth_column))
            for log in list(output.columns.values):
                las.add_curve(log, output[log], unit = logs_dict.get(log))
            las.write(lasfile, version=2, fmt = "%10.10g")

    return output_filenames

merge_las(root, wells, file_name, depth_column, merge_root, merge_name, merge_depth_col, out_root, null)







