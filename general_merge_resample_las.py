import lasio
import pandas as pd
import numpy as np
from scipy.interpolate import interp1d
from units_dict import units_dict
import sys
from load_well_table import read_well_table



# root for main las file
root = r'F:\Namibia\Data\las'
well = "Shark-1" #
file_name = "Shark-1_merge"
depth_column = "DEPTH" #this will be the column that the files are merged on. These files should be sampled at the same intervals.
well_table = read_well_table(r'F:\Namibia\well_database.xlsx', well)
print (well_table)

# root for las file to be merged with main las
merge = False
merge_root = root
merge_name = "2012_13-1_9"
merge_depth_col = "DEPTH" #this will be the column that the files are merged on. These files should be sampled at the same intervals.

#root for output las file
out_root = root
out_suffix = "_merge_merge_edit"

#output las interval
#output_int = 0.1524
output_int = 0.1524

# datum info. Required to create vert TVDml
kb = well_table.get("kb") # [28.05]
wd = well_table.get("wd") #[1387.5]
easting = well_table.get("coord_e")
northing = well_table.get("coord_n")

# set to none to use maximum range
start_depth = None
end_depth = None

create_vert_TVDml = True
# functionality for non vertical wells not yet created
vertical = True
# Null value
null = -999.25

# set to none where not required
rename = {"Rho_FINAL": "RHOB", "DT_FINAL": "DTCO", "DTS_FINAL": "DTSM", "VCLAY": "VSH"}

# set to none where not required
set_units = {"PHIT": "frac", "SW": "frac", "VP": "m/s", "VS": "m/s", "VCL": "frac", "PHIE": "frac", "RHOB": "g/cc"}

remove = False
remove_logs = ["TVD", "TVDSS", "X_Coord", "Y_Coord", "TWT", "HTEN", "PEF"]

# set to none where not required
select = True
select_logs_main = [depth_column, "AI", "CALI_x", "DTCO_x", "GR_x", "NPHI", "PR", "RES_D", "RES_S", "RHOB_x", "SW_y", "SO", "VClay", "VQuartz", "VCalcite", "VP", "VS", "PHIE", "PHIT", "VCL", "RT", "N_P", "BADHOLE", "DTSM", "HLLD", "HLLS"]
select_logs_merge = ["AT90", "AT10", "BADHOLE", "N_P"]

#interpolate = False # not yet built in

def merge_las(root, well, kb, wd, create_vert_TVDml, file_name, depth_column, merge_root, merge_name, merge_depth_col, out_root, null, start_depth, end_depth, rename):

    filepaths = []
    output_filenames = []



    # create list/dictionary of logs and units
    logs_dict = units_dict

    if merge == True:
        # load merge file
        merge_path = merge_root + "\\" + well + "\\" + merge_name + ".las"
        data = lasio.read(merge_path)

        logs_merge = []
        units_merge = []
        for curve in data.curves:
            logs_merge.append(curve.mnemonic)
            units_merge.append(curve.unit)
            if curve.unit != "UNKNOWN":
                if curve.unit != "":
                    logs_dict[curve.mnemonic] = curve.unit

        # select desired logs for output file from file to be merged in
        if select == True:
            if select_logs_merge != None:
                logs_merge = select_logs_merge
                logs_merge.append(merge_depth_col)

        merge_df = pd.DataFrame()


        # write nan to null values in file to be merged in
        for log in logs_merge:
            merge_df[str(log)] = data[str(log)]
            merge_df[str(log)] = np.where(merge_df[str(log)] == null, np.nan, merge_df[str(log)])

        # rename depth column to match main las file to ensure merge
        merge_df.rename(columns={merge_depth_col: depth_column}, inplace=True)

    # load main las file
    path = root + "\\" + well + "\\" + file_name + ".las"
    #filepaths.append(path)

    # create list/dictionary of logs and units
    data_main = lasio.read(path)
    logs_main = []
    units_main = []
    for curve in data_main.curves:
        logs_main.append(curve.mnemonic)
        units_main.append(curve.unit)
        if curve.unit != "UNKNOWN":
            if curve.unit != "":
                logs_dict[curve.mnemonic] = curve.unit

    if set_units != None:
        logs_dict.update(set_units)

    # to select logs
    if select == True:
        if select_logs_main != None:
            logs_main = select_logs_main
            logs_main.append(depth_column)

    # write nan to null values in main las file
    main_df = pd.DataFrame()

    for log in logs_main:
        main_df[str(log)] = data_main[str(log)]
        main_df[str(log)] = np.where(main_df[str(log)] == null, np.nan, main_df[str(log)])

    # unless specified, min/max depth extracted from range of main las file
    if start_depth == None:
        start_depth = main_df[depth_column].iloc[0]
    if end_depth == None:
        end_depth = main_df[depth_column].iloc[-1]

    if merge == True:
    # merge input las files into single dataframe
        output = pd.merge(main_df, merge_df, on = depth_column, how = "outer")
    else:
        output = main_df

    # remove unwanted logs
    if remove == True:
        for rem in remove_logs:
            print (list(output.columns.values))
            print (output.head())
            if rem in list(output.columns.values):
                output.drop(rem, axis = 1, inplace = True)

    # rename logs
    if rename != None:
        output.rename(columns = rename, inplace=True)

    output.dropna(inplace = True, how = "all")

    # interpolate logs to common depth range/sample interval
    logs_list = list(output.columns.values)
    logs_list.remove(depth_column)
    output_interp = pd.DataFrame()
    output_interp[depth_column] = np.arange(start_depth, end_depth + output_int, output_int)
    print (logs_list)
    print (output)

    for log in logs_list:
        filt = pd.DataFrame(data = output[[depth_column, log]], copy = True)
        filt.dropna(how = "any", inplace = True)
        log_start = filt[depth_column].min()
        log_end = filt[depth_column].max()
        log_depth = output_interp[depth_column].where((output_interp[depth_column] >= log_start) & (output_interp[depth_column] <= log_end))
        #print(filt)
        f = interp1d(filt[depth_column], filt[log], kind = "linear", fill_value = "extrapolate", assume_sorted=True)
        log_data = f(log_depth)
        log_out = pd.DataFrame()
        log_out[depth_column] = log_depth
        log_out[log] = log_data
        output_interp = output_interp.merge(log_out, on = depth_column, how = "left")

    # create alternative depth logs for vertical wells
    if create_vert_TVDml == True:
        if vertical == True:
            output_interp["TVDss"] = output_interp[depth_column] - kb
            output_interp["TVDml"] = output_interp["TVDss"] - wd
        else:
            print("function does not yet exist for non vertical wells!")

    # write and output las file
    output_file = out_root + "\\" + well + "\\" + well + out_suffix + ".las"
    output_filenames.append(output_file)

    output.rename(columns={depth_column: "DEPTH"}, inplace=True)
    output.set_index("DEPTH", inplace = True)
    output.dropna(how = "all", axis = 0, inplace = True)#, subset = ss)
    output.sort_index(axis = 1, ascending = True, inplace = True)


    with open(output_file, mode = 'w') as lasfile:
        las = lasio.LASFile()
        las.depth = ["DEPTH"]
        las.well["LOC"].value = str(easting) + "," + str(northing)
        las.well["WELL"].value = str(well)
        las.well["NULL"].value = null
        las.add_curve("DEPTH", output.index.values, unit = logs_dict.get(depth_column))
        for log in list(output.columns.values):
            las.add_curve(log, output[log], unit = logs_dict.get(log))
        las.write(lasfile, version=2, fmt = "%10.10g")

    return output_filenames

merge_las(root, well, kb, wd, create_vert_TVDml, file_name, depth_column, merge_root, merge_name, merge_depth_col, out_root, null, start_depth, end_depth, rename)







