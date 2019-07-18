import pandas as pd
import numpy as np
import lasio
from scipy.interpolate import interp1d
root = r'F:\Data\Namibia\las'
wells = ["2012_13-1"] # ["Murombe-1"] #"["Shark-1"] #["Wingat-1"] #
filename = "2012_13-1_run2" #""Murombe-1_merge" #"Shark-1_merge" #"Wingat-1_merge"
depth_column = "MD" #""DEPTH"

out_root = root
out_filename = "_interpolate"

null = -999.25
output_int = 0.1524
kb = [12.5] #[28.05] #[27] #[27.1] #
wd = [688]#[1387.5] #[274.7] #[1004.5]#

start_depth = "null"
end_depth = "null"

create_vert_TVDml = True


def interpolate_log(root, wells, filename, depth_column, start_depth, end_depth, null, output_int, out_root, out_filename, create_vert_TVDml, kb, wd):

    for well, kb_well, wd_well in zip(wells, kb, wd):

        filepath = root + "\\" + well + "\\" + filename + ".las"
        data = lasio.read(filepath)#, null_subs = True)

        logs_list = []
        units_list = []
        logs_dict = {}

        for curve in data.curves:
            logs_list.append(curve.mnemonic)
            units_list.append(curve.unit)
            logs_dict[curve.mnemonic] = curve.unit


        las_df = pd.DataFrame()

        for log in logs_list:
            las_df[str(log)] = data[str(log)]
            las_df[str(log)] = np.where(las_df[str(log)] == null, np.nan, las_df[str(log)])

        logs_list.remove(depth_column)

        if depth_column != "DEPTH":
            las_df.rename(columns={depth_column: "DEPTH"}, inplace=True)
            depth_column = "DEPTH"

        if start_depth == "null":
            start_depth = las_df[depth_column].iloc[0]
        if end_depth == "null":
            end_depth = las_df[depth_column].iloc[-1]


        #print (logs_list)
        #print (start_depth)
        #print (end_depth)
        #print (las_df.head())
        output_df = pd.DataFrame()
        output_df[depth_column] = np.arange(start_depth, end_depth + output_int, output_int)
        #print (output_df.head())

        for log in logs_list:
            filt = pd.DataFrame(data = las_df[[depth_column, log]], copy = True)
            filt.dropna(how = "any", inplace = True)
            log_start = filt[depth_column].min()
            log_end = filt[depth_column].max()
            log_depth = output_df[depth_column].where((output_df[depth_column] >= log_start) & (output_df[depth_column] <= log_end))
            #print(filt)
            f = interp1d(filt[depth_column], filt[log], kind = "linear", fill_value = "extrapolate", assume_sorted=True)
            log_data = f(log_depth)
            log_out = pd.DataFrame()
            log_out[depth_column] = log_depth
            log_out[log] = log_data
            output_df = output_df.merge(log_out, on = depth_column, how = "left")
            #print (output_df.head())


        if create_vert_TVDml == True:
            output_df["TVDss"] = output_df[depth_column] - kb_well
            output_df["TVDml"] = output_df["TVDss"] - wd_well

        #print (output_df.head)

        output_file = out_root + "\\" + well + "\\" + well + out_filename + ".las"
            
        with open(output_file, mode = 'w') as lasfile:
            las = lasio.LASFile()
            las.depth = ["DEPTH"]
            las.well["WELL"].value = str(well)
            las.well["NULL"].value = null
            for log in list(output_df.columns.values):
                las.add_curve(log, output_df[log], unit = logs_dict.get(log))
            las.write(lasfile, version=2, fmt = "%10.10g")
            


    return


interpolate_log(root, wells, filename, depth_column, start_depth, end_depth, null, output_int, out_root, out_filename, create_vert_TVDml, kb, wd)


