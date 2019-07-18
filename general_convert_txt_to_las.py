import lasio
import pandas as pd

## takes tab delimited text files with filename in following format: root\well_scenario.txt
#outputs a las file with same name to same directory
#column names and units must be specified in the order they appear in the text file
#header lines must be specified
#requires pandas and lasio libraries


root = r"\\DST-CHTECH2\Data\Country Folders\Morocco\RABAT DEEP\Technical\Drilling\Wireline_Final_Data\Reprocessed Sonic"
scenarios = [""]
wells = ["RSD-1_vels"]
header_lines = 1

column_names = ["MD", "BS", "DTCO_FINAL", "DTSH_FINAL", "GR_EDTC", "PR_FINAL", "TENS", "VPVS_FINAL", "Vp", "Vs"] #"["DEPTH","VP", "VS", "RHOB", "PHIE", "VSH", "SW"]
column_units = ["m", "in", "us/ft", "us/ft", "API", "ratio", "lbf", "ratio", "m/s", "ms"] #["m", "m/s", "m/s", "g/cc", "frac", "frac", "frac"]

#root = r'T:\Country Folders\Morocco\USER WORKING\JW\Morocco\Data\LBS-1'
#scenarios = [""]
#wells = ["LBS-1_INI"]
#header_lines = 2

#column_names = ["DEPTH" ,"CALI" ,"PEF", "DRHO", "RHOB", "GR", "NPHI", "TNPH", "VCL", "SW", "PHIE", "PHIT", "ILM", "ILD", "LLD", "SFLU", "DT", "SP", "VP", "DT_edit", "VP_ed", "P_HYDRO", "RHOB (Gardner_shale)", "RHOB (Gardner_sand)", "RHOB (Gardner_std)", "RHOB (Gardner_lsq)", "RHOB gardner splice", "RHOB gardner-data splice", "Vs_Castagna_shale", "Vs_Castagna_sand", "Vs_Castagna_splice", "DTSM_Castgna_splice"]
#column_units = ["m", "in", "-", "g/cc", "g/cc", "API", "frac", "v/v", "frac", "frac", "frac", "frac", "ohm.m", "ohm.m", "ohm.m", "ohm.m", "us/ft", "mv", "m/s", "us/ft", "m/s", "MPa", "g/cc", "g/cc", "g/cc", "g/cc", "g/cc", "g/cc", "m/s", "m/s", "m/s", "us/ft"]



#root = r"T:\Country Folders\Morocco\USER WORKING\JW\Morocco\Data\RSD-1\Logs\180423"
#scenarios = [""]
#wells = ["RSD-1"]
#header_lines = 1

#column_names = ["DEPTH", "TVD","GRCFM", "DTC", "SVC", "TTC", "DTCMH", "DTCML", "Vp_DTC", "RhoB_gardner"]
#column_units = ["m", "m", "API", "us/ft", "", "us", "us/ft", "us/ft", "m/s", "g/cc"]

#column_names = ["DEPTH", "GR_Run6b_230418", "DTC_Run6b_230418", "Vp_Run6b_230418", "RHOB_Gardner_Run6b_230418"]
#column_units = ["m", "API", "us/ft", "m/s", "g/cc"]


def data_load(filepath,header_lines=None, column_names=None, null_value=None):
    input = pd.read_table(filepath, index_col=False, sep=('\t'), delim_whitespace=True, header=header_lines, names=column_names, lineterminator='\n', na_values=null_value)
    return input

filepaths = []
for well in wells:
    for scenario in scenarios:
        path = root + "\\" + well + scenario + ".txt"
        output = root + "\\" + well + scenario + ".las"
        filepaths.append(path)
        #print (len(column_names))
        #print (len(column_units))
        logs = data_load(path, header_lines = header_lines, column_names = column_names)

        with open(output, mode = 'w') as lasfile:
            las = lasio.LASFile()
            las.depth = ["DEPTH"]
            las.well["WELL"].value = str(well)
            las.well["NULL"].value = -999.25
            for log, name, unit in zip(logs, column_names, column_units):
                las.add_curve(name, logs[name].tolist(), unit = unit)
            las.write(lasfile, version=2, fmt='%10.9g')


