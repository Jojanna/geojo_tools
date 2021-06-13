import lasio
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

## takes las files with filename in following format: root\well_scenario.las and displays as a chart for quick viewing
# if chart = True, the output graph will be saved with filename: root\well_scenario.png
# use the null value to prevent these values being plotted

root = r"\\dst-chtech2\data\Country Folders\Morocco\LIXUS\Subsurface\Petrophysics\RF IP Work\Petrophysics .Las File per Well"
#scenarios = [""]
wells = ["Anchois-1_ST_Final_RF", "Deep Thon-1_Final_RF", "Merou-1_Final_RF", "Deep Thon-1_Final_RF_DTS", "Merou-1_Final_RF_DTS", "LAR-A-1_Final_RF"]
null = -999.0000
depth_column = "DEPTH"

depth_range = None #[[3520, 4460.0032]] #"null" #[[6775, 6975]]

chart = True

filepaths = []

for well in wells:
    path = root + "\\" + well + ".las"
    filepaths.append(path)

    data = lasio.read(path)
    logs = []
    units = []
    for curve in data.curves:
        logs.append(curve.mnemonic)
        units.append(curve.unit)

    logs_df = pd.DataFrame()

    for log in logs:
        logs_df[str(log)] = data[str(log)]
        logs_df[str(log)] = np.where(logs_df[str(log)] == null, np.nan, logs_df[str(log)])
        filt = (logs_df[str(log)] <=0)
        logs_df[str(log)].where(~filt, np.nan, inplace=True)
    logs_df.dropna(how = "all", inplace = True)
    print(logs)
    logs.remove(str(depth_column))
    md_unit = units[0]
    del units[0]

    md = logs_df[str(depth_column)].tolist()
    ax = np.arange(1, (len(logs)), 1)
    ax_ref = list(map(lambda x: "ax" + "_" + str(x), ax))

    if depth_range is None:
        min = np.min(md)
        max = np.max(md)
        range = [min, max]
        print (range)

    #idx = wells.index(well)
    #range = depth_range[idx]

    fig1, ((ax_ref)) = plt.subplots(1, len(logs))
    fig1.set_size_inches(len(logs) * 3, 12.3)
    size = 8
    plt.title(str(well), ha="center")
    for log, ax_n, unit in zip(logs, ax_ref, units):
        #print (log)
        #print (logs_df[log])
        df = pd.Series(index = md, data = logs_df[log].tolist(), ).dropna(axis = 0)
        #print (df)
        if unit in ["ohm.m", "ohmm", "Ohm.m", "Ohmm", "OHM.M", "OHMM"]:
            ax_n.semilogx(df.tolist(), df.index.tolist(), linewidth=0.5)
        else:
            ax_n.plot(df.tolist(), df.index.tolist(), linewidth=0.5)
        ax_n.set_xlabel(str(log) + " " + str(unit), fontsize = size)
        ax_n.set_ylabel("Depth" + " " + md_unit, fontsize = size)
        ax_n.set_axisbelow(True)
        ax_n.grid()
        ax_n.tick_params(axis='both', which='major', labelsize=6)
        ax_n.set_ylim(range[0], range[1])


        ax_n.invert_yaxis()


    plt.tight_layout()

    if chart == True:
        chart_name = root + "\\" + well + ".png"
        fig1.savefig(chart_name, dpi = 400)
    plt.show()

