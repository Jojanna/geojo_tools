import pandas as pd
import sys

root = r'F:\Namibia\Exports\Polygons'


def polygon_convert(root, file_in, file_out):
    path = root + "\\" + file_in
    data = pd.read_table(path, skiprows = 8, skipfooter = 3, delim_whitespace= True, skip_blank_lines=True, header = None)
    data = data.transpose()

    polygon = pd.DataFrame()
    x = data.iloc[[0, 2, 4, 6, 8]] #pd.DataFrame()
    y = data.iloc[[1, 3, 5, 7, 9]]

    c = []
    d = []

    for a in x.columns.values:
        c = c +  x[a].tolist()
    for b in y.columns.values:
        d = d + y[b].tolist()

    polygon["x"] = c
    polygon["y"] = d
    polygon.dropna(how = "any", axis = 0, inplace = True)
    print(polygon)

    polygon.to_csv(root + "\\" + file_out, sep = "\t", index = False)

    return

#file_in = "RJWF_W94_deep_4160_contact.dat"
#file_out = "RJWF_W94_deep_4160_contact_edit_ordering.dat"
#polygon_convert(root, file_in, file_out)

files = ["RJWF_S94_Top_Reservoir_3620m.dat", "RJWF_S94_Top_Reservoir_3660m.dat", "RJWF_W_W94_Deep_4140m_contact.dat", "RJWF_W_W94_Deep_W06_4295m_contact.dat", "RJWF_W94_4045_contact.dat"]
files_out = ["RJWF_S94_Top_Reservoir_3620m_edit_ordering.dat", "RJWF_S94_Top_Reservoir_3660m_edit_ordering.dat", "RJWF_W_W94_Deep_4140m_contact_edit_ordering.dat", "RJWF_W_W94_Deep_W06_4295m_contact_edit_ordering.dat", "RJWF_W94_4045_contact_edit_ordering.dat"]
for x, y in zip(files, files_out):
    polygon_convert(root, x, y)

