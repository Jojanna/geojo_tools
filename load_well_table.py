import pandas as pd

def read_well_table(filename, well_name):
    well_table = pd.read_excel(filename)
    well_data = well_table.where(well_table["Name"] == well_name).dropna(axis = 0, how = "all")
    well_data.rename(columns = {"WD (m MSL)": "wd", "KB (m MSL)": "kb","TD (m KB)": "td", "E": "coord_e", "N": "coord_n"}, inplace = True)
    well_dict = well_data.to_dict("records")[0]
    return well_dict



#filename = r'F:\Namibia\well_database.xlsx'
#read_well_table(filename, well_name = "Wingat-1")
