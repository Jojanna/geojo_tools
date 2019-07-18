from matplotlib.colors import ListedColormap, LinearSegmentedColormap, Normalize
import numpy as np
import re


def cmap_load(path, file, sep="\s+", cols=None, header=None):
    file = open(path + "\\" + file)
    data = file.read()
    data = data.strip().split("\n")
    length = len(data) - header
    all_lines = []
    for line in data:
        line = line.strip()
        line = re.split(sep, line)
        if cols == None:
            cols = len(line)
        all_lines.append(line)
        for x in line:
            float_line = []
            if type(x) != str:
                x = float(x)
                float_line.append(x)
            else:
                float_line = line

    if header != None:
        if header == 1:
            del all_lines[0]
        else:
            head = np.arange(0, header - 1, 1)
            for h in head:
                del all_lines[h]

    out = np.array(all_lines, dtype=float)
    out = np.transpose(out)
    out = out.reshape(cols, length)

    out = out.transpose()

    cmap = ListedColormap(out)

    return cmap


def load_seismic_petrel():
    cb = cmap_load(r"C:\Users\Joannaw\Documents\C Drive Chariot\Resources\ColourBars\Geoview", "Seismic_Petrel_JW_256_strip.txt", header = 1)
    return cb
