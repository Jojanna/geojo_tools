import pandas as pd

def full_stack(df, theta1, synth_prefix = "Syn_", z = "MDKB", min_angle=0, max_angle=45, write_csv = True, outpath= None, file = None):
    if file == None:
        file =  "Full_stack_Synthetic_%s_%.0f-%.0fdeg" % (z, min_angle, max_angle)
    list_syn = []
    # [x for x in log_list if x not in syn_list]
    angles = [x for x in theta1 if (x <= max_angle and x >= min_angle)]

    for angle in angles:
        list_syn.append(synth_prefix + str(angle))

    fstk = pd.DataFrame()
    fstk[z] = df[z]
    fstk["%.0f-%.0f deg FULL STACK" % (min_angle, max_angle)] = df[list_syn].mean(axis=1, skipna=True)
    if write_csv == True:
        fstk.to_csv(outpath + "\\" + file + ".txt", sep="\t", mode="w", header=True, index=False)

    return fstk

