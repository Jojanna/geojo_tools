import numpy as np
import pandas as pd

def create_blocky_logs(df, vp_mod, vs_mod, rho_mod, depth_bounds):
    df["VP_Blocky"] = np.nan
    df["VS_Blocky"] = np.nan
    df["RHOB_Blocky"] = np.nan

    for i in range(len(depth_bounds)-1):
        df.loc[(df["MDKB"] > depth_bounds[i]) & (df["MDKB"] <= depth_bounds[i+1]), "VP"] = vp_mod[i]
        df.loc[(df["MDKB"] > depth_bounds[i]) & (df["MDKB"] <= depth_bounds[i + 1]), "VS"] = vs_mod[i]
        df.loc[(df["MDKB"] > depth_bounds[i]) & (df["MDKB"] <= depth_bounds[i + 1]), "RHOB"] = rho_mod[i]

    return df


