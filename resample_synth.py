import pandas as pd
from scipy.interpolate import interp1d
from scipy.signal import resample

import numpy as np
import math
from .las_handling import write_las


def resample_depth_in_twt(df, log_list, depth="MDKB", twt="TWT", output_sample_int=0.5,
                  tmin=None, tmax=None, las_out = False, out_path = None, out_file = None):

    if tmin == None:
        tmin = 0
    if tmax == None:
        tmax = math.ceil(float(df[twt].max()))

    print ("tmin = %i, tmax = %i" % (tmin, tmax))

    nsamp_out = int((tmax - tmin) / output_sample_int)  + 1
    print("nsamp = %i" % nsamp_out)

    df_resamp = pd.DataFrame()
    df_resamp[twt] = np.linspace(tmin, tmax, nsamp_out)
    units_dict = {}
    units_dict[twt] = "ms"
    units_dict[depth] = "m"
    # print (log_list)

    f = interp1d(df[twt], df[depth], kind="linear", fill_value="extrapolate", assume_sorted=True)
    df_resamp[depth] = f(df_resamp[twt])


    for log in log_list:
        units_dict[log] = None
        f = interp1d(df[twt], df[log], kind="linear", fill_value="extrapolate", assume_sorted=True)
        df_resamp[log] = f(df_resamp[twt])

    log_list.append(depth)

    if las_out == True:
        write_las(out_path, out_file, df_resamp, log_list, twt, units_dict)

    return df_resamp



def resample_synth(df, syn_list, log_list, depth="MDKB", twt="TWT", output_sample_int=2,
                  output_start=0, output_end=3000):
    import pandas as pd
    from scipy.interpolate import interp1d
    from scipy.signal import resample
    # check the mod stuff - may not work
    mod = (output_end - output_start) % output_sample_int


    nsamp_out = int((output_end - output_start - mod) / output_sample_int)
    df = df.loc[df[twt] <=output_end-mod]
    print("nsamp = %f" % nsamp_out)

    df_resamp = pd.DataFrame()
    t = df[twt].tolist()

    for trace in syn_list:
        x = df[trace].tolist()

        df_resamp[trace], df_resamp[twt] = resample(x[:-1], nsamp_out, t=t[:-1])  # resample trace

    f = interp1d(df[twt], df[depth], kind="linear", fill_value="extrapolate", assume_sorted=True)

    df_resamp[depth] = f(df_resamp[twt])


    for log in log_list:
        #fill_value = (float(df[log][0]), float(df[log][-1]))
        fill_value = (0, 0)
        f = interp1d(df[twt], df[log], kind="linear", fill_value=fill_value, assume_sorted=True)
        df_resamp[log] = f(df_resamp[twt])


    return df_resamp

"""
def resample_data(df, syn_list, log_list, depth="MDKB", twt="TWT", output_sample_int=2,
                  output_start=0, output_end=3000):
    nsamp_out = int((output_end - output_start) / output_sample_int)
    print("nsamp = %f" % nsamp_out)

    df_resamp = pd.DataFrame()

    t = df[twt].tolist()
    for trace in syn_list:
        units_dict[trace] = None

        x = df[trace].tolist()

        df_resamp[trace], df_resamp[twt] = resample(x[:-1], nsamp_out, t=t[:-1])  # resample trace

    f = interp1d(df[twt], df[depth], kind="cubic", fill_value="extrapolate", assume_sorted=True)
    df_resamp[depth] = f(df_resamp[twt])


    for log in log_list:
        units_dict[log] = None
        f = interp1d(df[twt], df[log], kind="linear", fill_value="extrapolate", assume_sorted=True)
        df_resamp[log] = f(df_resamp[twt])



    return df_resamp

"""

def resample_volume(sgy, nmodel, nangles, dt_in, dt_out, twt_in = None):

    from scipy.signal import resample
    import numpy as np
    nsamp_in = int(len(sgy[0, 0, 0, :]))

    mod = nsamp_in % (dt_out / dt_in)
    print (mod)
    print("No. samples in: %i" % nsamp_in)

    nsamp_out = int((nsamp_in-mod) / (dt_out / dt_in))

    print("No. samples out: %i" % nsamp_out)
    print("Sample reduction factor (dt_out/dt_in): %f" % (dt_out / dt_in))
    sgy_out = np.ndarray(shape=(1, len(sgy[0, :, 0, 0]), len(sgy[0, 0, :, 0]), nsamp_out))
    #print (np.shape(sgy))
    if mod != 0:
        sgy = sgy[:, :, :, :-int(mod)]
        twt_in = twt_in[:-int(mod)]
    #print(np.shape(sgy))
    for i in range(nmodel):
        if nangles == 1:
            syn_resamp, twt_out = resample(sgy[0, i, 0, :], num=nsamp_out, t=twt_in)
            sgy_out[0, i, 0, :] = syn_resamp

        else:
            for j in range(nangles):
                #print (j)
                if twt_in != None:
                    #max_amp = max(sgy[0, i, j, :])
                    #min_amp = min(sgy[0, i, j, :])
                    syn_resamp, twt_out = resample(sgy[0, i, j, :], num=nsamp_out, t = twt_in)
                    #print (np.shape(syn_resamp))
                    #syn_resamp[~np.isfinite(syn_resamp)] = 0
                    #syn_resamp[syn_resamp == np.inf] = 0
                    #syn_resamp[syn_resamp == np.NINF] = 0

                    sgy_out[0, i, j, :] = syn_resamp


                else:

                    syn_resamp = resample(sgy[0, i, j, :], num=nsamp_out)
                    #syn_resamp[~np.isfinite(syn_resamp)] = 0
                    sgy_out[0, i, j, :] = syn_resamp
                    twt_out = None
    print (np.inf in sgy_out.flatten())
    if twt_out is None:
        return sgy_out
    else:
        return sgy_out, twt_out


"""
option to add las output?

def resample_synth(df, syn_list, log_list, depth="MDKB", twt="TWT", output_sample_int=2,
                  output_start=0, output_end=3000, las_out = True, out_path = None, out_file = None):

    nsamp_out = int((output_end - output_start) / output_sample_int)
    df = df.loc[df[twt] <=output_end]
    print("nsamp = %f" % nsamp_out)

    df_resamp = pd.DataFrame()

    units_dict = {}
    units_dict[twt] = "ms"
    units_dict[depth] = "m"
    # print (log_list)

    t = df[twt].tolist()
    # print (df[twt].max())
    for trace in syn_list:
        units_dict[trace] = None

        x = df[trace].tolist()

        df_resamp[trace], df_resamp[twt] = resample(x[:-1], nsamp_out, t=t[:-1])  # resample trace
        #print (len(df_resamp[twt]))
        #print (df_resamp[twt])

        # print (df_resamp[twt])

    f = interp1d(df[twt], df[depth], kind="linear", fill_value="extrapolate", assume_sorted=True)
    df_resamp[depth] = f(df_resamp[twt])


    for log in log_list:
        units_dict[log] = None
        f = interp1d(df[twt], df[log], kind="linear", fill_value="extrapolate", assume_sorted=True)
        df_resamp[log] = f(df_resamp[twt])

    if las_out == True:
        write_las(out_path, out_file, df_resamp, syn_list, twt, units_dict)

    return df_resamp





def resample_data(df, syn_list, log_list, out_path, out_file, depth="MDKB", twt="TWT", output_sample_int=2,
                  output_start=0, output_end=3000, las_out = False):
    nsamp_out = int((output_end - output_start) / output_sample_int)
    print("nsamp = %f" % nsamp_out)

    df_resamp = pd.DataFrame()

    units_dict = {}
    units_dict[twt] = "ms"
    units_dict[depth] = "m"
    log_list.remove(twt)
    log_list.remove(depth)

    t = df[twt].tolist()
    for trace in syn_list:
        units_dict[trace] = None

        x = df[trace].tolist()

        df_resamp[trace], df_resamp[twt] = resample(x[:-1], nsamp_out, t=t[:-1])  # resample trace

    f = interp1d(df[twt], df[depth], kind="cubic", fill_value="extrapolate", assume_sorted=True)
    df_resamp[depth] = f(df_resamp[twt])


    for log in log_list:
        units_dict[log] = None
        f = interp1d(df[twt], df[log], kind="linear", fill_value="extrapolate", assume_sorted=True)
        df_resamp[log] = f(df_resamp[twt])

    # write_las(out_path, out_file, data_time_out, syn_list, twt, units_dict)

    return df_resamp
"""