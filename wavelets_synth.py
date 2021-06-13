import numpy as np
import scipy.signal as signal
import pandas as pd

"""
## not sure this is worth it

def import_wavelet(wavelet_type = "statistical", wvlt_export = "dtect", dt = 0.5, wvlt_path = None, wvlt_file = None, wvlt_cfreq = 30, wvlt_phase = 0, wvlt_length = 200):
    
    if wvlt_type == 'ricker':
        wvlt_t, wvlt_amp = ricker(wvlt_cfreq, wvlt_phase, dt, wvlt_length)
        a1, t1 = wvlt_amp, wvlt_t

    elif wvlt_type == 'bandpass':
        wvlt_t, wvlt_amp = wvlt_bpass(f1, f2, f3, f4, wvlt_phase, dt, wvlt_length)
        a1, t1 = wvlt_amp, wvlt_t

    elif wvlt_type == 'statistical':
        if wvlt_export == "dtect":
            wvlt_t, wvlt_amp, a1, t1 = wvlt_dtect(wvlt_path, wvlt_file, dt)
        elif wvlt_export == "geoview":
            wvlt_t, wvlt_amp, a1, t1 = wvlt_geoview(wvlt_path, wvlt_file, dt)
        elif wvlt_export == "kingdom":
            wvlt_t, wvlt_amp, a1, t1 = wvlt_kingdom(wvlt_path, wvlt_file, dt)
        else:
            print("build me!")
    
"""




# create wavelet
def ricker(cfreq, phase, dt, wvlt_length):
    '''
    Calculate a zero-phase ricker wavelet

    Usage:
    ------
    t, wvlt = wvlt_ricker(cfreq, dt, wvlt_length)

    cfreq: central frequency of wavelet in Hz
    phase: wavelet phase in degrees
    dt: sample rate in ms
    wvlt_length: length of wavelet in ms
    '''

    import numpy as np
    import scipy.signal as signal

    nsamp = int(wvlt_length / dt + 1)
    t_max = wvlt_length * 0.5
    t_min = -t_max

    t = np.arange(t_min, t_max, dt)
    t = np.linspace(-int(wvlt_length / 2), int((wvlt_length) / 2), int(wvlt_length / dt)+1)
    t = t/1000
    wvlt = (1.0 - 2.0 * (np.pi ** 2) * (cfreq ** 2) * (t ** 2)) * np.exp(-(np.pi ** 2) * (cfreq ** 2) * (t ** 2))

    if phase != 0:
        phase = phase * np.pi / 180.0
        wvlth = signal.hilbert(wvlt)
        wvlth = np.imag(wvlth)
        wvlt = np.cos(phase) * wvlt - np.sin(phase) * wvlth
    t = t * 1000
    return t, wvlt


def wvlt_bpass(f1, f2, f3, f4, phase, dt, wvlt_length):
    '''
    Calculate a trapezoidal bandpass wavelet

    Usage:
    ------
    t, wvlt = wvlt_ricker(f1, f2, f3, f4, phase, dt, wvlt_length)

    f1: Low truncation frequency of wavelet in Hz
    f2: Low cut frequency of wavelet in Hz
    f3: High cut frequency of wavelet in Hz
    f4: High truncation frequency of wavelet in Hz
    phase: wavelet phase in degrees
    dt: sample rate in seconds
    wvlt_length: length of wavelet in seconds
    '''

    from numpy.fft import fft, ifft, fftfreq, fftshift, ifftshift
    wvlt_length = wvlt_length / 1000
    dt = dt / 1000

    nsamp = int(wvlt_length / dt + 1)

    freq = fftfreq(nsamp, dt)
    freq = fftshift(freq)
    aspec = freq * 0.0
    pspec = freq * 0.0

    # Calculate slope and y-int for low frequency ramp
    M1 = 1 / (f2 - f1)
    b1 = -M1 * f1

    # Calculate slop and y-int for high frequency ramp
    M2 = -1 / (f4 - f3)
    b2 = -M2 * f4

    # Build initial frequency and filter arrays
    freq = fftfreq(nsamp, dt)
    freq = fftshift(freq)
    filt = np.zeros(nsamp)

    # Build LF ramp
    idx = np.nonzero((np.abs(freq) >= f1) & (np.abs(freq) < f2))
    filt[idx] = M1 * np.abs(freq)[idx] + b1

    # Build central filter flat
    idx = np.nonzero((np.abs(freq) >= f2) & (np.abs(freq) <= f3))
    filt[idx] = 1.0

    # Build HF ramp
    idx = np.nonzero((np.abs(freq) > f3) & (np.abs(freq) <= f4))
    filt[idx] = M2 * np.abs(freq)[idx] + b2

    # Unshift the frequencies and convert filter to fourier coefficients
    filt2 = ifftshift(filt)
    Af = filt2 * np.exp(np.zeros(filt2.shape) * 1j)

    # Convert filter to time-domain wavelet
    wvlt = fftshift(ifft(Af))
    wvlt = np.real(wvlt)
    wvlt = wvlt / np.max(np.abs(wvlt))  # normalize wavelet by peak amplitude

    # Generate array of wavelet times
    t = np.linspace(-wvlt_length * 0.5, wvlt_length * 0.5, nsamp)

    # Apply phase rotation if desired
    if phase != 0:
        phase = phase * np.pi / 180.0
        wvlth = signal.hilbert(wvlt)
        wvlth = np.imag(wvlth)
        wvlt = np.cos(phase) * wvlt - np.sin(phase) * wvlth

    t = t * 1000
    return t, wvlt


def wvlt_dtect(wvlt_path, wvlt_file, dt):
    wvlt_df = pd.read_csv(filepath_or_buffer=wvlt_path + "\\" + wvlt_file, delim_whitespace=True, header=None,
                            dtype=np.float64)
    wvlt_df.columns = ["Time", "Amp"]
    upper_fill, lower_fill = 0, 0
    # f = interpolate.interp1d(wvlt_df["Time"],wvlt_df["Amp"], kind = "cubic", bounds_error  = False, fill_value = (upper_fill, lower_fill))
    samp_freq = int(wvlt_df["Time"].iloc[1] - wvlt_df["Time"].iloc[0])
    print ("input sample freq = %f" % samp_freq)
    wvlt_length = (len(wvlt_df.index)) * samp_freq

    nsamp = int(wvlt_length / dt)
    print("output sample freq = %f" % dt)

    wvlt, t = signal.resample(wvlt_df["Amp"].tolist(), num=nsamp, t=wvlt_df["Time"].tolist())

    return t, wvlt, wvlt_df["Amp"].tolist(), wvlt_df["Time"].tolist()


# wvlt_path = r"T:\Country Folders\Morocco\USER WORKING\JW\Morocco\2018 QI Project\RSD_1 Revisit\wavelets"
# wvlt_file = "5-35deg_5-30deg_RSD-1_200x200lines_stat.txt"

def wvlt_geoview(wvlt_path, wvlt_file, dt):
    wvlet_df = None
    params_df = pd.read_table(wvlt_path + "\\" + wvlt_file, skiprows=26, nrows=4, delim_whitespace=True)
    params_df.index = ["Sample Rate", "Sample no. T = zero", "No. Samples", "Phase Rotation"]
    wvlt_df = pd.read_table(wvlt_path + "\\" + wvlt_file, skiprows=30, delim_whitespace=True, usecols=[0])
    wvlt_df.columns = ["Amp"]
    samp_freq = params_df.iloc[0]
    wvlt_length = float(params_df.iloc[0] * params_df.iloc[2])

    wvlt_df["Time"] = np.linspace((wvlt_length / 2) * -1, (wvlt_length / 2 - samp_freq), params_df.iloc[2])
    nsamp = int(wvlt_length / dt)

    wvlt, t = signal.resample(wvlt_df["Amp"].tolist(), num=nsamp, t=wvlt_df["Time"].tolist())

    return t, wvlt, wvlt_df["Amp"].tolist(), wvlt_df["Time"].tolist()


def wvlt_kingdom(wvlt_path, wvlt_file, dt):
    wvlt_df = pd.read_csv(filepath_or_buffer=wvlt_path + "\\" + wvlt_file, skiprows=4, delim_whitespace=True)
    wvlt_df.columns = ["Time", "Amp"]
    wvlt_df["Time"] = wvlt_df["Time"] * 1000

    samp_freq = int(abs((wvlt_df["Time"].iloc[1] - wvlt_df["Time"].iloc[0])))
    wvlt_length = (len(wvlt_df.index)) * samp_freq
    nsamp = int(wvlt_length / dt)

    wvlt, t = signal.resample(wvlt_df["Amp"].tolist(), num=nsamp, t=wvlt_df["Time"].tolist())
    return t, wvlt, wvlt_df["Amp"].tolist(), wvlt_df["Time"].tolist()
