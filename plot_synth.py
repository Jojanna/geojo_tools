import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.patches as mpatches
import matplotlib.lines as mlines

from geojo.colourmaps import load_seismic_petrel


def upper_lower_ref(df_resamp, upper_ref_z, lower_ref_z, angles, synth_prefix = "Syn_", mdkb = "MDKB", twt = "TWT"):
    # upper_df = ((data_time["Phi_Scenario"] == phi_scenarios[0]) & (data_time["Fluid_Scenario"] == fluid_scenarios[0]) & (data_time["Ovb"] == overburdens[0]) & (data_time["MDKB"] >= upper_ref_z))
    # lower_df = ((data_time["Phi_Scenario"] == phi_scenarios[0]) & (data_time["Fluid_Scenario"] == fluid_scenarios[0]) & (data_time["Ovb"] == overburdens[0]) & (data_time["MDKB"] >= lower_ref_z))
    upper_df = df_resamp.loc[(df_resamp[mdkb] >= upper_ref_z), ["TWT"]]
    lower_df = df_resamp.loc[(df_resamp[mdkb] >= lower_ref_z), ["TWT"]]

    upper_ref_t = float(upper_df.iloc[0].values)
    lower_ref_t = float(lower_df.iloc[0].values)
    # print (upper_ref_t, lower_ref_t)

    # print ("upper_ref_z = %f m MDKB, upper_ref_t = %f ms TWT, lower_ref_z = %f m MDKB, lower_ref_t = %f ms TWT" % (upper_ref_z, upper_ref_t, lower_ref_z, lower_ref_t))

    #   Copy convolved top/base reflectivity values to Lists for easier plotting
    line1 = []
    line2 = []
    for angle, i in zip(angles, range(0, len(angles))):
        # print (angle)
        # print (data_time["Syn_" + str(angle)].loc[(data_time["TWT"] == upper_ref_t)])
        line1.append(
            [i, float(df_resamp[synth_prefix + str(angle)].loc[(df_resamp[twt] == upper_ref_t)].values)])  # scale?
        line2.append([i, float(df_resamp[synth_prefix + str(angle)].loc[(df_resamp[twt] == lower_ref_t)].values)])

    line1 = list(map(list, zip(*line1)))
    line2 = list(map(list, zip(*line2)))

    return line1, line2, upper_ref_t, lower_ref_t

def plot_vawig_wiggle(axhdl, df_resamp, excursion, angles, synth_prefix = "Syn_", twt = "TWT", vminmax=[-5, 5], fill=True, colours="seismic_petrel"):

    if colours == "seismic_petrel":
        seismic_petrel = load_seismic_petrel()
        colours = seismic_petrel

    col_list = list(map(lambda x: synth_prefix + str(x), angles))
    print(col_list)
    max_amp = df_resamp[col_list].abs().values.max()

    print("max amplitude = %f" % max_amp)

    # [ntrc, nsamp] = data.shape

    # t = np.hstack([0, t, t.max()])
    cmap = plt.get_cmap(colours)

    for angle, i in zip(angles, range(0, len(angles))):

        df_resamp[synth_prefix + str(angle)] = excursion * df_resamp[synth_prefix + str(angle)] / max_amp + i

        norm = plt.Normalize(vmin=vminmax[0], vmax=vminmax[1])

        # tbuf = np.hstack([i, tbuf, i])
        # print (angle, i)
        axhdl.plot(df_resamp["Syn_scale_" + str(angle)], df_resamp["TWT"], color='black', linewidth=0.5)

        # plt.fill_betweenx(data_time["TWT"], data_time["Syn_scale_" + str(angle)], i, where=data_time["Syn_scale_" + str(angle)] > i, facecolor=[0.6, 0.6, 1.0], linewidth=0)
        # plt.fill_betweenx(data_time["TWT"], data_time["Syn_scale_" + str(angle)], i, where=data_time["Syn_scale_" + str(angle)] < i, facecolor=[1.0, 0.7, 0.7], linewidth=0)
        if fill == True:

            for j in range(len(df_resamp["TWT"]) - 1):
                a = df_resamp["TWT"].tolist()[j]
                b = df_resamp["TWT"].tolist()[j + 1]
                c = df_resamp["Syn_scale_" + str(angle)].tolist()[j]
                c2 = df_resamp[synth_prefix + str(angle)].tolist()[j]
                d = df_resamp["Syn_scale_" + str(angle)].tolist()[j + 1]

                axhdl.fill_betweenx([a, b], [c, d], i, color=cmap(norm(c2)))

    axhdl.set_xticklabels(labels=angles)
    axhdl.xaxis.set_ticks(range(0, len(angles)))
    axhdl.set_xlim((-excursion, len(angles) + excursion))
    axhdl.xaxis.tick_top()
    axhdl.xaxis.set_label_position('top')
    axhdl.set_xlabel("Theta, deg")
    axhdl.invert_yaxis()

    return


def plot_vawig_mesh(axhdl, df_resamp, excursion, angles, synth_prefix = "Syn_", twt = "TWT", vminmax=[-5, 5], fill=True, colours="seismic_petrel"):

    if colours == "seismic_petrel":
        seismic_petrel = load_seismic_petrel()
        colours = seismic_petrel

    col_list = list(map(lambda x: synth_prefix + str(x), angles))
    print(col_list)
    max_amp = df_resamp[col_list].abs().values.max()

    print("max amplitude = %f" % max_amp)

    # [ntrc, nsamp] = data.shape

    # t = np.hstack([0, t, t.max()])
    cmap = plt.get_cmap(colours)

    syn_list = []
    for angle, i in zip(angles, range(0, len(angles))):
        df_resamp["Syn_scale_" + str(angle)] = excursion * df_resamp[synth_prefix + str(angle)] / max_amp + angle
        axhdl.plot(df_resamp["Syn_scale_" + str(angle)], df_resamp["TWT"], color='black', linewidth=0.5)
        #print(df_resamp[synth_prefix + str(angle)].max())
        syn_list.append(synth_prefix + str(angle))
    if fill == True:
        syn_vol = df_resamp[syn_list].to_numpy()

        cb1 = axhdl.pcolormesh(angles, df_resamp["TWT"].tolist(), syn_vol, cmap=colours,
                               norm=plt.Normalize(vmin=vminmax[0], vmax=vminmax[1]))

    axhdl.xaxis.tick_top()

    # norm = plt.Normalize(vmin = vminmax[0] + i, vmax = vminmax[1] + i)

    # axhdl.set_xticklabels(labels = theta1)
    # axhdl.xaxis.set_ticks(range(0, len(theta1)))
    # axhdl.set_xlim((-excursion, len(theta1) + excursion))
    # axhdl.xaxis.tick_top()
    # axhdl.xaxis.set_label_position('top')
    axhdl.set_xlabel("Theta, deg")
    axhdl.invert_yaxis()

    return cb1


def synthetic_plot(df_resamp, angles, wvlt_t, wvlt_amp, vp = "VP", vs = "VS", rhob = "RHOB", twt = "TWT", mdkb = "MDKB", synth_prefix = "Syn_", dt = 0.5, upper_ref_z = 0, lower_ref_z = 100, excursion=2, vminmax=[-10, 10], fill=True, style="mesh", colours = "seismic_petrel"):



    if colours == "seismic_petrel":
        seismic_petrel = load_seismic_petrel()
        colours = seismic_petrel


    # Calculate twt of sampled interfaces
    line1, line2, upper_ref_t, lower_ref_t = upper_lower_ref(df_resamp, upper_ref_z, lower_ref_z, angles, synth_prefix = "Syn_", mdkb = "MDKB", twt = "TWT")
    print("upper ref TWT = %f, lower ref TWT = %f" % (upper_ref_t, lower_ref_t))

    min_plot_time = upper_ref_t - 50
    max_plot_time = lower_ref_t + 50

    #   Plotting Display Parameters
    df_resamp = df_resamp.loc[(df_resamp["TWT"] > min_plot_time) & (df_resamp["TWT"] < max_plot_time)].dropna(how="all")

    img = plt.imshow(np.array([vminmax]), cmap=seismic_petrel, norm=plt.Normalize(vmin=vminmax[0], vmax=vminmax[1]))
    img.set_visible(False)

    #   Create the plot figure
    fig = plt.figure(figsize=(16, 12))
    fig.set_facecolor('white')

    #   Plot log curves in two-way time
    ax0a = fig.add_subplot(261)  # 161 1x6, 1st grid
    l_vp_dig, = ax0a.plot(df_resamp[vp], df_resamp[twt], 'k', lw=1)
    ax0a.set_ylim((min_plot_time, max_plot_time))
    ax0a.set_xlim(1000, 6000)
    ax0a.invert_yaxis()
    ax0a.set_ylabel('TWT (ms)')
    ax0a.xaxis.tick_top()
    ax0a.xaxis.set_label_position('top')
    ax0a.set_xlabel('Vp (m/s)')
    ax0a.axhline(upper_ref_t, c="b")
    ax0a.axhline(lower_ref_t, c="r")
    # ax0a.axhline(lyr_times[0, 0], color='blue', lw=2, alpha=0.5)
    # ax0a.axhline(lyr_times[0, 1], color='red', lw=2, alpha=0.5)
    ax0a.grid()

    ax0b = fig.add_subplot(262)  # 163 1x6, 2nd grid
    l_vs_dig, = ax0b.plot(df_resamp[vs], df_resamp[twt], 'k', lw=1)
    ax0b.set_ylim((min_plot_time, max_plot_time))
    ax0b.set_xlim((0, 2500))
    ax0b.invert_yaxis()
    ax0b.xaxis.tick_top()
    ax0b.xaxis.set_label_position('top')
    ax0b.set_xlabel('Vs (m/s)')
    ax0b.set_yticklabels('')
    ax0b.axhline(upper_ref_t, c="b")
    ax0b.axhline(lower_ref_t, c="r")
    # ax0b.axhline(lyr_times[0, 0], color='blue', lw=2, alpha=0.5)
    # ax0b.axhline(lyr_times[0, 1], color='red', lw=2, alpha=0.5)
    ax0b.grid()

    ax0c = fig.add_subplot(263)  # 163 1x6, 3rd grid
    l_rho_dig, = ax0c.plot(df_resamp[rhob], df_resamp[twt], 'k', lw=1)
    ax0c.set_ylim((min_plot_time, max_plot_time))
    ax0c.set_xlim((1, 3))
    ax0c.invert_yaxis()
    ax0c.xaxis.tick_top()
    ax0c.xaxis.set_label_position('top')
    ax0c.set_xlabel('RHOB, g/cc')
    ax0c.set_yticklabels('')

    ax0c.axhline(upper_ref_t, c="b")
    ax0c.axhline(lower_ref_t, c="r")
    # ax0c.axhline(lyr_times[0, 0], color='blue', lw=2, alpha=0.5)
    # ax0c.axhline(lyr_times[0, 1], color='red', lw=2, alpha=0.5)
    ax0c.grid()

    # syn_zoep_pp =
    #   Plot synthetic gather and model top & base interfaces in two-way time
    ax1 = fig.add_subplot(222)
    if style == "wiggle":
        plot_vawig_wiggle(ax1, df_resamp, excursion, angles, synth_prefix = "Syn_", twt = "TWT", vminmax=vminmax, fill=fill, colours = colours)
    elif style == "mesh":
        plot_vawig_mesh(ax1, df_resamp, excursion, angles, synth_prefix = "Syn_", twt = "TWT",  vminmax=vminmax, fill=fill, colours = colours)

    ax1.set_ylim((min_plot_time, max_plot_time))
    # l_int1, = ax1.plot(lyr_times[:, 0], color='blue', lw=2)
    # l_int2, = ax1.plot(lyr_times[:, 1], color='red', lw=2)

    # plt.legend([l_int1, l_int2], ['Interface 1', 'Interface 2'], loc=4)
    ax1.invert_yaxis()
    # label_str = 'Synthetic angle gather\nLayer 2 thickness = %4.1fm' % thickness
    # ax1.set_xlabel(label_str, fontsize=14)
    ax1.grid()
    ax1.set_ylabel('TWT (ms)')
    plt.colorbar(img, ax=ax1, orientation="vertical")

    #   Plot Zoeppritz and convolved reflectivity curves

    ax1.axhline(upper_ref_t, c="b")
    ax1.axhline(lower_ref_t, c="r")

    # ax2 = fig.add_subplot(2, 2, 3)

    ax3 = fig.add_subplot(2, 2, 4)
    l_syn1, = ax3.plot(line1[0], line1[1], color='blue', linewidth=2, label="Upper Interface Reflectivity")
    # l_rc1, = ax2.plot(rc_zoep_pp[:, 0], '--', color='blue', lw=2)
    # plt.legend([l_syn1, l_rc1], ['Convolved', 'Zoepprtiz'], loc=0)

    l_syn2, = ax3.plot(line2[0], line2[1], color='red', linewidth=2, label="Lower Interface Reflectivity")
    # l_rc2, = ax3.plot(rc_zoep_pp[:, 1], '--', color='red', lw=2)
    ax3.set_xlim(-excursion, len(angles) + excursion)
    # ax3.set_ylim(-1, 1)
    # ax3.set_ylim(-0.5, 0.5)
    ax3.grid()
    ax3.set_xlabel('Angle of incidence (deg)')
    ax3.set_ylabel('Reflection coefficient')
    ax3.set_title('Lower interface reflectivity')
    ax3.axhline(0, c="k")

    ax3.set_xticklabels(labels=angles)
    ax3.xaxis.set_ticks(range(0, len(angles)))
    ax3.xaxis.tick_top()
    ax3.legend()
    ax3.set_ylim(vminmax[0] - 5, vminmax[1] + 5)

    ax4 = fig.add_subplot(4, 2, 5)
    ax5 = fig.add_subplot(4, 2, 7)
    ## check input wavelet before/after resampling

    # Fs = The sampling frequency (samples per time unit)
    ax4.plot(wvlt_t, wvlt_amp, label=None)

    velocity = df_resamp[vp].mean()
    mag_spectrum, mag_freqs, mag_line = ax5.magnitude_spectrum(wvlt_amp, Fs=1 / (dt / 1000))
    i = np.ndarray.argmax(mag_spectrum)
    wavelength = velocity / mag_freqs[i]

    ax5.axvline(mag_freqs[i], c="k")
    ax5.annotate(s="Peak Frequency = %.2f Hz" % mag_freqs[i], xy=(mag_freqs[i] + 2, 0), xycoords="data")
    legend_wv = []
    legend_wv.append(mpatches.Rectangle((0, 0), 1, 1, fc="w", edgecolor="none", fill=False,
                                        label="Peak freq = %.2f Hz. At velocity = %.0f m/s, %s = %.0f m" % (
                                        mag_freqs[i], velocity, chr(955), wavelength)))
    legend_wv.append(mpatches.Rectangle((0, 0), 1, 1, fc="w", edgecolor="none", fill=False,
                                        label="Estimated tuning thickness = %s/4 = %.2fm" % (chr(955), wavelength / 4)))

    plt.legend(handles=legend_wv, loc=3, bbox_to_anchor=(0, -0.5), ncol=1)

    ax4.grid()
    ax4.set_xlabel("Time, ms")
    ax4.set_ylabel("Amplitude")
    ax5.set_xlim(0, 125)

    ax5.grid()

    # plt.legend([l_syn2, l_rc2], ['Convolved', 'Zoepprtiz'], loc=0)

    return fig
