from segpy.reader import create_reader
from segpy.writer import write_segy
import numpy as np

#https://github.com/kwinkunks/notebooks/blob/master/Seismic_read_and_write.ipynb
#https://github.com/kwinkunks/ascii2segy/blob/master/ascii2segy.py

import segpy
from segpy import toolkit
from segpy import encoding


def write_stack_sgy_trace(sgy, outfile, output_sample_int = 2, header = None, il = 1000, xl = 1000, loc_x = 100000, loc_y = 1000000):
    dtype = 5  # IEE float
    segy_type = segpy.datatypes.DATA_SAMPLE_FORMAT_TO_SEG_Y_TYPE[dtype]

    if header == None:
        import datetime
        now = datetime.datetime.now()
        now_str = now.strftime("%Y-%m-%d %H:%M")
        header = ["Synthetic gather generated %s by J Wallis in segpy" % (now_str)]

    with open(outfile, 'wb') as fo:
        #n = create_reader(sgy)

        # text header (textural reel header
        info_header = header
        trh = [s[:80]+(80-len(s))*' ' for s in info_header]
        toolkit.write_textual_reel_header(fo, trh, encoding.EBCDIC) #ASCII

        # binary reel header
        brh = segpy.binary_reel_header.BinaryReelHeader()

        # MANDATORY FIELDS
        brh.data_traces_per_ensemble =1 # Pre-stack data
        brh.auxiliary_traces_per_ensemble = 1 # Pre-stack data
        brh.sample_interval = output_sample_int * 1000 #*1000 # microseconds? not millisec?
        brh.num_samples = sgy.shape[-1]
        brh.fixed_length_trace_flag = 1 # Default is 0 = varying trace length
        brh.format_revision_num = 0x100 # Default is 1
        brh.data_sample_format = dtype  # Default is 5 = IEEE float
        brh.num_extended_textual_headers = 0 # Default is 0
        # RECOMMENDED FIELDS
        brh.ensemble_fold = 0 #1 # Default = 0 offset=3227, default=0
        brh.trace_sorting = 4 # 2 = cdp, 4 = stack, default 0 = unknown COMMON_MIDPOINT = 8
        brh.measurement_system = 1 # m, default 0 = unknown, 2 = ft

        # Write the binary header.
        toolkit.write_binary_reel_header(fo, brh)

        # Pre-format trace header format.
        trace_header_packer = toolkit.make_header_packer(segpy.trace_header.TraceHeaderRev1)

        for i, xline, y in zip(range(0, len(xl)), xl, loc_y):
            samples = sgy[0, i, 0, :]
            #print(np.shape(samples))
            trace_header = segpy.trace_header.TraceHeaderRev1()

            trace_header.line_sequence_num = 1
            trace_header.file_sequence_num = 1

            trace_header.ensemble_num = 1 # cdp num --> 1 for synthetic, no spatial range
            trace_header.ensemble_trace_num = 1
            trace_header.trace_num = 1
            trace_header.group_x = loc_x
            trace_header.group_y = y
            trace_header.field_record_num = 1
            trace_header.num_samples = sgy.size
            trace_header.sample_interval = output_sample_int
            trace_header.doc_nreceiver_offset = 0

            trace_header.inline_number = il
            trace_header.crossline_number = xline
            trace_header.xy_scalar = 100
            trace_header.cdp_x = loc_x
            trace_header.cdp_y = y

            toolkit.write_trace_header(fo, trace_header, trace_header_packer)
            toolkit.write_trace_samples(fo, samples, seg_y_type=segy_type)

    return outfile

def write_gather_sgy_trace(sgy, outfile, theta1, output_sample_int = 2, header = None, il = 1000, xl = 1000, loc_x = 100000, loc_y = 1000000):
    dtype = 5  # IEE float
    segy_type = segpy.datatypes.DATA_SAMPLE_FORMAT_TO_SEG_Y_TYPE[dtype]

    if header == None:
        import datetime
        now = datetime.datetime.now()
        now_str = now.strftime("%Y-%m-%d %H:%M")
        header = ["Synthetic gather generated %s by J Wallis in segpy" % (now_str)]

    with open(outfile, 'wb') as fo:
        #n = create_reader(sgy)

        # text header (textural reel header
        info_header = header
        trh = [s[:80]+(80-len(s))*' ' for s in info_header]
        toolkit.write_textual_reel_header(fo, trh, encoding.EBCDIC) #ASCII

        # binary reel header
        brh = segpy.binary_reel_header.BinaryReelHeader()

        # MANDATORY FIELDS
        brh.data_traces_per_ensemble = len(theta1) # Pre-stack data
        brh.auxiliary_traces_per_ensemble = len(theta1) # Pre-stack data
        brh.sample_interval = output_sample_int * 1000 #*1000 # microseconds? not millisec?
        brh.num_samples = sgy.shape[-1]
        brh.fixed_length_trace_flag = 1 # Default is 0 = varying trace length
        brh.format_revision_num = 0x100 # Default is 1
        brh.data_sample_format = dtype  # Default is 5 = IEEE float
        brh.num_extended_textual_headers = 0 # Default is 0
        # RECOMMENDED FIELDS
        brh.ensemble_fold = len(theta1) #1 # Default = 0 offset=3227, default=0
        brh.trace_sorting = 2 # 2 = cdp, 4 = stack, default 0 = unknown COMMON_MIDPOINT = 8
        brh.measurement_system = 1 # m, default 0 = unknown, 2 = ft

        # Write the binary header.
        toolkit.write_binary_reel_header(fo, brh)

        # Pre-format trace header format.
        trace_header_packer = toolkit.make_header_packer(segpy.trace_header.TraceHeaderRev1)

        i = 0 # counting no. traces
        for angle, angle_no in zip(theta1, range(0, len(theta1))):
            i += 1
            #pint (i)
            samples = sgy[0, 0, angle_no, :]
            #print(np.shape(samples))
            trace_header = segpy.trace_header.TraceHeaderRev1()

            trace_header.line_sequence_num = i
            trace_header.file_sequence_num = i

            trace_header.ensemble_num = 1 # cdp num --> 1 for synthetic, no spatial range
            trace_header.ensemble_trace_num = angle_no
            trace_header.trace_num = i
            trace_header.group_x = loc_x
            trace_header.group_y = loc_y
            trace_header.field_record_num = i
            trace_header.num_samples = sgy.size
            trace_header.sample_interval = output_sample_int
            trace_header.doc_nreceiver_offset = angle

            trace_header.inline_number = il
            trace_header.crossline_number = xl
            trace_header.xy_scalar = 100
            trace_header.cdp_x = loc_x
            trace_header.cdp_y = loc_y

            toolkit.write_trace_header(fo, trace_header, trace_header_packer)
            toolkit.write_trace_samples(fo, samples, seg_y_type=segy_type)

    return outfile


def write_gather_sgy_wedge(sgy, outfile, angles,  dz_min = 0, dz_max = 100, dz_step = 10, theta_min = 0, theta_max = 45, theta_step = 5, output_sample_int = 2, header = None, loc_x = 1000, loc_y = None, il =  1000, xl = None):

    nmodel = int((dz_max - dz_min) / dz_step + 1)
    thickness = np.linspace(dz_min, dz_max, nmodel)

    nangles = int((theta_max - theta_min) / theta_step + 1)

    angles = []
    for i in range(0, nangles):
        theta = i * theta_step + theta_min
        angles.append(theta)

    if loc_y is None:
        loc_y = thickness * 1000

    if xl is None:
        xl = thickness
        xl = [int(x) for x in xl]


    i_size = 1  # segy_dataset.num_inlines()
    x_size = len(sgy[0, :, 0, 0])  # thicknesses/no models, segy_dataset.num_xlines()
    t_size = len(sgy[0, 0, 0, :])  # trace length/no. time samples, segy_dataset.max_num_trace_samples()
    a_size = len(angles)

    xl = [int(x) for x in thickness]

    if header == None:
        import datetime
        now = datetime.datetime.now()
        now_str = now.strftime("%Y-%m-%d %H:%M")
        header = ["Synthetic wedge/2D line generated %s by J Wallis in segpy" % (now_str)]

    dtype = 5  # IEE float
    segy_type = segpy.datatypes.DATA_SAMPLE_FORMAT_TO_SEG_Y_TYPE[dtype]


    with open(outfile, 'wb') as fo:
        # n = create_reader(sgy)

        # text header (textural reel header
        info_header = header
        trh = [s[:80] + (80 - len(s)) * ' ' for s in info_header]
        toolkit.write_textual_reel_header(fo, trh, encoding.EBCDIC)  # ASCII

        # binary reel header
        brh = segpy.binary_reel_header.BinaryReelHeader()

        # MANDATORY FIELDS
        brh.data_traces_per_ensemble = len(angles)  # Pre-stack data
        brh.auxiliary_traces_per_ensemble = len(angles)  # Pre-stack data
        brh.sample_interval = output_sample_int * 1000  # *1000 # microseconds? not millisec?
        brh.num_samples = sgy.shape[-1]
        brh.fixed_length_trace_flag = 1  # Default is 0 = varying trace length
        brh.format_revision_num = 0x100  # Default is 1
        brh.data_sample_format = dtype  # Default is 5 = IEEE float
        brh.num_extended_textual_headers = 0  # Default is 0
        # RECOMMENDED FIELDS
        brh.ensemble_fold = len(angles)  # 1 # Default = 0 offset=3227, default=0
        brh.trace_sorting = 2  # 2 = cdp, 4 = stack, default 0 = unknown COMMON_MIDPOINT = 8
        brh.measurement_system = 1  # m, default 0 = unknown, 2 = ft

        # Write the binary header.
        toolkit.write_binary_reel_header(fo, brh)

        # Pre-format trace header format.
        trace_header_packer = toolkit.make_header_packer(segpy.trace_header.TraceHeaderRev1)

        i = 0  # counting no. traces
        for n, thick, xline, y in zip(range(0, nmodel), thickness, xl, loc_y):
            # print (n, thick, xline, y) # = n = 0-11, index for thickness
            #print("X: %f, Y: %f, IL: %i, XL: %i, nmodel: %i" % (loc_x, y, il, xline, n))
            trace_header = segpy.trace_header.TraceHeaderRev1()
            trace_header.group_x = loc_x
            trace_header.group_y = loc_y[n]

            for angle, angle_no in zip(angles, range(0, len(angles))):
                i += 1
                #print("angle: %i" % angle)
                samples = sgy[0, n, angle_no, :]
                # print (max(samples))
                # print (max(samples))

                trace_header.line_sequence_num = i
                trace_header.file_sequence_num = i

                trace_header.ensemble_num = 1  # cdp num --> 1 for synthetic, no spatial range
                trace_header.ensemble_trace_num = angle_no
                trace_header.trace_num = i

                trace_header.field_record_num = i
                trace_header.num_samples = t_size  # in a full stack trace = no. samples in whole file
                trace_header.sample_interval = output_sample_int
                trace_header.doc_nreceiver_offset = angle

                trace_header.inline_number = il
                trace_header.crossline_number = xline
                trace_header.xy_scalar = 100
                trace_header.cdp_x = loc_x
                trace_header.cdp_y = y
                #print(loc_x, y, angle)

                toolkit.write_trace_header(fo, trace_header, trace_header_packer)
                toolkit.write_trace_samples(fo, samples, seg_y_type=segy_type)

    return outfile


def write_gather_sgy_2d(sgy, outfile, xl_list, theta_min = 0, theta_max = 45, theta_step = 5, output_sample_int = 2, header = None, loc_x = 1000, loc_y = None, il =  1000, xl = None):

    nmodel = len(xl_list) #int((dz_max - dz_min) / dz_step + 1)
    #thickness = np.linspace(dz_min, dz_max, nmodel)

    nangles = int((theta_max - theta_min) / theta_step + 1)

    angles = []
    for i in range(0, nangles):
        theta = i * theta_step + theta_min
        angles.append(theta)


    i_size = 1  # segy_dataset.num_inlines()
    x_size = len(sgy[0, :, 0, 0])  # thicknesses/no models, segy_dataset.num_xlines()
    t_size = len(sgy[0, 0, 0, :])  # trace length/no. time samples, segy_dataset.max_num_trace_samples()
    a_size = len(angles)

    xl_list = [int(x) for x in xl_list]

    if header == None:
        import datetime
        now = datetime.datetime.now()
        now_str = now.strftime("%Y-%m-%d %H:%M")
        header = ["Synthetic wedge/2D line generated %s by J Wallis in segpy" % (now_str)]

    dtype = 5  # IEE float
    segy_type = segpy.datatypes.DATA_SAMPLE_FORMAT_TO_SEG_Y_TYPE[dtype]


    with open(outfile, 'wb') as fo:
        # n = create_reader(sgy)

        # text header (textural reel header
        info_header = header
        trh = [s[:80] + (80 - len(s)) * ' ' for s in info_header]
        toolkit.write_textual_reel_header(fo, trh, encoding.EBCDIC)  # ASCII

        # binary reel header
        brh = segpy.binary_reel_header.BinaryReelHeader()

        # MANDATORY FIELDS
        brh.data_traces_per_ensemble = len(angles)  # Pre-stack data
        brh.auxiliary_traces_per_ensemble = len(angles)  # Pre-stack data
        brh.sample_interval = output_sample_int * 1000  # *1000 # microseconds? not millisec?
        brh.num_samples = sgy.shape[-1]
        brh.fixed_length_trace_flag = 1  # Default is 0 = varying trace length
        brh.format_revision_num = 0x100  # Default is 1
        brh.data_sample_format = dtype  # Default is 5 = IEEE float
        brh.num_extended_textual_headers = 0  # Default is 0
        # RECOMMENDED FIELDS
        brh.ensemble_fold = len(angles)  # 1 # Default = 0 offset=3227, default=0
        brh.trace_sorting = 2  # 2 = cdp, 4 = stack, default 0 = unknown COMMON_MIDPOINT = 8
        brh.measurement_system = 1  # m, default 0 = unknown, 2 = ft

        # Write the binary header.
        toolkit.write_binary_reel_header(fo, brh)

        # Pre-format trace header format.
        trace_header_packer = toolkit.make_header_packer(segpy.trace_header.TraceHeaderRev1)

        i = 0  # counting no. traces
        for n, xline, y in zip(range(0, nmodel), xl_list, loc_y):
            # print (n, thick, xline, y) # = n = 0-11, index for thickness
            #print (loc_x, y, il, xline, n)
            #print("X: %f, Y: %f, IL: %i, XL: %i, nmodel: %i" % (loc_x, y, il, xline, n))
            trace_header = segpy.trace_header.TraceHeaderRev1()
            trace_header.group_x = loc_x
            trace_header.group_y = loc_y[n]

            for angle, angle_no in zip(angles, range(0, len(angles))):
                i += 1
                #print("angle: %i" % angle)
                samples = sgy[0, n, angle_no, :]
                # print (max(samples))
                # print (max(samples))

                trace_header.line_sequence_num = i
                trace_header.file_sequence_num = i

                trace_header.ensemble_num = 1  # cdp num --> 1 for synthetic, no spatial range
                trace_header.ensemble_trace_num = angle_no
                trace_header.trace_num = i

                trace_header.field_record_num = i
                trace_header.num_samples = t_size  # in a full stack trace = no. samples in whole file
                trace_header.sample_interval = output_sample_int
                trace_header.doc_nreceiver_offset = angle

                trace_header.inline_number = il
                trace_header.crossline_number = xline
                trace_header.xy_scalar = 100
                trace_header.cdp_x = loc_x
                trace_header.cdp_y = y
                #print(loc_x, y, angle)

                toolkit.write_trace_header(fo, trace_header, trace_header_packer)
                toolkit.write_trace_samples(fo, samples, seg_y_type=segy_type)

    return outfile



def write_stack_sgy_wedge(sgy, outfile,  dz_min = 0, dz_max = 100, dz_step = 10, output_sample_int = 2, header = None, loc_x = 1000, loc_y = None, il =  1000, xl = None):

    nmodel = int((dz_max - dz_min) / dz_step + 1)
    thickness = np.linspace(dz_min, dz_max, nmodel)

    if loc_y is None:
        loc_y = thickness * 1000

    if xl is None:
        xl = thickness
        xl = [int(x) for x in xl]


    i_size = 1  # segy_dataset.num_inlines()
    x_size = len(sgy[0, :, 0, 0])  # thicknesses/no models, segy_dataset.num_xlines()
    t_size = len(sgy[0, 0, 0, :])  # trace length/no. time samples, segy_dataset.max_num_trace_samples()
    a_size = 1

    xl = [int(x) for x in thickness]

    if header == None:
        import datetime
        now = datetime.datetime.now()
        now_str = now.strftime("%Y-%m-%d %H:%M")
        header = ["Synthetic wedge/2D line generated %s by J Wallis in segpy" % (now_str)]

    dtype = 5  # IEE float
    segy_type = segpy.datatypes.DATA_SAMPLE_FORMAT_TO_SEG_Y_TYPE[dtype]


    with open(outfile, 'wb') as fo:
        # n = create_reader(sgy)

        # text header (textural reel header
        info_header = header
        trh = [s[:80] + (80 - len(s)) * ' ' for s in info_header]
        toolkit.write_textual_reel_header(fo, trh, encoding.EBCDIC)  # ASCII

        # binary reel header
        brh = segpy.binary_reel_header.BinaryReelHeader()

        # MANDATORY FIELDS
        brh.data_traces_per_ensemble = 1  # Pre-stack data
        brh.auxiliary_traces_per_ensemble = 1  # Pre-stack data
        brh.sample_interval = output_sample_int * 1000  # *1000 # microseconds? not millisec?
        brh.num_samples = sgy.shape[-1]
        brh.fixed_length_trace_flag = 1  # Default is 0 = varying trace length
        brh.format_revision_num = 0x100  # Default is 1
        brh.data_sample_format = dtype  # Default is 5 = IEEE float
        brh.num_extended_textual_headers = 0  # Default is 0
        # RECOMMENDED FIELDS
        brh.ensemble_fold = 0  # 1 # Default = 0 offset=3227, default=0
        brh.trace_sorting = 4  # 2 = cdp, 4 = stack, default 0 = unknown COMMON_MIDPOINT = 8
        brh.measurement_system = 1  # m, default 0 = unknown, 2 = ft

        # Write the binary header.
        toolkit.write_binary_reel_header(fo, brh)

        # Pre-format trace header format.
        trace_header_packer = toolkit.make_header_packer(segpy.trace_header.TraceHeaderRev1)

        i = 0  # counting no. traces
        for n, thick, xline, y in zip(range(0, nmodel), thickness, xl, loc_y):
            #print (n, thick, xline, y) # = n = 0-11, index for thickness
            #print("X: %f, Y: %f, IL: %i, XL: %i, nmodel: %i" % (loc_x, y, il, xline, n))
            trace_header = segpy.trace_header.TraceHeaderRev1()
            trace_header.group_x = loc_x
            trace_header.group_y = loc_y[n]


            samples = sgy[0, n, 0, :]
            # print (max(samples))
            # print (max(samples))

            trace_header.line_sequence_num = i
            trace_header.file_sequence_num = i

            trace_header.ensemble_num = 1  # cdp num --> 1 for synthetic, no spatial range
            trace_header.ensemble_trace_num = 0
            trace_header.trace_num = i

            trace_header.field_record_num = i
            trace_header.num_samples = t_size  # in a full stack trace = no. samples in whole file
            trace_header.sample_interval = output_sample_int
            trace_header.doc_nreceiver_offset = 0

            trace_header.inline_number = il
            trace_header.crossline_number = xline
            trace_header.xy_scalar = 100
            trace_header.cdp_x = il
            trace_header.cdp_y = xline
            #print(loc_x, y, angle)

            toolkit.write_trace_header(fo, trace_header, trace_header_packer)
            toolkit.write_trace_samples(fo, samples, seg_y_type=segy_type)

    return outfile


def write_stack_sgy_2d(sgy, outfile, xl_list, output_sample_int = 2, header = None, loc_x = 1000, loc_y = None, il =  1000, xl = None):

    nmodel = len(xl_list) #int((dz_max - dz_min) / dz_step + 1)
    #thickness = np.linspace(dz_min, dz_max, nmodel)

    i_size = 1  # segy_dataset.num_inlines()
    x_size = len(sgy[0, :, 0, 0])  # thicknesses/no models, segy_dataset.num_xlines()
    t_size = len(sgy[0, 0, 0, :])  # trace length/no. time samples, segy_dataset.max_num_trace_samples()
    a_size = 1

    xl_list = [int(x) for x in xl_list]

    if header == None:
        import datetime
        now = datetime.datetime.now()
        now_str = now.strftime("%Y-%m-%d %H:%M")
        header = ["Synthetic wedge/2D line generated %s by J Wallis in segpy" % (now_str)]

    dtype = 5  # IEE float
    segy_type = segpy.datatypes.DATA_SAMPLE_FORMAT_TO_SEG_Y_TYPE[dtype]


    with open(outfile, 'wb') as fo:
        # n = create_reader(sgy)

        # text header (textural reel header
        info_header = header
        trh = [s[:80] + (80 - len(s)) * ' ' for s in info_header]
        toolkit.write_textual_reel_header(fo, trh, encoding.EBCDIC)  # ASCII

        # binary reel header
        brh = segpy.binary_reel_header.BinaryReelHeader()

        # MANDATORY FIELDS
        brh.data_traces_per_ensemble = 1  # Pre-stack data
        brh.auxiliary_traces_per_ensemble = 1 # Pre-stack data
        brh.sample_interval = output_sample_int * 1000  # *1000 # microseconds? not millisec?
        brh.num_samples = sgy.shape[-1]
        brh.fixed_length_trace_flag = 1  # Default is 0 = varying trace length
        brh.format_revision_num = 0x100  # Default is 1
        brh.data_sample_format = dtype  # Default is 5 = IEEE float
        brh.num_extended_textual_headers = 0  # Default is 0
        # RECOMMENDED FIELDS
        brh.ensemble_fold = 1  # 1 # Default = 0 offset=3227, default=0
        brh.trace_sorting = 4  # 2 = cdp, 4 = stack, default 0 = unknown COMMON_MIDPOINT = 8
        brh.measurement_system = 1  # m, default 0 = unknown, 2 = ft

        # Write the binary header.
        toolkit.write_binary_reel_header(fo, brh)

        # Pre-format trace header format.
        trace_header_packer = toolkit.make_header_packer(segpy.trace_header.TraceHeaderRev1)

        i = 0  # counting no. traces
        for n, xline, y in zip(range(0, nmodel), xl_list, loc_y):
            # print (n, thick, xline, y) # = n = 0-11, index for thickness
            #print (loc_x, y, il, xline, n)
            #print("X: %f, Y: %f, IL: %i, XL: %i, nmodel: %i" % (loc_x, y, il, xline, n))
            trace_header = segpy.trace_header.TraceHeaderRev1()
            trace_header.group_x = loc_x
            trace_header.group_y = loc_y[n]

            samples = sgy[0, n, 0, :]
            # print (max(samples))
            # print (max(samples))

            trace_header.line_sequence_num = 1
            trace_header.file_sequence_num = 1

            trace_header.ensemble_num = 1  # cdp num --> 1 for synthetic, no spatial range
            trace_header.ensemble_trace_num = 1
            trace_header.trace_num = i

            trace_header.field_record_num = i
            trace_header.num_samples = t_size  # in a full stack trace = no. samples in whole file
            trace_header.sample_interval = output_sample_int
            trace_header.doc_nreceiver_offset = 0

            trace_header.inline_number = il
            trace_header.crossline_number = xline
            trace_header.xy_scalar = 100
            trace_header.cdp_x = loc_x
            trace_header.cdp_y = y
            #print(loc_x, y, angle)

            toolkit.write_trace_header(fo, trace_header, trace_header_packer)
            toolkit.write_trace_samples(fo, samples, seg_y_type=segy_type)

    return outfile
