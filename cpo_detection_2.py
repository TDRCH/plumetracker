import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import numpy as np
from mpl_toolkits.basemap import Basemap
import datetime
import shelve
import multiprocessing
import os
import copy
import os.path
from scipy import ndimage as ndi
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
import matplotlib.cm as cm
import scipy.signal
import tables
from copy import deepcopy
import cv2
from skimage import io, morphology, img_as_bool, segmentation
from scipy.ndimage import measurements

import utilities
import plumes
import pinkdust
import plotting
import get_llj_prob_model

def cloud_mask_mw(i, datetimes, oneday_datetimes):
    """
    Moving window cloud masking to be used by multiprocessing
    :return:
    """

    # Get lats and lons
    sdf_test = Dataset(
        '/soge-home/data_not_backed_up/satellite/meteosat/seviri'
        '/15-min/0.03x0.03/sdf/nc/JUNE2010/SDF_v2/SDF_v2.'
        '201006031500.nc')

    lons, lats = np.meshgrid(sdf_test.variables['longitude'][:],
                             sdf_test.variables['latitude'][:])

    date = oneday_datetimes[i]
    window_datetime_lower = datetime.datetime(datetimes[0].year,
                                              datetimes[0].month,
                                              datetimes[0].day,
                                              datetimes[0].hour,
                                              datetimes[0].minute) \
                            - datetime.timedelta(days=7)
    window_datetime_upper = datetime.datetime(datetimes[-1].year,
                                              datetimes[-1].month,
                                              datetimes[-1].day,
                                              datetimes[-1].hour,
                                              datetimes[-1].minute) \
                            + datetime.timedelta(days=7)

    # Get datetime objects between the above bounds
    time_params_7dayw = np.array([window_datetime_lower.year,
                                  window_datetime_upper.year,
                                  window_datetime_lower.month,
                                  window_datetime_upper.month,
                                  window_datetime_lower.day,
                                  window_datetime_upper.day,
                                  date.hour,
                                  date.hour,
                                  date.minute,
                                  date.minute])
    datetimes_7dayw = utilities.get_daily_datetime_objects(
        time_params_7dayw)

    bt_15day = np.zeros((datetimes_7dayw.shape[0], 3,
                         lats.shape[0],
                         lons.shape[1]))

    print str(oneday_datetimes[i].hour) + '_' + str(
            oneday_datetimes[i].minute)

    g = tables.open_file(
        '/soge-home/projects/seviri_dust/sdf/intermediary_files/bt_15d_' +
        str(oneday_datetimes[i].year) + '_' + str(oneday_datetimes[i].month) +
        '_' + str(oneday_datetimes[i].hour) +
        '_' + str(
            oneday_datetimes[i].minute) + '.hdf', 'w')
    atom = tables.Atom.from_dtype(bt_15day.dtype)
    filters = tables.Filters(complib='blosc', complevel=5)
    bts = g.create_carray(g.root, 'data', atom,
                         bt_15day.shape,
                         filters=filters)

    # Loop through each day of the time window for this time of day
    for j in np.arange(0, len(datetimes_7dayw)):
        date_w = datetimes_7dayw[j]

        # Extract BT data for this timestep
        filename = '/ouce-home/data/satellite/meteosat/seviri/15-min/' \
                   '0.03x0.03/bt/nc' \
                   '/' + str(date_w.strftime(
            '%B').upper()) + str(date_w.year) + \
                   '/H-000-MSG2__-MSG2________-' \
                   'IR_BrightnessTemperatures___-000005___-' + str(
            date_w.strftime('%Y')) + str(date_w.strftime('%m')) + \
                   str(date_w.strftime('%d')) + str(date_w.strftime('%H')) \
                   + str(date_w.strftime('%M')) + '-__.nc'
        if os.path.isfile(filename):
            btdata = Dataset(filename, 'r')
        else:
            filename = '/ouce-home/data/satellite/meteosat/seviri/15-min/' \
                       '0.03x0.03/bt/nc' \
                       '/' + str(date_w.strftime(
                '%B').upper()) + str(date_w.year) + \
                       '/H-000-MSG1__-MSG1________-' \
                       'IR_BrightnessTemperatures___-000005___-' + str(
                date_w.strftime('%Y')) + str(date_w.strftime('%m')) + \
                       str(date_w.strftime('%d')) + str(date_w.strftime('%H')) \
                       + str(date_w.strftime('%M')) + '-__.nc'
            if os.path.isfile(filename):
                btdata = Dataset(filename, 'r')
            else:
                print 'Found no BT data for ' + filename
                continue

        bt087 = btdata.variables['bt087'][:][0]
        bt108 = btdata.variables['bt108'][:][0]
        bt120 = btdata.variables['bt120'][:][0]

        bts[j, 0] = bt087
        bts[j, 1] = bt108
        bts[j, 2] = bt120
        btdata.close()

    # Save cloud masked data for this time of day to file
    g.close()

def detect_cpo(btdiff_2_anom_prev, btdiff_2_anom_prev_2,
               btdiff_2_anom_prev_3, datetimes, datestrings, date_i):

    # Get lats and lons
    sdf_test = Dataset(
        '/soge-home/data_not_backed_up/satellite/meteosat/seviri'
        '/15-min/0.03x0.03/sdf/nc/JUNE2010/SDF_v2/SDF_v2.'
        '201006031500.nc')

    lons, lats = np.meshgrid(sdf_test.variables['longitude'][:],
                             sdf_test.variables['latitude'][:])

    # Get cloud lats and lons
    cloud_test = Dataset(
            '/soge-home/data/satellite/meteosat/seviri/15-min/'
            '0.03x0.03/cloudmask'
            '/nc/'
            +
            'JUNE2010_CLOUDS/eumetsat.cloud.'
            + '201006031500.nc')
    cloud_lons = cloud_test.variables['lon'][:]
    cloud_lats = cloud_test.variables['lat'][:]
    cloud_lons, cloud_lats = np.meshgrid(cloud_lons, cloud_lats)
    lonmask = lons > 360
    latmask = lats > 90
    lons = np.ma.array(lons, mask=lonmask)
    lats = np.ma.array(lats, mask=latmask)
    used_ids = []

    runtime = datetimes[date_i] - datetimes[0]
    #print '\n' + datestrings[date_i] + '\n'
    totaltest = datetime.datetime.now()

    found_file = True

    if os.path.isfile('/ouce-home/data/satellite/meteosat/seviri/15-min/'
        '0.03x0.03/bt'
        '/nc/'
        +
        datetimes[date_i].strftime("%B").upper(
        ) + str(datetimes[date_i].year) + '/H-000-MSG2__'
                                          '-MSG2________-'
                                          'IR_BrightnessTemperatures___'
                                          '-000005___-'
        + datestrings[date_i] +
        '-__.nc'):
        bt = Dataset(
            '/ouce-home/data/satellite/meteosat/seviri/15-min/'
            '0.03x0.03/bt'
            '/nc/'
            +
            datetimes[date_i].strftime("%B").upper(
            ) + str(datetimes[date_i].year) + '/H-000-MSG2__'
                                              '-MSG2________-'
                                              'IR_BrightnessTemperatures___'
                                              '-000005___-'
            + datestrings[date_i] +
            '-__.nc')
    elif os.path.isfile('/ouce-home/data/satellite/meteosat/seviri/15-min/'
        '0.03x0.03/bt'
        '/nc/'
        +
        datetimes[date_i].strftime("%B").upper(
        ) + str(datetimes[date_i].year) + '/H-000-MSG1__'
                                          '-MSG1________-'
                                          'IR_BrightnessTemperatures___'
                                          '-000005___-'
        + datestrings[date_i] +
        '-__.nc'):
        bt = Dataset(
            '/ouce-home/data/satellite/meteosat/seviri/15-min/'
            '0.03x0.03/bt'
            '/nc/'
            +
            datetimes[date_i].strftime("%B").upper(
            ) + str(datetimes[date_i].year) + '/H-000-MSG1__'
                                              '-MSG1________-'
                                              'IR_BrightnessTemperatures___'
                                              '-000005___-'
            + datestrings[date_i] +
            '-__.nc')

    try:
        cloudmask = Dataset(
            '/soge-home/data/satellite/meteosat/seviri/15-min/'
            '0.03x0.03/cloudmask'
            '/nc/'
            +
            datetimes[date_i].strftime("%B").upper(
            ) + str(datetimes[date_i].year) + '_CLOUDS/eumetsat.cloud.'
            + datestrings[date_i] + '.nc')
        clouds_now = cloudmask.variables['cmask'][:][0]
        cloud_lons = cloudmask.variables['lon'][:]
        cloud_lats = cloudmask.variables['lat'][:]
    except:
        clouds_now = np.zeros(cloud_lons.shape)
        found_file = False

    if found_file:
        # Produce 12-10.8 imagery
        bt087 = bt.variables['bt087'][:][0]
        bt12 = bt.variables['bt120'][:][0]
        bt108 = bt.variables['bt108'][:][0]

        orig_lons = bt.variables['longitude'][:]
        orig_lats = bt.variables['latitude'][:]

        # print bt12.shape
        # print clouds_now.shape

        f = tables.open_file(
            '/soge-home/projects/seviri_dust/sdf/intermediary_files/bt_15d_' +
            str(datetimes[date_i].year) + '_' + str(
                datetimes[date_i].month) + '_' +
            str(datetimes[date_i].hour) + '_' + str(
                datetimes[date_i].minute) + '.hdf')
        arrobj = f.get_node('/data')
        bt_15day = arrobj.read()
        f.close()

        window_datetime_lower = datetime.datetime(datetimes[0].year,
                                                  datetimes[0].month,
                                                  datetimes[0].day,
                                                  datetimes[0].hour,
                                                  datetimes[0].minute) \
                                - datetime.timedelta(days=7)
        window_datetime_upper = datetime.datetime(datetimes[-1].year,
                                                  datetimes[-1].month,
                                                  datetimes[-1].day,
                                                  datetimes[-1].hour,
                                                  datetimes[-1].minute) \
                                + datetime.timedelta(days=7)

        BT_15_day_lower_bound = datetimes[date_i] - datetime.timedelta(
            days=7)
        BT_15_day_upper_bound = datetimes[date_i] + datetime.timedelta(
            days=7)

        # Get datetime objects between the above bounds
        time_params_7dayw = np.array([window_datetime_lower.year,
                                      window_datetime_upper.year,
                                      window_datetime_lower.month,
                                      window_datetime_upper.month,
                                      window_datetime_lower.day,
                                      window_datetime_upper.day,
                                      datetimes[date_i].hour,
                                      datetimes[date_i].hour,
                                      datetimes[date_i].minute,
                                      datetimes[date_i].minute])
        datetimes_7dayw = utilities.get_daily_datetime_objects(
            time_params_7dayw)

        BT_15_days = \
            bt_15day[np.asarray([j >= BT_15_day_lower_bound
                                 and j <= BT_15_day_upper_bound
                                 for j in datetimes_7dayw])]

        bt_15day_087 = BT_15_days[:, 0]
        #bt_15day_108 = BT_15_days[:, 1]
        bt_15day_120 = BT_15_days[:, 2]

        bt_15day_087_mean = np.nanmean(bt_15day_087, axis=0)
        #bt_15day_108_mean = np.nanmean(bt_15day_108, axis=0)
        bt_15day_120_mean = np.nanmean(bt_15day_120, axis=0)

        btdiff_2_15daymean = bt_15day_120_mean - bt_15day_087_mean

        orig_lons, orig_lats = np.meshgrid(orig_lons, orig_lats)

        bt087_regridded = pinkdust.regrid_data(orig_lons, orig_lats,
                                               cloud_lons,
                                               cloud_lats, bt087)

        bt12_regridded = pinkdust.regrid_data(orig_lons, orig_lats,
                                              cloud_lons,
                                              cloud_lats, bt12)

        btdiff_2_15daymean_regridded = pinkdust.regrid_data(orig_lons,
                                                            orig_lats,
                                                            cloud_lons,
                                                            cloud_lats,
                                                            btdiff_2_15daymean)

        btdiff_2 = bt12_regridded - bt087_regridded
        btdiff_2_anom = btdiff_2 - btdiff_2_15daymean_regridded

        btdiff_2_anom[clouds_now > 0] = np.nan

        if btdiff_2_anom_prev_3 != None:
            # Get the difference between this timestep and the one two
            # timesteps before
            btdiff_2_anom_diff = btdiff_2_anom - btdiff_2_anom_prev_3
            orig_btdiff_2_anom_diff = deepcopy(btdiff_2_anom_diff)
            btdiff_2_anom_diff += \
                orig_btdiff_2_anom_diff - btdiff_2_anom_prev_2
            btdiff_2_anom_diff += \
                orig_btdiff_2_anom_diff - btdiff_2_anom_prev
        else:
            btdiff_2_anom_diff = np.zeros((btdiff_2_anom.shape))
        if date_i == 0:
            btdiff_2_anom_prev = btdiff_2_anom
        elif date_i == 1:
            btdiff_2_anom_prev_2 = deepcopy(btdiff_2_anom_prev)
            btdiff_2_anom_prev = btdiff_2_anom
        elif date_i == 2:
            btdiff_2_anom_prev_3 = deepcopy(btdiff_2_anom_prev_2)
            btdiff_2_anom_prev_2 = deepcopy(btdiff_2_anom_prev)
            btdiff_2_anom_prev = deepcopy(btdiff_2_anom)
        elif date_i > 2:
            btdiff_2_anom_prev_3 = deepcopy(btdiff_2_anom_prev_2)
            btdiff_2_anom_prev_2 = deepcopy(btdiff_2_anom_prev)
            btdiff_2_anom_prev = deepcopy(btdiff_2_anom)

        lat_grad, lon_grad = np.gradient(btdiff_2_anom)
        total_grad = np.sqrt(lat_grad ** 2 + lon_grad ** 2)
        convolution = scipy.signal.convolve2d(clouds_now,
                                              np.ones((5,
                                                       5)),
                                              mode='same')
        clouds_now = convolution > 0
        total_grad[clouds_now == 1] = np.nan

        ### PASS I ###
        # In the FIRST PASS the LORD sayeth unto the image, 'Let all
        # whose BTD is below -10K be classified as CPOs, and remove the
        # tiny ones'
        # And those who fulfilled this condition were classified,
        # and it was good

        convolution = scipy.signal.convolve2d(clouds_now,
                                              np.ones((5,
                                                       5)),
                                              mode='same')
        clouds_now = convolution > 0

        btdiff_2_anom_diff_um = deepcopy(btdiff_2_anom_diff)
        btdiff_2_anom_diff[clouds_now > 0] = np.nan

        cpo_mask_pass_1 = btdiff_2_anom_diff < -7
        label_objects, nb_labels = ndi.label(cpo_mask_pass_1)
        sizes = np.bincount(label_objects.ravel())
        # Set clusters smaller than size 20 to zero
        mask_sizes = sizes > 20
        mask_sizes[0] = 0
        cpo_mask_pass_1 = mask_sizes[label_objects]

        ### PASS II ###
        # In the SECOND PASS the LORD sayeth unto the image, 'Let all
        # those included in the first pass which contain pixels which
        # are below -15K be classified'
        # And those who fulfilled this condition were classified,
        # and it was better

        cpo_mask_pass_2 = deepcopy(cpo_mask_pass_1)

        # Label the image and get all connected elements
        cpo_mask_pass_1, num = measurements.label(cpo_mask_pass_1)

        # Then loop through each labeled blob and find if one of the
        # pixels has a -15K. If it does, the whole blob is a freaking CPO.
        blob_ids = np.unique(cpo_mask_pass_1)
        blob_ids = blob_ids[blob_ids != 0]

        for i in np.arange(0, len(blob_ids)):
            target_region = cpo_mask_pass_1 == blob_ids[i]
            # First check if this region is within the original
            # time-gradient identified region (i.e. not one introduced
            # with the new generous gradient checker)
            if np.any(btdiff_2_anom_diff[target_region == 1] < -15):
                # Next check if there is a generous gradient contained
                # within this region somewhere
                pass
            else:
                cpo_mask_pass_2[target_region == 1] = 0

        # For identified CPO regions, undo the convolution on the cloud
        # mask

        cpo_mask_um = btdiff_2_anom_diff_um < -7
        # Label the image and get all connected elements
        cpo_mask_um, num = measurements.label(cpo_mask_um)

        # Then loop through each labeled blob and find if one of the
        # pixels has a -15K. If it does, the whole blob is a freaking CPO.
        blob_ids = np.unique(cpo_mask_um)
        blob_ids = blob_ids[blob_ids != 0]

        for i in np.arange(0, len(blob_ids)):
            target_region = cpo_mask_um == blob_ids[i]
            if np.any(cpo_mask_pass_2[target_region == 1] == 1):
                cpo_mask_pass_2[target_region == 1] = 1

        return cpo_mask_pass_2, btdiff_2_anom_prev, btdiff_2_anom_prev_2,\
               btdiff_2_anom_prev_3

if __name__ == '__main__':

    moving_window = False

    year_lower = 2010
    year_upper = 2010
    month_lower = 7
    month_upper = 7
    day_lower = 30
    day_upper = 31
    hour_lower = 18
    hour_upper = 23
    minute_lower = 0
    minute_upper = 45

    time_params_oneday = np.array([year_lower, year_lower, month_lower,
                                   month_lower, day_lower, day_lower,
                                   hour_lower, hour_upper,
                                   minute_lower, minute_upper])

    oneday_datetimes = utilities.get_datetime_objects(time_params_oneday)

    # Get lats and lons
    sdf_test = Dataset(
        '/soge-home/data_not_backed_up/satellite/meteosat/seviri'
        '/15-min/0.03x0.03/sdf/nc/JUNE2010/SDF_v2/SDF_v2.'
        '201006031500.nc')

    lons, lats = np.meshgrid(sdf_test.variables['longitude'][:],
                             sdf_test.variables['latitude'][:])

    # Here specify a set of years and months for which you will run.
    # Then the wrapper function is just all the below, taking month_lower,
    # month_upper, year_lower, year_upper as an argument and specifying the
    # other time bounds

    time_params = np.array([year_lower, year_upper, month_lower,
                            month_upper, day_lower, day_upper,
                            hour_lower, hour_upper, minute_lower,
                            minute_upper])

    datetimes = utilities.get_datetime_objects(time_params)
    datestrings = [j.strftime("%Y%m%d%H%M") for j in datetimes]

    # Get lats and lons
    sdf_test = Dataset(
        '/soge-home/data_not_backed_up/satellite/meteosat/seviri'
        '/15-min/0.03x0.03/sdf/nc/JUNE2010/SDF_v2/SDF_v2.'
        '201006031500.nc')

    lons, lats = np.meshgrid(sdf_test.variables['longitude'][:],
                             sdf_test.variables['latitude'][:])

    # Get cloud lats and lons
    cloud_test = Dataset(
            '/soge-home/data/satellite/meteosat/seviri/15-min/'
            '0.03x0.03/cloudmask'
            '/nc/'
            +
            'JUNE2010_CLOUDS/eumetsat.cloud.'
            + '201006031500.nc')
    cloud_lons = cloud_test.variables['lon'][:]
    cloud_lats = cloud_test.variables['lat'][:]
    cloud_lons, cloud_lats = np.meshgrid(cloud_lons, cloud_lats)
    lonmask = lons > 360
    latmask = lats > 90
    lons = np.ma.array(lons, mask=lonmask)
    lats = np.ma.array(lats, mask=latmask)
    used_ids = []
    found_file = True

    last_10_ids = []

    ### BACKGROUND PHASE ###

    if moving_window:
        print 'Generating cloud masks for all timesteps'
        # Loop through each time in a single day

        # Loop through blocks of 15 - what happens if we don't have 15 separate
        # dates to do?
        processes = [multiprocessing.Process(target=cloud_mask_mw, args=(i,
                                    datetimes, oneday_datetimes))
                     for i in np.arange(0, len(oneday_datetimes))]
        for p in processes:
            p.start()
        for p in processes:
            p.join()

    ### DETECTION PHASE ###
    print 'Detecting CPOs'

    btdiff_2_anom_prev = None
    btdiff_2_anom_prev_2 = None
    btdiff_2_anom_prev_3 = None

    for date_i in np.arange(0, len(datestrings)):
        runtime = datetimes[date_i] - datetimes[0]
        print '\n' + datestrings[date_i] + '\n'
        totaltest = datetime.datetime.now()

        found_file = True

        if os.path.isfile('/ouce-home/data/satellite/meteosat/seviri/15-min/'
                          '0.03x0.03/bt'
                          '/nc/'
                                  +
                                  datetimes[date_i].strftime("%B").upper(
                                  ) + str(
            datetimes[date_i].year) + '/H-000-MSG2__'
                                      '-MSG2________-'
                                      'IR_BrightnessTemperatures___'
                                      '-000005___-'
                                  + datestrings[date_i] +
                                  '-__.nc'):
            bt = Dataset(
                '/ouce-home/data/satellite/meteosat/seviri/15-min/'
                '0.03x0.03/bt'
                '/nc/'
                +
                datetimes[date_i].strftime("%B").upper(
                ) + str(datetimes[date_i].year) + '/H-000-MSG2__'
                                                  '-MSG2________-'
                                                  'IR_BrightnessTemperatures___'
                                                  '-000005___-'
                + datestrings[date_i] +
                '-__.nc')
            found_file = True
        elif os.path.isfile('/ouce-home/data/satellite/meteosat/seviri/15-min/'
                            '0.03x0.03/bt'
                            '/nc/'
                                    +
                                    datetimes[date_i].strftime("%B").upper(
                                    ) + str(
            datetimes[date_i].year) + '/H-000-MSG1__'
                                      '-MSG1________-'
                                      'IR_BrightnessTemperatures___'
                                      '-000005___-'
                                    + datestrings[date_i] +
                                    '-__.nc'):
            bt = Dataset(
                '/ouce-home/data/satellite/meteosat/seviri/15-min/'
                '0.03x0.03/bt'
                '/nc/'+datetimes[date_i].strftime("%B").upper() + str(
                    datetimes[date_i].year) + '/H-000-MSG1__-MSG1________-IR_BrightnessTemperatures___-000005___-'+ datestrings[date_i] + '-__.nc')
            found_file = True
        else:
            found_file = False
            continue

        try:
            cloudmask = Dataset(
                '/soge-home/data/satellite/meteosat/seviri/15-min/'
                '0.03x0.03/cloudmask'
                '/nc/'
                +
                datetimes[date_i].strftime("%B").upper(
                ) + str(datetimes[date_i].year) + '_CLOUDS/eumetsat.cloud.'
                + datestrings[date_i] + '.nc')
            clouds_now = cloudmask.variables['cmask'][:][0]
            cloud_lons = cloudmask.variables['lon'][:]
            cloud_lats = cloudmask.variables['lat'][:]
        except:
            clouds_now = np.zeros(cloud_lons.shape)
            continue

        if found_file:
            # Produce 12-10.8 imagery
            bt087 = bt.variables['bt087'][:][0]
            bt12 = bt.variables['bt120'][:][0]
            bt108 = bt.variables['bt108'][:][0]


            orig_lons = bt.variables['longitude'][:]
            orig_lats = bt.variables['latitude'][:]


            #print bt12.shape
            #print clouds_now.shape

            f = tables.open_file(
                '/soge-home/projects/seviri_dust/sdf/intermediary_files/bt_15d_' +
                str(datetimes[date_i].year) + '_' + str(
                    datetimes[date_i].month) + '_' +
                str(datetimes[date_i].hour) + '_' + str(
                    datetimes[date_i].minute) + '.hdf')
            arrobj = f.get_node('/data')
            bt_15day = arrobj.read()
            f.close()

            window_datetime_lower = datetime.datetime(year_lower, month_lower,
                                                      day_lower, hour_lower,
                                                      minute_lower) \
                                    - datetime.timedelta(days=7)
            window_datetime_upper = datetime.datetime(year_upper, month_upper,
                                                      day_upper, hour_upper,
                                                      minute_upper) \
                                    + datetime.timedelta(days=7)

            BT_15_day_lower_bound = datetimes[date_i] - datetime.timedelta(
                days=7)
            BT_15_day_upper_bound = datetimes[date_i] + datetime.timedelta(
                days=7)

            # Get datetime objects between the above bounds
            time_params_7dayw = np.array([window_datetime_lower.year,
                                          window_datetime_upper.year,
                                          window_datetime_lower.month,
                                          window_datetime_upper.month,
                                          window_datetime_lower.day,
                                          window_datetime_upper.day,
                                          datetimes[date_i].hour,
                                          datetimes[date_i].hour,
                                          datetimes[date_i].minute,
                                          datetimes[date_i].minute])
            datetimes_7dayw = utilities.get_daily_datetime_objects(
                time_params_7dayw)

            BT_15_days = \
                bt_15day[np.asarray([j >= BT_15_day_lower_bound
                                                  and j <= BT_15_day_upper_bound
                                                  for j in datetimes_7dayw])]

            bt_15day_087 = BT_15_days[:, 0]
            bt_15day_108 = BT_15_days[:, 1]
            bt_15day_120 = BT_15_days[:, 2]

            bt_15day_087_mean = np.nanmean(bt_15day_087, axis=0)
            bt_15day_108_mean = np.nanmean(bt_15day_108, axis=0)
            bt_15day_120_mean = np.nanmean(bt_15day_120, axis=0)

            btdiff_2_15daymean = bt_15day_120_mean - bt_15day_087_mean

            orig_lons, orig_lats = np.meshgrid(orig_lons, orig_lats)

            bt087_regridded  = pinkdust.regrid_data(orig_lons, orig_lats,
                                                  cloud_lons,
                                          cloud_lats, bt087)

            bt12_regridded = pinkdust.regrid_data(orig_lons, orig_lats,
                                                  cloud_lons,
                                          cloud_lats, bt12)

            btdiff_2_15daymean_regridded = pinkdust.regrid_data(orig_lons,
                                                               orig_lats,
                                                              cloud_lons,
                                                              cloud_lats,
                                                              btdiff_2_15daymean)

            btdiff_2 = bt12_regridded-bt087_regridded
            btdiff_2_anom = btdiff_2 - btdiff_2_15daymean_regridded

            if btdiff_2_anom_prev_3 != None:
                # Get the difference between this timestep and the one two
                # timesteps before
                # But if it takes the diff between a field with a cloud mask
                #  and one now, the cloud mask will always expand through
                # time - so say the first field has no cloud, second field
                # has tons. First field has no mask, second has tons. Diff
                # for second will have tons. Third has none. Diff for third
                # to second is tons, third to first is zero, so add them and
                #  third has tons. If there's then zero for the next three
                # timesteps

                # So the error is presumably from where the difference with
                # the previous is taken.

                # Question is, should it be possible to be rid of the cloud
                # in the first timestep
                # anom is screened and has loads of cloud
                # then in the next timestep there's none - diff field has loads
                # in the next timestep, still none - diff field has loads
                # next timestep, still none - diff field has loads

                btdiff_2_anom_diff = btdiff_2_anom-btdiff_2_anom_prev_3
                btdiff_2_anom_diff += \
                    (btdiff_2_anom-btdiff_2_anom_prev_2)
                btdiff_2_anom_diff += \
                    (btdiff_2_anom-btdiff_2_anom_prev)

            else:
                btdiff_2_anom_diff = np.zeros((btdiff_2_anom.shape))

            if date_i == 0:
                btdiff_2_anom_prev = deepcopy(btdiff_2_anom)
            elif date_i == 1:
                btdiff_2_anom_prev_2 = deepcopy(btdiff_2_anom_prev)
                btdiff_2_anom_prev = deepcopy(btdiff_2_anom)
            elif date_i == 2:
                btdiff_2_anom_prev_3 = deepcopy(btdiff_2_anom_prev_2)
                btdiff_2_anom_prev_2 = deepcopy(btdiff_2_anom_prev)
                btdiff_2_anom_prev = deepcopy(btdiff_2_anom)
            elif date_i > 2:
                btdiff_2_anom_prev_3 = deepcopy(btdiff_2_anom_prev_2)
                btdiff_2_anom_prev_2 = deepcopy(btdiff_2_anom_prev)
                btdiff_2_anom_prev = deepcopy(btdiff_2_anom)

            lat_grad, lon_grad = np.gradient(btdiff_2_anom)
            total_grad = np.sqrt(lat_grad ** 2 + lon_grad ** 2)
            convolution = scipy.signal.convolve2d(clouds_now,
                                                  np.ones((5,
                                                           5)),
                                                  mode='same')
            clouds_now = convolution > 0
            total_grad[clouds_now == 1] = np.nan

            plt.close()

            ### PASS I ###
            # In the FIRST PASS the LORD sayeth unto the image, 'Let all
            # whose BTD is below -10K be classified as CPOs, and remove the
            # tiny ones'
            # And those who fulfilled this condition were classified,
            # and it was good

            convolution = scipy.signal.convolve2d(clouds_now,
                                                  np.ones((5,
                                                           5)),
                                                  mode='same')
            clouds_now = convolution > 0

            btdiff_2_anom_diff_um = deepcopy(btdiff_2_anom_diff)
            btdiff_2_anom_diff[clouds_now > 0] = np.nan

            cpo_mask_pass_1 = btdiff_2_anom_diff < -7
            label_objects, nb_labels = ndi.label(cpo_mask_pass_1)
            sizes = np.bincount(label_objects.ravel())
            # Set clusters smaller than size 20 to zero
            mask_sizes = sizes > 20
            mask_sizes[0] = 0
            cpo_mask_pass_1 = mask_sizes[label_objects]

            ### PASS II ###
            # In the SECOND PASS the LORD sayeth unto the image, 'Let all
            # those included in the first pass which contain pixels which
            # are below -15K be classified'
            # And those who fulfilled this condition were classified,
            # and it was better

            cpo_mask_pass_2 = deepcopy(cpo_mask_pass_1)

            # Label the image and get all connected elements
            cpo_mask_pass_1, num = measurements.label(cpo_mask_pass_1)

            # Then loop through each labeled blob and find if one of the
            # pixels has a -15K. If it does, the whole blob is a freaking CPO.
            blob_ids = np.unique(cpo_mask_pass_1)
            blob_ids = blob_ids[blob_ids != 0]

            for i in np.arange(0, len(blob_ids)):
                target_region = cpo_mask_pass_1 == blob_ids[i]
                # First check if this region is within the original
                # time-gradient identified region (i.e. not one introduced
                # with the new generous gradient checker)
                if np.any(btdiff_2_anom_diff[target_region == 1] < -15):
                    # Next check if there is a generous gradient contained
                    # within this region somewhere
                    pass
                else:
                    cpo_mask_pass_2[target_region == 1] = 0

            # For identified CPO regions, undo the convolution on the cloud
            # mask

            cpo_mask_um = btdiff_2_anom_diff_um < -7
            # Label the image and get all connected elements
            cpo_mask_um, num = measurements.label(cpo_mask_um)

            # Then loop through each labeled blob and find if one of the
            # pixels has a -15K. If it does, the whole blob is a freaking CPO.
            blob_ids = np.unique(cpo_mask_um)
            blob_ids = blob_ids[blob_ids != 0]

            for i in np.arange(0, len(blob_ids)):
                target_region = cpo_mask_um == blob_ids[i]
                if np.any(cpo_mask_pass_2[target_region == 1] == 1):
                    cpo_mask_pass_2[target_region == 1] = 1

            plt.close()

            #cpo_mask_pass_2 = deepcopy(cpo_mask_pass_1)

            extent = (
                np.min(cloud_lons), np.max(cloud_lons), np.min(cloud_lats),
                np.max(
                    cloud_lats))
            m = Basemap(projection='cyl', llcrnrlon=extent[0],
                        urcrnrlon=extent[1],
                        llcrnrlat=extent[2], urcrnrlat=extent[3],
                        resolution='i')

            m.drawcoastlines(linewidth=0.5)
            m.drawcountries(linewidth=0.5)
            parallels = np.arange(10., 40, 2.)
            # labels = [left,right,top,bottom]
            m.drawparallels(parallels, labels=[False, True, True, False],
                            linewidth=0.5)
            meridians = np.arange(-20., 17., 4.)
            m.drawmeridians(meridians, labels=[True, False, False, True],
                            linewidth=0.5)

            min = 0
            max = 1

            levels = MaxNLocator(nbins=15).tick_values(min, max)

            cmap = cm.get_cmap('Blues_r')
            norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)

            m.pcolormesh(cloud_lons, cloud_lats, cpo_mask_pass_2,
                         cmap=cmap, vmin=min, vmax=max, norm=norm)

            cbar = plt.colorbar(orientation='horizontal', fraction=0.056,
                                pad=0.06)
            cbar.ax.set_xlabel('CPO mask pass 2')
            plt.tight_layout()
            plt.savefig('CPO_mask_pass_2_' + datestrings[
                date_i] + '.png',
                        bbox_inches='tight')

            plt.close()


            ### PASS IV ###

            # In the FOURTH PASS the LORD said 'let's shake things up a
            # little - I didn't just create static blobs - my CPOs are
            # LIVING, baby.' And so the LORD created dynamic CPO objects
            # and with these did he determine whether the CPOs had a long
            # enough duration to be detected, as well as whether they
            # originated close to deep convection and whether the centroid
            # had travelled far enough for detection
            # and it was Excellent

            # So how about a blob is allowed for a short duration if it merges
            #  with a bigger one

            # Put the SDF imagery through the plumetracker detection ID
            # assigner, but with a provision that it doesn't throw out
            # sub-250px ones

            # Basically we use exactly the same mechanics as plumetracker,
            # but with all this extra detection mumbo jumbo to start with (
            # which would be SDFs for plumetracker), and then we'll be able
            # to label all the properties, like distance from convection at
            # emission etc










