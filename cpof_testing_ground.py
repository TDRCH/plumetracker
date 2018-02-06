import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from netCDF4 import num2date
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
from pyresample import utils
from netCDF4 import date2num
from skimage.morphology import skeletonize, skeletonize_3d

import utilities
import plumes
import pinkdust
import plotting
import get_llj_prob_model

def plot(array, min=None, max=None, title=None):
    """

    :param array:
    :return:
    """

    if title == None:
        title = 'CPO_detection'

    extent = (
        np.min(orig_lons), np.max(orig_lons), np.min(orig_lats),
        np.max(
            orig_lats))

    m = Basemap(projection='cyl', llcrnrlon=extent[0],
                urcrnrlon=extent[1],
                llcrnrlat=extent[2], urcrnrlat=extent[3],
                resolution='i')

    m.drawcoastlines(linewidth=0.5)
    m.drawcountries(linewidth=0.5)
    parallels = np.arange(10., 40, 2.)
    m.drawparallels(parallels, labels=[False, True, True, False],
                    linewidth=0.5)
    meridians = np.arange(-20., 17., 4.)
    m.drawmeridians(meridians, labels=[True, False, False, True],
                    linewidth=0.5)

    if min != None:
        levels = MaxNLocator(nbins=15).tick_values(min, max)
        vmin = min
        vmax = max

    cmap = cm.get_cmap('Blues_r')

    if min != None:
        norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
        plot = m.pcolormesh(orig_lons, orig_lats, array,
                            cmap=cmap, vmin=min, vmax=max, norm=norm)
    else:
        plot = m.pcolormesh(orig_lons, orig_lats, array,
                            cmap=cmap)

    cbar = plt.colorbar(orientation='horizontal',
                        fraction=0.056,
                        pad=0.06)

    #cbar.ax.set_xlabel(labels[i])
    plt.tight_layout()
    plt.savefig(title + '_' + datestrings[
        date_i] + '.png',
                bbox_inches='tight')

    plt.close()

#def detect_speed(arr):
    # First get skeletonized arrays from each of the timesteps
    # Then mark the nearest point between the new skeleton and the old one

def detect_cpo(arr):

    # There will often be times when the thing breaks up, but because the
    # main signal is another segment, the breakaway gets undetected. So,
    # we should allow for detection in places where there has been a
    # detection in the previous timestep, as long as it's still over the
    # threshold. So the question is, what is lost between pass 1 and pass 2?

    threshold_1 = -6
    threshold_2 = -20

    arr_um = deepcopy(arr)

    cpo_mask_pass_1 = arr < threshold_1
    label_objects, nb_labels = ndi.label(cpo_mask_pass_1)
    sizes = np.bincount(label_objects.ravel())
    mask_sizes = sizes > 20
    mask_sizes[0] = 0
    cpo_mask_pass_1 = mask_sizes[label_objects]

    cpo_mask_p1 = deepcopy(cpo_mask_pass_1)

    cpo_mask_pass_2 = deepcopy(cpo_mask_pass_1)
    cpo_mask_pass_1, num = measurements.label(cpo_mask_pass_1)
    blob_ids = np.unique(cpo_mask_pass_1)
    blob_ids = blob_ids[blob_ids != 0]

    for i in np.arange(0, len(blob_ids)):
        target_region = cpo_mask_pass_1 == blob_ids[i]
        # First check if this region is within the original
        # time-gradient identified region (i.e. not one introduced
        # with the new generous gradient checker)
        if np.any(arr[target_region == 1] < threshold_2):
            # Next check if there is a generous gradient contained
            # within this region somewhere
            pass
        else:
            cpo_mask_pass_2[target_region == 1] = 0

    cpo_mask_um = arr_um < threshold_1
    # Label the image and get all connected elements
    cpo_mask_um, num = measurements.label(cpo_mask_um)

    # Then loop through each labeled blob and find if one of the
    # pixels has a -15K. If it does, the whole blob is a freaking CPO.
    blob_ids = np.unique(cpo_mask_um)
    blob_ids = blob_ids[blob_ids != 0]

    if 1 in cpo_mask_pass_2:
        for i in np.arange(0, len(blob_ids)):
            target_region = cpo_mask_um == blob_ids[i]
            if 1 in cpo_mask_pass_2[target_region]:
                # if np.any(cpo_mask_pass_2[target_region == 1] == 1):
                cpo_mask_pass_2[target_region] = 1

    return cpo_mask_pass_2, cpo_mask_p1

def wrapper(datetimes, datestrings, date_i, lons,
               lats,
               cloud_lons, cloud_lats,
               daily_clouds=False, double_digits=False, mesh=False,
               daily_bt=False):

    used_ids = []

    runtime = datetimes[date_i] - datetimes[0]
    # print '\n' + datestrings[date_i] + '\n'
    totaltest = datetime.datetime.now()

    found_file = True

    if daily_bt == False:
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
                '/nc/'
                +
                datetimes[date_i].strftime("%B").upper(
                ) + str(datetimes[date_i].year) + '/H-000-MSG1__'
                                                  '-MSG1________-'
                                                  'IR_BrightnessTemperatures___'
                                                  '-000005___-'
                + datestrings[date_i] +
                '-__.nc')
            found_file = True
        else:
            found_file = False

    if daily_clouds:
        try:
            cloudmask = Dataset('/soge-home/projects/seviri_dust/raw_seviri_'
                                'data/cloudmask_nc/' + datetimes[
                                    date_i].strftime("%B%Y") + '/cloudmask_' +
                                datetimes[
                                    date_i].strftime("%Y%m%d") + '.nc')
            cloudmask_times = num2date(cloudmask.variables['time'][:],
                                       cloudmask.variables['time'].units)
            cloudmask_times = np.asarray([datetime.datetime(j.year, j.month,
                                                            j.day, j.hour,
                                                            j.minute) for j
                                          in cloudmask_times])

            cloudmask_bool = cloudmask_times == datetimes[date_i]
            clouds_now = cloudmask.variables['cloud_mask'][cloudmask_bool][0]
            found_file = True
        except:
            print 'Found no cloud mask file!'
            clouds_now = np.zeros(cloud_lons.shape)
            found_file = False

        sdf_root = '/soge-home/projects/seviri_dust/sdf/' \
                   + datetimes[date_i].strftime('%B') \
                   + str(datetimes[date_i].year) + '/'

        if os.path.isfile(sdf_root + 'SDF_v2.' + \
                                  datestrings[date_i] + '.nc'):

            sdf = Dataset(
                sdf_root + 'SDF_v2.' + \
                datestrings[date_i] + '.nc')
            found_file = True
            # print sdf
            sdf_now = sdf.variables['SDF'][:]
        else:
            print 'No SDF file found for this date'
            found_file = False

    else:
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
            print 'Found no cloud mask file!'

        sdf_root = '/soge-home/data_not_backed_up/satellite/meteosat' \
                   '/seviri/15' \
                   '-min/0.03x0.03/sdf/nc/' + datetimes[date_i].strftime(
            '%B').upper() + str(datetimes[date_i].year) + '/SDF_v2/'
        if os.path.isfile(sdf_root + 'SDF_v2.' + \
                                  datestrings[date_i] + '.nc'):

            sdf = Dataset(
                sdf_root + 'SDF_v2.' + \
                datestrings[date_i] + '.nc')
            found_file = True
            # print sdf
            if 'time' in sdf.variables:
                sdf_now = sdf.variables['bt108'][0]
            else:
                sdf_now = sdf.variables['bt108'][:]
        else:
            print 'No SDF file found for this date'
            found_file = False


    if found_file:

        if daily_bt == False:
            bt087 = bt.variables['bt087'][:][0]
            bt108 = bt.variables['bt108'][:][0]
            bt12 = bt.variables['bt120'][:][0]
            orig_lons = bt.variables['longitude'][:]
            orig_lats = bt.variables['latitude'][:]
            orig_lons, orig_lats = np.meshgrid(orig_lons, orig_lats)
        else:
            orig_lons = lons
            orig_lats = lats

        # print bt12.shape
        # print clouds_now.shape

        window_datetime_lower = datetime.datetime(datetimes[0].year,
                                                  6,
                                                  1,
                                                  0,
                                                  45) \
                                - datetime.timedelta(days=7)
        window_datetime_upper = datetime.datetime(datetimes[-1].year,
                                                  8,
                                                  31,
                                                  23,
                                                  45) \
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

        indices = np.arange(0, len(datetimes_7dayw))
        lower_ind = datetimes_7dayw == BT_15_day_lower_bound
        lower_ind = indices[lower_ind][0]
        upper_ind = datetimes_7dayw == BT_15_day_upper_bound
        upper_ind = indices[upper_ind][0]
        current_ind = datetimes_7dayw == datetimes[date_i]
        current_ind = indices[current_ind][0]

        if double_digits:
            f = tables.open_file(
                '/soge-home/projects/seviri_dust/sdf/intermediary_files'
                '/bt_15d_' + datetimes[date_i].strftime('%Y_%H_%M') + '.hdf')
            BT_15_days = f.root.data[lower_ind:upper_ind]
            bt_data = f.root.data[current_ind]
            f.close()
        else:
            f = tables.open_file(
                '/soge-home/projects/seviri_dust/sdf/intermediary_files'
                '/bt_15d_' + str(datetimes[date_i].year) +
                '_' + str(datetimes[date_i].hour) +
                '_' + str(
                    datetimes[date_i].minute) + '.hdf')
            BT_15_days = f.root.data[lower_ind:upper_ind]
            bt_data = f.root.data[current_ind]
            f.close()

        if daily_bt:
            bt087 = bt_data[0]
            bt108 = bt_data[1]
            bt12 = bt_data[2]

        bt_15day_087 = BT_15_days[:, 0]
        bt_15day_108 = BT_15_days[:, 1]
        bt_15day_120 = BT_15_days[:, 2]

        bt_15day_087_mean = np.nanmean(bt_15day_087, axis=0)
        bt_15day_108_mean = np.nanmean(bt_15day_108, axis=0)
        bt_15day_120_mean = np.nanmean(bt_15day_120, axis=0)

        if mesh:
            cloud_lons, cloud_lats = np.meshgrid(cloud_lons, cloud_lats)

        clouds_now_regridded = pinkdust.regrid_data(cloud_lons, cloud_lats,
                                                    orig_lons, orig_lats,
                                                    clouds_now, mesh=True)


        if daily_clouds:
            clouds_now_regridded = clouds_now_regridded > 1

        if show_all_indicators:

            bt_108_087 = bt108-bt087
            bt_120_108 = bt12-bt108
            bt_120_087 = bt12-bt087

            bt_108_087_mean = bt_15day_108_mean-bt_15day_087_mean
            bt_120_108_mean = bt_15day_120_mean-bt_15day_108_mean
            bt_120_087_mean = bt_15day_120_mean - bt_15day_087_mean

            bt_108_087_anom = bt_108_087-bt_108_087_mean
            bt_120_108_anom = bt_120_108-bt_120_108_mean
            bt_120_087_anom = bt_120_087-bt_120_087_mean

            bt_087_anom = bt087 - bt_15day_087_mean
            bt_108_anom = bt108 - bt_15day_108_mean
            bt_120_anom = bt12 - bt_15day_120_mean

            lat_grad, lon_grad = np.gradient(bt_108_087_anom)
            total_grad = np.sqrt(lat_grad ** 2 + lon_grad ** 2)
            grad_108_087_anom = deepcopy(total_grad)

            lat_grad, lon_grad = np.gradient(bt_120_108_anom)
            total_grad = np.sqrt(lat_grad ** 2 + lon_grad ** 2)
            grad_120_108_anom = deepcopy(total_grad)

            lat_grad, lon_grad = np.gradient(bt_120_087_anom)
            total_grad = np.sqrt(lat_grad ** 2 + lon_grad ** 2)
            grad_120_087_anom = deepcopy(total_grad)

            lat_grad, lon_grad = np.gradient(bt_087_anom)
            total_grad = np.sqrt(lat_grad ** 2 + lon_grad ** 2)
            grad_087_anom = deepcopy(total_grad)

            lat_grad, lon_grad = np.gradient(bt_108_anom)
            total_grad = np.sqrt(lat_grad ** 2 + lon_grad ** 2)
            grad_108_anom = deepcopy(total_grad)

            lat_grad, lon_grad = np.gradient(bt_120_anom)
            total_grad = np.sqrt(lat_grad ** 2 + lon_grad ** 2)
            grad_120_anom = deepcopy(total_grad)

            bt_108_087_anom_m = deepcopy(bt_108_087_anom)
            bt_120_087_anom_m = deepcopy(bt_120_087_anom)
            bt_120_108_anom_m = deepcopy(bt_120_108_anom)

            clouds_now_regridded[sdf_now == 1] = 0

            bt_108_087_anom_m[clouds_now_regridded == 1] = np.nan
            bt_120_087_anom_m[clouds_now_regridded == 1] = np.nan
            bt_120_108_anom_m[clouds_now_regridded == 1] = np.nan

        return bt087, bt108, bt12, bt_108_087_anom, bt_120_087_anom, \
               bt_120_108_anom, grad_087_anom, grad_108_anom, grad_120_anom,\
               grad_108_087_anom, grad_120_108_anom, grad_120_087_anom, \
               bt_108_087_anom_m, bt_120_087_anom_m, bt_120_108_anom_m, \
               orig_lons, orig_lats, clouds_now_regridded, sdf_now

if __name__ == '__main__':

    daily_bt = False
    daily_cloudmask = False
    bt_root = ''
    daily_cloudmask_root = '/soge-home/projects/seviri_dust/raw_seviri_data' \
                           '/cloudmask_nc/'
    cloud_mask = False
    debug = False
    check_duration = False
    show_all_indicators = True
    bulk = True

    if daily_cloudmask:
        double_digits = True  # Set to True if using daily cloudmask
        mesh = False # Set this to False if using daily cloudmask
    else:
        double_digits = False
        mesh = True

    year_lower = 2010
    year_upper = 2010
    month_lower = 7
    month_upper = 7
    day_lower = 13
    day_upper = 13
    hour_lower = 1
    hour_upper = 2
    minute_lower = 15
    minute_upper = 0

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

    target_area = utils.load_area(
        '/soge-home/projects/seviri_dust/areas.def',
        'NorthAfrica')
    cloud_lons, cloud_lats = target_area.get_lonlats()

    if bulk == False:

        for date_i in np.arange(0, len(datetimes)):
            print datetimes[date_i]
            bt087_4444, bt108_4444, bt12_4444, bt_108_087_anom_4444, \
            bt_120_087_anom_4444, \
            bt_120_108_anom_4444, grad_087_anom_4444, grad_108_anom_4444, \
            grad_120_anom_4444, \
            grad_108_087_anom_4444, grad_120_108_anom_4444, \
            grad_120_087_anom_4444, bt_108_087_anom_m_4444, bt_120_087_anom_m_4444, \
            bt_120_108_anom_m_4444, \
            orig_lons, orig_lats, clouds_4444, sdf_now_4444 = \
                wrapper(
                datetimes, datestrings, date_i,
                lons, lats,
                cloud_lons, cloud_lats,
                daily_cloudmask, double_digits, mesh, daily_bt)

    indices = np.arange(0, len(datetimes))

    if bulk:
        bt087_1111, bt108_1111, bt12_1111, bt_108_087_anom_1111, \
        bt_120_087_anom_1111, \
        bt_120_108_anom_1111, grad_087_anom_1111, grad_108_anom_1111, \
        grad_120_anom_1111, \
        grad_108_087_anom_1111, grad_120_108_anom_1111, \
        grad_120_087_anom_1111, bt_108_087_anom_m_1111, bt_120_087_anom_m_1111, \
        bt_120_108_anom_m_1111, \
        orig_lons, orig_lats, clouds_1111, sdf_now_1111 = \
            wrapper(
                datetimes, datestrings, indices[0],
                lons, lats,
                cloud_lons, cloud_lats,
                daily_cloudmask, double_digits, mesh, daily_bt)

        bt087_2222, bt108_2222, bt12_2222, bt_108_087_anom_2222, \
        bt_120_087_anom_2222, \
        bt_120_108_anom_2222, grad_087_anom_2222, grad_108_anom_2222, \
        grad_120_anom_2222, \
        grad_108_087_anom_2222, grad_120_108_anom_2222, \
        grad_120_087_anom_2222, bt_108_087_anom_m_2222, bt_120_087_anom_m_2222, \
        bt_120_108_anom_m_2222, \
        orig_lons, orig_lats, clouds_2222, sdf_now_2222 = \
            wrapper(
                datetimes, datestrings, indices[1],
                lons, lats,
                cloud_lons, cloud_lats,
                daily_cloudmask, double_digits, mesh, daily_bt)

        bt087_3333, bt108_3333, bt12_3333, bt_108_087_anom_3333, \
        bt_120_087_anom_3333, \
        bt_120_108_anom_3333, grad_087_anom_3333, grad_108_anom_3333, \
        grad_120_anom_3333, \
        grad_108_087_anom_3333, grad_120_108_anom_3333, \
        grad_120_087_anom_3333, bt_108_087_anom_m_3333, bt_120_087_anom_m_3333, \
        bt_120_108_anom_m_3333, \
        orig_lons, orig_lats, clouds_3333, sdf_now_3333 = \
            wrapper(
                datetimes, datestrings, indices[2],
                lons, lats,
                cloud_lons, cloud_lats,
                daily_cloudmask, double_digits, mesh, daily_bt)

        bt087_4444, bt108_4444, bt12_4444, bt_108_087_anom_4444, \
        bt_120_087_anom_4444, \
        bt_120_108_anom_4444, grad_087_anom_4444, grad_108_anom_4444, \
        grad_120_anom_4444, \
        grad_108_087_anom_4444, grad_120_108_anom_4444, \
        grad_120_087_anom_4444, bt_108_087_anom_m_4444, bt_120_087_anom_m_4444, \
        bt_120_108_anom_m_4444, \
        orig_lons, orig_lats, clouds_4444, sdf_now_4444 = \
            wrapper(
                datetimes, datestrings, indices[3],
                lons, lats,
                cloud_lons, cloud_lats,
                daily_cloudmask, double_digits, mesh, daily_bt)

        arra = (bt_120_087_anom_m_4444 - bt_120_087_anom_m_1111) + (
            bt_120_087_anom_m_4444 - bt_120_087_anom_m_2222) + (
                   bt_120_087_anom_m_4444 - bt_120_087_anom_m_3333)

        arrb = (bt_120_108_anom_m_4444 - bt_120_108_anom_m_1111) + (
        bt_120_108_anom_m_4444 - bt_120_108_anom_m_2222) + (
        bt_120_108_anom_m_4444 - bt_120_108_anom_m_3333)

        arrc = (bt_108_087_anom_m_4444 - bt_108_087_anom_m_1111) + (
        bt_108_087_anom_m_4444 - bt_108_087_anom_m_2222) + (
        bt_108_087_anom_m_4444 - bt_108_087_anom_m_3333)

        arr = (arrc-arrb) + (arra-arrb)
        #arr = (arra-arrb)
        #arr = (arrc-arrb)

        arrm = np.ma.masked_where(np.isnan(arr), arr)
        arra = np.ma.masked_where(np.isnan(arra), arra)
        arrb = np.ma.masked_where(np.isnan(arrb), arrb)
        arrc = np.ma.masked_where(np.isnan(arrc), arrc)

        date_i = indices[-1]

        cpo, cpo_p1 = detect_cpo(arrm)
        plot(cpo, title='CPO_boolean_v3')

        cposk = skeletonize(cpo)

        plot(cposk, title='CPO_boolean_skeleton')


