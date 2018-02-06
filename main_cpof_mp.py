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
from numba import jit

import utilities
import plumes
import pinkdust
import plotting
import get_llj_prob_model

def wrapper(bt_120_108_anom_m_prev_1, bt_120_108_anom_m_prev_2,
            bt_120_108_anom_m_prev_3, bt_108_087_anom_m_prev_1,
            bt_108_087_anom_m_prev_2, bt_108_087_anom_m_prev_3,
            bt_120_087_anom_m_prev_1, bt_120_087_anom_m_prev_2,
            bt_120_087_anom_m_prev_3,
            datetimes,
            datestrings, date_i, lons,
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

        btdiff_2_15daymean = bt_15day_120_mean - bt_15day_087_mean

        bt_108_087 = bt108 - bt087
        bt_120_108 = bt12 - bt108
        bt_120_087 = bt12 - bt087

        bt_108_087_mean = bt_15day_108_mean - bt_15day_087_mean
        bt_120_108_mean = bt_15day_120_mean - bt_15day_108_mean
        bt_120_087_mean = bt_15day_120_mean - bt_15day_087_mean

        bt_108_087_anom = bt_108_087 - bt_108_087_mean
        bt_120_108_anom = bt_120_108 - bt_120_108_mean
        bt_120_087_anom = bt_120_087 - bt_120_087_mean


        if mesh:
            cloud_lons, cloud_lats = np.meshgrid(cloud_lons, cloud_lats)

        clouds_now_regridded = pinkdust.regrid_data(cloud_lons, cloud_lats,
                                                    orig_lons, orig_lats,
                                                    clouds_now, mesh=True)

        bt_108_087_anom_m = deepcopy(bt_108_087_anom)
        bt_120_087_anom_m = deepcopy(bt_120_087_anom)
        bt_120_108_anom_m = deepcopy(bt_120_108_anom)

        clouds_now_regridded[sdf_now == 1] = 0

        bt_108_087_anom_m[clouds_now_regridded == 1] = np.nan
        bt_120_087_anom_m[clouds_now_regridded == 1] = np.nan
        bt_120_108_anom_m[clouds_now_regridded == 1] = np.nan

        #btdiff_2 = bt12 - bt087
        #btdiff_2_anom = btdiff_2 - btdiff_2_15daymean

        if bt_108_087_anom_m_prev_1 != None:

            arra = (bt_120_087_anom_m - bt_120_087_anom_m_prev_3) + (
                bt_120_087_anom_m - bt_120_087_anom_m_prev_2) + (
                       bt_120_087_anom_m - bt_120_087_anom_m_prev_1)

            arrb = (bt_120_108_anom_m - bt_120_108_anom_m_prev_3) + (
                bt_120_108_anom_m- bt_120_108_anom_m_prev_2) + (
                       bt_120_108_anom_m - bt_120_108_anom_m_prev_1)

            arrc = (bt_108_087_anom_m - bt_108_087_anom_m_prev_3) + (
                bt_108_087_anom_m - bt_108_087_anom_m_prev_2) + (
                       bt_108_087_anom_m - bt_108_087_anom_m_prev_1)

            detected_bt = (arrc - arrb) + (arra - arrb)

        else:
            detected_bt = np.zeros((bt_108_087_anom.shape))

        if date_i == 0:
            bt_120_087_anom_m_prev_1 = deepcopy(bt_120_087_anom_m)
            bt_120_108_anom_m_prev_1 = deepcopy(bt_120_108_anom_m)
            bt_108_087_anom_m_prev_1 = deepcopy(bt_108_087_anom_m)
        elif date_i == 1:
            bt_120_087_anom_m_prev_2 = deepcopy(bt_120_087_anom_m_prev_1)
            bt_120_108_anom_m_prev_2 = deepcopy(bt_120_108_anom_m_prev_1)
            bt_108_087_anom_m_prev_2 = deepcopy(bt_108_087_anom_m_prev_1)
            bt_120_087_anom_m_prev_1 = deepcopy(bt_120_087_anom_m)
            bt_120_108_anom_m_prev_1 = deepcopy(bt_120_108_anom_m)
            bt_108_087_anom_m_prev_1 = deepcopy(bt_108_087_anom_m)
        elif date_i == 2:
            bt_120_087_anom_m_prev_3 = deepcopy(bt_120_087_anom_m_prev_2)
            bt_120_108_anom_m_prev_3 = deepcopy(bt_120_108_anom_m_prev_2)
            bt_108_087_anom_m_prev_3 = deepcopy(bt_108_087_anom_m_prev_2)
            bt_120_087_anom_m_prev_2 = deepcopy(bt_120_087_anom_m_prev_1)
            bt_120_108_anom_m_prev_2 = deepcopy(bt_120_108_anom_m_prev_1)
            bt_108_087_anom_m_prev_2 = deepcopy(bt_108_087_anom_m_prev_1)
            bt_120_087_anom_m_prev_1 = deepcopy(bt_120_087_anom_m)
            bt_120_108_anom_m_prev_1 = deepcopy(bt_120_108_anom_m)
            bt_108_087_anom_m_prev_1 = deepcopy(bt_108_087_anom_m)
        elif date_i > 2:
            bt_120_087_anom_m_prev_3 = deepcopy(bt_120_087_anom_m_prev_2)
            bt_120_108_anom_m_prev_3 = deepcopy(bt_120_108_anom_m_prev_2)
            bt_108_087_anom_m_prev_3 = deepcopy(bt_108_087_anom_m_prev_2)
            bt_120_087_anom_m_prev_2 = deepcopy(bt_120_087_anom_m_prev_1)
            bt_120_108_anom_m_prev_2 = deepcopy(bt_120_108_anom_m_prev_1)
            bt_108_087_anom_m_prev_2 = deepcopy(bt_108_087_anom_m_prev_1)
            bt_120_087_anom_m_prev_1 = deepcopy(bt_120_087_anom_m)
            bt_120_108_anom_m_prev_1 = deepcopy(bt_120_108_anom_m)
            bt_108_087_anom_m_prev_1 = deepcopy(bt_108_087_anom_m)

        if daily_clouds:
            clouds_now_regridded = clouds_now_regridded > 1

        ### PASS I ###
        # In the FIRST PASS the LORD sayeth unto the image, 'Let all
        # whose BTD is below -10K be classified as CPOs, and remove the
        # tiny ones'
        # And those who fulfilled this condition were classified,
        # and it was good

        # NOTE: Why cloud mask here? What if there was a strong cloud
        # gradient in the previous three timesteps which disappeared in this
        #  one but still pushed us over the threshold? It wouldn't be cloud
        # masked
        detected_bt_um = deepcopy(detected_bt)
        #if cloud_mask:
        #   detected_bt[clouds_now > 0] = np.nan

        levels = np.arange(-24, 40, 8)
        cpo_mask_pass_1 = detected_bt < -6

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
            if np.any(detected_bt[target_region == 1] < -20):
                # Next check if there is a generous gradient contained
                # within this region somewhere
                pass
            else:
                cpo_mask_pass_2[target_region == 1] = 0

        # For identified CPO regions, undo the convolution on the cloud
        # mask

        cpo_mask_um = detected_bt_um < -6
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

        if debug:
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
            # labels = [left,right,top,bottom]
            m.drawparallels(parallels, labels=[False, True, True, False],
                            linewidth=0.5)
            meridians = np.arange(-20., 17., 4.)
            m.drawmeridians(meridians, labels=[True, False, False, True],
                            linewidth=0.5)

            min = -20
            max = 5

            levels = MaxNLocator(nbins=15).tick_values(min, max)

            cmap = cm.get_cmap('Blues_r')
            norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)

            m.pcolormesh(orig_lons, orig_lats, detected_bt,
                         cmap=cmap, vmin=min, vmax=max, norm=norm)

            cbar = plt.colorbar(orientation='horizontal', fraction=0.056,
                                pad=0.06)
            cbar.ax.set_xlabel('BTdiff 2 anom diff')
            plt.tight_layout()
            plt.savefig('BTdiff_2_anom_diff' + datestrings[
                date_i] + '.png',
                        bbox_inches='tight')

            plt.close()

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
            # labels = [left,right,top,bottom]
            m.drawparallels(parallels, labels=[False, True, True, False],
                            linewidth=0.5)
            meridians = np.arange(-20., 17., 4.)
            m.drawmeridians(meridians, labels=[True, False, False, True],
                            linewidth=0.5)

            min = 180
            max = 320

            levels = MaxNLocator(nbins=15).tick_values(min, max)

            cmap = cm.get_cmap('Blues_r')
            norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)

            m.pcolormesh(orig_lons, orig_lats, bt108,
                         cmap=cmap, vmin=min, vmax=max, norm=norm)

            cbar = plt.colorbar(orientation='horizontal', fraction=0.056,
                                pad=0.06)
            cbar.ax.set_xlabel('BT 10.8')
            plt.tight_layout()
            plt.savefig('BT108_' + datestrings[
                date_i] + '.png',
                        bbox_inches='tight')

            plt.close()

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

            indicators = [bt087, bt108, bt12,
                          bt_108_087_anom, bt_120_108_anom,
                          bt_120_087_anom, grad_108_087_anom,
                          grad_120_108_anom, grad_120_087_anom,
                          grad_087_anom, grad_108_anom, grad_120_anom]

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

            mins = [180, 180, 180, -13, -13, -13, 0, 0, 0, 0, 0, 0]
            maxs = [320, 320, 320, 6, 6, 6, 9, 9, 9, 50, 50, 50]
            labels = ['bt087', 'bt108', 'bt120', 'bt_108_087_anom',
                      'bt_120_108_anom', 'bt_120_087_anom','grad_108_087_anom',
                      'grad_120_108_anom', 'grad_120_087_anom',
                      'grad_087_anom', 'grad_108_anom', 'grad_120_anom']

            for i in np.arange(0, len(indicators)):

                min = mins[i]
                max = maxs[i]

                levels = MaxNLocator(nbins=15).tick_values(min, max)

                cmap = cm.get_cmap('Blues_r')
                norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)

                plot = m.pcolormesh(orig_lons, orig_lats, indicators[i],
                             cmap=cmap, vmin=min, vmax=max, norm=norm)

                cbar = plt.colorbar(orientation='horizontal',
                                    fraction=0.056,
                                    pad=0.06)

                cbar.ax.set_xlabel(labels[i])
                plt.tight_layout()
                plt.savefig(labels[i] + '_' + datestrings[
                    date_i] + '.png',
                            bbox_inches='tight')

                cbar.remove()
                plot.remove()

        return cpo_mask_pass_2, bt_120_108_anom_m, \
               bt_120_108_anom_m_prev_1, bt_120_108_anom_m_prev_2, \
               bt_120_108_anom_m_prev_3, bt_108_087_anom_m,\
               bt_108_087_anom_m_prev_1, bt_108_087_anom_m_prev_2, \
               bt_108_087_anom_m_prev_3, bt_120_087_anom_m,\
               bt_120_087_anom_m_prev_1, bt_120_087_anom_m_prev_2, \
               bt_120_087_anom_m_prev_3
    else:
        if mesh:
            empty_arr = np.zeros((lats.shape[0], lons.shape[1]))
        else:
            empty_arr  = np.zeros((lats.shape[0], lons.shape[0]))
        empty_arr[:] = np.nan
        return empty_arr, None, None, None

if __name__ == '__main__':

    daily_bt = False
    daily_cloudmask = False
    bt_root = ''
    daily_cloudmask_root = '/soge-home/projects/seviri_dust/raw_seviri_data' \
                           '/cloudmask_nc/'
    cloud_mask = False
    debug = True
    check_duration = False
    show_all_indicators = True

    if daily_cloudmask:
        double_digits = True  # Set to True if using daily cloudmask
        mesh = False # Set this to False if using daily cloudmask
    else:
        double_digits = False
        mesh = True

    year_lower = 2010
    year_upper = 2010
    month_lower = 7
    month_upper = 8
    day_lower = 2
    day_upper = 31
    hour_lower = 15
    hour_upper = 23
    minute_lower = 0
    minute_upper = 45

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

    bt_120_108_anom_m_prev_1 = None
    bt_120_108_anom_m_prev_2 = None
    bt_120_108_anom_m_prev_3 = None
    bt_108_087_anom_m_prev_1 = None
    bt_108_087_anom_m_prev_2 = None
    bt_108_087_anom_m_prev_3 = None
    bt_120_087_anom_m_prev_1 = None
    bt_120_087_anom_m_prev_2 = None
    bt_120_087_anom_m_prev_3 = None

    cpo_prev_1 = None
    cpo_prev_2 = None
    cpo_prev_3 = None
    cpo_prev_4 = None

    for date_i in np.arange(0, len(datetimes)):
        print datetimes[date_i]
        cpo_now, bt_120_108_anom_m, \
               bt_120_108_anom_m_prev_1, bt_120_108_anom_m_prev_2, \
               bt_120_108_anom_m_prev_3, bt_108_087_anom_m,\
               bt_108_087_anom_m_prev_1, bt_108_087_anom_m_prev_2, \
               bt_108_087_anom_m_prev_3, bt_120_087_anom_m,\
               bt_120_087_anom_m_prev_1, bt_120_087_anom_m_prev_2, \
               bt_120_087_anom_m_prev_3 = wrapper(
            bt_120_108_anom_m_prev_1, bt_120_108_anom_m_prev_2,
            bt_120_108_anom_m_prev_3, bt_108_087_anom_m_prev_1,
            bt_108_087_anom_m_prev_2, bt_108_087_anom_m_prev_3,
            bt_120_087_anom_m_prev_1, bt_120_087_anom_m_prev_2,
            bt_120_087_anom_m_prev_3, datetimes, datestrings, date_i,
            lons, lats,
            cloud_lons, cloud_lats,
            daily_cloudmask, double_digits, mesh, daily_bt)
        cpo_now[cpo_now > 0] = 1

        if np.any(np.isnan(cpo_now)):
            cpo_now[:] = False
            cpo_now_raw = deepcopy(cpo_now)
            cpo_now = cpo_now == 1
            cpo_now_raw = cpo_now_raw == 1

        else:
            cpo_now_raw = deepcopy(cpo_now)
            if check_duration:
                if cpo_prev_4 != None:
                    # Get overlap of 3 and 4
                    overlap_4 = cpo_prev_4 & cpo_prev_3
                    #print np.unique(overlap_4)

                    # Label all blobs in 3
                    cpo_prev_3_lab, num = measurements.label(cpo_prev_3)

                    # Go through each blob in 3 and check if an overlap True exists. If
                    blob_ids = np.unique(cpo_prev_3_lab)
                    blob_ids = blob_ids[blob_ids != 0]

                    for i in np.arange(0, len(blob_ids)):
                        target_region = cpo_prev_3_lab == blob_ids[i]
                        if True in overlap_4[target_region]:
                            continue
                        else:
                            cpo_prev_3_lab[target_region] = 0

                    cpo_prev_3_lab[cpo_prev_3_lab > 0] = 1

                    # Get overlap of 2 and 3
                    overlap_3 = cpo_prev_3 & cpo_prev_2
                    #print np.unique(overlap_3)

                    # Label all blobs in 2
                    cpo_prev_2_lab, num = measurements.label(cpo_prev_2)

                    # Go through each blob in 3 and check if an overlap True exists. If
                    blob_ids = np.unique(cpo_prev_2_lab)
                    blob_ids = blob_ids[blob_ids != 0]

                    for i in np.arange(0, len(blob_ids)):
                        target_region = cpo_prev_2_lab == blob_ids[i]
                        if True in overlap_3[target_region]:
                            continue
                        else:
                            cpo_prev_2_lab[target_region] = 0

                    cpo_prev_2_lab[cpo_prev_2_lab > 0] = 1

                    # Get overlap of 1 and 2
                    overlap_2 = cpo_prev_2_lab & cpo_prev_1
                    #print np.unique(overlap_2)

                    # Label all blobs in 1
                    cpo_prev_1_lab, num = measurements.label(cpo_prev_1)

                    # Go through each blob in 1 and check if an overlap True exists. If
                    blob_ids = np.unique(cpo_prev_1_lab)
                    blob_ids = blob_ids[blob_ids != 0]

                    for i in np.arange(0, len(blob_ids)):
                        target_region = cpo_prev_1_lab == blob_ids[i]
                        if True in overlap_2[target_region]:
                            continue
                        else:
                            cpo_prev_1_lab[target_region] = 0

                    cpo_prev_1_lab[cpo_prev_1_lab > 0] = 1

                    # Get overlap of now and 1
                    overlap_1 = cpo_prev_1 & cpo_now
                    #print np.unique(overlap_1)

                    # Label all blobs in now
                    cpo_now_lab, num = measurements.label(cpo_now)

                    # Go through each blob in now and check if an overlap True exists. If
                    blob_ids = np.unique(cpo_now_lab)
                    blob_ids = blob_ids[blob_ids != 0]

                    for i in np.arange(0, len(blob_ids)):
                        target_region = cpo_now_lab == blob_ids[i]
                        if True in overlap_1[target_region]:
                            continue
                        else:
                            cpo_now_lab[target_region] = 0

                    cpo_now_lab[cpo_now_lab > 0] = 1

                    cpo_now = deepcopy(cpo_now_lab)

        #print np.unique(cpo_now)

        cpo_prev_4 = deepcopy(cpo_prev_3)
        cpo_prev_3 = deepcopy(cpo_prev_2)
        cpo_prev_2 = deepcopy(cpo_prev_1)
        cpo_prev_1 = deepcopy(cpo_now_raw)

        filename = '/soge-home/projects/seviri_dust/cpof/'+datetimes[
            date_i].strftime("%B%Y")+'/CPOF_'+datetimes[date_i].strftime(
            "%Y%m%d%H%M")+'.nc'

        if debug:
            filename = 'test.nc'

        cpoffile = Dataset(filename, 'w', format='NETCDF4_CLASSIC')
        cpoffile.description = 'CPOF values from EUMETSAT SEVIRI on a ' \
                              'regular lat/lon grid.'
        cpoffile.createDimension('time', None)
        cpoffile.createDimension('lat', lats.shape[0])
        cpoffile.createDimension('lon', lons.shape[1])

        # Variables
        times = cpoffile.createVariable('time', np.float64, ('time'))
        cpofvar = cpoffile.createVariable('CPOF', np.int32, ('lat', 'lon'))

        # Attributes
        times.units = 'hours since 0001-01-01 00:00:00'
        times.calendar = 'gregorian'

        # Writing data to file and closing
        times[:] = date2num(datetimes[date_i], units=times.units,
                            calendar=times.calendar)
        if cpo_now.shape[0] != lats.shape[0] or cpo_now.shape[1] != \
                lons.shape[1]:
            print 'Wrong shape of array for file!'
            cpo_now = np.zeros((lats.shape[0], lons.shape[1]))
            cpo_now[:] = False
            cpo_now = cpo_now == 1

            # Reset past CPOs
            cpo_prev_4 = None
            cpo_prev_3 = None
            cpo_prev_2 = None
            cpo_prev_1 = None
        cpofvar[:] = cpo_now
        cpoffile.close()