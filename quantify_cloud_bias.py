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
from pyresample import utils
import datetime as dt
from copy import deepcopy

import utilities
import plumes
import pinkdust
import plotting
import get_llj_prob_model

def count_cloud_pixel_duration():
    """
    Calculates and maps the mean duration of cloud cover at each pixel
    :return:
    """

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
    target_lons, target_lats = np.meshgrid(cloud_lons, cloud_lats)

    """

    target_area = utils.load_area(
        '/soge-home/projects/seviri_dust/areas.def',
        'NorthAfrica')
    root2_lons, root2_lats = target_area.get_lonlats()

    years = np.arange(2004, 2011)
    years = [2010]

    datetimes = []

    for i in years:
        year_lower = i
        year_upper = i
        month_lower = 6
        month_upper = 8
        day_lower = 1
        day_upper = 31
        hour_lower = 0
        hour_upper = 23
        minute_lower = 0
        minute_upper = 45

        time_params = np.array([year_lower, year_upper, month_lower,
                                month_upper, day_lower, day_upper,
                                hour_lower, hour_upper, minute_lower,
                                minute_upper])

        year_datetimes = utilities.get_datetime_objects(time_params)
        datetimes = np.append(datetimes, year_datetimes)

    root1 = '/soge-home/data/satellite/meteosat/seviri/15-min/0.03x0.03/cloudmask/nc/'
    root2 = '/soge-home/projects/seviri_dust/raw_seviri_data/cloudmask_nc/'

    year_cm_roots = [root1, root1, root1, root1, root1, root1, root1]
    year_cm_roots = [root1]

    duration_array = np.zeros(target_lons.shape)
    total_duration_array = np.zeros(target_lons.shape)
    counts_array = np.zeros(target_lons.shape)

    clouds_previous = None

    for date_i in np.arange(0, len(datetimes)):
        date = datetimes[date_i]
        print date
        year_bool = np.asarray([j == date.year for j in years])
        root = np.asarray(year_cm_roots)[year_bool][0]
        if root == root2:
            # We're using TCH's daily cloudmask files. These are values from
            #  0-2 (I think), and need selecting and regridding.
            try:
                clouddata = Dataset(root+date.strftime("%B%Y")+'/cloudmask_' + datetimes[
                        date_i].strftime("%Y%m%d") + '.nc')
            except:
                print 'No cloud data for', date
                continue
            cloudmask_times = num2date(clouddata.variables['time'][:],
                                       clouddata.variables['time'].units)
            cloudmask_times = np.asarray([dt.datetime(j.year, j.month,
                                                            j.day, j.hour,
                                                            j.minute) for j
                                          in cloudmask_times])
            cloudmask_bool = cloudmask_times == datetimes[date_i]
            clouds_now = clouddata.variables['cloud_mask'][cloudmask_bool]
            #clouds_now = clouds_data[cloudmask_bool]

            clouds_now = clouds_now > 1

            clouds_now_regridded = pinkdust.regrid_data(root2_lons, root2_lats,
                                                        target_lons,
                                                        target_lats,
                                                        clouds_now, mesh=True)

            if clouds_previous == None:
                clouds_previous = deepcopy(clouds_now_regridded)
            else:
                clouds_previous = clouds_previous == 1
                clouds_now = clouds_now_regridded == 1
                cloud_overlaps = clouds_previous & clouds_now
                duration_array[clouds_previous == 1] += 1
                duration_array[cloud_overlaps == 1] += 1
                total_duration_array[cloud_overlaps == 0] += duration_array[
                    cloud_overlaps == 0]
                c1 = cloud_overlaps == 0
                c1 = c1 == 1
                c2 = duration_array > 0
                c2 = c2 == 1
                combo_bool = c1 & c2
                counts_array[combo_bool] += 1
                duration_array[cloud_overlaps == 0] = 0
                clouds_previous = deepcopy(clouds_now_regridded)

            if date.day == 31 and date.hour == 23 and date.minute == 45:
                np.save('cloud_duration_array_'+str(date.year),
                        total_duration_array)
                np.save('cloud_events_array_' + str(date.year),
                        counts_array)

            print total_duration_array
            print counts_array


        else:
            # We're using Ian's original cloud mask files. These are just a
            # binary one or zero and won't be regridded.
            try:
                clouddata = Dataset(root+
                    datetimes[date_i].strftime("%B").upper(
                    ) + str(datetimes[date_i].year) + '_CLOUDS/eumetsat.cloud.'
                    + datetimes[
                        date_i].strftime("%Y%m%d%H%M") + '.nc')
            except:
                print 'No cloud data for', date
                continue
            clouds_now = clouddata.variables['cmask'][:][0]
            clouds_now = np.ma.filled(clouds_now, 0)

            if clouds_previous == None:
                clouds_previous = deepcopy(clouds_now)
            else:
                clouds_previous = clouds_previous == 1
                clouds_now = clouds_now == 1
                cloud_overlaps = clouds_previous & clouds_now
                duration_array[cloud_overlaps == 1] += 1
                total_duration_array[cloud_overlaps == 0] += duration_array[
                    cloud_overlaps == 0]
                c1 = cloud_overlaps == 0
                c1 = c1 == 1
                c2 = duration_array > 0
                c2 = c2 == 1
                combo_bool = c1 & c2
                counts_array[combo_bool] += 1
                duration_array[cloud_overlaps == 0] = 0
                clouds_previous = deepcopy(clouds_now)

            if date.day == 31 and date.hour == 23 and date.minute == 45:
                np.save('cloud_duration_array_'+str(date.year),
                        total_duration_array)
                np.save('cloud_events_array_' + str(date.year),
                        counts_array)

            print total_duration_array
            print counts_array

    np.save('cloud_duration_array', total_duration_array)
    np.save('cloud_events_array', counts_array)

    """

    data_array = np.load('cloud_duration_array.npy')
    counts_array = np.load('cloud_events_array.npy')
    data_array_mean = data_array/counts_array
    extent = (
        np.min(target_lons), np.max(target_lons), np.min(target_lats),
        np.max(target_lats))
    m = Basemap(projection='cyl', llcrnrlon=extent[0], urcrnrlon=extent[1],
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
    max = 6

    levels = MaxNLocator(nbins=15).tick_values(min, max)

    # discrete_cmap = utilities.cmap_discretize(cm.RdYlBu_r, 10)
    cmap = cm.get_cmap('RdYlBu_r')
    norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)

    # m.pcolormesh(lons, lats, correlation_array, extent=extent,
    #             origin='lower',
    #         interpolation='none',
    #         cmap=cmap, vmin=min, vmax=max+1)

    #data_array_masked = np.ma.masked_where(np.isnan(data_array),
    #                                       data_array)

    m.pcolormesh(target_lons, target_lats, (data_array_mean*15)/60,
                 cmap=cmap, norm=norm, vmin=min, vmax=max)

    cbar = plt.colorbar(orientation='horizontal', fraction=0.056, pad=0.06)
    cbar.ax.set_xlabel('Mean cloud duration (hours)')
    plt.tight_layout()
    plt.savefig('Cloud_duration_2004.png',
                bbox_inches='tight')

    plt.close()

    print 'Mean cloud duration:'
    print np.nanmean(data_array_mean), 'timesteps'
    print 'Or', (np.nanmean(data_array_mean)*15)/60, 'hours'

if __name__ == '__main__':
    count_cloud_pixel_duration()