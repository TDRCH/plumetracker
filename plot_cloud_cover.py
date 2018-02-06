import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import numpy as np
from netCDF4 import Dataset
from mpop.satout import netcdf4
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.cbook as cbook
import matplotlib.cm as cm
import os
import pymc3 as pm
from scipy import ndimage as ndi
import cmocean
from matplotlib import colors
import shelve
from netCDF4 import num2date
import datetime as dt
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
import pandas as pd
import seaborn
from copy import deepcopy
from pyresample import utils
import multiprocessing

import pinkdust
import utilities

def wrapper(i):

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

    datetimes = utilities.get_datetime_objects(time_params)

    for date_i in np.arange(0, len(datetimes)):
        date = datetimes[date_i]
        print date
        year_bool = np.asarray(years == date.year)
        root = np.asarray(year_cm_roots)[year_bool][0]
        if root == root2:
            # We're using TCH's daily cloudmask files. These are values from
            #  0-2 (I think), and need selecting and regridding.
            try:
                clouddata = Dataset(root+date.strftime("%B%Y")+'/cloudmask_' + datetimes[
                        date_i].strftime("%Y%m%d") + '.nc')
                cloudmask_times = num2date(clouddata.variables['time'][:],
                                           clouddata.variables['time'].units)
                cloudmask_times = np.asarray([dt.datetime(j.year, j.month,
                                                          j.day, j.hour,
                                                          j.minute) for j
                                              in cloudmask_times])
                cloudmask_bool = cloudmask_times == date
                clouds_data = clouddata.variables['cloud_mask'][cloudmask_bool]
                clouds_now = clouds_data#[cloudmask_bool]

                clouds_now_regridded = pinkdust.regrid_data(root2_lons, root2_lats,
                                                            target_lons,
                                                            target_lats,
                                                            clouds_now, mesh=True)
                data_array[clouds_now_regridded > 1] += 1
            except:
                print 'No cloud data for', date

        else:
            # We're using Ian's original cloud mask files. These are just a
            # binary one or zero and won't be regridded.
            try:
                clouddata = Dataset(root+
                    datetimes[date_i].strftime("%B").upper(
                    ) + str(datetimes[date_i].year) + '_CLOUDS/eumetsat.cloud.'
                    + datetimes[
                        date_i].strftime("%Y%m%d%H%M") + '.nc')
                clouds_now = clouddata.variables['cmask'][:][0]
                data_array[clouds_now == 1] += 1
            except:
                print 'No cloud data for', date

        if date_i == 100:
            print data_array
            np.save('cloudcover_array', data_array)

            extent = (
            np.min(target_lons), np.max(target_lons), np.min(target_lats),
            np.max(target_lats))
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
            max = 70000

            levels = MaxNLocator(nbins=15).tick_values(min, max)

            # discrete_cmap = utilities.cmap_discretize(cm.RdYlBu_r, 10)
            cmap = cm.get_cmap('RdYlBu_r')
            norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)

            data_array = np.ma.masked_where(np.isnan(data_array),
                                            data_array)

            m.pcolormesh(target_lons, target_lats, data_array,
                         cmap=cmap, norm=norm, vmin=min, vmax=max)

            cbar = plt.colorbar(orientation='horizontal', fraction=0.056,
                                pad=0.06)
            cbar.ax.set_xlabel('Counts of cloud occurrence')
            plt.tight_layout()
            plt.savefig('CPO_multiyear_cloud_frequency_2004_2012.png',
                        bbox_inches='tight')

            plt.close()


if __name__ == '__main__':

    pre_calculated = False

    """
    Goes through every cloud mask file and plots the count of cloud coverage in each region
    Note that daily cloudmask files are not binary -> values greater than 1
    need to be selected, while for Ian's cloudmask files it's binary.
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

    target_area = utils.load_area(
        '/soge-home/projects/seviri_dust/areas.def',
        'NorthAfrica')
    root2_lons, root2_lats = target_area.get_lonlats()

    if pre_calculated:
        data_array = np.load('cloudcover_array.npy')
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
        max = 70000

        levels = MaxNLocator(nbins=15).tick_values(min, max)

        # discrete_cmap = utilities.cmap_discretize(cm.RdYlBu_r, 10)
        cmap = cm.get_cmap('RdYlBu_r')
        norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)

        # m.pcolormesh(lons, lats, correlation_array, extent=extent,
        #             origin='lower',
        #         interpolation='none',
        #         cmap=cmap, vmin=min, vmax=max+1)

        data_array = np.ma.masked_where(np.isnan(data_array),
                                        data_array)

        m.pcolormesh(target_lons, target_lats, data_array,
                     cmap=cmap, norm=norm, vmin=min, vmax=max)

        cbar = plt.colorbar(orientation='horizontal', fraction=0.056, pad=0.06)
        cbar.ax.set_xlabel('Counts of cloud occurrence')
        plt.tight_layout()
        plt.savefig('CPO_multiyear_cloud_frequency_2004_2012.png',
                    bbox_inches='tight')

        plt.close()

    else:

        years = np.arange(2004, 2013)
        #years = np.asarray([2010])
        #years = np.arange(2004, 2011)

        root1 = '/soge-home/data/satellite/meteosat/seviri/15-min/0.03x0.03/cloudmask/nc/'
        root2 = '/soge-home/projects/seviri_dust/raw_seviri_data/cloudmask_nc/'

        year_cm_roots = [root1, root1, root1, root1, root1, root1, root1, root2,
                         root2]

        year_cm_roots = [root1, root1, root1, root1, root1, root1, root1]
        #year_cm_roots = [root2, root2]
        #year_cm_roots = [root1]
        data_array = np.zeros(target_lons.shape)

        pool = multiprocessing.Pool()
        for i in years:
            pool.apply_async(wrapper, args=(i,))
        pool.close()
        pool.join()

