import numpy as np
import multiprocessing
from netCDF4 import Dataset
import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpop.satellites import GeostationaryFactory
from mpop.satout import netcdf4
# from mpop.utils import debug_on

# debug_on()
import datetime
import numpy as np
import iris
import iris.quickplot as qplt
from netCDF4 import Dataset
from netCDF4 import date2num
from netCDF4 import num2date
from pycoast import ContourWriterAGG
from PIL import Image
from pyproj import Proj
from mpl_toolkits.basemap import Basemap
from pyresample.geometry import SwathDefinition
from pyresample.kd_tree import resample_nearest
from pyresample import utils
import pygrib
import glob
# import multiprocessing
# from multiprocessing import Process, Value, Array
import multiprocess
from multiprocess import Process
import tables
from pyresample import image, geometry
import email
import smtplib

import utilities
import pinkdust

def sub_wrapper(day_datetimes, j):
    print day_datetimes[j]
    try:
        clouddata = Dataset(
            '/soge-home/data/satellite/meteosat/seviri/15-min/'
            'native/cloudmask'
            '/nc/'
            +
            day_datetimes[j].strftime("%B").upper(
            ) + str(day_datetimes[
                        j].year) + '_CLOUDS/eumetsat.cloud.'
            + day_datetimes[j].strftime("%Y%m%d%H%M") + '.nc')
        cloudmaskarray = clouddata.variables['cmask'][:]
        cloudmask_regridded = pinkdust.regrid_data(
            cloudlons,
            cloudlats,
            lons,
            lats,
            cloudmaskarray,
            mesh=True)
    except:
        print 'Missing cloud mask data for', \
            day_datetimes[j]
        cloudmask_regridded = np.zeros((lats.shape[0],
                                        lons.shape[1]))
        cloudmask_regridded[:] = np.nan

    # Now, instead of writing to day array, you write to hdf
    f = tables.open_file('/soge-home/projects/seviri_dust/raw_seviri_data'
                         '/intermediary_files/cloudmask_' +
                         day_datetimes[j].strftime(
        '%Y%m%d%H%M%S.hdf'),
                         'w')
    atom = tables.Atom.from_dtype(cloudmask_regridded.dtype)
    filters = tables.Filters(complib='blosc', complevel=5)
    ds = f.create_carray(f.root, 'data', atom,
                         cloudmask_regridded.shape,
                         filters=filters)
    ds[:] = cloudmask_regridded
    f.close()

    print 'Wrote', day_datetimes[j]

def wrapper(y):

    Year_lower = years_list[y]
    Year_upper = years_list[y]
    Month_lower = 5
    Month_upper = 9
    Day_lower = 24
    Day_upper = 8
    Hour_lower = 0
    Hour_upper = 23
    Minute_lower = 0
    Minute_upper = 45

    # Generate datetime objects corresponding to these time bounds
    time_params = np.array(
        [Year_lower, Year_upper, Month_lower, Month_upper,
         Day_lower, Day_upper, Hour_lower, Hour_upper,
         Minute_lower, Minute_upper])

    datetimes = utilities.get_datetime_objects(time_params)

    years = np.unique(np.asarray(([j.year for j in datetimes])))

    months = np.unique(np.asarray(([j.month for j in datetimes])))

    days = np.unique(np.asarray(([j.day for j in datetimes])))
    for m in np.arange(0, len(years)):
        for k in np.arange(0, len(months)):

            for i in np.arange(0, len(days)):
                day_datetimes_bool = np.asarray(
                    [j.day == days[i] and j.month ==
                     months[k] and j.year ==
                     years[m] for j in datetimes])
                day_datetimes = datetimes[day_datetimes_bool]

                if len(day_datetimes) == 0:
                    continue

                ncfile = pinkdust.create_time_cloudmask_nc_file(
                    '/soge-home/projects/seviri_dust/raw_seviri_data'
                    '/cloudmask_nc/'+day_datetimes[0].strftime('%B%Y')
                    +'/cloudmask_'+day_datetimes[0].strftime(
                        '%Y%m%d')+'.nc',
                    day_datetimes, lats, lons)

                #sub_wrapper(day_datetimes, 0)

                pool = multiprocessing.Pool()
                for j in np.arange(0, len(day_datetimes)):
                    pool.apply_async(sub_wrapper, args=(day_datetimes, j,))
                pool.close()
                pool.join()

                for j in np.arange(0, len(day_datetimes)):
                    f = tables.open_file(
                        '/soge-home/projects/seviri_dust/raw_seviri_data'
                        '/intermediary_files/cloudmask_' + day_datetimes[
                            j].strftime(
                            '%Y%m%d%H%M%S.hdf'))
                    arrobj = f.get_node('/data')
                    cloudmask = arrobj.read()
                    f.close()

                    ncfile = pinkdust.save_cloudmask_to_existing_nc(ncfile,
                                                                    cloudmask,
                                                                    day_datetimes[
                                                                        j])

                print 'Day ' + str(days[i]) + ' done'
                ncfile.close()
            print 'Month ' + str(months[k]) + ' done'


if __name__ == '__main__':

    years_list = [2011]
    days_list = [1]
    target_area = utils.load_area(
        '/soge-home/projects/seviri_dust/areas.def',
        'NorthAfrica')
    lons, lats = target_area.get_lonlats()

    # Make a lonlats file
    lonlatsnc = pinkdust. \
        create_lonlats_file('/soge-home/projects/seviri_dust/'
                            'raw_seviri_data/cloudmask_nc/cloudmask_lonlats.nc',
                            lats, lons)
    lonlatsnc.close()

    # Get cloud lats and lons
    # Get cloud lats and lons
    cloud_test = Dataset(
        '/soge-home/data/satellite/meteosat/seviri/15-min/'
        'native/lonlats.NA_MiddleEast.nc')
    cloudlons = cloud_test.variables['longitude'][:]
    cloudlats = cloud_test.variables['latitude'][:]
    lonmask = cloudlons > 360
    latmask = cloudlats > 90
    cloudlons = np.ma.array(cloudlons, mask=lonmask)
    cloudlats = np.ma.array(cloudlats, mask=latmask)

    day_array = np.zeros((96, lons.shape[0], lats.shape[1]))

    wrapper(0)
    print 1/0

    pool = multiprocessing.Pool()
    for y in np.arange(0, len(years_list)):
        pool.apply_async(wrapper, args=(y,))
        # Pull a netCDF dataset object out so the projection coordinates can
        # be obtained
    pool.close()
    pool.join()





