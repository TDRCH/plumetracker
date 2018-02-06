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
from copy import deepcopy

import utilities
import pinkdust

def sub_wrapper(day_datetimes, j):
    print day_datetimes[j]
    date = day_datetimes[j]
    try:
        date_1 = date + datetime.timedelta(minutes=13)
        global_data, time_slot = pinkdust.load_channels(date,
                                                        date_1)
        data_087 = global_data[8.7].data
        data_108 = global_data[10.8].data
        data_120 = global_data[12.0].data
        data = np.zeros(
            (data_087.shape[0], data_087.shape[1], 3))

        data[:, :, 0] = data_087
        data[:, :, 1] = data_108
        data[:, :, 2] = data_120

    except:
        print 'Adding date to list of missing dates'
        with open('/soge-home/projects/seviri_dust'
                  '/raw_seviri_data/bt_native/missing_msgnative_dates'
                  '.txt', 'a') as my_file:
            my_file.write('\n')
            my_file.write(date.strftime('%Y%m%d%H%M%S'))
        data = np.zeros(
            (msg_area.shape[0], msg_area.shape[1], 3))
        data[:] = np.nan

    if (data.shape[0] != msg_area.shape[0]) or (
                data.shape[1] !=
                msg_area.shape[1]):
        print 'Native file data has wrong dimensions - using the second half of the array'
        data_copy = deepcopy(data)
        data = np.zeros(
            (msg_area.shape[0], msg_area.shape[1], 3))
        data[:] = data_copy[msg_area.shape[0]:, :]

    msg_con_nn = image.ImageContainerNearest(data, msg_area,
                                             radius_of_influence=50000)
    area_con_nn = msg_con_nn.resample(target_area)
    result_data_nn = area_con_nn.image_data

    bt_data = np.zeros((3, lons.shape[0], lons.shape[1]))
    bt_data[0] = result_data_nn[:, :, 0]
    bt_data[1] = result_data_nn[:, :, 1]
    bt_data[2] = result_data_nn[:, :, 2]

    cloudmask = None

    # Now, instead of writing to day array, you write to hdf
    f = tables.open_file('/soge-home/projects/seviri_dust/raw_seviri_data'
                         '/intermediary_files/BT_087_' + date.strftime(
        '%Y%m%d%H%M%S.hdf'),
                         'w')
    atom = tables.Atom.from_dtype(bt_data[0].dtype)
    filters = tables.Filters(complib='blosc', complevel=5)
    ds = f.create_carray(f.root, 'data', atom,
                         bt_data[0].shape,
                         filters=filters)
    ds[:] = bt_data[0]
    f.close()

    f = tables.open_file('/soge-home/projects/seviri_dust/raw_seviri_data'
                         '/intermediary_files/BT_108_' + date.strftime(
        '%Y%m%d%H%M%S.hdf'), 'w')
    atom = tables.Atom.from_dtype(bt_data[1].dtype)
    filters = tables.Filters(complib='blosc', complevel=5)
    ds = f.create_carray(f.root, 'data', atom,
                         bt_data[0].shape,
                         filters=filters)
    ds[:] = bt_data[1]
    f.close()

    f = tables.open_file('/soge-home/projects/seviri_dust/raw_seviri_data'
                         '/intermediary_files/BT_120_' + date.strftime(
        '%Y%m%d%H%M%S.hdf'), 'w')
    atom = tables.Atom.from_dtype(bt_data[2].dtype)
    filters = tables.Filters(complib='blosc', complevel=5)
    ds = f.create_carray(f.root, 'data', atom,
                         bt_data[0].shape,
                         filters=filters)
    ds[:] = bt_data[2]
    f.close()

    print 'Wrote', day_datetimes[j]

def wrapper(y):

    Year_lower = years_list[y]
    Year_upper = years_list[y]
    Month_lower = 9
    Month_upper = 9
    Day_lower = 1
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

                ncfile = pinkdust.create_time_nc_file(
                    '/soge-home/projects/seviri_dust/raw_seviri_data/bt_nc/'
                    + day_datetimes[0].strftime(
                        "%B%Y") + '/BT_' + day_datetimes[
                        0].strftime(
                        '%Y%m%d') + '.nc', day_datetimes, lats, lons)

                day_array[:] = 0

                pool = multiprocessing.Pool()
                for j in np.arange(0, len(day_datetimes)):
                    pool.apply_async(sub_wrapper, args=(day_datetimes, j,))
                pool.close()
                pool.join()

                # Now read back in all the intermediaries and write to nc

                for j in np.arange(0, len(day_datetimes)):

                    f = tables.open_file(
                        '/soge-home/projects/seviri_dust/raw_seviri_data'
                        '/intermediary_files/BT_087_' + day_datetimes[
                            j].strftime(
                            '%Y%m%d%H%M%S.hdf'))
                    arrobj = f.get_node('/data')
                    bt_087 = arrobj.read()
                    f.close()

                    f = tables.open_file(
                        '/soge-home/projects/seviri_dust/raw_seviri_data'
                        '/intermediary_files/BT_108_' + day_datetimes[
                            j].strftime(
                            '%Y%m%d%H%M%S.hdf'))
                    arrobj = f.get_node('/data')
                    bt_108 = arrobj.read()
                    f.close()

                    f = tables.open_file(
                        '/soge-home/projects/seviri_dust/raw_seviri_data'
                        '/intermediary_files/BT_120_' + day_datetimes[
                            j].strftime(
                            '%Y%m%d%H%M%S.hdf'))
                    arrobj = f.get_node('/data')
                    bt_120 = arrobj.read()
                    f.close()

                    ncfile = pinkdust.save_to_existing_nc(ncfile,
                                                          bt_087,
                                                bt_108,
                                            bt_120,
                                            None, day_datetimes[j])

                print 'Day ' + str(days[i]) + ' done'
                ncfile.close()
            print 'Month ' + str(months[k]) + ' done'


if __name__ == '__main__':

    years_list = [2016]
    days_list = [1]
    target_area = utils.load_area(
        '/soge-home/projects/seviri_dust/areas.def',
        'NorthAfrica')
    lons, lats = target_area.get_lonlats()

    msg_area = utils.load_area('/soge-home/projects/seviri_dust/areas.def',
                               'met09globeFull')

    day_array = np.zeros((96, target_area.shape[0], target_area.shape[1], 3))

    # Make a lonlats file
    lonlatsnc = pinkdust. \
        create_lonlats_file('/soge-home/projects/seviri_dust/'
                            'raw_seviri_data/cloudmask_nc/cloudmask_lonlats.nc',
                            lats, lons)
    lonlatsnc.close()

    # Get cloud lats and lons
    testclouds = pygrib.open(glob.glob(
                            '/soge-home/projects/seviri_dust/raw_seviri_data/cloudmask_grib/MSG*-'
                            'SEVI-MSGCLMK-0100-0100-20130601000000' + '*')[0])
    grb = testclouds.select()[0]
    cloudmaskarray = grb.values[:, ::-1]
    cloudlats, cloudlons = grb.latlons()

    cloudmaskarray[cloudmaskarray >= 3] = np.nan

    cloudlats[cloudlats > 90] = np.nan
    cloudlons[cloudlons > 180] = np.nan

    wrapper(0)

    print 1/0

    pool = multiprocessing.Pool()
    for y in np.arange(0, len(years_list)):
        pool.apply_async(wrapper, args=(y,))
        # Pull a netCDF dataset object out so the projection coordinates can
        # be obtained
    pool.close()
    pool.join()





