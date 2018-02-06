import numpy as np
import multiprocessing
from netCDF4 import Dataset
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpop.satellites import GeostationaryFactory
from mpop.satout import netcdf4
from mpop.utils import debug_on
debug_on()
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
import multiprocessing
from multiprocessing import Process, Value, Array
import tables
from pyresample import image, geometry

import utilities
import pinkdust

def wrapper(date):
    print date
    test1 = datetime.datetime.now()
    global_data, time_slot = pinkdust.load_channels(date)
    print 'Loading data to a scene object:', datetime.datetime.now() - test1

    # Reproject it to North Africa - note that this requires a custom
    # projection in areas.def config file
    test2 = datetime.datetime.now()
    #projected_data = global_data.project("NorthAfrica")
    msg_area = utils.load_area('/soge-home/projects/seviri_dust/areas.def',
                           'met09globeFull')
    target_area = utils.load_area('/soge-home/projects/seviri_dust/areas.def',
                           'NorthAfrica')

    #print global_data.data
    data = global_data[8.7].data
    msg_con_nn = image.ImageContainerNearest(data, msg_area,
                                             radius_of_influence=50000)
    area_con_nn = msg_con_nn.resample(target_area)
    result_data_nn = area_con_nn.image_data
    lons, lats = target_area.get_lonlats()
    print np.min(lons)
    print np.max(lons)
    print np.min(lats)
    print np.max(lats)
    """

    print 'Constraining the scene to North Africa:', datetime.datetime.now() \
                                                     - test2

    test3 = datetime.datetime.now()
    netcdf4.save(projected_data,
                 '/soge-home/projects/seviri_dust/raw_seviri_data/bt_nc/'
                 + date.strftime(
                     "%B%Y") + '/North_Africa_SEVIRI_BTs_' + str(
                     date.strftime("%Y%m%d%H%M")) + '.nc',
                 compression=True,
                 dtype=np.int16,
                 band_axis=2) #area_aggregation=True
    # #time_dimension=True
    print 'Saving to netCDF:', datetime.datetime.now() - test3
    test4 = datetime.datetime.now()

    # Pull a netCDF dataset object out so the projection coordinates can
    # be obtained
    ncdata = Dataset(
        '/soge-home/projects/seviri_dust/raw_seviri_data/bt_nc/'
                 + date.strftime(
                     "%B%Y") + '/North_Africa_SEVIRI_BTs_' + str(
                     date.strftime("%Y%m%d%H%M")) + '.nc')
    print 'Loading from netCDF:', datetime.datetime.now() - test4


    # Reproject the original geos projection coordinates
    test5 = datetime.datetime.now()
    lons, lats = pinkdust.reproject_to_latlon(ncdata)
    print 'Reprojecting to lat/lon:', datetime.datetime.now() - test5



    # Regrid the irregular to regular lat lon
    test7 = datetime.datetime.now()
    data_array = np.zeros((lons.shape[0], lats.shape[1], 3))
    data_array_ordered = np.zeros((lons.shape[0], lats.shape[1], 3))

    data_labels = np.array(
        [ncdata.variables['band0'][:][0], ncdata.variables['band1'][:][0],
         ncdata.variables['band2'][:][0]])

    data_array[:, :, 0] = ncdata.variables['Image0'][:]
    data_array[:, :, 1] = ncdata.variables['Image1'][:]
    data_array[:, :, 2] = ncdata.variables['Image2'][:]

    IR_indices = np.array([0, 1, 2])

    data_array_ordered[:, :, 0] = data_array[:, :,
                                  IR_indices[data_labels == 'IR_087'][0]]
    data_array_ordered[:, :, 1] = data_array[:, :,
                                  IR_indices[data_labels == 'IR_108'][0]]
    data_array_ordered[:, :, 2] = data_array[:, :,
                                  IR_indices[data_labels == 'IR_120'][0]]

    data_regridded = pinkdust.regrid_data_to_regular(lons, lats,
                                                 data_array_ordered)
    print 'Regridding to a regular lat/lon:', datetime.datetime.now() - test7



    # Pull out cloud mask data to be regridded and added to the nc file
    # glob is used for a wildcard
    #clouddata = pygrib.open(glob.glob(
    #    '/soge-home/projects/seviri_dust/raw_seviri_data/cloudmask_grib/MSG*-'
    #    'SEVI-MSGCLMK-0100-0100-' + date.strftime(
    #        "%Y%m%d%H%M%S") + '*')[0])
    #grb = clouddata.select()[0]
    #cloudmaskarray = grb.values[:, ::-1]
    #cloudlats, cloudlons = grb.latlons()

    #cloudmaskarray[cloudmaskarray >= 3] = np.nan
    #cloudlats[cloudlats > 90] = np.nan
    #cloudlons[cloudlons > 180] = np.nan

    # Generate a regular lat/lon grid for the cloud mask
    regular_lons = np.linspace(np.min(lons), np.max(lons), lons.shape[1])
    regular_lats = np.linspace(np.min(lats), np.max(lats), lats.shape[0])

    # Regrid the cloud mask to the above regular grid (note a new
    # function was defined as it was needed for a previous version of
    # the code...)
    #cloudmask_regridded = pinkdust.regrid_data(cloudlons, cloudlats,
    #                                         regular_lons,
    #                                  regular_lats, cloudmaskarray)
    time = num2date(ncdata.variables['time'][:],
                    ncdata.variables['time'].units)

    """

    f = tables.open_file('/soge-home/projects/seviri_dust/raw_seviri_data'
                         '/intermediary_files/BT_087_'+date.strftime(
        '%Y%m%d%H%M%S.hdf'),
                         'w')
    atom = tables.Atom.from_dtype(data_regridded[:, :, 0].dtype)
    filters = tables.Filters(complib='blosc', complevel=5)
    ds = f.create_carray(f.root, 'data', atom,
                         data_regridded[:, :, 0].shape,
                         filters=filters)
    ds[:] = data_regridded[:, :, 0]
    f.close()

    f = tables.open_file('/soge-home/projects/seviri_dust/raw_seviri_data'
                         '/intermediary_files/BT_108_'+date.strftime(
        '%Y%m%d%H%M%S.hdf'), 'w')
    atom = tables.Atom.from_dtype(data_regridded[:, :, 1].dtype)
    filters = tables.Filters(complib='blosc', complevel=5)
    ds = f.create_carray(f.root, 'data', atom,
                         data_regridded[:, :, 0].shape,
                         filters=filters)
    ds[:] = data_regridded[:, :, 1]
    f.close()

    f = tables.open_file('/soge-home/projects/seviri_dust/raw_seviri_data'
                         '/intermediary_files/BT_120_'+date.strftime(
        '%Y%m%d%H%M%S.hdf'), 'w')
    atom = tables.Atom.from_dtype(data_regridded[:, :, 2].dtype)
    filters = tables.Filters(complib='blosc', complevel=5)
    ds = f.create_carray(f.root, 'data', atom,
                         data_regridded[:, :, 0].shape,
                         filters=filters)
    ds[:] = data_regridded[:, :, 2]
    f.close()

    #f = tables.open_file('/soge-home/projects/seviri_dust/raw_seviri_data'
    #                     '/intermediary_files/cloudmask_
    # '+date.strftime(
    #    '%Y%m%d%H%M%S.hdf'), 'w')
    #atom = tables.Atom.from_dtype(cloudmask_regridded.dtype)
    #filters = tables.Filters(complib='blosc', complevel=5)
    #ds = f.create_carray(f.root, 'data', atom,
    #                     data_regridded[:, :, 0].shape,
    #                     filters=filters)
    #ds[:] = cloudmask_regridded
    f.close()

    ncdata.close()

if __name__ == '__main__':

    Year_lower = 2013
    Year_upper = 2013
    Month_lower = 6
    Month_upper = 6
    Day_lower = 23
    Day_upper = 23
    Hour_lower = 23
    Hour_upper = 23
    Minute_lower = 45
    Minute_upper = 45

    # Generate datetime objects corresponding to these time bounds
    time_params = np.array([Year_lower, Year_upper, Month_lower, Month_upper,
                            Day_lower, Day_upper, Hour_lower, Hour_upper,
                            Minute_lower, Minute_upper])

    datetimes = utilities.get_datetime_objects(time_params)

    wrapper(datetimes[0])

    """
    processes = [multiprocessing.Process(target=wrapper, args=(i,))
                 for i in datetimes]

    for p in processes:
        p.start()
    for p in processes:
        p.join()

    """

    # Run once to get the lats and lons
    test1 = datetime.datetime.now()
    global_data, time_slot = pinkdust.load_channels(datetimes[0])

    # Reproject it to North Africa - note that this requires a custom
    # projection in areas.def config file
    test2 = datetime.datetime.now()
    projected_data = global_data
    #projected_data = global_data.project("NorthAfrica")

    test3 = datetime.datetime.now()
    netcdf4.save(projected_data,
                 '/soge-home/projects/seviri_dust/raw_seviri_data/bt_nc/'
                 +datetimes[0].strftime(
                     "%B%Y")+'/North_Africa_SEVIRI_BTs_' + str(
        datetimes[0].strftime("%Y%m%d%H%M")) + '.nc', compression=True,
                 dtype=np.int16,
                 band_axis=2, area_aggregation=True, time_dimension=True)
    test4 = datetime.datetime.now()

    # Pull a netCDF dataset object out so the projection coordinates can
    # be obtained
    ncdata = Dataset(
        '/soge-home/projects/seviri_dust/raw_seviri_data/bt_nc/'
                 +datetimes[0].strftime(
                     "%B%Y")+'/North_Africa_SEVIRI_BTs_' + str(
        datetimes[0].strftime("%Y%m%d%H%M"))
        + '.nc')

    # Reproject the original geos projection coordinates
    test5 = datetime.datetime.now()
    lons, lats = pinkdust.reproject_to_latlon(ncdata)
    # Generate a regular lat/lon grid for the cloud mask
    regular_lons = np.linspace(np.min(lons), np.max(lons), lons.shape[1])
    regular_lats = np.linspace(np.min(lats), np.max(lats), lats.shape[0])
    print 'Obtaining lats and lons for netCDF file:', datetime.datetime.now() \
                                                     - test1

    ncdata.close()

    ncfile = pinkdust.create_time_nc_file('/soge-home/projects/seviri_dust/raw_seviri_data/bt_nc/'
                 +datetimes[0].strftime(
                     "%B%Y")+'/BT_regridded_'+datetimes[
        0].strftime(
        '%Y%m%d%H%M%S')+'_'+datetimes[
        -1].strftime(
        '%Y%m%d%H%M%S')+'.nc', datetimes, regular_lats, regular_lons)

    for i in np.arange(0, len(datetimes)):
        f = tables.open_file(
            '/soge-home/projects/seviri_dust/raw_seviri_data'
            '/intermediary_files/'
                 +datetimes[0].strftime(
                     "%B%Y")+'/BT_087_'+datetimes[i].strftime(
        '%Y%m%d%H%M%S.hdf'))
        arrobj = f.get_node('/data')
        bt_087 = arrobj.read()
        f.close()

        f = tables.open_file('/soge-home/projects/seviri_dust/raw_seviri_data'
            '/intermediary_files/'
                 +datetimes[0].strftime(
                     "%B%Y")+'/BT_108_'+datetimes[i].strftime(
        '%Y%m%d%H%M%S.hdf'))
        arrobj = f.get_node('/data')
        bt_108 = arrobj.read()
        f.close()

        f = tables.open_file('/soge-home/projects/seviri_dust/raw_seviri_data'
            '/intermediary_files/'
                 +datetimes[0].strftime(
                     "%B%Y")+'/BT_120_'+datetimes[i].strftime(
        '%Y%m%d%H%M%S.hdf'))
        arrobj = f.get_node('/data')
        bt_120 = arrobj.read()
        f.close()

        #f = tables.open_file('/soge-home/projects/seviri_dust/raw_seviri_data'
        #    '/intermediary_files/'
        #         +datetimes[0].strftime(
        #             "%B%Y")+'/cloudmask_' + datetimes[i].strftime(
        #    '%Y%m%d%H%M%S.hdf'))
        #arrobj = f.get_node('/data')
        #cloudmask = arrobj.read()
        #f.close()

        cloudmask = None

        ncfile = pinkdust.save_to_existing_nc(ncfile, bt_087, bt_108,
                                 bt_120,
                            cloudmask, datetimes[i])

    print 'Done'
    ncfile.close()


