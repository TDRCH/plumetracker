import matplotlib
matplotlib.use('Agg')
from mpop.satellites import GeostationaryFactory
import datetime
import numpy as np
from netCDF4 import Dataset
from netCDF4 import date2num
from pycoast import ContourWriterAGG
from PIL import Image
from pyproj import Proj
from pyresample.geometry import SwathDefinition
from pyresample.kd_tree import resample_nearest

def load_channels(datetime):
    """
    Load channel data into an mpop scene object
    :param year: desired year integer
    :param month: desired month integer
    :param day: desired day integer
    :param hour: desired day integer
    :param minute: desired minute integer
    :param second: desired second integer
    :return global_data: an mpop scene object with data from IR channels
    12.0, 10.8, 8.7 for the desired time
    :return time_slot: a datetime object corresponding to the desired timestep
    """
    # Create a datetime object
    time_slot = datetime

    # Create a scene object and load channels
    global_data = GeostationaryFactory.create_scene("Meteosat-9", "", "seviri",
                                                    time_slot)
    global_data.load([12.0, 10.8, 8.7])
    return global_data, time_slot


def generate_image(global_data, filename):
    """
    This function generates an RGB composite directly from a scene with
    loaded channels and saves to a GeoTiff
    :param global_data: an mpop scene object with data from IR channels
    12.0, 10.8, 8.7
    :param filename: a string with the name of the image file to be generated
    """

    # Generate a dust composite with Pytroll inbuilt function
    img = global_data.image.dust()

    # Save the image initially so it can be read in as an image object
    img.save("/ouce-home/students/hert4173/SEVIRI_imagery/" + filename)

    # Read it back in as an image object
    img = Image.open("/ouce-home/students/hert4173/SEVIRI_imagery/" + filename)

    # Define projection parameters
    proj4_string = '+proj=geos +lon_0=0.0 +a=6378169.00 +b=6356583.80 ' \
                   '+h=35785831.0'
    area_extent = (-2000000, 0, 3000000, 4000000)
    area_def = (proj4_string, area_extent)

    # ContourWriterAGG here requires the 'aggdraw' package
    cw = ContourWriterAGG(
        '/ouce-home/students/hert4173/.conda_envs/virtual_env/lib/python2.7/'
        'pycoast/GSHHS_DATA_ROOT')
    cw.add_coastlines(img, area_def, resolution='i', width=3)
    cw.add_borders(img, area_def, outline=(255, 255, 255), resolution='i',
                   width=3)
    img.save("/ouce-home/students/hert4173/SEVIRI_imagery/" + filename)


def generate_image_from_array(data_array):
    """
    Generate a pink dust image from a SEVIRI netCDF file with channels 12.0,
    10.8, 8.7
    :param ncdata: a netCDF4 dataset object containing SEVIRI data
    :param filename: string containing the name of the output jpg image
    :return img: a GeoImage object
    """

    # Read in brightness temperature channels
    IR_087 = data_array[:, :, 0]
    IR_108 = data_array[:, :, 1]
    IR_120 = data_array[:, :, 2]

    # Parameters for pink dust formula from Brindley et al. (2012)
    MinR = -4.0
    MaxR = 2.0
    MinG = 0.0
    MaxG = 15.0
    MinB = 261.0
    MaxB = 289.0
    GammaR = 1.0
    GammaG = 2.5
    GammaB = 1.0

    # Masked array elements are replaced with zero
    imgR = np.ma.filled(
        (255 * (((IR_120 - IR_108) - MinR) / (MaxR - MinR)) ** (1.0 / GammaR)),
        fill_value=0)
    imgG = np.ma.filled(
        (255 * (((IR_108 - IR_087) - MinG) / (MaxG - MinG)) ** (1.0 / GammaG)),
        fill_value=0)
    imgB = np.ma.filled(
        (255 * ((IR_108 - MinB) / (MaxB - MinB)) ** (1.0 / GammaB)),
        fill_value=0)

    # Elements outside the RGB range are set to zero or 255
    imgR[imgR < 0] = 0
    imgR[imgR > 255] = 255
    imgG[imgG < 0] = 0
    imgG[imgG > 255] = 255
    imgB[imgB < 0] = 0
    imgB[imgB > 255] = 255

    # All values are set to 'uint8' form from which an image can be generated
    rgbArray = np.zeros((imgR.shape[0], imgR.shape[1], 3), 'uint8')
    rgbArray[:, :, 0] = imgR
    rgbArray[:, :, 1] = imgG
    rgbArray[:, :, 2] = imgB

    return rgbArray


def get_datetime_objects(time_params):
    """
    Generates an array of datetime objects at 15 minute intervals within the
    input time bounds
    :param time_params: an array of integers corresponding to time bounds
    :return datetimes: array of datetime objects
    """

    # Get the lower and upper bound of the datetime objects
    datetime_lower = datetime.datetime(time_params[0], time_params[2],
                                       time_params[4], time_params[6],
                                       time_params[8])
    datetime_upper = datetime.datetime(time_params[1], time_params[3],
                                       time_params[5], time_params[7],
                                       time_params[9])

    # Get the difference between the two dates in minutes using a timedelta
    # object
    td = datetime_upper - datetime_lower
    mins_difference = divmod(td.days * 86400 + td.seconds, 60)[0]

    # Get an array of minute values with 15 minute intervals
    minutes_range = np.arange(0, mins_difference + 15, 15)

    # Get an array of datetime objects for each 15 minute interval
    datetimes = np.array(
        [datetime_lower + datetime.timedelta(minutes=i) for i in
         minutes_range])

    return datetimes


def reproject_to_latlon(nc):
    """
    Reprojects geostationary projection coordinates to latitude and longitude
    :param nc: a netCDF dataset object
    :return: longitude and latitude values
    """

    # Get the satellite viewing height
    h = nc.variables['grid_mapping_0'].perspective_point_height

    # Define a geostationary projection
    p = Proj(proj='geos', h=h, lon_0=0)

    # Longitude values are duplicated along the y-axis to create a 2D array
    # (note, should use np.meshgrid)
    x = np.resize(nc.variables['x0'][:] * h, (2400, 3000))
    # Latitude values are duplicated along the x-axis to create a 2D array
    y = np.array([nc.variables['y0'][:] * h] * 3000).transpose()

    # Reproject
    lons, lats = p(x, y, inverse=True)

    return lons, lats


def regrid_data(lons, lats, target_lons, target_lats, array):
    """
    Regrids the irregular lons and lats produced by the reprojection so that
    imshow can handle it
    :param lons: an array of irregular longitude values produced by the
    reprojection function
    :param lats: an array of irregular latitude values produced by the
    reprojection function
    :param array: an array BTs or RGB values for each pixel
    :return: the array on an interpolated regular grid
    """

    # Define a target grid
    XI = target_lons
    YI = target_lats
    XI, YI = np.meshgrid(XI, YI)

    # Resample BT data
    def_a = SwathDefinition(lons=XI, lats=YI)
    def_b = SwathDefinition(lons=lons, lats=lats)
    interp_dat = resample_nearest(def_b, array, def_a,
                                  radius_of_influence=7000)

    return interp_dat

def regrid_data_to_regular(lons, lats, array):
    """
    Regrids the irregular lons and lats produced by the reprojection so that
    imshow can handle it
    :param lons: an array of irregular longitude values produced by the
    reprojection function
    :param lats: an array of irregular latitude values produced by the
    reprojection function
    :param array: an array BTs or RGB values for each pixel
    :return: the array on an interpolated regular grid
    """

    # Define the average interval in the existing lons and lats
    intervals_lon = (np.max(lons) - np.min(lons)) / len(lons[0])
    intervals_lat = (np.max(lats) - np.min(lats)) / len(lats)

    # Define a regular grid
    XI = np.arange(np.min(lons), np.max(lons), intervals_lon)
    YI = np.arange(np.min(lats), np.max(lats), intervals_lat)
    XI, YI = np.meshgrid(XI, YI)

    # Resample BT data
    def_a = SwathDefinition(lons=XI, lats=YI)
    def_b = SwathDefinition(lons=lons, lats=lats)
    interp_dat_0 = resample_nearest(def_b, array[:, :, 0], def_a,
                                    radius_of_influence=7000)
    interp_dat_1 = resample_nearest(def_b, array[:, :, 1], def_a,
                                    radius_of_influence=7000)
    interp_dat_2 = resample_nearest(def_b, array[:, :, 2], def_a,
                                    radius_of_influence=7000)

    # Recombine into a single array
    interp_dat_BTs = array
    interp_dat_BTs[:, :, 0] = interp_dat_0
    interp_dat_BTs[:, :, 1] = interp_dat_1
    interp_dat_BTs[:, :, 2] = interp_dat_2

    return interp_dat_BTs

def save_regridded_data_to_nc(filename, array_87, array_108, array_120,
                              cloudmask_regridded, lons, lats, date):
    """
    Takes BT values on a regular grid and writes them to an ncfile.
    :param filename: string containing the file name to write to
    :param array_87: array with data for BT channel 8.7 on a regular grid
    :param array_108: array with data for BT channel 10.8 on a regular grid
    :param array_120: array with data for BT channel 12.0 on a regular grid
    :param lons: 1D array with regular lons
    :param lats: 1D array with regular lats
    """
    ncfile = Dataset(filename, 'w', format='NETCDF4_CLASSIC')
    ncfile.description = 'Brightness temperature values from Meteosat ' \
                         'SEVIRI, reprojected to lat/lon from Geostationary.'
    time = ncfile.createDimension('time', None)
    lat = ncfile.createDimension('lat', len(lats))
    lon = ncfile.createDimension('lon', len(lons))
    times = ncfile.createVariable('time', np.float64, ('time'))
    latitudes = ncfile.createVariable('latitude', np.float32, ('lat'))
    longitudes = ncfile.createVariable('longitude', np.float32, ('lon'))
    channel_087 = ncfile.createVariable('channel_087', np.float32,
                                        ('lat', 'lon'))
    channel_108 = ncfile.createVariable('channel_108', np.float32,
                                        ('lat', 'lon'))
    channel_120 = ncfile.createVariable('channel_120', np.float32,
                                        ('lat', 'lon'))
    cloudmask = ncfile.createVariable('cloud_mask', np.float32, ('lat', 'lon'))
    latitudes.units = 'degrees_north'
    longitudes.units = 'degrees_east'
    times.units = 'hours since 0001-01-01 00:00:00'
    times.calendar = 'gregorian'
    latitudes[:] = lats
    longitudes[:] = lons
    times[:] = date2num(date, times.units)
    channel_087[:] = array_87
    channel_108[:] = array_108
    channel_120[:] = array_120
    cloudmask[:] = cloudmask_regridded
    # Forgot to write to time coordinate variable
    ncfile.close()
