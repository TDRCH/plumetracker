"""
A python version of the SEVIRI Dust Flag algorithm v2 originally written in shellscript by Ian Ashpole (ported to python by Thomas Caton Harrison)
"""

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import datetime
from netCDF4 import Dataset
from netCDF4 import date2num
from mpl_toolkits.basemap import Basemap
import os
import tables

import pinkdust
import utilities

def cloud_screen(btdata, clouddata):
    """
    Cloud screens SEVIRI data using EUMETSAT cloud screen
    :param btdata: an nc Dataset object with the three SEVIRI bands for 10.8,
    12.0, 8.7 BTs
    :param cloudscreendata: an nc Dataset object from EUMETSAT's SEVIRI cloud
    screen
    :return cloudscreenedbt: cloud screened SEVIRI band array for the given
    day
    """
    cloudscreen = clouddata.variables['cmask'][:]
    cloudscreen_bool = np.asarray(cloudscreen != 2)
    # Note that the below assumes that the lat/long coord variables are not
    # in mesh grid form
    seviribands = np.zeros((3, cloudscreen.shape[0],
                            cloudscreen.shape[1]))
    channel_087 = btdata.variables['bt087'][:]
    channel_108 = btdata.variables['bt108'][:]
    channel_120 = btdata.variables['bt120'][:]

    seviribands[0] = channel_087
    seviribands[1] = channel_108
    seviribands[2] = channel_120
    seviribands[0][cloudscreen_bool == False] = np.nan
    seviribands[1][cloudscreen_bool == False] = np.nan
    seviribands[2][cloudscreen_bool == False] = np.nan

    return seviribands, channel_087, channel_108, channel_120

def extract_15day_mean(cloudscreenedbt_15day):
    """
    Finds the 15 day time mean of the BT10.8-8.7 channels around a given day
    on which SDF is being calculated
    :param cloudscreenedbt_15day: cloud screened array with three BT channels
    for the 15 day window
    :return: 2D array with a 15 day time mean of the BT10.8-8.7 channels
    """
    BT_108_87 = cloudscreenedbt_15day[:, 1] - cloudscreenedbt_15day[:, 0]
    cloudscreenedbt_15daymean_108_87 = np.nanmean(BT_108_87, axis=0)
    return cloudscreenedbt_15daymean_108_87

def generate_dust_mask(noncloudscreenedarray, cloudscreenedarray,
                       mean_15day_108_87):
    """
    Finds the dust mask value for each of the R, G and B channels
    :param noncloudscreenedarray: array of BT channel values with no cloud
    screen applied
    :param cloudscreenedarray: array of BT channel values with a cloud screen
    applied
    :param mean_15day_108_87: array with the 15-day time mean around the
    chosen day
    :return: arrays with dust masks for each of the channels
    """
    BT_87 = noncloudscreenedarray[0, 0]
    BT_108 = noncloudscreenedarray[0, 1]
    BT_120 = noncloudscreenedarray[0, 2]

    BT_array = np.zeros((BT_87.shape[0], BT_87.shape[1], 3))
    BT_array[:, :, 0] = BT_87
    BT_array[:, :, 1] = BT_108
    BT_array[:, :, 2] = BT_120

    rgb_array = pinkdust.generate_image_from_array(BT_array)
    """
    extent = (-26, 44, 1, 43)
    m = Basemap(projection='cyl', llcrnrlon=extent[0], urcrnrlon=extent[1],
                llcrnrlat=extent[2], urcrnrlat=extent[3],
                resolution='i')
    m.drawcoastlines(linewidth=0.5)
    m.drawcountries(linewidth=0.5)
    m.imshow(rgb_array, extent=extent, origin='lower')
    plt.savefig('/ouce-home/students/hert4173/SEVIRI_imagery/'
                'SDF_test_pinkdust.png',
                dpi=200)
    """
    BT_108_screened = cloudscreenedarray[0, 1]
    BT_87_screened = cloudscreenedarray[0, 0]

    # CORRECTION - ONLY THE ANOMALY BTs SHOULD BE CLOUD SCREENED

    # levels = np.arange(200, 360, 10)

    # plt.close()
    # plt.contourf(lons, lats, BT_87, levels)
    # plt.colorbar()
    # plt.savefig('BT_87.png')

    # plt.close()
    # plt.contourf(lons, lats, BT_108, levels)
    # plt.colorbar()
    # plt.savefig('BT_108.png')

    # plt.close()
    # plt.contourf(lons, lats, BT_120, levels)
    # plt.colorbar()
    # plt.savefig('BT_120.png')

    BT_108_87 = BT_108 - BT_87
    BT_108_87_screened = BT_108_screened - BT_87_screened
    BT_120_108 = BT_120 - BT_108

    # plt.close()
    # plt.contourf(lons, lats, BT_108_87)
    # plt.colorbar()
    # plt.savefig('BT_108_87.png')

    # plt.close()
    # plt.contourf(lons, lats, BT_120_108)
    # plt.colorbar()
    # plt.savefig('BT_120_108.png')

    # Calculate dust masks
    dust_mask_108 = 1 * (BT_108 >= 285)
    dust_mask_120_108 = 1 * (BT_120_108 >= 0)
    dust_mask_108_87 = 1 * (BT_108_87 <= 10)

    # Calculate anomalies between BT108(t) and BT108(15-day mean), and screen
    # them
    dust_mask_108_87_anom = BT_108_87 - mean_15day_108_87
    dust_mask_108_87_anom_screened = 1 * (dust_mask_108_87_anom <= -2)

    # plt.close()
    # plt.contourf(lons, lats, dust_mask_108_87_anom)
    # plt.colorbar()
    # plt.savefig('BT_mean_15day_108_87_anom.png')

    return dust_mask_108, dust_mask_108_87, dust_mask_120_108, \
           dust_mask_108_87_anom_screened


def generate_SDF(dust_mask_108, dust_mask_108_87, dust_mask_120_108,
                 dust_mask_108_87_anom_screened, filename, datetime,
                 lats, lons):
    """
    Calculates the SDF using dust masks for the R, G and B channels as well
    as 108-87 15 day BT anomalies, then saves to an nc file
    :param dust_mask_108: an array containing the dust mask values for the
    10.8 BT channel
    :param dust_mask_108_87: an array containing dust mask values for the
    10.8-8.7 channel
    :param dust_mask_120_108: an array containing dust mask values for the
    12.0-10.8 channel
    :param dust_mask_108_87_anoms: an array containing the anomalies in the
    dust mask values for 10.8-8.7 channel
    :return: SDF array, which has been written to an nc file
    """
    # Calculate SDF from the sum of dust mask values
    dust_mask_sum = dust_mask_108 + dust_mask_108_87 + dust_mask_120_108 \
                    + dust_mask_108_87_anom_screened
    # plt.contourf(lons, lats, dust_mask_sum)
    # plt.colorbar()
    # plt.savefig('dust_mask_sum.png')
    SDF = 1 * (dust_mask_sum >= 4)

    """
    # Write an output nc file with an SDF variable

    # File description and dimensions
    SDFfile = Dataset(filename, 'w', format='NETCDF4_CLASSIC')
    SDFfile.description = 'SDF values from EUMETSAT SEVIRI on a ' \
                          'regular lat/lon grid.'
    SDFfile.createDimension('time', None)
    SDFfile.createDimension('lat', lats.shape[0])
    SDFfile.createDimension('lon', lons.shape[1])

    # Variables
    times = SDFfile.createVariable('time', np.float64, ('time'))
    SDFvar = SDFfile.createVariable('SDF', np.int32, ('lat', 'lon'))

    # Attributes
    times.units = 'hours since 0001-01-01 00:00:00'
    times.calendar = 'gregorian'

    # Writing data to file and closing
    times[:] = date2num(datetime, units=times.units, calendar=times.calendar)
    SDFvar[:] = SDF
    SDFfile.close()
    """

    return SDF