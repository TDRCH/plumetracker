import numpy as np
import datetime
from netCDF4 import Dataset

import utilities
import sdf
import plotting

if __name__ == '__main__':

    # Select the time parameters within which to parse input nc files and
    # obtain SDF values
    year_lower = 2012
    year_upper = 2012
    month_lower = 7
    month_upper = 7
    day_lower = 1
    day_upper = 30
    hour_lower = 0
    hour_upper = 23
    minute_lower = 0
    minute_upper = 45

    test1 = datetime.datetime.now()

    # Calculates the cloud screened BT between the earliest time -7 and the
    # latest time +7 (so we only read in all the data we need once)
    window_datetime_lower = datetime.datetime(year_lower, month_lower,
                                              day_lower, hour_lower,
                                              minute_lower) \
                            - datetime.timedelta(
        days=7)
    window_datetime_upper = datetime.datetime(year_upper, month_upper,
                                              day_upper, hour_upper,
                                              minute_upper) \
                            + datetime.timedelta(
        days=7)

    # Get datetime objects between the above bounds
    time_params = np.array([window_datetime_lower.year,
                            window_datetime_upper.year,
                            window_datetime_lower.month,
                            window_datetime_upper.month,
                            window_datetime_lower.day,
                            window_datetime_upper.day,
                            window_datetime_lower.hour,
                            window_datetime_upper.hour,
                            window_datetime_lower.minute,
                            window_datetime_upper.minute])
    datetimes = utilities.get_daily_datetime_objects(time_params)
    datestrings = [j.strftime("%Y%m%d%H%M") for j in datetimes]

    # Get datetime objects without the 7 days on the ends
    time_params = np.array([year_lower, year_upper, month_lower,
                            month_upper, day_lower, day_upper,
                            hour_lower, hour_upper, minute_lower,
                            minute_upper])

    SDF_datetimes = utilities.get_datetime_objects(time_params)
    SDF_datestrings = [j.strftime("%Y%m%d%H%M") for j in SDF_datetimes]

    print 'Obtaining datetimes for the chosen time window:',\
        datetime.datetime.now() - test1

    test2 = datetime.datetime.now()
    # Read in the SEVIRI file one time to extract some lats and lons

    lonlats = Dataset(
        '/ouce-home/data/satellite/meteosat/seviri/15-min/native/'
        'lonlats.NA_MiddleEast.nc')

    # These need to be regridded to regular for consistency with cloud mask
    lons = lonlats.variables['longitude'][:]
    lats = lonlats.variables['latitude'][:]
    lonmask = lons > 360
    latmask = lats > 90
    lons = np.ma.array(lons, mask=lonmask)
    lats = np.ma.array(lats, mask=latmask)

    cloudscreenedbt_15day = np.zeros((datetimes.shape[0], 3,
                                      lats.shape[0],
                                      lons.shape[1]))

    bt_15day = np.zeros((datetimes.shape[0], 3,
                         lats.shape[0],
                         lons.shape[1]))

    cloudscreenedbt_15day[:] = np.nan
    bt_15day[:] = np.nan

    print 'Reading in SEVIRI file for lats/lons and initialising arrays:', \
        datetime.datetime.now() - test2
    # Set an initial value for BT data
    test3 = datetime.datetime.now()
    # Generate cloud masked BTs for each day, including +- 7 days on the ends
    # k is initiated here at -1 so it'll go to 0 on the first round of the
    # for loop
    k = -1
    for i in datetimes:
        k += 1
        filename = '/ouce-home/data/satellite/meteosat/seviri/15-min/' \
                   'native/bt/nc' \
                   '/'+str(i.strftime(
            '%B').upper())+str(i.year)+'/H-000-MSG2__-MSG2________-' \
                   'IR_BrightnessTemperatures___-000005___-' + str(
            i.strftime('%Y')) \
                   + str(i.strftime('%m')) \
                   + str(i.strftime('%d')) + str(
            i.strftime('%H')) \
                   + str(i.strftime('%M')) + '-__.nc'
        try:
            btdata = Dataset(filename, 'r')
        except:
            print 'Found no BT data for ' + filename
            continue
        cloudfilename = '/ouce-home/data/satellite/meteosat/seviri/15-min/' \
                   'native/cloudmask/nc' \
                   '/' + str(i.strftime(
            '%B').upper()) + str(i.year) + \
                        '_CLOUDS/eumetsat.cloud.' + str(
            i.strftime('%Y')) \
                   + str(i.strftime('%m')) \
                   + str(i.strftime('%d')) + str(
            i.strftime('%H')) \
                   + str(i.strftime('%M')) + '.nc'
        try:
            clouddata = Dataset(cloudfilename, 'r')
        except:
            print 'Found no cloud data for ' + cloudfilename
            continue
        # Actually, would like to be getting the cloud screen data from
        # the original BTs nc file, and it needs to be on the same lat
        # lon grid as BTs!!
        # How to do that: read in the original cloud mask dataset in
        # pink_dust and find a way of matching up the cloud mask
        # (without interpolating it since it's discrete values, i.e.
        # nearest neighbour)
        cloudscreenedbt_15day[k], bt_15day[k, 0, :, :], bt_15day[k, 1, :, :], \
        bt_15day[k, 2, :, :] = sdf.cloud_screen(btdata, clouddata)
        btdata.close()
    print 'Cloud screening 15 days of data:', datetime.datetime.now() - test3

    # Loop through each day in the window except the +-7 days on the ends
    test4 = datetime.datetime.now()
    for i in SDF_datetimes:
        BT_15_day_lower_bound = i - datetime.timedelta(days=7)
        BT_15_day_upper_bound = i + datetime.timedelta(days=7)
        BT_15_days_screened = \
            cloudscreenedbt_15day[np.asarray([j >= BT_15_day_lower_bound
                                              and j <= BT_15_day_upper_bound
                                              for j in datetimes])]
        cloudscreenedbt_15daymean_108_87 = \
            sdf.extract_15day_mean(BT_15_days_screened)
        dust_mask_108, \
        dust_mask_108_87, \
        dust_mask_120_108, \
        dust_mask_108_87_anom_screened = sdf.generate_dust_mask(
            bt_15day[datetimes == i],
            cloudscreenedbt_15day[datetimes == i],
            cloudscreenedbt_15daymean_108_87)
        ofilename = 'SDF_' + str(i.strftime('%Y')) + \
                    str(i.strftime('%m')) + \
                    str(i.strftime('%d')) + str(
            i.strftime('%H')) + str(i.strftime('%M')) + '.nc'
        SDF = sdf.generate_SDF(dust_mask_108,
                           dust_mask_108_87,
                           dust_mask_120_108,
                           dust_mask_108_87_anom_screened,
                           ofilename, i, lats, lons)
    print 'Running SDF algorithm for chosen date:', \
        datetime.datetime.now() - test4

    test5 = datetime.datetime.now()
    """
    plt.close()
    extent = (np.min(lons), np.max(lons), np.min(lats), np.max(lats))
    m = Basemap(projection='cyl', llcrnrlon=extent[0], urcrnrlon=extent[1],
                llcrnrlat=extent[2], urcrnrlat=extent[3],
                resolution='i')
    lons, lats = np.meshgrid(lons, lats)
    m.drawcoastlines(linewidth=0.5)
    m.drawcountries(linewidth=0.5)
    m.contourf(lons, lats, SDF)
    # m.contourf(cloudlons, cloudlats, cloudmaskarray)
    plt.colorbar()
    plt.savefig('SDF_test_25_06_2016.png')
    plt.close()
    print 'Plotting:', datetime.datetime.now() - test5
    """
    # SDF tests:

    # 11/06/2010 1200
    # 08/07/2016 1500
    # 25/06/2016 1700

    plotting.plot_existing_SDF()