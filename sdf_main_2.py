import numpy as np
import datetime
from netCDF4 import Dataset
import tables

import utilities
import sdf
import plotting

if __name__ == '__main__':

    # Select the time parameters within which to parse input nc files and
    # obtain SDF values

    # This needs to be a period for which there is a 15 day time window
    # around it
    year_lower = 2012
    year_upper = 2012
    month_lower = 6
    month_upper = 6
    day_lower = 10
    day_upper = 13
    hour_lower = 0
    hour_upper = 0
    minute_lower = 0
    minute_upper = 15

    # Generate the full set of datetimes for the entire period
    time_params = np.array([year_lower, year_upper, month_lower,
                            month_upper, day_lower, day_upper,
                            hour_lower, hour_upper, minute_lower,
                            minute_upper])

    SDF_datetimes = utilities.get_datetime_objects(time_params)
    SDF_datestrings = [j.strftime("%Y%m%d%H%M") for j in SDF_datetimes]

    ### CLOUD MASKING PHASE ###

    # Generate a set of datetimes for one day, and loop through just one day
    time_params_oneday = np.array([year_lower, year_lower, month_lower,
                                   month_lower, day_lower, day_lower,
                                   hour_lower, hour_upper,
                                   minute_lower, minute_upper])

    # Reads in lats and lons
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

    oneday_datetimes = utilities.get_datetime_objects(time_params_oneday)

    print 'Generating cloud masks for all timesteps'
    # Loop through each time in a single day
    for i in np.arange(0, len(oneday_datetimes)):
        print oneday_datetimes[i]
        window_datetime_lower = datetime.datetime(year_lower, month_lower,
                                                  day_lower, hour_lower,
                                                  minute_lower) \
                                                  - datetime.timedelta(days=7)
        window_datetime_upper = datetime.datetime(year_upper, month_upper,
                                                  day_upper, hour_upper,
                                                  minute_upper) \
                                                  + datetime.timedelta(days=7)

        # Get datetime objects between the above bounds
        time_params_7dayw = np.array([window_datetime_lower.year,
                                window_datetime_upper.year,
                                window_datetime_lower.month,
                                window_datetime_upper.month,
                                window_datetime_lower.day,
                                window_datetime_upper.day,
                                oneday_datetimes[i].hour,
                                oneday_datetimes[i].hour,
                                oneday_datetimes[i].minute,
                                oneday_datetimes[i].minute])
        datetimes_7dayw = utilities.get_daily_datetime_objects(
            time_params_7dayw)

        bt_15day = np.zeros((datetimes_7dayw.shape[0], 3,
                         lats.shape[0],
                         lons.shape[1]))

        cloudscreenedbt_15day = np.zeros((datetimes_7dayw.shape[0], 3,
                         lats.shape[0],
                         lons.shape[1]))

        # Loop through each day of the time window for this time of day
        for j in np.arange(0, len(datetimes_7dayw)):
            date = datetimes_7dayw[j]

            # Extract BT data for this timestep
            filename = '/ouce-home/data/satellite/meteosat/seviri/15-min/' \
                       'native/bt/nc' \
                       '/' + str(date.strftime(
                '%B').upper()) + str(date.year) + \
                       '/H-000-MSG2__-MSG2________-' \
                            'IR_BrightnessTemperatures___-000005___-' + str(
                date.strftime('%Y')) + str(date.strftime('%m')) + \
                       str(date.strftime('%d')) + str(date.strftime('%H')) \
                       + str(date.strftime('%M')) + '-__.nc'

            try:
                btdata = Dataset(filename, 'r')
            except:
                print 'Found no BT data for ' + filename
                continue

            # Extract cloud mask data for this timestep
            cloudfilename = '/ouce-home/data/satellite/meteosat/' \
                            'seviri/15-min/'\
                            'native/cloudmask/nc' \
                            '/' + str(date.strftime('%B').upper()) \
                            + str(date.year) + '_CLOUDS/eumetsat.cloud.' \
                            + str(date.strftime('%Y')) \
                            + str(date.strftime('%m')) \
                            + str(date.strftime('%d')) + str(
                date.strftime('%H')) \
                            + str(date.strftime('%M')) + '.nc'
            try:
                clouddata = Dataset(cloudfilename, 'r')
            except:
                print 'Found no cloud data for ' + cloudfilename
                continue

            # Apply cloud screening
            cloudscreenedbt_15day[i], bt_15day[i, 0, :, :], bt_15day[i, 1,
                                                                :, :], \
                bt_15day[i, 2, :, :] = sdf.cloud_screen(btdata, clouddata)
            btdata.close()

            f = tables.open_file('cloud_masked_bt_15d_'+
                                str(oneday_datetimes[i].hour)+str(
                oneday_datetimes[i].minute)+'.hdf', 'w')
            atom = tables.Atom.from_dtype(cloudscreenedbt_15day.dtype)
            filters = tables.Filters(complib='blosc', complevel=5)
            ds = f.create_carray(f.root, 'somename', atom,
                                cloudscreenedbt_15day.shape,
                                filters=filters)
            ds[:] = cloudscreenedbt_15day
            f.close()

            # Save cloud masked data for this time of day to file
            f = tables.open_file('bt_15d_' +
                                str(oneday_datetimes[i].hour) + str(
                oneday_datetimes[i].minute) + '.hdf', 'w')
            atom = tables.Atom.from_dtype(bt_15day.dtype)
            filters = tables.Filters(complib='blosc', complevel=5)
            ds = f.create_carray(f.root, 'somename2', atom,
                                bt_15day.shape,
                                filters=filters)
            ds[:] = bt_15day
            f.close()


    ### SDF PHASE ###

    # Open an nc file in which to store all SDF data for the run period

    print 'Generating SDFs'
    # Loop through the full set of datetimes for the entire period
    for i in np.arange(0, len(SDF_datetimes)):
        print SDF_datetimes[i]

        # Read in the cloud mask data for that time of day from file
        cloudscreenedbt_15day = np.load('cloud_masked_bt_15d_'+
                                        str(SDF_datetimes[i].hour)+str(
                                        SDF_datetimes[i].minute)+'.npy')

        bt_15day = np.load('bt_15d_'+str(SDF_datetimes[i].hour)+str(
                                        SDF_datetimes[i].minute)+'.npy')

        window_datetime_lower = datetime.datetime(year_lower, month_lower,
                                                  day_lower, hour_lower,
                                                  minute_lower) \
                                - datetime.timedelta(days=7)
        window_datetime_upper = datetime.datetime(year_upper, month_upper,
                                                  day_upper, hour_upper,
                                                  minute_upper) \
                                + datetime.timedelta(days=7)

        BT_15_day_lower_bound = SDF_datetimes[i] - datetime.timedelta(days=7)
        BT_15_day_upper_bound = SDF_datetimes[i] + datetime.timedelta(days=7)

        # Get datetime objects between the above bounds
        time_params_7dayw = np.array([window_datetime_lower.year,
                                      window_datetime_upper.year,
                                      window_datetime_lower.month,
                                      window_datetime_upper.month,
                                      window_datetime_lower.day,
                                      window_datetime_upper.day,
                                      SDF_datetimes[i].hour,
                                      SDF_datetimes[i].hour,
                                      SDF_datetimes[i].minute,
                                      SDF_datetimes[i].minute])
        datetimes_7dayw = utilities.get_daily_datetime_objects(
            time_params_7dayw)

        BT_15_days_screened = \
            cloudscreenedbt_15day[np.asarray([j >= BT_15_day_lower_bound
                                              and j <= BT_15_day_upper_bound
                                              for j in datetimes_7dayw])]

        # Get the mean of the 15 day window of cloud screened data
        cloudscreenedbt_15daymean_108_87 = \
            sdf.extract_15day_mean(BT_15_days_screened)

        dust_mask_108, \
        dust_mask_108_87, \
        dust_mask_120_108, \
        dust_mask_108_87_anom_screened = \
            sdf.generate_dust_mask(
                bt_15day[datetimes_7dayw == SDF_datetimes[i]],
                cloudscreenedbt_15day[datetimes_7dayw == SDF_datetimes[i]],
                cloudscreenedbt_15daymean_108_87)

        ofilename = 'SDF_' + str(SDF_datetimes[i].strftime('%Y')) + \
                    str(SDF_datetimes[i].strftime('%m')) + \
                    str(SDF_datetimes[i].strftime('%d')) + str(
            SDF_datetimes[i].strftime('%H')) + \
                    str(SDF_datetimes[i].strftime('%M')) + '.nc'

        SDF = sdf.generate_SDF(dust_mask_108,
                               dust_mask_108_87,
                               dust_mask_120_108,
                               dust_mask_108_87_anom_screened,
                               ofilename, SDF_datetimes[i], lats, lons)

    # Subset the cloud mask data so it is only the 15 day time window around
    # the given timestep

    # Do dust masking

    # Save the SDF to file