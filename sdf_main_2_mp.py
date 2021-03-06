import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import datetime
from netCDF4 import Dataset
import multiprocessing
import tables
from netCDF4 import date2num
import os
import smtplib

import utilities
import sdf
import plotting
import pinkdust

def cloud_mask_mw(i):
    """
    Moving window cloud masking to be used by multiprocessing
    :return:
    """

    date = oneday_datetimes[i]
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
                                  date.hour,
                                  date.hour,
                                  date.minute,
                                  date.minute])
    datetimes_7dayw = utilities.get_daily_datetime_objects(
        time_params_7dayw)

    bt_15day = np.zeros((datetimes_7dayw.shape[0], 3,
                         ianlats.shape[0],
                         ianlons.shape[0]))

    #bt_15day[:] = np.nan

    cloudscreenedbt_15day = np.zeros((datetimes_7dayw.shape[0], 3,
                                      ianlats.shape[0],
                                      ianlons.shape[0]))

    #cloudscreenedbt_15day[:] = np.nan

    f = tables.open_file(
        '/soge-home/projects/seviri_dust/sdf/intermediary_files/'
        'cloud_masked_bt_15d_' + oneday_datetimes[i].strftime(
                "%Y_%H_%M") + '.hdf', 'w')
    atom = tables.Atom.from_dtype(cloudscreenedbt_15day.dtype)
    filters = tables.Filters(complib='blosc', complevel=5)
    cs = f.create_carray(f.root, 'data', atom,
                         cloudscreenedbt_15day.shape,
                         filters=filters)

    g = tables.open_file(
        '/soge-home/projects/seviri_dust/sdf/intermediary_files/bt_15d_' + oneday_datetimes[i].strftime(
                "%Y_%H_%M") + '.hdf', 'w')
    atom = tables.Atom.from_dtype(bt_15day.dtype)
    filters = tables.Filters(complib='blosc', complevel=5)
    bts = g.create_carray(g.root, 'data', atom,
                          bt_15day.shape,
                          filters=filters)

    # Loop through each day of the time window for this time of day
    for j in np.arange(0, len(datetimes_7dayw)):
        print str(oneday_datetimes[i].hour) + str(
            oneday_datetimes[i].minute), 'day ', j + 1
        date_w = datetimes_7dayw[j]

        # Extract BT data for this timestep
        filename = '/ouce-home/data/satellite/meteosat/seviri/15-min/' \
                   'native/bt/nc' \
                   '/' + str(date_w.strftime(
            '%B').upper()) + str(date_w.year) + \
                   '/H-000-MSG2__-MSG2________-' \
                   'IR_BrightnessTemperatures___-000005___-' + str(
            date_w.strftime('%Y')) + str(date_w.strftime('%m')) + \
                   str(date_w.strftime('%d')) + str(date_w.strftime('%H')) \
                   + str(date_w.strftime('%M')) + '-__.nc'

        try:
            btdata = Dataset(filename, 'r')
        except:
            print 'Found no BT data for ' + filename
            continue

        # Extract cloud mask data for this timestep
        cloudfilename = '/ouce-home/data/satellite/meteosat/' \
                        'seviri/15-min/' \
                        'native/cloudmask/nc' \
                        '/' + str(date_w.strftime('%B').upper()) \
                        + str(date_w.year) + '_CLOUDS/eumetsat.cloud.' \
                        + str(date_w.strftime('%Y')) \
                        + str(date_w.strftime('%m')) \
                        + str(date_w.strftime('%d')) + str(
            date_w.strftime('%H')) \
                        + str(date_w.strftime('%M')) + '.nc'
        try:
            clouddata = Dataset(cloudfilename, 'r')
        except:
            print 'Found no cloud data for ' + cloudfilename
            continue

        # Apply cloud screening
        cloudscreenedbt_15day_array, bt087, bt108, \
        bt120 = sdf.cloud_screen(btdata, clouddata)

        regcs_array = np.zeros((len(cloudscreenedbt_15day_array),
                                ianlats.shape[
                                    0], ianlons.shape[0]))

        for k in np.arange(0, len(cloudscreenedbt_15day_array)):
            regcs = pinkdust.regrid_data(lons, lats, ianlons, ianlats,
                                         cloudscreenedbt_15day_array[k],
                                                        mesh=False)
            regcs_array[k] = regcs

        regbt087 = pinkdust.regrid_data(lons, lats, ianlons, ianlats,
                                        bt087, mesh=False)
        regbt108 = pinkdust.regrid_data(lons, lats, ianlons, ianlats,
                                        bt108, mesh=False)
        regbt120 = pinkdust.regrid_data(lons, lats, ianlons, ianlats,
                                        bt120, mesh=False)

        cs[j] = regcs_array
        bts[j, 0] = regbt087
        bts[j, 1] = regbt108
        bts[j, 2] = regbt120
        btdata.close()

    # Save cloud masked data for this time of day to file

    f.close()

    # Save cloud masked data for this time of day to file

    g.close()

if __name__ == '__main__':

    # Select the time parameters within which to parse input nc files and
    # obtain SDF values

    # This needs to be a period for which there is a 15 day time window
    # around it
    year_lower = 2011
    year_upper = 2011
    month_lower = 6
    month_upper = 8
    day_lower = 1
    day_upper = 31
    hour_lower = 0
    hour_upper = 23
    minute_lower = 0
    minute_upper = 45
    pickup = True
    pickupdate = datetime.datetime(2011, 8, 6, 0, 0)

    # Generate the full set of datetimes for the entire period
    time_params = np.array([year_lower, year_upper, month_lower,
                            month_upper, day_lower, day_upper,
                            hour_lower, hour_upper, minute_lower,
                            minute_upper])

    SDF_datetimes = utilities.get_datetime_objects(time_params)
    SDF_datestrings = [j.strftime("%Y%m%d%H%M") for j in SDF_datetimes]

    # Get Ian's lats and lons
    sdf_test = Dataset(
        '/soge-home/data_not_backed_up/satellite/meteosat/seviri'
        '/15-min/0.03x0.03/sdf/nc/JUNE2010/SDF_v2/SDF_v2.'
        '201006031500.nc')

    ianlons = sdf_test.variables['longitude'][:]
    ianlats = sdf_test.variables['latitude'][:]

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

    # Get Ian's lats and lons
    sdf_test = Dataset(
        '/soge-home/data_not_backed_up/satellite/meteosat/seviri'
        '/15-min/0.03x0.03/sdf/nc/JUNE2010/SDF_v2/SDF_v2.'
        '201006031500.nc')

    ianlons = sdf_test.variables['longitude'][:]
    ianlats = sdf_test.variables['latitude'][:]

    """

    oneday_datetimes = utilities.get_datetime_objects(time_params_oneday)

    print 'Generating cloud masks for all timesteps'
    # Loop through each time in a single day

    # Loop through blocks of 15 - what happens if we don't have 15 separate
    # dates to do?
    pool = multiprocessing.Pool()
    for i in np.arange(0, len(oneday_datetimes)):
        pool.apply_async(cloud_mask_mw, args=(i,))
    pool.close()
    pool.join()

    gmail_user = 'tcharrison4@gmail.com'
    gmail_password = None

    sent_from = gmail_user
    to = ['tcharrison4@gmail.com', 'thomas.catonharrison@keble.ox.ac.uk']
    subject = 'SDF moving window:  complete'
    body = 'Dear Thomas,\n\nMoving window processing of selected ' \
           'SDFs had been completed.\n\nLots ' \
           'of love, ' \
           '\n\nPython'

    email_text = 'Subject: {}\n\n{}'.format(subject, body)

    s = smtplib.SMTP("smtp.gmail.com", 587)
    s.ehlo()
    s.starttls()
    s.login(gmail_user, gmail_password)
    s.sendmail(sent_from, to, email_text)
    s.close()

    """
    
    ### SDF PHASE ###

    # Open an nc file in which to store all SDF data for the run period

    """

    print 'Generating SDF nc file'
    # File description and dimensions
    SDFfile = Dataset('/soge-home/projects/seviri_dust/sdf/SDF_'
                      ''+SDF_datetimes[0].strftime(
        '%Y_%m_%d')+'_to_'+SDF_datetimes[-1].strftime('%Y_%m_%d')+'.nc', 'w',
        format='NETCDF4_CLASSIC')
    SDFfile.description = 'SDF values from EUMETSAT SEVIRI on a ' \
                          'regular lat/lon grid.'
    SDFfile.createDimension('time', len(SDF_datetimes))
    SDFfile.createDimension('lat', lats.shape[0])
    SDFfile.createDimension('lon', lons.shape[1])

    # Variables
    times = SDFfile.createVariable('time', np.float64, ('time'))
    SDFvar = SDFfile.createVariable('SDF', np.int32, ('time', 'lat', 'lon'))

    # Attributes
    times.units = 'hours since 0001-01-01 00:00:00'
    times.calendar = 'gregorian'

    # Writing data to file and closing
    times[:] = date2num(SDF_datetimes, units=times.units,
                        calendar=times.calendar)

    """

    print 'Generating SDFs'
    if pickup:
        date_indices = np.arange(0, len(SDF_datetimes))
        pickup_index = SDF_datetimes == pickupdate
        pickup_index = date_indices[pickup_index][0]
    else:
        pickup_index = 0

    # Loop through the full set of datetimes for the entire period
    for i in np.arange(pickup_index, len(SDF_datetimes)):
        print SDF_datetimes[i]

        # Read in the cloud mask data for that time of day from file
        # Remove the file after use as it is no longer needed
        #print 'here 1'
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

        f = tables.open_file(
            '/soge-home/projects/seviri_dust/sdf/intermediary_files/'
            'cloud_masked_bt_15d_' + SDF_datetimes[i].strftime(
                "%Y_%H_%M") + '.hdf')
        #arrobj = f.get_node('/data')

        # Take out indices since pytables doesn't seem to support fancy
        # indexing with booleans...
        indices = np.arange(0, len(datetimes_7dayw))
        lower_ind = datetimes_7dayw == BT_15_day_lower_bound
        lower_ind = indices[lower_ind][0]
        upper_ind = datetimes_7dayw == BT_15_day_upper_bound
        upper_ind = indices[upper_ind][0]
        current_ind = datetimes_7dayw == SDF_datetimes[i]
        current_ind = indices[current_ind][0]

        BT_15_days_screened = f.root.data[lower_ind:upper_ind]
        BT_screened = f.root.data[current_ind]
        #cloudscreenedbt_15day = arrobj.read()
        # print 'here 3'
        f.close()

        f = tables.open_file(
            '/soge-home/projects/seviri_dust/sdf/intermediary_files/bt_15d_' + SDF_datetimes[i].strftime(
                "%Y_%H_%M") + '.hdf')
        arrobj = f.get_node('/data')
        bt_15day = f.root.data[current_ind]
        #bt_15day = arrobj.read()
        f.close()

        #print 'here 4'

        #BT_15_days_screened = \
         #   cloudscreenedbt_15day[np.asarray([j >= BT_15_day_lower_bound
          #                                    and j <= BT_15_day_upper_bound
           #                                   for j in datetimes_7dayw])]
        #print 'here 5'
        # Get the mean of the 15 day window of cloud screened data
        cloudscreenedbt_15daymean_108_87 = \
            sdf.extract_15day_mean(BT_15_days_screened)

        dust_mask_108, \
        dust_mask_108_87, \
        dust_mask_120_108, \
        dust_mask_108_87_anom_screened = \
            sdf.generate_dust_mask(
                bt_15day,
                BT_screened,
                cloudscreenedbt_15daymean_108_87)

        ofilename = '/soge-home/projects/seviri_dust/sdf/SDF_v2.' + \
                    SDF_datetimes[i].strftime(
                "%Y%m%d%H%M") + '.nc'

        SDF = sdf.generate_SDF(dust_mask_108,
                               dust_mask_108_87,
                               dust_mask_120_108,
                               dust_mask_108_87_anom_screened,
                               ofilename, SDF_datetimes[i], ianlats, ianlons)

        #print 'here 6'
        # THIS BIT IS A MAJOR TIME SINK
        #SDF_data = SDFvar[:]
        #print 'here 7'
        #time_bool = np.asarray([j == SDF_datetimes[i] for j in SDF_datetimes])
        #print 'here 8'
        #SDF_data[time_bool] = SDF
        #print 'here 9'
        #SDFvar[:] = SDF_data
        #print 'here 10'
    #SDFfile.close()

