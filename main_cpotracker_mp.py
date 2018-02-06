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
from pyresample import utils

import utilities
import plumes
import cpo_detection
import cpo
import pinkdust
import plotting
import get_llj_prob_model

def wrapper(yearmonth):

    year_lower = yearmonth[0]
    year_upper = yearmonth[0]
    month_lower = yearmonth[-1][0]
    #month_lower = 7
    month_upper = yearmonth[-1][-1]
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
    datestrings = [j.strftime("%Y%m%d%H%M") for j in datetimes]

    # Get lats and lons
    sdf_test = Dataset(
        '/soge-home/data_not_backed_up/satellite/meteosat/seviri'
        '/15-min/0.03x0.03/sdf/nc/JUNE2010/SDF_v2/SDF_v2.'
        '201006031500.nc')

    lons, lats = np.meshgrid(sdf_test.variables['longitude'][:],
                             sdf_test.variables['latitude'][:])
    lonmask = lons > 360
    latmask = lats > 90
    lons = np.ma.array(lons, mask=lonmask)
    lats = np.ma.array(lats, mask=latmask)
    cpo_previous = None
    raw_cpo_prev = []
    clouds_previous = None
    ids_previous = []
    cloud_ids_previous = []
    deep_conv_IDs_prev = None
    LLJ_plumes_IDs_prev = []
    k = 0
    available_colours = np.arange(0, 41)
    # To pick back up where you left off, simply add all the keys in plume
    # archive to used IDs here - perhaps include an option to do that
    used_ids = []
    used_cloud_ids = []
    used_colour_IDs = {}
    flicker_ids = []
    reintroduced_ids = []
    last_10_ids = []
    last_10_ids = np.asarray(last_10_ids)
    prev_dust_assoc_clouds = None
    found_file = True

    # Restrict the lons and lats to the CWS alone
    lonbool = np.asarray([j >= -20 and j <= 10 for j in lons[0]])


    cpo_objects = {}

    cpo_archive = shelve.open(
        '/soge-home/projects/seviri_dust/plumetracker/'
                                'cpo_archive_v4'
        ''+str(yearmonth[0]))

    if pickup:
        archived_ids = np.asarray([j for j in cpo_archive])
        for i in archived_ids:
            used_ids.append(int(i))
        with open('date_i_'+str(yearmonth[0])+'.txt', 'r') as f:
            pickup_date_i = f.read()
        print 'Picking up from timestep', pickup_date_i
        datestrings = datestrings[int(pickup_date_i):]
        datetimes = datetimes[int(pickup_date_i):]

    cpo_archive.close()

    year_lower = datetimes[0].year
    year_upper = datetimes[-1].year
    month_lower = datetimes[0].month
    month_upper = datetimes[-1].month
    day_lower = datetimes[0].day
    day_upper = datetimes[-1].day
    hour_lower = 0
    hour_upper = 23
    minute_lower = 0
    minute_upper = 45

    time_params_oneday = np.array([year_lower, year_lower, month_lower,
                                   month_lower, day_lower, day_lower,
                                   hour_lower, hour_upper,
                                   minute_lower, minute_upper])

    oneday_datetimes = utilities.get_datetime_objects(time_params_oneday)

    if moving_window:
        # Here run the moving window for our selected time period
        print 'Generating cloud masks for all timesteps'
        # Loop through each time in a single day

        pool = multiprocessing.Pool()
        for i in np.arange(0, len(oneday_datetimes)):
            pool.apply_async(cpo_detection.cloud_mask_mw, args=(i, datetimes,
                                                 oneday_datetimes,
                                            ianlons, ianlats))
            # Pull a netCDF dataset object out so the projection coordinates can
            # be obtained
        pool.close()
        pool.join()

    btdiff_2_anom_prev = None
    btdiff_2_anom_prev_2 = None
    btdiff_2_anom_prev_3 = None

    for date_i in np.arange(0, len(datestrings)):
        cpo_archive = shelve.open(
            '/soge-home/projects/seviri_dust/plumetracker/'
            'cpo_archive_v4_'
            '' + str(yearmonth[0]))
        runtime = datetimes[date_i] - datetimes[0]
        print '\n' + datestrings[date_i] + '\n'
        totaltest = datetime.datetime.now()

        if found_file:

            #if 'time' in bt.variables:
            #    bt_108 = bt.variables['bt108'][:][0]
            #else:
            #    bt_108 = bt.variables['bt108'][:]
            #clouds = bt_108 < 270

            # Here run CPO detection
            # IF A RAW DETECTED CPO FILE EXISTS, USE THAT (assuming the
            # algorithm hasn't changed)
            if calculate_cpof:
                cpo_now, btdiff_2_anom_prev, btdiff_2_anom_prev_2, \
                btdiff_2_anom_prev_3 = cpo_detection.detect_cpo(
                    btdiff_2_anom_prev, btdiff_2_anom_prev_2,
                           btdiff_2_anom_prev_3, datetimes, datestrings, date_i,
                    lons, lats,
                    cloud_lons, cloud_lats,
                    daily_cloudmask, double_digits, mesh, daily_bt)
                cpo_now[cpo_now > 0] = 1
            else:
                cpo_data = Dataset(
                    '/ouce-home/projects/seviri_dust/cpof/'+datetimes[
                        date_i].strftime("%B%Y")+'/CPOF_'+datetimes[
                        date_i].strftime("%Y%m%d%H%M")+'.nc')
                cpo_now = cpo_data.variables['CPOF'][:]

            # Save the raw CPO detection to file
            # These should be stored in the project directory, not home
            #np.save('raw_detected_cpo_'+datestrings[date_i], cpo_now)

            # Add the reintroduced plumes back into the running
            for i in np.arange(0, len(reintroduced_ids)):
                #print 'Reintroducing CPO', reintroduced_ids[i]
                cpo_previous[cpo_previous == flicker_ids[i]] = \
                    reintroduced_ids[i]
                ids_previous = np.append(ids_previous, reintroduced_ids[i])

            # Get plumes first by scanning for them
            cpo_blobs, new_ids, blob_ids, merge_ids = cpo.\
                scan_for_blobs(
                cpo_now,
                cpo_previous,
                used_ids)

            #print 'New IDs', new_ids

            for i in np.arange(0, len(reintroduced_ids)):
                if reintroduced_ids[i] not in blob_ids:
                    #print 'But plume '+str(reintroduced_ids[i])+' sadly ' \
                    #                                            'died.'
                    pass

            # Clear flicker IDs and reintroduced IDs
            flicker_ids = []
            reintroduced_ids = []

            for i in new_ids:
                used_ids.append(i)

            old_bool = np.asarray([j in ids_previous for j in blob_ids])
            if len(old_bool) > 0:
                old_ids = blob_ids[old_bool]
            else:
                old_ids = []

            # Then, for each new ID, we initialise plume objects
            for i in np.arange(0, len(new_ids)):
                if debug:
                    print 'Creating new CPO', new_ids[i]
                cold_pool = cpo.CPO(new_ids[i], datetimes[date_i])
                cold_pool.update_position(lats, lons, cpo_blobs, new_ids[i])
                cold_pool.update_duration(datetimes[date_i])
                #if plume.plume_id == 18:
                #    print plume.dates_observed
                #    print plume.dates_observed
                cold_pool.update_bbox()
                cold_pool.update_majorminor_axes(lons, lats)
                cold_pool.update_area()
                cold_pool.update_max_extent()
                cold_pool.update_centroid_speed()
                cold_pool.update_centroid_direction()
                #cold_pool.check_conv_distance(lats, lons, clouds)
                cold_pool.update_most_likely_source()
                # plume.update_leading_edge_4(sdf_plumes, lons, lats)
                cpo_objects[str(new_ids[i])] = cold_pool
                missing_cold_pool, missing_id, flickered = \
                    cold_pool.chain_flickerchecker(
                    raw_cpo_prev)
                if flickered:
                    raise ValueError('Found an overlapping CPO in the '
                                     'previous timestep larger than size '
                                     '50 - a new CPO should not be '
                                     'initiated here')
                steps_back = 1
                # As long as there is an overlapping previous plume,
                #  keep updating it back in time
                while len(missing_cold_pool) > 0:
                    if debug:
                        print 'Rolling back CPO', new_ids[i]
                    # We can only step back to the first timestep and no
                    # earlier
                    if (date_i - steps_back) < 0:
                        missing_cold_pool = []
                        break
                    else:
                        missing_date = datetimes[date_i-steps_back]
                        missing_cpo_blobs = np.zeros(missing_cold_pool.shape)
                        missing_cold_pool = missing_cold_pool == 1
                        missing_cpo_blobs[missing_cold_pool] = missing_id
                        # Run all the updates that would be used for a
                        # new plume
                        cold_pool.update_position(lats, lons,
                                              missing_cpo_blobs,
                                              missing_id)
                        cold_pool.update_duration(missing_date)
                        #if plume.plume_id == 18:
                        #    print plume.dates_observed
                        cold_pool.update_bbox()
                        cold_pool.update_majorminor_axes(lons, lats)
                        cold_pool.update_area()
                        cold_pool.update_max_extent()
                        cold_pool.update_centroid_speed()
                        cold_pool.update_centroid_direction()
                        #cold_pool.check_conv_distance(lats, lons, clouds)
                        cold_pool.update_most_likely_source()
                        cold_pool.process_missing_plume()
                        #print 'Updated missing plume back '+str(
                        #    steps_back)+' steps'
                        steps_back += 1
                        if (date_i - steps_back) < 0:
                            missing_cpo = []
                            break
                        # Pull out data from the timestep before to
                        # continue the chain
                        try:
                            raw_cpo_prev_prev = \
                                np.load('raw_detected_cpo_'+datestrings[
                                date_i]+'.npy')
                        except:
                            #print 'Adding date to list of missing dates'
                            with open('missing_dates.txt', 'a') as my_file:
                                my_file.write('\n'+datestrings[date_i-
                                                          steps_back])
                            break

                        missing_cold_pool, missing_id, flickered = \
                            cold_pool.chain_flickerchecker(
                        raw_cpo_prev_prev)

                        if flickered:
                            #print 'Found a flickered plume. Searching ' \
                            #      'for the corresponding archived plume.'
                            # We have a plume in a previous timestep
                            # which flickered
                            cold_pool_archive_keys = last_10_ids.astype(int)
                            # Sort the keys in reverse order as the
                            # plume we want is most likely to have a
                            # high ID
                            cold_pool_archive_keys[::-1].sort()

                            missing_cold_pool = missing_cold_pool == 1
                            missing_cpo_blobs = np.zeros(
                                missing_cold_pool.shape)
                            missing_cpo_blobs[missing_cold_pool] = missing_id
                            cold_pool_bool = missing_cpo_blobs == missing_id

                            search_centroid_lon = \
                                np.nanmean(lons[cold_pool_bool])

                            search_centroid_lat = \
                                np.nanmean(lats[cold_pool_bool])

                            cold_pool_lons = lons[cold_pool_bool]
                            cold_pool_lats = lats[cold_pool_bool]

                            search_date = datetimes[date_i-steps_back]
                            #print search_date

                            found_cold_pool = False
                            for key in cold_pool_archive_keys:
                                #print 'Searching in plume archive'
                                if search_centroid_lon == \
                                        cpo_archive[\
                                        str(key)].centroid_lon and \
                                        search_centroid_lat ==\
                                        cpo_archive[str(
                                            key)].centroid_lat and \
                                        cpo_archive[str(
                                            key)].dates_observed[-1] == \
                                                search_date:
                                    #print 'Found it in plume archive. ' \
                                    #      'ID is', \
                                    #    plume_archive[str(key)].plume_id
                                    found_cold_pool = True

                                    correct_cold_pool = cpo_archive[str(key)]
                                    cold_pool_to_append = cold_pool
                                    # Append the flickered plume to the
                                    # old one which was archived
                                    correct_cold_pool.append_missing_cold_pool(
                                        cold_pool_to_append)

                                    # Add it to plume objects and remove
                                    # it from archives
                                    cpo_objects[str(key)] = correct_cold_pool
                                    del cpo_archive[str(key)]

                                    # Add it to old IDs, replacing the
                                    # ID of the plume which was found to be
                                    # flickered
                                    flicker_ids.append(cold_pool.plume_id)
                                    reintroduced_ids.append(key)
                                    missing_cold_pool = []

                                    # Reintroduced plumes also get removed
                                    # from the record of the last 10 ids
                                    index = np.argwhere(last_10_ids
                                                        == key)
                                    last_10_ids = np.delete(last_10_ids,
                                                            index)

                                    break
                            # If we didn't find the plume in the plume
                            # archive, it must still be active
                            if found_cold_pool == False:
                                cold_pool_object_keys = np.asarray([int(i) for
                                                                 i in
                                                            cpo_objects])
                                # Sort the keys in reverse order as the
                                # plume we want is most likely to have a
                                # high ID
                                cold_pool_object_keys[::-1].sort()
                                for key in cold_pool_object_keys:
                                    #print 'Searching in plume objects'

                                    if search_centroid_lon == \
                                            cpo_objects[ \
                                                    str(
                                                key)].centroid_lon and \
                                                search_centroid_lat == \
                                                cpo_objects[str(
                                                key)].centroid_lat and \
                                                cpo_objects[str(
                                                key)].dates_observed[
                                                -1] == \
                                                search_date:

                                        found_cold_pool = True

                                        correct_cold_pool = cpo_archive[
                                            str(key)]
                                        cold_pool_to_append = cold_pool
                                        # Append the flickered plume to the
                                        # old one which was archived
                                        correct_cold_pool.append_missing_cold_pool(
                                            cold_pool_to_append)

                                        # Add it to plume objects and
                                        # remove it from archives
                                        cpo_objects[
                                            str(key)] = correct_cold_pool
                                        del cpo_archive[str(key)]

                                        # Add it to old IDs, replacing the
                                        # ID of the plume which was found
                                        # to be flickered
                                        flicker_ids.append(cold_pool.plume_id)
                                        reintroduced_ids.append(key)
                                        missing_cold_pool = []

                                        index = np.argwhere(last_10_ids
                                                            == key)
                                        last_10_ids = np.delete(
                                            last_10_ids,
                                            index)

                                        break

                            break

            # Remove any IDs which were actually flickers
            for i in np.arange(0, len(flicker_ids)):
                index = np.argwhere(new_ids==flicker_ids[i])
                new_ids = np.delete(new_ids, index)
                index = np.argwhere(ids_previous == flicker_ids[i])
                ids_previous = np.delete(ids_previous, index)
                index = np.argwhere(blob_ids == flicker_ids[i])
                blob_ids = np.delete(blob_ids, index)
                del cpo_objects[str(flicker_ids[i])]

            #For merged IDs, we move the tracks to pre-merge tracks
            for i in np.arange(0, len(merge_ids)):
                cold_pool = cpo_objects[str(merge_ids[i])]
                cold_pool.merge()

            # For old IDs, we just run an update.
            for i in np.arange(0, len(old_ids)):
                if debug:
                    print 'Updating CPO', old_ids[i]
                cold_pool = cpo_objects[str(old_ids[i])]
                #if plume.plume_id == 2:
                #    print plume.dates_observed
                cold_pool.update_position(lats, lons, cpo_blobs, old_ids[i])
                cold_pool.update_duration(datetimes[date_i])
                #if plume.plume_id == 18:
                #    print plume.dates_observed
                cold_pool.update_bbox()
                cold_pool.update_majorminor_axes(lons, lats)
                cold_pool.update_area()
                cold_pool.update_centroid_speed()
                cold_pool.update_centroid_direction()
                cold_pool.update_max_extent()
                #plume_assoc_clouds, plume_check_bool \
                #    = plume.flag_dust_associated_convection(
                #    clouds, prev_dust_assoc_clouds)
                #dust_assoc_clouds += plume_assoc_clouds
                #check_bool += plume_check_bool.astype(int)
                # plume.update_leading_edge_4(sdf_plumes, lons, lats)
                cpo_objects[str(old_ids[i])] = cold_pool

            # Plumes which no longer exist are removed and archived
            if len(ids_previous) == 0:
                removed_ids = []
            else:
                removed_bool = np.asarray([j not in blob_ids for j in
                                           ids_previous])
                removed_ids = ids_previous[removed_bool]

            for i in np.arange(0, len(removed_ids)):
                if debug:
                    print 'Archiving CPO', removed_ids[i]
                cold_pool = cpo_objects[str(removed_ids[i])]
                cold_pool.update_GPE_speed()
                # plume.update_mechanism_likelihood()
                cold_pool.update_mean_axis_offset()
                if cold_pool.duration.total_seconds() >= 3600:
                    cpo_archive[str(removed_ids[i])] = cold_pool
                    del cpo_objects[str(removed_ids[i])]
                    print '== CPO detected =='
                    print 'CPO emission time:', cold_pool.emission_time
                    print 'CPO emission lat:', cold_pool.track_centroid_lat[0]
                    print 'CPO emission lon:', cold_pool.track_centroid_lon[0]
                    last_10_ids = np.append(last_10_ids, removed_ids[i])
                else:
                    #print 'Too short duration - not archived'
                    del cpo_objects[str(removed_ids[i])]
                #print 'Plume', plume.plume_id,'removed. Final lons and \
                #                                             lats:'
                #print plume.centroid_lon
                #print plume.centroid_lat

            if len(last_10_ids) > 10:
                last_10_ids = last_10_ids[0:10]

            if len(np.unique(cpo_blobs)) < 2:
                cpo_previous = None
                ids_previous = []
                raw_cpo_prev = []
            else:
                cpo_previous = cpo_blobs
                ids_previous = blob_ids
                raw_cpo_prev = cpo_now

        else:
            print 'Adding date to list of missing dates'
            with open('missing_dates.txt', 'a') as my_file:
                my_file.write('\n'+datestrings[date_i])

        with open('date_i_'+str(yearmonth[0])+'.txt', 'w') as f:
            f.write('%d' % date_i)
        cpo_archive.close()

    # There seems to be a residual massive file - remove this

if __name__ == '__main__':

    pickup = False
    moving_window = False
    calculate_cpof = False
    daily_bt = False
    daily_cloudmask = False
    bt_root = ''
    daily_cloudmask_root = '/soge-home/projects/seviri_dust/raw_seviri_data' \
                           '/cloudmask_nc/'
    debug = True

    if daily_cloudmask:
        double_digits = True # Set to True if using daily cloudmask
        mesh = False# Set this to False if using daily cloudmask
    else:
        double_digits = False
        mesh = True

    # Get Ian's lats and lons
    sdf_test = Dataset(
        '/soge-home/data_not_backed_up/satellite/meteosat/seviri'
        '/15-min/0.03x0.03/sdf/nc/JUNE2010/SDF_v2/SDF_v2.'
        '201006031500.nc')

    ianlons = sdf_test.variables['longitude'][:]
    ianlats = sdf_test.variables['latitude'][:]

    target_area = utils.load_area(
            '/soge-home/projects/seviri_dust/areas.def',
            'NorthAfrica')
    cloud_lons, cloud_lats = target_area.get_lonlats()

    # Get lats and lons
    sdf_test = Dataset(
        '/soge-home/data_not_backed_up/satellite/meteosat/seviri'
        '/15-min/0.03x0.03/sdf/nc/JUNE2010/SDF_v2/SDF_v2.'
        '201006031500.nc')

    lons, lats = np.meshgrid(sdf_test.variables['longitude'][:],
                             sdf_test.variables['latitude'][:])

    years = [2008, 2009, 2010]
    months = [[6, 7, 8], [6, 7, 8], [6, 7, 8]]
    years = [2009]
    months = [[6, 7, 8]]
    #months = [[7]]
    #months = [[7], [7], [7], [7], [7], [7], [7]]
    yearmonths = zip(years, months)

    wrapper(yearmonths[0])
    print 1/0

    pool = multiprocessing.Pool()
    for i in yearmonths:
        pool.apply_async(wrapper, args=(i,))
        # Pull a netCDF dataset object out so the projection coordinates can
        # be obtained
    pool.close()
    pool.join()

    print 1/0

    # 2007, 2010, 2009, 2006, 2008

