import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import numpy as np
from mpl_toolkits.basemap import Basemap
import datetime
import shelve
import multiprocessing
import os
import copy
import os.path

import utilities
import plumes
import pinkdust
import plotting
import get_llj_prob_model

def wrapper(yearmonth):

    year_lower = yearmonth[0]
    year_upper = yearmonth[0]
    month_lower = yearmonth[1][0]
    month_upper = yearmonth[1][-1]
    month_upper = 6
    day_lower = 1
    day_upper = 1
    hour_lower = 0
    hour_upper = 10
    minute_lower = 0
    minute_upper = 0

    time_params = np.array([year_lower, year_upper, month_lower,
                            month_upper, day_lower, day_upper,
                            hour_lower, hour_upper, minute_lower,
                            minute_upper])

    datetimes = utilities.get_datetime_objects(time_params)
    datestrings = [j.strftime("%Y%m%d%H%M") for j in datetimes]

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
    sdf_previous = None
    raw_sdf_prev = []
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
    plume_objects = []
    flicker_ids = []
    reintroduced_ids = []
    last_50_ids = []
    last_50_ids = np.asarray(last_50_ids)

    # Restrict the lons and lats to the CWS alone
    lonbool = np.asarray([j >= -20 and j <= 10 for j in lons[0]])

    if run:

        plume_objects = {}

        plume_archive = shelve.open(
            '/soge-home/projects/seviri_dust/plumetracker/'
                                    'plume_archive_flicker_v2_debug_'+str(
                yearmonth[
                                                                        0]))

        if pickup:
            archived_ids = np.asarray([j for j in plume_archive])
            for i in archived_ids:
                used_ids.append(int(i))

        for date_i in np.arange(0, len(datestrings)):
            runtime = datetimes[date_i] - datetimes[0]
            print '\n' + datestrings[date_i] + '\n'
            totaltest = datetime.datetime.now()

            if os.path.isfile('/soge-home/data_not_backed_up/satellite/'
                              'meteosat/seviri/15-min/0.03x0.03/sdf/nc/' +
                datetimes[date_i].strftime("%B").upper(
                ) + str(datetimes[date_i].year) + '/SDF_v2/SDF_v2.' + \
                datestrings[date_i] + '.nc'):

                sdf = Dataset(
                    '/soge-home/data_not_backed_up/satellite/meteosat/seviri/'
                    '15-min/0.03x0.03/sdf/nc/' +
                    datetimes[date_i].strftime("%B").upper(
                    ) + str(datetimes[date_i].year) + '/SDF_v2/SDF_v2.' + \
                    datestrings[date_i] + '.nc')
                found_file = True
            else:
                'No file found for this date'
                found_file = False

            try:
                bt = Dataset(
                    '/ouce-home/data/satellite/meteosat/seviri/15-min/'
                    '0.03x0.03/bt'
                    '/nc/'
                    +
                    datetimes[date_i].strftime("%B").upper(
                    ) + str(datetimes[date_i].year) + '/H-000-MSG2__'
                                                      '-MSG2________-'
                                        'IR_BrightnessTemperatures___'
                                                      '-000005___-'
                    + datestrings[date_i] +
                    '-__.nc')
            except:

                if os.path.isfile('/ouce-home/data/satellite/meteosat/seviri/'
                                  '15-min/'
                    '0.03x0.03/bt'
                    '/nc/'
                    +
                    datetimes[date_i].strftime("%B").upper(
                    ) + str(datetimes[date_i].year) + '/H-000-MSG1__'
                                                      '-MSG1________-'
                                        'IR_BrightnessTemperatures___'
                                                      '-000005___-'
                    + datestrings[date_i] +
                    '-__.nc'):
                    bt = Dataset(
                        '/ouce-home/data/satellite/meteosat/seviri/15-min/'
                        '0.03x0.03/bt'
                        '/nc/'
                        +
                        datetimes[date_i].strftime("%B").upper(
                        ) + str(datetimes[date_i].year) + '/H-000-MSG1__'
                                                          '-MSG1________-'
                                            'IR_BrightnessTemperatures___'
                                                          '-000005___-'
                        + datestrings[date_i] +
                        '-__.nc')
                else:
                    found_file = False

            if found_file:

                # If the SDF file has its own longitudes and latitudes,
                # use these rather than the preset
                if 'longitude' in sdf.variables:
                    lons = sdf.variables['longitude'][:]
                    lonmask = lons > 360
                    lons = np.ma.array(lons, mask=lonmask)
                if 'latitude' in sdf.variables:
                    lats = sdf.variables['latitude'][:]
                    latmask = lats > 90
                    lats = np.ma.array(lats, mask=latmask)
                    lons, lats = np.meshgrid(lons, lats)

                # clouds_now = np.load('/ouce-home/students/hert4173/'
                #                        'cloud_mask_numpy_files/cloudmask'
                #    ''+datetimes[date_i]
                #    .strftime(
                #    "%Y%m%d%H%M%S")+'.npy')

                # Some SDF files have a time dimension. For these index it out.
                if 'time' in sdf.variables:
                    sdf_now = sdf.variables['bt108'][0]
                else:
                    sdf_now = sdf.variables['bt108'][:]

                if 'time' in bt.variables:
                    bt_108 = bt.variables['bt108'][:][0]
                else:
                    bt_108 = bt.variables['bt108'][:]
                clouds = bt_108 < 270

                # Add the reintroduced plumes back into the running
                for i in np.arange(0, len(reintroduced_ids)):
                    #print 'Reintroducing plume', reintroduced_ids[i]
                    sdf_previous[sdf_previous == flicker_ids[i]] = \
                        reintroduced_ids[i]
                    ids_previous = np.append(ids_previous, reintroduced_ids[i])

                # Get plumes first by scanning for them
                sdf_plumes, new_ids, plume_ids, merge_ids = plumes.\
                    scan_for_plumes(
                    sdf_now,
                    sdf_previous,
                    used_ids,
                    clouds)

                print plume_ids

                for i in np.arange(0, len(reintroduced_ids)):
                    if reintroduced_ids[i] not in plume_ids:
                        #print 'But plume '+str(reintroduced_ids[i])+' sadly ' \
                        #                                            'died.'
                        pass

                # Clear flicker IDs and reintroduced IDs
                flicker_ids = []
                reintroduced_ids = []

                # We could here do infilling, then at least you'll have it for
                # the next iteration. But if you've labelled them already you
                # don't have a binary.

                # Get clouds by scanning for them
                # clouds, new_cloud_ids, cloud_ids, merge_cloud_ids = \
                #    plumes.scan_for_plumes(clouds_now,
                #                           clouds_previous,
                #                           used_cloud_ids)

                for i in new_ids:
                    used_ids.append(i)

                # for i in new_cloud_ids:
                #    used_cloud_ids.append(i)

                old_bool = np.asarray([j in ids_previous for j in plume_ids])
                if len(old_bool) > 0:
                    old_ids = plume_ids[old_bool]
                else:
                    old_ids = []

                # old_cloud_bool = np.asarray([j in clouds_previous for
                # j in cloud_ids])
                # old_cloud_ids = cloud_ids[old_cloud_bool]

                # Then, for each new ID, we initialise plume objects
                for i in np.arange(0, len(new_ids)):
                    print 'Creating new plume', new_ids[i]
                    plume = plumes.Plume(new_ids[i], datetimes[date_i])
                    plume.update_position(lats, lons, sdf_plumes, new_ids[i])
                    plume.update_duration(datetimes[date_i])
                    if plume.plume_id == 18:
                        print plume.dates_observed
                    #    print plume.dates_observed
                    plume.update_bbox()
                    plume.update_majorminor_axes()
                    plume.update_area()
                    plume.update_max_extent()
                    plume.update_centroid_speed()
                    plume.update_centroid_direction()
                    plume.check_conv_distance(lats, lons, clouds)
                    plume.update_most_likely_source()
                    # plume.update_leading_edge_4(sdf_plumes, lons, lats)
                    plume_objects[str(new_ids[i])] = plume
                    missing_plume, missing_id, flickered = \
                        plume.chain_flickerchecker(
                        raw_sdf_prev)
                    if flickered:
                        raise ValueError('Found an overlaping plume in the '
                                         'previous timestep larger than size '
                                         '250 - a new plume should not be '
                                         'initiated here')
                    steps_back = 1
                    # As long as there is an overlapping previous plume,
                    #  keep updating it back in time
                    while len(missing_plume) > 0:
                        print 'Rolling back plume', new_ids[i]
                        # We can only step back to the first timestep and no
                        # earlier
                        if (date_i - steps_back) < 0:
                            missing_plume = []
                            break
                        else:
                            missing_date = datetimes[date_i-steps_back]
                            missing_sdf_plumes = np.zeros(missing_plume.shape)
                            missing_plume = missing_plume == 1
                            missing_sdf_plumes[missing_plume] = missing_id
                            # Run all the updates that would be used for a
                            # new plume
                            plume.update_position(lats, lons,
                                                  missing_sdf_plumes,
                                                  missing_id)
                            plume.update_duration(missing_date)
                            if plume.plume_id == 18:
                                print plume.dates_observed
                            plume.update_bbox()
                            plume.update_majorminor_axes()
                            plume.update_area()
                            plume.update_max_extent()
                            plume.update_centroid_speed()
                            plume.update_centroid_direction()
                            plume.check_conv_distance(lats, lons, clouds)
                            plume.update_most_likely_source()
                            plume.process_missing_plume()
                            #print 'Updated missing plume back '+str(
                            #    steps_back)+' steps'
                            steps_back += 1
                            if (date_i - steps_back) < 0:
                                missing_plume = []
                                break
                            # Pull out data from the timestep before to
                            # continue the chain
                            try:
                                raw_sdf_prev_prev_data = Dataset(
                                '/soge-home/data_not_backed_up/satellite/'
                                'meteosat/'
                                'seviri/'
                                '15-min/0.03x0.03/sdf/nc/' +
                                datetimes[date_i-steps_back].strftime(
                                    "%B").upper(
                                ) + str(datetimes[date_i-steps_back].year) +
                                '/SDF_v2/SDF_v2.' + \
                                datestrings[date_i-steps_back] + '.nc')
                            except:
                                print 'Adding date to list of missing dates'
                                with open('missing_dates.txt', 'a') as my_file:
                                    my_file.write('\n'+datestrings[date_i-
                                                              steps_back])
                                break
                            if 'time' in raw_sdf_prev_prev_data.variables:
                                raw_sdf_prev_prev = \
                                    raw_sdf_prev_prev_data.variables[
                                        'bt108'][0]
                            else:
                                raw_sdf_prev_prev = \
                                    raw_sdf_prev_prev_data.variables[
                                'bt108'][:]
                            missing_plume, missing_id, flickered = \
                                plume.chain_flickerchecker(
                            raw_sdf_prev_prev)

                            if flickered:
                                #print 'Found a flickered plume. Searching ' \
                                #      'for the corresponding archived plume.'
                                # We have a plume in a previous timestep
                                # which flickered
                                plume_archive_keys = last_50_ids.astype(int)
                                # Sort the keys in reverse order as the
                                # plume we want is most likely to have a
                                # high ID
                                plume_archive_keys[::-1].sort()

                                missing_plume = missing_plume == 1
                                missing_sdf_plumes = np.zeros(
                                    missing_plume.shape)
                                missing_sdf_plumes[missing_plume] = missing_id
                                plume_bool = missing_sdf_plumes == missing_id

                                search_centroid_lon = \
                                    np.nanmean(lons[plume_bool])

                                search_centroid_lat = \
                                    np.nanmean(lats[plume_bool])

                                plume_lons = lons[plume_bool]
                                plume_lats = lats[plume_bool]

                                search_date = datetimes[date_i-steps_back]
                                #print search_date

                                found_plume = False
                                for key in plume_archive_keys:
                                    #print 'Searching in plume archive'
                                    if search_centroid_lon == \
                                            plume_archive[\
                                            str(key)].centroid_lon and \
                                            search_centroid_lat ==\
                                            plume_archive[str(
                                                key)].centroid_lat and \
                                            plume_archive[str(
                                                key)].dates_observed[-1] == \
                                                    search_date:
                                        #print 'Found it in plume archive. ' \
                                        #      'ID is', \
                                        #    plume_archive[str(key)].plume_id
                                        found_plume = True

                                        correct_plume = plume_archive[str(key)]
                                        plume_to_append = plume
                                        # Append the flickered plume to the
                                        # old one which was archived
                                        correct_plume.append_missing_plume(
                                            plume_to_append)

                                        # Add it to plume objects and remove
                                        # it from archives
                                        plume_objects[str(key)] = correct_plume
                                        del plume_archive[str(key)]

                                        # Add it to old IDs, replacing the
                                        # ID of the plume which was found to be
                                        # flickered
                                        flicker_ids.append(plume.plume_id)
                                        reintroduced_ids.append(key)
                                        missing_plume = []

                                        # Reintroduced plumes also get removed
                                        # from the record of the last 50 ids
                                        index = np.argwhere(last_50_ids
                                                            == key)
                                        last_50_ids = np.delete(last_50_ids,
                                                                index)

                                        break
                                # If we didn't find the plume in the plume
                                # archive, it must still be active
                                if found_plume == False:
                                    plume_object_keys = np.asarray([int(i) for
                                                                     i in
                                                                plume_objects])
                                    # Sort the keys in reverse order as the
                                    # plume we want is most likely to have a
                                    # high ID
                                    plume_object_keys[::-1].sort()
                                    for key in plume_object_keys:
                                        #print 'Searching in plume objects'

                                        if search_centroid_lon == \
                                                plume_objects[ \
                                                        str(
                                                    key)].centroid_lon and \
                                                    search_centroid_lat == \
                                                    plume_objects[str(
                                                    key)].centroid_lat and \
                                                    plume_objects[str(
                                                    key)].dates_observed[
                                                    -1] == \
                                                    search_date:
                                            #print 'Found it in plume ' \
                                            #      'objects. ID is', \
                                            #    plume_archive[
                                            #        str(key)].plume_id


                                            found_plume = True

                                            correct_plume = plume_archive[
                                                str(key)]
                                            plume_to_append = plume
                                            # Append the flickered plume to the
                                            # old one which was archived
                                            correct_plume.append_missing_plume(
                                                plume_to_append)

                                            # Add it to plume objects and
                                            # remove it from archives
                                            plume_objects[
                                                str(key)] = correct_plume
                                            del plume_archive[str(key)]

                                            # Add it to old IDs, replacing the
                                            # ID of the plume which was found
                                            # to be flickered
                                            flicker_ids.append(plume.plume_id)
                                            reintroduced_ids.append(key)
                                            missing_plume = []

                                            index = np.argwhere(last_50_ids
                                                                == key)
                                            last_50_ids = np.delete(
                                                last_50_ids,
                                                index)

                                            break

                                break


                # Remove any new IDs which were actually flickers
                for i in np.arange(0, len(flicker_ids)):
                    index = np.argwhere(new_ids==flicker_ids[i])
                    new_ids = np.delete(new_ids, index)

                print merge_ids

                # For merged IDs, we move the tracks to pre-merge tracks
                for i in np.arange(0, len(merge_ids)):
                    plume = plume_objects[str(merge_ids[i])]
                    plume.merge()

                # For old IDs, we just run an update.
                for i in np.arange(0, len(old_ids)):
                    print 'Updating plume', old_ids[i]
                    plume = plume_objects[str(old_ids[i])]
                    #if plume.plume_id == 2:
                    #    print plume.dates_observed
                    plume.update_position(lats, lons, sdf_plumes, old_ids[i])
                    plume.update_duration(datetimes[date_i])
                    if plume.plume_id == 18:
                        print plume.dates_observed
                    plume.update_bbox()
                    plume.update_majorminor_axes()
                    plume.update_area()
                    plume.update_centroid_speed()
                    plume.update_centroid_direction()
                    plume.update_max_extent()
                    # plume.update_leading_edge_4(sdf_plumes, lons, lats)
                    plume_objects[str(old_ids[i])] = plume

                # Plumes which no longer exist are removed and archived
                if len(ids_previous) == 0:
                    removed_ids = []
                else:
                    removed_bool = np.asarray([j not in plume_ids for j in
                                               ids_previous])
                    removed_ids = ids_previous[removed_bool]

                for i in np.arange(0, len(removed_ids)):
                    print 'Archiving plume', removed_ids[i]
                    plume = plume_objects[str(removed_ids[i])]
                    plume.update_GPE_speed()
                    # plume.update_mechanism_likelihood()
                    plume.update_llj_probability(trace)
                    plume_archive[str(removed_ids[i])] = plume
                    #print 'Plume', plume.plume_id,'removed. Final lons and \
                    #                                             lats:'
                    #print plume.centroid_lon
                    #print plume.centroid_lat
                    del plume_objects[str(removed_ids[i])]
                    last_50_ids = np.append(last_50_ids, removed_ids[i])

                if len(last_50_ids) > 50:
                    last_50_ids = last_50_ids[0:50]

                if len(np.unique(sdf_plumes)) < 2:
                    sdf_previous = None
                    ids_previous = []
                    raw_sdf_prev = []
                else:
                    sdf_previous = sdf_plumes
                    ids_previous = plume_ids
                    raw_sdf_prev = sdf_now

                # if runtime > datetime.timedelta(hours=0):
                #    plotting.plot_plumes(plume_objects, sdf_plumes, lats, lons, bt,
                #                         datetimes[date_i], datestrings[date_i])

            else:
                print 'Adding date to list of missing dates'
                with open('missing_dates.txt', 'a') as my_file:
                    my_file.write('\n'+datestrings[date_i])
        # After the script has finished, add remaining plumes to the plume archive
        for i in plume_objects:
            plume_archive[i] = plume_objects[i]

        plume_archive.close()

        # There seems to be a residual massive file - remove this

if __name__ == '__main__':

    run = True
    run_mcmc = False
    no_trace = True
    pickup = False

    if run_mcmc:
        data = get_llj_prob_model.create_plume_dataframe(
            'LLJ_manual_ID.csv', 'plume_archive_LLJ')

        trace = get_llj_prob_model.logistic_llj_model(data, 'LLJ',
                                                      'conv_distance',
                                                      'time_to_09',
                                                      'duration')

    if no_trace:
        trace = None

    # Here specify a set of years and months for which you will run.
    # Then the wrapper function is just all the below, taking month_lower,
    # month_upper, year_lower, year_upper as an argument and specifying the
    # other time bounds

    years = [2004, 2005, 2006, 2007, 2008, 2009, 2010]
    months = [[6, 7, 8], [6, 7, 8], [6, 7, 8], [6, 7, 8], [6, 7, 8], [6, 7, 8],
              [6, 7, 8]]
    yearmonths = zip(years, months)

    wrapper(yearmonths[6])

    """
    processes = [multiprocessing.Process(target=wrapper, args=(i,))
                 for i in yearmonths]
    for p in processes:
        p.start()
    for p in processes:
        p.join()
    """



