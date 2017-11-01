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
    clouds_previous = None
    ids_previous = []
    cloud_ids_previous = []
    deep_conv_IDs_prev = None
    LLJ_plumes_IDs_prev = []
    k = 0
    available_colours = np.arange(0, 41)
    used_ids = []
    used_cloud_ids = []
    used_colour_IDs = {}
    plume_objects = []

    # Restrict the lons and lats to the CWS alone
    lonbool = np.asarray([j >= -20 and j <= 10 for j in lons[0]])

    if run:

        plume_objects = {}

        plume_archive = shelve.open('plume_archive_'+str(yearmonth[0]))

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

                if os.path.isfile('/ouce-home/data/satellite/meteosat/seviri/15-min/'
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

                # Get plumes first by scanning for them
                sdf_plumes, new_ids, plume_ids, merge_ids = plumes.\
                    scan_for_plumes(
                    sdf_now,
                    sdf_previous,
                    used_ids,
                    clouds)

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
                    # print 'Creating new plume', new_ids[i]
                    plume = plumes.Plume(new_ids[i], datetimes[date_i])
                    plume.update_position(lats, lons, sdf_plumes, new_ids[i])
                    plume.update_duration(datetimes[date_i])
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

                # For merged IDs, we move the tracks to pre-merge tracks
                for i in np.arange(0, len(merge_ids)):
                    plume = plume_objects[str(merge_ids[i])]
                    plume.merge()

                # For old IDs, we just run an update.
                for i in np.arange(0, len(old_ids)):
                    # print 'Updating plume', old_ids[i]
                    plume = plume_objects[str(old_ids[i])]
                    plume.update_position(lats, lons, sdf_plumes, old_ids[i])
                    plume.update_duration(datetimes[date_i])
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
                    # print 'Archiving plume', removed_ids[i]
                    plume = plume_objects[str(removed_ids[i])]
                    plume.update_GPE_speed()
                    # plume.update_mechanism_likelihood()
                    plume.update_llj_probability(trace)
                    plume_archive[str(removed_ids[i])] = plume
                    del plume_objects[str(removed_ids[i])]

                if len(np.unique(sdf_plumes)) < 2:
                    sdf_previous = None
                    ids_previous = []
                else:
                    sdf_previous = sdf_plumes
                    ids_previous = plume_ids

                # if runtime > datetime.timedelta(hours=0):
                #    plotting.plot_plumes(plume_objects, sdf_plumes, lats, lons, bt,
                #                         datetimes[date_i], datestrings[date_i])

            else:
                print 'Adding date to list of missing dates'
                with open('missing_dates.txt', 'a') as my_file:
                    my_file.write(datestrings[date_i])
        # After the script has finished, add remaining plumes to the plume archive
        for i in plume_objects:
            plume_archive[i] = plume_objects[i]

        plume_archive.close()

        # There seems to be a residual massive file - remove this

if __name__ == '__main__':

    run = True
    run_mcmc = False
    no_trace = True

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

    #wrapper(yearmonths[6])

    processes = [multiprocessing.Process(target=wrapper, args=(i,))
                 for i in yearmonths]
    for p in processes:
        p.start()
    for p in processes:
        p.join()




