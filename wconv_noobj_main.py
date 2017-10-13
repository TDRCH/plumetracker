import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import numpy as np
from mpl_toolkits.basemap import Basemap
import datetime
import shelve

import utilities
import plumes
import pinkdust
import plotting

if __name__ == '__main__':
    #### NOTES ####
    # These parameters should be handled in a spreadsheet
    # Plotting should be handled in a separate function
    # If you have to give the sdf_plumes to the plotting function you won't
    # be able to run it as a standalone after having generated all the plumes
    # Well you can't do that anyway without having a date associated with
    # each plume - actually you could work it out from the duration
    # So you can store the valid lons and valid lats for each plume,
    # then you'd need to use a boolean on the lons and lats, but in that
    # case they'd have to match exactly. So you could look for the nearest
    # matching lat and lon in the plotting lats and lons?
    # If a plume merges, its centroid track should be cleared

    # Cloud mask is to be read in from grib and regridded with a function in
    #  pinkdust. Then you feed it to the scan_for_plumes function, which'll
    # return a whole pile of IDs. Those in turn can go to the convection
    # object generator.

    year_lower = 2012
    year_upper = 2012
    month_lower = 6
    month_upper = 6
    day_lower = 1
    day_upper = 30
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
    cloud_ids_previous = []
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

    plume_objects = shelve.open('plume_objects_clouds')
    plume_objects.clear()

    convection_objects = shelve.open('convection_cloud_objects')
    convection_objects.clear()

    plume_archive = shelve.open('plume_archive_clouds')
    cloud_archive = shelve.open('cloud_archive')

    for date_i in np.arange(0, len(datestrings)):
        runtime = datetimes[date_i] - datetimes[0]
        print '\n' + datestrings[date_i] + '\n'
        totaltest = datetime.datetime.now()
        sdf = Dataset(
            '/ouce-home/data/satellite/meteosat/seviri/15-min/native/sdf/nc/'
            'JUNE2012/SDF_v2/SDF_v2.' + datestrings[date_i] + '.nc')
        bt = Dataset(
            '/ouce-home/data/satellite/meteosat/seviri/15-min/native/bt/nc/'
            'JUNE2012/H-000-MSG2__-MSG2________-'
            'IR_BrightnessTemperatures___-000005___-' + datestrings[date_i] +
            '-__.nc')

        cloudnc = Dataset('/soge-home/data/satellite/meteosat/seviri/15-min/'
                        'native/cloudmask/nc/JUNE2012_CLOUDS/eumetsat.cloud'
                             '.'+datetimes[date_i]
            .strftime(
            "%Y%m%d%H%M")+'.nc')

        sdf_now = sdf.variables['bt108'][:]
        clouds_now = cloudnc.variables['cmask'][:]
        clouds_now[clouds_now == 1] = 0
        clouds_now[clouds_now == 3] = 0
        clouds_now[clouds_now == 2] = 1

        # Get plumes first by scanning for them
        sdf_plumes, new_ids, plume_ids, merge_ids = plumes.scan_for_plumes(
            sdf_now,
            sdf_previous,
            used_ids)

        # Get clouds by scanning for them
        clouds, new_cloud_ids, cloud_ids, merge_cloud_ids = \
            plumes.scan_for_plumes(clouds_now,
                                   clouds_previous,
                                   used_cloud_ids)

        # Now just need a way whereby the plume objects can interact with
        # cloud objects, e.g. check if they are near them or check if clouds
        #  have been initiated on top of plumes, or check their propagaton
        # relative to that of the cloud idk

        # Then, for each new ID, we initialise plume objects
        for i in np.arange(0, len(new_ids)):
            #print 'Creating new plume', new_ids[i]
            plume = plumes.Plume(new_ids[i], datetimes[date_i])
            plume.update_position(lats, lons, sdf_plumes, new_ids[i])
            plume.update_duration(datetimes[date_i])
            plume.update_bbox()
            plume.update_majorminor_axes()
            plume.update_area()
            plume.update_max_extent()
            plume.update_centroid_speed()
            plume.update_centroid_direction()
            #plume.update_leading_edge_4(sdf_plumes, lons, lats)
            plume_objects[str(new_ids[i])] = plume

        # For merged IDs, we move the tracks to pre-merge tracks
        for i in np.arange(0, len(merge_ids)):
            plume = plume_objects[str(merge_ids[i])]
            plume.merge()

        # For old IDs, we just run an update.
        for i in np.arange(0, len(old_ids)):
            #print 'Updating plume', old_ids[i]
            plume = plume_objects[str(old_ids[i])]
            plume.update_position(lats,lons, sdf_plumes, old_ids[i])
            plume.update_duration(datetimes[date_i])
            plume.update_bbox()
            plume.update_majorminor_axes()
            plume.update_area()
            plume.update_centroid_speed()
            plume.update_centroid_direction()
            plume.update_max_extent()
            #plume.update_leading_edge_4(sdf_plumes, lons, lats)
            plume_objects[str(old_ids[i])] = plume

        # Plumes which no longer exist are removed and archived
        if len(ids_previous) == 0:
            removed_ids = []
        else:
            removed_bool = np.asarray([j not in plume_ids for j in
                                       ids_previous])
            removed_ids = ids_previous[removed_bool]

        for i in np.arange(0, len(removed_ids)):
            #print 'Archiving plume', removed_ids[i]
            plume = plume_objects[str(removed_ids[i])]
            plume.update_GPE_speed()
            plume_archive[str(removed_ids[i])] = plume
            del plume_objects[str(removed_ids[i])]

        sdf_previous = sdf_plumes
        ids_previous = plume_ids

        if runtime > datetime.timedelta(hours=0):
            plotting.plot_plumes(convection_objects, clouds, lats, lons, bt,
                                 datetimes[date_i], datestrings[date_i])

    # After the script has finished, add remaining plumes to the plume archive
    for i in plume_objects:
        plume_archive[i] = plume_objects[i]

    for i in convection_objects:
        cloud_archive[i] = convection_objects[i]

    # Summary plots

    plume_objects.close()
    plume_archive.close()

    convection_objects.close()
    cloud_archive.close()


