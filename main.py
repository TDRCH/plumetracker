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
    #

    SDF_0 = Dataset(
        '/ouce-home/data/satellite/meteosat/seviri/15-min/native/sdf/nc'
        '/JUNE2012/SDF_v2/SDF_v2.201206221800.nc')

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
    lons = lonlats.variables['longitude'][:]
    lats = lonlats.variables['latitude'][:]
    lonmask = lons > 360
    latmask = lats > 90
    lons = np.ma.array(lons, mask=lonmask)
    lats = np.ma.array(lats, mask=latmask)
    fig, ax = plt.subplots()
    extent = (-21, 31, 10, 41)
    m = Basemap(projection='cyl', llcrnrlon=extent[0], urcrnrlon=extent[1],
                llcrnrlat=extent[2], urcrnrlat=extent[3],
                resolution='i')

    m.drawcoastlines(linewidth=0.5)
    m.drawcountries(linewidth=0.5)

    sdf_previous = None
    ids_previous = []
    deep_conv_IDs_prev = None
    LLJ_plumes_IDs_prev = []
    k = 0
    available_colours = np.arange(0, 41)
    used_ids = []
    used_colour_IDs = {}
    plume_objects = []

    plume_objects = shelve.open('plume_objects')
    plume_objects.clear()

    convection_objects = shelve.open('convection_objects')
    convection_objects.clear()

    plume_archive = shelve.open('plume_archive')

    for date_i in np.arange(0, len(datestrings)):
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

        sdf_now = sdf.variables['bt108'][:]

        sdf_plumes, new_ids, plume_ids = plumes.scan_for_plumes(sdf_now,
                                                              sdf_previous,
                                                                used_ids)

        for i in new_ids:
            used_ids.append(i)

        old_bool = np.asarray([j in ids_previous for j in plume_ids])
        old_ids = plume_ids[old_bool]

        # Then, for each new ID, we initialise plume objects
        for i in np.arange(0, len(new_ids)):
            print 'Creating new plume', new_ids[i]
            plume = plumes.Plume(new_ids[i], datetimes[date_i])
            plume.update_position(lats, lons, sdf_plumes, new_ids[i])
            plume.update_duration(datetimes[date_i])
            plume.update_bbox()
            plume.update_majorminor_axes()
            plume.update_leading_edge(sdf_plumes)
            plume_objects[str(new_ids[i])] = plume

        # For old IDs, we just run an update.
        for i in np.arange(0, len(old_ids)):
            print 'Updating plume', old_ids[i]
            plume = plume_objects[str(old_ids[i])]
            plume.update_position(lats,lons, sdf_plumes, old_ids[i])
            plume.update_duration(datetimes[date_i])
            plume.update_bbox()
            plume.update_majorminor_axes()
            plume.update_leading_edge(sdf_plumes)
            plume_objects[str(old_ids[i])] = plume

        # Plumes which no longer exist are removed
        if len(ids_previous) == 0:
            removed_ids = []
        else:
            removed_bool = np.asarray([j not in plume_ids for j in
                                       ids_previous])
            removed_ids = ids_previous[removed_bool]

        for i in np.arange(0, len(removed_ids)):
            print 'Archiving plume', removed_ids[i]
            plume = plume_objects[str(removed_ids[i])]
            plume_archive[str(removed_ids[i])] = plume
            del plume_objects[str(removed_ids[i])]

        sdf_previous = sdf_plumes
        ids_previous = plume_ids



        plotting.plot_plumes(plume_objects, sdf_plumes, lats, lons, bt,
                             datetimes[date_i], datestrings[date_i])

    plume_objects.close()


