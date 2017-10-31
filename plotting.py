import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import numpy as np
from netCDF4 import Dataset
from mpop.satout import netcdf4
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.cbook as cbook
import matplotlib.cm as cm
import os
import pymc3 as pm

import pinkdust
import utilities

def plot_plume_overviews(plume_archive, lats, lons):
    """
    Plots plume overview plots for each timestep used in the plume archive
    :param plume_archive:
    :return:
    """

    fig, ax = plt.subplots()
    extent = (np.min(lons), np.max(lons), np.min(lats), np.max(lats))
    m = Basemap(projection='cyl', llcrnrlon=extent[0], urcrnrlon=extent[1],
                llcrnrlat=extent[2], urcrnrlat=extent[3],
                resolution='i')

    m.drawcoastlines(linewidth=0.5)
    m.drawcountries(linewidth=0.5)

    # All you need is a list of IDs active for each datetime
    date_dictionary = {}

    for i in plume_archive:
        dates_observed = plume_archive[i].dates_observed
        for j in dates_observed:
            if j in date_dictionary:
                date_dictionary[j].append(i)
            else:
                date_dictionary[j] = []
                date_dictionary[j].append(i)

    total_dates = len(plume_archive)
    plots_done = 0

    # Then make one plot per date
    for i in date_dictionary:
        print '\n\nProgress:', str((float(plots_done)/float(
            total_dates))*100)+'%\n\n'
        plume_data = np.zeros((lats.shape[0], lons.shape[1]))
        plume_data[:] = np.nan
        centroid_lats = []
        centroid_lons = []
        plume_ids = []

        active_plumes = date_dictionary[i]

        for j in active_plumes:

            dates_observed = plume_archive[j].dates_observed
            date_bool = np.asarray(dates_observed) == i

            if i.hour == 9 and i.minute == 30:
                print plume_archive[j].plume_id
                print plume_archive[j].dates_observed

            plume_bools = np.asarray(plume_archive[j].track_plume_bool)

            plume_bool = plume_bools[date_bool]

            plume_bool = plume_bool[0].toarray()

            if plume_archive[j].clear_LLJ:
                plume_data[plume_bool] = 2
            else:
                plume_data[plume_bool] = 1

            track_centroid_lats = plume_archive[j].track_centroid_lat
            track_centroid_lons = plume_archive[j].track_centroid_lon
            pre_merge_centroid_lats = plume_archive[
                j].pre_merge_track_centroid_lat
            pre_merge_centroid_lons = plume_archive[
                j].pre_merge_track_centroid_lon

            if plume_archive[j].merged == False:
                centroid_lat = np.asarray(track_centroid_lats)[date_bool][0]
                centroid_lon = np.asarray(track_centroid_lons)[date_bool][0]
            else:
                combined_track_lat = np.append(pre_merge_centroid_lats,
                                               track_centroid_lats)
                combined_track_lon = np.append(pre_merge_centroid_lons,
                                               track_centroid_lons)
                centroid_lat = np.asarray(combined_track_lat)[date_bool][0]
                centroid_lon = np.asarray(combined_track_lon)[date_bool][0]

            centroid_lats.append(centroid_lat)
            centroid_lons.append(centroid_lon)
            plume_ids.append(plume_archive[j].plume_id)

        centroid_x, centroid_y = m(centroid_lons, centroid_lats)

        levels = [1, 2]

        contourplot = m.contourf(lons, lats, plume_data, levels, colors=(
            'r', 'g'))
        centroidplot = m.scatter(centroid_x, centroid_y, s=1)

        anns = []

        for label_i, txt in enumerate(plume_ids):
            ann = ax.annotate(txt, (centroid_x[label_i], centroid_y[
                label_i]), fontsize=8)
            anns.append(ann)

        plt.savefig('Plume_overview_'+i.strftime("%Y%m%d%H%M")+'.png')

        centroidplot.remove()
        for coll in contourplot.collections:
            coll.remove()

        for ann in anns:
            ann.remove()

        plots_done += 1

def plot_plumes(plume_objects, sdf_plumes, lats, lons, bt, datetime, \
                                          datestring):
    """
    Plots a pinkdust RGB composite with SDF plumes overlain
    :param sdf_plumes:
    :param lats:
    :param lons:
    :param pinkdust_rgb:
    :return:
    """

    fig, ax = plt.subplots()
    extent = (np.min(lons), np.max(lons), np.min(lats), np.max(lats))
    m = Basemap(projection='cyl', llcrnrlon=extent[0], urcrnrlon=extent[1],
                llcrnrlat=extent[2], urcrnrlat=extent[3],
                resolution='i')

    m.drawcoastlines(linewidth=0.5)
    m.drawcountries(linewidth=0.5)

    """
    data_array = np.zeros((lons.shape[0], lats.shape[1], 3))


    data_array[:, :, 0] = bt.variables['bt087'][:]
    data_array[:, :, 1] = bt.variables['bt108'][:]
    data_array[:, :, 2] = bt.variables['bt120'][:]

    data_regridded = pinkdust.regrid_data_to_regular(lons, lats,
                                                 data_array)

    rgb_array = pinkdust.generate_image_from_array(data_regridded)

    m.imshow(rgb_array, extent=extent, origin='lower')
    """

    for id in plume_objects:
        plume_bool = sdf_plumes == plume_objects[id].plume_id
        plume_data = np.zeros((lats.shape[0], lons.shape[1]))
        plume_data[:] = np.nan
        plume_data[plume_bool] = plume_objects[id].plume_id
        m.contourf(lons, lats, plume_data)

        if plume_objects[id].leading_edge_lon != None:
            leading_edge_x, leading_edge_y = m(plume_objects[
                                                     id].leading_edge_lon,
                                                 plume_objects[
                                                     id].leading_edge_lat)
            m.plot(leading_edge_x, leading_edge_y,
                   linewidth=1)

        if len(plume_objects[id].track_lons) > 0:
            print plume_objects[id].track_centroid_lat
            print plume_objects[id].track_centroid_lon
            track_x, track_y = m(plume_objects[id].track_centroid_lon,
                                 plume_objects[id].track_centroid_lat)

            m.plot(track_x, track_y, linewidth=1)

    plt.savefig('pinkdust_'+datestring+'.png')
    plt.close()

def plot_existing_SDF():
    root = '/ouce-home/data/satellite/meteosat/seviri/15-min/native/sdf/' \
           'nc/JUNE2012/SDF_v2/'
    lonlats = Dataset('/ouce-home/data/satellite/meteosat/seviri/15-min/'
                      'native/lonlats.NA_MiddleEast.nc')
    lons = lonlats.variables['longitude'][:]
    lats = lonlats.variables['latitude'][:]
    lonmask = lons > 360
    latmask = lats > 90
    lons = np.ma.array(lons, mask=lonmask)
    lats = np.ma.array(lats, mask=latmask)

    extent = (np.min(lons), np.max(lons), np.min(lats), np.max(lats))
    m = Basemap(projection='cyl', llcrnrlon=extent[0], urcrnrlon=extent[1],
                llcrnrlat=extent[2], urcrnrlat=extent[3],
                resolution='i')

    m.drawcoastlines(linewidth=0.5)
    m.drawcountries(linewidth=0.5)

    for item in os.listdir(root):
        if os.path.isfile(os.path.join(root, item)):
            SDF = Dataset(root + str(item))
            c = m.contourf(lons, lats, SDF.variables['bt108'][:])
            plt.savefig(str(item)[:-3] + '.png')
            print str(item)[:-3]
            for coll in c.collections:
                plt.gca().collections.remove(coll)

# Function to plot active plume number through time
def plot_plume_count(plume_archive, min_lat=None, min_lon=None,
                    max_lat=None, max_lon=None):
    """
    Plots the number of active plumes at each timestep in the plume archive
    :param plume_archive:
    :return:
    """

    plt.close()
    fig, ax = plt.subplots()

    if min_lat == None:
        # Extract all the dates in use in the plume archive
        all_dates = []
        for i in plume_archive:
            dates = plume_archive[i].dates_observed
            for j in dates:
                all_dates.append(j)
    else:
        all_dates = []
        for i in plume_archive:
            plume_centroids_lat = plume_archive[i].track_centroid_lat
            plume_centroids_lon = plume_archive[i].track_centroid_lon
            include_bool_lat = np.asarray([j >= min_lat and j <= max_lat for j
                                           in
                                           plume_centroids_lat])
            include_bool_lon = np.asarray([j >= min_lon and j <= max_lon for j
                                           in plume_centroids_lon])
            combined_bool = np.asarray([a and b for a, b in zip(
            include_bool_lat, include_bool_lon)])
            dates = np.asarray(plume_archive[i].dates_observed)[combined_bool]
            for j in dates:
                all_dates.append(j)

    unique_dates = np.unique(all_dates)
    unique_dates = np.sort(unique_dates)
    date_count = []

    for i in unique_dates:
        date_count.append(all_dates.count(i))

    ax.plot(unique_dates, date_count, color='black')

    plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%d/%m/%Y'))
    plt.gca().xaxis.set_major_locator(mdates.DayLocator(interval=2))
    #plt.gca().xaxis.set_minor_locator(mdates.HourLocator())

    plt.gcf().autofmt_xdate()

    ax.set_ylabel('Plume count')
    ax.set_xlabel('Time')
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    ax.set_ylim(0, 25)
    plt.grid()

    plt.savefig('plume_count_through_time.png')
    plt.close()

def plot_plume_area(plume_archive, min_lat=None, min_lon=None,
                    max_lat=None, max_lon=None):
    """
    Plots the mean area of emitted plumes through time
    Note that if there are many plumes they will simply outweight the large
    ones
    :param plume_archive:
    :return:
    """

    plt.close()
    fig, ax = plt.subplots()

    date_dictionary = {}

    for i in plume_archive:
        if min_lat != None:
            plume_centroids_lat = plume_archive[i].track_centroid_lat
            plume_centroids_lon = plume_archive[i].track_centroid_lon
            include_bool_lat = np.asarray([j >= min_lat and j <= max_lat for j
                                           in
                                       plume_centroids_lat])
            include_bool_lon = np.asarray([j >= min_lon and j <= max_lon for j
                                           in plume_centroids_lon])
            combined_bool = np.asarray([a and b for a, b in zip(
                include_bool_lat, include_bool_lon)])

            areas = np.asarray(plume_archive[i].track_area)[combined_bool]
            dates = np.asarray(plume_archive[i].dates_observed)[combined_bool]
        else:
            areas = np.asarray(plume_archive[i].track_area)
            dates = np.asarray(plume_archive[i].dates_observed)

        # If this date is not present in the dictionary, create a list and
        # append areas to it. Otherwise just add entries to the existing key
        for j in np.arange(0, len(dates)):
            if dates[j] in date_dictionary:
                date_dictionary[dates[j]].append(areas[j])
            else:
                date_dictionary[dates[j]] = []
                date_dictionary[dates[j]].append(areas[j])

    # Average each entry
    for i in date_dictionary:
        area_list = date_dictionary[i]
        averaged_area = np.nanmean(area_list)
        date_dictionary[i] = averaged_area

    # Plot
    lists = sorted(date_dictionary.items())
    x, y = zip(*lists)

    ax.plot(x, y, color='black')

    plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%d/%m/%Y'))
    plt.gca().xaxis.set_major_locator(mdates.DayLocator())
    plt.gca().xaxis.set_minor_locator(mdates.HourLocator())

    plt.gcf().autofmt_xdate()

    ax.set_ylabel('Mean plume area (m^2)')
    ax.set_xlabel('Time')
    plt.grid()

    plt.savefig('plume_area_through_time.png')
    plt.close()

def plot_plume_total_area(plume_archive, min_lat=None, min_lon=None,
                    max_lat=None, max_lon=None):
    """
    Plots the total area of emitted plumes through time
    Note that if there are many plumes they will simply outweigh the large
    ones
    :param plume_archive:
    :return:
    """

    plt.close()
    fig, ax = plt.subplots()

    date_dictionary = {}

    for i in plume_archive:
        if min_lat != None:
            plume_centroids_lat = plume_archive[i].track_centroid_lat
            plume_centroids_lon = plume_archive[i].track_centroid_lon
            include_bool_lat = np.asarray([j >= min_lat and j <= max_lat for j
                                           in
                                       plume_centroids_lat])
            include_bool_lon = np.asarray([j >= min_lon and j <= max_lon for j
                                           in plume_centroids_lon])
            combined_bool = np.asarray([a and b for a, b in zip(
                include_bool_lat, include_bool_lon)])

            areas = np.asarray(plume_archive[i].track_area)[combined_bool]
            dates = np.asarray(plume_archive[i].dates_observed)[combined_bool]
        else:
            areas = np.asarray(plume_archive[i].track_area)
            dates = np.asarray(plume_archive[i].dates_observed)

        # If this date is not present in the dictionary, create a list and
        # append areas to it. Otherwise just add entries to the existing key
        for j in np.arange(0, len(dates)):
            if dates[j] in date_dictionary:
                date_dictionary[dates[j]].append(areas[j])
            else:
                date_dictionary[dates[j]] = []
                date_dictionary[dates[j]].append(areas[j])

    # Average each entry
    for i in date_dictionary:
        area_list = date_dictionary[i]
        averaged_area = np.nansum(area_list)
        date_dictionary[i] = averaged_area

    # Plot
    lists = sorted(date_dictionary.items())
    x, y = zip(*lists)

    ax.plot(x, y, color='black')

    plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%d/%m/%Y'))
    plt.gca().xaxis.set_major_locator(mdates.DayLocator(interval=2))
    #plt.gca().xaxis.set_minor_locator(mdates.HourLocator())

    #plt.gcf().autofmt_xdate()

    ax.set_ylabel('Total plume area (m^2)')
    ax.set_xlabel('Time')
    plt.xticks(rotation=45, ha='right')
    plt.grid()
    plt.tight_layout()

    plt.savefig('plume_total_area_through_time.png')
    plt.close()

def plot_plume_centroid_speed(plume_archive, min_lat=None, min_lon=None,
                    max_lat=None, max_lon=None):
    """
    Plots the mean speed of centroids through time
    :param plume_archive:
    :return:
    """

    plt.close()
    fig, ax = plt.subplots()

    date_dictionary = {}

    for i in plume_archive:
        if min_lat != None:
            plume_centroids_lat = plume_archive[i].track_centroid_lat
            plume_centroids_lon = plume_archive[i].track_centroid_lon
            include_bool_lat = np.asarray([j >= min_lat and j <= max_lat for j
                                           in
                                       plume_centroids_lat])
            include_bool_lon = np.asarray([j >= min_lon and j <= max_lon for j
                                           in plume_centroids_lon])
            combined_bool = np.asarray([a and b for a, b in zip(
                include_bool_lat, include_bool_lon)])

            speeds = np.asarray(plume_archive[i].track_speed_centroid)[
                combined_bool]
            dates = np.asarray(plume_archive[i].dates_observed)[combined_bool]
        else:
            speeds = np.asarray(plume_archive[i].track_speed_centroid)
            dates = np.asarray(plume_archive[i].dates_observed)

        # If this date is not present in the dictionary, create a list and
        # append areas to it. Otherwise just add entries to the existing key
        for j in np.arange(0, len(dates)):
            if dates[j] in date_dictionary:
                date_dictionary[dates[j]].append(speeds[j])
            else:
                date_dictionary[dates[j]] = []
                date_dictionary[dates[j]].append(speeds[j])

    # Average each entry
    for i in date_dictionary:
        speed_list = date_dictionary[i]
        averaged_speed = np.nanmean(speed_list)
        date_dictionary[i] = averaged_speed

    # Plot
    lists = sorted(date_dictionary.items())
    x, y = zip(*lists)

    ax.plot(x, y, color='black')

    plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%d/%m/%Y'))
    plt.gca().xaxis.set_major_locator(mdates.DayLocator(interval=2))

    # plt.gcf().autofmt_xdate()

    ax.set_ylabel('Plume speed (m/s)')
    ax.yaxis.set_tick_params()
    ax.set_xlabel('Time')
    plt.grid()
    plt.ylim(0, 100)
    plt.xticks(rotation=45, ha='right')
    ax.xaxis.set_tick_params()
    plt.tight_layout()

    plt.savefig('plume_speed_through_time.png')
    plt.close()

def plot_plume_centroid_direction(plume_archive, min_lat=None, min_lon=None,
                    max_lat=None, max_lon=None):
    """
    Plots the mean direction of centroids through time
    :param plume_archive:
    :return:
    """

    plt.close()
    fig, ax = plt.subplots()

    date_dictionary = {}

    for i in plume_archive:
        if min_lat != None:
            plume_centroids_lat = plume_archive[i].track_centroid_lat
            plume_centroids_lon = plume_archive[i].track_centroid_lon
            include_bool_lat = np.asarray([j >= min_lat and j <= max_lat for j
                                           in
                                       plume_centroids_lat])
            include_bool_lon = np.asarray([j >= min_lon and j <= max_lon for j
                                           in plume_centroids_lon])
            combined_bool = np.asarray([a and b for a, b in zip(
                include_bool_lat, include_bool_lon)])

            directions = np.asarray(plume_archive[i].track_centroid_direction)[
                combined_bool]
            dates = np.asarray(plume_archive[i].dates_observed)[combined_bool]
        else:
            directions = np.asarray(plume_archive[i].track_centroid_direction)
            dates = np.asarray(plume_archive[i].dates_observed)

        # If this date is not present in the dictionary, create a list and
        # append areas to it. Otherwise just add entries to the existing key
        for j in np.arange(0, len(dates)):
            if dates[j] in date_dictionary:
                date_dictionary[dates[j]].append(directions[j])
            else:
                date_dictionary[dates[j]] = []
                date_dictionary[dates[j]].append(directions[j])

    # Average each entry
    for i in date_dictionary:
        direction_list = date_dictionary[i]
        averaged_area = np.nanmean(direction_list)
        date_dictionary[i] = averaged_area

    # Plot
    lists = sorted(date_dictionary.items())
    x, y = zip(*lists)

    ax.plot(x, y, color='black')

    plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%d/%m/%Y'))
    plt.gca().xaxis.set_major_locator(mdates.DayLocator(interval=2))
    #plt.gca().xaxis.set_minor_locator(mdates.HourLocator())

    # plt.gcf().autofmt_xdate()

    ax.set_ylabel('Plume direction (degrees)')
    #ax.yaxis.set_tick_params(labelsize=20)
    ax.set_xlabel('Time')
    plt.grid()
    plt.ylim(0, 360)
    plt.xticks(rotation=45, ha='right')
    #ax.xaxis.set_tick_params()
    plt.tight_layout()

    plt.savefig('plume_direction_through_time.png')
    plt.close()

def plot_emission_speed_map(plume_archive, lats, lons, min_lat=None,
                            min_lon=None,
                    max_lat=None, max_lon=None):

    """
    For each gridbox in the SDF region of choice, plots the mean
    emission speed of plumes emitted from that gridbox
    :param plume_archive:
    :param lats:
    :param lons:
    :param min_lat:
    :param min_lon:
    :param max_lat:
    :param max_lon:
    :return:
    """

    plt.close()

    fig, ax = plt.subplots()
    extent = (np.min(lons), np.max(lons), np.min(lats), np.max(lats))
    m = Basemap(projection='cyl', llcrnrlon=extent[0], urcrnrlon=extent[1],
                llcrnrlat=extent[2], urcrnrlat=extent[3],
                resolution='i')

    m.drawcoastlines(linewidth=0.5)
    m.drawcountries(linewidth=0.5)

    latlon_dictionary = {}

    for i in plume_archive:
        # If the plume has merged, the plume source is found in the
        # pre-merge track
        if plume_archive[i].merged == True:
            if len(plume_archive[i].pre_merge_track_speed_centroid) < 4:
                continue
            emission_lat = plume_archive[i].pre_merge_track_centroid_lat[0]
            emission_lon = plume_archive[i].pre_merge_track_centroid_lon[0]
            emission_speeds = plume_archive[i].pre_merge_track_speed_centroid[
                              3]
            if emission_speeds >= 50:
                emission_speeds = np.nan
            print emission_speeds
        else:
            if len(plume_archive[i].track_speed_centroid) < 4:
                continue
            emission_lat = plume_archive[i].track_centroid_lat[0]
            emission_lon = plume_archive[i].track_centroid_lon[0]
            emission_speeds = plume_archive[i].track_speed_centroid[3]
            if emission_speeds >= 50:
                emission_speeds = np.nan

        # If this entry already exists in the dictionary we just append to it
        if (emission_lat, emission_lon) in latlon_dictionary:
            latlon_dictionary[(emission_lat, emission_lon)].append(np.nanmean(
                emission_speeds))
        else:
            latlon_dictionary[(emission_lat, emission_lon)] = []
            latlon_dictionary[(emission_lat, emission_lon)].append(
                np.nanmean(emission_speeds))

    # Average together lists for each entry in the dictionary
    for i in latlon_dictionary:
        latlon_dictionary[i] = np.nanmean(latlon_dictionary[i])

    latlons = np.asarray([i for i in latlon_dictionary])
    lats = np.asarray([j[0] for j in latlons])
    lons = np.asarray([j[1] for j in latlons])
    data = np.asarray([latlon_dictionary[i] for i in latlon_dictionary])

    #lons, lats = np.meshgrid(lons, lats)

    # What you actually need is to bin the data into discrete gridboxes,
    # because the centroid lat/lons will not align with the SDF lat/lons
    # So you round each one to the nearest degree, for instance,

    # So say you just have the right lats and the right lons - how would you
    #  then get a meshgrid?

    # Generate a data array which will be plotted
    #data_array[:] = np.nan

    m.scatter(lons, lats, c=data, cmap=cm.RdYlBu_r)

    plt.savefig('emission_speed_map.png')

    plt.close()

# Function to plot mean plume size through time
# Function to plot plume size distribution by region
# Function to plot mean plume direction for given regions
# Function to plot source region density map
# Function to plot LLJ likelihood and CPO likelihood by region

def plot_emission_speed_map_2(plume_archive, lats, lons, min_lat=None,
                            min_lon=None,
                    max_lat=None, max_lon=None):
    # First we need a grid on which to map all our speeds
    # This can come from lats and lons
    # Second, you set that grid to a whole pile of nans
    # Then you write all your speed values to a nice dictionary
    # Then when you have your dictionary, average each entry
    # Then for each latlon tuple, bin it to the nearest latlon on the grid
    # and write that to the array of nans
    # Then you do like imshow or something, but you'd have to get regular
    # latlons probably, unless there's some other grid plotting method which
    # allows regular latlons

    plt.close()

    fig, ax = plt.subplots()

    latlon_dictionary = {}

    lats = np.linspace(np.min(lats), np.max(lats), 100)
    lons = np.linspace(np.min(lons), np.max(lons), 100)

    lons, lats = np.meshgrid(lons, lats)

    # Get your zip grid of latlon coordinates
    latlon_zip = np.array(zip(lats.ravel(), lons.ravel()),
                          dtype=('f8,'
                                 'f8')).reshape(lats.shape)

    latlon_zip = latlon_zip.flatten()

    for i in plume_archive:
        # If the plume has merged, the plume source is found in the
        # pre-merge track
        if plume_archive[i].merged == True:
            if len(plume_archive[i].pre_merge_track_speed_centroid) < 4:
                continue
            emission_lat = plume_archive[i].pre_merge_track_centroid_lat[0]
            emission_lon = plume_archive[i].pre_merge_track_centroid_lon[0]
            emission_speeds = plume_archive[i].pre_merge_track_speed_centroid[
                3]
            if emission_speeds >= 50:
                emission_speeds = np.nan
        else:
            if len(plume_archive[i].track_speed_centroid) < 4:
                continue
            emission_lat = plume_archive[i].track_centroid_lat[0]
            emission_lon = plume_archive[i].track_centroid_lon[0]
            emission_speeds = plume_archive[i].track_speed_centroid[3]
            if emission_speeds >= 50:
                emission_speeds = np.nan

        nearest = min(latlon_zip, key=lambda x: utilities.
                      haversine(x[1],
                                x[0],
                                emission_lon,
                                emission_lat))

        # If this entry already exists in the dictionary we just append to it
        if (nearest[0], nearest[1]) in latlon_dictionary:
            latlon_dictionary[(nearest[0], nearest[1])].append(np.nanmean(
                emission_speeds))
        else:
            latlon_dictionary[(nearest[0], nearest[1])] = []
            latlon_dictionary[(nearest[0], nearest[1])].append(
                np.nanmean(emission_speeds))

    # Average together lists for each entry in the dictionary
    for i in latlon_dictionary:
        latlon_dictionary[i] = np.nanmean(latlon_dictionary[i])

    # Data array
    data_array = np.zeros((lons.shape))

    # Find the nearest latlon for each dictionary key, and extract the index
    # for it

    for i in latlon_dictionary:
        lat_bool = lats == i[0]
        lon_bool = lons == i[1]
        union_bool = lat_bool & lon_bool
        data_array[union_bool] = latlon_dictionary[i]

    data_array[data_array == 0] = np.nan

    extent = (np.min(lons), np.max(lons), np.min(lats), np.max(lats))
    m = Basemap(projection='cyl', llcrnrlon=extent[0], urcrnrlon=extent[1],
                llcrnrlat=extent[2], urcrnrlat=extent[3],
                resolution='i')

    m.drawcoastlines(linewidth=0.5)
    m.drawcountries(linewidth=0.5)

    #m.contourf(lons, lats, data_array)
    discrete_cmap = utilities.cmap_discretize(cm.RdYlBu_r, 10)
    m.imshow(data_array, extent=extent, origin='lower', interpolation='none',
             cmap=discrete_cmap)
    plt.colorbar(orientation='horizontal', fraction=0.056, pad=0.04)
    plt.tight_layout()
    plt.savefig('plume_speed_map.png', bbox_inches='tight')

    plt.close()

def plot_emission_direction_map(plume_archive, lats, lons, min_lat=None,
                            min_lon=None,
                    max_lat=None, max_lon=None):
    # First we need a grid on which to map all our speeds
    # This can come from lats and lons
    # Second, you set that grid to a whole pile of nans
    # Then you write all your speed values to a nice dictionary
    # Then when you have your dictionary, average each entry
    # Then for each latlon tuple, bin it to the nearest latlon on the grid
    # and write that to the array of nans
    # Then you do like imshow or something, but you'd have to get regular
    # latlons probably, unless there's some other grid plotting method which
    # allows regular latlons

    plt.close()

    fig, ax = plt.subplots()

    latlon_dictionary = {}
    lats = np.linspace(np.min(lats), np.max(lats), 100)
    lons = np.linspace(np.min(lons), np.max(lons), 100)

    lons, lats = np.meshgrid(lons, lats)

    # Get your zip grid of latlon coordinates
    latlon_zip = np.array(zip(lats.ravel(), lons.ravel()),
                          dtype=('f8,'
                          'f8')).reshape(
        lats.shape)
    latlon_zip = latlon_zip.flatten()

    for i in plume_archive:
        # If the plume has merged, the plume source is found in the
        # pre-merge track
        if plume_archive[i].merged == True:
            if len(plume_archive[i].pre_merge_track_centroid_direction) < 4:
                continue
            emission_lat = plume_archive[i].pre_merge_track_centroid_lat[0]
            emission_lon = plume_archive[i].pre_merge_track_centroid_lon[0]
            emission_directions = plume_archive[
                i].pre_merge_track_centroid_direction[
                3]
        else:
            if len(plume_archive[i].track_centroid_direction) < 4:
                continue
            emission_lat = plume_archive[i].track_centroid_lat[0]
            emission_lon = plume_archive[i].track_centroid_lon[0]
            emission_directions = plume_archive[i].track_centroid_direction[3]

        nearest = min(latlon_zip, key=lambda x: utilities.
                      haversine(x[1],
                                x[0],
                                emission_lon,
                                emission_lat))

        # If this entry already exists in the dictionary we just append to it
        if (nearest[0], nearest[1]) in latlon_dictionary:
            latlon_dictionary[(nearest[0], nearest[1])].append(np.nanmean(
                emission_directions))
        else:
            latlon_dictionary[(nearest[0], nearest[1])] = []
            latlon_dictionary[(nearest[0], nearest[1])].append(
                np.nanmean(emission_directions))

    # Average together lists for each entry in the dictionary
    for i in latlon_dictionary:
        latlon_dictionary[i] = np.nanmean(latlon_dictionary[i])

    # Data array
    data_array = np.zeros((lons.shape))
    data_array[:] = -999

    # Find the nearest latlon for each dictionary key, and extract the index
    # for it

    for i in latlon_dictionary:
        lat_bool = lats == i[0]
        lon_bool = lons == i[1]
        union_bool = lat_bool & lon_bool
        data_array[union_bool] = latlon_dictionary[i]

    data_array[data_array == -999] = np.nan

    extent = (np.min(lons), np.max(lons), np.min(lats), np.max(lats))
    m = Basemap(projection='cyl', llcrnrlon=extent[0], urcrnrlon=extent[1],
                llcrnrlat=extent[2], urcrnrlat=extent[3],
                resolution='i')

    m.drawcoastlines(linewidth=0.5)
    m.drawcountries(linewidth=0.5)

    # m.contourf(lons, lats, data_array)
    discrete_cmap = utilities.cmap_discretize(cm.RdYlBu_r, 10)
    m.imshow(data_array, extent=extent, origin='lower', interpolation='none',
             cmap=discrete_cmap, vmax=np.max(data_array[
                                                 np.isfinite(data_array)]) +
                                      40)

    plt.colorbar(orientation='horizontal', fraction=0.056, pad=0.04)
    plt.tight_layout()
    plt.savefig('plume_direction_map.png', bbox_inches='tight')

    plt.close()

def plot_emission_count_map(plume_archive, lats, lons, res=100, min_lat=None,
                            min_lon=None,
                    max_lat=None, max_lon=None):
    """
    Plots a map of emission count, akin to Ashpole and Washington (2012)
    key figure (to be updated to frequency)
    :param plume_archive:
    :param lats:
    :param lons:
    :param min_lat:
    :param min_lon:
    :param max_lat:
    :param max_lon:
    :return:
    """

    plt.close()

    fig, ax = plt.subplots()

    latlon_dictionary = {}

    lats = np.linspace(np.min(lats), np.max(lats), res)
    lons = np.linspace(np.min(lons), np.max(lons), res)
    lons, lats = np.meshgrid(lons, lats)

    # Get your zip grid of latlon coordinates
    latlon_zip = np.array(zip(lats.ravel(), lons.ravel()), dtype=('f8,'
                                                                  'f8')). \
        reshape(
        lats.shape)

    latlon_zip = latlon_zip.flatten()

    for i in plume_archive:
        # If the plume has merged, the plume source is found in the
        # pre-merge track
        if plume_archive[i].merged == True:
            emission_lat = plume_archive[i].pre_merge_track_centroid_lat[0]
            emission_lon = plume_archive[i].pre_merge_track_centroid_lon[0]
        else:
            emission_lat = plume_archive[i].track_centroid_lat[0]
            emission_lon = plume_archive[i].track_centroid_lon[0]

        nearest = min(latlon_zip, key=lambda x: utilities.
                      haversine(x[1],
                      x[0],
                      emission_lon,
                      emission_lat))

        # If this entry already exists in the dictionary we just add to it
        if (nearest[0], nearest[1]) in latlon_dictionary:
            latlon_dictionary[(nearest[0], nearest[1])] += 1
        else:
            latlon_dictionary[(nearest[0], nearest[1])] = 0
            latlon_dictionary[(nearest[0], nearest[1])] += 1

    print latlon_dictionary

    # Data array
    data_array = np.zeros((lons.shape))

    # Find the nearest latlon for each dictionary key, and extract the index
    # for it

    for i in latlon_dictionary:
        lat_bool = lats == i[0]
        lon_bool = lons == i[1]
        union_bool = lat_bool & lon_bool
        data_array[union_bool] = latlon_dictionary[i]

    data_array[data_array == 0] = np.nan

    extent = (np.min(lons), np.max(lons), np.min(lats), np.max(lats))
    m = Basemap(projection='cyl', llcrnrlon=extent[0], urcrnrlon=extent[1],
                llcrnrlat=extent[2], urcrnrlat=extent[3],
                resolution='i')

    m.drawcoastlines(linewidth=0.5)
    m.drawcountries(linewidth=0.5)

    # m.contourf(lons, lats, data_array)
    discrete_cmap = utilities.cmap_discretize(cm.RdYlBu_r, 10)
    m.imshow(data_array, extent=extent, origin='lower', interpolation='none',
             cmap=discrete_cmap)

    plt.colorbar(orientation='horizontal', fraction=0.056, pad=0.04)
    plt.tight_layout()
    plt.savefig('plume_count_map.png', bbox_inches='tight')

    plt.close()

def plot_plume_tracks(plume_archive, lats, lons, min_lat=None,
                            min_lon=None,
                    max_lat=None, max_lon=None):
    # First we need a grid on which to map all our speeds
    # This can come from lats and lons
    # Second, you set that grid to a whole pile of nans
    # Then you write all your speed values to a nice dictionary
    # Then when you have your dictionary, average each entry
    # Then for each latlon tuple, bin it to the nearest latlon on the grid
    # and write that to the array of nans
    # Then you do like imshow or something, but you'd have to get regular
    # latlons probably, unless there's some other grid plotting method which
    # allows regular latlons

    plt.close()

    fig, ax = plt.subplots()

    latlon_dictionary = {}

    extent = (np.min(lons), np.max(lons), np.min(lats), np.max(lats))

    m = Basemap(projection='cyl', llcrnrlon=extent[0], urcrnrlon=extent[1],
                llcrnrlat=extent[2], urcrnrlat=extent[3],
                resolution='i')

    m.drawcoastlines(linewidth=0.5)
    m.drawcountries(linewidth=0.5)

    for i in plume_archive:
        if plume_archive[i].merged == True:
            track_lats = plume_archive[i].pre_merge_track_centroid_lat
            track_lons = plume_archive[i].pre_merge_track_centroid_lon
        else:
            track_lats = plume_archive[i].track_centroid_lat
            track_lons = plume_archive[i].track_centroid_lon

        track_x, track_y = m(track_lons, track_lats)
        m.plot(track_x, track_y)

    plt.savefig('plume_tracks.png')

    plt.close()

def plot_traces(traces, retain=1000):
    """
    Convenience function:
    Plot traces with overlaid means and values
    """

    ax = pm.traceplot(traces[-retain:],
                      figsize=(12, len(traces.varnames) * 1.5),
                      lines={k: v['mean'] for k, v in
                             pm.df_summary(traces[-retain:]).
                      iterrows()})

    for i, mn in enumerate(pm.df_summary(traces[-retain:])['mean']):
        ax[i, 0].annotate('{:.2f}'.format(mn), xy=(mn, 0), xycoords='data'
                          , xytext=(5, 10), textcoords='offset points',
                          rotation=90
                          , va='bottom', fontsize='large', color='#AA0022')

# Assign it to that index in the empty data array of nans, and Bob's
# your uncle



