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

import pinkdust

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

# Function to plot active plume number through time
def plot_plume_count(plume_archive):
    """
    Plots the number of active plumes at each timestep in the plume archive
    :param plume_archive:
    :return:
    """

    plt.close()
    fig, ax = plt.subplots()

   # Extract all the dates in use in the plume archive
    all_dates = []
    for i in plume_archive:
        dates = plume_archive[i].dates_observed
        for j in dates:
            all_dates.append(j)

    unique_dates = np.unique(all_dates)
    unique_dates = np.sort(unique_dates)
    date_count = []

    for i in unique_dates:
        date_count.append(all_dates.count(i))

    ax.plot(unique_dates, date_count, color='black')

    plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%m/%d/%Y'))
    plt.gca().xaxis.set_major_locator(mdates.DayLocator())
    plt.gca().xaxis.set_minor_locator(mdates.HourLocator())

    #plt.gcf().autofmt_xdate()

    ax.set_ylabel('Plume count')
    ax.set_xlabel('Time')
    ax.set_ylim(0, 25)
    plt.grid()

    plt.savefig('plume_count_through_time_3.png')
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

    plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%m/%d/%Y'))
    plt.gca().xaxis.set_major_locator(mdates.DayLocator())
    plt.gca().xaxis.set_minor_locator(mdates.HourLocator())

    #plt.gcf().autofmt_xdate()

    ax.set_ylabel('Mean plume area')
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

    plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%m/%d/%Y'))
    plt.gca().xaxis.set_major_locator(mdates.DayLocator())
    plt.gca().xaxis.set_minor_locator(mdates.HourLocator())

    #plt.gcf().autofmt_xdate()

    ax.set_ylabel('Total plume area')
    ax.set_xlabel('Time')
    plt.grid()

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
            speeds = np.asarray(plume_archive[i].track_area)
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
        averaged_area = np.nanmean(speed_list)
        date_dictionary[i] = averaged_area

    # Plot
    lists = sorted(date_dictionary.items())
    x, y = zip(*lists)

    ax.plot(x, y, color='black')

    plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%m/%d/%Y'))
    plt.gca().xaxis.set_major_locator(mdates.DayLocator())
    plt.gca().xaxis.set_minor_locator(mdates.HourLocator())

    # plt.gcf().autofmt_xdate()

    ax.set_ylabel('Plume speed (m/s)')
    ax.set_xlabel('Time')
    plt.grid()

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

    plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%m/%d/%Y'))
    plt.gca().xaxis.set_major_locator(mdates.DayLocator())
    plt.gca().xaxis.set_minor_locator(mdates.HourLocator())

    # plt.gcf().autofmt_xdate()

    ax.set_ylabel('Plume direction (deg)')
    ax.set_xlabel('Time')
    plt.grid()

    plt.savefig('plume_direction_through_time.png')
    plt.close()

# Function to plot mean plume size through time
# Function to plot plume size distribution by region
# Function to plot mean plume direction for given regions
# Function to plot source region density map
# Function to plot LLJ likelihood and CPO likelihood by region


