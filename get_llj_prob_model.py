import matplotlib
matplotlib.use('Agg')
import pandas as pd
import numpy as np
import pymc3 as pm
import matplotlib.pyplot as plt
import seaborn
import warnings
warnings.filterwarnings('ignore')
from collections import OrderedDict
from time import time
from scipy.optimize import fmin_powell
from scipy import integrate
import theano as thno
import theano.tensor as T
import plotting
import datetime as dt
import shelve
from netCDF4 import Dataset
import math
from skimage.measure import label, regionprops

import utilities

def create_plume_dataframe(manual_filename, plume_archive_filename):
    """
    Creates a pandas dataframe from manually identified dust plumes,
    assigning them attributes found by plumetracker
    :param manual_filename:
    :param plume_archive_filename:
    :return:
    """

    # Read in the csv
    data = pd.read_csv('/ouce-home/students/hert4173/workspace/plumetracker/'
                       +manual_filename,
        header=None, names=['plume_ID', 'LLJ'])

    # Read in the plume archive
    plume_archive = shelve.open(plume_archive_filename)

    conv_distances = []
    axis_direction_offsets = []
    durations = []
    emission_speeds = []
    distances_from_09 = []
    max_extents = []

    for index, row in data.iterrows():
        plume = plume_archive[str(row['plume_ID'])]
        conv_distances.append(plume.conv_distance)
        #axis_direction_offsets.append(plume.mean_axis_offset)
        durations.append(plume.duration.total_seconds())
        emission_speeds.append(plume.track_speed_centroid[1])

        greatest_extent_ind = np.where(plume.track_area == np.max(
            plume.track_area))[0][0]

        """
        if greatest_extent_ind > 0 and greatest_extent_ind < len(
                plume.track_centroid_lat)-1:
            # Extract centroid positions around the greatest extent time
            first_position_lat = plume.track_centroid_lat[
                greatest_extent_ind-1]
            first_position_lon = plume.track_centroid_lon[
                greatest_extent_ind-1]
            last_position_lat = plume.track_centroid_lat[
                greatest_extent_ind+1]
            last_position_lon = plume.track_centroid_lon[
                greatest_extent_ind+1]
        """

        # Get the direction from the total distance travelled by the plume
        first_position_lat = plume.track_centroid_lat[0]
        first_position_lon = plume.track_centroid_lon[0]
        last_position_lat = plume.track_centroid_lat[-1]
        last_position_lon = plume.track_centroid_lon[-1]

        # Get the direction
        final_direction = utilities.calculate_initial_compass_bearing((
            first_position_lat, first_position_lon),(last_position_lat,
            last_position_lon))

        # Do the same as before for the offset
        mean_primary_axis_direction = np.nanmean(
            plume.track_primary_axis_direction)

        #plt.close()

        #plt.contourf(lons, lats, plume.track_plume_bool[
        #    greatest_extent_ind].toarray())

        image = plume.track_plume_bool[
            greatest_extent_ind].toarray()
        label_img = label(image)
        regions = regionprops(label_img)

        for props in regions:
            orientation = props.orientation
            orientation = math.degrees(orientation)+90

        u = 1 * np.sin(math.radians(orientation))
        v = 1 * np.cos(math.radians(orientation))

        u1 = 1 * np.sin(math.radians(final_direction))
        v1 = 1 * np.cos(math.radians(final_direction))

        #plt.imshow(plume.track_plume_bool[greatest_extent_ind].toarray(),
        #           origin='lower')
        #plt.quiver(plume.track_centroid_lon[greatest_extent_ind],
        #           plume.track_centroid_lat[greatest_extent_ind], u, v,
        #           color='r')
        #plt.quiver(plume.track_centroid_lon[greatest_extent_ind],
        #           plume.track_centroid_lat[greatest_extent_ind], u1, v1,
        #           color='g')

        #plt.savefig('Orientation_check_plot_'+str(plume.plume_id)+'.png')
        #plt.close()

        mean_axis_difference = abs(utilities.get_distance(
            final_direction, orientation))

        #print 'Direction travelled', final_direction
        #print 'Orientation of axis', orientation

        diff_180 = abs(mean_axis_difference - 180)
        diff_0 = abs(mean_axis_difference)

        mean_axis_offset = np.min([diff_180, diff_0])
        axis_direction_offsets.append(mean_axis_offset)

        #if str(row['plume_ID']) == '375':
        #    print 'Orientation', orientation
        #    print 'Final direction', final_direction

        #print mean_axis_offset

        emission_time = plume.dates_observed[0]

        # Take the smallest time to the nearest 0900UTC
        nine = emission_time.replace(hour=9, minute=0, second=0)
        nine_1 = nine+dt.timedelta(days=1)
        distance_nine = abs(emission_time-nine)
        distance_nine_1 = abs(nine_1-emission_time)
        distance_from_09 = np.min([distance_nine.total_seconds(),
                                   distance_nine_1.total_seconds()])
        distances_from_09.append(distance_from_09)
        max_extents.append(np.max(plume.track_area))

    data['conv_distance'] = conv_distances
    data['axis_direction_offset'] = axis_direction_offsets
    data['duration'] = durations
    #data['emission_speed'] = emission_speeds
    data['time_to_09'] = distances_from_09
    data['max_extent'] = max_extents

    # Remove the plume ID field
    del data['plume_ID']

    # Close shelf
    plume_archive.close()

    # For each identified ID, create a row in a pandas data frame and
    # extract the equivalent attributes from the plume archive
    # Save the data frame as a csv
    return data

def logistic_llj_model(data, dpdtvarname, idpdtvarname1, idpdtvarname2,
                       idpdtvarname3):
    """
    Samples the probability space for parameters in a multivariate logistic
    regression model
    :param data:
    :param dpdtvarname:
    :param idpdtvarname1:
    :param idpdtvarname2:
    :param idpdtvarname3:
    :return:
    """

    with pm.Model() as logistic_model:
        # String is R syntax: want income as a function of the following
        # variables
        pm.glm.GLM.from_formula(dpdtvarname +' ~ '+ idpdtvarname1 + ' + ' +
                                idpdtvarname2 + ' + ' + idpdtvarname3,
                                data, family = pm.glm.families.Binomial())
        trace_logistic_model = pm.sample(4000)
        trace = trace_logistic_model[1000:]

    return trace

def plot_comparison_histograms(data):
    """
    Plots histograms comparing the distributions of independent variables
    for LLJ and CPO cases
    :param data:
    :return:
    """

    llj_data = data.loc[data['LLJ']==1]
    cpo_data = data.loc[data['LLJ']==0]

    independent_variables = ['conv_distance', 'time_to_09',
                             'axis_direction_offset', 'duration']
    labels = ['Distance from deep convection (km)', 'Time away from 0900UTC '
              '(seconds)', 'Direction offset from major axis (degrees)',
              'Duration (seconds)']
    bins = [np.arange(0, 1100, 100), np.arange(0, 50000, 5000), np.arange(
        0, 100, 10), np.arange(0, 800000, 100000)]

    for i in np.arange(0, len(independent_variables)):
        plt.close()
        llj_data_var = llj_data[independent_variables[i]]
        cpo_data_var = cpo_data[independent_variables[i]]
        hist_bins = bins[i]
        xlabel = labels[i]
        llj_data_var.plot.hist(alpha=0.5, bins=hist_bins, label='Low-level '
                                                                'jet')
        cpo_data_var.plot.hist(alpha=0.5, bins=hist_bins, label='Cold pool '
                                                                'outflow')
        plt.xlabel(xlabel)
        plt.legend()
        plt.savefig('llj_vs_cpo_'+independent_variables[i]+'.png')
        plt.close()

if __name__ == '__main__':


    sdf_test = Dataset('/soge-home/data_not_backed_up/satellite/meteosat/'
                       'seviri/15-min/0.03x0.03/sdf/nc/JUNE2010/SDF_v2/'
                       'SDF_v2.201006031500.nc')
    lons, lats = np.meshgrid(sdf_test.variables['longitude'][:],
                             sdf_test.variables['latitude'][:])

    data = create_plume_dataframe('LLJ_manual_ID_2010.csv',
                                  '/soge-home/projects/seviri_dust/'
                                  'plumetracker/completed_archives'
                                  '/plume_archive_flicker_v3_2010', lons, lats)
    data = data.dropna(0)

    g = seaborn.pairplot(data)
    plt.savefig('LLJ_pairplot.png')
    plt.close()

    # Compute the correlation matrix
    corr = data.corr()

    # Generate a mask for the upper triangle
    mask = np.zeros_like(corr, dtype=np.bool)
    mask[np.triu_indices_from(mask)] = True

    # Set up the matplotlib figure
    f, ax = plt.subplots(figsize=(11, 9))

    # Generate a custom diverging colormap
    cmap = seaborn.diverging_palette(220, 10, as_cmap=True)

    # Draw the heatmap with the mask and correct aspect ratio
    seaborn.heatmap(corr, mask=mask, cmap=cmap, vmax=.3,
                    linewidths=.5, cbar_kws={"shrink": .5}, ax=ax)
    plt.savefig('LLJ_correlheatmap.png')
    plt.close()

    with pm.Model() as logistic_model:
        # String is R syntax: want income as a function of the following
        # variables
        pm.glm.GLM.from_formula('LLJ ~ conv_distance + time_to_09 + '
                                'duration',
                                data,
                                family = pm.glm.families.Binomial())
        trace_logistic_model = pm.sample(4000)

        plotting.plot_traces(trace_logistic_model, retain=1000)

    plt.figure(figsize=(9, 7))
    trace = trace_logistic_model[1000:]
    print np.nanmean(trace['conv_distance'])
    seaborn.jointplot(trace['conv_distance'], trace['time_to_09'], kind="hex",
                      color="#4CB391")
    plt.xlabel("beta_conv_distance")
    plt.ylabel("beta_time_to_09")
    plt.savefig('LLJ_probability_jointplot.png')
    plt.close()

    prob = np.nanmean(1 / (1 + np.exp(-(trace['Intercept'] +
                              trace['conv_distance'] * 500 +
                              trace['time_to_09'] * 80 +
                              trace['duration'] * 10000))))

    # This should be a distribution of values across our parameter space

    plot_comparison_histograms(data)