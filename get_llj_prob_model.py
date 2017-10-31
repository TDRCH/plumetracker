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
        mean_direction = np.nanmean(plume.track_centroid_direction)
        mean_primary_axis_direction = np.nanmean(
            plume.track_primary_axis_direction)
        axis_direction_offsets.append(
            abs(mean_direction-mean_primary_axis_direction))
        durations.append(plume.duration.total_seconds())
        emission_speeds.append(plume.track_speed_centroid[1])
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
    data['emission_speed'] = emission_speeds
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

if __name__ == '__main__':

    data = create_plume_dataframe('LLJ_manual_ID.csv', 'plume_archive_LLJ')

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
