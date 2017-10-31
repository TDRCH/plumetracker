import matplotlib
matplotlib.use('Agg')
import pandas as pd
import numpy as np
import shelve
import seaborn
import matplotlib.pyplot as plt

"""
Assigns detected attributes to manually detected plumes, creating a data
table for the events
"""

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

    data['conv_distance'] = conv_distances
    data['axis_direction_offset'] = axis_direction_offsets
    data['duration'] = durations
    data['emission_speed'] = emission_speeds

    # Remove the plume ID field
    del data['plume_ID']

    # For each identified ID, create a row in a pandas data frame and
    # extract the equivalent attributes from the plume archive
    # Save the data frame as a csv
    return data

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