import numpy as np
from scipy.ndimage import measurements
from scipy import ndimage as ndi

"""
Handling of plumes objects, including updating of attributes and testing for
mechanism type
"""

# Global function to scan the SDFs for unique plumes
def scan_for_plumes(sdf_now, sdf_prev):
    """
    Scans a set of SDFs for plumes and labels them
    :param SDF_now:
    :return:
    """

    if sdf_prev is None:
        label_objects, nb_labels = ndi.label(sdf_now)
        sizes = np.bincount(label_objects.ravel())

        # Set clusters smaller than size 250 to zero
        mask_sizes = sizes > 250
        mask_sizes[0] = 0
        sdf_now = mask_sizes[label_objects]

        sdf_clusters, num = measurements.label(sdf_now)
        plume_ids = np.unique(sdf_clusters)
        large_plume_ids = np.unique(sdf_clusters[sdf_clusters != 0])
        new_ids = large_plume_ids
    else:
        label_objects, nb_labels = ndi.label(sdf_now)
        sizes = np.bincount(label_objects.ravel())

        # Set clusters smaller than size 250 to zero
        mask_sizes = sizes > 250
        mask_sizes[0] = 0
        sdf_now = mask_sizes[label_objects]

        sdf_clusters, num = measurements.label(sdf_now)
        plume_ids = np.unique(sdf_clusters)
        large_plume_ids = np.unique(sdf_clusters[sdf_clusters != 0])

        # Increase the plume_ID so that they are all new
        old_id_max = np.max(sdf_prev)
        sdf_clusters[sdf_clusters != 0] += old_id_max

        # Get an array of plumes which are overlapping
        overlaps = (sdf_clusters > 0) & (sdf_prev > 0)
        overlapping_ids = np.unique(sdf_clusters[overlaps])
        large_plume_ids = np.unique(sdf_clusters[sdf_clusters != 0])
        new_ids = [j for j in large_plume_ids if j not in overlapping_ids]
        if new_ids == [0]:
            new_ids = []
        else:
            new_ids = np.asarray(new_ids)
            new_id_bool = np.asarray([i != 0 for i in new_ids])
            if new_ids.shape[0] > 0:
                new_ids = np.unique(new_ids[new_id_bool])

        old_id_array = np.zeros(np.shape(sdf_clusters))
        for i in overlapping_ids:
            prev_ids = sdf_prev[sdf_clusters == i]
            # Take the most common of the previous IDs as the one which should
            # be applied to the new plume
            counts = np.bincount(prev_ids)
            prev_id = np.argmax(counts)
            # Set prev_ID to whatever that is
            sdf_clusters[sdf_clusters == i] = prev_id
            old_id_array[sdf_clusters == i] = prev_id
        large_plume_ids = np.unique(sdf_clusters[sdf_clusters != 0])
    return sdf_clusters, new_ids, large_plume_ids

# This returns a set of labeled plumes
# We then loop through those at the beginning and make objects
# You're doing label at every timestep, and where there is a plume which has
#  absolutely no overlap with a previous plume, you're calling a Plume object
# Then, for each existing plume, you're updating it. If the plume is gone,
# it dies (as determined by update_position)

# Ok so you start with a set of clusters. Loop through the clusters and get
# a whole pile of objects

# Then in the next timestep, you do the image processing again. Then loop
# through each individual plume which is already active and call the update
# methods (which may involve killing the plume)

# But you need a global overlap function to extract new IDs, then in the
# next timestep you will just get objects for each new plume

## Class of plumes objects
class Plume:

    plume_id = 0
    area = None
    centroid_lat = None
    centroid_lon = None
    source_lat = None
    source_lon = None
    centroid_speed = None
    centroid_direction = 0
    duration = 0
    emission_time = 0
    major_axis_position = 0
    minor_axis_position = 0
    LLJ_likelihood = 0
    CPO_likelihood = 0

    def __init__(self, plume_id, emission_time):
        self.plume_id = plume_id
        self.emission_time = emission_time

        # Now this function needs to work out for us where the new plume is
        # located by using the SDF map and finding where there is an ID
        # equal to the plume ID, then working out the centroid
        # Here we need a way to get the centroid for this object




    # So, at every timestep you have a whole load of SDFs
    # Then, when a new pixel becomes 1, an instance of this class is called
    # Can get a lat and lon for that, pass it to the object
    # Calling the update functions gets a new snapshot of the SDF for the
    # Sahara, looks for the overlapping plume and picks out the centroid. If
    #  none exists, the die method is called

    # So at every timestep you get the plumes and label them, then get
    # instances

    def update_position(self, lats, lons, sdf_now, plume_id):
        """
        Takes an SDF map and updates the centroid lat and lon
        :return:
        """
        valid_lons,


    def update_duration(self):
        """
        Takes an SDF map and updates the centroid lat and lon
        :return:
        """

    def update_speed(self):
        """
        Takes an SDF map and updates the centroid lat and lon
        :return:
        """

    def update_direction(self):
        """
        Takes an SDF map and updates the centroid lat and lon
        :return:
        """



    def move(self):
        pass

    def merge(self):
        pass

    def die(self):
        pass

    def update_position(self):
        pass

    def update_duration(self):
        pass

    def update_speed(self):
        pass

    def update_direction(self):
        pass

    def update_axes(self):
        pass

    def update_mechanism_likelihood(self):
        pass


    # Ok so there could be a method to update various parameters
    # Then call these each time an object instance is created

    # Ok so if you make an __init__(self) function, python will run it every
    #  time the class is invoked.

# Area attribute
# Leading edge speed attribute
# Centroid speed attribute
# Centroid position attribute
# Source latitude attribute
# Source longitude attribute
# Duration attribute
# Emission timing attribute
# Major axis attribute
# Minor axis attribute
# LLJ likelihood attribute
# CPO likelihood attribute
# Other likelihood attribute

# Move method
# Merge method
# Die method
# Implement LLJ checks method
# Implement CPO checks method

## Class of convection objects
# Area attribute
# Centroid position attribute

# Move method
# Merge method
# Die method

# Functions to cloud screen, generate SDFs and categorise to get instances
# of objects, so main should only have to run these at each timestep,
# and you'll get a whole bunch of objects generated within this function
# which you don't see in main

# Functions to plot a snapshot of all the plumes at a given moment in time (
# this comes from an aggregate of plumes - doesn't make too much sense as a
# part of a class)
