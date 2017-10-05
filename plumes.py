import numpy as np
from scipy.ndimage import measurements
from scipy import ndimage as ndi
from skimage.measure import label, regionprops
import datetime

import utilities

"""
Handling of plumes objects, including updating of attributes and testing for
mechanism type
"""

# Global function to scan the SDFs for unique plumes
def scan_for_plumes(sdf_now, sdf_prev, used_ids):
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
        print large_plume_ids
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
            old_id_array[sdf_clusters == i] = prev_id
            sdf_clusters[sdf_clusters == i] = prev_id
        large_plume_ids = np.unique(sdf_clusters[sdf_clusters != 0])

        # Fix the new IDs so that they follow directly from the previous
        # largest one
        if len(new_ids) > 0:
            old_ids = np.unique(sdf_prev[sdf_prev != 0])
            replacement_new_ids = np.arange(np.max(used_ids)+1, np.max(
                used_ids)+1+len(new_ids))
            for i in np.arange(0, len(new_ids)):
                sdf_clusters[sdf_clusters == new_ids[i]] = (
                    replacement_new_ids[i])
                large_plume_ids[large_plume_ids == new_ids[i]] = (
                replacement_new_ids[i])
                new_ids[i] = replacement_new_ids[i]

    return sdf_clusters, new_ids, large_plume_ids

class Plume:

    """
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
    track_lat = []
    track_lon = []
    major_axis_position = 0
    minor_axis_position = 0
    LLJ_likelihood = 0
    CPO_likelihood = 0
    """
    def __init__(self, plume_id, emission_time):
        self.plume_id = plume_id
        self.emission_time = emission_time
        # Attributes to track the centroid position through time
        self.track_centroid_lat = []
        self.track_centroid_lon = []
        # The previous two timesteps are recorded
        self.track_lons = []
        self.track_lats = []

    # So, at every timestep you have a whole load of SDFs
    # Then, when a new pixel becomes 1, an instance of this class is called
    # Can get a lat and lon for that, pass it to the object
    # Calling the update functions gets a new snapshot of the SDF for the
    # Sahara, looks for the overlapping plume and picks out the centroid. If
    #  none exists, the die method is called

    # So at every timestep you get the plumes and label them, then get
    # instances

    def update_position(self, lats, lons, sdf_plumes, plume_id):
        """
        Updates the position of the plume, including its centroid and point
        coordinates, archiving them in a track
        :return:
        """

        plume_bool = sdf_plumes == plume_id

        plume_lons = lons[plume_bool]
        plume_lats = lats[plume_bool]

        # Calculate the centroid
        sum_x = np.sum(plume_lons)
        sum_y = np.sum(plume_lats)

        centroid_lon = sum_x/plume_lons.shape[0]
        centroid_lat = sum_y/plume_lons.shape[0]

        self.centroid_lon = centroid_lon
        self.centroid_lat = centroid_lat

        self.track_centroid_lat.append(centroid_lat)
        self.track_centroid_lon.append(centroid_lon)

        self.plume_lons = plume_lons
        self.plume_lats = plume_lats

        self.track_lons.append(plume_lons)
        self.track_lats.append(plume_lats)

        # Only record the previous three timesteps (including this one)
        if len(self.track_lons) > 3:
            del self.track_lons[0]

        if len(self.track_lats) > 3:
            del self.track_lats[0]

        print 'Centroid lat', centroid_lat
        print 'Centroid lon', centroid_lon

        print '\n'

    def update_leading_edge(self, sdf_plumes):
        """
        Updates the identification of the pixels associated with a leading
        edge on the plume
        :return:
        """

        # Ok, very simple, just compare it to the last centroid

        # Assuming this function is run AFTER the update_position,
        # the previous centroid is two steps back in the track

        if self.duration < datetime.timedelta(hours=0.75):
            self.leading_edge_lon = None
            self.leading_edge_lat = None
        else:
            previous_centroid_lat = self.track_centroid_lat[-2]
            previous_centroid_lon = self.track_centroid_lon[-2]

            plume_bool = sdf_plumes == self.plume_id

            # Calculate distances of all coordinates from the previous centroid
            distances = np.zeros((len(self.plume_lons)))
            distances = np.asarray([utilities.haversine(self.plume_lons[j],
                                                        self.plume_lats[j],
                                                        previous_centroid_lon,
                                                        previous_centroid_lat)
                                    for j in np.arange(0,
                                                       len(self.plume_lons))])

            furthest_index = np.where(distances == np.max(distances))[0][0]

            self.leading_edge_lon = self.plume_lons[furthest_index]
            self.leading_edge_lat = self.plume_lats[furthest_index]

    def update_leading_edge_2(self, sdf_plumes):
        """
        Attempt #2 looks for the edge point with the greatest distance from
        the nearest edge of the previous plume
        :param sdf_plumes:
        :return:
        """

        # Identify the edge of the plume by scanning for an edge point then
        # checking its neighbours

        plume_bool = sdf_plumes == self.plume_id
        for idx, val in enumerate(plume_bool):
            # Keep going until we find a value
            if val:
                # Check if there is an adjacent zero
                left_x = idx[0]-1
                left_y = idx[1]
                above_x = idx[0]
                above_y = idx[1]+1
                right_x = idx[0]+1
                right_y = idx[1]
                below_x = idx[0]
                below_y = idx[1]-1

                adjacent_values = [plume_bool[left_x, left_y], plume_bool[
                    above_x, above_y], plume_bool[right_x, right_y],
                                   plume_bool[below_x, below_y]]

                # Now need to check if either False or 0 is in these
                # adjacent values. If it is, then we know we're at an edge
                # point. In that case we need to initiate a while loop which
                #  searches for another edge point
                # Otherwise, you need to continue the loop until we do find
                # an edge point
                # The while condition is that our current coordinates are
                # not equal to the first edge point which you do find




    def update_bbox(self):
        """
        Updates the bounding box coordinates associated with the plume
        :param valid_lons:
        :param valid_lats:
        :return:
        """

        self.bbox_leftlon = np.min(self.plume_lons)
        self.bbox_rightlon = np.max(self.plume_lons)

        self.bbox_bottomlat = np.min(self.plume_lats)
        self.bbox_toplat = np.max(self.plume_lats)

        print self.bbox_leftlon
        print self.bbox_rightlon
        print self.bbox_bottomlat
        print self.bbox_toplat
        print '\n'

    def update_majorminor_axes(self):
        """
        Updates the major and minor axes of a plume with PCA
        :return:
        """

        # Subtract the mean from each dimension
        lons = self.plume_lons - np.mean(self.plume_lons)
        lats = self.plume_lats - np.mean(self.plume_lats)

        # Calculate the covariance matrix of the coordinates
        coords = np.vstack([lons, lats])
        cov = np.cov(coords)
        evals, evecs = np.linalg.eig(cov)

        # Sort eigenvalues
        sort_indices = np.argsort(evals)[::-1]
        evec1, evec2 = evecs[:, sort_indices]

        # Add the mean back in
        self.lon_major = evec1[0]+np.mean(self.plume_lons)
        self.lat_major = evec1[1]+np.mean(self.plume_lats)
        self.lon_minor = evec2[0]+np.mean(self.plume_lons)
        self.lat_minor = evec2[1]+np.mean(self.plume_lats)

        print self.lon_major
        print self.lat_major
        print self.lon_minor
        print self.lat_minor
        print '\n'

    def update_duration(self, date):
        """
        Updates the duration of the dust plume with the time since its emission
        :return:
        """

        self.duration = date - self.emission_time

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
