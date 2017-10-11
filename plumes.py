"""
Handling of plumes objects, including updating of attributes and testing for
mechanism type
"""

import numpy as np
from scipy.ndimage import measurements
from scipy import ndimage as ndi
from skimage.measure import label, regionprops
from skimage import feature
import datetime
import random
from copy import deepcopy
from shapely.geometry import shape
from shapely.ops import transform
import pyproj
import shapely.ops as ops
from shapely.geometry.polygon import Polygon
from functools import partial

import utilities

# Global function to scan the SDFs for unique plumes
def scan_for_plumes(sdf_now, sdf_prev, used_ids):
    """
    Scans a set of SDFs for plumes and labels them
    :param SDF_now:
    :return:
    """

    merge_ids = []

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
        #print large_plume_ids
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
        old_id_max = np.max(used_ids)
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
            # Remove zeros
            prev_ids = prev_ids[prev_ids != 0]
            # Take the most common of the previous IDs as the one which should
            # be applied to the new plume
            counts = np.bincount(prev_ids)
            prev_id = np.argmax(counts)
            # Set prev_ID to whatever that is
            old_id_array[sdf_clusters == i] = prev_id
            sdf_clusters[sdf_clusters == i] = prev_id
            if len(np.unique(prev_ids)) > 1:
                merge_ids.append(prev_id)

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

    return sdf_clusters, new_ids, large_plume_ids, merge_ids

class Plume:

    def __init__(self, plume_id, emission_time):
        self.plume_id = plume_id
        self.emission_time = emission_time
        # Attributes to track the centroid position through time
        self.track_centroid_lat = []
        self.track_centroid_lon = []
        # The previous two timesteps are recorded
        self.track_lons = []
        self.track_lats = []
        self.track_edges_lat = []
        self.track_edges_lon = []
        self.dates_observed = []
        self.track_area = []
        self.track_centroid_direction = []
        self.leading_edge_lon = None
        self.leading_edge_lat = None
        self.merged = False
        self.track_speed_centroid = []
        self.pre_merge_track_lons = []
        self.pre_merge_track_lats = []
        self.pre_merge_track_centroid_lat = []
        self.pre_merge_track_centroid_lon = []
        self.pre_merge_track_area = []
        self.pre_merge_track_speed_centroid = []
        self.pre_merge_dates_observed = []
        self.pre_merge_track_direction = []

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

        #print 'Centroid lat', centroid_lat
        #print 'Centroid lon', centroid_lon

        #print '\n'

    def update_area(self):

        cop = {"type": "Polygon", "coordinates": [zip(self.plume_lons,
                                                      self.plume_lats)]}
        s = shape(cop)

        geom = Polygon([(a, b) for a, b in zip(self.plume_lons,
                                               self.plume_lats)])
        geom_area = ops.transform(
            partial(
                pyproj.transform,
                pyproj.Proj(init='EPSG:4326'),
                pyproj.Proj(
                    proj='aea',
                    lat1=geom.bounds[1],
                    lat2=geom.bounds[3])),
            geom)

        # Store area in m^2
        self.area = geom_area.area
        self.track_area.append(geom_area.area)

    def update_leading_edge_4(self, sdf_plumes, lons, lats):
        """
        All those previous version were joke versions - this is the true
        version, the only true version
        :param sdf_plumes:
        :return:
        """

        plume_bool = sdf_plumes == self.plume_id
        edge_plume_bool = feature.canny(plume_bool, sigma=0.2)

        #print np.unique(plume_bool)
        #print np.unique(edge_plume_bool)

        edge_lons = lons[edge_plume_bool]
        edge_lats = lats[edge_plume_bool]

        # No leading edge if the duration of the plume is too low
        if self.duration < datetime.timedelta(hours=0.75):
            self.leading_edge_lon = None
            self.leading_edge_lat = None
            self.track_edges_lon.append(edge_lons)
            self.track_edges_lat.append(edge_lats)

        # Same applies for duration since a merge
        if self.merged == True:
            duration_since_merge = self.dates_observed[-1] - self.merge_date
            if duration_since_merge < datetime.timedelta(hours=0.75):
                self.leading_edge_lon = None
                self.leading_edge_lat = None
                self.track_edges_lon.append(edge_lons)
                self.track_edges_lat.append(edge_lats)

        else:
            previous_edge_lons = self.track_edges_lon[-1]
            previous_edge_lats = self.track_edges_lat[-1]

            # We now need the edge with the greatest distance from its
            # nearest neighbour

            # So for each edge pixel, calculate the distance to all previous
            # edge pixels
            edge_distances = []

            for i in np.arange(0, len(edge_lats)):
                distances = []
                for j in np.arange(0, len(previous_edge_lats)):
                    distances.append(utilities.haversine(edge_lons[i],
                                                        edge_lats[
                        i], previous_edge_lons[j], previous_edge_lats[j]))
                edge_distances.append(np.min(distances))

            #print np.where(edge_distances == np.max(edge_distances))[0]
            greatest_edge_distance_indices = np.where(edge_distances ==
                                                      np.max(
                                                          edge_distances))[0]

            # Take the 5 largest distances as the leading edge
            sorted_distances = np.asarray(deepcopy(edge_distances))
            sorted_distance_indices = sorted_distances.argsort()[-5:][::-1]

            self.leading_edge_lon = [edge_lons[j]
                                for j in sorted_distance_indices]

            self.leading_edge_lat = [edge_lats[j]
                                for j in sorted_distance_indices]

            self.track_edges_lon.append(edge_lons)
            self.track_edges_lat.append(edge_lats)

            # Remove edges from the track which are older than 2 timesteps
            if len(self.track_edges_lon) > 3:
                del self.track_edges_lon[0]

            if len(self.track_edges_lat) > 3:
                del self.track_edges_lat[0]

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

        #print self.bbox_leftlon
        #print self.bbox_rightlon
        #print self.bbox_bottomlat
        #print self.bbox_toplat
        #print '\n'

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

    def update_duration(self, date):
        """
        Updates the duration of the dust plume with the time since its emission
        :return:
        """

        self.duration = date - self.emission_time
        self.dates_observed.append(date)

    def update_max_extent(self):
        """
        If the plume is at its greatest extent so far, updates the maximum
        extent
        :return:
        """

        # Note: if this is too costly, just store the date and use that to
        # index the lats and lons later
        if len(self.track_area[:]) <= 1:
            self.maximum_extent_lons = self.plume_lons
            self.maximum_extent_lats = self.plume_lats
            self.maximum_extent_time = self.dates_observed[-1]

        elif self.area >= np.max(self.track_area[:-1]):
            #print 'True for ', self.plume_id
            self.maximum_extent_lons = self.plume_lons
            self.maximum_extent_lats = self.plume_lats
            self.maximum_extent_time = self.dates_observed[-1]

    def update_GPE_speed(self):
        """
        Updates the greatest plume extent speed of the plume at the end of
        its lifecycle
        :return:
        """

        # Distance of each pixel in the maximum extent from the source region
        plume_distances = [utilities.haversine(self.maximum_extent_lons[j],
                                               self.maximum_extent_lats[j],
                                               self.track_centroid_lon[0],
                                               self.track_centroid_lat[0])
                           for j in
                           np.arange(0, len(self.maximum_extent_lons))]

        greatest_plume_distance = np.max(plume_distances)

        date_of_max_extent = self.maximum_extent_time
        secs_to_max_extent = (date_of_max_extent-self.dates_observed[
            0]).seconds
        if secs_to_max_extent == 0:
            self.speed_gpe = np.nan
        else:
            self.speed_gpe = (greatest_plume_distance*1000)/secs_to_max_extent
        #print 'Plume', self.plume_id
        #print 'GPE speed', self.speed_gpe

    def update_centroid_speed(self):
        """
        Updates the centroid speed of the plume and its track. This must be
        run after updating the duration.
        :return:
        """

        # Only get a centroid speed if we have a previous centroid and we
        # haven't merged recently
        if self.duration <= datetime.timedelta(hours=0.5):
            self.speed_centroid = np.nan
            self.track_speed_centroid.append(self.speed_centroid)
        elif self.merged == True:
            if (self.dates_observed[-1] - self.merge_date) <= \
                    datetime.timedelta(hours=0.5):
                self.speed_centroid = np.nan
                self.track_speed_centroid.append(self.speed_centroid)

        else:
            centroid_distance = utilities.haversine(self.centroid_lon,
                                                    self.centroid_lat,
                                                    self.track_centroid_lon[
                                                        -2],
                                                    self.track_centroid_lat[
                                                        -2])
            secs_since_previous = (self.dates_observed[-1] -
                                   self.dates_observed[
                -2]).seconds

            self.speed_centroid = (centroid_distance*1000/secs_since_previous)
            self.track_speed_centroid.append(self.speed_centroid)
        #print 'Plume', self.plume_id
        #print 'Centroid speed', self.speed_centroid
        #print 'Mean centroid speed', np.nanmean(self.track_speed_centroid)

    def merge(self):
        """
        Used when the plume has merged - moves all track data to a pre-merge
        track, so that velocities etc are not calculated across a merge
        :return:
        """

        self.merged = True
        self.merge_date = self.dates_observed[-1]

        self.pre_merge_track_lons = self.track_lons
        self.pre_merge_track_lats = self.track_lats
        self.pre_merge_track_centroid_lat = self.track_centroid_lat
        self.pre_merge_track_centroid_lon = self.track_centroid_lon
        self.pre_merge_track_area = self.track_area
        self.pre_merge_track_speed_centroid = self.track_speed_centroid
        self.pre_merge_dates_observed = self.dates_observed
        self.pre_merge_track_centroid_direction = self.track_centroid_direction

        self.track_lons = []
        self.track_lats = []
        self.track_centroid_lat = []
        self.track_centroid_lon = []
        self.track_area = []
        self.track_speed_centroid = []
        self.dates_observed = []
        self.track_centroid_direction = []

    def update_centroid_direction(self):
        """
        Identify the direction in which the centroid has travelled since the
        previous timestep
        :return:
        """

        # Only get a centroid direction if we have a previous centroid and
        # haven't merged recently
        if self.duration <= datetime.timedelta(hours=0.5):
            self.centroid_direction = np.nan
            self.track_centroid_direction.append(self.centroid_direction)
        elif self.merged == True:
            if (self.dates_observed[-1] - self.merge_date) <= \
                    datetime.timedelta(hours=0.5):
                self.centroid_direction = np.nan
                self.track_centroid_direction.append(self.centroid_direction)
        else:
            self.centroid_direction = utilities.\
                calculate_initial_compass_bearing((self.track_centroid_lat[-2],
                                                   self.track_centroid_lon[-2]),
                                                  (self.centroid_lat,
                                                   self.centroid_lon))

            #print self.centroid_direction

            self.track_centroid_direction.append(self.centroid_direction)

    def move(self):
        pass


    def update_axes(self):
        pass

    def update_mechanism_likelihood(self):
        pass


    # Ok so there could be a method to update various parameters
    # Then call these each time an object instance is created

    # Ok so if you make an __init__(self) function, python will run it every
    #  time the class is invoked.

class Convection:

    def __init__(self, cloud_id, initiation_time):
        self.cloud_id = cloud_id
        self.initiation_time = initiation_time
        # Attributes to track the centroid position through time
        self.track_centroid_lat = []
        self.track_centroid_lon = []
        # The previous two timesteps are recorded
        self.track_lons = []
        self.track_lats = []
        self.track_edges_lat = []
        self.track_edges_lon = []
        self.leading_edge_lon = None
        self.leading_edge_lat = None

    def update_position(self, lats, lons, conv_clouds, cloud_id):
        """
        Updates the position of the plume, including its centroid and point
        coordinates, archiving them in a track
        :return:
        """

        cloud_bool = conv_clouds == cloud_id

        cloud_lons = lons[cloud_bool]
        cloud_lats = lats[cloud_bool]

        # Calculate the centroid
        sum_x = np.sum(cloud_lons)
        sum_y = np.sum(cloud_lats)

        centroid_lon = sum_x/cloud_lons.shape[0]
        centroid_lat = sum_y/cloud_lons.shape[0]

        self.centroid_lon = centroid_lon
        self.centroid_lat = centroid_lat

        self.track_centroid_lat.append(centroid_lat)
        self.track_centroid_lon.append(centroid_lon)

        self.plume_lons = cloud_lons
        self.plume_lats = cloud_lats

        self.track_lons.append(cloud_lons)
        self.track_lats.append(cloud_lats)

        # Only record the previous three timesteps (including this one)
        if len(self.track_lons) > 3:
            del self.track_lons[0]

        if len(self.track_lats) > 3:
            del self.track_lats[0]

        #print 'Centroid lat', centroid_lat
        #print 'Centroid lon', centroid_lon

        #print '\n'

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