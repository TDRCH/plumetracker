"""
Handling of plumes objects, including updating of attributes and testing for
mechanism type
"""

import numpy as np
from scipy.ndimage import measurements
from scipy import ndimage as ndi
from skimage import feature
import datetime
from copy import deepcopy
from shapely.geometry import shape
import pyproj
import shapely.ops as ops
from shapely.geometry.polygon import Polygon
from functools import partial
import scipy.signal
from scipy import spatial
from scipy import sparse
import cv2
import matplotlib.pyplot as plt
import datetime as dt
import math
import shelve
from mpl_toolkits.basemap import Basemap
from skimage.measure import label, regionprops
from netCDF4 import Dataset

import utilities
import pinkdust

# Global function to scan the SDFs for unique plumes
def scan_for_plumes(sdf_now, sdf_prev, used_ids, clouds):
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

        # Plume infilling is here - switch off for now
        #sdf_now, infilled = plume_infilling(sdf_now, clouds)

        sdf_clusters, num = measurements.label(sdf_now)
        if len(used_ids) > 0:
            # Increase the plume_ID so that they are all new
            old_id_max = np.max(used_ids)
            sdf_clusters[sdf_clusters != 0] += old_id_max
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

        # Plume infilling is here - switch off for now
        #sdf_now, infilled = plume_infilling(sdf_now, clouds)

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
        # An array of previous IDs so we can check if the same ID has been
        # reassigned more than once (i.e. a split)
        overlap_prev_ids = {}
        split_ids = []
        for i in overlapping_ids:
            current = sdf_clusters[sdf_clusters == i]
            prev_ids = sdf_prev[sdf_clusters == i]
            # Remove zeros
            prev_ids = prev_ids[prev_ids != 0]
            # Take the most common of the previous IDs as the one which should
            # be applied to the new plume
            counts = np.bincount(prev_ids)
            current_size = np.bincount(current)
            prev_id = np.argmax(counts)
            if len(np.unique(prev_ids)) > 1:
                merge_ids.append(prev_id)

            old_id_array[sdf_clusters == i] = prev_id
            sdf_clusters[sdf_clusters == i] = prev_id

            """

            # Check if the ID is already a key in the dictionary
            if prev_id in overlap_prev_ids:
                # This previous ID was involved in a split
                split_ids.append(prev_id)
                if overlap_prev_ids[prev_id] < current_size:
                    # This is the largest sized plume for this previous ID
                    overlap_prev_ids[prev_id] = current_size
                    # Set prev_ID to whatever that is
                    old_id_array[sdf_clusters == i] = prev_id
                    sdf_clusters[sdf_clusters == i] = prev_id
            else:
                overlap_prev_ids[prev_id] = current_size
            """
            # The IDs we want now are the new plume IDs which are not
                # continuations of the old, so that we don't misidentify
                # them as a source. We want the plume IDs which were not
                # converted to a prev ID because another plume ID had a
                # larger area
            # We also want to identify the prev ID which split so that we
                # don't get the propagation wrong on it

            # So as for sub-250 plumes - they need to be assigned temporary
                # IDs, then as plume tracker rolls through, it gets to the
                # end of these temp plumes' lifetime, checks if they ever
                # got to 250 size and, if they did, grants them status as a
                # true plume. Otherwise they get cast into the flames of
                # oblivion

            # The above is quite a computationally intensive method,
                # however. It requires you to basically be tracking a LOT of
                #  random little baby plumes, making plume bools for them
                # and all

                # How about instead, if a plume is close to the 250
                # limit, it checks for flickering by looking at the raw
                # SDF for the previous timestep and seeing if there's
                # any overlapping SDF which are not already accounted for
                # if there are, it adds a previous timestep and does all
                #  the bells and whistles, then it checks the timestep
                # before that...HENCE IT IS THE CHAIN FLICKER CHECKER

                # If you were to do that with the plumescanner part,
                # you'd have to look at the previous timestep again,
                # and then the one before..all rather tricky

                # With plumetracker, you at least have access to all the
                #  timesteps from the tracks.

        """
        prev_id_counts = np.unique(overlap_prev_ids, return_counts=True)[1]
        for i in np.arange(0, len(prev_id_counts)):
            if prev_id_counts[i] > 1:
                split_ids.append(overlap_prev_ids[i])
                # Find the size of the plume
                # If the plume is the largest, it keeps its ID
                # Otherwise it gets assigned a new one
                # Or, as you're assigning prev ids above, add the size of
                # the overlapping id plume to a list. Then, if it turns out
                # a plume is within overlap prev ids, add the prev id to
                # split ids, then test if the size of this plume is larger
                # than the size of the previous plume wit the same assigned
                # prev id. If it's larger, it gets the prev id, and the
                # previous plume gets assigned a different one
                # So for this you need a dictionary of overlapping
                # prev ids and the size of their respective plumes

                ## WATCH OUT ## This current code creates plumes from
                # splits, which will be identified as new plumes
        """
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

    return sdf_clusters, new_ids, large_plume_ids, merge_ids#, split_ids

def update_plume_probs():
    """
    Updates the probabilities of all plumes with any conditions or new
    predictors (without needing to run plumetracker in its entirety)
    :return:
    """

    # Get cpo archives
    ca2004 = shelve.open('/soge-home/projects/seviri_dust/plumetracker/'
                         'cpo_archive_2004')
    ca2005 = shelve.open('/soge-home/projects/seviri_dust/plumetracker/'
                         'cpo_archive_2005')
    ca2006 = shelve.open('/soge-home/projects/seviri_dust/plumetracker/'
                         'cpo_archive_2006')
    ca2007 = shelve.open('/soge-home/projects/seviri_dust/plumetracker/'
                         'cpo_archive_2007')
    ca2008 = shelve.open('/soge-home/projects/seviri_dust/plumetracker/'
                         'cpo_archive_2008')
    ca2009 = shelve.open('/soge-home/projects/seviri_dust/plumetracker/'
                         'cpo_archive_2009')
    ca2010 = shelve.open('/soge-home/projects/seviri_dust/plumetracker/'
                         'cpo_archive_2010')
    cpo_archives = [ca2004, ca2005, ca2006, ca2007, ca2008, ca2009, ca2010]

    # Get plume archives
    pa2004 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                         'plume_archive_flicker_v3_prob_v2_2004')
    pa2005 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                         'plume_archive_flicker_v3_prob_v2_2005')
    pa2006 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                         'plume_archive_flicker_v3_prob_v2_2006')
    pa2007 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                         'plume_archive_flicker_v3_prob_v2_2007')
    pa2008 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                         'plume_archive_flicker_v3_prob_v2_2008')
    pa2009 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                         'plume_archive_flicker_v3_prob_v2_2009')
    pa2010 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                         'plume_archive_flicker_v3_prob_v2_2010')
    plume_archives = [pa2004, pa2005, pa2006, pa2007, pa2008, pa2009, pa2010]

    sdf_test = Dataset(
        '/soge-home/data_not_backed_up/satellite/meteosat/seviri'
        '/15-min/0.03x0.03/sdf/nc/JUNE2010/SDF_v2/SDF_v2.'
        '201006031500.nc')

    lons, lats = np.meshgrid(sdf_test.variables['longitude'][:],
                             sdf_test.variables['latitude'][:])

    # Get lats and lons
    cloud_test = Dataset(
        '/soge-home/data/satellite/meteosat/seviri/15-min/'
        '0.03x0.03/cloudmask'
        '/nc/'
        +
        'JUNE2010_CLOUDS/eumetsat.cloud.'
        + '201006031500.nc')
    cloud_lons = cloud_test.variables['lon'][:]
    cloud_lats = cloud_test.variables['lat'][:]
    cloud_lons, cloud_lats = np.meshgrid(cloud_lons, cloud_lats)

    # Loop through all plume archives
    # Loop through all plumes
    # Check the emission time of the plume
    # If the emission time is outside the LLJ emission window, set the LLJ
    # prob to 0 (emission time should be more generous after than before)
    # If there is an active CPO at this time, check if the plume is emitted
    # within its sweep. If so, set LLJ probability to 0

    print 'Processing CPO dates'
    for archive_i in np.arange(0, len(cpo_archives)):
        # All you need is a list of IDs active for each datetime
        date_dictionary = {}

        for i in cpo_archives[archive_i]:
            if cpo_archives[archive_i][i].merged == True:
                dates_observed = cpo_archives[archive_i][i].pre_merge_dates_observed
                post_merge_dates = cpo_archives[archive_i][i].dates_observed
                for j in post_merge_dates:
                    dates_observed.append(j)
            else:
                dates_observed = cpo_archives[archive_i][i].dates_observed
            for j in dates_observed:
                if j in date_dictionary:
                    date_dictionary[j].append(i)
                else:
                    date_dictionary[j] = []
                    date_dictionary[j].append(i)

    archive_idx = 1

    print 'Processing plume probabilities'
    for plume_archive in plume_archives:
        print 'Plume archive', archive_idx
        archive_idx += 1
        archive_size = len(plume_archive)
        plume_idx = 1
        used_percentages = []
        for i in plume_archive:
            if int((float(plume_idx) / float(
                    archive_size)) * 100) % 10 == 0 and int((float(
                plume_idx) / float(archive_size)) * 100) not in \
                    used_percentages:
                print str(int((float(plume_idx) / float(
                    archive_size)) * 100)) + "%"
                # This percentage has been printed already
                used_percentages.append(int((float(plume_idx) / float(
                    archive_size)) * 100))
            plume_idx += 1
            # If the plume has merged, the plume source is found in the
            # pre-merge track
            plume = plume_archive[i]

            # Emission times - no LLJs can be emitted outside of these
            # generous time bounds
            emission_time = plume.emission_time
            if emission_time.hour < 3 or emission_time.hour > 16:
                plume.archive_LLJ_prob = deepcopy(plume.LLJ_prob)
                plume.LLJ_prob = 0

            emission_bool = plume.track_plume_bool[0]

            for date in plume.dates_observed:
                cpo_sweeps = np.zeros((cloud_lons.shape))
                if date in date_dictionary:
                    active_cpos = date_dictionary[date]
                    for j in np.arange(0, len(active_cpos)):
                        cpo_i = active_cpos[j]
                        cpo = cpo_archives[archive_idx-1][cpo_i]
                        cpo_dates = cpo.dates_observed
                        cpo_index_bool = np.asarray([k <= date for k in
                                                     cpo_dates])
                        cpo_bools_to_date = cpo.track_plume_bool[cpo_index_bool]
                        cpo_bools_to_date = np.asarray([k.toarray() for k in
                                                        cpo_bools_to_date])
                        # Find the union of all CPO bools to date
                        cpo_sweep = np.sum(cpo_bools_to_date, axis=0)
                        cpo_sweep = cpo_sweep > 0
                        cpo_sweeps += cpo_sweep

                    emission_bool_regridded = pinkdust.regrid_data(lons,
                                                               lats,
                                                       cloud_lons,
                                                       cloud_lats, emission_bool)
                    if np.any(emission_bool_regridded & cpo_sweeps):
                        # An overlap between the emitted plume and a CPO sweep
                        plume.archive_LLJ_prob = deepcopy(plume.LLJ_prob)
                        plume.LLJ_prob = 0
                        break
            print i
            plume_archive[i] = plume

    for i in plume_archives:
        i.close()
    for i in cpo_archives:
        i.close()

def plume_infilling(sdf_now, clouds):
    """
    Fills in small areas of dust surrounding convection. These plumes are
    assigned a lower confidence. But then how would we maintain the
    identification of what's been infilled...So maybe this should be done
    at the end of the plume scanning phase, or just simply afterwards in a
    separate function called by main? But what about merging. We need to run a
    comparison with and without this infilling. Also wouldn't this make sense
    as something which only runs if the plume is close to convection in the
    first place?
    :param sdf_now:
    :param clouds:
    :return:
    """

    sdf_copy = sdf_now.astype(np.uint8).copy()
    plt.close()
    plt.contourf(sdf_copy)
    plt.savefig('sdf_copy.png')
    plt.close()
    #sdf_copy[sdf_copy!=0] = 255
    # Get a version of the SDF with no holes in it
    holes = cv2.floodFill(sdf_copy, None, (0, 0), 1)
    holes = ~holes[1]
    hole_conv_colocation = (clouds.astype(np.uint8) & holes)
    plt.close()
    plt.contourf(hole_conv_colocation)
    plt.savefig('hole_conv_colocation.png')
    plt.close()
    sdf_infilled = sdf_now + hole_conv_colocation
    plt.contourf(sdf_infilled)
    plt.savefig('sdf_copy_filled.png')
    plt.close()

    return sdf_infilled, hole_conv_colocation

class Plume:

    def __init__(self, plume_id, emission_time):
        self.plume_id = plume_id
        self.emission_time = emission_time
        # Attributes to track the centroid position through time
        self.track_centroid_lat = []
        self.track_centroid_lon = []
        # The previous two timesteps are recorded
        self.track_plume_bool = []
        self.track_edges_lat = []
        self.track_edges_lon = []
        self.dates_observed = []
        self.track_area = []
        self.track_centroid_direction = []
        self.track_primary_axis_direction = []
        self.leading_edge_lon = None
        self.leading_edge_lat = None
        self.merged = False
        self.clear_LLJ = False
        self.track_speed_centroid = []
        self.in_filled = False
        self.pre_merge_track_centroid_lat = []
        self.pre_merge_track_centroid_lon = []
        self.pre_merge_track_area = []
        self.pre_merge_track_speed_centroid = []
        self.pre_merge_dates_observed = []
        self.pre_merge_track_direction = []
        self.LLJ_prob = np.nan
        self.LLJ_prob_std = np.nan
        self.mean_axis_offset = None

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

        # Convert the plume bool to a sparse matrix for storage
        sparse_plume_bool = sparse.csr_matrix(plume_bool)
        self.plume_bool = sparse_plume_bool
        self.track_plume_bool.append(self.plume_bool)

        plume_lons = lons[plume_bool]
        plume_lats = lats[plume_bool]

        # Calculate the centroid
        sum_x = np.sum(plume_lons)
        sum_y = np.sum(plume_lats)

        centroid_lon = sum_x/plume_lons.shape[0]
        centroid_lat = sum_y/plume_lons.shape[0]
        self.centroid_lon = centroid_lon
        self.centroid_lat = centroid_lat

        #print 'Lons and lats of plume', self.plume_id
        #print centroid_lon
        #print centroid_lat

        self.track_centroid_lat.append(centroid_lat)
        self.track_centroid_lon.append(centroid_lon)

        self.plume_lons = plume_lons
        self.plume_lats = plume_lats

        #print 'Centroid lat', centroid_lat
        #print 'Centroid lon', centroid_lon

        #print '\n'

    def update_area(self):

        if len(self.plume_lons) > 2 and len(self.plume_lats) > 2:

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
        #if self.merged == True:
        #    duration_since_merge = self.dates_observed[-1] - self.merge_date
        #    if duration_since_merge < datetime.timedelta(hours=0.75):
        #        self.leading_edge_lon = None
        #        self.leading_edge_lat = None
        #        self.track_edges_lon.append(edge_lons)
        #        self.track_edges_lat.append(edge_lats)

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

    def update_majorminor_axes(self, lons, lats):
        """
        Updates the major and minor axes of a plume with PCA
        :return:
        """

        if len(self.plume_lons) > 1 and len(self.plume_lats) > 1:

            image = self.plume_bool.toarray()
            label_img = label(image)
            regions = regionprops(label_img)

            for props in regions:
                orientation = props.orientation
                # Add 90 degrees so it's the offset from the vertical axis
                # for consistency with the centroid direction
                orientation = math.degrees(orientation) + 90

            self.track_primary_axis_direction.append(orientation)
            self.primary_axis_direction = orientation

    def update_duration(self, date):
        """
        Updates the duration of the dust plume with the time since its emission

        NOTE: if the plume only lasts one timestep but has backwards
        flickerchecking, it will record a duration of zero (i.e. this is a
        bug ya big goon)
        :return:
        """

        if date >= self.emission_time:
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
        if self.duration < datetime.timedelta(hours=0.25):
            self.speed_centroid = np.nan
            self.track_speed_centroid.append(self.speed_centroid)
        #elif self.merged == True:
        #    if (self.dates_observed[-1] - self.merge_date) < \
        #           datetime.timedelta(hours=0.25):
        #        self.speed_centroid = np.nan
        #        self.track_speed_centroid.append(self.speed_centroid)

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

            if secs_since_previous == 0:
                print 'WARNING: Found a timestep value of zero seconds at',
                self.dates_observed[-1]
                self.track_speed_centroid.append(np.nan)
            else:
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

        if len(self.dates_observed) > 0:
            self.merged = True
            self.merge_idx = len(self.track_plume_bool)-1

            """
            self.merge_date = self.dates_observed[-1]
            self.pre_merge_track_centroid_lat = self.track_centroid_lat
            self.pre_merge_track_centroid_lon = self.track_centroid_lon
            self.pre_merge_track_area = self.track_area
            self.pre_merge_track_speed_centroid = self.track_speed_centroid
            self.pre_merge_dates_observed = self.dates_observed
            self.pre_merge_track_centroid_direction = self.track_centroid_direction

            self.track_centroid_lat = []
            self.track_centroid_lon = []
            self.track_area = []
            self.track_speed_centroid = []
            self.dates_observed = []
            self.track_centroid_direction = []
            """

    def update_centroid_direction(self):
        """
        Identify the direction in which the centroid has travelled since the
        previous timestep
        :return:
        """

        # Only get a centroid direction if we have a previous centroid and
        # haven't merged recently
        if self.duration < datetime.timedelta(hours=0.25):
            self.centroid_direction = np.nan
            self.track_centroid_direction.append(self.centroid_direction)
       # elif self.merged == True:
        #    if (self.dates_observed[-1] - self.merge_date) < \
        #            datetime.timedelta(hours=0.25):
        #        self.centroid_direction = np.nan
        #        self.track_centroid_direction.append(self.centroid_direction)
        else:
            self.centroid_direction = utilities.\
                calculate_initial_compass_bearing((self.track_centroid_lat[-2],
                                                   self.track_centroid_lon[-2]),
                                                  (self.centroid_lat,
                                                   self.centroid_lon))


            #print self.centroid_direction

            self.track_centroid_direction.append(self.centroid_direction)

            axis_offset = utilities.get_distance(self.centroid_direction,
                                                 self.primary_axis_direction)

            #print 'Centroid direction:', self.centroid_direction
            #print 'Primary axis direction:', self.primary_axis_direction

    def check_conv_distance(self, lats, lons, clouds):
        """
        Checks how far the plume is from a deep convection cloud.
        This should be done at the beginning of the plume life cycle to
        assess its association with deep convection.
        :param clouds:
        :return:
        """

        # Get clouds from the BT field rather than the cloud mask - then we
        # can also use the BT field to get the 'pinkness' value of the plume

        edge_cloud_bool = feature.canny(clouds, sigma=0.2)
        cloud_lats = lats[edge_cloud_bool == 1]
        cloud_lons = lons[edge_cloud_bool == 1]
        centroid_lat = self.centroid_lat
        centroid_lon = self.centroid_lon
        cloud_latlons_array = np.dstack((cloud_lats, cloud_lons))[0]

        distance, index = spatial.cKDTree(cloud_latlons_array).query([
            centroid_lat, centroid_lon])

        nearest_coord = cloud_latlons_array[index]
        self.conv_distance = utilities.haversine(centroid_lon, centroid_lat,
                                       nearest_coord[1], nearest_coord[0])

    def check_conv_distance_2(self, sdf_plumes, lats, lons, clouds):
        """
        Checks how far the plume is from a deep convection cloud.
        This should be done at the beginning of the plume life cycle to
        assess its association with deep convection.
        Attempt #2: I draw a buffer around the plume, iteratively increasing it
        :param self:
        :param sdf_plumes:
        :param lats:
        :param lons:
        :param clouds:
        :return:
        """

        found_convection = False
        buffered_plume = sdf_plumes == self.plume_id

        while found_convection == False:
            # Create a 2-deep buffer around the plume
            print 'here'
            convolution = scipy.signal.convolve2d(buffered_plume, np.ones((10,
                                                                         10)),
                                                  mode='same')
            check_grid = convolution > 0
            print 'here 1'
            # Get only the edges of the clouds
            edge_cloud_bool = feature.canny(clouds, sigma=0.2)
            cloud_lats = lats[edge_cloud_bool == 1]
            cloud_lons = lons[edge_cloud_bool == 1]
            cloud_latlons_array = np.dstack((cloud_lats, cloud_lons))[0]
            print 'here 2'
            # Check if any cloud lats/lons are within the buffer
            check_lats = lats[check_grid]
            check_lons = lons[check_grid]
            check_edge_cloud_bool = edge_cloud_bool[check_grid]
            check_edge_cloud_lats = check_lats[check_edge_cloud_bool]
            check_edge_cloud_lons = check_lons[check_edge_cloud_bool]
            print 'here 3'
            # If cloud has been found, stop the search
            if len(check_edge_cloud_lats) > 0:
                found_convection = True
            print 'here 4'
            # Otherwise we build a larger buffer
            buffered_plume = convolution

        smallest_distance = \
            np.min([utilities.haversine(self.centroid_lon,
                                        self.centroid_lat,
                                        check_edge_cloud_lons[j],
                                        check_edge_cloud_lats[j])
                    for j in np.arange(0, len(check_edge_cloud_lons))])

        print smallest_distance

    def update_most_likely_source(self):
        """
        Identifies the source to which the plume was nearest at the point of
        emission. Uses sources defined by AW13. Only assigns sources to
        plumes within the CWS region. Plumes with a source not within 50km
        of a source are not assigned any most likely source. Note: a better
        approach to this would be to assign box regions for each source.
        :return:
        """

        # Hard coded emission sources taken from AW13
        source_dict = {'A': (26, -6.5), 'B': (30, -6), 'C': (23, -5),
                       'D': (24.5, -5), 'E': (21, 1.75), 'F': (34, 0),
                       'G': (20, 0.5), 'H': (21.5, 1), 'I': (26.5, 1.5),
                       'J': (22.5, 2), 'K': (24, 3), 'L': (19.5, 3),
                       'M': (21, 3), 'N': (20, 5), 'O': (20, 7.5)}

        CWS_min_lat = 15
        CWS_max_lat = 35
        CWS_min_lon = -17
        CWS_max_lon = 15

        #if self.merged == True:
        #    source_lat = self.pre_merge_track_centroid_lat[0]
        #    source_lon = self.pre_merge_track_centroid_lon[0]
        #else:
        source_lat = self.track_centroid_lat[0]
        source_lon = self.track_centroid_lon[0]

        # Find the distance of the detected source to known source regions
        source_distances = np.asarray([utilities.haversine(source_lon,
                                                           source_lat,
                                                           source_dict[j][1],
                                                           source_dict[j][0])
                                       for
                                       j in source_dict])

        dict_indices = np.asarray([j for j in source_dict])

        if np.min(source_distances) > 50:
            self.most_likely_source = None
        else:
            smallest_distance = source_distances == np.min(source_distances)
            self.most_likely_source = dict_indices[smallest_distance][0]

    def update_mechanism_likelihood(self):
        ## Criteria: ##
        # Emission time
        # Distance from deep convection
        # Shape of the plume
        # Depth and duration of nearby convection
        # ERAI 6-hourly fields
        # Proximity to an identified CPO event
        # Size of the final plume (note this risks skewing results)

        total_criteria = 3
        criteria_fulfilled = 0

        if self.dates_observed[0].hour >= 9 and self.dates_observed[0].hour \
                <= 12:
            criteria_fulfilled += 1

        if self.conv_distance > 100:
            criteria_fulfilled += 1

        mean_direction = np.nanmean(self.track_centroid_direction)
        mean_primary_axis_direction = np.nanmean(
            self.track_primary_axis_direction)
        if abs(mean_direction-mean_primary_axis_direction < 20):
            criteria_fulfilled += 1

        if criteria_fulfilled == 3:
            print 'Likely an LLJ:', self.plume_id
            print self.track_centroid_lat
            print self.track_centroid_lon
            self.clear_LLJ = True
        else:
            self.clear_LLJ = False

    def update_plume_confidence(self):
        """
        Tests whether the plume has been affected by deep convection,
        proximity to the plume size detection threshold, exceedance of the
        BT threshold, merging with other plumes, splitting
        :return:
        """

    # Ok so there could be a method to update various parameters
    # Then call these each time an object instance is created

    # Ok so if you make an __init__(self) function, python will run it every
    #  time the class is invoked.

    def update_llj_probability(self, trace):
        """
        Updates the probability of being an LLJ based off a pre-existing
        empirical multivariate logistic regression model which links
        predictor variables to the probability of a plume being an LLJ. The
        MCMC method to sample the parameter landscape could be run once
        every time plumetracker is run (not that that's necessary, since the
        training dataset is not going to change from month to month). This
        update function should only be run at the end of a plume's lifecycle
        so we can retrospectively update its LLJ probability

        NOTE: It would be good to also include the Bayesian aspect of this,
        which is that we have a distribution of probabilities from the
        sampled parameter space - in some regions this landscape may be
        undersampled so our probability may be very uncertain

        #NOTE2: You need to read up on what to actually do with the prior,
        amigo! Just using whatever the default is is not really good enough,
        and hardly Bayesian!
        :return:
        """

        if trace == None:
            self.LLJ_prob = np.nan
            self.LLJ_prob_std = np.nan
        else:
            emission_time = self.dates_observed[0]

            # Take the smallest time to the nearest 0900UTC
            nine = emission_time.replace(hour=9, minute=0, second=0)
            nine_1 = nine + dt.timedelta(days=1)
            distance_nine = abs(emission_time - nine)
            distance_nine_1 = abs(nine_1 - emission_time)
            distance_from_09 = np.min([distance_nine.total_seconds(),
                                       distance_nine_1.total_seconds()])

            # The probability is calculated as the mean of the distribution of
            # probabilities from our logistic function
            if self.mean_axis_offset != None:
                LLJ_probs = (1 / (1 + np.exp(-(trace['Intercept'] +
                                                    trace['conv_distance'] *
                                                    self.conv_distance +
                                                    trace['time_to_09'] *
                                                    distance_from_09 +
                                                    trace['axis_direction_offset'] *
                                                    self.mean_axis_offset)
                                                           )))
                mean_LLJ_prob = np.nanmean(LLJ_probs)
                stdev_LLJ_prob = np.std(LLJ_probs)

                self.LLJ_prob = mean_LLJ_prob
                self.LLJ_prob_std = stdev_LLJ_prob

                #print 'conv distance', self.conv_distance
                #print 'distance from 09', distance_from_09
                #print 'mean axis offset', self.mean_axis_offset
                #print 'LLJ prob', self.LLJ_prob
                #print 'LLJ std', self.LLJ_prob_std


    def convection_infilling(self, lats, lons, clouds):
        """
        Fills in small gaps in plumes from spots of convection, saving coords
        which have been infilled. Only run when plume is in close proximity
        to convection. Updates the plume's positional variables.
        :return:
        """

        lat_min_bool = lats >= self.bbox_bottomlat
        lat_max_bool = lats <= self.bbox_toplat
        lon_min_bool = lons >= self.bbox_leftlon
        lon_max_bool = lons <= self.bbox_rightlon
        plume_bbox_bool = (lat_min_bool & lat_max_bool & lon_min_bool &
                           lon_max_bool)

        # So now you have your boolean, loop through that, and test each
        # point for a cloud True - easy
        cloud_box_all = np.where(plume_bbox_bool, clouds, False)
        mask = cloud_box_all == 0
        rows = np.flatnonzero((~mask).sum(axis=1))
        cols = np.flatnonzero((~mask).sum(axis=0))
        cloud_box = cloud_box_all[rows.min():rows.max()+1, cols.min(
        ):cols.max()+1]

        # Need here to initialise an array of the right size
        # Count the number of True elements along a row

        sdf_plumes = self.plume_bool.toarray()

        tpb = sdf_plumes[plume_bbox_bool]
        print tpb.shape

        # So that's a boolean of the plume itself, with a whole pile of both
        #  Trues and Falses therein

        # The below is a selection of sdf plumes using the plume bbox bool,
        # leaving zeros elsewhere, so a big array of mostly zeros, leaving
        # sdf plumes, which has zeros in it too...

        # Temporarily assign zero values in sdf_plumes to 2 so we can use
        # the method below - need to do the same for the cloud box

        sdf_plumes[sdf_plumes == 0] = 2

        plume_box_all = np.where(plume_bbox_bool, sdf_plumes, 0)
        mask = plume_box_all == 0
        rows = np.flatnonzero((~mask).sum(axis=1))
        cols = np.flatnonzero((~mask).sum(axis=0))
        plume_box = plume_box_all[rows.min():rows.max() + 1, cols.min(
        ):cols.max() + 1]

        print plume_box.shape[0]*plume_box.shape[1]

        # Need to reshape these arrays so they're 2D squares rather than 1D

        if np.any(cloud_box):
            for idx, val in np.ndenumerate(cloud_box):
                if val == True:
                    # Get indices adjacent to this cloud point
                    neighbours = utilities.get_neighbours(idx,
                                                         shape=cloud_box.shape)
                    total_neighbours = len(neighbours)

                    # Count the amount of neighbouring coordinates which have
                    # dust
                    adjacent_dust = 0
                    for j in neighbours:
                        if plume_box[j[0], j[1]] == True:
                            adjacent_dust += 1

                    # Determine the percentage of adjacent coordinates which
                    # are
                    # dusty
                    adjacent_dust_percentage = adjacent_dust/total_neighbours

                    # If enough adjacent dust exists, fill in this point as an
                    # SDF
                    if adjacent_dust_percentage > 0.8:
                        plume_box[idx] = True

            # Can't we just reverse the thing where we cut down the array,
            # since sdf plumes isn't going to have anything in it except
            # plume box


            sdf_plumes[plume_bbox_bool] = plume_box.flatten()
            selection_bool = deepcopy(sdf_plumes)

            sdf_plumes = sdf_plumes.astype(int)

            sdf_plumes[selection_bool] = self.plume_id

            self.in_filled = True
            print 'Did infilling'

            # Update the plume's position using the newly calculated SDF plumes
            self.update_position(lats, lons, sdf_plumes, self.plume_id)

    def convection_infilling_2(self, lats, lons, clouds):
        """
        Fills in small gaps in plumes from spots of convection, saving coords
        which have been infilled. Only run when plume is in close proximity
        to convection. Updates the plume's positional variables.
        :return:
        """

        lat_min_bool = lats >= self.bbox_bottomlat
        lat_max_bool = lats <= self.bbox_toplat
        lon_min_bool = lons >= self.bbox_leftlon
        lon_max_bool = lons <= self.bbox_rightlon
        plume_bbox_bool = (lat_min_bool & lat_max_bool & lon_min_bool &
                           lon_max_bool)

        # So now you have your boolean, loop through that, and test each
        # point for a cloud True - easy
        cloud_box_all = (clouds & plume_bbox_bool)
        cloud_box_indices = np.indices(cloud_box_all.shape)
        cloud_box_indices = np.dstack((cloud_box_indices[0],
                                       cloud_box_indices[1]))
        cloud_box_indices_plume = cloud_box_indices[plume_bbox_bool]

        if True in cloud_box_all:

            sdf_plumes = self.plume_bool.toarray()

            for i in cloud_box_indices_plume:
                cloudval = cloud_box_all[i[0], i[1]]
                if cloudval == True:
                    # Get indices adjacent to this cloud point
                    # So this is actually slightly more complex, because you
                    #  want the neighbours of the whole cloud, not just one
                    # point
                    # So if one of the point's neighbours is also cloud,
                    # it needs to add it to a list of that cloud's points.
                    # Then, at the end, it needs to compare the number of
                    # adjacent dust plume points found with the number of
                    # edge cloud points found. OK so it needs to keep adding
                    #  edge points. Or we use canny and get the edge cloud out.

                    # So say you get the edge cloud out with canny. You are
                    # looping through all cloud points. You find one,
                    # then you get its neighbours. You go through the
                    # neighbours and find which one are cloud to get the
                    # total...but then we're stuck when it comes to little
                    # jutting out clouds surrounded by dust, no?

                    # So say you've got a wee baby cloud, which is a few
                    # dozen pixels. Test any point in that cloud and it
                    # won't have more than 2 or 3 adjacent dust points.

                    # Now say you got the centroid of that cloud

                    neighbours = utilities.get_neighbours(i,
                                                shape=cloud_box_all.shape)
                    total_neighbours = len(neighbours)

                    # Count the amount of neighbouring coordinates which have
                    # dust
                    adjacent_dust = 0
                    for j in neighbours:
                        if sdf_plumes[j[0], j[1]]:
                            print 'Found dust'
                            adjacent_dust += 1

                    # Determine the percentage of adjacent coordinates which
                    # are dusty
                    adjacent_dust_percentage = float(adjacent_dust) / \
                                               float(total_neighbours)

                    if adjacent_dust_percentage != 0:
                        print adjacent_dust_percentage

                    # If enough adjacent dust exists, fill in this point as an
                    # SDF
                    if adjacent_dust_percentage > 0.8:
                        sdf_plumes[i[0], i[1]] = True

            # Can't we just reverse the thing where we cut down the array,
            # since sdf plumes isn't going to have anything in it except
            # plume box

            selection_bool = deepcopy(sdf_plumes)
            sdf_plumes = sdf_plumes.astype(int)
            sdf_plumes[selection_bool] = self.plume_id

            self.in_filled = True
            print 'Did infilling'

            # Update the plume's position using the newly calculated SDF plumes
            self.update_position(lats, lons, sdf_plumes, self.plume_id)

    def chain_flickerchecker(self, raw_sdf_prev):
        """
        If the plume is new, checks if any plume actually exists in the
        previous timestep which was filtered out by the 250 minimum size
        :return:
        """

        if len(raw_sdf_prev) > 0:
            plume_bool = self.track_plume_bool[0].toarray()
            # Only proceed if we have more than one pixel activated
            plumes_prev = deepcopy(raw_sdf_prev)
            label_objects, nb_labels = ndi.label(plumes_prev)
            local_label_objects = deepcopy(label_objects)
            local_label_objects[~plume_bool] = 0
            sizes = np.bincount(local_label_objects.ravel())
            # Only proceed if we have plumes larger than size 14
            if len(sizes) > 1:
                if np.max(sizes[1:] > 14):
                    sizes = sizes[sizes != 0]
                    sizes_sorted = deepcopy(sizes)
                    sizes_sorted.sort()
                    second_largest = sizes_sorted[-2]
                    # We want the second largest as the largest is the
                    # background
                    mask_sizes = sizes == second_largest
                    largest_id = np.unique(local_label_objects)[mask_sizes][0]
                    mask_sizes[0] = 0
                    #plumes_prev = mask_sizes[label_objects]
                    plumes_prev[label_objects != largest_id] = False
                    plumes_prev = plumes_prev.astype(int)
                    #print np.bincount(plumes_prev.ravel())
                    plumes_prev = plumes_prev == 1
                    if np.bincount(plumes_prev.ravel())[1] > 250:
                        # We found a plume which went underneath the plume
                        # size threshold - this is a flickering plume
                        return plumes_prev, self.plume_id, True
                    return plumes_prev, self.plume_id, False
                else:
                    return [], None, False
            else:
                return [], None, False
        else:
            return [], None, False

    def process_missing_plume(self):
        """
        After running an update needed due to a missing plume, moves the
        track values to the start of the track rather than the end (missing
        plumes are always in the timestep before).
        :param sdf_plumes:
        :param datetime:
        :return:
        """

        self.track_centroid_direction = list(np.roll(
            self.track_centroid_direction, 1))
        self.centroid_direction = self.track_centroid_direction[-1]
        self.track_centroid_lat = list(np.roll(self.track_centroid_lat, 1))
        self.track_centroid_lon = list(np.roll(self.track_centroid_lon, 1))
        self.centroid_lat = self.track_centroid_lat[-1]
        self.centroid_lon = self.track_centroid_lon[-1]
        self.track_plume_bool = list(np.roll(self.track_plume_bool, 1))
        self.plume_bool = self.track_plume_bool[-1]
        #self.track_edges_lat = list(np.roll(self.track_edges_lat, 1))
        #self.track_edges_lon = list(np.roll(self.track_edges_lat, 1))
        self.dates_observed = list(np.roll(self.dates_observed, 1))
        self.track_area = list(np.roll(self.track_area, 1))
        self.area = self.track_area[-1]
        self.track_centroid_direction = list(np.roll(
            self.track_centroid_direction, 1))
        self.centroid_direction = self.track_centroid_direction[-1]
        self.track_primary_axis_direction = list(np.roll(
            self.track_primary_axis_direction, 1))
        self.primary_axis_direction = self.track_primary_axis_direction
        self.track_speed_centroid = list(np.roll(self.track_speed_centroid, 1))
        self.speed_centroid = self.track_speed_centroid[-1]

        self.emission_time = self.dates_observed[0]
        self.duration = self.dates_observed[-1]-self.emission_time

    def append_missing_plume(self, plume_to_append):
        """
        When a plume has been found to be flickered, append the new one to it
        :return:
        """

        for i in np.arange(0, len(plume_to_append.dates_observed)):
            self.track_centroid_direction.append(
                plume_to_append.track_centroid_direction[i])
            self.centroid_direction = self.track_centroid_direction[-1]
            self.track_centroid_lat.append(
                plume_to_append.track_centroid_lat[i])
            self.track_centroid_lon.append(
                plume_to_append.track_centroid_lon[i])
            self.centroid_lat = self.track_centroid_lat[-1]
            self.centroid_lon = self.track_centroid_lon[-1]
            self.track_plume_bool.append(plume_to_append.track_plume_bool[i])
            self.plume_bool = self.track_plume_bool[-1]
            #self.track_edges_lat.append(plume_to_append.track_edges_lat[i])
            #self.track_edges_lon.append(plume_to_append.track_edges_lon[i])
            self.dates_observed.append(plume_to_append.dates_observed[i])
            self.track_area.append(plume_to_append.track_area[i])
            self.area = self.track_area[-1]
            self.track_centroid_direction.append(
                plume_to_append.track_centroid_direction[i])
            self.track_primary_axis_direction.append(
                plume_to_append.track_primary_axis_direction[i])
            self.track_speed_centroid.append(
                plume_to_append.track_speed_centroid[i])
            self.speed_centroid = self.track_speed_centroid[-1]
            self.duration = self.dates_observed[-1]-self.emission_time
        self.plume_lats = plume_to_append.plume_lats
        self.plume_lons = plume_to_append.plume_lons

    def flag_dust_associated_convection(self, new_clouds, prev_dust_assoc):
        """
        From an array of clouds with unique IDs associated with each
        separate object, identifies IDs of new clouds which either overlap
        with this plume in the previous timestep, or which are directly
        adjacent to this plume in the current timestep
        :param clouds:
        :return:
        """

        prev_dust_assoc = prev_dust_assoc == 1

        # First, check for an overlap with the previous timestep's dust
        # Take the previous plume bool, and simply bool out any new cloud.
        # Store all IDs which are associated with this plume and return them.
        previous_plume_bool = self.track_plume_bool[-2]
        previous_plume_bool = previous_plume_bool.toarray()
        convolution = scipy.signal.convolve2d(previous_plume_bool, np.ones((
            9,9)),mode='same')
        check_grid = convolution > 0
        if prev_dust_assoc != None:
            check_bool = check_grid | prev_dust_assoc
        else:
            check_bool = check_grid

        dust_assoc_new_clouds = deepcopy(new_clouds)
        dust_assoc_new_clouds[~check_bool] = 0

        return dust_assoc_new_clouds, check_bool

        # Second, do the convolution thing to draw a little buffer around
        # all the new cloud, then do exactly the same thing but with the
        # present plume bool. Store all IDs which are associated with the
        # plume and return them

        # Don't you think the plume should have to be of a certain age
        # before it's allowed to create convection with an association?
        # So for example you can't just get a new plume with a big
        # convection object appearing next to it at the same time

        # Or maybe a lag should be allowed? I.e. you could look back
        #  even more timesteps

        # Before coding this, the first thing to do would be to just
        #  look at the deep convection field

        # A cloud is also dust assoc if it overlaps with previous dust
        # assoc, no?

    def update_mean_axis_offset(self):
        """
        Calculates the mean offset from parallel between the centroid
        direction of the plume and its primary axis

        parallel is either a difference closer to zero or closer to 180
        perpendicular is a difference closer to 90 or to 270

        How about including the standard deviation of the axis offset and
        centroid direction? So if there's a really high one we don't bother?
        :return:
        """

        # Take the smaller of the difference between 180 and the offset and
        # 0 and the offset

        mean_centroid_direction = np.nanmean(self.track_centroid_direction)
        mean_primary_axis_direction = np.nanmean(
            self.track_primary_axis_direction)

        # Get the direction from the total distance travelled by the plume
        first_position_lat = self.track_centroid_lat[0]
        first_position_lon = self.track_centroid_lon[0]
        last_position_lat = self.track_centroid_lat[-1]
        last_position_lon = self.track_centroid_lon[-1]

        # Get the direction
        final_direction = utilities.calculate_initial_compass_bearing((
            first_position_lat, first_position_lon), (last_position_lat,
                                                      last_position_lon))

        std_centroid_direction = np.nanstd(self.track_centroid_direction)
        std_primary_axis_direction = np.nanstd(
            self.track_primary_axis_direction)

        #print '\n'
        #print std_centroid_direction
        #print std_primary_axis_direction

        mean_axis_difference = abs(utilities.get_distance(
            final_direction, mean_primary_axis_direction))

        diff_180 = abs(mean_axis_difference-180)
        diff_0 = abs(mean_axis_difference)

        self.mean_axis_offset = np.min([diff_180, diff_0])

        if np.isnan(self.mean_axis_offset):
            self.mean_axis_offset = None
        #print '\n'

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

    def update_axes(self):
        pass

    def update_mechanism_likelihood(self):
        pass