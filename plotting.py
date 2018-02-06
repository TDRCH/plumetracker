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
from scipy import ndimage as ndi
import cmocean
from matplotlib import colors
import shelve
from netCDF4 import num2date
import datetime as dt
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
import pandas as pd
import seaborn
from copy import deepcopy
from pyresample import utils
from skimage.morphology import skeletonize

import pinkdust
import utilities

def plot_cpofs():
    """
    Plots all CPOFs which contain a CPO
    :return:
    """
    root = '/ouce-home/projects/seviri_dust/cpof/'

    year_lower = 2011
    year_upper = 2011
    month_lower = 6
    month_upper = 8
    day_lower = 1
    day_upper = 31
    hour_lower = 0
    hour_upper = 23
    minute_lower = 0
    minute_upper = 45

    time_params = np.array([year_lower, year_upper, month_lower,
                            month_upper, day_lower, day_upper,
                            hour_lower, hour_upper, minute_lower,
                            minute_upper])

    datetimes = utilities.get_datetime_objects(time_params)

    # Get lats and lons
    sdf_test = Dataset(
        '/soge-home/data_not_backed_up/satellite/meteosat/seviri'
        '/15-min/0.03x0.03/sdf/nc/JUNE2010/SDF_v2/SDF_v2.'
        '201006031500.nc')

    lons, lats = np.meshgrid(sdf_test.variables['longitude'][:],
                             sdf_test.variables['latitude'][:])

    fig, ax = plt.subplots()
    extent = (np.min(lons), np.max(lons), np.min(lats), np.max(lats))
    m = Basemap(projection='cyl', llcrnrlon=extent[0], urcrnrlon=extent[1],
                llcrnrlat=extent[2], urcrnrlat=extent[3],
                resolution='i')

    m.drawcoastlines(linewidth=0.5)
    m.drawcountries(linewidth=0.5)

    parallels = np.arange(10., 40, 2.)
    # labels = [left,right,top,bottom]
    m.drawparallels(parallels, labels=[False, True, True, False],
                    linewidth=0.5)
    meridians = np.arange(-20., 17., 4.)
    m.drawmeridians(meridians, labels=[True, False, False, True],
                    linewidth=0.5)

    for i in np.arange(0, len(datetimes)):
        print datetimes[i]
        cpof_dataset = Dataset(root+datetimes[i].strftime(
            "%B%Y")+"/CPOF_"+datetimes[i].strftime("%Y%m%d%H%M")+".nc")
        cpof_data = cpof_dataset.variables['CPOF'][:]
        if 1 in cpof_data:
            contourplot = m.contourf(lons, lats, cpof_data)
            plt.savefig("CPOF_"+datetimes[i].strftime("%Y%m%d%H%M")+".png")
            for coll in contourplot.collections:
                coll.remove()
        else:
            continue
    plt.close()

def plot_cpof_plume_cloud_overviews():
    """
    Plots CPOs, as well as dust and clouds
    :return:
    """
    cpof_root = '/ouce-home/projects/seviri_dust/cpof/'
    cloud_root = '/soge-home/data/satellite/meteosat/seviri/15-min/' \
                 '0.03x0.03/cloudmask/nc/'
    dust_root = '/soge-home/data_not_backed_up/satellite/meteosat/' \
                'seviri/15-min/0.03x0.03/sdf/nc/'

    year_lower = 2010
    year_upper = 2010
    month_lower = 6
    month_upper = 8
    day_lower = 1
    day_upper = 31
    hour_lower = 0
    hour_upper = 23
    minute_lower = 0
    minute_upper = 45

    time_params = np.array([year_lower, year_upper, month_lower,
                            month_upper, day_lower, day_upper,
                            hour_lower, hour_upper, minute_lower,
                            minute_upper])

    datetimes = utilities.get_datetime_objects(time_params)

    # Get lats and lons
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

    fig, ax = plt.subplots()
    extent = (np.min(lons), np.max(lons), np.min(lats), np.max(lats))
    m = Basemap(projection='cyl', llcrnrlon=extent[0], urcrnrlon=extent[1],
                llcrnrlat=extent[2], urcrnrlat=extent[3],
                resolution='i')

    m.drawcoastlines(linewidth=0.5)
    m.drawcountries(linewidth=0.5)

    for i in np.arange(0, len(datetimes)):
        print datetimes[i]
        cpof_dataset = Dataset(cpof_root + datetimes[i].strftime(
            "%B%Y") + "/CPOF_" + datetimes[i].strftime("%Y%m%d%H%M") + ".nc")
        cpof_data = cpof_dataset.variables['CPOF'][:]

        sdf_dataset = Dataset(dust_root + datetimes[i].strftime("%B").upper(
        )+datetimes[i].strftime("%Y")+"/SDF_v2/" + "SDF_v2."
                              + datetimes[i].strftime("%Y%m%d%H%M") + ".nc")
        if 'time' in sdf_dataset.variables:
            sdf_data = sdf_dataset.variables['bt108'][0]
        else:
            sdf_data = sdf_dataset.variables['bt108'][:]

        cloud_dataset = Dataset(cloud_root + datetimes[i].strftime("%B").upper(
        )+datetimes[i].strftime("%Y")+"_CLOUDS/eumetsat.cloud."+
                                datetimes[i].strftime("%Y%m%d%H%M")+".nc")
        if 'time' in cloud_dataset.variables:
            cloud_data = cloud_dataset.variables['cmask'][0]
        else:
            cloud_data = cloud_dataset.variables['cmask'][:]

        cloud_data_regridded = pinkdust.regrid_data(cloud_lons, cloud_lats,
                                                    lons, lats, cloud_data,
                                                                mesh=True)
        cpof_data = np.ma.masked_equal(cpof_data, 0)
        sdf_data = np.ma.masked_equal(sdf_data, 0)
        cloud_data_regridded = np.ma.masked_equal(cloud_data_regridded, 0)
        cpof_contour = m.pcolormesh(lons, lats, cpof_data, cmap='Blues_r')
        dust_contour = m.pcolormesh(lons, lats, sdf_data, cmap='Oranges_r')
        cloud_contour = m.pcolormesh(lons, lats, cloud_data_regridded,
                                   cmap='Greys')
        plt.savefig('cpof_dust_cloud_overview_'+datetimes[i].strftime(
            "%Y%m%d%H%M")+".png")

        cpof_contour.remove()
        dust_contour.remove()
        cloud_contour.remove()

        """
        for coll in cpof_contour.collections:
            coll.remove()

        for coll in dust_contour.collections:
            coll.remove()

        for coll in cloud_contour.collections:
            coll.remove()
        """

    plt.close()

def plot_plume_overviews():
    """
    Plots plume overview plots for each timestep used in the plume archive
    :param plume_archive:
    :return:
    """

    # Get plume archives
    pa2004 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                         'plume_archive_flicker_v3_prob_v3_2004')
    pa2005 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                         'plume_archive_flicker_v3_prob_v3_2005')
    pa2006 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                         'plume_archive_flicker_v3_prob_v3_2006')
    pa2007 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                         'plume_archive_flicker_v3_prob_v3_2007')
    pa2008 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                         'plume_archive_flicker_v3_prob_v3_2008')
    pa2009 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                         'plume_archive_flicker_v3_prob_v3_2009')
    pa2010 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                         'plume_archive_flicker_v3_prob_v3_2010')
    pa2011 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                         'plume_archive_flicker_v3_prob_v3_2011')
    pa2012 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                         'plume_archive_flicker_v3_prob_v3_2012')
    plume_archives = [pa2004, pa2005, pa2006, pa2007, pa2008, pa2009,
                      pa2010, pa2011, pa2012]

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
    lons, lats = np.meshgrid(cloud_lons, cloud_lats)

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

    fig, ax = plt.subplots()
    extent = (np.min(lons), np.max(lons), np.min(lats), np.max(lats))
    m = Basemap(projection='cyl', llcrnrlon=extent[0], urcrnrlon=extent[1],
                llcrnrlat=extent[2], urcrnrlat=extent[3],
                resolution='i')

    m.drawcoastlines(linewidth=0.5)
    m.drawcountries(linewidth=0.5)

    for plume_archive in plume_archives:
        # All you need is a list of IDs active for each datetime
        date_dictionary = {}

        for i in plume_archive:
            if plume_archive[i].merged == True:
                dates_observed = plume_archive[i].pre_merge_dates_observed
                post_merge_dates = plume_archive[i].dates_observed
                for j in post_merge_dates:
                    dates_observed.append(j)
            else:
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
        for i in sorted(date_dictionary.keys()):
            print i
            print '\n\nProgress:', str(np.round((float(plots_done)/float(
                total_dates))*100, 2))+'%\n\n'
            plume_data = np.zeros((lats.shape[0], lons.shape[1]))
            plume_data[:] = np.nan
            centroid_lats = []
            centroid_lons = []
            plume_ids = []

            active_plumes = date_dictionary[i]

            for j in active_plumes:

                print j

                if plume_archive[j].merged == True:
                    dates_observed = plume_archive[j].pre_merge_dates_observed
                    post_merge_dates = plume_archive[j].dates_observed
                    for k in post_merge_dates:
                        dates_observed.append(k)
                else:
                    dates_observed = plume_archive[j].dates_observed

                date_bool = np.asarray(dates_observed) == i

                plume_bools = np.asarray(plume_archive[j].track_plume_bool)

                plume_bool = plume_bools[date_bool]

                plume_bool = plume_bool[0].toarray()

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

            levels = np.arange(0, 1.1, 0.1)
            contourplot = m.contourf(lons, lats, plume_data, levels=levels)
            centroidplot = m.scatter(centroid_x, centroid_y, s=1)

            anns = []

            for label_i, txt in enumerate(plume_ids):
                ann = ax.annotate(txt, (centroid_x[label_i], centroid_y[
                    label_i]), fontsize=8)
                anns.append(ann)

            plt.savefig('plume_overview_'+i.strftime(
                "%Y%m%d%H%M")+'.png')

            centroidplot.remove()
            for coll in contourplot.collections:
                coll.remove()

            for ann in anns:
                ann.remove()

            plots_done += 1

def plot_cpo_skeleton_overviews():
    """
    Plots plume overview plots for each timestep used in the plume archive
    :param plume_archive:
    :return:
    """

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
    lons, lats = np.meshgrid(cloud_lons, cloud_lats)

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
    ca2011 = shelve.open('/soge-home/projects/seviri_dust/plumetracker/'
                         'cpo_archive_2011')
    ca2012 = shelve.open('/soge-home/projects/seviri_dust/plumetracker/'
                         'cpo_archive_2012')
    cpo_archives = [ca2004, ca2005, ca2006, ca2007, ca2008, ca2009, ca2010,
                    ca2011, ca2012]

    fig, ax = plt.subplots()
    extent = (np.min(lons), np.max(lons), np.min(lats), np.max(lats))
    m = Basemap(projection='cyl', llcrnrlon=extent[0], urcrnrlon=extent[1],
                llcrnrlat=extent[2], urcrnrlat=extent[3],
                resolution='i')

    m.drawcoastlines(linewidth=0.5)
    m.drawcountries(linewidth=0.5)

    for plume_archive in cpo_archives:
        # All you need is a list of IDs active for each datetime
        date_dictionary = {}

        for i in plume_archive:
            if plume_archive[i].merged == True:
                dates_observed = plume_archive[i].pre_merge_dates_observed
                post_merge_dates = plume_archive[i].dates_observed
                for j in post_merge_dates:
                    dates_observed.append(j)
            else:
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
        for i in sorted(date_dictionary.keys()):
            print i
            print '\n\nProgress:', str(np.round((float(plots_done)/float(
                total_dates))*100, 2))+'%\n\n'
            plume_data = np.zeros((lats.shape[0], lons.shape[1]))
            plume_data[:] = np.nan
            centroid_lats = []
            centroid_lons = []
            plume_ids = []

            active_plumes = date_dictionary[i]

            for j in active_plumes:

                print j

                if plume_archive[j].merged == True:
                    dates_observed = plume_archive[j].pre_merge_dates_observed
                    post_merge_dates = plume_archive[j].dates_observed
                    for k in post_merge_dates:
                        dates_observed.append(k)
                else:
                    dates_observed = plume_archive[j].dates_observed

                date_bool = np.asarray(dates_observed) == i

                plume_bools = np.asarray(plume_archive[j].track_plume_bool)

                plume_bool = plume_bools[date_bool]

                plume_bool = plume_bool[0].toarray()

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

            levels = np.arange(0, 1.1, 0.1)

            # Skeletonise plume data
            skeleton_data = skeletonize(plume_data)

            contourplot = m.contourf(lons, lats, skeleton_data, levels=levels)
            centroidplot = m.scatter(centroid_x, centroid_y, s=1)

            anns = []

            for label_i, txt in enumerate(plume_ids):
                ann = ax.annotate(txt, (centroid_x[label_i], centroid_y[
                    label_i]), fontsize=8)
                anns.append(ann)

            plt.savefig('CPO_skeleton_overview_'+i.strftime(
                "%Y%m%d%H%M")+'.png')

            centroidplot.remove()
            for coll in contourplot.collections:
                coll.remove()

            for ann in anns:
                ann.remove()

            plots_done += 1

def plot_cpo_dust_overviews():
    """
    Plots CPOs and contemporaneous dust along with its LLJ likelihood
    :param plume_archive:
    :return:
    """

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

    fig, ax = plt.subplots()
    extent = (np.min(lons), np.max(lons), np.min(lats), np.max(lats))
    m = Basemap(projection='cyl', llcrnrlon=extent[0], urcrnrlon=extent[1],
                llcrnrlat=extent[2], urcrnrlat=extent[3],
                resolution='i')

    m.drawcoastlines(linewidth=0.5)
    m.drawcountries(linewidth=0.5)

    for archive_i in np.arange(0, len(cpo_archives)):
        # All you need is a list of IDs active for each datetime
        date_dictionary = {}
        plume_date_dictionary = {}

        for i in plume_archives[archive_i]:
            if plume_archives[archive_i][i].merged == True:
                dates_observed = plume_archives[archive_i][i].pre_merge_dates_observed
                post_merge_dates = plume_archives[archive_i][i].dates_observed
                for j in post_merge_dates:
                    dates_observed.append(j)
            else:
                dates_observed = plume_archives[archive_i][i].dates_observed
            for j in dates_observed:
                if j in plume_date_dictionary:
                    plume_date_dictionary[j].append(i)
                else:
                    plume_date_dictionary[j] = []
                    plume_date_dictionary[j].append(i)

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

        total_dates = len(plume_archives[archive_i])
        plots_done = 0

        # Then make one plot per date
        for i in sorted(date_dictionary.keys()):
            plume_data = np.zeros((lats.shape[0], lons.shape[1]))
            plume_data[:] = np.nan
            centroid_lats = []
            centroid_lons = []
            plume_ids = []

            cpo_data = np.zeros((cloud_lats.shape[0], cloud_lons.shape[1]))
            cpo_data[:] = np.nan
            cpo_centroid_lats = []
            cpo_centroid_lons = []
            cpo_ids = []

            active_cpos = date_dictionary[i]
            if i in plume_date_dictionary:
                active_plumes = plume_date_dictionary[i]
            else:
                active_plumes = []

            for j in active_cpos:

                print j

                cpo = cpo_archives[archive_i][j]

                if cpo.merged == True:
                    dates_observed = cpo.pre_merge_dates_observed
                    post_merge_dates = cpo.dates_observed
                    for k in post_merge_dates:
                        dates_observed.append(k)
                else:
                    dates_observed = cpo.dates_observed

                date_bool = np.asarray(dates_observed) == i

                plume_bools = np.asarray(cpo.track_plume_bool)

                plume_bool = plume_bools[date_bool]

                plume_bool = plume_bool[0].toarray()

                cpo_data[plume_bool] = 1

                track_centroid_lats = cpo.track_centroid_lat
                track_centroid_lons = cpo.track_centroid_lon
                pre_merge_centroid_lats = cpo.pre_merge_track_centroid_lat
                pre_merge_centroid_lons = cpo.pre_merge_track_centroid_lon

                if cpo.merged == False:
                    centroid_lat = np.asarray(track_centroid_lats)[date_bool][0]
                    centroid_lon = np.asarray(track_centroid_lons)[date_bool][0]
                else:
                    combined_track_lat = np.append(pre_merge_centroid_lats,
                                                   track_centroid_lats)
                    combined_track_lon = np.append(pre_merge_centroid_lons,
                                                   track_centroid_lons)
                    centroid_lat = np.asarray(combined_track_lat)[date_bool][0]
                    centroid_lon = np.asarray(combined_track_lon)[date_bool][0]

                cpo_centroid_lats.append(centroid_lat)
                cpo_centroid_lons.append(centroid_lon)
                cpo_ids.append(cpo.plume_id)

            for j in active_plumes:

                print j

                plume = plume_archives[archive_i][j]

                if plume.merged == True:
                    dates_observed = plume.pre_merge_dates_observed
                    post_merge_dates = plume.dates_observed
                    for k in post_merge_dates:
                        dates_observed.append(k)
                else:
                    dates_observed = plume.dates_observed

                date_bool = np.asarray(dates_observed) == i

                plume_bools = np.asarray(plume.track_plume_bool)

                plume_bool = plume_bools[date_bool]

                plume_bool = plume_bool[0].toarray()

                if plume.LLJ_prob != None:
                    plume_data[plume_bool] = plume.LLJ_prob
                else:
                    plume_data[plume_bool] = 1

                track_centroid_lats = plume.track_centroid_lat
                track_centroid_lons = plume.track_centroid_lon
                pre_merge_centroid_lats = plume.pre_merge_track_centroid_lat
                pre_merge_centroid_lons = plume.pre_merge_track_centroid_lon

                if plume.merged == False:
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
                plume_ids.append(plume.plume_id)

            centroid_x, centroid_y = m(cpo_centroid_lons, cpo_centroid_lats)

            levels = np.arange(0, 1.1, 0.1)
            contourplot1 = m.contourf(cloud_lons, cloud_lats, cpo_data, \
                                            levels=levels,
                                      cmap='Blues')
            centroidplot1 = m.scatter(centroid_x, centroid_y, s=1)

            anns1 = []

            for label_i, txt in enumerate(cpo_ids):
                ann = ax.annotate(txt, (centroid_x[label_i], centroid_y[
                    label_i]), fontsize=8)
                anns1.append(ann)

            centroid_x, centroid_y = m(centroid_lons, centroid_lats)

            levels = np.arange(0, 1.1, 0.1)
            contourplot = m.contourf(lons, lats, plume_data, levels=levels)
            centroidplot = m.scatter(centroid_x, centroid_y, s=1)

            anns = []

            for label_i, txt in enumerate(plume_ids):
                ann = ax.annotate(txt, (centroid_x[label_i], centroid_y[
                    label_i]), fontsize=8)
                anns.append(ann)

            plt.savefig('CPO_dust_overview_'+i.strftime(
                "%Y%m%d%H%M")+'.png')

            centroidplot.remove()
            for coll in contourplot.collections:
                coll.remove()

            centroidplot1.remove()
            for coll in contourplot1.collections:
                coll.remove()

            for ann in anns:
                ann.remove()

            for ann in anns1:
                ann.remove()

            plots_done += 1

def plot_cpo_dust_cloud_overviews():
    """
    Plots CPOs and contemporaneous dust along with its LLJ likelihood
    :param plume_archive:
    :return:
    """

    # This one needs to loop through active dust plumes rather than active CPOs
    # Then extract the cloud from file for that date. It would be nice if we
    #  had this in an intermediary file rather than having to regrid...
    # We could work out the cloud mask from the cloudmasked_BT couldn't we?
    # Just look for nans and bob's your uncle - yeah there shouln't be any
    # nans in the unmasked stuff

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

    fig, ax = plt.subplots()
    extent = (np.min(lons), np.max(lons), np.min(lats), np.max(lats))
    m = Basemap(projection='cyl', llcrnrlon=extent[0], urcrnrlon=extent[1],
                llcrnrlat=extent[2], urcrnrlat=extent[3],
                resolution='i')

    m.drawcoastlines(linewidth=0.5)
    m.drawcountries(linewidth=0.5)

    for archive_i in np.arange(0, len(cpo_archives)):
        # All you need is a list of IDs active for each datetime
        date_dictionary = {}
        plume_date_dictionary = {}

        for i in plume_archives[archive_i]:
            if plume_archives[archive_i][i].merged == True:
                dates_observed = plume_archives[archive_i][i].pre_merge_dates_observed
                post_merge_dates = plume_archives[archive_i][i].dates_observed
                for j in post_merge_dates:
                    dates_observed.append(j)
            else:
                dates_observed = plume_archives[archive_i][i].dates_observed
            for j in dates_observed:
                if j in plume_date_dictionary:
                    plume_date_dictionary[j].append(i)
                else:
                    plume_date_dictionary[j] = []
                    plume_date_dictionary[j].append(i)

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

        total_dates = len(plume_archives[archive_i])
        plots_done = 0

        # Then make one plot per date
        for i in sorted(date_dictionary.keys()):
            plume_data = np.zeros((lats.shape[0], lons.shape[1]))
            plume_data[:] = np.nan
            centroid_lats = []
            centroid_lons = []
            plume_ids = []

            cpo_data = np.zeros((cloud_lats.shape[0], cloud_lons.shape[1]))
            cpo_data[:] = np.nan
            cpo_centroid_lats = []
            cpo_centroid_lons = []
            cpo_ids = []

            active_cpos = date_dictionary[i]
            if i in plume_date_dictionary:
                active_plumes = plume_date_dictionary[i]
            else:
                active_plumes = []

            for j in active_cpos:

                cpo = cpo_archives[archive_i][j]

                if cpo.merged == True:
                    dates_observed = cpo.pre_merge_dates_observed
                    post_merge_dates = cpo.dates_observed
                    for k in post_merge_dates:
                        dates_observed.append(k)
                else:
                    dates_observed = cpo.dates_observed

                date_bool = np.asarray(dates_observed) == i

                plume_bools = np.asarray(cpo.track_plume_bool)

                plume_bool = plume_bools[date_bool]

                plume_bool = plume_bool[0].toarray()

                cpo_data[plume_bool] = 1

                track_centroid_lats = cpo.track_centroid_lat
                track_centroid_lons = cpo.track_centroid_lon
                pre_merge_centroid_lats = cpo.pre_merge_track_centroid_lat
                pre_merge_centroid_lons = cpo.pre_merge_track_centroid_lon

                if cpo.merged == False:
                    centroid_lat = np.asarray(track_centroid_lats)[date_bool][0]
                    centroid_lon = np.asarray(track_centroid_lons)[date_bool][0]
                else:
                    combined_track_lat = np.append(pre_merge_centroid_lats,
                                                   track_centroid_lats)
                    combined_track_lon = np.append(pre_merge_centroid_lons,
                                                   track_centroid_lons)
                    centroid_lat = np.asarray(combined_track_lat)[date_bool][0]
                    centroid_lon = np.asarray(combined_track_lon)[date_bool][0]

                cpo_centroid_lats.append(centroid_lat)
                cpo_centroid_lons.append(centroid_lon)
                cpo_ids.append(cpo.plume_id)

            for j in active_plumes:

                plume = plume_archives[archive_i][j]

                if plume.merged == True:
                    dates_observed = plume.pre_merge_dates_observed
                    post_merge_dates = plume.dates_observed
                    for k in post_merge_dates:
                        dates_observed.append(k)
                else:
                    dates_observed = plume.dates_observed

                date_bool = np.asarray(dates_observed) == i

                plume_bools = np.asarray(plume.track_plume_bool)

                plume_bool = plume_bools[date_bool]

                plume_bool = plume_bool[0].toarray()

                if plume.LLJ_prob != None:
                    plume_data[plume_bool] = plume.LLJ_prob
                else:
                    plume_data[plume_bool] = 1

                track_centroid_lats = plume.track_centroid_lat
                track_centroid_lons = plume.track_centroid_lon
                pre_merge_centroid_lats = plume.pre_merge_track_centroid_lat
                pre_merge_centroid_lons = plume.pre_merge_track_centroid_lon

                if plume.merged == False:
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
                plume_ids.append(plume.plume_id)

            centroid_x, centroid_y = m(cpo_centroid_lons, cpo_centroid_lats)

            levels = np.arange(0, 1.1, 0.1)
            contourplot1 = m.contourf(cloud_lons, cloud_lats, cpo_data, \
                                            levels=levels,
                                      cmap='Blues')
            centroidplot1 = m.scatter(centroid_x, centroid_y, s=1)

            anns1 = []

            for label_i, txt in enumerate(cpo_ids):
                ann = ax.annotate(txt, (centroid_x[label_i], centroid_y[
                    label_i]), fontsize=8)
                anns1.append(ann)

            centroid_x, centroid_y = m(centroid_lons, centroid_lats)

            levels = np.arange(0, 1.1, 0.1)
            contourplot = m.contourf(lons, lats, plume_data, levels=levels)
            centroidplot = m.scatter(centroid_x, centroid_y, s=1)

            anns = []

            for label_i, txt in enumerate(plume_ids):
                ann = ax.annotate(txt, (centroid_x[label_i], centroid_y[
                    label_i]), fontsize=8)
                anns.append(ann)

            plt.savefig('CPO_dust_overview_'+i.strftime(
                "%Y%m%d%H%M")+'.png')

            centroidplot.remove()
            for coll in contourplot.collections:
                coll.remove()

            centroidplot1.remove()
            for coll in contourplot1.collections:
                coll.remove()

            for ann in anns:
                ann.remove()

            for ann in anns1:
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
    root = '/soge-home/projects/seviri_dust/sdf/June2014/'
    cloud_test = Dataset(
        '/soge-home/data/satellite/meteosat/seviri/15-min/'
        '0.03x0.03/cloudmask'
        '/nc/'
        +
        'JUNE2010_CLOUDS/eumetsat.cloud.'
        + '201006031500.nc')
    lons = cloud_test.variables['lon'][:]
    lats = cloud_test.variables['lat'][:]
    lons, lats = np.meshgrid(lons, lats)

    # Get Ian's lats and lons
    sdf_test = Dataset(
        '/soge-home/data_not_backed_up/satellite/meteosat/seviri'
        '/15-min/0.03x0.03/sdf/nc/JUNE2010/SDF_v2/SDF_v2.'
        '201006031500.nc')

    ianlons = sdf_test.variables['longitude'][:]
    ianlats = sdf_test.variables['latitude'][:]
    lons, lats = np.meshgrid(ianlons, ianlats)

    extent = (np.min(lons), np.max(lons), np.min(lats), np.max(lats))
    m = Basemap(projection='cyl', llcrnrlon=extent[0], urcrnrlon=extent[1],
                llcrnrlat=extent[2], urcrnrlat=extent[3],
                resolution='i')

    m.drawcoastlines(linewidth=0.5)
    m.drawcountries(linewidth=0.5)

    filelist = os.listdir(root)
    filelist.sort()

    for item in filelist:
        if os.path.isfile(os.path.join(root, item)):
            SDF = Dataset(root + str(item))
            #cpo = np.load(root+str(item))
            c = m.contourf(lons, lats, SDF.variables['SDF'][:])
            #c = m.contourf(lons, lats, cpo)
            plt.savefig(str(item)[:-3] + '.png')
            print str(item)[:-3]
            for coll in c.collections:
                plt.gca().collections.remove(coll)

def plot_SDF_overviews(datetimes, lats=None, lons=None):
    """
    Plots raw SDF, without 250 pixel filtering, for comparison with plume plots
    :param datetimes:
    :return:
    """

    datestrings = [j.strftime("%Y%m%d%H%M") for j in datetimes]

    extent = (np.min(lons), np.max(lons), np.min(lats), np.max(lats))
    m = Basemap(projection='cyl', llcrnrlon=extent[0],
                urcrnrlon=extent[1],
                llcrnrlat=extent[2], urcrnrlat=extent[3],
                resolution='i')

    m.drawcoastlines(linewidth=0.5)
    m.drawcountries(linewidth=0.5)

    for date_i in np.arange(0, len(datetimes)):
        print datetimes[date_i]
        try:
            sdf = Dataset(
                '/soge-home/data_not_backed_up/satellite/meteosat/seviri/'
                '15-min/0.03x0.03/sdf/nc/' +
                datetimes[date_i].strftime("%B").upper(
                ) + str(datetimes[date_i].year) + '/SDF_v2/SDF_v2.' + \
                datestrings[date_i] + '.nc')

            if 'time' in sdf.variables:
                sdf_now = sdf.variables['bt108'][0]
            else:
                sdf_now = sdf.variables['bt108'][:]

            if lats == None:
                if 'latitude' in sdf.variables:
                    lats = sdf.variables['latitude'][:]
                    lons = sdf.variables['longitude'][:]
                else:
                    raise ValueError('No lat/lon variable - please specify '
                                     'lats and lons')

            title = 'SDF_overview_map_'+datestrings[date_i]+'.png'

            # m.contourf(lons, lats, data_array)
            SDF_img = m.imshow(sdf_now, extent=extent, origin='lower',
                     interpolation='none')

            plt.tight_layout()
            plt.savefig(title, bbox_inches='tight')

            SDF_img.remove()

            plt.close()

        except:
            'No SDF data found for', datestrings[date_i]

def plot_SDF_250px_overviews(datetimes, lats=None, lons=None):
    """
    Plots raw SDF, without 250 pixel filtering, for comparison with plume plots
    :param datetimes:
    :return:
    """

    datestrings = [j.strftime("%Y%m%d%H%M") for j in datetimes]

    sdf = Dataset(
        '/soge-home/data_not_backed_up/satellite/meteosat/seviri/'
        '15-min/0.03x0.03/sdf/nc/' +
        datetimes[0].strftime("%B").upper(
        ) + str(datetimes[0].year) + '/SDF_v2/SDF_v2.' + \
        datestrings[0] + '.nc')

    if lats == None:
        if 'latitude' in sdf.variables:
            lats = sdf.variables['latitude'][:]
            lons = sdf.variables['longitude'][:]
        else:
            raise ValueError(
                'No lat/lon variable - please specify '
                'lats and lons')

    extent = (np.min(lons), np.max(lons), np.min(lats), np.max(lats))
    m = Basemap(projection='cyl', llcrnrlon=extent[0],
                urcrnrlon=extent[1],
                llcrnrlat=extent[2], urcrnrlat=extent[3],
                resolution='i')

    m.drawcoastlines(linewidth=0.5)
    m.drawcountries(linewidth=0.5)

    for date_i in np.arange(0, len(datetimes)):
        print datetimes[date_i]
        try:
            sdf = Dataset(
                '/soge-home/data_not_backed_up/satellite/meteosat/seviri/'
                '15-min/0.03x0.03/sdf/nc/' +
                datetimes[date_i].strftime("%B").upper(
                ) + str(datetimes[date_i].year) + '/SDF_v2/SDF_v2.' + \
                datestrings[date_i] + '.nc')

            if 'time' in sdf.variables:
                sdf_now = sdf.variables['bt108'][0]
            else:
                sdf_now = sdf.variables['bt108'][:]

            label_objects, nb_labels = ndi.label(sdf_now)
            sizes = np.bincount(label_objects.ravel())

            # Set clusters smaller than size 250 to zero
            mask_sizes = sizes > 250
            mask_sizes[0] = 0
            sdf_now = mask_sizes[label_objects]

            title = 'SDF_overview_map_' + datestrings[date_i] + '_250px.png'

            # m.contourf(lons, lats, data_array)
            SDF_img = m.imshow(sdf_now, extent=extent, origin='lower',
                               interpolation='none')

            plt.tight_layout()
            plt.savefig(title, bbox_inches='tight')

            SDF_img.remove()

            plt.close()

        except:
            'No SDF data found for', datestrings[date_i]


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

def plot_llj_prob_map(plume_archive, lats, lons, min_lat=None,
                            min_lon=None,
                    max_lat=None, max_lon=None, res=100):
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

    # Get lats and lons
    sdf_test = Dataset(
        '/soge-home/data_not_backed_up/satellite/meteosat/seviri'
        '/15-min/0.03x0.03/sdf/nc/JUNE2010/SDF_v2/SDF_v2.'
        '201006031500.nc')

    lons, lats = np.meshgrid(sdf_test.variables['longitude'][:],
                             sdf_test.variables['latitude'][:])

    linlats = np.linspace(np.min(lats), np.max(lats), res)
    linlons = np.linspace(np.min(lons), np.max(lons), res)

    lons, lats = np.meshgrid(linlons, linlats)

    # Get plume archives
    pa2004 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                         'plume_archive_flicker_v3_prob_2004')
    pa2005 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                         'plume_archive_flicker_v3_prob_2005')
    pa2006 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                         'plume_archive_flicker_v3_prob_2006')
    pa2007 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                         'plume_archive_flicker_v3_prob_2007')
    pa2008 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                         'plume_archive_flicker_v3_prob_2008')
    pa2009 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                         'plume_archive_flicker_v3_prob_2009')
    pa2010 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                         'plume_archive_flicker_v3_prob_2010')
    plume_archives = [pa2004, pa2005, pa2006, pa2007, pa2008, pa2009, pa2010]

    plt.close()

    fig, ax = plt.subplots()

    latlon_dictionary = {}

    # Get your zip grid of latlon coordinates
    latlon_zip = np.array(zip(lats.ravel(), lons.ravel()),
                          dtype=('f8,'
                                 'f8')).reshape(lats.shape)

    latlon_zip = latlon_zip.flatten()

    for i in plume_archive:
        # If the plume has merged, the plume source is found in the
        # pre-merge track

        emission_lat = plume_archive[i].track_centroid_lat[0]
        emission_lon = plume_archive[i].track_centroid_lon[0]
        llj_prob = plume_archive[i].LLJ_prob

        if plume_archive[i].LLJ_prob == None:
            continue

        nearest = min(latlon_zip, key=lambda x: utilities.
                      haversine(x[1],
                                x[0],
                                emission_lon,
                                emission_lat))

        # If this entry already exists in the dictionary we just append to it
        if (nearest[0], nearest[1]) in latlon_dictionary:
            latlon_dictionary[(nearest[0], nearest[1])].append(np.nanmean(
                llj_prob))
        else:
            latlon_dictionary[(nearest[0], nearest[1])] = []
            latlon_dictionary[(nearest[0], nearest[1])].append(
                np.nanmean(llj_prob))

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
    parallels = np.arange(10., 40, 2.)
    # labels = [left,right,top,bottom]
    m.drawparallels(parallels, labels=[False, True, True, False],
                    linewidth=0.5)
    meridians = np.arange(-20., 17., 2.)
    m.drawmeridians(meridians, labels=[True, False, False, True],
                    linewidth=0.5)

    #m.contourf(lons, lats, data_array)
    discrete_cmap = utilities.cmap_discretize(cm.YlGnBu, 7)
    m.imshow(data_array, extent=extent, origin='lower', interpolation='none',
             cmap=discrete_cmap, vmin=0.3, vmax=1)

    ticks = np.arange(0.3, 1.1, 0.1)
    plt.colorbar(orientation='horizontal', fraction=0.056, pad=0.06,
                 ticks=ticks)
    plt.tight_layout()
    plt.savefig('LLJ_prob_map_2010.png', bbox_inches='tight')

    np.save('prob_array_2010', data_array)

    plt.close()

def plot_multiyear_llj_prob_map(min_lat=None,
                            min_lon=None,
                    max_lat=None, max_lon=None, res=100,
                                pre_calculated=False):
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

    if pre_calculated:
        # Get lats and lons
        sdf_test = Dataset(
            '/soge-home/data_not_backed_up/satellite/meteosat/seviri'
            '/15-min/0.03x0.03/sdf/nc/JUNE2010/SDF_v2/SDF_v2.'
            '201006031500.nc')

        lons, lats = np.meshgrid(sdf_test.variables['longitude'][:],
                                 sdf_test.variables['latitude'][:])

        linlats = np.linspace(np.min(lats), np.max(lats), res)
        linlons = np.linspace(np.min(lons), np.max(lons), res)

        lons, lats = np.meshgrid(linlons, linlats)
        extent = (np.min(lons), np.max(lons), np.min(lats), np.max(lats))
        m = Basemap(projection='cyl', llcrnrlon=extent[0], urcrnrlon=extent[1],
                    llcrnrlat=extent[2], urcrnrlat=extent[3],
                    resolution='i')

        m.drawcoastlines(linewidth=0.5)
        m.drawcountries(linewidth=0.5)
        parallels = np.arange(10., 40, 2.)
        # labels = [left,right,top,bottom]
        m.drawparallels(parallels, labels=[False, True, True, False],
                        linewidth=0.5)
        meridians = np.arange(-20., 17., 4.)
        m.drawmeridians(meridians, labels=[True, False, False, True],
                        linewidth=0.5)

        prob_array = np.load('llj_prob_array_v4.npy')
        counts_array = np.load('multiyear_dust_frequency_array_noflicker.npy')

        prob_array[counts_array < 3] = np.nan

        discrete_cmap = utilities.cmap_discretize(cm.RdYlBu_r, 10)
        discrete_cmap = utilities.cmap_discretize(cm.Blues, 10)
        m.imshow(prob_array, extent=extent, origin='lower',
                 interpolation='none',
                 cmap=discrete_cmap, vmin=0, vmax=1)

        #m.contourf(lons, lats, prob_array, cmap='YlGnBu', levels=np.arange(0,
        #           1.1, 0.1))

        ticks = np.arange(0, 1.1, 0.1)
        cbar = plt.colorbar(orientation='horizontal', fraction=0.056, pad=0.06,
                            ticks=ticks)
        cbar.ax.set_xlabel('Mean LLJ probability')
        plt.tight_layout()
        plt.savefig('LLJ_prob_map_v6_altcm_masked_1.png', bbox_inches='tight')

        plt.close()

    else:

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
        #cpo_archives = [ca2005, ca2006, ca2007, ca2008, ca2009, ca2010]

        # Get plume archives
        pa2004 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                             'plume_archive_flicker_v4_prob_v4_2004')
        pa2005 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                             'plume_archive_flicker_v4_prob_v4_2005')
        pa2006 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                             'plume_archive_flicker_v4_prob_v4_2006')
        pa2007 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                             'plume_archive_flicker_v4_prob_v4_2007')
        pa2008 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                             'plume_archive_flicker_v4_prob_v4_2008')
        pa2009 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                             'plume_archive_flicker_v4_prob_v4_2009')
        pa2010 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                             'plume_archive_flicker_v4_prob_v4_2010')
        pa2011 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                             'plume_archive_flicker_v4_prob_v4_2011')
        pa2012 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                             'plume_archive_flicker_v4_prob_v4_2012')
        pa2013 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                             'plume_archive_flicker_v4_prob_v4_2013')
        pa2014 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                             'plume_archive_flicker_v4_prob_v4_2014')
        pa2015 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                             'plume_archive_flicker_v4_prob_v4_2015')
        pa2016 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                             'plume_archive_flicker_v4_prob_v4_2016')
        pa2017 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                             'plume_archive_flicker_v4_prob_v4_2017')
        plume_archives = [pa2004, pa2005, pa2006, pa2007, pa2008, pa2009,
                          pa2010, pa2011, pa2012, pa2013, pa2014, pa2015,
                          pa2016, pa2017]

        sdf_test = Dataset(
            '/soge-home/data_not_backed_up/satellite/meteosat/seviri'
            '/15-min/0.03x0.03/sdf/nc/JUNE2010/SDF_v2/SDF_v2.'
            '201006031500.nc')

        lons, lats = np.meshgrid(sdf_test.variables['longitude'][:],
                                 sdf_test.variables['latitude'][:])

        linlats = np.linspace(np.min(lats), np.max(lats), res)
        linlons = np.linspace(np.min(lons), np.max(lons), res)

        linlons, linlats = np.meshgrid(linlons, linlats)

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
        #cloud_lons, cloud_lats = np.meshgrid(cloud_lons, cloud_lats)

        plt.close()

        fig, ax = plt.subplots()

        #lons, lats = np.meshgrid(lons, lats)
        """
        print 'Processing CPO dates'
        for archive_i in np.arange(0, len(cpo_archives)):
            # All you need is a list of IDs active for each datetime
            date_dictionary = {}

            for i in cpo_archives[archive_i]:
                if cpo_archives[archive_i][i].merged == True:
                    dates_observed = cpo_archives[archive_i][
                        i].pre_merge_dates_observed
                    post_merge_dates = cpo_archives[archive_i][
                        i].dates_observed
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

        """

        archive_idx = 1

        data_array = np.zeros((linlons.shape))
        prob_array = np.zeros((linlons.shape))

        print cloud_lons.shape

        for plume_archive in plume_archives:
            print 'Plume archive', archive_idx
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

                emission_lat = plume.track_centroid_lat[0]
                emission_lon = plume.track_centroid_lon[0]
                llj_prob = plume.LLJ_prob

                if plume.LLJ_prob == None:
                    continue

                emission_time = plume.emission_time
                if emission_time.hour < 3 or emission_time.hour > 16:
                    llj_prob = 0

                emission_bool = plume.track_plume_bool[0].toarray()
                """
                for date in plume.dates_observed:
                    cpo_sweeps = np.zeros((cloud_lats.shape[0],
                                           cloud_lons.shape[0]))
                    if date in date_dictionary:
                        active_cpos = date_dictionary[date]
                        for j in np.arange(0, len(active_cpos)):
                            cpo_i = active_cpos[j]
                            cpo = cpo_archives[archive_idx - 1][cpo_i]
                            cpo_dates = cpo.dates_observed
                            cpo_index_bool = np.asarray([k <= date for k in
                                                         cpo_dates])
                            cpo_bools_to_date = np.asarray(
                                cpo.track_plume_bool)[
                                cpo_index_bool]
                            cpo_bools_to_date = np.asarray(
                                [k.toarray() for k in
                                 cpo_bools_to_date])
                            # Find the union of all CPO bools to date
                            cpo_sweep = np.sum(cpo_bools_to_date, axis=0)
                            cpo_sweep = cpo_sweep > 0
                            cpo_sweeps += cpo_sweep

                        emission_bool_regridded = pinkdust.regrid_data(lons,
                                                                       lats,
                                                                       cloud_lons,
                                                                       cloud_lats,
                                                                       emission_bool)

                        if np.any(emission_bool_regridded & (cpo_sweeps==1)):
                            # An overlap between the emitted plume and a CPO sweep
                            plume.archive_LLJ_prob = deepcopy(plume.LLJ_prob)
                            llj_prob = 0
                            break

                """

                lat_diff = np.abs(linlats - emission_lat)
                lon_diff = np.abs(linlons - emission_lon)
                lat_bool = lat_diff == lat_diff.min()
                lon_bool = lon_diff == lon_diff.min()
                union_bool = lat_bool & lon_bool

                data_array[union_bool] += 1
                prob_array[union_bool] += llj_prob

            archive_idx += 1

        # Find the nearest latlon for each dictionary key, and extract the index
        # for it

        #data_array[data_array == 0] = np.nan
        prob_array = prob_array/data_array

        np.save('llj_prob_array_v7', prob_array)

        extent = (np.min(lons), np.max(lons), np.min(lats), np.max(lats))
        m = Basemap(projection='cyl', llcrnrlon=extent[0], urcrnrlon=extent[1],
                    llcrnrlat=extent[2], urcrnrlat=extent[3],
                    resolution='i')

        m.drawcoastlines(linewidth=0.5)
        m.drawcountries(linewidth=0.5)
        parallels = np.arange(10., 40, 2.)
        # labels = [left,right,top,bottom]
        m.drawparallels(parallels, labels=[False, True, True, False],
                        linewidth=0.5)
        meridians = np.arange(-20., 17., 4.)
        m.drawmeridians(meridians, labels=[True, False, False, True],
                        linewidth=0.5)

        #m.contourf(lons, lats, data_array)
        discrete_cmap = utilities.cmap_discretize(cm.YlGnBu, 10)
        m.imshow(prob_array, extent=extent, origin='lower', interpolation='none',
                 cmap=discrete_cmap, vmin=0, vmax=1)

        ticks = np.arange(0, 1.1, 0.1)
        cbar = plt.colorbar(orientation='horizontal', fraction=0.056, pad=0.06,
                     ticks=ticks)
        cbar.ax.set_xlabel('Mean LLJ probability')
        plt.tight_layout()
        plt.savefig('LLJ_prob_map_v7_2004_2017.png', bbox_inches='tight')

        plt.close()

def plot_multiyear_llj_prob_map_std(min_lat=None,
                            min_lon=None,
                    max_lat=None, max_lon=None, res=100,
                                pre_calculated=False):
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

    if pre_calculated:
        # Get lats and lons
        sdf_test = Dataset(
            '/soge-home/data_not_backed_up/satellite/meteosat/seviri'
            '/15-min/0.03x0.03/sdf/nc/JUNE2010/SDF_v2/SDF_v2.'
            '201006031500.nc')

        lons, lats = np.meshgrid(sdf_test.variables['longitude'][:],
                                 sdf_test.variables['latitude'][:])

        linlats = np.linspace(np.min(lats), np.max(lats), res)
        linlons = np.linspace(np.min(lons), np.max(lons), res)

        lons, lats = np.meshgrid(linlons, linlats)
        extent = (np.min(lons), np.max(lons), np.min(lats), np.max(lats))
        m = Basemap(projection='cyl', llcrnrlon=extent[0], urcrnrlon=extent[1],
                    llcrnrlat=extent[2], urcrnrlat=extent[3],
                    resolution='i')

        m.drawcoastlines(linewidth=0.5)
        m.drawcountries(linewidth=0.5)
        parallels = np.arange(10., 40, 2.)
        # labels = [left,right,top,bottom]
        m.drawparallels(parallels, labels=[False, True, True, False],
                        linewidth=0.5)
        meridians = np.arange(-20., 17., 4.)
        m.drawmeridians(meridians, labels=[True, False, False, True],
                        linewidth=0.5)

        prob_array = np.load('llj_prob_array_v4.npy')
        counts_array = np.load('multiyear_dust_frequency_array.npy')

        #prob_array[counts_array < 2] = np.nan

        discrete_cmap = utilities.cmap_discretize(cm.RdYlBu_r, 10)
        m.imshow(prob_array, extent=extent, origin='lower',
                 interpolation='none',
                 cmap=discrete_cmap)#, vmin=0, vmax=1)

        #m.contourf(lons, lats, prob_array, cmap='YlGnBu', levels=np.arange(0,
        #           1.1, 0.1))

        #ticks = np.arange(0, 1.1, 0.1)
        cbar = plt.colorbar(orientation='horizontal', fraction=0.056,
                            pad=0.06)#,ticks=ticks)
        cbar.ax.set_xlabel('Mean LLJ probability StD')
        plt.tight_layout()
        plt.savefig('LLJ_prob_map_v5_std.png', bbox_inches='tight')

        plt.close()

    else:

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
        #cpo_archives = [ca2005, ca2006, ca2007, ca2008, ca2009, ca2010]

        # Get plume archives
        pa2004 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                             'plume_archive_flicker_v3_prob_v4_2004')
        pa2005 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                             'plume_archive_flicker_v3_prob_v4_2005')
        pa2006 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                             'plume_archive_flicker_v3_prob_v4_2006')
        pa2007 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                             'plume_archive_flicker_v3_prob_v4_2007')
        pa2008 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                             'plume_archive_flicker_v3_prob_v4_2008')
        pa2009 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                             'plume_archive_flicker_v3_prob_v4_2009')
        pa2010 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                             'plume_archive_flicker_v3_prob_v4_2010')
        pa2011 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                             'plume_archive_flicker_v3_prob_v4_2011')
        pa2012 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                             'plume_archive_flicker_v3_prob_v4_2012')
        plume_archives = [pa2004, pa2005, pa2006, pa2007, pa2008, pa2009,
                          pa2010, pa2011, pa2012]

        sdf_test = Dataset(
            '/soge-home/data_not_backed_up/satellite/meteosat/seviri'
            '/15-min/0.03x0.03/sdf/nc/JUNE2010/SDF_v2/SDF_v2.'
            '201006031500.nc')

        lons, lats = np.meshgrid(sdf_test.variables['longitude'][:],
                                 sdf_test.variables['latitude'][:])

        linlats = np.linspace(np.min(lats), np.max(lats), res)
        linlons = np.linspace(np.min(lons), np.max(lons), res)

        linlons, linlats = np.meshgrid(linlons, linlats)

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
        #cloud_lons, cloud_lats = np.meshgrid(cloud_lons, cloud_lats)

        plt.close()

        fig, ax = plt.subplots()

        #lons, lats = np.meshgrid(lons, lats)
        """
        print 'Processing CPO dates'
        for archive_i in np.arange(0, len(cpo_archives)):
            # All you need is a list of IDs active for each datetime
            date_dictionary = {}

            for i in cpo_archives[archive_i]:
                if cpo_archives[archive_i][i].merged == True:
                    dates_observed = cpo_archives[archive_i][
                        i].pre_merge_dates_observed
                    post_merge_dates = cpo_archives[archive_i][
                        i].dates_observed
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

        """

        archive_idx = 1

        data_array = np.zeros((linlons.shape))
        prob_array = np.zeros((linlons.shape))

        print cloud_lons.shape

        for plume_archive in plume_archives:
            print 'Plume archive', archive_idx
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

                emission_lat = plume.track_centroid_lat[0]
                emission_lon = plume.track_centroid_lon[0]
                llj_prob = plume.LLJ_prob_std

                if plume.LLJ_prob == np.nan:
                    continue

                emission_time = plume.emission_time
                if emission_time.hour < 3 or emission_time.hour > 16:
                    llj_prob = 0

                emission_bool = plume.track_plume_bool[0].toarray()
                """
                for date in plume.dates_observed:
                    cpo_sweeps = np.zeros((cloud_lats.shape[0],
                                           cloud_lons.shape[0]))
                    if date in date_dictionary:
                        active_cpos = date_dictionary[date]
                        for j in np.arange(0, len(active_cpos)):
                            cpo_i = active_cpos[j]
                            cpo = cpo_archives[archive_idx - 1][cpo_i]
                            cpo_dates = cpo.dates_observed
                            cpo_index_bool = np.asarray([k <= date for k in
                                                         cpo_dates])
                            cpo_bools_to_date = np.asarray(
                                cpo.track_plume_bool)[
                                cpo_index_bool]
                            cpo_bools_to_date = np.asarray(
                                [k.toarray() for k in
                                 cpo_bools_to_date])
                            # Find the union of all CPO bools to date
                            cpo_sweep = np.sum(cpo_bools_to_date, axis=0)
                            cpo_sweep = cpo_sweep > 0
                            cpo_sweeps += cpo_sweep

                        emission_bool_regridded = pinkdust.regrid_data(lons,
                                                                       lats,
                                                                       cloud_lons,
                                                                       cloud_lats,
                                                                       emission_bool)

                        if np.any(emission_bool_regridded & (cpo_sweeps==1)):
                            # An overlap between the emitted plume and a CPO sweep
                            plume.archive_LLJ_prob = deepcopy(plume.LLJ_prob)
                            llj_prob = 0
                            break

                """

                lat_diff = np.abs(linlats - emission_lat)
                lon_diff = np.abs(linlons - emission_lon)
                lat_bool = lat_diff == lat_diff.min()
                lon_bool = lon_diff == lon_diff.min()
                union_bool = lat_bool & lon_bool

                data_array[union_bool] += 1
                prob_array[union_bool] += llj_prob

            archive_idx += 1

        # Find the nearest latlon for each dictionary key, and extract the index
        # for it

        #data_array[data_array == 0] = np.nan
        prob_array = prob_array/data_array

        np.save('llj_prob_array_v4_std', prob_array)

        extent = (np.min(lons), np.max(lons), np.min(lats), np.max(lats))
        m = Basemap(projection='cyl', llcrnrlon=extent[0], urcrnrlon=extent[1],
                    llcrnrlat=extent[2], urcrnrlat=extent[3],
                    resolution='i')

        m.drawcoastlines(linewidth=0.5)
        m.drawcountries(linewidth=0.5)
        parallels = np.arange(10., 40, 2.)
        # labels = [left,right,top,bottom]
        m.drawparallels(parallels, labels=[False, True, True, False],
                        linewidth=0.5)
        meridians = np.arange(-20., 17., 4.)
        m.drawmeridians(meridians, labels=[True, False, False, True],
                        linewidth=0.5)

        #m.contourf(lons, lats, data_array)
        discrete_cmap = utilities.cmap_discretize(cm.YlGnBu, 10)
        m.imshow(prob_array, extent=extent, origin='lower', interpolation='none',
                 cmap=discrete_cmap, vmin=0, vmax=1)

        ticks = np.arange(0, 1.1, 0.1)
        cbar = plt.colorbar(orientation='horizontal', fraction=0.056, pad=0.06,
                     ticks=ticks)
        cbar.ax.set_xlabel('Mean LLJ probability StDev')
        plt.tight_layout()
        plt.savefig('LLJ_prob_map_v5_std.png', bbox_inches='tight')

        plt.close()

def plot_bool_multiyear_llj_prob_map(min_lat=None,
                            min_lon=None,
                    max_lat=None, max_lon=None, res=100, pre_calculated=False):
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

    if pre_calculated:
        # Get lats and lons
        sdf_test = Dataset(
            '/soge-home/data_not_backed_up/satellite/meteosat/seviri'
            '/15-min/0.03x0.03/sdf/nc/JUNE2010/SDF_v2/SDF_v2.'
            '201006031500.nc')

        lons, lats = np.meshgrid(sdf_test.variables['longitude'][:],
                                 sdf_test.variables['latitude'][:])

        linlats = np.linspace(np.min(lats), np.max(lats), res)
        linlons = np.linspace(np.min(lons), np.max(lons), res)

        lons, lats = np.meshgrid(linlons, linlats)
        extent = (np.min(lons), np.max(lons), np.min(lats), np.max(lats))
        m = Basemap(projection='cyl', llcrnrlon=extent[0], urcrnrlon=extent[1],
                    llcrnrlat=extent[2], urcrnrlat=extent[3],
                    resolution='i')

        m.drawcoastlines(linewidth=0.5)
        m.drawcountries(linewidth=0.5)
        parallels = np.arange(10., 40, 2.)
        # labels = [left,right,top,bottom]
        m.drawparallels(parallels, labels=[False, True, True, False],
                        linewidth=0.5)
        meridians = np.arange(-20., 17., 4.)
        m.drawmeridians(meridians, labels=[True, False, False, True],
                        linewidth=0.5)

        prob_array = np.load('llj_prob_array_v3.npy')
        counts_array = np.load('multiyear_dust_frequency_array.npy')

        prob_array[counts_array < 2] = np.nan

        discrete_cmap = utilities.cmap_discretize(cm.YlGnBu, 10)
        m.imshow(prob_array, extent=extent, origin='lower',
                 interpolation='none',
                 cmap=discrete_cmap, vmin=0, vmax=1)

        #m.contourf(lons, lats, prob_array, cmap='YlGnBu', levels=np.arange(0,
        #           1.1, 0.1))

        ticks = np.arange(0, 1.1, 0.1)
        cbar = plt.colorbar(orientation='horizontal', fraction=0.056, pad=0.06,
                            ticks=ticks)
        cbar.ax.set_xlabel('Mean LLJ probability')
        plt.tight_layout()
        plt.savefig('LLJ_prob_map_v3_masked.png', bbox_inches='tight')

        plt.close()

    else:

        # Get lats and lons
        sdf_test = Dataset(
            '/soge-home/data_not_backed_up/satellite/meteosat/seviri'
            '/15-min/0.03x0.03/sdf/nc/JUNE2010/SDF_v2/SDF_v2.'
            '201006031500.nc')

        lons, lats = np.meshgrid(sdf_test.variables['longitude'][:],
                                 sdf_test.variables['latitude'][:])

        linlats = np.linspace(np.min(lats), np.max(lats), res)
        linlons = np.linspace(np.min(lons), np.max(lons), res)

        lons, lats = np.meshgrid(linlons, linlats)

        # Get plume archives
        pa2004 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                             'plume_archive_flicker_v3_prob_v3_2004')
        pa2005 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                             'plume_archive_flicker_v3_prob_v3_2005')
        pa2006 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                             'plume_archive_flicker_v3_prob_v3_2006')
        pa2007 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                             'plume_archive_flicker_v3_prob_v3_2007')
        pa2008 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                             'plume_archive_flicker_v3_prob_v3_2008')
        pa2009 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                             'plume_archive_flicker_v3_prob_v3_2009')
        pa2010 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                             'plume_archive_flicker_v3_prob_v3_2010')
        plume_archives = [pa2004, pa2005, pa2006, pa2007, pa2008, pa2009, pa2010]

        plt.close()

        fig, ax = plt.subplots()

        lats = np.linspace(np.min(lats), np.max(lats), 100)
        lons = np.linspace(np.min(lons), np.max(lons), 100)

        lons, lats = np.meshgrid(lons, lats)

        archive_idx = 1

        data_array = np.zeros((lons.shape))
        prob_array = np.zeros((lons.shape))

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

                emission_lat = plume.track_centroid_lat[0]
                emission_lon = plume.track_centroid_lon[0]
                llj_prob = plume.LLJ_prob

                if plume.LLJ_prob == None:
                    continue

                lat_diff = np.abs(lats - emission_lat)
                lon_diff = np.abs(lons - emission_lon)
                lat_bool = lat_diff == lat_diff.min()
                lon_bool = lon_diff == lon_diff.min()
                union_bool = lat_bool & lon_bool

                data_array[union_bool] += 1
                prob_array[union_bool] += llj_prob

        # Find the nearest latlon for each dictionary key, and extract the index
        # for it

        data_array[data_array == 0] = np.nan
        prob_array = prob_array/data_array

        np.save('llj_prob_array_v7', prob_array)

        extent = (np.min(lons), np.max(lons), np.min(lats), np.max(lats))
        m = Basemap(projection='cyl', llcrnrlon=extent[0], urcrnrlon=extent[1],
                    llcrnrlat=extent[2], urcrnrlat=extent[3],
                    resolution='i')

        m.drawcoastlines(linewidth=0.5)
        m.drawcountries(linewidth=0.5)
        parallels = np.arange(10., 40, 2.)
        # labels = [left,right,top,bottom]
        m.drawparallels(parallels, labels=[False, True, True, False],
                        linewidth=0.5)
        meridians = np.arange(-20., 17., 4.)
        m.drawmeridians(meridians, labels=[True, False, False, True],
                        linewidth=0.5)

        #m.contourf(lons, lats, data_array)
        discrete_cmap = utilities.cmap_discretize(cm.YlGnBu, 7)
        m.imshow(prob_array, extent=extent, origin='lower', interpolation='none',
                 cmap=discrete_cmap, vmin=0, vmax=0.8)

        ticks = np.arange(0, 0.9, 0.1)
        cbar = plt.colorbar(orientation='horizontal', fraction=0.056, pad=0.06,
                     ticks=ticks)
        cbar.ax.set_xlabel('Mean LLJ probability')
        plt.tight_layout()
        plt.savefig('LLJ_prob_map_v7_2004_2017.png', bbox_inches='tight')

        plt.close()

def plot_multiyear_llj_prob_map_wzeros(plume_archives, lats, lons,
                                       min_lat=None, min_lon=None,
                                       max_lat=None, max_lon=None,
                                       total_timesteps=61824,
                                       title='LLJ_prob_map_total_wzeros.png',
                                       pre_calculated = False):

    if pre_calculated:
        plt.close()
        fig, ax = plt.subplots()
        extent = (np.min(lons), np.max(lons), np.min(lats), np.max(lats))
        m = Basemap(projection='cyl', llcrnrlon=extent[0], urcrnrlon=extent[1],
                    llcrnrlat=extent[2], urcrnrlat=extent[3],
                    resolution='i')

        m.drawcoastlines(linewidth=0.5)
        m.drawcountries(linewidth=0.5)
        parallels = np.arange(10., 40, 2.)
        # labels = [left,right,top,bottom]
        m.drawparallels(parallels, labels=[False, True, True, False],
                        linewidth=0.5)
        meridians = np.arange(-20., 17., 2.)
        m.drawmeridians(meridians, labels=[True, False, False, True],
                        linewidth=0.5)

        data_array = np.load('prob_array_total_mask5.npy')
        data_array[data_array == 0] = np.nan
        data_array[data_array < 0.5] = np.nan

        # m.contourf(lons, lats, data_array)
        discrete_cmap = utilities.cmap_discretize(cm.YlGnBu, 7)
        m.imshow(data_array, extent=extent, origin='lower',
                 interpolation='none',
                 cmap=discrete_cmap, vmin=0.3, vmax=1.0)

        ticks = np.arange(0.3, 1.1, 0.1)

        plt.colorbar(orientation='horizontal', fraction=0.056, pad=0.06,
                     ticks=ticks)
        plt.tight_layout()
        plt.savefig(title, bbox_inches='tight')

        plt.close()

    else:

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

        # Data array
        data_array = np.zeros((lons.shape))

        archive_idx = 1

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

                emission_lat = plume_archive[i].track_centroid_lat[0]
                emission_lon = plume_archive[i].track_centroid_lon[0]
                llj_prob = plume_archive[i].LLJ_prob

                if plume_archive[i].LLJ_prob == None:
                    continue

                nearest = min(latlon_zip, key=lambda x: utilities.
                              haversine(x[1],
                                        x[0],
                                        emission_lon,
                                        emission_lat))

                lat_bool = lats == nearest[0]
                lon_bool = lons == nearest[1]
                union_bool = lat_bool & lon_bool

                data_array[union_bool] += np.nanmean(llj_prob)

        data_array = data_array/total_timesteps

        extent = (np.min(lons), np.max(lons), np.min(lats), np.max(lats))
        m = Basemap(projection='cyl', llcrnrlon=extent[0], urcrnrlon=extent[1],
                    llcrnrlat=extent[2], urcrnrlat=extent[3],
                    resolution='i')

        m.drawcoastlines(linewidth=0.5)
        m.drawcountries(linewidth=0.5)
        parallels = np.arange(10., 40, 2.)
        # labels = [left,right,top,bottom]
        m.drawparallels(parallels, labels=[False, True, True, False],
                        linewidth=0.5)
        meridians = np.arange(-20., 17., 2.)
        m.drawmeridians(meridians, labels=[True, False, False, True],
                        linewidth=0.5)

        #m.contourf(lons, lats, data_array)
        discrete_cmap = utilities.cmap_discretize(cm.YlGnBu, 7)
        m.imshow(data_array, extent=extent, origin='lower', interpolation='none',
                 cmap=discrete_cmap)#, vmin=0.3, vmax=1)

        #ticks = np.arange(0.3, 1.1, 0.1)
        plt.colorbar(orientation='horizontal', fraction=0.056, pad=0.06)#,
                     #ticks=ticks)
        plt.tight_layout()
        plt.savefig('LLJ_prob_map_total_wzeros.png', bbox_inches='tight')

        np.save('prob_array_total_wzeros', data_array)

        plt.close()

def plot_emission_count_map(plume_archive,
                            lats, lons, title='plume_count_map.png', res=100,
                            min_lat=None, min_lon=None,
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
    plt.savefig(title, bbox_inches='tight')

    plt.close()

def plot_multiyear_emission_count_map(plume_archives,
                            lats, lons, title='plume_count_map_multiyear.png',
                                      res=500,
                            min_lat=None, min_lon=None,
                             max_lat=None, max_lon=None, pre_calculated=False):
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

    if pre_calculated:
        plt.close()
        fig, ax = plt.subplots()
        extent = (np.min(lons), np.max(lons), np.min(lats), np.max(lats))
        m = Basemap(projection='cyl', llcrnrlon=extent[0], urcrnrlon=extent[1],
                    llcrnrlat=extent[2], urcrnrlat=extent[3],
                    resolution='i')

        m.drawcoastlines(linewidth=0.5)
        m.drawcountries(linewidth=0.5)
        parallels = np.arange(10., 40, 2.)
        # labels = [left,right,top,bottom]
        m.drawparallels(parallels, labels=[False, True, True, False],
                        linewidth=0.5)
        meridians = np.arange(-20., 17., 2.)
        m.drawmeridians(meridians, labels=[True, False, False, True],
                        linewidth=0.5)

        data_array = np.load('multiyear_dust_frequency_array_v3.npy')

        # m.contourf(lons, lats, data_array)
        discrete_cmap = utilities.cmap_discretize(cm.YlGnBu, 10)
        m.imshow(data_array, extent=extent, origin='lower',
                 interpolation='none',
                 cmap=discrete_cmap, vmin=0.3, vmax=1.0)

        ticks = np.arange(0.3, 1.1, 0.1)

        plt.colorbar(orientation='horizontal', fraction=0.056, pad=0.06,
                     ticks=ticks)
        plt.tight_layout()
        plt.savefig(title, bbox_inches='tight')

        plt.close()

    else:

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

        for plume_archive in plume_archives:

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
                 cmap=discrete_cmap, vmin=0, vmax=25)

        plt.colorbar(orientation='horizontal', fraction=0.056, pad=0.04)
        plt.tight_layout()
        plt.savefig(title, bbox_inches='tight')

        plt.close()

        np.save('multiyear_dust_frequency_array', data_array)

def plot_multiyear_emission_count_map_2(
        title='plume_count_map_2004_2017_noflicker15.png',
                                      res=100,
                            min_lat=None, min_lon=None,
                             max_lat=None, max_lon=None, pre_calculated=False):
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

    if pre_calculated:

        sdf_test = Dataset(
            '/soge-home/data_not_backed_up/satellite/meteosat/seviri'
            '/15-min/0.03x0.03/sdf/nc/JUNE2010/SDF_v2/SDF_v2.'
            '201006031500.nc')

        lons, lats = np.meshgrid(sdf_test.variables['longitude'][:],
                                 sdf_test.variables['latitude'][:])

        plt.close()
        fig, ax = plt.subplots()
        extent = (np.min(lons), np.max(lons), np.min(lats), np.max(lats))
        m = Basemap(projection='cyl', llcrnrlon=extent[0], urcrnrlon=extent[1],
                    llcrnrlat=extent[2], urcrnrlat=extent[3],
                    resolution='i')

        m.drawcoastlines(linewidth=0.5)
        m.drawcountries(linewidth=0.5)
        parallels = np.arange(10., 40, 2.)
        # labels = [left,right,top,bottom]
        m.drawparallels(parallels, labels=[False, True, True, False],
                        linewidth=0.5)
        meridians = np.arange(-20., 17., 4.)
        m.drawmeridians(meridians, labels=[True, False, False, True],
                        linewidth=0.5)

        lats = np.linspace(np.min(lats), np.max(lats), res)
        lons = np.linspace(np.min(lons), np.max(lons), res)
        lons, lats = np.meshgrid(lons, lats)
        data_array = np.load('multiyear_dust_frequency_array_noflicker15.npy')
        centroids_lon = np.load('centroids_lon.npy')
        centroids_lat = np.load('centroids_lat.npy')

        # m.contourf(lons, lats, data_array)
        discrete_cmap = utilities.cmap_discretize(cm.YlOrBr, 10)
        #bounds = np.arange(0.3, 1.1, 0.1)
        #norm = colors.BoundaryNorm(bounds, discrete_cmap)
        m.imshow(data_array, extent=extent, origin='lower',
                 interpolation='none',
                 cmap=discrete_cmap, vmin=0, vmax=60)
        plt.colorbar(orientation='horizontal', fraction=0.056, pad=0.06)
        #m.scatter(centroids_lon, centroids_lat, s=2)
        vmin = 0.0
        vmax = 60.0
        step = vmax/10
        levels = np.arange(vmin, vmax, step)
        #m.contourf(lons, lats, data_array, cmap=cm.YlOrBr, levels=levels)
        print 'hello'
        #m.contourf(lons, lats, data_array, cmap=cmocean.cm.solar,
        #           levels=levels)

        #ticks = np.arange(0.3, 1.1, 0.1)

        plt.tight_layout()
        plt.savefig(title, bbox_inches='tight')

        plt.close()

    else:

        pa2004 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                             'plume_archive_flicker_v4_prob_v4_2004')
        pa2005 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                             'plume_archive_flicker_v4_prob_v4_2005')
        pa2006 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                             'plume_archive_flicker_v4_prob_v4_2006')
        pa2007 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                             'plume_archive_flicker_v4_prob_v4_2007')
        pa2008 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                             'plume_archive_flicker_v4_prob_v4_2008')
        pa2009 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                             'plume_archive_flicker_v4_prob_v4_2009')
        pa2010 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                             'plume_archive_flicker_v4_prob_v4_2010')
        pa2011 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                             'plume_archive_flicker_v4_prob_v4_2011')
        pa2012 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                             'plume_archive_flicker_v4_prob_v4_2012')
        pa2013 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                             'plume_archive_flicker_v4_prob_v4_2013')
        pa2014 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                             'plume_archive_flicker_v4_prob_v4_2014')
        pa2015 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                             'plume_archive_flicker_v4_prob_v4_2015')
        pa2016 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                             'plume_archive_flicker_v4_prob_v4_2016')
        pa2017 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                             'plume_archive_flicker_v4_prob_v4_2017')
        plume_archives = [pa2004, pa2005, pa2006, pa2007, pa2008, pa2009,
                          pa2010, pa2011, pa2012, pa2013, pa2014, pa2015,
                          pa2016, pa2017]

        sdf_test = Dataset(
            '/soge-home/data_not_backed_up/satellite/meteosat/seviri'
            '/15-min/0.03x0.03/sdf/nc/JUNE2010/SDF_v2/SDF_v2.'
            '201006031500.nc')

        lons, lats = np.meshgrid(sdf_test.variables['longitude'][:],
                                 sdf_test.variables['latitude'][:])
        ianlons = deepcopy(lons)
        ianlats = deepcopy(lats)

        plt.close()

        fig, ax = plt.subplots()

        sourced_array = np.zeros((lons.shape))

        latlon_dictionary = {}

        lats = np.linspace(np.min(lats), np.max(lats), res)
        lons = np.linspace(np.min(lons), np.max(lons), res)
        lons, lats = np.meshgrid(lons, lats)

        archive_idx = 1

        # Data array
        data_array = np.zeros((lons.shape))
        centroids_lon = []
        centroids_lat = []

        for plume_archive in plume_archives:
            year_data_array = np.zeros((lons.shape))
            print 'Plume archive', archive_idx
            archive_idx += 1
            archive_size = len(plume_archive)
            plume_idx = 1
            used_percentages = []
            for i in plume_archive:
                if int((float(plume_idx)/float(
                    archive_size))*100)%10 == 0 and int((float(
                    plume_idx)/float(archive_size))*100) not in \
                        used_percentages:
                    print str(int((float(plume_idx)/float(
                    archive_size))*100))+"%"
                    # This percentage has been printed already
                    used_percentages.append(int((float(plume_idx)/float(
                    archive_size))*100))
                plume_idx += 1
                # If the plume has merged, the plume source is found in the
                # pre-merge track
                #if plume_archive[i].merged == True:
                #    emission_lat = np.round(plume_archive[
                #        i].pre_merge_track_centroid_lat[0], 1)
                #    emission_lon = np.round(plume_archive[
                #        i].pre_merge_track_centroid_lon[0], 1)
                #else:
                emission_lat = np.round(plume_archive[
                    i].track_centroid_lat[0], 1)
                emission_lon = np.round(plume_archive[
                    i].track_centroid_lon[0], 1)

                lat_diff = np.abs(lats-emission_lat)
                lon_diff = np.abs(lons-emission_lon)
                lat_bool = lat_diff == lat_diff.min()
                lon_bool = lon_diff == lon_diff.min()
                union_bool = lat_bool & lon_bool
                data_array[union_bool] += 1
                year_data_array[union_bool] += 1
                centroids_lon.append(emission_lon)
                centroids_lat.append(emission_lat)

                """

                # If this entry already exists in the dictionary we just add to it
                if (emission_lat, emission_lon) in latlon_dictionary:
                    latlon_dictionary[(emission_lat, emission_lon)] += 1
                else:
                    latlon_dictionary[(emission_lat, emission_lon)] = 0
                    latlon_dictionary[(emission_lat, emission_lon)] += 1

                """
            year_data_array[year_data_array == 0] = np.nan
            np.save('dust_frequency_array_'+str(archive_idx-1),
                    year_data_array)
            """
            extent = (np.min(lons), np.max(lons), np.min(lats), np.max(lats))
            m = Basemap(projection='cyl', llcrnrlon=extent[0], urcrnrlon=extent[1],
                        llcrnrlat=extent[2], urcrnrlat=extent[3],
                        resolution='i')

            m.drawcoastlines(linewidth=0.5)
            m.drawcountries(linewidth=0.5)
            parallels = np.arange(10., 40, 2.)
            # labels = [left,right,top,bottom]
            m.drawparallels(parallels, labels=[False, True, True, False],
                            linewidth=0.5)
            meridians = np.arange(-20., 17., 4.)
            m.drawmeridians(meridians, labels=[True, False, False, True],
                            linewidth=0.5)

            # m.contourf(lons, lats, data_array)
            discrete_cmap = utilities.cmap_discretize(cm.YlOrBr, 10)
            m.imshow(year_data_array, extent=extent, origin='lower',
                     interpolation='none',
                     cmap=discrete_cmap, vmin=0, vmax=10)
            m.scatter(centroids_lon, centroids_lat, marker='+')

            plt.colorbar(orientation='horizontal', fraction=0.056, pad=0.06)
            plt.tight_layout()
            plt.savefig('dust_emission_count_'+str(archive_idx-1)+'.png',
                        bbox_inches='tight')

            plt.close()

            """

        # Find the nearest latlon for each dictionary key, and extract the index
        # for it

        """

        ll_dict_idx = 1
        for i in latlon_dictionary:
            print float(ll_dict_idx)/float(len(latlon_dictionary))
            ll_dict_idx += 1
            # HERE YOU NEED TO BE LOOKING FOR THE NEAREST ONE, NOT THE EXACT
            #  EQUAL - ALSO HOLY HELL THIS IS GONNA TAKE A LONG TIME
            lat_diff = np.abs(lats-i[0])
            lon_diff = np.abs(lons-i[1])
            lat_bool = lat_diff == lat_diff.min()
            lon_bool = lon_diff == lon_diff.min()
            #lat_bool = np.round(lats, 2) == i[0]
            #lon_bool = np.round(lons, 2) == i[1]
            union_bool = lat_bool & lon_bool
            data_array[union_bool] = latlon_dictionary[i]

        """

        data_array[data_array == 0] = np.nan

        np.save('centroids_lon', centroids_lon)
        np.save('centroids_lat', centroids_lat)

        extent = (np.min(lons), np.max(lons), np.min(lats), np.max(lats))
        m = Basemap(projection='cyl', llcrnrlon=extent[0], urcrnrlon=extent[1],
                    llcrnrlat=extent[2], urcrnrlat=extent[3],
                    resolution='i')

        m.drawcoastlines(linewidth=0.5)
        m.drawcountries(linewidth=0.5)
        parallels = np.arange(10., 40, 2.)
        # labels = [left,right,top,bottom]
        m.drawparallels(parallels, labels=[False, True, True, False],
                        linewidth=0.5)
        meridians = np.arange(-20., 17., 4.)
        m.drawmeridians(meridians, labels=[True, False, False, True],
                        linewidth=0.5)

        # m.contourf(lons, lats, data_array)
        discrete_cmap = utilities.cmap_discretize(cm.YlOrBr, 10)
        m.imshow(data_array, extent=extent, origin='lower', interpolation='none',
                 cmap=discrete_cmap, vmin=0, vmax=40)
        #m.scatter(centroids_lon, centroids_lat, marker='+')

        plt.colorbar(orientation='horizontal', fraction=0.056, pad=0.06)
        plt.tight_layout()
        plt.savefig(title, bbox_inches='tight')

        plt.close()

        np.save('multiyear_dust_frequency_array_noflicker15', data_array)
        np.save('sourced_array_noflicker15', sourced_array)

def plot_multiyear_emission_area_map(plume_archives,
                            lats, lons, title='plume_area_map_multiyear.png',
                                      res=200,
                            min_lat=None, min_lon=None,
                             max_lat=None, max_lon=None, pre_calculated=False):
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

    if pre_calculated:
        plt.close()
        fig, ax = plt.subplots()
        extent = (np.min(lons), np.max(lons), np.min(lats), np.max(lats))
        m = Basemap(projection='cyl', llcrnrlon=extent[0], urcrnrlon=extent[1],
                    llcrnrlat=extent[2], urcrnrlat=extent[3],
                    resolution='i')

        m.drawcoastlines(linewidth=0.5)
        m.drawcountries(linewidth=0.5)
        parallels = np.arange(10., 40, 2.)
        # labels = [left,right,top,bottom]
        m.drawparallels(parallels, labels=[False, True, True, False],
                        linewidth=0.5)
        meridians = np.arange(-20., 17., 2.)
        m.drawmeridians(meridians, labels=[True, False, False, True],
                        linewidth=0.5)

        lats = np.linspace(np.min(lats), np.max(lats), res)
        lons = np.linspace(np.min(lons), np.max(lons), res)
        lons, lats = np.meshgrid(lons, lats)
        data_array = np.load('multiyear_dust_frequency_array_v3.npy')
        area_data_array = np.load('multiyear_dust_area_array.npy')
        area_data_array[area_data_array > 1000000000] = np.nan
        print np.max(area_data_array[np.isfinite(area_data_array)])
        # m.contourf(lons, lats, data_array)
        discrete_cmap = utilities.cmap_discretize(cm.YlOrBr, 10)
        m.imshow(area_data_array, extent=extent, origin='lower',
                 interpolation='none',
                 cmap=discrete_cmap)
        vmin = 0.0
        vmax = 40.0
        step = vmax/10
        levels = np.arange(vmin, vmax, step)
        #m.contourf(lons, lats, data_array, cmap=cm.YlOrBr, levels=levels)
        print 'hello'
        #m.contourf(lons, lats, data_array, cmap=cmocean.cm.solar,
        #           levels=levels)

        plt.colorbar(orientation='horizontal', fraction=0.056, pad=0.06)
        plt.tight_layout()
        plt.savefig(title, bbox_inches='tight')

        plt.close()

    else:

        plt.close()

        fig, ax = plt.subplots()

        latlon_dictionary = {}

        lats = np.linspace(np.min(lats), np.max(lats), res)
        lons = np.linspace(np.min(lons), np.max(lons), res)
        lons, lats = np.meshgrid(lons, lats)

        archive_idx = 1

        # Data array
        data_array = np.zeros((lons.shape))
        area_data_array = np.zeros((lons.shape))

        for plume_archive in plume_archives:
            print 'Plume archive', archive_idx
            archive_idx += 1
            archive_size = len(plume_archive)
            plume_idx = 1
            used_percentages = []
            for i in plume_archive:
                if int((float(plume_idx)/float(
                    archive_size))*100)%10 == 0 and int((float(
                    plume_idx)/float(archive_size))*100) not in \
                        used_percentages:
                    print str(int((float(plume_idx)/float(
                    archive_size))*100))+"%"
                    # This percentage has been printed already
                    used_percentages.append(int((float(plume_idx)/float(
                    archive_size))*100))
                plume_idx += 1
                # If the plume has merged, the plume source is found in the
                # pre-merge track
                if plume_archive[i].merged == True:
                    emission_lat = np.round(plume_archive[
                        i].pre_merge_track_centroid_lat[0], 1)
                    emission_lon = np.round(plume_archive[
                        i].pre_merge_track_centroid_lon[0], 1)
                else:
                    emission_lat = np.round(plume_archive[
                        i].track_centroid_lat[0], 1)
                    emission_lon = np.round(plume_archive[
                        i].track_centroid_lon[0], 1)
                    emission_area = plume_archive[i].track_area[0]

                lat_diff = np.abs(lats-emission_lat)
                lon_diff = np.abs(lons-emission_lon)
                lat_bool = lat_diff == lat_diff.min()
                lon_bool = lon_diff == lon_diff.min()
                union_bool = lat_bool & lon_bool
                data_array[union_bool] += 1
                area_data_array[union_bool] += emission_area

        # Find the nearest latlon for each dictionary key, and extract the index
        # for it

        data_array[data_array == 0] = np.nan
        area_data_array = area_data_array/data_array

        extent = (np.min(lons), np.max(lons), np.min(lats), np.max(lats))
        m = Basemap(projection='cyl', llcrnrlon=extent[0], urcrnrlon=extent[1],
                    llcrnrlat=extent[2], urcrnrlat=extent[3],
                    resolution='i')

        m.drawcoastlines(linewidth=0.5)
        m.drawcountries(linewidth=0.5)
        parallels = np.arange(10., 40, 2.)
        # labels = [left,right,top,bottom]
        m.drawparallels(parallels, labels=[False, True, True, False],
                        linewidth=0.5)
        meridians = np.arange(-20., 17., 2.)
        m.drawmeridians(meridians, labels=[True, False, False, True],
                        linewidth=0.5)

        # m.contourf(lons, lats, data_array)
        discrete_cmap = utilities.cmap_discretize(cm.YlOrBr, 10)
        m.imshow(area_data_array, extent=extent, origin='lower',
                 interpolation='none',
                 cmap=discrete_cmap)

        plt.colorbar(orientation='horizontal', fraction=0.056, pad=0.06)
        plt.tight_layout()
        plt.savefig(title, bbox_inches='tight')

        plt.close()

        np.save('multiyear_dust_frequency_array', data_array)
        np.save('multiyear_dust_area_array', area_data_array)

def plot_multiyear_emission_time_map(title='plume_time_map_multiyear.png',
                                      res=100,
                            min_lat=None, min_lon=None,
                             max_lat=None, max_lon=None, pre_calculated=False):
    """
    Plots a map of mean emission time across all plume archives, quantified
    as distance (positive or negative) from 1200UTC (i.e. first and second
    halves of the day) on a divergent colour map
    :param plume_archive:
    :param lats:
    :param lons:
    :param min_lat:
    :param min_lon:
    :param max_lat:
    :param max_lon:
    :return:
    """

    # Get lats and lons
    sdf_test = Dataset(
        '/soge-home/data_not_backed_up/satellite/meteosat/seviri'
        '/15-min/0.03x0.03/sdf/nc/JUNE2010/SDF_v2/SDF_v2.'
        '201006031500.nc')

    lons, lats = np.meshgrid(sdf_test.variables['longitude'][:],
                             sdf_test.variables['latitude'][:])

    linlats = np.linspace(np.min(lats), np.max(lats), res)
    linlons = np.linspace(np.min(lons), np.max(lons), res)

    lons, lats = np.meshgrid(linlons, linlats)

    # Get plume archives
    pa2004 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                         'plume_archive_flicker_v3_prob_2004')
    pa2005 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                         'plume_archive_flicker_v3_prob_2005')
    pa2006 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                         'plume_archive_flicker_v3_prob_2006')
    pa2007 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                         'plume_archive_flicker_v3_prob_2007')
    pa2008 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                         'plume_archive_flicker_v3_prob_2008')
    pa2009 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                         'plume_archive_flicker_v3_prob_2009')
    pa2010 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                         'plume_archive_flicker_v3_prob_2010')
    plume_archives = [pa2004, pa2005, pa2006, pa2007, pa2008, pa2009, pa2010]

    if pre_calculated:
        plt.close()
        fig, ax = plt.subplots()
        extent = (np.min(lons), np.max(lons), np.min(lats), np.max(lats))
        m = Basemap(projection='cyl', llcrnrlon=extent[0], urcrnrlon=extent[1],
                    llcrnrlat=extent[2], urcrnrlat=extent[3],
                    resolution='i')

        m.drawcoastlines(linewidth=0.5)
        m.drawcountries(linewidth=0.5)
        parallels = np.arange(10., 40, 2.)
        # labels = [left,right,top,bottom]
        m.drawparallels(parallels, labels=[False, True, True, False],
                        linewidth=0.5)
        meridians = np.arange(-20., 17., 2.)
        m.drawmeridians(meridians, labels=[True, False, False, True],
                        linewidth=0.5)

        lats = np.linspace(np.min(lats), np.max(lats), res)
        lons = np.linspace(np.min(lons), np.max(lons), res)
        lons, lats = np.meshgrid(lons, lats)
        counts_data_array = np.load('multiyear_dust_frequency_array.npy')
        time_data_array = np.load('multiyear_dust_time_array.npy')
        time_data_array[counts_data_array < 5] = np.nan
        inv = utilities.cmap_map(lambda x: 1 - x, cm.RdYlBu_r)
        discrete_cmap = utilities.cmap_discretize(cm.RdYlBu_r, 8)
        m.imshow(time_data_array, extent=extent, origin='lower',
                 interpolation='none',
                 cmap=discrete_cmap)

        plt.colorbar(orientation='horizontal', fraction=0.056, pad=0.06)
        plt.tight_layout()
        plt.savefig(title, bbox_inches='tight')
        plt.close()

    else:

        plt.close()

        fig, ax = plt.subplots()

        latlon_dictionary = {}

        lats = np.linspace(np.min(lats), np.max(lats), res)
        lons = np.linspace(np.min(lons), np.max(lons), res)
        lons, lats = np.meshgrid(lons, lats)

        archive_idx = 1

        # Data array
        data_array = np.zeros((lons.shape))
        area_data_array = np.zeros((lons.shape))

        for plume_archive in plume_archives:
            print 'Plume archive', archive_idx
            archive_idx += 1
            archive_size = len(plume_archive)
            plume_idx = 1
            used_percentages = []
            for i in plume_archive:
                if int((float(plume_idx)/float(
                    archive_size))*100)%10 == 0 and int((float(
                    plume_idx)/float(archive_size))*100) not in \
                        used_percentages:
                    print str(int((float(plume_idx)/float(
                    archive_size))*100))+"%"
                    # This percentage has been printed already
                    used_percentages.append(int((float(plume_idx)/float(
                    archive_size))*100))
                plume_idx += 1
                # If the plume has merged, the plume source is found in the
                # pre-merge track
                plume = plume_archive[i]
                emission_lat = np.round(plume.track_centroid_lat[0], 1)
                emission_lon = np.round(plume.track_centroid_lon[0], 1)
                emission_time = plume.emission_time
                # Take the smallest time to the nearest 0900UTC
                twelve = emission_time.replace(hour=12, minute=0, second=0)
                twelve_1 = twelve + dt.timedelta(days=1)
                distance_twelve = emission_time - twelve
                distance_twelve_1 = twelve_1 - emission_time
                distances = np.asarray([distance_twelve,
                                        distance_twelve_1])
                distances_absolute = np.asarray([abs(distance_twelve),
                                        abs(distance_twelve_1)])
                selection = distances_absolute == np.min(
                    distances_absolute)
                final_distance = distances[selection][0].total_seconds()

                lat_diff = np.abs(lats-emission_lat)
                lon_diff = np.abs(lons-emission_lon)
                lat_bool = lat_diff == lat_diff.min()
                lon_bool = lon_diff == lon_diff.min()
                union_bool = lat_bool & lon_bool
                data_array[union_bool] += 1
                area_data_array[union_bool] += final_distance

        # Find the nearest latlon for each dictionary key, and extract the index
        # for it

        data_array[data_array == 0] = np.nan
        time_data_array = area_data_array/data_array

        np.save('multiyear_dust_frequency_array', data_array)
        np.save('multiyear_dust_time_array', time_data_array)

        extent = (np.min(lons), np.max(lons), np.min(lats), np.max(lats))
        m = Basemap(projection='cyl', llcrnrlon=extent[0], urcrnrlon=extent[1],
                    llcrnrlat=extent[2], urcrnrlat=extent[3],
                    resolution='i')

        m.drawcoastlines(linewidth=0.5)
        m.drawcountries(linewidth=0.5)
        parallels = np.arange(10., 40, 2.)
        # labels = [left,right,top,bottom]
        m.drawparallels(parallels, labels=[False, True, True, False],
                        linewidth=0.5)
        meridians = np.arange(-20., 17., 2.)
        m.drawmeridians(meridians, labels=[True, False, False, True],
                        linewidth=0.5)

        discrete_cmap = utilities.cmap_discretize(cm.RdYlBu_r, 8)
        time_data_array[data_array < 5] = np.nan
        m.imshow(area_data_array, extent=extent, origin='lower',
                 interpolation='none',
                 cmap=discrete_cmap)

        plt.colorbar(orientation='horizontal', fraction=0.056, pad=0.06)
        plt.tight_layout()
        plt.savefig(title, bbox_inches='tight')

        plt.close()

def plot_multiyear_emission_count_map_3(
        title='plume_count_map_multiyear_bool_v3_2004_2017_flicker15.png',
                                        res=100,
                                        min_lat=None, min_lon=None,
                                        max_lat=None, max_lon=None,
                                        pre_calculated=False):
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

    if pre_calculated:
        sdf_test = Dataset(
            '/soge-home/data_not_backed_up/satellite/meteosat/seviri'
            '/15-min/0.03x0.03/sdf/nc/JUNE2010/SDF_v2/SDF_v2.'
            '201006031500.nc')

        lons, lats = np.meshgrid(sdf_test.variables['longitude'][:],
                                 sdf_test.variables['latitude'][:])
        plt.close()
        fig, ax = plt.subplots()
        extent = (np.min(lons), np.max(lons), np.min(lats), np.max(lats))
        m = Basemap(projection='cyl', llcrnrlon=extent[0], urcrnrlon=extent[1],
                    llcrnrlat=extent[2], urcrnrlat=extent[3],
                    resolution='i')

        m.drawcoastlines(linewidth=0.5)
        m.drawcountries(linewidth=0.5)
        parallels = np.arange(10., 40, 2.)
        # labels = [left,right,top,bottom]
        m.drawparallels(parallels, labels=[False, True, True, False],
                        linewidth=0.5)
        meridians = np.arange(-20., 17., 4.)
        m.drawmeridians(meridians, labels=[True, False, False, True],
                        linewidth=0.5)

        #lats = np.linspace(np.min(lats), np.max(lats), res)
        #lons = np.linspace(np.min(lons), np.max(lons), res)
        #lons, lats = np.meshgrid(lons, lats)
        data_array = np.load(
            'multiyear_dust_frequency_array_bool_noflickercheck_2004_2012'
            '.npy')
        #data_array[data_array < 5] = np.nan

        # m.contourf(lons, lats, data_array)
        discrete_cmap = utilities.cmap_discretize(cm.YlOrBr, 10)
        m.imshow(data_array, extent=extent, origin='lower',
                 interpolation='none',
                 cmap=discrete_cmap, vmin=0, vmax=100)
        vmin = 0.0
        vmax = 40.0
        step = vmax / 10
        levels = np.arange(vmin, vmax, step)
        # m.contourf(lons, lats, data_array, cmap=cm.YlOrBr, levels=levels)
        # m.contourf(lons, lats, data_array, cmap=cmocean.cm.solar,
        #           levels=levels)

        cbar = plt.colorbar(orientation='horizontal', fraction=0.056, pad=0.06)
        plt.tight_layout()
        plt.savefig(title, bbox_inches='tight')

        plt.close()

    else:
        """
        pa2004 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                             'plume_archive_flicker_v4_prob_v4_2004')
        pa2005 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                             'plume_archive_flicker_v4_prob_v4_2005')
        pa2006 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                             'plume_archive_flicker_v4_prob_v4_2006')
        pa2007 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                             'plume_archive_flicker_v4_prob_v4_2007')
        pa2008 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                             'plume_archive_flicker_v4_prob_v4_2008')
        pa2009 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                             'plume_archive_flicker_v4_prob_v4_2009')
        pa2010 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                             'plume_archive_flicker_v4_prob_v4_2010')
        pa2011 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                             'plume_archive_flicker_v4_prob_v4_2011')
        pa2012 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                             'plume_archive_flicker_v4_prob_v4_2012')
        pa2013 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                             'plume_archive_flicker_v4_prob_v4_2013')
        pa2014 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                             'plume_archive_flicker_v4_prob_v4_2014')
        pa2015 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                             'plume_archive_flicker_v4_prob_v4_2015')
        pa2016 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                             'plume_archive_flicker_v4_prob_v4_2016')
        pa2017 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                             'plume_archive_flicker_v4_prob_v4_2017')
        plume_archives = [pa2004, pa2005, pa2006, pa2007, pa2008, pa2009,
                          pa2010, pa2011, pa2012, pa2013, pa2014, pa2015,
                          pa2016, pa2017]
        """

        pa2004 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                             'plume_archive_v3_prob_v3_2004')
        pa2005 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                             'plume_archive_v3_prob_v3_2005')
        pa2006 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                             'plume_archive_v3_prob_v3_2006')
        pa2007 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                             'plume_archive_v3_prob_v3_2007')
        pa2008 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                             'plume_archive_v3_prob_v3_2008')
        pa2009 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                             'plume_archive_v3_prob_v3_2009')
        pa2010 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                             'plume_archive_v3_prob_v3_2010')
        pa2011 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                             'plume_archive_v3_prob_v3_2011')
        pa2012 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                             'plume_archive_v3_prob_v3_2012')

        plume_archives = [pa2004, pa2005, pa2006, pa2007, pa2008, pa2009,
                          pa2010, pa2011, pa2012]

        sdf_test = Dataset(
            '/soge-home/data_not_backed_up/satellite/meteosat/seviri'
            '/15-min/0.03x0.03/sdf/nc/JUNE2010/SDF_v2/SDF_v2.'
            '201006031500.nc')

        lons, lats = np.meshgrid(sdf_test.variables['longitude'][:],
                                 sdf_test.variables['latitude'][:])

        plt.close()

        fig, ax = plt.subplots()

        latlon_dictionary = {}

        #lats = np.linspace(np.min(lats), np.max(lats), res)
        #lons = np.linspace(np.min(lons), np.max(lons), res)
        #lons, lats = np.meshgrid(lons, lats)

        archive_idx = 1

        # Data array
        data_array = np.zeros((lons.shape))
        tidihelt_array = np.zeros((lons.shape))

        # Tidihelt year dictionary
        tidihelt_dict = shelve.open('tidihelt_dict')
        tidihelt_llj_dict = shelve.open('tidihelt_llj_dict')
        tidihelt_dict.clear()
        tidihelt_llj_dict.clear()

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
                plume = plume_archive[i]
                # If the plume has merged, the plume source is found in the
                # pre-merge track
                """
                if plume_archive[i].merged == True:
                    emission_lat = np.round(plume_archive[
                                                i].pre_merge_track_centroid_lat[
                                                0], 1)
                    emission_lon = np.round(plume_archive[
                                                i].pre_merge_track_centroid_lon[
                                                0], 1)
                """

                emission_bool = plume.track_plume_bool[0]
                emission_bool = emission_bool.toarray()

                em_lats = np.round(lats[emission_bool], 0)
                em_lons = np.round(lons[emission_bool], 0)

                emission_lat = np.round(plume.track_centroid_lat[0], 0)
                emission_lon = np.round(plume.track_centroid_lon[0], 0)

                """

                if (emission_lat == 26 or emission_lat == 27) and (
                        emission_lon == 1 or emission_lon == 2):

                    print 'Found a Tidihelt with centroid lat', \
                        emission_lat
                    print 'and centroid lon', \
                        emission_lon
                    # union_bool = lat_bool & lon_bool
                    print 'Found these emission lats', np.nanmean(lats[emission_bool])
                    print 'Count', len(lats[emission_bool])
                    print 'and these emission lons', np.nanmean(lons[emission_bool])
                    print 'Count', len(lons[emission_bool])

                    tidihelt_array[emission_bool] += 1
                    print 'Tidihelt array at these points:'
                    print tidihelt_array[emission_bool]

                    print 'Plotting Tidihelt emission'
                    plt.close()
                    plt.contourf(lons, lats, emission_bool)
                    plt.savefig('Tidihelt_emission'+i+'.png')
                    plt.close()

                    emission_time = plume_archive[i].emission_time
                    emission_year = emission_time.year

                    if str(emission_year) in tidihelt_dict:
                        tidihelt_dict[str(emission_year)] += 1
                    else:
                        tidihelt_dict[str(emission_year)] = 1
                    if str(emission_year) in tidihelt_llj_dict:
                        a = tidihelt_llj_dict[str(emission_year)]
                        p = plume_archive[i].LLJ_prob
                        a.append(p)
                        tidihelt_llj_dict[str(emission_year)] = a

                    else:
                        tidihelt_llj_dict[str(emission_year)] = \
                            [plume_archive[i].LLJ_prob]

                """

                data_array[emission_bool] += 1

        print np.unique(data_array)

        data_array[data_array == 0] = np.nan

        extent = (np.min(lons), np.max(lons), np.min(lats), np.max(lats))
        m = Basemap(projection='cyl', llcrnrlon=extent[0], urcrnrlon=extent[1],
                    llcrnrlat=extent[2], urcrnrlat=extent[3],
                    resolution='i')

        m.drawcoastlines(linewidth=0.5)
        m.drawcountries(linewidth=0.5)
        parallels = np.arange(10., 40, 2.)
        # labels = [left,right,top,bottom]
        m.drawparallels(parallels, labels=[False, True, True, False],
                        linewidth=0.5)
        meridians = np.arange(-20., 17., 4.)
        m.drawmeridians(meridians, labels=[True, False, False, True],
                        linewidth=0.5)

        # m.contourf(lons, lats, data_array)
        discrete_cmap = utilities.cmap_discretize(cm.YlOrBr, 10)
        m.imshow(data_array, extent=extent, origin='lower',
                 interpolation='none',
                 cmap=discrete_cmap, vmin=0, vmax=50)

        cbar = plt.colorbar(orientation='horizontal', fraction=0.056, pad=0.06)
        plt.tight_layout()
        plt.savefig(title, bbox_inches='tight')

        plt.close()

        np.save('multiyear_dust_frequency_array_bool_noflickercheck_2004_2012',
                data_array)

        """
        print tidihelt_dict

        yearly_tidihelt = np.asarray([tidihelt_dict[j] for j in tidihelt_dict])
        yearly_tidihelt[0] = tidihelt_dict['2004']
        yearly_tidihelt[1] = tidihelt_dict['2005']
        yearly_tidihelt[2] = tidihelt_dict['2006']
        yearly_tidihelt[3] = tidihelt_dict['2007']
        yearly_tidihelt[4] = tidihelt_dict['2008']
        yearly_tidihelt[5] = tidihelt_dict['2009']
        yearly_tidihelt[6] = tidihelt_dict['2010']

        width = 0.35
        N = 7
        ind = np.arange(N)
        fig, ax = plt.subplots()
        rects1 = ax.bar(ind, yearly_tidihelt, width, color='r')
        ax.set_ylabel('Emission events')
        ax.set_xticks(ind+width/2)
        ax.set_xticklabels(('2004', '2005', '2006', '2007', '2008', '2009',
                            '2010'))
        plt.savefig('yearly_tidihelt_frequency.png')
        plt.close()

        yearly_tidihelt_llj = np.asarray([tidihelt_dict[j] for j in
                                       tidihelt_dict])
        yearly_tidihelt_llj[0] = sum(j>0.5 for j in tidihelt_llj_dict['2004'])
        yearly_tidihelt_llj[1] = sum(j>0.5 for j in tidihelt_llj_dict['2005'])
        yearly_tidihelt_llj[2] = sum(j>0.5 for j in tidihelt_llj_dict['2006'])
        yearly_tidihelt_llj[3] = sum(j>0.5 for j in tidihelt_llj_dict['2007'])
        yearly_tidihelt_llj[4] = sum(j>0.5 for j in tidihelt_llj_dict['2008'])
        yearly_tidihelt_llj[5] = sum(j>0.5 for j in tidihelt_llj_dict['2009'])
        yearly_tidihelt_llj[6] = sum(j>0.5 for j in tidihelt_llj_dict['2010'])

        width = 0.35
        N = 7
        ind = np.arange(N)
        fig, ax = plt.subplots()
        rects1 = ax.bar(ind, yearly_tidihelt_llj, width, color='b')
        ax.set_ylabel('LLJ events')
        ax.set_xticks(ind + width / 2)
        ax.set_xticklabels(('2004', '2005', '2006', '2007', '2008', '2009',
                            '2010'))
        ax.set_ylim(0, np.max(yearly_tidihelt_llj))
        plt.savefig('yearly_tidihelt_llj_frequency.png')
        plt.close()

        tidihelt_dict.close()
        tidihelt_llj_dict.close()

        """

def plot_multiyear_emission_bool_area_map(
        title='plume_count_map_multiyear_bool_area_2004_2017.png',
                                        res=100,
                                        min_lat=None, min_lon=None,
                                        max_lat=None, max_lon=None,
                                        pre_calculated=False):
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

    if pre_calculated:
        sdf_test = Dataset(
            '/soge-home/data_not_backed_up/satellite/meteosat/seviri'
            '/15-min/0.03x0.03/sdf/nc/JUNE2010/SDF_v2/SDF_v2.'
            '201006031500.nc')

        lons, lats = np.meshgrid(sdf_test.variables['longitude'][:],
                                 sdf_test.variables['latitude'][:])
        plt.close()
        fig, ax = plt.subplots()
        extent = (np.min(lons), np.max(lons), np.min(lats), np.max(lats))
        m = Basemap(projection='cyl', llcrnrlon=extent[0], urcrnrlon=extent[1],
                    llcrnrlat=extent[2], urcrnrlat=extent[3],
                    resolution='i')

        m.drawcoastlines(linewidth=0.5)
        m.drawcountries(linewidth=0.5)
        parallels = np.arange(10., 40, 2.)
        # labels = [left,right,top,bottom]
        m.drawparallels(parallels, labels=[False, True, True, False],
                        linewidth=0.5)
        meridians = np.arange(-20., 17., 4.)
        m.drawmeridians(meridians, labels=[True, False, False, True],
                        linewidth=0.5)

        #lats = np.linspace(np.min(lats), np.max(lats), res)
        #lons = np.linspace(np.min(lons), np.max(lons), res)
        #lons, lats = np.meshgrid(lons, lats)
        data_array = np.load('multiyear_dust_area_array_bool_noflicker.npy')
        #data_array[data_array < 5] = np.nan

        # m.contourf(lons, lats, data_array)
        discrete_cmap = utilities.cmap_discretize(cm.YlOrBr, 10)
        m.imshow(data_array, extent=extent, origin='lower',
                 interpolation='none',
                 cmap=discrete_cmap, vmin=0, vmax=100000000000)
        vmin = 0.0
        vmax = 40.0
        step = vmax / 10
        levels = np.arange(vmin, vmax, step)
        # m.contourf(lons, lats, data_array, cmap=cm.YlOrBr, levels=levels)
        # m.contourf(lons, lats, data_array, cmap=cmocean.cm.solar,
        #           levels=levels)

        cbar = plt.colorbar(orientation='horizontal', fraction=0.056, pad=0.06)
        plt.tight_layout()
        plt.savefig(title, bbox_inches='tight')

        plt.close()

    else:

        pa2004 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                             'plume_archive_flicker_v3_prob_v4_2004')
        pa2005 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                             'plume_archive_flicker_v3_prob_v4_2005')
        pa2006 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                             'plume_archive_flicker_v3_prob_v4_2006')
        pa2007 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                             'plume_archive_flicker_v3_prob_v4_2007')
        pa2008 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                             'plume_archive_flicker_v3_prob_v4_2008')
        pa2009 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                             'plume_archive_flicker_v3_prob_v4_2009')
        pa2010 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                             'plume_archive_flicker_v3_prob_v4_2010')
        pa2011 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                             'plume_archive_flicker_v3_prob_v3_2011')
        pa2012 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                             'plume_archive_flicker_v3_prob_v3_2012')
        pa2013 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                             'plume_archive_flicker_v3_prob_v4_2013')
        pa2014 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                             'plume_archive_flicker_v3_prob_v4_2014')
        pa2015 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                             'plume_archive_flicker_v3_prob_v4_2015')
        pa2016 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                             'plume_archive_flicker_v3_prob_v4_2016')
        pa2017 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                             'plume_archive_flicker_v3_prob_v4_2017')
        plume_archives = [pa2004, pa2005, pa2006, pa2007, pa2008, pa2009,
                          pa2010, pa2011, pa2012, pa2013, pa2014, pa2015,
                          pa2016, pa2017]

        sdf_test = Dataset(
            '/soge-home/data_not_backed_up/satellite/meteosat/seviri'
            '/15-min/0.03x0.03/sdf/nc/JUNE2010/SDF_v2/SDF_v2.'
            '201006031500.nc')

        lons, lats = np.meshgrid(sdf_test.variables['longitude'][:],
                                 sdf_test.variables['latitude'][:])

        plt.close()

        fig, ax = plt.subplots()

        latlon_dictionary = {}

        #lats = np.linspace(np.min(lats), np.max(lats), res)
        #lons = np.linspace(np.min(lons), np.max(lons), res)
        #lons, lats = np.meshgrid(lons, lats)

        archive_idx = 1

        # Data array
        data_array = np.zeros((lons.shape))
        tidihelt_array = np.zeros((lons.shape))

        # Tidihelt year dictionary
        tidihelt_dict = shelve.open('tidihelt_dict')
        tidihelt_llj_dict = shelve.open('tidihelt_llj_dict')
        tidihelt_dict.clear()
        tidihelt_llj_dict.clear()

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
                plume = plume_archive[i]
                # If the plume has merged, the plume source is found in the
                # pre-merge track
                emission_bool = plume.track_plume_bool[0]
                emission_bool = emission_bool.toarray()

                emission_area = plume.track_area[0]

                em_lats = np.round(lats[emission_bool], 0)
                em_lons = np.round(lons[emission_bool], 0)

                emission_lat = np.round(plume.track_centroid_lat[0], 0)
                emission_lon = np.round(plume.track_centroid_lon[0], 0)

                data_array[emission_bool] += emission_area

        data_array[data_array == 0] = np.nan

        counts_array = np.load(
            'multiyear_dust_frequency_array_bool_noflicker.npy')

        data_array = data_array/counts_array

        extent = (np.min(lons), np.max(lons), np.min(lats), np.max(lats))
        m = Basemap(projection='cyl', llcrnrlon=extent[0], urcrnrlon=extent[1],
                    llcrnrlat=extent[2], urcrnrlat=extent[3],
                    resolution='i')

        m.drawcoastlines(linewidth=0.5)
        m.drawcountries(linewidth=0.5)
        parallels = np.arange(10., 40, 2.)
        # labels = [left,right,top,bottom]
        m.drawparallels(parallels, labels=[False, True, True, False],
                        linewidth=0.5)
        meridians = np.arange(-20., 17., 4.)
        m.drawmeridians(meridians, labels=[True, False, False, True],
                        linewidth=0.5)

        # m.contourf(lons, lats, data_array)
        discrete_cmap = utilities.cmap_discretize(cm.YlOrBr, 10)
        m.imshow(data_array, extent=extent, origin='lower',
                 interpolation='none',
                 cmap=discrete_cmap, vmin=0, vmax=10000000000)

        cbar = plt.colorbar(orientation='horizontal', fraction=0.056, pad=0.06)
        plt.tight_layout()
        plt.savefig(title, bbox_inches='tight')

        plt.close()

        np.save('multiyear_dust_area_array_bool_noflicker', data_array)

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

def plot_midas_comparison(pre_calculated=False):
    """
    Plots a multiyear comparison of MIDAS 0900UTC wind speed data and Tidihelt
    depression dust activation, including distributions of wind speed for
    LLJ and no-LLJ and a scatterplot of yearly mean wind speed against
    activation frequency
    dust source activation
    :param midas_filename:
    :param plume_archives:
    :return:
    """

    midas_filename = '17667.60630.SYNOP.midas.gl.nc'

    insalah = Dataset(midas_filename)

    pa2004 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                         'plume_archive_v3_prob_v3_2004')
    pa2005 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                         'plume_archive_v3_prob_v3_2005')
    pa2006 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                         'plume_archive_v3_prob_v3_2006')
    pa2007 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                         'plume_archive_v3_prob_v3_2007')
    pa2008 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                         'plume_archive_v3_prob_v3_2008')
    pa2009 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                         'plume_archive_v3_prob_v3_2009')
    pa2010 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                         'plume_archive_v3_prob_v3_2010')
    pa2011 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                         'plume_archive_v3_prob_v3_2011')
    pa2012 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                         'plume_archive_v3_prob_v3_2012')
    plume_archives = [pa2004, pa2005, pa2006, pa2007, pa2008, pa2009,
                      pa2010, pa2011, pa2012]

    emission_times = []

    # Go through every plume in the archive and take out its emission lat
    # and lon

    total_plumes = np.sum([len(plume_archive) for plume_archive in
                           plume_archives])

    archive_count = 1

    if pre_calculated:
        emission_times = np.load('tidihelt_emission_dates.npy')
        # Get 0900UTC In Salah values and their datetime
        insalah_wspd = insalah.variables['WIND_SPEED'][:]
        insalah_times = num2date(insalah.variables['TIME'][:],
                                 insalah.variables['TIME'].units)
        insalah_wdir = insalah.variables['WIND_DIRECTION'][:]

        bool_09 = np.asarray(
            [(j.hour == 9) and (j.year >= 2004 and j.year <=
                                               2012
                              and j.month in [6, 7, 8]) for j in
             insalah_times])

        insalah_wspd_09 = insalah_wspd[bool_09]
        insalah_times_09 = insalah_times[bool_09]
        insalah_wdir_09 = insalah_wdir[bool_09]

        #Include 1200UTC values - set these to 0900 so they are matched
        #with emission events
        #for i in np.arange(0, len(insalah_times_09)):
        #    if insalah_times_09[i].hour == 12:
        #        insalah_times_09[i] = insalah_times_09[i].\
        #            replace(hour=9, minute=0)

        # Go through every day in the 0900UTC values, and try to find a dust
        # activation in the dictionary. If one does not exist, we assign that
        # day as a zero, otherwise we put it in the one bin
        event_wspds = []
        no_event_wspds = []

        for i in np.arange(0, len(insalah_times_09)):
            if insalah_wspd_09[i] >= 20:
                print insalah_times_09[i]
                if insalah_times_09[i] in emission_times:
                    print 'Flagged as dusty'
            if np.ma.is_masked(insalah_wspd_09[i]):
                continue
            if insalah_times_09[i] in emission_times:
                event_wspds.append(insalah_wspd_09[i])
            else:
                no_event_wspds.append(insalah_wspd_09[i])

        # NEerly only
        """
        for i in np.arange(0, len(insalah_times_09)):
            if np.ma.is_masked(insalah_wspd_09[i]):
                continue
            if insalah_times_09[i] in emission_times and insalah_wdir_09[i] \
                    > 0 \
                    and insalah_wdir_09[i] < 90:
                event_wspds.append(insalah_wspd_09[i])
            elif insalah_wdir_09[i] > 0 and insalah_wdir_09[i] < 90:
                no_event_wspds.append(insalah_wspd_09[i])
        """

        bins = np.arange(0, 45, 5)

        # Plot two distributions, LLJ and no LLJ
        fig, ax = plt.subplots()
        ax.hist(no_event_wspds, bins=bins, color='lightcoral', label='No dust '
                                                                   'activation '
                                                                   'event',
                alpha=0.5)
        ax.hist(event_wspds, bins=bins, color='cornflowerblue', label='Dust '
                                                                'activation '
                                                                    'event',
                alpha=0.5)
        plt.legend()
        ax.set_xlabel('0900UTC In Salah MIDAS wind speed')
        ax.set_ylabel('Frequency')
        plt.savefig('midas_comparison_distributions_NEerly.png')

    else:

        for plume_archive in plume_archives:
            plume_idx = 0
            archive_size = len(plume_archive)
            used_percentages = []
            print 'Plume archive', archive_count
            archive_count += 1
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
                #emission_lat = plume_archive[i].track_centroid_lat[0]
                #emission_lon = plume_archive[i].track_centroid_lon[0]
                emission_lat = np.round(plume_archive[i].track_centroid_lat[0], 0)
                emission_lon = np.round(plume_archive[i].track_centroid_lon[0], 0)
                if (emission_lat == 25 or emission_lat == 26 or emission_lat\
                        == 27) and (
                                emission_lon == -1 or emission_lon == 0 or
                                emission_lon == 1 or emission_lon == 2):
                    emission_time = plume_archive[i].emission_time
                    # Only allow emission times which are close to 0900UTC
                    if emission_time.hour >= 8 and emission_time.hour <= 12:
                        emission_time = emission_time.replace(hour=9, minute=0)
                        emission_times.append(emission_time)
                #print (float(plumes_done)/float(total_plumes))*100

        np.save('tidihelt_emission_dates', emission_times)

        # Get 0900UTC In Salah values and their datetime
        insalah_wspd = insalah.variables['WIND_SPEED'][:]
        insalah_times = num2date(insalah.variables['TIME'][:],
                                 insalah.variables['TIME'].units)

        bool_09 = np.asarray([j.hour == 9 and (j.year >= 2004 and j.year <= 2010
                              and j.month in [6,7,8]) for j in insalah_times])

        insalah_wspd_09 = insalah_wspd[bool_09]
        insalah_times_09 = insalah_times[bool_09]

        # Go through every day in the 0900UTC values, and try to find a dust
        # activation in the dictionary. If one does not exist, we assign that
        # day as a zero, otherwise we put it in the one bin
        event_wspds = []
        no_event_wspds = []

        for i in np.arange(0, len(insalah_times_09)):
            if np.ma.is_masked(insalah_wspd_09[i]):
                continue
            if insalah_times_09[i] in emission_times:
                event_wspds.append(insalah_wspd_09[i])
            else:
                no_event_wspds.append(insalah_wspd_09[i])

        # Plot two distributions, LLJ and no LLJ
        fig, ax = plt.subplots()
        bins = np.arange(0, 45, 5)
        ax.hist(event_wspds, bins=bins, color='cornflowerblue', label='Dust '
                                                                    'activation '
                                                                    'event',
                alpha=0.5)
        ax.hist(no_event_wspds, bins=bins, color='lightcoral', label='No dust '
                                                                   'activation '
                                                                   'event',
                alpha=0.5)
        plt.legend()
        ax.set_xlabel('0900UTC In Salah MIDAS wind speed')
        ax.set_ylabel('Frequency')
        plt.savefig('midas_comparison_distributions.png')

def plot_insalah_yearly_timeseries():

    midas_filename = '17667.60630.SYNOP.midas.gl.nc'
    insalah = Dataset(midas_filename)
    insalah_times = num2date(insalah.variables['TIME'][:],
                             insalah.variables['TIME'].units)
    tidihelt_emission_dates = np.load('tidihelt_emission_dates.npy')

    pa2004 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                         'plume_archive_v3_prob_v3_2004')
    pa2005 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                         'plume_archive_v3_prob_v3_2005')
    pa2006 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                         'plume_archive_v3_prob_v3_2006')
    pa2007 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                         'plume_archive_v3_prob_v3_2007')
    pa2008 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                         'plume_archive_v3_prob_v3_2008')
    pa2009 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                         'plume_archive_v3_prob_v3_2009')
    pa2010 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                         'plume_archive_v3_prob_v3_2010')
    pa2011 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                         'plume_archive_v3_prob_v3_2011')
    pa2012 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                         'plume_archive_v3_prob_v3_2012')
    plume_archives = [pa2004, pa2005, pa2006, pa2007, pa2008, pa2009,
                      pa2010, pa2011, pa2012]

    years = np.arange(2004, 2013)

    for i in np.arange(0, len(years)):
        year_bool = np.asarray([j.year == years[i] and j.month in [6, 7,
                                                                   8] and
                                j.hour == 9
                                                                    for j in
                                insalah_times])
        year_emission_bool = np.asarray([j.year == years[i] and j.month in [6, 7,
                                                                   8] and
                                j.hour == 9
                                                                    for j in
                                tidihelt_emission_dates])
        year_emissions = tidihelt_emission_dates[year_emission_bool]
        year_times = insalah_times[year_bool]
        year_data = insalah.variables['WIND_SPEED'][year_bool]

        fig, ax = plt.subplots(figsize=(20, 5))

        plt.plot(year_times, year_data)

        ax.xaxis.set_major_locator(mdates.DayLocator(
            interval=2))  # to get a tick every day
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%d/%m/%Y'))

        plt.xticks(rotation=45, ha='right')

        for emission in year_emissions:
            plt.axvline(x=emission, color='r')

        plt.tight_layout()

        plt.savefig('Insalah_wind_speed_'+str(years[i])+'.png')

        plt.close()

def plot_tidihelt_spatial_midas(res=20, pre_calculated=False):
    """
    For the region of the Tidihelt Depression, plots a high resolution dust
    source map with pixels coloured according to the wind speed reading at
    0900UTC at In Salah. I want this to look super cool, so really high res
    with some nifty orography plotted. We can also use it for just a
    straight dust emission map
    :param midas_filename:
    :param plume_archives:
    :return:
    """

    tidihelt_emission_dates = np.load('tidihelt_emission_dates.npy')
    for i in np.arange(0, len(tidihelt_emission_dates)):
        emission_date = tidihelt_emission_dates[i]
    midas_filename = '17667.60630.SYNOP.midas.gl.nc'
    insalah = Dataset(midas_filename)

    # Get lats and lons for the extremely generous Tidihelt region and get a
    # linearly spaced 10 by 10 regular grid
    linlats = np.linspace(25, 28, res)
    linlons = np.linspace(-1, 3, res)

    lons, lats = np.meshgrid(linlons, linlats)
    data_array = np.zeros((lats.shape))
    counts_array = np.zeros((lats.shape))
    areas_array = np.zeros((lats.shape))

    insalah_wspd = insalah.variables['WIND_SPEED'][:]
    insalah_times = num2date(insalah.variables['TIME'][:],
                             insalah.variables['TIME'].units)
    insalah_wdir = insalah.variables['WIND_DIRECTION'][:]


    # Lines to constrain to northeasterly
    neerly_bool = np.asarray([j >= 0 and j <= 90 for j in insalah_wdir])
    neerly_bool = neerly_bool == 1
    insalah_wspd = insalah_wspd[neerly_bool]
    insalah_times = insalah_times[neerly_bool]

    archive_count = 1

    bool_period = np.asarray(
        [(j.year >= 2004 and j.year <=
                            2012
                            and j.month in [6, 7, 8]) for j in
         insalah_times])
    insalah_wspd = insalah_wspd[bool_period]
    insalah_times = insalah_times[bool_period]
    emission_wspds = []
    no_emission_wspds = insalah_wspd
    emission_times = []

    if pre_calculated:

        data_array = np.load('tidihelt_totals.npy')
        counts_array = np.load('tidihelt_counts.npy')

        data_array_mean = data_array / counts_array
        #data_array_mean[counts_array < 5] = np.nan

        plt.close()

        extent = (np.min(lons), np.max(lons), np.min(lats), np.max(lats))
        m = Basemap(projection='cyl', llcrnrlon=extent[0], urcrnrlon=extent[1],
                    llcrnrlat=extent[2], urcrnrlat=extent[3],
                    resolution='i')

        m.drawcoastlines(linewidth=0.5)
        m.drawcountries(linewidth=0.5)
        parallels = np.arange(25, 28, 0.5)
        # labels = [left,right,top,bottom]
        m.drawparallels(parallels, labels=[False, True, True, False],
                        linewidth=0.5)
        meridians = np.arange(-1, 3, 0.5)
        m.drawmeridians(meridians, labels=[True, False, False, True],
                        linewidth=0.5)

        # m.contourf(lons, lats, data_array)
        discrete_cmap = utilities.cmap_discretize(cm.YlOrBr, 10)
        m.imshow(counts_array, extent=extent, origin='lower',
                 interpolation='none',
                 cmap=discrete_cmap)

        plt.colorbar(orientation='horizontal', fraction=0.056, pad=0.06)

        m.scatter(insalah.variables['LON'][:], insalah.variables['LAT'][:])

        m.contour(lons, lats, orog_regridded)

        plt.tight_layout()
        plt.savefig('Tidihelt_midas_counts.png', bbox_inches='tight')

        plt.close()

    else:

        pa2004 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                             'plume_archive_v3_prob_v3_2004')
        pa2005 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                             'plume_archive_v3_prob_v3_2005')
        pa2006 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                             'plume_archive_v3_prob_v3_2006')
        pa2007 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                             'plume_archive_v3_prob_v3_2007')
        pa2008 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                             'plume_archive_v3_prob_v3_2008')
        pa2009 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                             'plume_archive_v3_prob_v3_2009')
        pa2010 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                             'plume_archive_v3_prob_v3_2010')
        pa2011 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                             'plume_archive_v3_prob_v3_2011')
        pa2012 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                             'plume_archive_v3_prob_v3_2012')
        plume_archives = [pa2004, pa2005, pa2006, pa2007, pa2008, pa2009,
                          pa2010, pa2011, pa2012]

        for plume_archive in plume_archives:
            print 'Plume archive', archive_count
            archive_count += 1
            plume_idx = 0
            archive_size = len(plume_archive)
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
                plume = plume_archive[i]
                emission_lat = plume.track_centroid_lat[0]
                emission_lon = plume.track_centroid_lon[0]
                if np.round(emission_lat, 0) in [25, 26, 27, 28
                                                 ] and np.round(
                    emission_lon, 0) in [-1, 0, 1, 2, 3]:
                    emission_time = plume.emission_time
                    # Find the nearest emission time from MIDAS 3-hourly
                    available_hours = [0, 3, 6, 9, 12, 15, 18, 21]
                    possible_times = [plume.emission_time.
                                          replace(hour=available_hours[k],
                                                  minute=0) for
                                      k in
                                                np.arange(0,
                                                          len(available_hours))]
                    time_differences = [abs(emission_time-j).total_seconds() for j
                                        in
                                        possible_times]
                    nearest_time = np.asarray(possible_times)[
                        time_differences==np.min(
                        time_differences)]

                    emission_area = plume.track_area[0]

                    emission_time = nearest_time[0]
                    time_bool = insalah_times == emission_time
                    if np.all(time_bool == False):
                        continue
                    insalah_wspd_matched = insalah_wspd[time_bool][0]
                    if np.ma.is_masked(insalah_wspd_matched):
                        continue
                    emission_wspds.append(insalah_wspd_matched)
                    emission_times.append(emission_time)
                    diff_array_lat = abs(lats - emission_lat)
                    diff_array_lon = abs(lons - emission_lon)
                    latbool = diff_array_lat == np.min(diff_array_lat)
                    lonbool = diff_array_lon == np.min(diff_array_lon)
                    union_bool = latbool & lonbool
                    data_array[union_bool] += insalah_wspd_matched
                    counts_array[union_bool] += 1
                    areas_array[union_bool] += emission_area

        remove_bool = np.asarray([j in emission_times for j in insalah_times])
        no_emission_wspds = np.asarray(no_emission_wspds)
        no_emission_wspds[remove_bool] = -999

        areas_array = areas_array/counts_array

        no_emission_wspds = no_emission_wspds[no_emission_wspds!=-999]

        #emission_wspds = emission_wspds[emission_wspds < 100]
        no_emission_wspds = no_emission_wspds[no_emission_wspds < 100]

        np.save('tidihelt_totals', data_array)
        np.save('tidihelt_counts', counts_array)
        np.save('tidihelt_areas', areas_array)

        data_array_mean = data_array/counts_array

        plt.close()

        extent = (np.min(lons), np.max(lons), np.min(lats), np.max(lats))
        m = Basemap(projection='cyl', llcrnrlon=extent[0], urcrnrlon=extent[1],
                    llcrnrlat=extent[2], urcrnrlat=extent[3],
                    resolution='i')

        m.drawcoastlines(linewidth=0.5)
        m.drawcountries(linewidth=0.5)
        parallels = np.arange(25, 28, 0.5)
        # labels = [left,right,top,bottom]
        m.drawparallels(parallels, labels=[False, True, True, False],
                        linewidth=0.5)
        meridians = np.arange(-1, 3, 0.5)
        m.drawmeridians(meridians, labels=[True, False, False, True],
                        linewidth=0.5)

        # m.contourf(lons, lats, data_array)
        discrete_cmap = utilities.cmap_discretize(cm.YlOrBr, 10)
        m.imshow(data_array_mean, extent=extent, origin='lower',
                 interpolation='none',
                 cmap=discrete_cmap)

        plt.colorbar(orientation='horizontal', fraction=0.056, pad=0.06)
        m.scatter(insalah.variables['LON'][:], insalah.variables['LAT'][:])
        plt.tight_layout()
        plt.savefig('Tidihelt_midas_wspd_time_matched_northeasterly.png',
                    bbox_inches='tight')

        plt.close()

        extent = (np.min(lons), np.max(lons), np.min(lats), np.max(lats))
        m = Basemap(projection='cyl', llcrnrlon=extent[0], urcrnrlon=extent[1],
                    llcrnrlat=extent[2], urcrnrlat=extent[3],
                    resolution='i')

        m.drawcoastlines(linewidth=0.5)
        m.drawcountries(linewidth=0.5)
        parallels = np.arange(25, 28, 0.5)
        # labels = [left,right,top,bottom]
        m.drawparallels(parallels, labels=[False, True, True, False],
                        linewidth=0.5)
        meridians = np.arange(-1, 3, 0.5)
        m.drawmeridians(meridians, labels=[True, False, False, True],
                        linewidth=0.5)

        discrete_cmap = utilities.cmap_discretize(cm.YlOrBr, 10)
        m.imshow(areas_array, extent=extent, origin='lower',
                 interpolation='none',
                 cmap=discrete_cmap)

        plt.colorbar(orientation='horizontal', fraction=0.056, pad=0.06)
        plt.tight_layout()
        plt.savefig('Tidihelt_areas.png',
                    bbox_inches='tight')

        plt.close()

        fig, ax = plt.subplots()
        bins = np.arange(0, 45, 5)
        ax.hist(emission_wspds, bins=bins, color='cornflowerblue',
                label='Dust '
                                                                      'activation '
                                                                      'event',
                alpha=0.5)
        ax.hist(no_emission_wspds, bins=bins, color='lightcoral', label='No '
                                                                      'dust '
                                                                     'activation '
                                                                     'event',
                alpha=0.5)
        plt.legend()
        ax.set_xlabel('0900UTC In Salah MIDAS wind speed')
        ax.set_ylabel('Frequency')
        plt.savefig('midas_comparison_distributions_northeasterly.png')

        plt.close()

def plot_midas_correlation_maps(res=100):
    """
    For each selected station in the MIDAS database, plots the correlation
    between surface wind speed at that station and dust emission from SEVIRI
    on a 100 by 100 grid
    :return:
    """

    # Read in MIDAS stations
    timimoun = Dataset('18225.60607.SYNOP.midas.gl.nc')
    adrar = Dataset('17666.60620.SYNOP.midas.gl.nc')
    insalah = Dataset('17667.60630.SYNOP.midas.gl.nc')
    beni_abbes = Dataset('18224.60602.SYNOP.midas.gl.nc')
    inguezzam = Dataset('20159.60690.SYNOP.midas.gl.nc')
    bbm = Dataset('23891.60686.SYNOP.midas.gl.nc')
    rhourd_nouss = Dataset('56367.60603.SYNOP.midas.gl.nc')
    elgolea = Dataset('17665.60590.SYNOP.midas.gl.nc')
    bechar = Dataset('17664.60571.SYNOP.midas.gl.nc')
    tindouf = Dataset('18227.60656.SYNOP.midas.gl.nc')
    tamanrasset = Dataset('17668.60680.SYNOP.midas.gl.nc')
    illizi = Dataset('20188.60640.SYNOP.midas.gl.nc')
    tessalit = Dataset('18239.61202.SYNOP.midas.gl.nc')

    labels = ['insalah', 'adrar', 'timimoun', 'beni_abbes', 'inguezzam', 'bbm',
              'rhourd_nouss', 'elgolea', 'bechar', 'tindouf', 'tamanrasset',
              'illizi', 'tessalit']

    midas_stations = [insalah, adrar, timimoun, beni_abbes, inguezzam, bbm,
                      rhourd_nouss, elgolea, bechar, tindouf, tamanrasset,
                      illizi, tessalit]

    # Get MIDAS datetimes from any one of the MIDAS stations
    insalah_times = num2date(insalah.variables['TIME'][:],
                             insalah.variables['TIME'].units)

    # Constrain the dates to the period for which we have plumetracker
    bool_period = np.asarray(
        [(j.year >= 2004 and j.year <=
          2017
          and j.month in [6, 7, 8]) for j in
         insalah_times])

    midas_dates = insalah_times[bool_period]

    # Get lats and lons
    sdf_test = Dataset('/soge-home/data_not_backed_up/satellite/meteosat/seviri'
                       '/15-min/0.03x0.03/sdf/nc/JUNE2010/SDF_v2/SDF_v2.'
                       '201006031500.nc')

    lons, lats = np.meshgrid(sdf_test.variables['longitude'][:],
                             sdf_test.variables['latitude'][:])

    linlats = np.linspace(np.min(lats), np.max(lats), res)
    linlons = np.linspace(np.min(lons), np.max(lons), res)

    lons, lats = np.meshgrid(linlons, linlats)

    # Generate a massive lat lon array
    plume_emission_array = np.zeros((res, res, len(midas_dates)))

    # Get plume archives
    pa2004 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                         'plume_archive_flicker_v3_prob_v4_2004')
    pa2005 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                         'plume_archive_flicker_v3_prob_v4_2005')
    pa2006 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                         'plume_archive_flicker_v3_prob_v4_2006')
    pa2007 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                         'plume_archive_flicker_v3_prob_v4_2007')
    pa2008 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                         'plume_archive_flicker_v3_prob_v4_2008')
    pa2009 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                         'plume_archive_flicker_v3_prob_v4_2009')
    pa2010 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                         'plume_archive_flicker_v3_prob_v4_2010')
    pa2011 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                         'plume_archive_flicker_v3_prob_v3_2011')
    pa2012 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                         'plume_archive_flicker_v3_prob_v3_2012')
    pa2013 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                         'plume_archive_flicker_v3_prob_v4_2013')
    pa2014 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                         'plume_archive_flicker_v3_prob_v4_2014')
    pa2015 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                         'plume_archive_flicker_v3_prob_v4_2015')
    pa2016 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                         'plume_archive_flicker_v3_prob_v4_2016')
    pa2017 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                         'plume_archive_flicker_v3_prob_v4_2017')
    plume_archives = [pa2004, pa2005, pa2006, pa2007, pa2008, pa2009,
                      pa2010, pa2011, pa2012, pa2013, pa2014, pa2015,
                      pa2016, pa2017]

    # Convert an array of datetimes corresponding to the plume archives into
    #  an index array
    """
    archive_count = 1

    for plume_archive in plume_archives:
        print 'Plume archive', archive_count
        archive_count += 1
        plume_idx = 0
        archive_size = len(plume_archive)
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

            # Extract the emission date of the plume
            plume = plume_archive[i]
            emission_time = plume.emission_time
            emission_lat = plume.track_centroid_lat[0]
            emission_lon = plume.track_centroid_lon[0]

            # Now need to assign this to the right time
            # Get the date then assign it to the nearest three hourly timeslot
            # Find the index in the dates array which matches that
            available_hours = [0, 3, 6, 9, 12, 15, 18, 21]
            possible_times = [plume.emission_time.
                                  replace(hour=available_hours[k],
                                          minute=0) for
                              k in
                              np.arange(0,
                                        len(available_hours))]
            time_differences = [abs(emission_time - j).total_seconds() for j
                                in
                                possible_times]
            nearest_time = np.asarray(possible_times)[
                time_differences == np.min(
                    time_differences)]

            # Assign that to the correct array index
            time_bool = midas_dates == nearest_time[0]
            diff_array_lat = abs(lats - emission_lat)
            diff_array_lon = abs(lons - emission_lon)
            latbool = diff_array_lat == np.min(diff_array_lat)
            lonbool = diff_array_lon == np.min(diff_array_lon)
            union_bool = latbool & lonbool

            # Finally, assign a 1 to this point in time and space
            if True in union_bool and True in time_bool:
                plume_emission_array[union_bool, time_bool] += 1

    # Save the final array to file
    np.save('plume_emission_array', plume_emission_array)

    """

    plume_emission_array = np.load('plume_emission_array.npy')

    print 'Finished processing tracked plumes. Calculating correlations.'

    # Loop through each MIDAS station
    for i in np.arange(0, len(midas_stations)):
        # Extract wind speed
        wspd = midas_stations[i].variables['WIND_SPEED'][:]
        wdir = midas_stations[i].variables['WIND_DIRECTION'][:]
        # Convert wind direction to offset from 0 degrees
        wdir_offsets = np.asarray([abs(utilities.get_distance(0, j)) for j in \
                wdir])
        rh = np.asarray(utilities.dpd_to_rh(midas_stations[i]))
        temp = midas_stations[i].variables['AIR_TEMPERATURE'][:]
        times = num2date(midas_stations[i].variables['TIME'][:],
                         midas_stations[i].variables['TIME'].units)

        # So we want to constrain the times of the midas wind speeds to when
        # SEVIRI plume emissions are available

        # We also want to constrain SEVIRI plume emissions to when midas
        # wind speeds are available

        # So one boolean to select only midas wind speeds which are in
        # SEVIRI times

        # So then our midas array could well be smaller than SEVIRI, but the
        #  second removes from SEVIRI according to what's not in midas,
        # so the two booleans should finally be the same shape

        print 'Matching up the times between the two datasets.'

        midas_bool_period = np.asarray(
            [j in midas_dates for j in
             times])
        times_constrained = times[midas_bool_period]
        seviri_wspds = wspd[midas_bool_period]
        seviri_wdir = wdir_offsets[midas_bool_period]
        seviri_rh = rh[midas_bool_period]
        seviri_temp = temp[midas_bool_period]

        seviri_bool_period = np.asarray(
            [j in times_constrained for j in
             midas_dates])
        plume_emission_array_constrained = plume_emission_array[:, :,
                                           seviri_bool_period]

        # Ideally get the overlap between the two date bools

        # Constrain SEVIRI activation data to the period for which MIDAS
        # data is available

        # Get an array of wind speeds out for the period

        correlation_array = np.zeros((res, res))

        station_lat = midas_stations[i].variables['LAT'][:]
        station_lon = midas_stations[i].variables['LON'][:]

        emission_wdirs = []
        no_emission_wdirs = []

        emission_rhs = []
        no_emission_rhs = []

        emission_temps = []
        no_emission_temps = []

        emission_wspds = []
        no_emission_wspds = []

        # Correlate the wind speed with every box in the grid
        for idx, val in np.ndenumerate(plume_emission_array[:, :, 0]):
            emission_array = plume_emission_array_constrained[idx]
            correlation = np.ma.corrcoef(seviri_wspds, emission_array)
            #if np.sum(emission_array) < 5:
            #    correlation_array[idx] = np.nan
            #else:
            correlation_array[idx] = correlation[0, 1]

            distance = utilities.haversine(lons[idx], lats[idx],
                                           station_lon, station_lat)

            if distance < 300:
                if correlation_array[idx] > 0.05:
                    emission_array_bool = emission_array == 1
                    wdir_emissions = seviri_wdir[emission_array_bool]
                    emission_wdirs.append(wdir_emissions)
                    wdir_no_emissions = seviri_wdir[~emission_array_bool]
                    no_emission_wdirs.append(wdir_no_emissions)
                    rh_emissions = seviri_rh[emission_array_bool]
                    emission_rhs.append(rh_emissions)
                    rh_no_emissions = seviri_rh[~emission_array_bool]
                    no_emission_rhs.append(rh_no_emissions)
                    temp_emissions = seviri_temp[emission_array_bool]
                    emission_temps.append(temp_emissions)
                    temp_no_emissions = seviri_temp[~emission_array_bool]
                    no_emission_temps.append(temp_no_emissions)
                    wspd_emissions = seviri_wspds[emission_array_bool]
                    emission_wspds.append(wspd_emissions)
                    wspd_no_emissions = seviri_wspds[~emission_array_bool]
                    no_emission_wspds.append(wspd_no_emissions)

        np.save('correlation_array_'+labels[i], correlation_array)

        emission_wdirs = np.concatenate(emission_wdirs).ravel()
        no_emission_wdirs = np.concatenate(no_emission_wdirs).ravel()

        emission_rhs = np.concatenate(emission_rhs).ravel()
        no_emission_rhs = np.concatenate(no_emission_rhs).ravel()

        emission_temps = np.concatenate(emission_temps).ravel()
        no_emission_temps = np.concatenate(no_emission_temps).ravel()

        emission_wspds = np.concatenate(emission_wspds).ravel()
        no_emission_wspds = np.concatenate(no_emission_wspds).ravel()

        emission_wdirs = emission_wdirs[np.isfinite(emission_wdirs)]
        no_emission_wdirs = no_emission_wdirs[np.isfinite(no_emission_wdirs)]
        emission_rhs = emission_rhs[np.isfinite(emission_rhs)]
        no_emission_rhs = no_emission_rhs[np.isfinite(no_emission_rhs)]
        emission_temps = emission_temps[np.isfinite(emission_temps)]
        no_emission_temps = no_emission_temps[np.isfinite(no_emission_temps)]
        emission_wspds = emission_wspds[np.isfinite(emission_wspds)]
        no_emission_wspds = no_emission_wspds[np.isfinite(no_emission_wspds)]

        plt.close()

        # Plot the correlation array along with the position of the station
        extent = (np.min(lons), np.max(lons), np.min(lats), np.max(lats))
        m = Basemap(projection='cyl', llcrnrlon=extent[0], urcrnrlon=extent[1],
                    llcrnrlat=extent[2], urcrnrlat=extent[3],
                    resolution='i')

        m.drawcoastlines(linewidth=0.5)
        m.drawcountries(linewidth=0.5)
        parallels = np.arange(10., 40, 2.)
        # labels = [left,right,top,bottom]
        m.drawparallels(parallels, labels=[False, True, True, False],
                        linewidth=0.5)
        meridians = np.arange(-20., 17., 4.)
        m.drawmeridians(meridians, labels=[True, False, False, True],
                        linewidth=0.5)

        correlation_array = np.ma.masked_where(np.isnan(correlation_array),
                                               correlation_array)

        max = 0.1
        min = -0.1

        levels = MaxNLocator(nbins=15).tick_values(min, max)

        #discrete_cmap = utilities.cmap_discretize(cm.RdYlBu_r, 10)
        cmap = cm.get_cmap('RdBu_r')
        norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)

        #m.pcolormesh(lons, lats, correlation_array, extent=extent,
        #             origin='lower',
        #         interpolation='none',
        #         cmap=cmap, vmin=min, vmax=max+1)

        m.pcolormesh(lons, lats, correlation_array,
                     cmap=cmap, norm=norm, vmin=min, vmax=max)

        cbar = plt.colorbar(orientation='horizontal', fraction=0.056, pad=0.06)
        cbar.ax.set_xlabel('Correlation')
        m.scatter(midas_stations[i].variables['LON'][:],
                  midas_stations[i].variables['LAT'][:])
        plt.tight_layout()
        plt.savefig('MIDAS_dust_emission_correlation_'+labels[
            i]+'_2004_2017.png',
                    bbox_inches='tight')

        plt.close()

        print len(emission_wspds)

        fig, ax = plt.subplots()
        bins = np.arange(0, 45, 5)
        ax.hist(no_emission_wspds, bins=bins, color='lightcoral',
                label='No dust activation event',
                alpha=0.5)
        ax2 = ax.twinx()
        ax2.hist(emission_wspds, bins=bins, color='cornflowerblue',
                label='Dust '
                      'activation '
                      'event',
                alpha=0.5)
        ax.legend(loc=2)
        ax2.legend(loc=1)
        #plt.legend()
        ax.set_xlabel('Wind speed (m/s)')
        ax.set_ylabel('Frequency (no dust)')
        ax2.set_ylabel('Frequency (dust)')
        plt.savefig('wspd_distributions_2004_2017_'+labels[i]+'.png')

        plt.close()

        print len(emission_wdirs)

        fig, ax = plt.subplots()
        bins = np.arange(0, 190, 10)
        ax.hist(no_emission_wdirs, bins=bins, color='lightcoral',
                label='No dust activation event',
                alpha=0.5)
        ax2 = ax.twinx()
        ax2.hist(emission_wdirs, bins=bins, color='cornflowerblue',
                label='Dust '
                      'activation '
                      'event',
                alpha=0.5)
        ax.legend()
        ax2.legend(loc=0)
        ax.set_xlabel('Wind direction (deg)')
        ax.set_ylabel('Frequency (no dust)')
        ax2.set_ylabel('Frequency (dust)')
        plt.savefig('wdir_distributions_2004_2017_' + labels[i] + '.png')

        plt.close()

        print len(emission_rhs)

        fig, ax = plt.subplots()
        bins = np.arange(0, 110, 10)
        ax.hist(no_emission_rhs, bins=bins, color='lightcoral',
                label='No dust activation event',
                alpha=0.5)
        ax2 = ax.twinx()
        ax2.hist(emission_rhs, bins=bins, color='cornflowerblue',
                label='Dust '
                      'activation '
                      'event',
                alpha=0.5)
        ax.legend()
        ax2.legend(loc=0)
        ax.set_xlabel('RH (%)')
        ax.set_ylabel('Frequency (no dust)')
        ax2.set_ylabel('Frequency (dust)')
        plt.savefig('rh_distributions_' + labels[i] + '.png')

        plt.close()

        print len(emission_temps)

        fig, ax = plt.subplots()
        bins = np.arange(0, 55, 5)
        ax.hist(no_emission_temps, bins=bins, color='lightcoral',
                label='No dust activation event',
                alpha=0.5)
        ax2 = ax.twinx()
        ax2.hist(emission_temps, bins=bins, color='cornflowerblue',
                label='Dust '
                      'activation '
                      'event',
                alpha=0.5)
        ax.legend()
        ax2.legend(loc=0)
        ax.set_xlabel('Temp (degC)')
        ax.set_ylabel('Frequency (no dust)')
        ax2.set_ylabel('Frequency (dust)')
        plt.savefig('temp_distributions_2004_2017' + labels[i] + '.png')

        plt.close()

def plume_midas_eda():
    """
    Plots a pandas EDA for dust emission events matched to the nearest MIDAS
    timestep
    :return:
    """

def plume_eda(pre_calculated=False):
    """
    Plots a pandas EDA for plume attributes
    :return:
    """

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

    total_plumes = np.sum([len(j) for j in plume_archives])

    if pre_calculated == False:

        # Create an array with columns for the number of plumes and rows for the
        #  number of attributes
        plume_arr = np.zeros((9, total_plumes))

        archive_idx = 1
        plume_count = 0

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
                plume = plume_archive[i]
                conv_distance = plume.conv_distance
                plume_arr[0, plume_count] = conv_distance

                emission_time = plume.dates_observed[0]
                # Take the smallest time to the nearest 0900UTC
                nine = emission_time.replace(hour=9, minute=0, second=0)
                nine_1 = nine + dt.timedelta(days=1)
                distance_nine = abs(emission_time - nine)
                distance_nine_1 = abs(nine_1 - emission_time)
                distance_from_09 = np.min([distance_nine.total_seconds(),
                                           distance_nine_1.total_seconds()])
                plume_arr[1, plume_count] = distance_from_09

                mean_axis_offset = plume.mean_axis_offset
                plume_arr[2, plume_count] = mean_axis_offset

                mean_centroid_speed = np.nanmean(plume.track_speed_centroid)
                plume_arr[3, plume_count] = mean_centroid_speed

                try:
                    mean_gpe_speed = plume.speed_gpe
                    plume_arr[4, plume_count] = mean_gpe_speed
                except:
                    plume_arr[4, plume_count] = np.nan

                maximum_extent = np.max(plume.track_area)
                plume_arr[5, plume_count] = maximum_extent

                duration = plume.duration.total_seconds()
                plume_arr[6, plume_count] = duration

                llj_prob = plume.LLJ_prob
                plume_arr[7, plume_count] = llj_prob

                # Get the direction from the total distance travelled by the plume
                first_position_lat = plume.track_centroid_lat[0]
                first_position_lon = plume.track_centroid_lon[0]
                last_position_lat = plume.track_centroid_lat[-1]
                last_position_lon = plume.track_centroid_lon[-1]

                # Get the direction
                final_direction = utilities.calculate_initial_compass_bearing((
                    first_position_lat, first_position_lon),
                    (last_position_lat,
                     last_position_lon))
                plume_arr[8, plume_count] = final_direction

                plume_count += 1

        np.save('plume_attributes', plume_arr)
        pre_calculated = True

    if pre_calculated == True:

        plume_arr = np.load('plume_attributes.npy')

        ids = np.arange(0, total_plumes)

        colnames = ['conv_distance', 'seconds_from_09', 'mean_axis_offset',
                    'mean_centroid_speed', 'mean_gpe_speed', 'maximum_extent',
                    'duration', 'llj_prob', 'final_centroid_direction']

        df = pd.DataFrame(plume_arr.T,
                          index=ids,
                          columns=colnames)

        df = df.dropna(0)

        plt.close()
        g = seaborn.pairplot(df)
        plt.savefig('plume_pairplot.png')
        plt.close()

        # Compute the correlation matrix
        corr = df.corr()

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
        plt.savefig('plume_correlheatmap.png')
        plt.close()



    # Attributes include convection distance, distance from 09, mean axis
    # offset, mean plume centroid speed, plume GPE speed, axis direction,
    # plume extent, plume duration, LLJ probability

def plot_llj_predictors(title_1='mean_convection_distance.png',
                        title_2='mean_emission_distance_from_09.png',
                        title_3='mean_axis_direction_offset.png',
                        title_4='mean_duration', res=100, pre_calculated=False):
    """
    Plots each of the constituent predictor variables for the LLJ
    probability maps to see if there are any strong controlling factors on
    the final distribution (e.g. I expect conv distance is a big one).
    :return:
    """

    # Get lats and lons
    sdf_test = Dataset(
        '/soge-home/data_not_backed_up/satellite/meteosat/seviri'
        '/15-min/0.03x0.03/sdf/nc/JUNE2010/SDF_v2/SDF_v2.'
        '201006031500.nc')

    lons, lats = np.meshgrid(sdf_test.variables['longitude'][:],
                             sdf_test.variables['latitude'][:])

    linlats = np.linspace(np.min(lats), np.max(lats), res)
    linlons = np.linspace(np.min(lons), np.max(lons), res)

    lons, lats = np.meshgrid(linlons, linlats)

    # Get plume archives
    pa2004 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                         'plume_archive_flicker_v3_prob_v3_2004')
    pa2005 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                         'plume_archive_flicker_v3_prob_v3_2005')
    pa2006 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                         'plume_archive_flicker_v3_prob_v3_2006')
    pa2007 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                         'plume_archive_flicker_v3_prob_v3_2007')
    pa2008 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                         'plume_archive_flicker_v3_prob_v3_2008')
    pa2009 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                         'plume_archive_flicker_v3_prob_v3_2009')
    pa2010 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                         'plume_archive_flicker_v3_prob_v3_2010')
    pa2011 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                         'plume_archive_flicker_v3_prob_v3_2010')
    pa2012 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                         'plume_archive_flicker_v3_prob_v3_2010')
    plume_archives = [pa2004, pa2005, pa2006, pa2007, pa2008, pa2009,
                      pa2010, pa2011, pa2012]

    if pre_calculated:
        plt.close()
        fig, ax = plt.subplots()
        extent = (np.min(lons), np.max(lons), np.min(lats), np.max(lats))
        m = Basemap(projection='cyl', llcrnrlon=extent[0], urcrnrlon=extent[1],
                    llcrnrlat=extent[2], urcrnrlat=extent[3],
                    resolution='i')

        m.drawcoastlines(linewidth=0.5)
        m.drawcountries(linewidth=0.5)
        parallels = np.arange(10., 40, 2.)
        # labels = [left,right,top,bottom]
        m.drawparallels(parallels, labels=[False, True, True, False],
                        linewidth=0.5)
        meridians = np.arange(-20., 17., 4.)
        m.drawmeridians(meridians, labels=[True, False, False, True],
                        linewidth=0.5)

        lats = np.linspace(np.min(lats), np.max(lats), res)
        lons = np.linspace(np.min(lons), np.max(lons), res)
        lons, lats = np.meshgrid(lons, lats)
        counts_data_array = np.load('multiyear_dust_frequency_array.npy')
        conv_data_array = np.load('multiyear_dust_conv_array.npy')
        time_data_array = np.load('multiyear_dust_time_array.npy')
        offset_data_array = np.load('multiyear_dust_offset_array.npy')

        plt.close()

        extent = (np.min(lons), np.max(lons), np.min(lats), np.max(lats))
        m = Basemap(projection='cyl', llcrnrlon=extent[0], urcrnrlon=extent[1],
                    llcrnrlat=extent[2], urcrnrlat=extent[3],
                    resolution='i')

        m.drawcoastlines(linewidth=0.5)
        m.drawcountries(linewidth=0.5)
        parallels = np.arange(10., 40, 2.)
        # labels = [left,right,top,bottom]
        m.drawparallels(parallels, labels=[False, True, True, False],
                        linewidth=0.5)
        meridians = np.arange(-20., 17., 4.)
        m.drawmeridians(meridians, labels=[True, False, False, True],
                        linewidth=0.5)

        min = 0
        max = 500

        levels = MaxNLocator(nbins=15).tick_values(min, max)

        # discrete_cmap = utilities.cmap_discretize(cm.RdYlBu_r, 10)
        cmap = cm.get_cmap('Greys')
        norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)

        # m.pcolormesh(lons, lats, correlation_array, extent=extent,
        #             origin='lower',
        #         interpolation='none',
        #         cmap=cmap, vmin=min, vmax=max+1)

        conv_data_array = np.ma.masked_where(np.isnan(conv_data_array),
                                             conv_data_array)

        m.pcolormesh(lons, lats, conv_data_array,
                     cmap=cmap, norm=norm, vmin=min, vmax=max)

        cbar = plt.colorbar(orientation='horizontal', fraction=0.056, pad=0.06)
        cbar.ax.set_xlabel('Distance to convection (km)')
        plt.tight_layout()
        plt.savefig(title_1, bbox_inches='tight')

        plt.close()

        extent = (np.min(lons), np.max(lons), np.min(lats), np.max(lats))
        m = Basemap(projection='cyl', llcrnrlon=extent[0], urcrnrlon=extent[1],
                    llcrnrlat=extent[2], urcrnrlat=extent[3],
                    resolution='i')

        m.drawcoastlines(linewidth=0.5)
        m.drawcountries(linewidth=0.5)
        parallels = np.arange(10., 40, 2.)
        # labels = [left,right,top,bottom]
        m.drawparallels(parallels, labels=[False, True, True, False],
                        linewidth=0.5)
        meridians = np.arange(-20., 17., 4.)
        m.drawmeridians(meridians, labels=[True, False, False, True],
                        linewidth=0.5)

        min = 0
        max = 40000

        levels = MaxNLocator(nbins=15).tick_values(min, max)

        # discrete_cmap = utilities.cmap_discretize(cm.RdYlBu_r, 10)
        cmap = cm.get_cmap('Greys')
        norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)

        # m.pcolormesh(lons, lats, correlation_array, extent=extent,
        #             origin='lower',
        #         interpolation='none',
        #         cmap=cmap, vmin=min, vmax=max+1)

        time_data_array = np.ma.masked_where(np.isnan(time_data_array),
                                             time_data_array)

        m.pcolormesh(lons, lats, time_data_array,
                     cmap=cmap, norm=norm, vmin=min, vmax=max)

        cbar = plt.colorbar(orientation='horizontal', fraction=0.056, pad=0.06)
        cbar.ax.set_xlabel('Time to 0900UTC (seconds)')
        plt.tight_layout()
        plt.savefig(title_2, bbox_inches='tight')

        plt.close()

        extent = (np.min(lons), np.max(lons), np.min(lats), np.max(lats))
        m = Basemap(projection='cyl', llcrnrlon=extent[0], urcrnrlon=extent[1],
                    llcrnrlat=extent[2], urcrnrlat=extent[3],
                    resolution='i')

        m.drawcoastlines(linewidth=0.5)
        m.drawcountries(linewidth=0.5)
        parallels = np.arange(10., 40, 2.)
        # labels = [left,right,top,bottom]
        m.drawparallels(parallels, labels=[False, True, True, False],
                        linewidth=0.5)
        meridians = np.arange(-20., 17., 4.)
        m.drawmeridians(meridians, labels=[True, False, False, True],
                        linewidth=0.5)

        min = 0
        max = 90

        levels = MaxNLocator(nbins=15).tick_values(min, max)

        # discrete_cmap = utilities.cmap_discretize(cm.RdYlBu_r, 10)
        cmap = cm.get_cmap('Greys')
        norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)

        # m.pcolormesh(lons, lats, correlation_array, extent=extent,
        #             origin='lower',
        #         interpolation='none',
        #         cmap=cmap, vmin=min, vmax=max+1)

        offset_data_array = np.ma.masked_where(np.isnan(offset_data_array),
                                             offset_data_array)

        m.pcolormesh(lons, lats, offset_data_array,
                     cmap=cmap, norm=norm, vmin=min, vmax=max)

        cbar = plt.colorbar(orientation='horizontal', fraction=0.056, pad=0.06)
        cbar.ax.set_xlabel('Offset from major axis (degrees)')
        plt.tight_layout()
        plt.savefig(title_3, bbox_inches='tight')

        plt.close()

    else:

        plt.close()

        fig, ax = plt.subplots()

        latlon_dictionary = {}

        lats = np.linspace(np.min(lats), np.max(lats), res)
        lons = np.linspace(np.min(lons), np.max(lons), res)
        lons, lats = np.meshgrid(lons, lats)

        archive_idx = 1

        # Data array
        data_array = np.zeros((lons.shape))
        offset_counts_array = np.zeros((lons.shape))
        conv_data_array = np.zeros((lons.shape))
        time_data_array = np.zeros((lons.shape))
        offset_data_array = np.zeros((lons.shape))

        data_array_jun = np.zeros((lons.shape))
        data_array_jul = np.zeros((lons.shape))
        data_array_aug = np.zeros((lons.shape))

        offset_counts_array_jun = np.zeros((lons.shape))
        offset_counts_array_jul = np.zeros((lons.shape))
        offset_counts_array_aug = np.zeros((lons.shape))

        conv_data_array_jun = np.zeros((lons.shape))
        conv_data_array_jul = np.zeros((lons.shape))
        conv_data_array_aug = np.zeros((lons.shape))

        time_data_array_jun = np.zeros((lons.shape))
        time_data_array_jul = np.zeros((lons.shape))
        time_data_array_aug = np.zeros((lons.shape))

        offset_data_array_jun = np.zeros((lons.shape))
        offset_data_array_jul = np.zeros((lons.shape))
        offset_data_array_aug = np.zeros((lons.shape))

        for plume_archive in plume_archives:
            print 'Plume archive', archive_idx
            archive_idx += 1
            archive_size = len(plume_archive)
            plume_idx = 1
            used_percentages = []
            for i in plume_archive:
                if int((float(plume_idx)/float(
                    archive_size))*100)%10 == 0 and int((float(
                    plume_idx)/float(archive_size))*100) not in \
                        used_percentages:
                    print str(int((float(plume_idx)/float(
                    archive_size))*100))+"%"
                    # This percentage has been printed already
                    used_percentages.append(int((float(plume_idx)/float(
                    archive_size))*100))
                plume_idx += 1
                # If the plume has merged, the plume source is found in the
                # pre-merge track
                plume = plume_archive[i]
                emission_time = plume.dates_observed[0]
                emission_lat = np.round(plume.track_centroid_lat[0], 1)
                emission_lon = np.round(plume.track_centroid_lon[0], 1)
                lat_diff = np.abs(lats-emission_lat)
                lon_diff = np.abs(lons-emission_lon)
                lat_bool = lat_diff == lat_diff.min()
                lon_bool = lon_diff == lon_diff.min()
                union_bool = lat_bool & lon_bool
                data_array[union_bool] += 1

                conv_distance = plume.conv_distance
                conv_data_array[union_bool] += conv_distance

                # Take the smallest time to the nearest 0900UTC
                nine = emission_time.replace(hour=9, minute=0, second=0)
                nine_1 = nine + dt.timedelta(days=1)
                distance_nine = abs(emission_time - nine)
                distance_nine_1 = abs(nine_1 - emission_time)
                distance_from_09 = np.min([distance_nine.total_seconds(),
                                           distance_nine_1.total_seconds()])
                time_data_array[union_bool] += distance_from_09
                offset = plume.mean_axis_offset
                # Offset needs a special N, as it has None sometimes
                if offset != None:
                    offset_data_array[union_bool] += offset
                    offset_counts_array[union_bool] += 1
                    if emission_time.month == 6:
                        offset_data_array_jun[union_bool] += offset
                        offset_counts_array_jun[union_bool] += 1
                    elif emission_time.month == 7:
                        offset_data_array_jul[union_bool] += offset
                        offset_counts_array_jul[union_bool] += 1
                    elif emission_time.month == 8:
                        offset_data_array_aug[union_bool] += offset
                        offset_counts_array_aug[union_bool] += 1

                if emission_time.month == 6:
                    conv_data_array_jun[union_bool] += conv_distance
                    time_data_array_jun[union_bool] += distance_from_09
                    data_array_jun[union_bool] += 1
                elif emission_time.month == 7:
                    conv_data_array_jul[union_bool] += conv_distance
                    time_data_array_jul[union_bool] += distance_from_09
                    data_array_jul[union_bool] += 1
                elif emission_time.month == 8:
                    conv_data_array_aug[union_bool] += conv_distance
                    time_data_array_aug[union_bool] += distance_from_09
                    data_array_aug[union_bool] += 1

        # Find the nearest latlon for each dictionary key, and extract the index
        # for it

        # Get means
        data_array[data_array == 0] = np.nan
        conv_data_array = conv_data_array/data_array
        time_data_array = time_data_array/data_array
        offset_data_array = offset_data_array/offset_counts_array

        data_array_jun[data_array_jun == 0] = np.nan
        data_array_jul[data_array_jul == 0] = np.nan
        data_array_aug[data_array_aug == 0] = np.nan

        conv_data_array_jun = conv_data_array_jun/data_array_jun
        conv_data_array_jul = conv_data_array_jul/data_array_jul
        conv_data_array_aug = conv_data_array_aug/data_array_aug

        time_data_array_jun = time_data_array_jun/data_array_jun
        time_data_array_jul = time_data_array_jul/data_array_jul
        time_data_array_aug = time_data_array_aug/data_array_aug

        offset_data_array_jun = time_data_array_jun/offset_counts_array_jun
        offset_data_array_jul = time_data_array_jul/offset_counts_array_jul
        offset_data_array_aug = time_data_array_aug/offset_counts_array_aug

        np.save('multiyear_dust_frequency_array', data_array)
        np.save('multiyear_dust_conv_array', conv_data_array)
        np.save('multiyear_dust_time_array', time_data_array)
        np.save('multiyear_dust_offset_array', offset_data_array)

        np.save('multiyear_dust_frequency_array', data_array_jun)
        np.save('multiyear_dust_conv_array', conv_data_array_jun)
        np.save('multiyear_dust_time_array', time_data_array_jun)
        np.save('multiyear_dust_offset_array', offset_data_array_jun)

        np.save('multiyear_dust_frequency_array', data_array_jul)
        np.save('multiyear_dust_conv_array', conv_data_array_jul)
        np.save('multiyear_dust_time_array', time_data_array_jul)
        np.save('multiyear_dust_offset_array', offset_data_array_jul)

        np.save('multiyear_dust_frequency_array', data_array_aug)
        np.save('multiyear_dust_conv_array', conv_data_array_aug)
        np.save('multiyear_dust_time_array', time_data_array_aug)
        np.save('multiyear_dust_offset_array', offset_data_array_aug)

        plt.close()

        extent = (np.min(lons), np.max(lons), np.min(lats), np.max(lats))
        m = Basemap(projection='cyl', llcrnrlon=extent[0], urcrnrlon=extent[1],
                    llcrnrlat=extent[2], urcrnrlat=extent[3],
                    resolution='i')

        m.drawcoastlines(linewidth=0.5)
        m.drawcountries(linewidth=0.5)
        parallels = np.arange(10., 40, 2.)
        # labels = [left,right,top,bottom]
        m.drawparallels(parallels, labels=[False, True, True, False],
                        linewidth=0.5)
        meridians = np.arange(-20., 17., 4.)
        m.drawmeridians(meridians, labels=[True, False, False, True],
                        linewidth=0.5)

        min = 0
        max = 500

        levels = MaxNLocator(nbins=15).tick_values(min, max)

        # discrete_cmap = utilities.cmap_discretize(cm.RdYlBu_r, 10)
        cmap = cm.get_cmap('Greys')
        norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)

        # m.pcolormesh(lons, lats, correlation_array, extent=extent,
        #             origin='lower',
        #         interpolation='none',
        #         cmap=cmap, vmin=min, vmax=max+1)

        conv_data_array = np.ma.masked_where(np.isnan(conv_data_array),
                                               conv_data_array)

        m.pcolormesh(lons, lats, conv_data_array,
                     cmap=cmap, norm=norm, vmin=min, vmax=max)

        cbar = plt.colorbar(orientation='horizontal', fraction=0.056, pad=0.06)
        cbar.ax.set_xlabel('Distance to convection (km)')
        plt.tight_layout()
        plt.savefig(title_1, bbox_inches='tight')

        plt.close()

        extent = (np.min(lons), np.max(lons), np.min(lats), np.max(lats))
        m = Basemap(projection='cyl', llcrnrlon=extent[0], urcrnrlon=extent[1],
                    llcrnrlat=extent[2], urcrnrlat=extent[3],
                    resolution='i')

        m.drawcoastlines(linewidth=0.5)
        m.drawcountries(linewidth=0.5)
        parallels = np.arange(10., 40, 2.)
        # labels = [left,right,top,bottom]
        m.drawparallels(parallels, labels=[False, True, True, False],
                        linewidth=0.5)
        meridians = np.arange(-20., 17., 4.)
        m.drawmeridians(meridians, labels=[True, False, False, True],
                        linewidth=0.5)

        min = 0
        max = 40000

        levels = MaxNLocator(nbins=15).tick_values(min, max)

        # discrete_cmap = utilities.cmap_discretize(cm.RdYlBu_r, 10)
        cmap = cm.get_cmap('Greys')
        norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)

        # m.pcolormesh(lons, lats, correlation_array, extent=extent,
        #             origin='lower',
        #         interpolation='none',
        #         cmap=cmap, vmin=min, vmax=max+1)

        time_data_array = np.ma.masked_where(np.isnan(time_data_array),
                                             time_data_array)

        m.pcolormesh(lons, lats, time_data_array,
                     cmap=cmap, norm=norm, vmin=min, vmax=max)

        cbar = plt.colorbar(orientation='horizontal', fraction=0.056, pad=0.06)
        cbar.ax.set_xlabel('Time to 0900UTC (seconds)')
        plt.tight_layout()
        plt.savefig(title_2, bbox_inches='tight')

        plt.close()

        extent = (np.min(lons), np.max(lons), np.min(lats), np.max(lats))
        m = Basemap(projection='cyl', llcrnrlon=extent[0], urcrnrlon=extent[1],
                    llcrnrlat=extent[2], urcrnrlat=extent[3],
                    resolution='i')

        m.drawcoastlines(linewidth=0.5)
        m.drawcountries(linewidth=0.5)
        parallels = np.arange(10., 40, 2.)
        # labels = [left,right,top,bottom]
        m.drawparallels(parallels, labels=[False, True, True, False],
                        linewidth=0.5)
        meridians = np.arange(-20., 17., 4.)
        m.drawmeridians(meridians, labels=[True, False, False, True],
                        linewidth=0.5)

        min = 0
        max = 90

        levels = MaxNLocator(nbins=15).tick_values(min, max)

        # discrete_cmap = utilities.cmap_discretize(cm.RdYlBu_r, 10)
        cmap = cm.get_cmap('Greys')
        norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)

        # m.pcolormesh(lons, lats, correlation_array, extent=extent,
        #             origin='lower',
        #         interpolation='none',
        #         cmap=cmap, vmin=min, vmax=max+1)

        offset_data_array = np.ma.masked_where(np.isnan(offset_data_array),
                                             offset_data_array)

        m.pcolormesh(lons, lats, offset_data_array,
                     cmap=cmap, norm=norm, vmin=min, vmax=max)

        cbar = plt.colorbar(orientation='horizontal', fraction=0.056, pad=0.06)
        cbar.ax.set_xlabel('Offset from major axis (degrees)')
        plt.tight_layout()
        plt.savefig(title_3, bbox_inches='tight')

        plt.close()

        plt.close()

        extent = (np.min(lons), np.max(lons), np.min(lats), np.max(lats))
        m = Basemap(projection='cyl', llcrnrlon=extent[0], urcrnrlon=extent[1],
                    llcrnrlat=extent[2], urcrnrlat=extent[3],
                    resolution='i')

        m.drawcoastlines(linewidth=0.5)
        m.drawcountries(linewidth=0.5)
        parallels = np.arange(10., 40, 2.)
        # labels = [left,right,top,bottom]
        m.drawparallels(parallels, labels=[False, True, True, False],
                        linewidth=0.5)
        meridians = np.arange(-20., 17., 4.)
        m.drawmeridians(meridians, labels=[True, False, False, True],
                        linewidth=0.5)

        min = 0
        max = 500

        levels = MaxNLocator(nbins=15).tick_values(min, max)

        # discrete_cmap = utilities.cmap_discretize(cm.RdYlBu_r, 10)
        cmap = cm.get_cmap('Greys')
        norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)

        # m.pcolormesh(lons, lats, correlation_array, extent=extent,
        #             origin='lower',
        #         interpolation='none',
        #         cmap=cmap, vmin=min, vmax=max+1)

        conv_data_array_jun = np.ma.masked_where(np.isnan(conv_data_array_jun),
                                             conv_data_array_jun)

        m.pcolormesh(lons, lats, conv_data_array_jun,
                     cmap=cmap, norm=norm, vmin=min, vmax=max)

        cbar = plt.colorbar(orientation='horizontal', fraction=0.056, pad=0.06)
        cbar.ax.set_xlabel('Distance to convection (km)')
        plt.tight_layout()
        plt.savefig('june_'+title_1, bbox_inches='tight')

        plt.close()

        extent = (np.min(lons), np.max(lons), np.min(lats), np.max(lats))
        m = Basemap(projection='cyl', llcrnrlon=extent[0], urcrnrlon=extent[1],
                    llcrnrlat=extent[2], urcrnrlat=extent[3],
                    resolution='i')

        m.drawcoastlines(linewidth=0.5)
        m.drawcountries(linewidth=0.5)
        parallels = np.arange(10., 40, 2.)
        # labels = [left,right,top,bottom]
        m.drawparallels(parallels, labels=[False, True, True, False],
                        linewidth=0.5)
        meridians = np.arange(-20., 17., 4.)
        m.drawmeridians(meridians, labels=[True, False, False, True],
                        linewidth=0.5)

        min = 0
        max = 40000

        levels = MaxNLocator(nbins=15).tick_values(min, max)

        # discrete_cmap = utilities.cmap_discretize(cm.RdYlBu_r, 10)
        cmap = cm.get_cmap('Greys')
        norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)

        # m.pcolormesh(lons, lats, correlation_array, extent=extent,
        #             origin='lower',
        #         interpolation='none',
        #         cmap=cmap, vmin=min, vmax=max+1)

        time_data_array_jun = np.ma.masked_where(np.isnan(time_data_array_jun),
                                             time_data_array_jun)

        m.pcolormesh(lons, lats, time_data_array_jun,
                     cmap=cmap, norm=norm, vmin=min, vmax=max)

        cbar = plt.colorbar(orientation='horizontal', fraction=0.056, pad=0.06)
        cbar.ax.set_xlabel('Time to 0900UTC (seconds)')
        plt.tight_layout()
        plt.savefig('june_'+title_2, bbox_inches='tight')

        plt.close()

        extent = (np.min(lons), np.max(lons), np.min(lats), np.max(lats))
        m = Basemap(projection='cyl', llcrnrlon=extent[0], urcrnrlon=extent[1],
                    llcrnrlat=extent[2], urcrnrlat=extent[3],
                    resolution='i')

        m.drawcoastlines(linewidth=0.5)
        m.drawcountries(linewidth=0.5)
        parallels = np.arange(10., 40, 2.)
        # labels = [left,right,top,bottom]
        m.drawparallels(parallels, labels=[False, True, True, False],
                        linewidth=0.5)
        meridians = np.arange(-20., 17., 4.)
        m.drawmeridians(meridians, labels=[True, False, False, True],
                        linewidth=0.5)

        min = 0
        max = 90

        levels = MaxNLocator(nbins=15).tick_values(min, max)

        # discrete_cmap = utilities.cmap_discretize(cm.RdYlBu_r, 10)
        cmap = cm.get_cmap('Greys')
        norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)

        # m.pcolormesh(lons, lats, correlation_array, extent=extent,
        #             origin='lower',
        #         interpolation='none',
        #         cmap=cmap, vmin=min, vmax=max+1)

        offset_data_array_jun = np.ma.masked_where(np.isnan(offset_data_array_jun),
                                             offset_data_array_jun)

        m.pcolormesh(lons, lats, offset_data_array_jun,
                     cmap=cmap, norm=norm, vmin=min, vmax=max)

        cbar = plt.colorbar(orientation='horizontal', fraction=0.056, pad=0.06)
        cbar.ax.set_xlabel('Offset from major axis (degrees)')
        plt.tight_layout()
        plt.savefig('june_'+title_3, bbox_inches='tight')

        plt.close()

        plt.close()

        extent = (np.min(lons), np.max(lons), np.min(lats), np.max(lats))
        m = Basemap(projection='cyl', llcrnrlon=extent[0], urcrnrlon=extent[1],
                    llcrnrlat=extent[2], urcrnrlat=extent[3],
                    resolution='i')

        m.drawcoastlines(linewidth=0.5)
        m.drawcountries(linewidth=0.5)
        parallels = np.arange(10., 40, 2.)
        # labels = [left,right,top,bottom]
        m.drawparallels(parallels, labels=[False, True, True, False],
                        linewidth=0.5)
        meridians = np.arange(-20., 17., 4.)
        m.drawmeridians(meridians, labels=[True, False, False, True],
                        linewidth=0.5)

        min = 0
        max = 500

        levels = MaxNLocator(nbins=15).tick_values(min, max)

        # discrete_cmap = utilities.cmap_discretize(cm.RdYlBu_r, 10)
        cmap = cm.get_cmap('Greys')
        norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)

        # m.pcolormesh(lons, lats, correlation_array, extent=extent,
        #             origin='lower',
        #         interpolation='none',
        #         cmap=cmap, vmin=min, vmax=max+1)

        conv_data_array_jul = np.ma.masked_where(np.isnan(conv_data_array_jul),
                                                 conv_data_array_jul)

        m.pcolormesh(lons, lats, conv_data_array_jul,
                     cmap=cmap, norm=norm, vmin=min, vmax=max)

        cbar = plt.colorbar(orientation='horizontal', fraction=0.056, pad=0.06)
        cbar.ax.set_xlabel('Distance to convection (km)')
        plt.tight_layout()
        plt.savefig('july_' + title_1, bbox_inches='tight')

        plt.close()

        extent = (np.min(lons), np.max(lons), np.min(lats), np.max(lats))
        m = Basemap(projection='cyl', llcrnrlon=extent[0], urcrnrlon=extent[1],
                    llcrnrlat=extent[2], urcrnrlat=extent[3],
                    resolution='i')

        m.drawcoastlines(linewidth=0.5)
        m.drawcountries(linewidth=0.5)
        parallels = np.arange(10., 40, 2.)
        # labels = [left,right,top,bottom]
        m.drawparallels(parallels, labels=[False, True, True, False],
                        linewidth=0.5)
        meridians = np.arange(-20., 17., 4.)
        m.drawmeridians(meridians, labels=[True, False, False, True],
                        linewidth=0.5)

        min = 0
        max = 40000

        levels = MaxNLocator(nbins=15).tick_values(min, max)

        # discrete_cmap = utilities.cmap_discretize(cm.RdYlBu_r, 10)
        cmap = cm.get_cmap('Greys')
        norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)

        # m.pcolormesh(lons, lats, correlation_array, extent=extent,
        #             origin='lower',
        #         interpolation='none',
        #         cmap=cmap, vmin=min, vmax=max+1)

        time_data_array_jul = np.ma.masked_where(np.isnan(time_data_array_jul),
                                                 time_data_array_jul)

        m.pcolormesh(lons, lats, time_data_array_jul,
                     cmap=cmap, norm=norm, vmin=min, vmax=max)

        cbar = plt.colorbar(orientation='horizontal', fraction=0.056, pad=0.06)
        cbar.ax.set_xlabel('Time to 0900UTC (seconds)')
        plt.tight_layout()
        plt.savefig('july_' + title_2, bbox_inches='tight')

        plt.close()

        extent = (np.min(lons), np.max(lons), np.min(lats), np.max(lats))
        m = Basemap(projection='cyl', llcrnrlon=extent[0], urcrnrlon=extent[1],
                    llcrnrlat=extent[2], urcrnrlat=extent[3],
                    resolution='i')

        m.drawcoastlines(linewidth=0.5)
        m.drawcountries(linewidth=0.5)
        parallels = np.arange(10., 40, 2.)
        # labels = [left,right,top,bottom]
        m.drawparallels(parallels, labels=[False, True, True, False],
                        linewidth=0.5)
        meridians = np.arange(-20., 17., 4.)
        m.drawmeridians(meridians, labels=[True, False, False, True],
                        linewidth=0.5)

        min = 0
        max = 90

        levels = MaxNLocator(nbins=15).tick_values(min, max)

        # discrete_cmap = utilities.cmap_discretize(cm.RdYlBu_r, 10)
        cmap = cm.get_cmap('Greys')
        norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)

        # m.pcolormesh(lons, lats, correlation_array, extent=extent,
        #             origin='lower',
        #         interpolation='none',
        #         cmap=cmap, vmin=min, vmax=max+1)

        offset_data_array_jul = np.ma.masked_where(
            np.isnan(offset_data_array_jul),
            offset_data_array_jul)

        m.pcolormesh(lons, lats, offset_data_array_jul,
                     cmap=cmap, norm=norm, vmin=min, vmax=max)

        cbar = plt.colorbar(orientation='horizontal', fraction=0.056, pad=0.06)
        cbar.ax.set_xlabel('Offset from major axis (degrees)')
        plt.tight_layout()
        plt.savefig('july_' + title_3, bbox_inches='tight')

        plt.close()

        plt.close()

        extent = (np.min(lons), np.max(lons), np.min(lats), np.max(lats))
        m = Basemap(projection='cyl', llcrnrlon=extent[0], urcrnrlon=extent[1],
                    llcrnrlat=extent[2], urcrnrlat=extent[3],
                    resolution='i')

        m.drawcoastlines(linewidth=0.5)
        m.drawcountries(linewidth=0.5)
        parallels = np.arange(10., 40, 2.)
        # labels = [left,right,top,bottom]
        m.drawparallels(parallels, labels=[False, True, True, False],
                        linewidth=0.5)
        meridians = np.arange(-20., 17., 4.)
        m.drawmeridians(meridians, labels=[True, False, False, True],
                        linewidth=0.5)

        min = 0
        max = 500

        levels = MaxNLocator(nbins=15).tick_values(min, max)

        # discrete_cmap = utilities.cmap_discretize(cm.RdYlBu_r, 10)
        cmap = cm.get_cmap('Greys')
        norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)

        # m.pcolormesh(lons, lats, correlation_array, extent=extent,
        #             origin='lower',
        #         interpolation='none',
        #         cmap=cmap, vmin=min, vmax=max+1)

        conv_data_array_aug = np.ma.masked_where(np.isnan(conv_data_array_aug),
                                                 conv_data_array_aug)

        m.pcolormesh(lons, lats, conv_data_array_aug,
                     cmap=cmap, norm=norm, vmin=min, vmax=max)

        cbar = plt.colorbar(orientation='horizontal', fraction=0.056, pad=0.06)
        cbar.ax.set_xlabel('Distance to convection (km)')
        plt.tight_layout()
        plt.savefig('august_' + title_1, bbox_inches='tight')

        plt.close()

        extent = (np.min(lons), np.max(lons), np.min(lats), np.max(lats))
        m = Basemap(projection='cyl', llcrnrlon=extent[0], urcrnrlon=extent[1],
                    llcrnrlat=extent[2], urcrnrlat=extent[3],
                    resolution='i')

        m.drawcoastlines(linewidth=0.5)
        m.drawcountries(linewidth=0.5)
        parallels = np.arange(10., 40, 2.)
        # labels = [left,right,top,bottom]
        m.drawparallels(parallels, labels=[False, True, True, False],
                        linewidth=0.5)
        meridians = np.arange(-20., 17., 4.)
        m.drawmeridians(meridians, labels=[True, False, False, True],
                        linewidth=0.5)

        min = 0
        max = 40000

        levels = MaxNLocator(nbins=15).tick_values(min, max)

        # discrete_cmap = utilities.cmap_discretize(cm.RdYlBu_r, 10)
        cmap = cm.get_cmap('Greys')
        norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)

        # m.pcolormesh(lons, lats, correlation_array, extent=extent,
        #             origin='lower',
        #         interpolation='none',
        #         cmap=cmap, vmin=min, vmax=max+1)

        time_data_array_aug = np.ma.masked_where(np.isnan(time_data_array_aug),
                                                 time_data_array_aug)

        m.pcolormesh(lons, lats, time_data_array_aug,
                     cmap=cmap, norm=norm, vmin=min, vmax=max)

        cbar = plt.colorbar(orientation='horizontal', fraction=0.056, pad=0.06)
        cbar.ax.set_xlabel('Time to 0900UTC (seconds)')
        plt.tight_layout()
        plt.savefig('august_' + title_2, bbox_inches='tight')

        plt.close()

        extent = (np.min(lons), np.max(lons), np.min(lats), np.max(lats))
        m = Basemap(projection='cyl', llcrnrlon=extent[0], urcrnrlon=extent[1],
                    llcrnrlat=extent[2], urcrnrlat=extent[3],
                    resolution='i')

        m.drawcoastlines(linewidth=0.5)
        m.drawcountries(linewidth=0.5)
        parallels = np.arange(10., 40, 2.)
        # labels = [left,right,top,bottom]
        m.drawparallels(parallels, labels=[False, True, True, False],
                        linewidth=0.5)
        meridians = np.arange(-20., 17., 4.)
        m.drawmeridians(meridians, labels=[True, False, False, True],
                        linewidth=0.5)

        min = 0
        max = 90

        levels = MaxNLocator(nbins=15).tick_values(min, max)

        # discrete_cmap = utilities.cmap_discretize(cm.RdYlBu_r, 10)
        cmap = cm.get_cmap('Greys')
        norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)

        # m.pcolormesh(lons, lats, correlation_array, extent=extent,
        #             origin='lower',
        #         interpolation='none',
        #         cmap=cmap, vmin=min, vmax=max+1)

        offset_data_array_aug = np.ma.masked_where(
            np.isnan(offset_data_array_aug),
            offset_data_array_aug)

        m.pcolormesh(lons, lats, offset_data_array_aug,
                     cmap=cmap, norm=norm, vmin=min, vmax=max)

        cbar = plt.colorbar(orientation='horizontal', fraction=0.056, pad=0.06)
        cbar.ax.set_xlabel('Offset from major axis (degrees)')
        plt.tight_layout()
        plt.savefig('august_' + title_3, bbox_inches='tight')

        plt.close()

def plot_multiyear_cpo_sweeps(pre_calculated=False):
    """
    Plots a map of the sum of CPO sweeps across all years of data availability
    :return:
    """

    if pre_calculated:
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
        lons, lats = np.meshgrid(cloud_lons, cloud_lats)

        # Get lats and lons
        sdf_test = Dataset(
            '/soge-home/data_not_backed_up/satellite/meteosat/seviri'
            '/15-min/0.03x0.03/sdf/nc/JUNE2010/SDF_v2/SDF_v2.'
            '201006031500.nc')

        lons, lats = np.meshgrid(sdf_test.variables['longitude'][:],
                                 sdf_test.variables['latitude'][:])

        plt.close()

        extent = (np.min(lons), np.max(lons), np.min(lats), np.max(lats))
        m = Basemap(projection='cyl', llcrnrlon=extent[0], urcrnrlon=extent[1],
                    llcrnrlat=extent[2], urcrnrlat=extent[3],
                    resolution='i')

        m.drawcoastlines(linewidth=0.5)
        m.drawcountries(linewidth=0.5)
        parallels = np.arange(10., 40, 2.)
        # labels = [left,right,top,bottom]
        m.drawparallels(parallels, labels=[False, True, True, False],
                        linewidth=0.5)
        meridians = np.arange(-20., 17., 4.)
        m.drawmeridians(meridians, labels=[True, False, False, True],
                        linewidth=0.5)

        min = 0
        max = 50

        levels = MaxNLocator(nbins=15).tick_values(min, max)

        # discrete_cmap = utilities.cmap_discretize(cm.RdYlBu_r, 10)
        cmap = cm.get_cmap('Blues')
        norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)

        data_array = np.load('cpo_array_v2.npy')

        data_array = np.ma.masked_where(np.isnan(data_array),
                                        data_array)

        m.pcolormesh(lons, lats, data_array,
                     cmap=cmap, norm=norm, vmin=min, vmax=max)

        cbar = plt.colorbar(orientation='horizontal', fraction=0.056, pad=0.06)
        cbar.ax.set_xlabel('Counts of CPO sweep occurrence')
        plt.tight_layout()

        # Read in cloud frequency
        cloud_frequency = np.load('cloudcover_array2010.npy')

        # Read in cloud lons, lats
        # Get cloud lats and lons
        cloud_test = Dataset(
            '/soge-home/data/satellite/meteosat/seviri/15-min/'
            '0.03x0.03/cloudmask'
            '/nc/'
            +
            'JUNE2010_CLOUDS/eumetsat.cloud.'
            + '201006031500.nc')
        cloud_lons = cloud_test.variables['lon'][:]
        cloud_lats = cloud_test.variables['lat'][:]
        orig_lons, orig_lats = np.meshgrid(cloud_lons, cloud_lats)

        #cloud_frequency[cloud_frequency < 24000] = np.nan
        #cloud_frequency[cloud_frequency > 25000] = np.nan

        clouds_frequency_regridded = pinkdust.regrid_data(orig_lons, orig_lats,
                                                          lons,
                                                          lats,
                                                          cloud_frequency,
                                                          mesh=True)

        clouds_frequency_regridded[lats > 22] = np.nan
        clouds_frequency_regridded[lons < -15.5] = np.nan

        m.contour(lons, lats, clouds_frequency_regridded, levels=[24000],
                  linestyles='dashed', colors='firebrick')

        plt.savefig('CPO_multiyear_sweep_frequency_2004_2017_clouds.png',
                    bbox_inches='tight')

        plt.close()

    else:

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
        lons, lats = np.meshgrid(cloud_lons, cloud_lats)

        # Get lats and lons
        sdf_test = Dataset(
            '/soge-home/data_not_backed_up/satellite/meteosat/seviri'
            '/15-min/0.03x0.03/sdf/nc/JUNE2010/SDF_v2/SDF_v2.'
            '201006031500.nc')

        lons, lats = np.meshgrid(sdf_test.variables['longitude'][:],
                                 sdf_test.variables['latitude'][:])

        # Get cpo archives
        ca2004 = shelve.open('/soge-home/projects/seviri_dust/plumetracker/'
                             'cpo_archive_v4_2004')
        ca2005 = shelve.open('/soge-home/projects/seviri_dust/plumetracker/'
                             'cpo_archive_v4_2005')
        ca2006 = shelve.open('/soge-home/projects/seviri_dust/plumetracker/'
                             'cpo_archive_v4_2006')
        ca2007 = shelve.open('/soge-home/projects/seviri_dust/plumetracker/'
                             'cpo_archive_v4_2007')
        ca2008 = shelve.open('/soge-home/projects/seviri_dust/plumetracker/'
                             'cpo_archive_v4_2008')
        ca2009 = shelve.open('/soge-home/projects/seviri_dust/plumetracker/'
                             'cpo_archive_v4_2009')
        ca2010 = shelve.open('/soge-home/projects/seviri_dust/plumetracker/'
                             'cpo_archive_v4_2010')
        ca2011 = shelve.open('/soge-home/projects/seviri_dust/plumetracker/'
                             'cpo_archive_v4_2011')
        ca2012 = shelve.open('/soge-home/projects/seviri_dust/plumetracker/'
                             'cpo_archive_v4_2012')
        ca2013 = shelve.open('/soge-home/projects/seviri_dust/plumetracker/'
                             'cpo_archive_v4_2013')
        ca2014 = shelve.open('/soge-home/projects/seviri_dust/plumetracker/'
                             'cpo_archive_v4_2014')
        ca2015 = shelve.open('/soge-home/projects/seviri_dust/plumetracker/'
                             'cpo_archive_v4_2015')
        ca2016 = shelve.open('/soge-home/projects/seviri_dust/plumetracker/'
                             'cpo_archive_v4_2016')
        ca2017 = shelve.open('/soge-home/projects/seviri_dust/plumetracker/'
                             'cpo_archive_v4_2017')
        cpo_archives = [ca2004, ca2005, ca2006, ca2007, ca2008, ca2009, ca2010,
                        ca2011, ca2012, ca2013, ca2014, ca2015, ca2016,
                        ca2017]

        archive_idx = 1

        data_array = np.zeros((lons.shape))

        count = 0

        for cpo_archive in cpo_archives:
            print 'CPO archive', archive_idx
            archive_idx += 1
            archive_size = len(cpo_archive)
            cpo_idx = 1
            used_percentages = []
            for i in cpo_archive:
                if int((float(cpo_idx) / float(
                        archive_size)) * 100) % 10 == 0 and int((float(
                    cpo_idx) / float(archive_size)) * 100) not in \
                        used_percentages:
                    print str(int((float(cpo_idx) / float(
                        archive_size)) * 100)) + "%"
                    # This percentage has been printed already
                    used_percentages.append(int((float(cpo_idx) / float(
                        archive_size)) * 100))
                cpo_idx += 1

                cpo = cpo_archive[i]
                cpo_bools = cpo.track_plume_bool
                cpo_lons = lons[cpo_bools[0].toarray()]
                # Remove artifact stretching across image
                count += 1

                #cpo_sweep = np.zeros((lons.shape)) # COMMENT THIS LINE FOR CPO
                # SWEEP V1
                min = 0
                max = 0
                for j in cpo_bools:
                    cpo_bool = j.toarray()
                    if archive_idx == 2:
                        if np.min(lons[cpo_bool]) < -17:
                    #        print cpo.plume_id
                            continue
                    #cpo_sweep[cpo_bool] = 1 # COMMENT THIS LINE FOR CPO SWEEP V1
                    if cpo_bool.shape == data_array.shape:
                        data_array += cpo_bool
                    else:
                        print 'Misshapen CPO found'
                        count -= 1
                #data_array += cpo_sweep # COMMENT THIS LINE FOR CPO SWEEP V1
            print 'Found', count, 'CPOs'
            print np.unique(data_array)

        plt.close()

        extent = (np.min(lons), np.max(lons), np.min(lats), np.max(lats))
        m = Basemap(projection='cyl', llcrnrlon=extent[0], urcrnrlon=extent[1],
                    llcrnrlat=extent[2], urcrnrlat=extent[3],
                    resolution='i')

        m.drawcoastlines(linewidth=0.5)
        m.drawcountries(linewidth=0.5)
        parallels = np.arange(10., 40, 2.)
        # labels = [left,right,top,bottom]
        m.drawparallels(parallels, labels=[False, True, True, False],
                        linewidth=0.5)
        meridians = np.arange(-20., 17., 4.)
        m.drawmeridians(meridians, labels=[True, False, False, True],
                        linewidth=0.5)

        min = 0
        max = 50

        levels = MaxNLocator(nbins=15).tick_values(min, max)

        # discrete_cmap = utilities.cmap_discretize(cm.RdYlBu_r, 10)
        cmap = cm.get_cmap('Blues')
        norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)

        # m.pcolormesh(lons, lats, correlation_array, extent=extent,
        #             origin='lower',
        #         interpolation='none',
        #         cmap=cmap, vmin=min, vmax=max+1)

        np.save('cpo_array_v2', data_array)

        data_array = np.ma.masked_where(np.isnan(data_array),
                                             data_array)

        m.pcolormesh(lons, lats, data_array,
                     cmap=cmap, norm=norm, vmin=min, vmax=max)

        cbar = plt.colorbar(orientation='horizontal', fraction=0.056, pad=0.06)
        cbar.ax.set_xlabel('Counts of CPO sweep occurrence')
        plt.tight_layout()

        # Read in cloud frequency
        cloud_frequency = np.load('cloudcover_array2010.npy')

        # Read in cloud lons, lats
        # Get cloud lats and lons
        cloud_test = Dataset(
            '/soge-home/data/satellite/meteosat/seviri/15-min/'
            '0.03x0.03/cloudmask'
            '/nc/'
            +
            'JUNE2010_CLOUDS/eumetsat.cloud.'
            + '201006031500.nc')
        cloud_lons = cloud_test.variables['lon'][:]
        cloud_lats = cloud_test.variables['lat'][:]
        orig_lons, orig_lats = np.meshgrid(cloud_lons, cloud_lats)

        clouds_frequency_regridded = pinkdust.regrid_data(orig_lons, orig_lats,
                                                    lons,
                                                    lats,
                                                    cloud_frequency, mesh=True)

        clouds_frequency_regridded[lats > 22] = np.nan
        clouds_frequency_regridded[lons < -15.7] = np.nan

        m.contour(lons, lats, clouds_frequency_regridded, levels=[24000],
                  linestyles='dashed', colors='firebrick')

        plt.savefig('CPO_multiyear_sweep_frequency_2004_2017_clouds_v2.png',
                    bbox_inches='tight')

        plt.close()

        for ca in cpo_archives:
            ca.close()

def plot_multiyear_cloud_cover(pre_calculated = False):
    """
    Goes through every cloud mask file and plots the count of cloud coverage in each region
    Note that daily cloudmask files are not binary -> values greater than 1
    need to be selected, while for Ian's cloudmask files it's binary.
    :return:
    """

    # Get cloud lats and lons
    cloud_test = Dataset(
        '/soge-home/data/satellite/meteosat/seviri/15-min/'
        '0.03x0.03/cloudmask'
        '/nc/'
        +
        'JUNE2010_CLOUDS/eumetsat.cloud.'
        + '201006031500.nc')
    cloud_lons = cloud_test.variables['lon'][:]
    cloud_lats = cloud_test.variables['lat'][:]
    target_lons, target_lats = np.meshgrid(cloud_lons, cloud_lats)

    target_area = utils.load_area(
        '/soge-home/projects/seviri_dust/areas.def',
        'NorthAfrica')
    root2_lons, root2_lats = target_area.get_lonlats()

    if pre_calculated:
        data_array = np.load('cloudcover_array.npy')
        extent = (
        np.min(target_lons), np.max(target_lons), np.min(target_lats),
        np.max(target_lats))
        m = Basemap(projection='cyl', llcrnrlon=extent[0], urcrnrlon=extent[1],
                    llcrnrlat=extent[2], urcrnrlat=extent[3],
                    resolution='i')

        m.drawcoastlines(linewidth=0.5)
        m.drawcountries(linewidth=0.5)
        parallels = np.arange(10., 40, 2.)
        # labels = [left,right,top,bottom]
        m.drawparallels(parallels, labels=[False, True, True, False],
                        linewidth=0.5)
        meridians = np.arange(-20., 17., 4.)
        m.drawmeridians(meridians, labels=[True, False, False, True],
                        linewidth=0.5)

        min = 0
        max = 70000

        levels = MaxNLocator(nbins=15).tick_values(min, max)

        # discrete_cmap = utilities.cmap_discretize(cm.RdYlBu_r, 10)
        cmap = cm.get_cmap('RdYlBu_r')
        norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)

        # m.pcolormesh(lons, lats, correlation_array, extent=extent,
        #             origin='lower',
        #         interpolation='none',
        #         cmap=cmap, vmin=min, vmax=max+1)

        data_array_masked = np.ma.masked_where(np.isnan(data_array),
                                        data_array)

        m.pcolormesh(target_lons, target_lats, data_array_masked,
                     cmap=cmap, norm=norm, vmin=min, vmax=max)

        cbar = plt.colorbar(orientation='horizontal', fraction=0.056, pad=0.06)
        cbar.ax.set_xlabel('Counts of cloud occurrence')
        plt.tight_layout()
        plt.savefig('CPO_multiyear_cloud_frequency_2004_2012.png',
                    bbox_inches='tight')

        plt.close()

    years = np.arange(2004, 2011)

    datetimes = []

    for i in years:

        year_lower = i
        year_upper = i
        month_lower = 6
        month_upper = 8
        day_lower = 1
        day_upper = 31
        hour_lower = 0
        hour_upper = 23
        minute_lower = 0
        minute_upper = 45

        time_params = np.array([year_lower, year_upper, month_lower,
                                month_upper, day_lower, day_upper,
                                hour_lower, hour_upper, minute_lower,
                                minute_upper])

        year_datetimes = utilities.get_datetime_objects(time_params)
        datetimes = np.append(datetimes, year_datetimes)

    root1 = '/soge-home/data/satellite/meteosat/seviri/15-min/0.03x0.03/cloudmask/nc/'
    root2 = '/soge-home/projects/seviri_dust/raw_seviri_data/cloudmask_nc/'

    year_cm_roots = [root1, root1, root1, root1, root1, root1, root1]

    data_array = np.zeros(target_lons.shape)

    for date_i in np.arange(0, len(datetimes)):
        date = datetimes[date_i]
        print date
        year_bool = np.asarray(years == date.year)
        root = np.asarray(year_cm_roots)[year_bool][0]
        if root == root2:
            # We're using TCH's daily cloudmask files. These are values from
            #  0-2 (I think), and need selecting and regridding.
            try:
                clouddata = Dataset(root+date.strftime("%B%Y")+'/cloudmask_' + datetimes[
                        date_i].strftime("%Y%m%d") + '.nc')
            except:
                print 'No cloud data for', date
                continue
            cloudmask_times = num2date(clouddata.variables['time'][:],
                                       clouddata.variables['time'].units)
            cloudmask_times = np.asarray([dt.datetime(j.year, j.month,
                                                            j.day, j.hour,
                                                            j.minute) for j
                                          in cloudmask_times])
            cloudmask_bool = cloudmask_times == datetimes[date_i]
            clouds_data = clouddata.variables['cloud_mask'][:]
            clouds_now = clouds_data[cloudmask_bool]

            clouds_now_regridded = pinkdust.regrid_data(root2_lons, root2_lats,
                                                        target_lons,
                                                        target_lats,
                                                        clouds_now, mesh=True)
            data_array[clouds_now_regridded > 1] += 1

        else:
            # We're using Ian's original cloud mask files. These are just a
            # binary one or zero and won't be regridded.
            try:
                clouddata = Dataset(root+
                    datetimes[date_i].strftime("%B").upper(
                    ) + str(datetimes[date_i].year) + '_CLOUDS/eumetsat.cloud.'
                    + datetimes[
                        date_i].strftime("%Y%m%d%H%M") + '.nc')
            except:
                print 'No cloud data for', date
                continue
            clouds_now = clouddata.variables['cmask'][:][0]
            data_array[clouds_now == 1] += 1
            if date.day == 31 and date.hour == 23 and date.minute == 45:
                extent = (
                np.min(target_lons), np.max(target_lons), np.min(target_lats),
                np.max(target_lats))
                m = Basemap(projection='cyl', llcrnrlon=extent[0],
                            urcrnrlon=extent[1],
                            llcrnrlat=extent[2], urcrnrlat=extent[3],
                            resolution='i')

                m.drawcoastlines(linewidth=0.5)
                m.drawcountries(linewidth=0.5)
                parallels = np.arange(10., 40, 2.)
                # labels = [left,right,top,bottom]
                m.drawparallels(parallels, labels=[False, True, True, False],
                                linewidth=0.5)
                meridians = np.arange(-20., 17., 4.)
                m.drawmeridians(meridians, labels=[True, False, False, True],
                                linewidth=0.5)

                #levels = MaxNLocator(nbins=15).tick_values(min, max)

                # discrete_cmap = utilities.cmap_discretize(cm.RdYlBu_r, 10)
                cmap = cm.get_cmap('RdYlBu_r')
                #norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)

                # m.pcolormesh(lons, lats, correlation_array, extent=extent,
                #             origin='lower',
                #         interpolation='none',
                #         cmap=cmap, vmin=min, vmax=max+1)

                data_array_masked = np.ma.masked_where(np.isnan(data_array),
                                                data_array)

                m.pcolormesh(target_lons, target_lats, data_array_masked,
                             cmap=cmap)

                cbar = plt.colorbar(orientation='horizontal', fraction=0.056,
                                    pad=0.06)
                cbar.ax.set_xlabel('Counts of cloud occurrence')
                plt.tight_layout()
                plt.savefig('CPO_multiyear_cloud_frequency_'+str(
                    date.year)+'.png',
                            bbox_inches='tight')

                plt.close()

                np.save('cloudcover_array'+str(date.year), data_array)

    np.save('cloudcover_array', data_array)

    extent = (np.min(target_lons), np.max(target_lons), np.min(target_lats),
              np.max(target_lats))
    m = Basemap(projection='cyl', llcrnrlon=extent[0], urcrnrlon=extent[1],
                llcrnrlat=extent[2], urcrnrlat=extent[3],
                resolution='i')

    m.drawcoastlines(linewidth=0.5)
    m.drawcountries(linewidth=0.5)
    parallels = np.arange(10., 40, 2.)
    # labels = [left,right,top,bottom]
    m.drawparallels(parallels, labels=[False, True, True, False],
                    linewidth=0.5)
    meridians = np.arange(-20., 17., 4.)
    m.drawmeridians(meridians, labels=[True, False, False, True],
                    linewidth=0.5)

    min = 0
    max = 70000

    levels = MaxNLocator(nbins=15).tick_values(min, max)

    # discrete_cmap = utilities.cmap_discretize(cm.RdYlBu_r, 10)
    cmap = cm.get_cmap('RdYlBu_r')
    norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)

    # m.pcolormesh(lons, lats, correlation_array, extent=extent,
    #             origin='lower',
    #         interpolation='none',
    #         cmap=cmap, vmin=min, vmax=max+1)

    data_array = np.ma.masked_where(np.isnan(data_array),
                                    data_array)

    m.pcolormesh(target_lons, target_lats, data_array,
                 cmap=cmap, norm=norm, vmin=min, vmax=max)

    cbar = plt.colorbar(orientation='horizontal', fraction=0.056, pad=0.06)
    cbar.ax.set_xlabel('Counts of cloud occurrence')
    plt.tight_layout()
    plt.savefig('CPO_multiyear_cloud_frequency_2004_2012.png',
                bbox_inches='tight')

    plt.close()


def plot_multiyear_cpo_dust_assoc():
    """
    Plots a map the percentage of dust emission events which are associated
    with a sweep from a CPO
    :return:
    """

    # A CPO sweep at any given point in time is the sum of the CPO booleans
    # (no greater than 1) - lasting for one hour after the end of the CPO
    # Dust is associated with the CPO if it is contained within 50km of the
    # CPO sweep

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
    lons, lats = np.meshgrid(cloud_lons, cloud_lats)

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

    archive_idx = 1

    data_array = np.zeros((lons.shape))

    for cpo_archive in cpo_archives:
        print 'CPO archive', archive_idx
        archive_idx += 1
        archive_size = len(cpo_archive)
        cpo_idx = 1
        used_percentages = []
        for i in cpo_archive:
            if int((float(cpo_idx) / float(
                    archive_size)) * 100) % 10 == 0 and int((float(
                cpo_idx) / float(archive_size)) * 100) not in \
                    used_percentages:
                print str(int((float(cpo_idx) / float(
                    archive_size)) * 100)) + "%"
                # This percentage has been printed already
                used_percentages.append(int((float(cpo_idx) / float(
                    archive_size)) * 100))
            cpo_idx += 1

            cpo = cpo_archive[i]
            cpo_bools = cpo.track_plume_bool

            for j in cpo_bools:
                cpo_bool = j.toarray()
                data_array += cpo_bool

    plt.close()

    extent = (np.min(lons), np.max(lons), np.min(lats), np.max(lats))
    m = Basemap(projection='cyl', llcrnrlon=extent[0], urcrnrlon=extent[1],
                llcrnrlat=extent[2], urcrnrlat=extent[3],
                resolution='i')

    m.drawcoastlines(linewidth=0.5)
    m.drawcountries(linewidth=0.5)
    parallels = np.arange(10., 40, 2.)
    # labels = [left,right,top,bottom]
    m.drawparallels(parallels, labels=[False, True, True, False],
                    linewidth=0.5)
    meridians = np.arange(-20., 17., 4.)
    m.drawmeridians(meridians, labels=[True, False, False, True],
                    linewidth=0.5)

    min = 0
    max = 20

    levels = MaxNLocator(nbins=15).tick_values(min, max)

    # discrete_cmap = utilities.cmap_discretize(cm.RdYlBu_r, 10)
    cmap = cm.get_cmap('Greys')
    norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)

    # m.pcolormesh(lons, lats, correlation_array, extent=extent,
    #             origin='lower',
    #         interpolation='none',
    #         cmap=cmap, vmin=min, vmax=max+1)

    data_array = np.ma.masked_where(np.isnan(data_array),
                                         data_array)

    m.pcolormesh(lons, lats, data_array,
                 cmap=cmap, norm=norm, vmin=min, vmax=max)

    cbar = plt.colorbar(orientation='horizontal', fraction=0.056, pad=0.06)
    cbar.ax.set_xlabel('Counts of CPO sweep occurrence')
    plt.tight_layout()
    plt.savefig('CPO_multiyear_sweep_frequency.png', bbox_inches='tight')

    plt.close()

def plot_multiyear_cpo_sweep_directions():
    """
    Plots a map of the sum of CPO sweeps across all years of data availability
    :return:
    """

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
    lons, lats = np.meshgrid(cloud_lons, cloud_lats)

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

    archive_idx = 1

    data_array = np.zeros((lons.shape))
    direction_data_array = np.zeros((lons.shape))
    counts_array = np.zeros((lons.shape))

    for cpo_archive in cpo_archives:
        print 'CPO archive', archive_idx
        archive_idx += 1
        archive_size = len(cpo_archive)
        cpo_idx = 1
        used_percentages = []
        for i in cpo_archive:
            if int((float(cpo_idx) / float(
                    archive_size)) * 100) % 10 == 0 and int((float(
                cpo_idx) / float(archive_size)) * 100) not in \
                    used_percentages:
                print str(int((float(cpo_idx) / float(
                    archive_size)) * 100)) + "%"
                # This percentage has been printed already
                used_percentages.append(int((float(cpo_idx) / float(
                    archive_size)) * 100))
            cpo_idx += 1

            cpo = cpo_archive[i]
            cpo_bools = cpo.track_plume_bool
            cpo_directions = cpo.track_speed_centroid

            for j in np.arange(0, len(cpo_bools)):
                cpo_bool = cpo_bools[j].toarray()
                data_array += cpo_bool
                if np.isfinite(cpo_directions[j]):
                    direction_data_array[cpo_bool == 1] += cpo_directions[j]
                    counts_array += cpo_bool

    plt.close()

    direction_data_array = direction_data_array/counts_array

    extent = (np.min(lons), np.max(lons), np.min(lats), np.max(lats))
    m = Basemap(projection='cyl', llcrnrlon=extent[0], urcrnrlon=extent[1],
                llcrnrlat=extent[2], urcrnrlat=extent[3],
                resolution='i')

    m.drawcoastlines(linewidth=0.5)
    m.drawcountries(linewidth=0.5)
    parallels = np.arange(10., 40, 2.)
    # labels = [left,right,top,bottom]
    m.drawparallels(parallels, labels=[False, True, True, False],
                    linewidth=0.5)
    meridians = np.arange(-20., 17., 4.)
    m.drawmeridians(meridians, labels=[True, False, False, True],
                    linewidth=0.5)

    min = 0
    max = 50

    levels = MaxNLocator(nbins=15).tick_values(min, max)

    # discrete_cmap = utilities.cmap_discretize(cm.RdYlBu_r, 10)
    cmap = cm.get_cmap('coolwarm')
    norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)

    # m.pcolormesh(lons, lats, correlation_array, extent=extent,
    #             origin='lower',
    #         interpolation='none',
    #         cmap=cmap, vmin=min, vmax=max+1)

    direction_data_array = np.ma.masked_where(np.isnan(direction_data_array),
                                         direction_data_array)

    m.pcolormesh(lons, lats, direction_data_array,
                 cmap=cmap, norm=norm, vmin=min, vmax=max)

    cbar = plt.colorbar(orientation='horizontal', fraction=0.056, pad=0.06)
    cbar.ax.set_xlabel('CPO sweep direction')
    plt.tight_layout()
    plt.savefig('CPO_multiyear_sweep_direction.png', bbox_inches='tight')

    plt.close()

def plot_multiyear_cpo_sweep_speed():
    """
    Plots a map of the sum of CPO sweeps across all years of data availability
    :return:
    """

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
    lons, lats = np.meshgrid(cloud_lons, cloud_lats)

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
    ca2011 = shelve.open('/soge-home/projects/seviri_dust/plumetracker/'
                         'cpo_archive_2011')
    ca2012 = shelve.open('/soge-home/projects/seviri_dust/plumetracker/'
                         'cpo_archive_2012')

    cpo_archives = [ca2004, ca2005, ca2006, ca2007, ca2008, ca2009, ca2010,
                    ca2011, ca2012]

    archive_idx = 1

    data_array = np.zeros((lons.shape))
    direction_data_array = np.zeros((lons.shape))
    counts_array = np.zeros((lons.shape))

    for cpo_archive in cpo_archives:
        print 'CPO archive', archive_idx
        archive_idx += 1
        archive_size = len(cpo_archive)
        cpo_idx = 1
        used_percentages = []
        for i in cpo_archive:
            if int((float(cpo_idx) / float(
                    archive_size)) * 100) % 10 == 0 and int((float(
                cpo_idx) / float(archive_size)) * 100) not in \
                    used_percentages:
                print str(int((float(cpo_idx) / float(
                    archive_size)) * 100)) + "%"
                # This percentage has been printed already
                used_percentages.append(int((float(cpo_idx) / float(
                    archive_size)) * 100))
            cpo_idx += 1

            cpo = cpo_archive[i]
            cpo_bools = cpo.track_plume_bool
            cpo_directions = cpo.track_speed_centroid

            for j in np.arange(0, len(cpo_bools)):
                cpo_bool = cpo_bools[j].toarray()
                data_array += cpo_bool
                if np.isfinite(cpo_directions[j]):
                    direction_data_array[cpo_bool == 1] += cpo_directions[j]
                    counts_array += cpo_bool

    plt.close()

    direction_data_array = direction_data_array/counts_array

    extent = (np.min(lons), np.max(lons), np.min(lats), np.max(lats))
    m = Basemap(projection='cyl', llcrnrlon=extent[0], urcrnrlon=extent[1],
                llcrnrlat=extent[2], urcrnrlat=extent[3],
                resolution='i')

    m.drawcoastlines(linewidth=0.5)
    m.drawcountries(linewidth=0.5)
    parallels = np.arange(10., 40, 2.)
    # labels = [left,right,top,bottom]
    m.drawparallels(parallels, labels=[False, True, True, False],
                    linewidth=0.5)
    meridians = np.arange(-20., 17., 4.)
    m.drawmeridians(meridians, labels=[True, False, False, True],
                    linewidth=0.5)

    min = 0
    max = 50

    levels = MaxNLocator(nbins=15).tick_values(min, max)

    # discrete_cmap = utilities.cmap_discretize(cm.RdYlBu_r, 10)
    cmap = cm.get_cmap('coolwarm')
    norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)

    # m.pcolormesh(lons, lats, correlation_array, extent=extent,
    #             origin='lower',
    #         interpolation='none',
    #         cmap=cmap, vmin=min, vmax=max+1)

    direction_data_array = np.ma.masked_where(np.isnan(direction_data_array),
                                         direction_data_array)

    m.pcolormesh(lons, lats, direction_data_array,
                 cmap=cmap, norm=norm, vmin=min, vmax=max)

    cbar = plt.colorbar(orientation='horizontal', fraction=0.056, pad=0.06)
    cbar.ax.set_xlabel('CPO sweep speed')
    plt.tight_layout()
    plt.savefig('CPO_multiyear_sweep_speed.png', bbox_inches='tight')

    plt.close()

def plot_mechanism_pixel_piechart(pre_calculated=False):
    """
    Plots a pie chart showing the proportion of pre-merge pixels covered at
    peak size associated with LLJs and CPOs
    :return:
    """

    if pre_calculated:
        counts = [2.41735814e+14, 7.81714255e+14]
        fig, ax1 = plt.subplots()
        labels = ['LLJ', 'Other (CPO)']
        colors = ['cornflowerblue', 'lightcoral']
        ax1.pie(counts, labels=labels, autopct='%1.1f%%', startangle=90,
                colors=colors)
        ax1.axis('equal')
        plt.savefig('Emission_area_mechanism_partition.png')
        plt.close()

    else:
        pa2004 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                             'plume_archive_flicker_v3_prob_v4_2004')
        pa2005 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                             'plume_archive_flicker_v3_prob_v4_2005')
        pa2006 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                             'plume_archive_flicker_v3_prob_v4_2006')
        pa2007 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                             'plume_archive_flicker_v3_prob_v4_2007')
        pa2008 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                             'plume_archive_flicker_v3_prob_v4_2008')
        pa2009 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                             'plume_archive_flicker_v3_prob_v4_2009')
        pa2010 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                             'plume_archive_flicker_v3_prob_v4_2010')
        pa2011 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                             'plume_archive_flicker_v3_prob_v3_2011')
        pa2012 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                             'plume_archive_flicker_v3_prob_v3_2012')
        pa2013 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                             'plume_archive_flicker_v3_prob_v4_2013')
        pa2014 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                             'plume_archive_flicker_v3_prob_v4_2014')
        pa2015 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                             'plume_archive_flicker_v3_prob_v4_2015')
        pa2016 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                             'plume_archive_flicker_v3_prob_v4_2016')
        pa2017 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                             'plume_archive_flicker_v3_prob_v4_2017')
        plume_archives = [pa2004, pa2005, pa2006, pa2007, pa2008, pa2009,
                          pa2010, pa2011, pa2012, pa2013, pa2014, pa2015,
                          pa2016, pa2017]

        archive_idx = 1

        # Data array
        counts = np.zeros((2))

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

                plume = plume_archive[i]

                if plume.merged:
                    merge_idx = plume.merge_idx
                    plume_areas = plume.track_area[0:merge_idx]
                    if len(plume_areas) > 0:
                        max_plume_area = np.max(plume_areas)
                    else:
                        continue
                else:
                    plume_areas = plume.track_area[:]
                    if len(plume_areas) > 0:
                        max_plume_area = np.max(plume_areas)
                    else:
                        continue

                llj_prob = plume.LLJ_prob
                if llj_prob == None:
                    continue
                if np.isnan(llj_prob):
                    continue
                elif np.isnan(max_plume_area):
                    continue
                else:
                    if llj_prob > 0.5:
                        counts[0] += max_plume_area
                    else:
                        counts[1] += max_plume_area

            print counts

        fig, ax1 = plt.subplots()
        labels = ['LLJ', 'Other (CPO)']
        print counts
        ax1.pie(counts, labels=labels, autopct='%1.1f%%',
                shadow=True, startangle=90)
        ax1.axis('equal')
        plt.savefig('Emission_area_emission_partition.png')
        plt.close()

def plot_mechanism_count_piechart(pre_calculated=False):
    """
    Plots a pie chart showing the proportion of pre-merge pixels covered at
    peak size associated with LLJs and CPOs
    :return:
    """

    if pre_calculated:
        counts = [12185, 31571]
        fig, ax1 = plt.subplots()
        labels = ['LLJ', 'Other (CPO)']
        colors = ['cornflowerblue', 'lightcoral']
        ax1.pie(counts, labels=labels, autopct='%1.1f%%', startangle=90,
                colors=colors)
        ax1.axis('equal')
        plt.savefig('Emission_count_partition.png')
        plt.close()

    else:
        pa2004 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                             'plume_archive_flicker_v3_prob_v4_2004')
        pa2005 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                             'plume_archive_flicker_v3_prob_v4_2005')
        pa2006 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                             'plume_archive_flicker_v3_prob_v4_2006')
        pa2007 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                             'plume_archive_flicker_v3_prob_v4_2007')
        pa2008 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                             'plume_archive_flicker_v3_prob_v4_2008')
        pa2009 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                             'plume_archive_flicker_v3_prob_v4_2009')
        pa2010 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                             'plume_archive_flicker_v3_prob_v4_2010')
        pa2011 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                             'plume_archive_flicker_v3_prob_v3_2011')
        pa2012 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                             'plume_archive_flicker_v3_prob_v3_2012')
        pa2013 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                             'plume_archive_flicker_v3_prob_v4_2013')
        pa2014 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                             'plume_archive_flicker_v3_prob_v4_2014')
        pa2015 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                             'plume_archive_flicker_v3_prob_v4_2015')
        pa2016 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                             'plume_archive_flicker_v3_prob_v4_2016')
        pa2017 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                             'plume_archive_flicker_v3_prob_v4_2017')
        plume_archives = [pa2004, pa2005, pa2006, pa2007, pa2008, pa2009,
                          pa2010, pa2011, pa2012, pa2013, pa2014, pa2015,
                          pa2016, pa2017]

        archive_idx = 1

        # Data array
        counts = np.zeros((2))

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

                plume = plume_archive[i]

                llj_prob = plume.LLJ_prob
                if llj_prob == None:
                    continue
                if np.isnan(llj_prob):
                    continue

                else:
                    if llj_prob > 0.5:
                        counts[0] += 1
                    else:
                        counts[1] += 1

            print counts

        fig, ax1 = plt.subplots()
        labels = ['LLJ', 'Other (CPO)']
        print counts
        ax1.pie(counts, labels=labels, autopct='%1.1f%%',
                shadow=True, startangle=90)
        ax1.axis('equal')
        plt.savefig('Emission_count_partition.png')
        plt.close()

def plot_cpo_diurnal_cycle():
    """
    Plots a diurnal cycle of emission frequency of cold pool outflows from
    the samples available in the CPO archive
    :return:
    """

    # Get cpo archives
    ca2004 = shelve.open('/soge-home/projects/seviri_dust/plumetracker/'
                         'cpo_archive_v3_2004')
    ca2005 = shelve.open('/soge-home/projects/seviri_dust/plumetracker/'
                         'cpo_archive_v3_2005')
    ca2006 = shelve.open('/soge-home/projects/seviri_dust/plumetracker/'
                         'cpo_archive_v3_2006')
    ca2007 = shelve.open('/soge-home/projects/seviri_dust/plumetracker/'
                         'cpo_archive_v3_2007')
    ca2008 = shelve.open('/soge-home/projects/seviri_dust/plumetracker/'
                         'cpo_archive_v3_2008')
    ca2009 = shelve.open('/soge-home/projects/seviri_dust/plumetracker/'
                         'cpo_archive_v3_2009')
    ca2010 = shelve.open('/soge-home/projects/seviri_dust/plumetracker/'
                         'cpo_archive_v3_2010')
    ca2011 = shelve.open('/soge-home/projects/seviri_dust/plumetracker/'
                         'cpo_archive_2011')
    ca2012 = shelve.open('/soge-home/projects/seviri_dust/plumetracker/'
                         'cpo_archive_2012')
    ca2013 = shelve.open('/soge-home/projects/seviri_dust/plumetracker/'
                         'cpo_archive_v3_2013')
    ca2014 = shelve.open('/soge-home/projects/seviri_dust/plumetracker/'
                         'cpo_archive_v3_2014')
    ca2015 = shelve.open('/soge-home/projects/seviri_dust/plumetracker/'
                         'cpo_archive_v3_2015')
    ca2016 = shelve.open('/soge-home/projects/seviri_dust/plumetracker/'
                         'cpo_archive_v3_2016')
    ca2017 = shelve.open('/soge-home/projects/seviri_dust/plumetracker/'
                         'cpo_archive_v3_2017')
    cpo_archives = [ca2004, ca2005, ca2006, ca2007, ca2008, ca2009, ca2010,
                    ca2011, ca2012, ca2013, ca2014, ca2015, ca2016, ca2017]

    # Get an array of datetimes at 15 minute resolution through the day
    # Get a bool by matching the emission time from each CPO to these datetimes
    # Then use that bool to increment a counter for each time of day
    # Plot that and bob's your uncle

    year_lower = 2013
    year_upper = 2013
    month_lower = 6
    month_upper = 6
    day_lower = 1
    day_upper = 1
    hour_lower = 0
    hour_upper = 23
    minute_lower = 0
    minute_upper = 45

    # Generate a set of datetimes for one day, and loop through just one day
    time_params_oneday = np.array([year_lower, year_lower, month_lower,
                                   month_lower, day_lower, day_lower,
                                   hour_lower, hour_upper,
                                   minute_lower, minute_upper])

    # It doesn't particularly matter what the year, month, day is, as we'll
    # just match the hour and minute


    oneday_datetimes = utilities.get_datetime_objects(time_params_oneday)

    diurnal_cycle = np.zeros((oneday_datetimes.shape))

    archive_idx = 1

    for plume_archive in cpo_archives:
        print 'CPO archive', archive_idx
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

            plume = plume_archive[i]

            emission_time = plume.dates_observed[0]
            print emission_time
            diurnal_cycle_bool = np.asarray([j.hour == emission_time.hour
                                             and j.minute == emission_time.minute
                                             for j in oneday_datetimes])
            diurnal_cycle[diurnal_cycle_bool] += 1

    plt.plot(oneday_datetimes, diurnal_cycle)
    plt.xlabel('Time')
    plt.ylabel('CPO emissions')
    plt.savefig('CPO_emission_diurnal_cycle.png')

def plot_interannual_dust_area():
    """
    Plots a bar chart with the total area of dust plumes emitted, quantified as
    the peak area which each plume reaches before it merges with another or
    dies
    :return:
    """

def plot_tidihelt_dustsourcing(pre_calculated=False, res=20):
    """
    For each emission event in the Tidihelt region, plots the boolean of the
    event and its centroid, as well as the gridded dust source map under a
    slightly transparent bools map (constructed up to that point)
    :return:
    """

    sdf_test = Dataset(
        '/soge-home/data_not_backed_up/satellite/meteosat/seviri'
        '/15-min/0.03x0.03/sdf/nc/JUNE2010/SDF_v2/SDF_v2.'
        '201006031500.nc')

    ianlons, ianlats = np.meshgrid(sdf_test.variables['longitude'][:],
                             sdf_test.variables['latitude'][:])

    tidihelt_emission_dates = np.load('tidihelt_emission_dates.npy')
    for i in np.arange(0, len(tidihelt_emission_dates)):
        emission_date = tidihelt_emission_dates[i]
    midas_filename = '17667.60630.SYNOP.midas.gl.nc'
    insalah = Dataset(midas_filename)

    # Get lats and lons for the extremely generous Tidihelt region and get a
    # linearly spaced 10 by 10 regular grid
    linlats = np.linspace(24, 29, res)
    linlons = np.linspace(-2, 4, res)

    lons, lats = np.meshgrid(linlons, linlats)
    data_array = np.zeros((lats.shape))
    counts_array = np.zeros((lats.shape))
    bool_data_array = np.zeros(ianlons.shape)

    insalah_wspd = insalah.variables['WIND_SPEED'][:]
    insalah_times = num2date(insalah.variables['TIME'][:],
                             insalah.variables['TIME'].units)
    insalah_wdir = insalah.variables['WIND_DIRECTION'][:]


    # Lines to constrain to northeasterly
    neerly_bool = np.asarray([j >= 0 and j <= 90 for j in insalah_wdir])
    neerly_bool = neerly_bool == 1
    insalah_wspd = insalah_wspd[neerly_bool]
    insalah_times = insalah_times[neerly_bool]

    archive_count = 1

    bool_period = np.asarray(
        [(j.year >= 2004 and j.year <=
                            2012
                            and j.month in [6, 7, 8]) for j in
         insalah_times])
    insalah_wspd = insalah_wspd[bool_period]
    insalah_times = insalah_times[bool_period]
    emission_wspds = []
    no_emission_wspds = insalah_wspd
    emission_times = []

    if pre_calculated:

        data_array = np.load('tidihelt_totals.npy')
        counts_array = np.load('tidihelt_counts.npy')

        data_array_mean = data_array / counts_array
        #data_array_mean[counts_array < 5] = np.nan

        plt.close()

        extent = (np.min(lons), np.max(lons), np.min(lats), np.max(lats))
        m = Basemap(projection='cyl', llcrnrlon=extent[0], urcrnrlon=extent[1],
                    llcrnrlat=extent[2], urcrnrlat=extent[3],
                    resolution='i')

        m.drawcoastlines(linewidth=0.5)
        m.drawcountries(linewidth=0.5)
        parallels = np.arange(25, 28, 0.5)
        # labels = [left,right,top,bottom]
        m.drawparallels(parallels, labels=[False, True, True, False],
                        linewidth=0.5)
        meridians = np.arange(-1, 3, 0.5)
        m.drawmeridians(meridians, labels=[True, False, False, True],
                        linewidth=0.5)

        # m.contourf(lons, lats, data_array)
        discrete_cmap = utilities.cmap_discretize(cm.YlOrBr, 10)
        m.imshow(counts_array, extent=extent, origin='lower',
                 interpolation='none',
                 cmap=discrete_cmap)

        plt.colorbar(orientation='horizontal', fraction=0.056, pad=0.06)

        m.scatter(insalah.variables['LON'][:], insalah.variables['LAT'][:])

        m.contour(lons, lats, orog_regridded)

        plt.tight_layout()
        plt.savefig('Tidihelt_midas_counts.png', bbox_inches='tight')

        plt.close()

    else:

        pa2004 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                             'plume_archive_flicker_v3_prob_v4_2004')
        pa2005 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                             'plume_archive_flicker_v3_prob_v4_2005')
        pa2006 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                             'plume_archive_flicker_v3_prob_v4_2006')
        pa2007 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                             'plume_archive_flicker_v3_prob_v4_2007')
        pa2008 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                             'plume_archive_flicker_v3_prob_v4_2008')
        pa2009 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                             'plume_archive_flicker_v3_prob_v4_2009')
        pa2010 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                             'plume_archive_flicker_v3_prob_v4_2010')
        pa2011 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                             'plume_archive_flicker_v3_prob_v3_2011')
        pa2012 = shelve.open('/ouce-home/projects/seviri_dust/plumetracker/'
                             'plume_archive_flicker_v3_prob_v3_2012')
        plume_archives = [pa2004, pa2005, pa2006, pa2007, pa2008, pa2009,
                          pa2010, pa2011, pa2012]

        fig, ax = plt.subplots()
        extent = (np.min(lons), np.max(lons), np.min(lats), np.max(lats))
        m = Basemap(projection='cyl', llcrnrlon=extent[0], urcrnrlon=extent[1],
                    llcrnrlat=extent[2], urcrnrlat=extent[3],
                    resolution='i')

        m.drawcoastlines(linewidth=0.5)
        m.drawcountries(linewidth=0.5)

        parallels = np.arange(10., 40, 2.)
        # labels = [left,right,top,bottom]
        m.drawparallels(parallels, labels=[False, True, True, False],
                        linewidth=0.5)
        meridians = np.arange(-20., 17., 4.)
        m.drawmeridians(meridians, labels=[True, False, False, True],
                        linewidth=0.5)

        for plume_archive in plume_archives:
            print 'Plume archive', archive_count
            archive_count += 1
            plume_idx = 0
            archive_size = len(plume_archive)
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
                plume = plume_archive[i]
                emission_lat = plume.track_centroid_lat[0]
                emission_lon = plume.track_centroid_lon[0]
                emission_bool = plume.track_plume_bool[0].toarray()
                if np.round(emission_lat, 0) in [25, 26, 27, 28
                                                 ] and np.round(
                    emission_lon, 0) in [-1, 0, 1, 2, 3]:
                    emission_time = plume.emission_time
                    # Find the nearest emission time from MIDAS 3-hourly
                    available_hours = [0, 3, 6, 9, 12, 15, 18, 21]
                    possible_times = [plume.emission_time.
                                          replace(hour=available_hours[k],
                                                  minute=0) for
                                      k in
                                                np.arange(0,
                                                          len(available_hours))]
                    time_differences = [abs(emission_time-j).total_seconds() for j
                                        in
                                        possible_times]
                    nearest_time = np.asarray(possible_times)[
                        time_differences==np.min(
                        time_differences)]

                    emission_time = nearest_time[0]
                    time_bool = insalah_times == emission_time
                    if np.all(time_bool == False):
                        continue
                    insalah_wspd_matched = insalah_wspd[time_bool][0]
                    if np.ma.is_masked(insalah_wspd_matched):
                        continue
                    emission_wspds.append(insalah_wspd_matched)
                    emission_times.append(emission_time)
                    diff_array_lat = abs(lats - emission_lat)
                    diff_array_lon = abs(lons - emission_lon)
                    latbool = diff_array_lat == np.min(diff_array_lat)
                    lonbool = diff_array_lon == np.min(diff_array_lon)
                    union_bool = latbool & lonbool
                    data_array[union_bool] += insalah_wspd_matched
                    counts_array[union_bool] += 1
                    bool_data_array[emission_bool] += 1

                    # Here mask the emission bool and centroid and plot
                    # Underneath should be data_array
                    emission_bool = np.ma.masked_where(emission_bool ==
                                                       False, emission_bool)

                    print np.unique(counts_array)
                    print np.unique(bool_data_array)

                    grid = m.pcolormesh(lons, lats, counts_array, vmin=0,
                                     vmax=70,
                                 cmap=cm.YlOrBr)
                    boolarr = m.pcolormesh(ianlons, ianlats, bool_data_array,
                                   vmin=0, vmax=70, alpha=0.5, cmap=cm.YlOrBr)
                    nowbool = m.pcolormesh(ianlons, ianlats, emission_bool,
                                 cmap=cm.YlOrBr_r)
                    nowgrid = m.scatter(emission_lon, emission_lat,
                                        marker='+', s=80,
                              color='blue')

                    plt.savefig('Tidihelt_emission_plume_'+str(
                        archive_count)+"_"+str(i)+".png")

                    grid.remove()
                    boolarr.remove()
                    nowbool.remove()
                    nowgrid.remove()

        plt.close()

        remove_bool = np.asarray([j in emission_times for j in insalah_times])
        no_emission_wspds = np.asarray(no_emission_wspds)
        no_emission_wspds[remove_bool] = -999

        no_emission_wspds = no_emission_wspds[no_emission_wspds!=-999]

        #emission_wspds = emission_wspds[emission_wspds < 100]
        no_emission_wspds = no_emission_wspds[no_emission_wspds < 100]

        np.save('tidihelt_totals', data_array)
        np.save('tidihelt_counts', counts_array)

        data_array_mean = data_array/counts_array

        plt.close()

        extent = (np.min(lons), np.max(lons), np.min(lats), np.max(lats))
        m = Basemap(projection='cyl', llcrnrlon=extent[0], urcrnrlon=extent[1],
                    llcrnrlat=extent[2], urcrnrlat=extent[3],
                    resolution='i')

        m.drawcoastlines(linewidth=0.5)
        m.drawcountries(linewidth=0.5)
        parallels = np.arange(25, 28, 0.5)
        # labels = [left,right,top,bottom]
        m.drawparallels(parallels, labels=[False, True, True, False],
                        linewidth=0.5)
        meridians = np.arange(-1, 3, 0.5)
        m.drawmeridians(meridians, labels=[True, False, False, True],
                        linewidth=0.5)

        # m.contourf(lons, lats, data_array)
        discrete_cmap = utilities.cmap_discretize(cm.YlOrBr, 10)
        m.imshow(data_array_mean, extent=extent, origin='lower',
                 interpolation='none',
                 cmap=discrete_cmap)

        plt.colorbar(orientation='horizontal', fraction=0.056, pad=0.06)
        m.scatter(insalah.variables['LON'][:], insalah.variables['LAT'][:])
        plt.tight_layout()
        plt.savefig('Tidihelt_midas_wspd_time_matched_northeasterly.png',
                    bbox_inches='tight')

        plt.close()

        fig, ax = plt.subplots()
        bins = np.arange(0, 45, 5)
        ax.hist(emission_wspds, bins=bins, color='cornflowerblue',
                label='Dust '
                                                                      'activation '
                                                                      'event',
                alpha=0.5)
        ax.hist(no_emission_wspds, bins=bins, color='lightcoral', label='No '
                                                                      'dust '
                                                                     'activation '
                                                                     'event',
                alpha=0.5)
        plt.legend()
        ax.set_xlabel('0900UTC In Salah MIDAS wind speed')
        ax.set_ylabel('Frequency')
        plt.savefig('midas_comparison_distributions_northeasterly.png')

        plt.close()







