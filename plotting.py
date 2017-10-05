import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import numpy as np
from netCDF4 import Dataset
from mpop.satout import netcdf4

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


