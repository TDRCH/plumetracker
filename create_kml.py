import matplotlib
matplotlib.use('Agg')
from mpl_toolkits.basemap import Basemap
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from copy import deepcopy
from simplekml import (Kml, OverlayXY, ScreenXY, Units, RotationXY,
                       AltitudeMode, Camera)
from matplotlib import colors as c
from scipy.stats import gaussian_kde


import pinkdust
import utilities

def make_kml(llcrnrlon, llcrnrlat, urcrnrlon, urcrnrlat,
             figs, colorbar=None, **kw):
    """TODO: LatLon bbox, list of figs, optional colorbar figure,
    and several simplekml kw..."""

    kml = Kml()
    altitude = kw.pop('altitude', 2e7)
    roll = kw.pop('roll', 0)
    tilt = kw.pop('tilt', 0)
    altitudemode = kw.pop('altitudemode', AltitudeMode.relativetoground)
    camera = Camera(latitude=np.mean([urcrnrlat, llcrnrlat]),
                    longitude=np.mean([urcrnrlon, llcrnrlon]),
                    altitude=altitude, roll=roll, tilt=tilt,
                    altitudemode=altitudemode)

    kml.document.camera = camera
    draworder = 0
    for fig in figs:  # NOTE: Overlays are limited to the same bbox.
        draworder += 1
        ground = kml.newgroundoverlay(name='GroundOverlay')
        ground.draworder = draworder
        ground.visibility = kw.pop('visibility', 1)
        ground.name = kw.pop('name', 'overlay')
        ground.color = kw.pop('color', '9effffff')
        ground.atomauthor = kw.pop('author', 'ocefpaf')
        ground.latlonbox.rotation = kw.pop('rotation', 0)
        ground.description = kw.pop('description', 'Matplotlib figure')
        ground.gxaltitudemode = kw.pop('gxaltitudemode',
                                       'clampToSeaFloor')
        ground.icon.href = fig
        ground.latlonbox.east = llcrnrlon
        ground.latlonbox.south = llcrnrlat
        ground.latlonbox.north = urcrnrlat
        ground.latlonbox.west = urcrnrlon

    if colorbar:  # Options for colorbar are hard-coded (to avoid a big mess).
        screen = kml.newscreenoverlay(name='ScreenOverlay')
        screen.icon.href = colorbar
        screen.overlayxy = OverlayXY(x=0, y=0,
                                     xunits=Units.fraction,
                                     yunits=Units.fraction)
        screen.screenxy = ScreenXY(x=0.015, y=0.075,
                                   xunits=Units.fraction,
                                   yunits=Units.fraction)
        screen.rotationXY = RotationXY(x=0.5, y=0.5,
                                       xunits=Units.fraction,
                                       yunits=Units.fraction)
        screen.size.x = 0
        screen.size.y = 0
        screen.size.xunits = Units.fraction
        screen.size.yunits = Units.fraction
        screen.visibility = 1

    kmzfile = kw.pop('kmzfile', 'overlay.kmz')
    kml.savekmz(kmzfile)

def gearth_fig(llcrnrlon, llcrnrlat, urcrnrlon, urcrnrlat, pixels=1024):
    """Return a Matplotlib `fig` and `ax` handles for a Google-Earth Image."""
    aspect = np.cos(np.mean([llcrnrlat, urcrnrlat]) * np.pi/180.0)
    xsize = np.ptp([urcrnrlon, llcrnrlon]) * aspect
    ysize = np.ptp([urcrnrlat, llcrnrlat])
    aspect = ysize / xsize

    if aspect > 1.0:
        figsize = (10.0 / aspect, 10.0)
    else:
        figsize = (10.0, 10.0 * aspect)

    if False:
        plt.ioff()  # Make `True` to prevent the KML components from poping-up.
    fig = plt.figure(figsize=figsize,
                     frameon=False,
                     dpi=pixels//10)
    # KML friendly image.  If using basemap try: `fix_aspect=False`.
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(llcrnrlon, urcrnrlon)
    ax.set_ylim(llcrnrlat, urcrnrlat)
    return fig, ax

def plot_dust_source_map():
    sdf_test = Dataset(
        '/soge-home/data_not_backed_up/satellite/meteosat/seviri'
        '/15-min/0.03x0.03/sdf/nc/JUNE2010/SDF_v2/SDF_v2.'
        '201006031500.nc')

    lons, lats = np.meshgrid(sdf_test.variables['longitude'][:],
                             sdf_test.variables['latitude'][:])
    lats = np.linspace(np.min(lats), np.max(lats), 100)
    lons = np.linspace(np.min(lons), np.max(lons), 100)
    lon, lat = np.meshgrid(lons, lats)

    pixels = 1024 * 10

    data_array = np.load('multiyear_dust_frequency_array_noflicker15.npy')
    data_array[np.isnan(data_array)] = 0
    data_array = np.ma.masked_equal(data_array, 0)

    fig, ax = gearth_fig(llcrnrlon=lon.min(),
                         llcrnrlat=lat.min(),
                         urcrnrlon=lon.max(),
                         urcrnrlat=lat.max(),
                         pixels=pixels)
    cs = ax.pcolormesh(lon, lat, data_array, cmap=cm.YlOrBr)
    ax.set_axis_off()
    plt.savefig('dust_sourcemap.png', transparent=False, format='png')
    plt.close()

    fig = plt.figure(figsize=(1.0, 4.0), facecolor=None, frameon=False)
    ax = fig.add_axes([0.0, 0.05, 0.2, 0.9])
    cb = fig.colorbar(cs, cax=ax)
    cb.set_label('Dust activation frequency', rotation=-90, color='k',
                 labelpad=20)
    fig.savefig('dust_legend.png', transparent=False, format='png')
    plt.close()

    fig, ax = gearth_fig(llcrnrlon=lon.min(),
                         llcrnrlat=lat.min(),
                         urcrnrlon=lon.max(),
                         urcrnrlat=lat.max(),
                         pixels=pixels)
    absence_array = deepcopy(data_array)
    absence_array[absence_array == 0] = -999
    absence_array = np.ma.masked_greater(absence_array, 0)
    absence_array[absence_array == -999] = 1
    cMap = c.ListedColormap(['cornflowerblue'])
    cs = ax.pcolormesh(lon, lat, absence_array, cmap=cMap)
    ax.set_axis_off()
    plt.savefig('dust_absence.png', transparent=False, format='png')
    plt.close()

    sdf_test = Dataset(
        '/soge-home/data_not_backed_up/satellite/meteosat/seviri'
        '/15-min/0.03x0.03/sdf/nc/JUNE2010/SDF_v2/SDF_v2.'
        '201006031500.nc')

    lons, lats = np.meshgrid(sdf_test.variables['longitude'][:],
                             sdf_test.variables['latitude'][:])

    fig, ax = gearth_fig(llcrnrlon=lons.min(),
                         llcrnrlat=lats.min(),
                         urcrnrlon=lons.max(),
                         urcrnrlat=lats.max(),
                         pixels=pixels)
    cpo_array = np.load('cpo_array_v2.npy')
    cpo_array = np.ma.masked_equal(cpo_array, 0)
    # Get lats and lons

    cs = ax.pcolormesh(lons, lats, cpo_array, cmap='Blues')
    ax.set_axis_off()
    plt.savefig('cpo_coverage.png', transparent=False, format='png')
    plt.close()

    fig = plt.figure(figsize=(1.0, 4.0), facecolor=None, frameon=False)
    ax = fig.add_axes([0.0, 0.05, 0.2, 0.9])
    cb = fig.colorbar(cs, cax=ax)
    cb.set_label('CPO coverage frequency', rotation=-90, color='k',
                 labelpad=20)
    fig.savefig('cpo_legend.png', transparent=False, format='png')
    plt.close()

    data_array = np.load(
        'multiyear_dust_frequency_array_bool_noflickercheck_2004_2012'
                         '.npy')
    data_array[np.isnan(data_array)] = 0
    data_array = np.ma.masked_equal(data_array, 0)

    fig, ax = gearth_fig(llcrnrlon=lons.min(),
                         llcrnrlat=lats.min(),
                         urcrnrlon=lons.max(),
                         urcrnrlat=lats.max(),
                         pixels=pixels)

    levels = np.arange(0, 160, 10)

    cs = ax.contour(lons, lats, data_array, cmap=cm.YlOrBr, levels=levels)
    ax.set_axis_off()
    plt.savefig('dust_emission_bool.png', transparent=False, format='png')
    plt.close()

    fig = plt.figure(figsize=(1.0, 4.0), facecolor=None, frameon=False)
    ax = fig.add_axes([0.0, 0.05, 0.2, 0.9])
    cb = fig.colorbar(cs, cax=ax)
    cb.set_label('Dust emission coverage frequency', rotation=-90, color='k',
                 labelpad=20)
    fig.savefig('dust_bool_legend.png', transparent=False, format='png')
    plt.close()

    cloud_test = Dataset(
        '/soge-home/data/satellite/meteosat/seviri/15-min/'
        '0.03x0.03/cloudmask'
        '/nc/'
        +
        'JUNE2010_CLOUDS/eumetsat.cloud.'
        + '201006031500.nc')
    cloud_lons = cloud_test.variables['lon'][:]
    cloud_lats = cloud_test.variables['lat'][:]

    fig, ax = gearth_fig(llcrnrlon=cloud_lons.min(),
                         llcrnrlat=cloud_lats.min(),
                         urcrnrlon=cloud_lons.max(),
                         urcrnrlat=cloud_lats.max(),
                         pixels=pixels)

    cloud_array = np.load('cloudcover_array2010.npy')
    cloud_array = np.ma.masked_equal(cloud_array, 0)

    cs = ax.pcolormesh(cloud_lons, cloud_lats, cloud_array, cmap='Reds',
                       vmin=0, vmax=40000)
    ax.set_axis_off()
    plt.savefig('cloud_coverage.png', transparent=False, format='png')
    plt.close()

    fig = plt.figure(figsize=(1.0, 4.0), facecolor=None, frameon=False)
    ax = fig.add_axes([0.0, 0.05, 0.2, 0.9])
    cb = fig.colorbar(cs, cax=ax)
    cb.set_label('Cloud coverage frequency', rotation=-90, color='k',
                 labelpad=20)
    fig.savefig('cloud_legend.png', transparent=False, format='png')
    plt.close()

    centroids_lon = np.load('centroids_lon.npy')
    centroids_lat = np.load('centroids_lat.npy')

    fig, ax = gearth_fig(llcrnrlon=lons.min(),
                         llcrnrlat=lats.min(),
                         urcrnrlon=lons.max(),
                         urcrnrlat=lats.max(),
                         pixels=pixels)

    ax.hist2d(centroids_lon, centroids_lat, (200, 200), cmap=plt.cm.YlOrBr)

    #sc = ax.scatter(centroids_lon, centroids_lat, s=1)

    plt.savefig('emission_hist.png', transparent=False, format='png')
    plt.close()

    make_kml(llcrnrlon=lon.min(), llcrnrlat=lat.min(),
             urcrnrlon=lon.max(), urcrnrlat=lat.max(),
             figs=['dust_sourcemap.png'],
             colorbar='dust_legend.png',
             kmzfile='dust_sourcemap.kmz', name='Dust source activation')

    make_kml(llcrnrlon=lon.min(), llcrnrlat=lat.min(),
             urcrnrlon=lon.max(), urcrnrlat=lat.max(),
             figs=['dust_absence.png'],
             kmzfile='dust_absence.kmz', name='Dust absence')

    make_kml(llcrnrlon=lon.min(), llcrnrlat=lat.min(),
             urcrnrlon=lon.max(), urcrnrlat=lat.max(),
             figs=['cpo_coverage.png'],
             colorbar='cpo_legend.png',
             kmzfile='cpo_coverage.kmz', name='CPO coverage')

    make_kml(llcrnrlon=lon.min(), llcrnrlat=lat.min(),
             urcrnrlon=lon.max(), urcrnrlat=lat.max(),
             figs=['cloud_coverage.png'],
             colorbar='cloud_legend.png',
             kmzfile='cloud_coverage.kmz', name='Cloud coverage')

    make_kml(llcrnrlon=lon.min(), llcrnrlat=lat.min(),
             urcrnrlon=lon.max(), urcrnrlat=lat.max(),
             figs=['dust_emission_bool.png'],
             colorbar='dust_bool_legend.png',
             kmzfile='dust_bool_coverage_noflickercheck.kmz', name='Dust '
                                                                'coverage at t0')

    make_kml(llcrnrlon=lon.min(), llcrnrlat=lat.min(),
             urcrnrlon=lon.max(), urcrnrlat=lat.max(),
             figs=['emission_hist.png'],
             kmzfile='emission_hist.kmz', name='Emission hist')

if __name__ == '__main__':
    plot_dust_source_map()