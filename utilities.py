import numpy as np
import datetime
import math
from math import radians, cos, sin, asin, sqrt
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.pyplot as plt

def get_datetime_objects(time_params):
    """
    Generates an array of datetime objects at 15 minute intervals within the
    input time bounds
    :param time_params: an array of integers corresponding to time bounds
    :return datetimes: array of datetime objects
    """

    # Get the lower and upper bound of the datetime objects
    datetime_lower = datetime.datetime(time_params[0], time_params[2],
                                       time_params[4], time_params[6],
                                       time_params[8])
    datetime_upper = datetime.datetime(time_params[1], time_params[3],
                                       time_params[5], time_params[7],
                                       time_params[9])

    # Get the difference between the two dates in minutes using a timedelta
    # object
    td = datetime_upper-datetime_lower
    mins_difference = divmod(td.days*86400+td.seconds, 60)[0]

    # Get an array of minute values with 15 minute intervals
    minutes_range = np.arange(0, mins_difference+15, 15)

    # Get an array of datetime objects for each 15 minute interval
    datetimes = np.array([datetime_lower + datetime.timedelta(minutes=i) for
                          i in minutes_range])
    return datetimes

def get_daily_datetime_objects(time_params):
    """
    Generates an array of datetime objects at daily intervals within the input
    time bounds
    :param time_params: an array of integers corresponding to time bounds
    :return datetimes: array of datetime objects
    """

    # Get the lower and upper bound of the datetime objects
    datetime_lower = datetime.datetime(time_params[0], time_params[2],
                                       time_params[4], time_params[6],
                                       time_params[8])
    datetime_upper = datetime.datetime(time_params[1], time_params[3],
                                       time_params[5], time_params[7],
                                       time_params[9])

    # Get the difference between the two dates in minutes using a timedelta
    # object
    td = datetime_upper - datetime_lower
    days_difference = td.days

    # Get an array of minute values with 15 minute intervals
    days_range = np.arange(0, td.days + 1)

    # Get an array of datetime objects for each 15 minute interval
    datetimes = np.array([datetime_lower + datetime.timedelta(days=i)
                          for i in days_range])

    return datetimes

def haversine(lon1, lat1, lon2, lat2):
    """
    Calculate the great circle distance between two points
    on the earth (specified in decimal degrees)
    """
    # convert decimal degrees to radians
    lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2])

    # haversine formula
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = sin(dlat / 2) ** 2 + cos(lat1) * cos(lat2) * sin(dlon / 2) ** 2
    c = 2 * asin(sqrt(a))
    r = 6371  # Radius of earth in kilometers. Use 3956 for miles
    return c * r

def get_neighbours(p, exclude_p=True, shape=None):

    ndim = len(p)

    # generate an (m, ndims) array containing all combinations of 0, 1, 2
    offset_idx = np.indices((3,) * ndim).reshape(ndim, -1).T

    # use these to index into np.array([-1, 0, 1]) to get offsets
    offsets = np.r_[-1, 0, 1].take(offset_idx)

    # optional: exclude offsets of 0, 0, ..., 0 (i.e. p itself)
    if exclude_p:
        offsets = offsets[np.any(offsets, 1)]

    neighbours = p + offsets  # apply offsets to p

    # optional: exclude out-of-bounds indices
    if shape is not None:
        valid = np.all((neighbours < np.array(shape)) & (neighbours >= 0),
                       axis=1)
        neighbours = neighbours[valid]

    return neighbours

def polyarea(x, y):
    return 0.5*np.abs(np.dot(x,np.roll(y,1))-np.dot(y,np.roll(x,1)))

def angle_between(p1, p2):
    ang1 = np.arctan2(*p1[::-1])
    ang2 = np.arctan2(*p2[::-1])
    return np.rad2deg((ang1 - ang2) % (2 * np.pi))

def calculate_initial_compass_bearing(pointA, pointB):
    """
    CREDIT: https://gist.github.com/jeromer/2005586
    Calculates the bearing between two points.
    :Parameters:
      - `pointA: The tuple representing the latitude/longitude for the
        first point. Latitude and longitude must be in decimal degrees
      - `pointB: The tuple representing the latitude/longitude for the
        second point. Latitude and longitude must be in decimal degrees
    :Returns:
      The bearing in degrees
    :Returns Type:
      float
    """

    if (type(pointA) != tuple) or (type(pointB) != tuple):
        raise TypeError("Only tuples are supported as arguments")

    lat1 = math.radians(pointA[0])
    lat2 = math.radians(pointB[0])

    diffLong = math.radians(pointB[1] - pointA[1])

    x = math.sin(diffLong) * math.cos(lat2)
    y = math.cos(lat1) * math.sin(lat2) - (math.sin(lat1)
            * math.cos(lat2) * math.cos(diffLong))

    initial_bearing = math.atan2(x, y)

    # Now we have the initial bearing but math.atan2 return values
    # from -180 to + 180 which is not what we want for a compass bearing
    # The solution is to normalize the initial bearing as shown below
    initial_bearing = math.degrees(initial_bearing)
    compass_bearing = (initial_bearing + 360) % 360

    return compass_bearing

def cmap_discretize(cmap, n_colors):
    """Return a discrete colormap from the continuous colormap cmap.

    Parameters
    ----------
    cmap : str or colormap object
        Colormap to discretize.
    n_colors : int
        Number of discrete colors to divide `cmap` into.

    Returns
    ----------
    discrete_cmap : LinearSegmentedColormap
        Discretized colormap.
    """

    if type(cmap) == basestring:
        cmap = plt.get_cmap(cmap)
    colors_i = np.concatenate((np.linspace(0, 1., n_colors), (0., 0., 0., 0.)))
    colors_rgba = cmap(colors_i)
    indices = np.linspace(0, 1., n_colors + 1)
    cdict = {}
    for ki, key in enumerate(('red', 'green', 'blue')):
        cdict[key] = [(indices[i], colors_rgba[i - 1, ki], colors_rgba[i, ki])
                      for i in range(n_colors + 1)]
    # Return colormap object.
    return LinearSegmentedColormap(cmap.name + "_%d" % n_colors, cdict, 1024)