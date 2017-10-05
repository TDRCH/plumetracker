import numpy as np
import datetime
from math import radians, cos, sin, asin, sqrt

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