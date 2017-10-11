import numpy as np
import multiprocessing
from netCDF4 import Dataset

import utilities
import pinkdust

def wrapper(datetime):
    lonlats = Dataset(
        '/ouce-home/data/satellite/meteosat/seviri/15-min/native/'
        'lonlats.NA_MiddleEast.nc')

    # These need to be regridded to regular for consistency with cloud mask
    lons = lonlats.variables['longitude'][:]
    lats = lonlats.variables['latitude'][:]
    lonmask = lons > 360
    latmask = lats > 90
    lons = np.ma.array(lons, mask=lonmask)
    lats = np.ma.array(lats, mask=latmask)

    cloudmask_now = pinkdust.process_cloud_mask_grib(datetime,
                                                     lats, lons)

    np.save('/ouce-home/students/hert4173/cloud_mask_numpy_files/cloudmask'
            ''+datetime
            .strftime(
            "%Y%m%d%H%M%S"), cloudmask_now)

    print datetime

if __name__ == '__main__':
    # Here testing running of multiprocessing
    #pool = multiprocessing.Pool(processes=4)
    #[pool.apply(wrapper, args=(i,)) for i in range(1, 16)]

    Year_lower = 2012
    Year_upper = 2012
    Month_lower = 6
    Month_upper = 6
    Day_lower = 23
    Day_upper = 26
    Hour_lower = 0
    Hour_upper = 23
    Minute_lower = 0
    Minute_upper = 0

    process_num = 15

    # Generate datetime objects corresponding to these time bounds
    time_params = np.array([Year_lower, Year_upper, Month_lower, Month_upper,
                            Day_lower, Day_upper, Hour_lower, Hour_upper,
                            Minute_lower, Minute_upper])

    datetimes = utilities.get_datetime_objects(time_params)

    for j in np.arange(0, len(datetimes), process_num):
        processes = [multiprocessing.Process(target=wrapper, args=(i,))
                     for i in datetimes[j:j+process_num]]
        for p in processes:
            p.start()
        for p in processes:
            p.join()