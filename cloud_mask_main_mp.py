import numpy as np
import multiprocessing
from netCDF4 import Dataset
import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpop.satellites import GeostationaryFactory
from mpop.satout import netcdf4
#from mpop.utils import debug_on

#debug_on()
import datetime
import numpy as np
import iris
import iris.quickplot as qplt
from netCDF4 import Dataset
from netCDF4 import date2num
from netCDF4 import num2date
from pycoast import ContourWriterAGG
from PIL import Image
from pyproj import Proj
from mpl_toolkits.basemap import Basemap
from pyresample.geometry import SwathDefinition
from pyresample.kd_tree import resample_nearest
from pyresample import utils
import pygrib
import glob
#import multiprocessing
#from multiprocessing import Process, Value, Array
import multiprocess
from multiprocess import Process
import tables
from pyresample import image, geometry
import email
import smtplib

import utilities
import pinkdust


def wrapper(date):
    print date
    test1 = datetime.datetime.now()

    # Pull out cloud mask data to be regridded and added to the nc file
    # glob is used for a wildcard
    try:
        clouddata = pygrib.open(glob.glob(
            '/soge-home/projects/seviri_dust/raw_seviri_data/cloudmask_grib/MSG*-'
            'SEVI-MSGCLMK-0100-0100-' + date.strftime(
                "%Y%m%d%H%M%S") + '*')[0])
        grb = clouddata.select()[0]
        cloudmaskarray = grb.values[:, ::-1]
        cloudlats, cloudlons = grb.latlons()

        cloudmaskarray[cloudmaskarray >= 3] = np.nan
        cloudlats[cloudlats > 90] = np.nan
        cloudlons[cloudlons > 180] = np.nan

        target_area = utils.load_area('/soge-home/projects/seviri_dust/areas.def',
                                      'NorthAfrica')
        lons, lats = target_area.get_lonlats()

        # Generate a regular lat/lon grid for the cloud mask
        #regular_lons = np.linspace(np.min(lons), np.max(lons), lons.shape[1])
        #regular_lats = np.linspace(np.min(lats), np.max(lats), lats.shape[0])

        # Regrid cloud mask data to the same grid as the BT data
        cloudmask_regridded = pinkdust.regrid_data(cloudlons,
                                                       cloudlats,
                                                       lons,
                                                       lats,
                                                       cloudmaskarray,
                                                   mesh=True)

        f = tables.open_file('/soge-home/projects/seviri_dust/raw_seviri_data'
                             '/intermediary_files/cloudmask/cloudmask_'+date
                             .strftime(
            '%Y%m%d%H%M%S.hdf'), 'w')
        atom = tables.Atom.from_dtype(cloudmask_regridded.dtype)
        filters = tables.Filters(complib='blosc', complevel=5)
        ds = f.create_carray(f.root, 'data', atom,
                             cloudmask_regridded.shape,
                             filters=filters)
        ds[:] = cloudmask_regridded
        f.close()
    except:
        with open('/soge-home/projects/seviri_dust'
                  '/raw_seviri_data/bt_native/missing_msgnative_dates'
                  '.txt', 'a') as my_file:
            my_file.write('\n')
            my_file.write(date.strftime('%Y%m%d%H%M%S')+' clouds')


if __name__ == '__main__':

    years_list = [2015]
    days_list = [1]

    gmail_user = 'tcharrison4@gmail.com'
    gmail_password = None

    sent_from = gmail_user
    to = ['tcharrison4@gmail.com', 'thomas.catonharrison@keble.ox.ac.uk']
    subject = 'Cloud mask processing: started processing'
    body = 'Dear Thomas,\n\nProcessing of cloud mask data files to HDF has ' \
           'been initiated. Lots of love, \n\nPython'

    email_text = 'Subject: {}\n\n{}'.format(subject, body)

    s = smtplib.SMTP("smtp.gmail.com", 587)
    s.ehlo()
    s.starttls()
    s.login(gmail_user, gmail_password)
    s.sendmail(sent_from, to, email_text)
    s.close()

    for y in np.arange(0, len(years_list)):

        Year_lower = years_list[y]
        Year_upper = years_list[y]
        Month_lower = 6
        Month_upper = 7
        Day_lower = 1
        Day_upper = 31
        Hour_lower = 0
        Hour_upper = 23
        Minute_lower = 0
        Minute_upper = 45

        # Generate datetime objects corresponding to these time bounds
        time_params = np.array([Year_lower, Year_upper, Month_lower, Month_upper,
                                Day_lower, Day_upper, Hour_lower, Hour_upper,
                                Minute_lower, Minute_upper])

        datetimes = utilities.get_datetime_objects(time_params)

        """

        pool = multiprocessing.Pool()
        for i in datetimes:
            pool.apply_async(wrapper, args=(i,))
        pool.close()
        pool.join()

        gmail_user = 'tcharrison4@gmail.com'
        gmail_password = None

        sent_from = gmail_user
        to = ['tcharrison4@gmail.com', 'thomas.catonharrison@keble.ox.ac.uk']
        subject = 'Cloud mask processing: year complete'
        body = 'Dear Thomas,\n\nA year of cloud mask data has been ' \
               'processed.\n\nLots ' \
               'of love, ' \
               '\n\nPython'

        email_text = 'Subject: {}\n\n{}'.format(subject, body)

        s = smtplib.SMTP("smtp.gmail.com", 587)
        s.ehlo()
        s.starttls()
        s.login(gmail_user, gmail_password)
        s.sendmail(sent_from, to, email_text)
        s.close()

        """

        years = np.unique(np.asarray(([j.year for j in datetimes])))

        months = np.unique(np.asarray(([j.month for j in datetimes])))

        days = np.unique(np.asarray(([j.day for j in datetimes])))

        # Pull a netCDF dataset object out so the projection coordinates can
        # be obtained
        target_area = utils.load_area('/soge-home/projects/seviri_dust/areas.def',
                                      'NorthAfrica')
        lons, lats = target_area.get_lonlats()

        # Make a lonlats file
        lonlatsnc = pinkdust.\
            create_lonlats_file('/soge-home/projects/seviri_dust/'
                            'raw_seviri_data/cloudmask_nc/cloudmask_lonlats.nc',
                                lats, lons)
        lonlatsnc.close()

        for m in np.arange(0, len(years)):

            for k in np.arange(0, len(months)):

                for i in np.arange(0, len(days)):
                    day_datetimes_bool = np.asarray([j.day == days[i] and j.month ==
                                                     months[k] and j.year ==
                                                     years[m] for j in datetimes])
                    day_datetimes = datetimes[day_datetimes_bool]

                    if len(day_datetimes) == 0:
                        continue

                    ncfile = pinkdust.create_time_cloudmask_nc_file(
                        '/soge-home/projects/seviri_dust/raw_seviri_data'
                        '/cloudmask_nc/'
                        + day_datetimes[0].strftime(
                            "%B%Y") + '/cloudmask_' + day_datetimes[
                            0].strftime(
                            '%Y%m%d') + '.nc', day_datetimes, lats, lons)

                    for j in np.arange(0, len(day_datetimes)):
                        f = tables.open_file('/soge-home/projects/seviri_dust/raw_seviri_data'
                            '/intermediary_files/cloudmask/cloudmask_' +
                                             datetimes[i].strftime(
                            '%Y%m%d%H%M%S.hdf'))
                        arrobj = f.get_node('/data')
                        cloudmask = arrobj.read()
                        f.close()

                        ncfile = pinkdust.save_cloudmask_to_existing_nc(ncfile,
                                                    cloudmask, day_datetimes[j])

                    print 'Day '+str(days[i])+' done'
                    ncfile.close()
                print 'Month ' + str(months[k]) + ' done'
            print 'Year ' + str(years[m]) + ' done'

            gmail_user = 'tcharrison4@gmail.com'
            gmail_password = None

            sent_from = gmail_user
            to = ['tcharrison4@gmail.com',
                  'thomas.catonharrison@keble.ox.ac.uk']
            subject = 'NetCDF cloud mask processing: complete'
            body = 'Dear Thomas,\n\nA year of cloud mask data has been written ' \
                   'to ' \
                   'netCDF.\n\nLots ' \
                   'of love, ' \
                   '\n\nPython'

            email_text = 'Subject: {}\n\n{}'.format(subject, body)

            s = smtplib.SMTP("smtp.gmail.com", 587)
            s.ehlo()
            s.starttls()
            s.login(gmail_user, gmail_password)
            s.sendmail(sent_from, to, email_text)
            s.close()

