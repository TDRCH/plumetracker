import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import numpy as np
from mpl_toolkits.basemap import Basemap
import datetime
import shelve

import utilities
import plumes

if __name__ == '__main__':
    # NOTES #
    # These parameters should be handled in a spreadsheet
    # Plotting should be handled in a separate function


    SDF_0 = Dataset(
        '/ouce-home/data/satellite/meteosat/seviri/15-min/native/sdf/nc'
        '/JUNE2012/SDF_v2/SDF_v2.201206221800.nc')

    year_lower = 2012
    year_upper = 2012
    month_lower = 6
    month_upper = 6
    day_lower = 1
    day_upper = 30
    hour_lower = 0
    hour_upper = 23
    minute_lower = 0
    minute_upper = 45

    time_params = np.array([year_lower, year_upper, month_lower,
                            month_upper, day_lower, day_upper,
                            hour_lower, hour_upper, minute_lower,
                            minute_upper])

    datetimes = utilities.get_datetime_objects(time_params)
    datestrings = [j.strftime("%Y%m%d%H%M") for j in datetimes]

    lonlats = Dataset(
        '/ouce-home/data/satellite/meteosat/seviri/15-min/native/'
        'lonlats.NA_MiddleEast.nc')
    lons = lonlats.variables['longitude'][:]
    lats = lonlats.variables['latitude'][:]
    lonmask = lons > 360
    latmask = lats > 90
    lons = np.ma.array(lons, mask=lonmask)
    lats = np.ma.array(lats, mask=latmask)
    fig, ax = plt.subplots()
    extent = (-21, 31, 10, 41)
    m = Basemap(projection='cyl', llcrnrlon=extent[0], urcrnrlon=extent[1],
                llcrnrlat=extent[2], urcrnrlat=extent[3],
                resolution='i')

    m.drawcoastlines(linewidth=0.5)
    m.drawcountries(linewidth=0.5)

    sdf_previous = None
    ids_previous = []
    deep_conv_IDs_prev = None
    LLJ_plumes_IDs_prev = []
    k = 0
    available_colours = np.arange(0, 41)
    used_ids = []
    used_colour_IDs = {}
    plume_objects = []

    for date_i in np.arange(0, len(datestrings)):
        print '\n' + datestrings[date_i] + '\n'
        totaltest = datetime.datetime.now()
        sdf = Dataset(
            '/ouce-home/data/satellite/meteosat/seviri/15-min/native/sdf/nc/'
            'JUNE2012/SDF_v2/SDF_v2.' + datestrings[date_i] + '.nc')
        bt = Dataset(
            '/ouce-home/data/satellite/meteosat/seviri/15-min/native/bt/nc/'
            'JUNE2012/H-000-MSG2__-MSG2________-'
            'IR_BrightnessTemperatures___-000005___-' + datestrings[date_i] +
            '-__.nc')

        sdf_now = sdf.variables['bt108'][:]

        sdf_plumes, new_ids, plume_ids = plumes.scan_for_plumes(sdf_now,
                                                              sdf_previous,
                                                                used_ids)

        for i in new_ids:
            used_ids.append(i)

        old_bool = np.asarray([j in ids_previous for j in plume_ids])
        old_ids = plume_ids[old_bool]

        plume_objects = shelve.open('plume_objects')

        # Then, for each new ID, we initialise plume objects
        for i in np.arange(0, len(new_ids)):
            print 'Creating new plume', new_ids[i]
            plume = plumes.Plume(new_ids[i], datetimes[date_i])
            plume.update_position(lats, lons, sdf_plumes, new_ids[i])
            plume.update_duration(datetimes[date_i])
            plume_objects[str(new_ids[i])] = plume

        # For old IDs, we just run an update.
        for i in np.arange(0, len(old_ids)):
            print 'Updating plume', old_ids[i]
            plume = plume_objects[str(old_ids[i])]
            plume.update_position(lats,lons, sdf_plumes, old_ids[i])
            plume.update_duration(datetimes[date_i])
            plume_objects[str(old_ids[i])] = plume

        # Plumes which no longer exist are removed
        if len(ids_previous) == 0:
            removed_ids = []
        else:
            removed_bool = np.asarray([j not in plume_ids for j in
                                       ids_previous])
            removed_ids = ids_previous[removed_bool]


        for i in np.arange(0, len(removed_ids)):
            del plume_objects[str(removed_ids[i])]
            print 'Deleted plume', removed_ids[i]
            # Now ideally here we would create a plume archive, which would
            # have a track for all past plumes and their duration and source
            #  location

        sdf_previous = sdf_plumes
        ids_previous = plume_ids

        print plume_objects['4'].track_lat
        print plume_objects['4'].track_lon
        plume_objects.close()

        """
        deep_conv_assoc, deep_conv_IDs = identify_convection(BT, SDF_plumes,
                                                             new_IDs,
                                                             deep_conv_IDs_prev)
        # print 'Identifying deep convection:', datetime.datetime.now() - test
        # test = datetime.datetime.now()
        LLJ_plumes, LLJ_IDs, new_plumes = classify_as_LLJ(SDF_plumes,
                                                          LLJ_plumes_IDs_prev,
                                                          new_IDs,
                                                          datetimes[k],
                                                          deep_conv_assoc)
        # print 'Identifying LLJs:', datetime.datetime.now() - test

        # test = datetime.datetime.now()

        # These throw an error at present - can't store sparse matrices in an nc file - could do possibly to save them to txt or something
        # Can then read them back in...
        SDF_plumes_s = sparse.csr_matrix(SDF_plumes)
        LLJ_plumes_s = sparse.csr_matrix(LLJ_plumes)
        LLJ_plumes_data = LLJ_plumes_s.data
        LLJ_plumes_indices = LLJ_plumes_s.indices
        LLJ_plumes_indptr = SDF_plumes_s.indptr
        SDF_plumes_data = SDF_plumes_s.data
        SDF_plumes_indices = SDF_plumes_s.indices
        SDF_plumes_indptr = SDF_plumes_s.indptr

        np.save('SDF_sparse_matrices/LLJ_plume_data_' + date, LLJ_plumes_data)
        np.save('SDF_sparse_matrices/LLJ_plume_indices_' + date,
                LLJ_plumes_indices)
        np.save('SDF_sparse_matrices/LLJ_plume_indptr_' + date,
                LLJ_plumes_indptr)
        np.save('SDF_sparse_matrices/SDF_plume_data_' + date, SDF_plumes_data)
        np.save('SDF_sparse_matrices/SDF_plume_indices_' + date,
                SDF_plumes_indices)
        np.save('SDF_sparse_matrices/SDF_plume_indptr_' + date,
                SDF_plumes_indptr)

        # generate_nc_file(SDF_plumes_s, LLJ_plumes_s, 'SDF_LLJ_plumes_'+date+'.nc', lats, lons, datetimes[k])
        # LLJs[k] = LLJ_plumes
        # SDF_plumes_array[k] = SDF_plumes
        # print 'Writing to storage arrays:', datetime.datetime.now() - test
        # np.save('LLJ_plumes_June_2012.npy', LLJs)

        ID_min = 0
        ID_max = 2

        # print plume_IDs

        # test = datetime.datetime.now()
        reavailable_IDs = []
        # Put colours of IDs back into availablity if that ID no longer exists
        for ID in used_colour_IDs:
            if ID not in plume_IDs:
                reavailable_IDs = np.append(reavailable_IDs, ID)
                available_colours = np.append(available_colours,
                                              used_colour_IDs[ID])
        for ID in reavailable_IDs:
            del used_colour_IDs[ID]

        # For IDs without an assigned colour, assign them a random colour from the available ones
        for i in np.arange(0, np.shape(plume_IDs)[0]):
            ID = plume_IDs[i]
            if ID not in used_colour_IDs:
                colour = np.random.choice(available_colours)
                index = np.where(available_colours == colour)[0][0]
                available_colours = np.delete(available_colours, index)
                used_colour_IDs[ID] = colour

        # Now create an equivalent array to the plumes, but where plume IDs are replaced with colour IDs
        coloured_plumes = np.zeros(np.shape(SDF_plumes))
        # print used_colour_IDs
        for ID in plume_IDs:
            coloured_plumes[SDF_plumes == ID] = used_colour_IDs[ID]

        # print 'Generating colour IDs:', datetime.datetime.now() - test
        clevs = np.arange(0, 41)

        # Need a way of dealing with the colouring of these plumes - the same ID should always be the same colour
        # But that's tricky because you're constantly switching the colour window to account for new IDs

        # So should really assign each ID to a unique colour, like a dictionary or something
        # When an ID is assigned to a colour the colour is removed from a list of available colours
        # When an ID is removed, the colour is re-added to the list of available colours

        # Note that if you just use the existing scheme of no colo
        # test = datetime.datetime.now()
        c1 = m.contourf(lons, lats, coloured_plumes, clevs)
        # print 'Contouring SDF plumes:', datetime.datetime.now() - test
        label_plume_properties(SDF_plumes, date, lons, lats)
        for coll in c1.collections:
            plt.gca().collections.remove(coll)
        # test = datetime.datetime.now()
        c2 = m.contourf(lons, lats, deep_conv_assoc)
        plt.savefig('SDF_plumes/Deep_conv_assoc_' + date + '.png')
        # print 'Contouring deep convection:', datetime.datetime.now() - test
        for coll in c2.collections:
            plt.gca().collections.remove(coll)
        # test = datetime.datetime.now()
        c3 = m.contourf(lons, lats, LLJ_plumes)
        plt.savefig('SDF_plumes/LLJ_plumes_' + date + '.png')
        # print 'Contouring LLJs:', datetime.datetime.now() - test
        for coll in c3.collections:
            plt.gca().collections.remove(coll)
        # test = datetime.datetime.now()
        # c4 = m.contourf(lons, lats, new_plumes)
        # plt.savefig('new_plumes_' + date + '.png')
        # print 'Contouring new plumes:', datetime.datetime.now() - test
        # for coll in c4.collections:
        #    plt.gca().collections.remove(coll)
        # test = datetime.datetime.now()
        SDF_previous = SDF_plumes
        deep_conv_IDs_prev = np.append(deep_conv_IDs_prev, deep_conv_IDs)
        LLJ_plumes_IDs_prev = np.append(LLJ_plumes_IDs_prev, LLJ_IDs)
        k += 1
        print '\nTotal time:', datetime.datetime.now() - totaltest
        """



