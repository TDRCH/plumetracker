from netCDF4 import Dataset
from netCDF4 import num2date
import pandas as pd
import numpy as np
import pymc3 as pm
import matplotlib.pyplot as plt
import seaborn
import plotting
import datetime as dt

"""
Works out the most likely CET frost day from a Bayesian linear regression model
"""

def extract_frost_variables(filename):
    """
    From CET, works out the climatological mean first frost day, last frost
    day and a bunch of associated variables
    :return:
    """

    cet = Dataset(filename)
    cetmin = cet.variables['min_CET'][:]
    dates = num2date(cet.variables['time'][:], cet.variables['time'].units)
    cetmean = cet.variables['mean_CET'][:]
    cetmax = cet.variables['max_CET'][:]

    # Calculate the climatological first frost date
    # This is the first day in the second half of the year with a CETmin
    # less than zero

    first_frost_days = []
    last_frost_days = []
    mean_summer_temps = []
    mean_spring_temps = []
    mean_sept_temps = []
    mean_3rd_week_Oct_temps = []
    years = []

    year_count = (dates[-1].year-dates[0].year)+1

    # Extract yearly statistics
    for i in np.arange(0, year_count):

        # No CETmin data exists prior to 1878
        if dates[0].year+i < 1878:
            continue

        # Get the first frost day for the year
        year_bool = np.asarray([j.year == dates[0].year+i for j in dates])
        years.append(dates[0].year+i)
        year_cetmin = cetmin[year_bool]
        year_dates = dates[year_bool]
        h2_bool = np.asarray([j.month >= 7 for j in year_dates])
        year_cetmin_h2 = year_cetmin[h2_bool]
        year_dates_h2 = year_dates[h2_bool]
        frost_bool = np.asarray([j <= 0 for j in year_cetmin_h2])
        first_frost_day = year_dates_h2[frost_bool][0]
        first_frost_days.append(first_frost_day)

        if dates[0].year+i < 2016:
            # Get the last frost day for the following year
            nextyear_bool = np.asarray([j.year == dates[0].year+i+1 for j in
                                        dates])
            nextyear_cetmin = cetmin[nextyear_bool]
            nextyear_dates = dates[nextyear_bool]
            h1_bool = np.asarray([j.month < 7 for j in nextyear_dates])
            nextyear_cetmin_h1 = nextyear_cetmin[h1_bool]
            nextyear_dates_h1 = nextyear_dates[h1_bool]
            nextfrost_bool = np.asarray([j <= 0 for j in nextyear_cetmin_h1])
            last_frost_day = nextyear_dates_h1[nextfrost_bool][-1]
            last_frost_days.append(last_frost_day)
        else:
            last_frost_days.append(np.nan)

        year_cetmean = cetmean[year_bool]
        # Try with CET min
        year_cetmean = cetmin[year_bool]

        # Get the mean summer temperature for the year
        summer_bool = np.asarray([j.month in [6, 7, 8] for j in year_dates])
        summer_cetmean = np.nanmean(year_cetmean[summer_bool])
        mean_summer_temps.append(summer_cetmean)

        # Get the mean spring temperature for the year
        spring_bool = np.asarray([j.month in [3, 4, 5] for j in year_dates])
        spring_cetmean = np.nanmean(year_cetmean[spring_bool])
        mean_spring_temps.append(spring_cetmean)

        # Get the mean September temperature for the year
        sept_bool = np.asarray([j.month == 9 for j in year_dates])
        september_cetmean = np.nanmean(year_cetmean[sept_bool])
        mean_sept_temps.append(september_cetmean)

        # Get the mean 3rd week October temperature for the year
        octw3_bool = np.asarray([j.month == 10 and j.day >= 14 and j.day <=
                                 21 for j in year_dates])
        octw3_cetmean = np.nanmean(year_cetmean[octw3_bool])
        mean_3rd_week_Oct_temps.append(octw3_cetmean)

    # Get first frost days of year
    first_frost_doys = np.asarray([j.timetuple().tm_yday for j in
                               first_frost_days])

    # Convert the last frost day in the last year to something invalid
    last_frost_days[-1] = dt.datetime.now()

    # Get last frost days of year
    last_frost_doys = np.asarray([j.timetuple().tm_yday for j in
                               last_frost_days]).astype(float)

    last_frost_doys[-1] = np.nan

    # The actual year count is the number of years for which we have data
    year_count_with_data = len(first_frost_doys)

    # Calculate the climatological mean first frost day
    mean_first_frost = np.round(np.nanmean(first_frost_doys))

    # Calculate the climatological mean last frost day
    mean_last_frost = np.round(np.nanmean(last_frost_doys))

    first_frost_diffs = []
    last_frost_diffs = []

    # Calculate the distance in days of each year's first frost day from the
    # climatological mean
    for i in np.arange(0, year_count_with_data):
        first_frost = first_frost_doys[i]
        last_frost = last_frost_doys[i]
        first_frost_diff = first_frost-mean_first_frost
        first_frost_diffs.append(first_frost_diff)
        last_frost_diff = last_frost-mean_last_frost
        last_frost_diffs.append(last_frost_diff)

    d = {'first_frost_diff': first_frost_diffs, 'last_frost_diff':
        last_frost_diffs,
                                                'mean_summer_temp':
        mean_summer_temps, 'mean_spring_temp': mean_spring_temps,
         'mean_september_temp': mean_sept_temps, 'mean_octw3_temps':
             mean_3rd_week_Oct_temps}

    index = years
    cetdf = pd.DataFrame(data=d, index=index)

    return cetdf, mean_first_frost, mean_last_frost

def linear_regression_model(cetdf, dpdtvarname, indpdtvarname1,
                            indpdtvarname2):

    with pm.Model() as model:
        # String is R syntax: want income as a function of the following
        # variables
        pm.glm.GLM.from_formula(dpdtvarname+' ~ '+indpdtvarname1 + ' + ' +
                                 indpdtvarname2,
                                cetdf)
        trace_model = pm.sample(4000)

        plotting.plot_traces(trace_model, retain=1000)

    plt.figure(figsize=(9, 7))
    trace = trace_model[1000:]

    """
    seaborn.jointplot(trace[indpdtvarname1], trace[indpdtvarname2],
                      kind="hex",
                      color="#4CB391")
    plt.xlabel("beta_mean_octw3_temps")
    plt.ylabel("beta_mean_summer_temp")
    plt.savefig('cet_probability_jointplot.png')
    plt.close()
    """

    return trace

def linreg(trace, indpdtvar1name, indpdtvar2name, indpdtvar1, indpdtvar2):
    return trace[indpdtvar1name] * indpdtvar1 \
           + trace[indpdtvar2name] * indpdtvar2 + trace['Intercept']

if __name__ == '__main__':

    cetdf, mean_first_frost, mean_last_frost = extract_frost_variables(
        'HadCET2.nc')
    cetdf = cetdf[~pd.isnull(cetdf['last_frost_diff'])]

    g = seaborn.pairplot(cetdf)
    plt.savefig('cet_pairplot.png')
    plt.close()

    # Compute the correlation matrix
    corr = cetdf.corr()

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
    plt.savefig('cet_correlheatmap.png')
    plt.close()



    trace = linear_regression_model(cetdf,
                                    'first_frost_diff', 'mean_octw3_temps',
                                    'mean_summer_temp')

    summer_cetmean_2017 = 16.1
    octw3_cetmean_2017 = np.nanmean([11.7, 11.1, 10.7, 11.8, 12.0, 10.5,
                                     12.9, 14.2])

    mlff_offset = linreg(trace,
                                    'mean_octw3_temps',
                                    'mean_summer_temp', octw3_cetmean_2017,
                                 summer_cetmean_2017)
    mlff_offset_mean = np.nanmean(mlff_offset)

    mlff_offset_p5 = np.percentile(mlff_offset, 5)
    mlff_offset_p95 = np.percentile(mlff_offset, 95)

    most_likely_ff = np.round(mean_first_frost + mlff_offset_mean)
    p5_ff = np.round(mean_first_frost + mlff_offset_p5)
    p95_ff = np.round(mean_first_frost + mlff_offset_p95)

    most_likely_ff_date = dt.datetime(2017, 1, 1) + dt.timedelta(
        most_likely_ff - 1)
    p5_ff_date = dt.datetime(2017, 1, 1) + dt.timedelta(
        p5_ff - 1)
    p95_ff_date = dt.datetime(2017, 1, 1) + dt.timedelta(
        p95_ff - 1)

    print 'Most likely first frost date:'
    print most_likely_ff_date
    print 'P5% first frost date:'
    print p5_ff_date
    print 'P95% first frost date:'
    print p95_ff_date

    trace2 = linear_regression_model(cetdf,
                                    'last_frost_diff', 'mean_summer_temp',
                                    'mean_september_temp')

    summer_cetmean_2017 = 16.1
    mean_september_cetmean_2017 = 13.5

    mllf_offset = linreg(trace2,
                                    'mean_summer_temp',
                                    'mean_september_temp', summer_cetmean_2017,
                         mean_september_cetmean_2017)
    mllf_offset_mean = np.nanmean(mllf_offset)

    mllf_offset_p5 = np.percentile(mllf_offset, 5)
    mllf_offset_p95 = np.percentile(mllf_offset, 95)

    most_likely_lf = np.round(mean_last_frost + mllf_offset_mean)
    p5_lf = np.round(mean_last_frost + mllf_offset_p5)
    p95_lf = np.round(mean_last_frost + mllf_offset_p95)

    most_likely_lf_date = dt.datetime(2018, 1, 1) + dt.timedelta(
        most_likely_lf - 1)
    p5_lf_date = dt.datetime(2018, 1, 1) + dt.timedelta(
        p5_lf - 1)
    p95_lf_date = dt.datetime(2018, 1, 1) + dt.timedelta(
        p95_lf - 1)

    print 'Most likely last frost date:'
    print most_likely_lf_date
    print 'P5% last frost date:'
    print p5_lf_date
    print 'P95% last frost date:'
    print p95_lf_date









