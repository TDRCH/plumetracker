import pandas as pd
import numpy as np
import pymc3 as pm
import matplotlib.pyplot as plt
import seaborn
import warnings
warnings.filterwarnings('ignore')
from collections import OrderedDict
from time import time
from scipy.optimize import fmin_powell
from scipy import integrate
import theano as thno
import theano.tensor as T

def run_models(df, upper_order=5):
    """
    Convenience function:
    Fit a range of pymc3 models of increasing polynomial complexity.
    Suggest limit to max order 5 since calculation time is exponential.
    """

    models, traces = OrderedDict(), OrderedDict()

    for k in range(1,upper_order+1):

        nm = 'k{}'.format(k)
        fml = create_poly_modelspec(k)

        with pm.Model() as models[nm]:

            print('\nRunning: {}'.format(nm))
            pm.glm.GLM.from_formula(fml, df, family=pm.glm.families.Normal())

            traces[nm] = pm.sample(2000, init=None)

    return models, traces

def plot_traces(traces, retain=1000):
    """
    Convenience function:
    Plot traces with overlaid means and values
    """

    ax = pm.traceplot(traces[-retain:], figsize=(12,len(traces.varnames)*1.5),
        lines={k: v['mean'] for k, v in pm.df_summary(traces[-retain:]).
                      iterrows()})

    for i, mn in enumerate(pm.df_summary(traces[-retain:])['mean']):
        ax[i,0].annotate('{:.2f}'.format(mn), xy=(mn,0), xycoords='data'
                    ,xytext=(5,10), textcoords='offset points', rotation=90
                    ,va='bottom', fontsize='large', color='#AA0022')

def create_poly_modelspec(k=1):
    """
    Convenience function:
    Create a polynomial modelspec string for patsy
    """
    return ('income ~ educ + hours + age ' + ' '.join(['+ np.power(age,{})'.
                                                      format(j)
                                     for j in range(2,k+1)])).strip()

if __name__ == '__main__':
    data = pd.read_csv(
        "https://archive.ics.uci.edu/ml/machine-learning-databases/adult/"
        "adult.data",
        header=None, names=['age', 'workclass', 'fnlwgt',
                            'education-categorical', 'educ',
                            'marital-status', 'occupation',
                            'relationship', 'race', 'sex',
                            'captial-gain', 'capital-loss',
                            'hours', 'native-country',
                            'income'])

    data.head(10)
    data = data[~pd.isnull(data['income'])]
    data[data['native-country'] == " United-States"]
    income = 1 * (data['income'] == " >50K")
    age2 = np.square(data['age'])

    data = data[['age', 'educ', 'hours']]
    data['age2'] = age2
    data['income'] = income

    g = seaborn.pairplot(data)
    plt.savefig('pairplot.png')
    plt.close()

    # Compute the correlation matrix
    corr = data.corr()

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
    plt.savefig('correlheatmap.png')
    plt.close()

    with pm.Model() as logistic_model:
        # String is R syntax: want income as a function of the following
        # variables
        pm.glm.GLM.from_formula('income ~ age + age2 + educ + hours', data,
                                family = pm.glm.families.Binomial())
        trace_logistic_model = pm.sample(4000)

        plot_traces(trace_logistic_model, retain=1000)

    plt.figure(figsize=(9, 7))
    trace = trace_logistic_model[1000:]
    seaborn.jointplot(trace['age'], trace['educ'], kind="hex", color="#4CB391")
    plt.xlabel("beta_age")
    plt.ylabel("beta_educ")
    plt.savefig('probability_jointplot.png')
    plt.close()

    # Linear model with hours == 50 and educ == 12 and age == 30
    lm = lambda x, samples: 1 / (1 + np.exp(-(samples['Intercept'] +
                                              samples['age'] * 30 +
                                              samples['age2'] * np.square(30) +
                                              samples['educ'] * 12 +
                                              samples['hours'] * 50)))

    # This should be a distribution of values across our parameter space

    # Ok so now there are two ways we could go: make some dummy LLJ data and
    #  test this method, or try and figure out what on earth it's doing



