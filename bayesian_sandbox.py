import matplotlib
matplotlib.use('Agg')
import numpy as np
import scipy.stats as stats
import pymc3 as pm
import theano.tensor as tt
from IPython.core.pylabtools import figsize
import matplotlib.pyplot as plt
from scipy.stats.mstats import mquantiles

def define_random_llj_data(n, varlower, varupper):
    """
    Defines a random set of LLJ occurence data with a predictor variable
    between certain bounds
    :param n:
    :return:
    """

    llj_data = np.zeros((n, 2))
    llj_data[:, 0] = np.random.uniform(50, 100, [n])
    llj_data[:, 1] = np.random.randint(0, 2, [n])

    return llj_data

def logistic(x, beta, alpha):
    return 1.0 / (1.0 + np.exp(np.dot(beta, x) + alpha))

def multivar_logistic(xs, betas, alpha, predictor_count):

    dot_products = np.zeros((betas.shape[0], xs.shape[-2]))
    for i in np.arange(0, predictor_count):
        dot_product = np.dot(betas[..., i], xs[..., i])
        dot_products += dot_product

    """
    print (np.exp(alpha+np.sum([np.dot(betas[:, :, j], xs[:, :, j]) for j in
                                np.arange(
        0, predictor_count)]))/(1. + np.exp(alpha +
        np.sum([np.dot(betas[:, :, j], xs[:, :, j]) for j in np.arange(0,
        predictor_count)])))).shape
    """

    return 1.0 / (1.0 + np.exp(
        alpha + dot_products))

    #return (np.exp(alpha+dot_products))/(1. + np.exp(
    #    alpha+dot_products))

def define_logistic_model(llj_data):
    """
    Defines a model with parameters to derive probabilities of an LLJ event
    for a given value of the predictor variable
    :param llj_data:
    :return:
    """

    predictor = llj_data[:, 0]
    llj_events = llj_data[:, 1]

    with pm.Model() as model:

        # Define probability distribution with a logistic function,
        # with parameters drawn from a normal distribution
        beta = pm.Normal("beta", mu=0, tau=0.001, testval=0)
        alpha = pm.Normal("alpha", mu=0, tau=0.001, testval=0)
        p = pm.Deterministic("p", 1.0/(1. + tt.exp(beta*predictor + alpha)))

        # Use a Bernoulli random variable to connect the probabilities to
        # the observations
        observed = pm.Bernoulli("bernoulli_obs", p, observed=llj_events)
        start = pm.find_MAP()
        step = pm.Metropolis()
        trace = pm.sample(120000, step=step, start=start)
        burned_trace = trace[100000::2]
        alpha_samples = burned_trace["alpha"][:, None]
        beta_samples = burned_trace["beta"][:, None]
        pred_vals = np.linspace(predictor.min(), predictor.max(), 50)[:,
                    None]
        p_pred = logistic(pred_vals.T, beta_samples, alpha_samples)
        mean_prob_pred = p_pred.mean(axis=0)

    return alpha_samples, beta_samples, mean_prob_pred, pred_vals, p_pred

def define_linreg_model(llj_data):
    """
    Defines a model with parameters to derive probabilities of an LLJ event
    for a given value of the predictor variable
    :param llj_data:
    :return:
    """

    predictor = llj_data[:, 0]
    llj_events = llj_data[:, 1]

    with pm.Model() as model:

        # Define probability distribution with a logistic function,
        # with parameters drawn from a normal distribution
        beta = pm.Normal("beta", mu=0, tau=0.001, testval=0)
        alpha = pm.Normal("alpha", mu=0, tau=0.001, testval=0)
        p = pm.Deterministic("p", 1.0/(1. + tt.exp(beta*predictor + alpha)))

        # Use a Bernoulli random variable to connect the probabilities to
        # the observations
        observed = pm.Bernoulli("bernoulli_obs", p, observed=llj_events)
        start = pm.find_MAP()
        step = pm.Metropolis()
        trace = pm.sample(120000, step=step, start=start)
        burned_trace = trace[100000::2]
        alpha_samples = burned_trace["alpha"][:, None]
        beta_samples = burned_trace["beta"][:, None]
        pred_vals = np.linspace(predictor.min()-5, predictor.max()+5, 50)[:,
                    None]
        p_pred = logistic(pred_vals.T, beta_samples, alpha_samples)
        mean_prob_pred = p_pred.mean(axis=0)

    return alpha_samples, beta_samples, mean_prob_pred, pred_vals, p_pred

def define_multivar_logistic_model(llj_data):
    """
    Defines a model with parameters to derive probabilities of an LLJ event
    for a given value of the predictor variable
    :param llj_data:
    :return:
    """

    # For each column, you'll extract alpha and beta normdists (it's
    # probably possible to combine these all into one variable in pm...)

    # Then you get a p value which is simply the same equation as below but
    # with a bunch of variables

    # Then you'll get alpha1 samples, beta1 samples, alpha2 samples,
    # beta2 samples. The point of this process is to infer a probability
    # distribution for the parameters given both your prior and the data you
    #  have.

    predictor_count = llj_data.shape[1]-1
    llj_events = llj_data[:, -1]

    with pm.Model() as model:

        betas = []

        alpha = pm.Normal("alpha", mu=0, tau=0.001, testval=0)

        for i in np.arange(0, predictor_count):
            betas.append(pm.Normal("beta"+str(i), mu=0, tau=0.001, testval=0))

        # Shoddy attempt at a multivariate logistic regression function
        p = pm.Deterministic("p", tt.exp(alpha + tt.sum([betas[j]*llj_data[
        :, j] for j in np.arange(0,  predictor_count)]))/(1. + tt.exp(alpha +
        tt.sum([betas[j]*llj_data[:, j] for j in np.arange(0,
        predictor_count)]))))

        # Use a Bernoulli random variable to connect the probabilities to
        # the observations
        observed = pm.Bernoulli("bernoulli_obs", p, observed=llj_events)
        start = pm.find_MAP()
        step = pm.Metropolis()
        trace = pm.sample(120000, step=step, start=start)
        burned_trace = trace[100000::2]
        alpha_samples = burned_trace["alpha"][:, None]
        beta_samples = burned_trace["beta0"][:, None]

        betas = np.zeros((10000, 1, predictor_count))
        for i in np.arange(0, predictor_count):
            betas[:, :, i] = burned_trace["beta"+str(i)][:, None]#[:, 0]

        predictor_vals = np.zeros((1, 50, predictor_count))
        for i in np.arange(0, predictor_count):
            predictor_vals[:, :, i] = np.linspace(llj_data[:, i].min(),
                                               llj_data[:, i].max(), 50)[:,
                                   None].T#[:, 0]
        p_pred = multivar_logistic(predictor_vals, betas, alpha_samples,
                                   predictor_count)

        mean_prob_pred = p_pred.mean(axis=0)

    return burned_trace, alpha_samples, betas, mean_prob_pred, \
           predictor_vals, p_pred

def define_multivar_logistic_model_2(llj_data)

def plot_posterior_histogram(samples, varname):
    """
    Plots a histogram of sampled data
    :param samples:
    :return:
    """

    figsize(12.5, 6)
    plt.subplot()
    plt.title("Posterior distributions of the variable "+varname)
    plt.hist(samples, histtype='stepfilled', bins=35, alpha=0.85,
             label=r"posterior of "+varname, color="#A60628", normed=True)
    plt.legend();
    plt.grid()
    plt.savefig('posterior_distribution_'+varname+'.png')
    plt.close()

def plot_posterior_expected_probability(mean_prob_pred, pred_vals, p_pred,
                                        predictor, llj_events, varname):
    """
    Plots a curve of the mean expected probability at each value of the
    predictor variable
    :param mean_prob_pred:
    :param pred_vals:
    :param predictor:
    :param llj_events:
    :param varname:
    :return:
    """

    qs = mquantiles(p_pred, [0.025, 0.975], axis=0)

    plt.fill_between(pred_vals[:, 0], *qs, alpha=0.7, color="#7A68A6")

    plt.plot(pred_vals[:, 0], qs[0], label="95% CI", color="#7A68A6",
             alpha=0.7)
    plt.plot(pred_vals, mean_prob_pred, lw=1, ls="--", color="k",
    label="average posterior \nprobability of event")
    plt.scatter(predictor, llj_events, color="k", s=50, alpha=0.5)
    plt.title("Posterior probability estimates")
    plt.legend(loc="lower left")
    plt.ylim(-0.02, 1.02)
    plt.xlim(pred_vals.min(), pred_vals.max())
    plt.ylabel("probability estimate")
    plt.xlabel(varname)
    plt.grid()
    plt.savefig("posterior_probability_estimates.png")
    plt.close()

if __name__ == '__main__':


    univariate_llj_data = np.zeros((25, 2))
    multivariate_llj_data = np.zeros((25, 3))

    llj_data_1 = define_random_llj_data(25, 50, 100)
    llj_data_2 = define_random_llj_data(25, 0.5, 1.5)

    multivariate_llj_data[:, 0] = llj_data_1[:, 0]
    multivariate_llj_data[:, 1] = llj_data_2[:, 0]
    multivariate_llj_data[:, 2] = llj_data_2[:, 1]

    multivariate_llj_data[:, 0] = np.arange(1, 26)
    multivariate_llj_data[:, 1] = np.arange(26, 51)

    multivariate_llj_data[1:10, 2] = 0
    multivariate_llj_data[10:, 2] = 1

    univariate_llj_data[1:10, 1] = 0
    univariate_llj_data[10:, 1] = 1
    univariate_llj_data[:, 0] = np.arange(1, 26)

    #alpha_samples, beta_samples, mean_prob_pred, pred_vals, p_pred = \
    #    define_logistic_model(univariate_llj_data)


    #plot_posterior_histogram(alpha_samples, "alpha")
    #plot_posterior_histogram(beta_samples, "beta")

    #plot_posterior_expected_probability(mean_prob_pred, pred_vals, p_pred,
    #                                    univariate_llj_data[:, 0],
    #                                    univariate_llj_data[:, 1],
    #                                   'var')

    burned_trace, alpha_samples, betas, mean_prob_pred, pred_vals, \
    p_pred \
        = define_multivar_logistic_model(multivariate_llj_data)

    """
    beta0 = burned_trace["beta0"][:, None]
    beta1 = burned_trace["beta1"][:, None]

    plot_posterior_histogram(alpha_samples, "alpha")
    plot_posterior_histogram(beta0, "beta0")
    plot_posterior_histogram(beta1, "beta1")

    print 'Probability of an LLJ'
    xs = np.zeros((1, 1, 2))
    xs[:, :, 0] = 0
    xs[:, :, 1] = 0

    ps = multivar_logistic(xs, betas, alpha_samples, 2)
    print np.nanmean(ps)

    ps1 = logistic(25, beta0, alpha_samples)
    print np.nanmean(ps1)

    plot_posterior_expected_probability(mean_prob_pred, pred_vals[:, :, 0].T,
                                        p_pred,
                                        multivariate_llj_data[:, 0],
                                        multivariate_llj_data[:, 2],
                                        "sample variable (one of two)")
    """




