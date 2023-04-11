import numpy as np
from numpy import exp, linspace, pi, random, sign, sin
from scipy import interpolate
from scipy.optimize import minimize
import emcee
import matplotlib.pyplot as plt
import pandas as pd


'''
p_true = Parameters()
p_true.add('amp', value=14.0)
p_true.add('period', value=5.46)
p_true.add('shift', value=0.123)
p_true.add('decay', value=0.032)


def residual(pars, x, data=None):
    """Model a decaying sine wave and subtract data."""
    vals = pars.valuesdict()
    print(vals)
    amp = vals['amp']
    per = vals['period']
    shift = vals['shift']
    decay = vals['decay']

    if abs(shift) > pi/2:
        shift = shift - sign(shift)*pi
    model = amp * sin(shift + x/per) * exp(-x*x*decay*decay)
    if data is None:
        return model
    return model - data


#random.seed(0)
x = linspace(0.0, 250., 1001)
noise = random.normal(scale=0.7215, size=x.size)
data = residual(p_true, x) + noise

fit_params = Parameters()
fit_params.add('amp', value=13.0)
fit_params.add('period', value=2)
fit_params.add('shift', value=0.0)
fit_params.add('decay', value=0.02)

out = minimize(residual, fit_params, method='basinhopping', args=(x,), kws={'data': data})

print(fit_report(out))

tgrid = np.linspace(4500, 6500, 21)
ggrid = np.linspace(4.0, 5.0, 6)
zgrid = np.linspace(-1., 1., 11)


params = Parameters()
params.add('Temp', value=5815)
params.add('g', value=4.45)
params.add('z', value=0.1)

def closest_param(param_list, param):
    index = (abs(param_list - param)).argmin()

    return param_list[index]

new_param = closest_param(tgrid, params['Temp'])
print(new_param)


# Choose the "true" parameters.
m_true = -0.9594
b_true = 4.294
f_true = 0.534

# Generate some synthetic data from the model.
N = 50
x = np.sort(10 * np.random.rand(N))
yerr = 0.1 + 0.5 * np.random.rand(N)
y = m_true * x + b_true
y += np.abs(f_true * y) * np.random.randn(N)
y += yerr * np.random.randn(N)

def log_likelihood(theta, x, y, yerr):

    m, b, log_f = theta
    model = m * x + b
    sigma2 = yerr**2 + model**2 * np.exp(2 * log_f)
    return -0.5 * np.sum((y - model) ** 2 / sigma2 + np.log(sigma2))

nll = lambda *args: -log_likelihood(*args)
initial = np.array([m_true, b_true, np.log(f_true)]) + 0.1 * np.random.randn(3)
soln = minimize(nll, initial, args=(x, y, yerr))
m_ml, b_ml, log_f_ml = soln.x

def log_prior(theta):
    print(theta)
    m, b, log_f = theta
    if -5.0 < m < 0.5 and 0.0 < b < 10.0 and -10.0 < log_f < 1.0:
        return 0.0
    return -np.inf

def log_probability(theta, x, y, yerr):
    lp = log_prior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + log_likelihood(theta, x, y, yerr)

pos = soln.x + 1e-4 * np.random.randn(10, 3)
nwalkers, ndim = 10, 3
#print(pos)

pos = np.array([np.array([-1.0, 4.4, 0.0]) + 1e-4 * np.random.randn(3) for i in range(10)])
print(pos)


sampler = emcee.EnsembleSampler(nwalkers, ndim, log_probability, args=(x, y, yerr))
sampler.run_mcmc(pos, 5, progress=True)

fig, axes = plt.subplots(3, figsize=(10, 7), sharex=True)
samples = sampler.get_chain()
labels = ["m", "b", "log(f)"]
for i in range(ndim):
    ax = axes[i]
    ax.plot(samples[:, :, i], "k", alpha=0.3)
    ax.set_xlim(0, len(samples))
    ax.set_ylabel(labels[i])
    ax.yaxis.set_label_coords(-0.1, 0.5)

axes[-1].set_xlabel("step number");
plt.show()
'''

spec_5000 = pd.read_csv('/home/fmendez/Turbospectrum2019/Utilities/model_5527g+4.5z+0.10vt01.TURBO_5000-5500_xit1.10.spec_rot', delimiter='\s+', names=['Wave','Flux','Int'], engine='python')
spec_5500 = pd.read_csv('/home/fmendez/Turbospectrum2019/Utilities/model_5527g+4.5z+0.10vt01.TURBO_5500-6000_xit1.10.spec_rot', delimiter='\s+', names=['Wave','Flux','Int'], engine='python')
interpolated_df = pd.read_csv('/home/fmendez/Desktop/interpol_T5527_g+4.5_z+0.10_vt1.10_vr-4.5_Na6.30_Mg7.15.spec', engine='python')

actual_spec = pd.concat([spec_5000, spec_5500],names=['Wave','Flux','Int']).reset_index(drop=True)
print(interpolated_df)
print(actual_spec)

plt.plot(actual_spec['Wave'], actual_spec['Flux'], 'r-', label='Actual Spectrum')
plt.plot(interpolated_df['Wave'], interpolated_df['Flux'], 'k-', label='Interpolated Spectrum')
plt.plot(actual_spec['Wave'], (interpolated_df['Flux'] - actual_spec['Flux']), 'b-', label='Residual')
plt.legend()
plt.show()
