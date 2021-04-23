"""
This script provides the computations necessary to replicate
the quantitative analysis that can be found in Section 5 of the paper
"Disaster Risks, Disaster Strikes, and Economic Growth: the Role of Preferences"
published in the Review of Economic Dynamics.

The sript is organized in four parts. If you want to reproduce a given table or
figure of the paper, you can simply go to part III and enter "True" in front
of this table or figure, and then run the code. Table I in part IV allows you
to play with the parameters to test alternative scenarios.
"""


"""
I. Import packages
"""

import numpy as np
from numpy import arange
from pylab import meshgrid

import matplotlib.pyplot as plt
from matplotlib.colors import Normalize

"""
II. Define functions:
"""

class MidpointNormalize(Normalize):
    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
        return np.ma.masked_array(np.interp(value, x, y))


norm = MidpointNormalize(midpoint=0)


def tau(e,g,l,w):

    return (((1-w**(1-g))*l*u)/(a*(1-g)))**(1/(1-u))


def psi(e,g,l,w):

    return e*r+(1-e)*(
        (1-tau(e,g,l,w))*a - (g*(s**2)/2) - l*(1+d-(tau(e,g,l,w))**u) * (1-w**(1-g))/(1-g)
        )


def trend_growth(e,g,l,w):

    return (1-tau(e,g,l,w))*a - psi(e,g,l,w)


def expected_growth(e,g,l,w):

    return trend_growth(e,g,l,w) - l*(1+d-(tau(e,g,l,w))**u) * (1-w)


def effect_disasters_expected_growth(e,g,l,w):

    return expected_growth(e,g,l,w) - expected_growth(e,g,0,w)


def lucas_measure(e,g,l,w):
    psi_optimal = psi(e,g,l,w)
    psi_bau = e*r+(1-e)*(
        (1-0)*a - (g*(s**2)/2) - l*(1+d-(0)**u) * (1-w**(1-g))/(1-g)
        )

    return (psi_optimal/psi_bau)**(1/(1-e)) - 1


def mrt_lambda_gdp(e,g,l,w):
    group_1 = (u**(u/(1-u))-u**(1/(1-u)))/(1-u)
    group_2 = ((1-w**(1-g))/((a**u)*(1-g)))**(1/(1-u))
    group_3 = (1+d)*(1-w**(1-g))/(1-g)

    return -(1/psi(e,g,l,w))*(
            l**(u/(1-u)) * group_1 * group_2 - group_3
            )


def mrt_omega_gdp(e,g,l,w):
    group_1 = (u**(u/(1-u))-u**(1/(1-u)))/(1-u)
    group_2 = ((1-w**(1-g))/(a*(1-g)))**(u/(1-u))

    return -(w**(-g))/psi(e,g,l,w)*(
            l*(1+d) - l**(1/(1-u)) * group_1 * group_2
            )


def rho_to_fit_growth(e,g,l,w,growth_target):

    group_1 = a*(1-tau(e,g,l,w))
    group_2 = l*(1+d-(tau(e,g,l,w))**u) * (1-w)
    group_3 = g*(s**2)/2
    group_4 = l*(1+d-(tau(e,g,l,w))**u) * (1-w**(1-g)) / (1-g)

    return 1/e * (group_1 - group_2 - (1-e)*(group_1 - group_3 - group_4) - growth_target)


"""
III. Choose "True" to display table or figure:
"""

table_2 = True # Main parameters
table_3 = False # MRS lambda
table_4 = False # MRS omega
table_5 = False # Tau (policy instrument)
table_6 = False # Lucas' measure
table_7 = False # Calibration rho
figure_1 = False # Heatmaps (choose scenario in the function)


"""
IV. Apply functions to obtain tables and figures:
"""

# Table I: Parameters used in the calibration (main specification).
g = 3.0
e = 1+1e-09
a = 0.069
w_1 = 0.948
w_2 = 0.85
w_3 = 0.60
l_1 = 0.0307
l_2 = 0.01064
l_3 = 0.003991
d = 1.0
s = 0.02
u = 0.25
growth_target = 0.0175


# Table II: Variables computed at parameters’ baseline value.
if table_2 == True:
    for (l,w,scenario) in [(l_1,w_1,'moderate'),(l_2,w_2,'large'),(l_3,w_3,'extreme')]:
        r = rho_to_fit_growth(e,g,l,w,growth_target)
        print('For', scenario, 'disasters:')
        print('Share of production consumed:', psi(e,g,l,w)/a*100, '%')
        print('Share of production in risk-mitigation:', tau(e,g,l,w)*100, '%')
        print('Reduction in prob. of an env. disaster:', tau(e,g,l,w)**u*100, '%')
        print('Expected growth:', expected_growth(e,g,l,w)*100, '%')
        print('Expected aggregate damages from env. dis. (per year):', l*(1-(tau(e,g,l,w))**u) * (1-w) * 100, '%')
        print('')
        del l,w


# Table III: Marginal rate of substitution between proportionate changes in GDP and in disaster probability.
#if table_3 == True:
#    for gamma in [1+1e-09, 3, 5, 10]:
#        r = rho_to_fit_growth(e,gamma,l_1,w_1,growth_target)
#        print 'with g=', gamma, ', w=', w_1, 'and l=', l_1, ':', mrt_lambda_gdp(e,gamma,l_1,w_1)
#        r = rho_to_fit_growth(e,gamma,l_2,w_2,growth_target)
#        print 'with g=', gamma, ', w=', w_2, 'and l=', l_2, ':', mrt_lambda_gdp(e,gamma,l_2,w_2)
#        r = rho_to_fit_growth(e,gamma,l_3,w_3,growth_target)
#        print 'with g=', gamma, ', w=', w_3, 'and l=', l_3, ':', mrt_lambda_gdp(e,gamma,l_3,w_3)
#
#
## Table IV: Marginal rate of substitution between proportionate changes in GDP and in disaster intensity.
#if table_4 == True:
#    for gamma in [1+1e-09, 3, 5, 10]:
#        r = rho_to_fit_growth(e,gamma,l_1,w_1,growth_target)
#        print 'with g=', gamma, ', w=', w_1, 'and l=', l_1, ':', mrt_omega_gdp(e,gamma,l_1,w_1)
#        r = rho_to_fit_growth(e,gamma,l_2,w_2,growth_target)
#        print 'with g=', gamma, ', w=', w_2, 'and l=', l_2, ':', mrt_omega_gdp(e,gamma,l_2,w_2)
#        r = rho_to_fit_growth(e,gamma,l_3,w_3,growth_target)
#        print 'with g=', gamma, ', w=', w_3, 'and l=', l_3, ':', mrt_omega_gdp(e,gamma,l_3,w_3)
#
#
## Table V: Optimal share of income spent in policy instrument.
#if table_5 == True:
#    for gamma in [1+1e-09, 3, 5, 10]:
#        r = rho_to_fit_growth(e,gamma,l_1,w_1,growth_target)
#        print 'with g=', gamma, ', w=', w_1, 'and l=', l_1, ':', tau(e,gamma,l_1,w_1)*100, '%'
#        r = rho_to_fit_growth(e,gamma,l_2,w_2,growth_target)
#        print 'with g=', gamma, ', w=', w_2, 'and l=', l_2, ':', tau(e,gamma,l_2,w_2)*100, '%'
#        r = rho_to_fit_growth(e,gamma,l_3,w_3,growth_target)
#        print 'with g=', gamma, ', w=', w_3, 'and l=', l_3, ':', tau(e,gamma,l_3,w_3)*100, '%'
#
#
## Table VI: Welfare benefits of the policy.
#if table_6 == True:
#    for gamma in [1+1e-09, 3, 5, 10]:
#        r = rho_to_fit_growth(e,gamma,l_1,w_1,growth_target)
#        print 'with g=', gamma, ', w=', w_1, 'and l=', l_1, ':', lucas_measure(e,gamma,l_1,w_1)*100, '%'
#        r = rho_to_fit_growth(e,gamma,l_2,w_2,growth_target)
#        print 'with g=', gamma, ', w=', w_2, 'and l=', l_2, ':', lucas_measure(e,gamma,l_2,w_2)*100, '%'
#        r = rho_to_fit_growth(e,gamma,l_3,w_3,growth_target)
#        print 'with g=', gamma, ', w=', w_3, 'and l=', l_3, ':', lucas_measure(e,gamma,l_3,w_3)*100, '%'
#
#
## Table VII : Calibration of time impatience to match a 1.75% expected growth rate.
#if table_7 == True:
#    for gamma in [1+1e-09, 3, 5, 10]:
#        for epsilon in [float(1)/3, 1+1e-09, 1.5]:
#            r = rho_to_fit_growth(epsilon,gamma,l_1,w_1,growth_target)
#            print 'with g=', gamma, 'and e=', epsilon, ', w=', w_1, 'and l=', l_1, ':', r
#            r = rho_to_fit_growth(epsilon,gamma,l_2,w_2,growth_target)
#            print 'with g=', gamma, 'and e=', epsilon, ', w=', w_2, 'and l=', l_2, ':', r
#            r = rho_to_fit_growth(epsilon,gamma,l_3,w_3,growth_target)
#            print 'with g=', gamma, 'and e=', epsilon, ', w=', w_3, 'and l=', l_3, ':', r
#
#
## Figure 1: Difference between long-run growth in a disaster vs. disaster free economy.
#if figure_1 == True:
#    # Create the grid
#    e_inverse = 1/(arange(4,1,-0.01) - 0.001)
#    e_normal = arange(1,3,0.01) + 0.001
#    e = np.concatenate([e_inverse,e_normal])
#    g = arange(6,1,-0.01) + 0.005
#    E,G = meshgrid(e, g)
#    r = 0 # r plays no role here. It is necessary to compute expected growth with and wihout disasters
#    # but it cancels out when we take the difference.
#
#    # Prepare axes
#    axe_g = ['6', '5', '4', '3', '2', '1']
#    axe_e = ['1/4', '1/3', '1/2', '1', '2', '3']
#
#    # Build heatmap - choose scenario by selection (l_1,w_1), (l_2,w_2) or (l_3,w_3)
#    values_heatmap = effect_disasters_expected_growth(E, G, l_1, w_1)
#
#    fig, ax = plt.subplots()
#    heatmap_expected_growth = ax.imshow(
#        values_heatmap, norm=norm, cmap=plt.cm.seismic, extent = [0,10,10,0], interpolation='none'
#        )
#    ax.set_xticklabels(axe_e)
#    ax.set_yticklabels(axe_g)
#    plt.xlabel('IES')
#    plt.ylabel('RRA')
#    cbar = plt.colorbar(heatmap_expected_growth)
#    cbar.set_label('percentage points')
#    plt.show()
#
