#! /usr/bin/env python

import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['text.usetex'] = True
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.serif'] = 'Computer Modern Roman'
#mpl.use('Agg')
fsize = 32

from numpy import pi, deg2rad, linspace, meshgrid, sqrt, arctan2
from ks_lib import omega


# material parameters
Ms = 139260.5752 #A/m
H = 55704.230082 #A/m
L = 5.1e-06 #m
gamma = 1.760859708e11 #rad/(sT)
alpha = 3.1e-16 #m^2
theta = deg2rad(90.0) #rad
trans_n = 0.0
mp = [Ms, H, L, gamma, alpha, theta, trans_n]

#k = np.logspace(0,7)
kmag = 2e6
pts = 1000
k = 2*pi*linspace(-kmag, kmag, pts) #meters


#summary figure number
##########
sfn= 1
##########
plt.figure(sfn, figsize = (12,9))
plt.title('Dispersion Surface')
plt.xlabel(r'$k_{\parallel}$, m$^{-1}$', fontsize = fsize, labelpad = 20)
plt.ylabel(r'$k_{\perp}$, m$^{-1}$', fontsize = fsize, labelpad = 20)

plt.minorticks_on()
plt.tick_params(labelsize = fsize, which = 'both', right = 'off', top = 'off', pad = -10)
plt.tick_params(which = 'major', direction = 'out', length = 8, width = 3)
plt.tick_params(which = 'minor', direction = 'out', length = 5, width = 2)

kpara, kperp = meshgrid(k, k)
kzeta = sqrt(kperp**2 + kpara**2)
phi = arctan2(kperp,kpara)
vo = 1e-9*omega(kzeta,phi,mp)/(2*pi)
disp_surface_limits = [ kpara.min(), kpara.max(), kperp.min(), kperp.max() ]
plt.imshow(vo, extent = disp_surface_limits)

##########
sfn=2
##########
plt.figure(sfn, figsize = (12,9))
plt.title('Slowness Curves')
plt.xlabel(r'$k_{\parallel}$, m$^{-1}$', fontsize = fsize, labelpad = 20)
plt.ylabel(r'$k_{\perp}$, m$^{-1}$', fontsize = fsize, labelpad = 20)

plt.minorticks_on()
plt.tick_params(labelsize = fsize, which = 'both', right = 'off', top = 'off', pad = -10)
plt.tick_params(which = 'major', direction = 'out', length = 8, width = 3)
plt.tick_params(which = 'minor', direction = 'out', length = 5, width = 2)

CS = plt.contour(k, k, vo)
plt.clabel(CS, inline=1, fontsize=20)

##########
sfn=3
##########
plt.figure(sfn, figsize = (12,9))
plt.title('BV Dispersion')
plt.xlabel(r'$|k|$, m$^{-1}$', fontsize = fsize, labelpad = 20)
plt.ylabel(r'$f$, GHz', fontsize = fsize, labelpad = 20)

plt.minorticks_on()
plt.tick_params(labelsize = fsize, which = 'both', right = 'off', top = 'off', pad = -10)
plt.tick_params(which = 'major', direction = 'out', length = 8, width = 3)
plt.tick_params(which = 'minor', direction = 'out', length = 5, width = 2)

phi=0
k_BV=2*pi*linspace(0,kmag,pts)
plt.plot(k_BV,1e-9*omega(k_BV,phi,mp)/(2*pi))

##########
sfn=4
##########
plt.figure(sfn, figsize = (12,9))
plt.figure(sfn, figsize = (12,9))
plt.title('DE Dispersion')
plt.xlabel(r'$|k|$, m$^{-1}$', fontsize = fsize, labelpad = 20)
plt.ylabel(r'$f$, GHz', fontsize = fsize, labelpad = 20)

plt.minorticks_on()
plt.tick_params(labelsize = fsize, which = 'both', right = 'off', top = 'off', pad = -10)
plt.tick_params(which = 'major', direction = 'out', length = 8, width = 3)
plt.tick_params(which = 'minor', direction = 'out', length = 5, width = 2)


phi=90
k_DE=2*pi*linspace(0,kmag,pts)
plt.plot(k_DE,1e-9*omega(k_DE,phi,mp)/(2*pi))

plt.tight_layout()
plt.show()
