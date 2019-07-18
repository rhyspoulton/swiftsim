###############################################################################
 # This file is part of SWIFT.
 # Copyright (c) 2015 Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 #                    Matthieu Schaller (matthieu.schaller@durham.ac.uk)
 # 
 # This program is free software: you can redistribute it and/or modify
 # it under the terms of the GNU Lesser General Public License as published
 # by the Free Software Foundation, either version 3 of the License, or
 # (at your option) any later version.
 # 
 # This program is distributed in the hope that it will be useful,
 # but WITHOUT ANY WARRANTY; without even the implied warranty of
 # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 # GNU General Public License for more details.
 # 
 # You should have received a copy of the GNU Lesser General Public License
 # along with this program.  If not, see <http://www.gnu.org/licenses/>.
 # 
 ##############################################################################

# Computes the analytical solution of the 3D Sedov blast wave.
# The script works for a given initial box and dumped energy and computes the solution at a later time t.

# Parameters
rho_0 = 1.          # Background Density
P_0 = 1.e-6         # Background Pressure
E_0 = 1.0           # Energy of the explosion
gas_gamma = 5./3.   # Gas polytropic index


# ---------------------------------------------------------------
# Don't touch anything after this.
# ---------------------------------------------------------------

import matplotlib
matplotlib.use("Agg")
from pylab import *
from scipy import stats
import h5py

from swiftsimio import load
import unyt

# Plot parameters
params = {'axes.labelsize': 10,
'axes.titlesize': 10,
'font.size': 12,
'legend.fontsize': 12,
'xtick.labelsize': 10,
'ytick.labelsize': 10,
'text.usetex': True,
 'figure.figsize' : (9.90,6.45),
'figure.subplot.left'    : 0.045,
'figure.subplot.right'   : 0.99,
'figure.subplot.bottom'  : 0.05,
'figure.subplot.top'     : 0.99,
'figure.subplot.wspace'  : 0.15,
'figure.subplot.hspace'  : 0.12,
'lines.markersize' : 6,
'lines.linewidth' : 3.,
'text.latex.unicode': True
}
rcParams.update(params)
rc('font',**{'family':'sans-serif','sans-serif':['Times']})


snap = int(sys.argv[1])


# Read the simulation data
sim = load("feedback_%04d.hdf5"%snap)
boxSize = sim.metadata.boxsize[0]
time = sim.metadata.t.to(unyt.Gyr)
scheme = sim.metadata.hydro_scheme["Scheme"]
kernel = sim.metadata.hydro_scheme["Kernel function"]
neighbours = sim.metadata.hydro_scheme["Kernel target N_ngb"]
eta = sim.metadata.hydro_scheme["Kernel eta"]
git = sim.metadata.code["Git Revision"]

pos = sim.gas.coordinates
x = pos[:,0] - boxSize / 2
y = pos[:,1] - boxSize / 2
z = pos[:,2] - boxSize / 2
vel = sim.gas.velocities
r = sqrt(x**2 + y**2 + z**2)
v_r = ((x * vel[:,0] + y * vel[:,1] + z * vel[:,2]) / r).to(unyt.km / unyt.s)
r.convert_to_units(unyt.kpc)
u = sim.gas.temperature.to(unyt.K)
S = sim.gas.entropy
P = sim.gas.pressure
rho = sim.gas.density.to(unyt.mh / (unyt.cm ** 3))

try:
    diffusion = sim.gas.diffusion 
    plot_diffusion = True
except:
    plot_diffusion = False

try:
    viscosity = sim.gas.viscosity
    plot_viscosity = True
except:
    plot_viscosity = False

# Bin te data
r_bin_edge = np.linspace(0., boxSize.to(unyt.kpc).value, 25)
r_bin = 0.5*(r_bin_edge[1:] + r_bin_edge[:-1])
rho_bin,_,_ = stats.binned_statistic(r, rho, statistic='mean', bins=r_bin_edge)
v_bin,_,_ = stats.binned_statistic(r, v_r, statistic='mean', bins=r_bin_edge)
P_bin,_,_ = stats.binned_statistic(r, P, statistic='mean', bins=r_bin_edge)
S_bin,_,_ = stats.binned_statistic(r, S, statistic='mean', bins=r_bin_edge)
u_bin,_,_ = stats.binned_statistic(r, u, statistic='mean', bins=r_bin_edge)
rho2_bin,_,_ = stats.binned_statistic(r, rho**2, statistic='mean', bins=r_bin_edge)
v2_bin,_,_ = stats.binned_statistic(r, v_r**2, statistic='mean', bins=r_bin_edge)
P2_bin,_,_ = stats.binned_statistic(r, P**2, statistic='mean', bins=r_bin_edge)
S2_bin,_,_ = stats.binned_statistic(r, S**2, statistic='mean', bins=r_bin_edge)
u2_bin,_,_ = stats.binned_statistic(r, u**2, statistic='mean', bins=r_bin_edge)
rho_sigma_bin = np.sqrt(rho2_bin - rho_bin**2)
v_sigma_bin = np.sqrt(v2_bin - v_bin**2)
P_sigma_bin = np.sqrt(P2_bin - P_bin**2)
S_sigma_bin = np.sqrt(S2_bin - S_bin**2)
u_sigma_bin = np.sqrt(u2_bin - u_bin**2)

if plot_diffusion:
    alpha_diff_bin,_,_ = stats.binned_statistic(r, diffusion, statistic='mean', bins=r_bin_edge)
    alpha2_diff_bin,_,_ = stats.binned_statistic(r, diffusion**2, statistic='mean', bins=r_bin_edge)
    alpha_diff_sigma_bin = np.sqrt(alpha2_diff_bin - alpha_diff_bin**2)

if plot_viscosity:
    alpha_visc_bin,_,_ = stats.binned_statistic(r, viscosity, statistic='mean', bins=r_bin_edge)
    alpha2_visc_bin,_,_ = stats.binned_statistic(r, viscosity**2, statistic='mean', bins=r_bin_edge)
    alpha_visc_sigma_bin = np.sqrt(alpha2_visc_bin - alpha_visc_bin**2)


# Now, work our the solution....

from scipy.special import gamma as Gamma
from numpy import *


# Plot the interesting quantities
figure()

# Velocity profile --------------------------------
subplot(231)
plot(r, v_r, '.', color='r', ms=0.5, alpha=0.2)
errorbar(r_bin, v_bin, yerr=v_sigma_bin, fmt='.', ms=8.0, color='b', lw=1.2)
xlabel("${\\rm{Radius}}~r$", labelpad=0)
ylabel("${\\rm{Radial~velocity}}~v_r$", labelpad=0)

# Density profile --------------------------------
subplot(232)
plot(r, rho, '.', color='r', ms=0.5, alpha=0.2)
errorbar(r_bin, rho_bin, yerr=rho_sigma_bin, fmt='.', ms=8.0, color='b', lw=1.2)
xlabel("${\\rm{Radius}}~r$", labelpad=0)
ylabel("${\\rm{Density}}~\\rho$", labelpad=2)

# Pressure profile --------------------------------
subplot(233)
plot(r, P, '.', color='r', ms=0.5, alpha=0.2)
errorbar(r_bin, P_bin, yerr=P_sigma_bin, fmt='.', ms=8.0, color='b', lw=1.2)
xlabel("${\\rm{Radius}}~r$", labelpad=0)
ylabel("${\\rm{Pressure}}~P$", labelpad=0)

# Internal energy profile -------------------------
subplot(234)
plot(r, u, '.', color='r', ms=0.5, alpha=0.2)
errorbar(r_bin, u_bin, yerr=u_sigma_bin, fmt='.', ms=8.0, color='b', lw=1.2)
xlabel("${\\rm{Radius}}~r$", labelpad=0)
ylabel("${\\rm{Temperature}}~u$", labelpad=0)

# Entropy profile ---------------------------------
subplot(235)
xlabel("${\\rm{Radius}}~r$", labelpad=0)


if plot_diffusion or plot_viscosity:
    if plot_diffusion:
        plot(r, diffusion, ".", color='r', ms=0.5, alpha=0.2)
        errorbar(r_bin, alpha_diff_bin, yerr=alpha_diff_sigma_bin, fmt=".", ms=8.0, color='b', lw=1.2, label="Diffusion")

    if plot_viscosity:
        plot(r, viscosity, ".", color='g', ms=0.5, alpha=0.2)
        errorbar(r_bin, alpha_visc_bin, yerr=alpha_visc_sigma_bin, fmt=".", ms=8.0, color='y', lw=1.2, label="Viscosity")

    ylabel("${\\rm{Rate~Coefficient}}~\\alpha$", labelpad=0)
    legend()
else:
    plot(r, S, '.', color='r', ms=0.5, alpha=0.2)
    errorbar(r_bin, S_bin, yerr=S_sigma_bin, fmt='.', ms=8.0, color='b', lw=1.2)
    ylabel("${\\rm{Entropy}}~S$", labelpad=0)

# Information -------------------------------------
subplot(236, frameon=False)

text(-0.49, 0.9, "Feedback with  $\\gamma=%.3f$ in 3D at $t=%.2f$"%(gas_gamma,time), fontsize=10)
text(-0.49, 0.8, "Background $\\rho_0=%.2f$"%(rho_0), fontsize=10)
text(-0.49, 0.7, "Energy injected $E_0=%.2f$"%(E_0), fontsize=10)
plot([-0.49, 0.1], [0.62, 0.62], 'k-', lw=1)
text(-0.49, 0.5, "$\\textsc{Swift}$ %s"%git, fontsize=10)
text(-0.49, 0.4, scheme, fontsize=10)
text(-0.49, 0.3, kernel, fontsize=10)
text(-0.49, 0.2, "$%.2f$ neighbours ($\\eta=%.3f$)"%(neighbours, eta), fontsize=10)
xlim(-0.5, 0.5)
ylim(0, 1)
xticks([])
yticks([])


savefig("Feedback.png", dpi=200)




