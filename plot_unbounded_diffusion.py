#!/usr/bin/env python
# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import numpy as np

plt.rc('text', usetex=True)

font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 10}

figsize = {'figsize': (3.5, 3.5)}

plt.rc('font', **font)
plt.rc('figure', **figsize)

def radial_bins(x_pos, y_pos, z_pos, nbins=120, r_max=1.2, offset=(0,0,0)):
    bins = np.zeros(nbins)
    diff = r_max/nbins
    centers = np.linspace(0, r_max-diff, nbins)
    for i in range(len(x_pos)):
        x = x_pos[i]-offset[0]
        y = y_pos[i]-offset[1]
        z = z_pos[i]-offset[2]
        r = np.sqrt(x**2 + y**2 + z**2)
        r_bin = int(r*nbins/r_max)
        if r_bin < nbins:
            bins[r_bin] += 1
    return bins, centers

def cartesian_bins(x_pos, nbins=200, x_min=-2, x_max=2):
    bins = np.zeros(nbins)
    diff = (x_max-x_min)/nbins
    centers = np.linspace(x_min, x_max-diff, nbins)
    for i in range(len(x_pos)):
        x = x_pos[i]
        x_bin = int((x-x_min)*nbins/(x_max-x_min))
        if x_bin < nbins and x_bin >= 0:
            bins[x_bin] += 1
    return bins, centers


def analytic_radial(l, dim=3, r_max=1.2, N=500):
    r = np.linspace(0,r_max,N)
    pr = 1e-4*4*np.pi*(1e-6*r)**2/(np.pi**(dim/2.)*l**dim)*np.exp(-(r*1e-6)**2/l**2)
    return pr, r

def analytic_cartesian(l, dim=3, x_min=-2., x_max=2., N=500, offset=0.):
    x = np.linspace(x_min, x_max, N)
    px = 100./(np.pi**(dim/2.)*l**dim)*np.exp(-((x-offset))**2/l**2)
    return px, x

def read_file(f_n):
    file = open(f_n, 'r')

    x_pos = []
    y_pos = []
    z_pos = []

    for line in file.readlines():
        inline = line.split()
        x_pos.append(float(inline[2]))
        y_pos.append(float(inline[3]))
        z_pos.append(float(inline[4]))

    return x_pos, y_pos, z_pos


D = 4e-10
t = 1e-4
l = np.sqrt(4*D*t)*1e6
R = 8.31446
F = 96485.
T = 300.
psi = R*T/F
z = 2.
ex = 1e5
ey = 3e5
ez = -2e5
offset_x = 1e6*ex*z*D*t/psi
offset_y = 1e6*ey*z*D*t/psi
offset_z = 1e6*ez*z*D*t/psi



# pr, r = analytic_radial(l)
# plt.figure(figsize=(6,6))
# plt.plot(r, pr, label='Analytic')
#
# f_n = 'output_1/viz_data/seed_0002/Scene.ascii.1.dat'
# x_pos, y_pos, z_pos = read_file(f_n)
#
# bins, centers = radial_bins(x_pos, y_pos, z_pos, nbins=120)
# plt.plot(centers, bins, 'o', markerfacecolor='none', label=r'$\Delta t = $1 us, x100')
#
# f_n = 'output_2/viz_data/seed_0002/Scene.ascii.100.dat'
# x_pos, y_pos, z_pos = read_file(f_n)
#
# bins, centers = radial_bins(x_pos, y_pos, z_pos, nbins=120)
# plt.plot(centers, bins, 'x', label=r'$\Delta t = $100 us')
# plt.legend()


f_n = 'output_1/viz_data/seed_0002/Scene.ascii.1.dat'
x_pos_1, y_pos_1, z_pos_1 = read_file(f_n)
f_n = 'output_2/viz_data/seed_0002/Scene.ascii.100.dat'
x_pos_2, y_pos_2, z_pos_2 = read_file(f_n)

px, x = analytic_cartesian(l, offset=offset_x)
plt.figure()
plt.plot(x, px, label='Analytic')
bins, centers = cartesian_bins(x_pos_1, nbins=200)
plt.plot(centers, bins, 'o', markerfacecolor='none', label=r'$\Delta t = $1 \textmu s, x100', markersize=5)
bins, centers = cartesian_bins(x_pos_2, nbins=200)
plt.plot(centers, bins, 'x', label=r'$\Delta t = $100 \textmu s', markersize=5)
plt.plot([offset_x, offset_x],[-10, 350], '--', color='k', label=r'$zDtE_x/\psi$')
plt.ylabel(r"Number of molecules")
plt.xlabel(r"Distance (\textmu m)")
plt.legend()
plt.xlim([-2, 2])
plt.ylim([-5, 320])
plt.tight_layout()
plt.savefig('unbounded_x.pdf')

px, x = analytic_cartesian(l, offset=offset_y)
plt.figure()
plt.plot(x, px, label='Analytic')
bins, centers = cartesian_bins(y_pos_1, nbins=200)
plt.plot(centers, bins, 'o', markerfacecolor='none', label=r'$\Delta t = $1 \textmu s, x100', markersize=5)
bins, centers = cartesian_bins(y_pos_2, nbins=200)
plt.plot(centers, bins, 'x', label=r'$\Delta t = $100 \textmu s', markersize=5)
plt.plot([offset_y, offset_y],[-10, 350], '--', color='k', label=r'$zD\Delta tE_y/\psi$')
plt.ylabel(r"Number of molecules")
plt.xlabel(r"Distance (\textmu m)")
plt.legend()
plt.xlim([-2, 2])
plt.ylim([-5, 320])
plt.tight_layout()
plt.savefig('unbounded_y.pdf')

px, x = analytic_cartesian(l, offset=offset_z)
plt.figure()
plt.plot(x, px, label='Analytic')
bins, centers = cartesian_bins(z_pos_1, nbins=200)
plt.plot(centers, bins, 'o', markerfacecolor='none', label=r'$\Delta t = $1 \textmu s, x100', markersize=5)
bins, centers = cartesian_bins(z_pos_2, nbins=200)
plt.plot(centers, bins, 'x', label=r'$\Delta t = $100 \textmu s', markersize=5)
plt.plot([offset_z, offset_z],[-10, 350], '--', color='k', label=r'$zD\Delta tE_z/\psi$')
plt.ylabel(r"Number of molecules")
plt.xlabel(r"Distance (\textmu m)")
plt.legend()
plt.xlim([-2, 2])
plt.ylim([-5, 320])
plt.tight_layout()
plt.savefig('unbounded_z.pdf')

plt.show()
