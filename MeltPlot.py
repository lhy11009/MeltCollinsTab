from mpl_toolkits.basemap import Basemap
from scipy.interpolate import griddata
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import os
import math
import argparse


nox = 64
sclMF_old = 1e6 / (4 * math.pi * 1740**2.0 / (12 * (nox - 1)**2.0))
sclMF = 1e6
sclMFTotal_old = 1.0 / (4 * math.pi * 1740**2.0 / (12 * (nox - 1)**2.0))
sclMFTotal = 1.0
jet_adjust_dict = {'red': ((0.0, 0, 0), (0.35, 0, 0), (0.66, 1, 1), (0.89, 1, 1),
                           (1, 0.5, 0.5)),
                   'green': ((0.0, 0, 0), (0.025, 0, 0), (0.275, 1, 1), (0.54, 1, 1),
                             (0.81, 0, 0), (1, 0, 0)),
                   'blue': ((0.0, 0.5, 0.5), (0.11, 1, 1), (0.18, 1, 1), (0.4, 0, 0),
                            (1, 0, 0))}
cName = ['All', 'Lower_Mantle', 'Upper_Mantle', 'IBC', 'Crust']
column = [0, 2, 3, 4, 5]


#   plot using python basemap
def BmPlot(filename, step, _pTotal=0, **kwargs):
    oroute = '../eps'
    try:
        _chemical = kwargs['chemical']
    except KeyError:
        _chemical = 0
    try:
        prj = kwargs['projection']
    except KeyError:
        prj = 'hammer'
    _cName = cName[_chemical]
    _column = column[_chemical]
    mfroute = oroute
    otype = 'eps'
    otype1 = 'jpg'
    if _pTotal == 1:
        outname = os.path.join(mfroute, "MF_%s_Total_%s" % (prj, _cName))
        sclField = sclMFTotal
        unit = "km"
    else:
        outname = os.path.join(mfroute, "MF_%s_%s_%d" % (prj, _cName, step))
        sclField = sclMF
        unit = "km/Myr"
    scl_ang = 180.0 / np.pi
    # read data
    data = np.genfromtxt(filename, skip_header=1)
    ilon = data[:, 1] * scl_ang
    ilat = data[:, 0] * scl_ang
    ifield = data[:, 2 + _column] * sclField
    # mesh data
    res = 1.0
    lon = np.arange(-180, 180 + res, 1.0)
    lat = np.arange(-90.0, 90.0 + res, 1.0)
    lon, lat = np.meshgrid(lon, lat)
    field = griddata((ilon, ilat), ifield, (lon, lat), method='linear')
    # initialize figure and axes
    cm = matplotlib.colors.LinearSegmentedColormap('my_colormap', jet_adjust_dict, 1024)
    if prj == 'hammer':
        fig = plt.figure(figsize=(8, 8))
        ax_map = fig.add_axes([0.1, 0.1, 0.8, 0.8])
        ax_cbr = fig.add_axes([0.91, 0.3, 0.01, 0.4])
        # set up map
        m = Basemap(projection='hammer', lon_0=0, resolution='c', ax=ax_map)
        # plot the model
        s = m.transform_scalar(field, lon[0, :], lat[:, 0], 1000, 500)
        im = m.imshow(s, cmap=cm)
        m.drawparallels(np.arange(-90, 90, 30))
        m.drawmeridians(np.arange(-180, 180, 30))
        cb = plt.colorbar(im, cax=ax_cbr)
    elif prj == 'cart':
        fig, ax = plt.subplots()
        im = ax.pcolormesh(lon, lat, field, cmap=cm)
        cb = plt.colorbar(im, ax=ax)
    else:
        print("Wrong type of projection: %s" % prj)
        exit()
    # color bar
    cb.set_label('Melt Rate [%s]' % unit)
    fig.savefig("%s.%s" % (outname, otype))
    fig.savefig("%s.%s" % (outname, otype1))


def main():
    parser = argparse.ArgumentParser(description='Plot Model Melt')
    parser.add_argument('-f', '--fname', type=str,
                        default='../MgH30v1e-2V5e20E100_hot/MgH30v1e-2V5e20E100_hot.MF0.1000',
                        help='File name')
    parser.add_argument('-s', '--step', type=int, default=1000,
                        help='Step')
    parser.add_argument('-c', '--chemical', type=int, default=0,
                        help='Chemical, 0: total(default), 1-4: chemical0 - chemical3')
    parser.add_argument('-t', '--total', type=int, default=0,
                        help='Plot Total Value, 0: single step(default), 1: Total')
    parser.add_argument('-p', '--projection', type=str, default='hammer',
                        help='Projecton type, hammer or cart')
    arg = parser.parse_args()
    BmPlot(arg.fname, arg.step, arg.total, chemical=arg.chemical,
           projection=arg.projection)


main()
