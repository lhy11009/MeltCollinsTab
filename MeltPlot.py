from mpl_toolkits.basemap import Basemap
from scipy.interpolate import griddata
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import os
import math
import argparse


nox = 33
sclMF = 1e6 / (4 * math.pi * 1740**2.0 / (12 * (nox - 1)**2.0))
sclMFTotal = 1.0 / (4 * math.pi * 1740**2.0 / (12 * (nox - 1)**2.0))
jet_adjust_dict = {'red': ((0.0, 0, 0), (0.35, 0, 0), (0.66, 1, 1), (0.89, 1, 1),
                           (1, 0.5, 0.5)),
                   'green': ((0.0, 0, 0), (0.025, 0, 0), (0.275, 1, 1), (0.54, 1, 1),
                             (0.81, 0, 0), (1, 0, 0)),
                   'blue': ((0.0, 0.5, 0.5), (0.11, 1, 1), (0.18, 1, 1), (0.4, 0, 0),
                            (1, 0, 0))}
cName = ['All', 'Lower_Mantle', 'Upper_Mantle', 'IBC', 'Crust']
column = [0, 2, 3, 4, 5]


#   plot using python basemap
def BmPlot(caseName, step, _pTotal=False, **kwargs):
    route = os.path.join('..', caseName)
    oroute = '../eps'
    try:
        _chemical = kwargs['chemical']
    except KeyError:
        _chemical = 0
    _cName = cName[_chemical]
    _column = column[_chemical]
    mfroute = oroute
    outname = os.path.join(mfroute, "MF_%s_%d" % (_cName, step))
    otype = 'eps'
    if _pTotal is True:
        filename = os.path.join(route, "%s.MF_Total" % caseName)
        outname = os.path.join(mfroute, "MF_Total_%s" % _cName)
        sclField = sclMFTotal
        unit = "km"
    else:
        filename = os.path.join(route, "%s.MF0.%d" % (caseName, step))
        sclField = sclMF
        unit = "km/Myr"
    scl_ang = 180.0 / np.pi
    # read data
    try:
        mm = kwargs['mm']
    except KeyError:
        mm = np.zeros(2)
        with open(filename, 'r') as f:
            inputs = f.readline().split()
        mm[0] = float(inputs[0])
        mm[1] = float(inputs[1])
    mm[0] = mm[0] * sclField
    mm[1] = mm[1] * sclField
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
    fig = plt.figure(figsize=(8, 8))
    ax_map = fig.add_axes([0.1, 0.1, 0.8, 0.8])
    ax_cbr = fig.add_axes([0.91, 0.3, 0.01, 0.4])
    # set up map
    m = Basemap(projection='hammer', lon_0=0, resolution='c', ax=ax_map)
    # plot the model
    s = m.transform_scalar(field, lon[0, :], lat[:, 0], 1000, 500)
    cm = matplotlib.colors.LinearSegmentedColormap('my_colormap', jet_adjust_dict, 1024)
    im = m.imshow(s, cmap=cm, vmin=0.0, vmax=mm[1])
    m.drawparallels(np.arange(-90, 90, 30))
    m.drawmeridians(np.arange(-180, 180, 30))
    # color bar
    cb = plt.colorbar(im, cax=ax_cbr)
    cb.set_label('Melt Rate [%s]' % unit)
    fig.savefig("%s.%s" % (outname, otype))


def main():
    parser = argparse.ArgumentParser(description='Plot Model Melt')
    parser.add_argument('-n', '--name', type=str, default='MgH30v1e-2V5e20E100_hot',
                        help='Case name')
    parser.add_argument('-s', '--step', type=int, default=3600,
                        help='Step')
    parser.add_argument('-c', '--chemical', type=int, default=0,
                        help='Chemical, 0: total(default), 1-4: chemical0 - chemical3')
    parser.add_argument('-t', '--total', type=int, default=False,
                        help='Plot Total Value, 0: single step(default), 1: Total')
    arg = parser.parse_args()
    BmPlot(arg.name, arg.step, arg.total, chemical=arg.chemical)


main()
