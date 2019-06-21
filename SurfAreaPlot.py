from mpl_toolkits.basemap import Basemap
from scipy.interpolate import griddata
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import os
from FGenLib import SurfArea
from ReadFiles import read_input


sclS = 1740**2.0
jet_adjust_dict = {'red': ((0.0, 0, 0), (0.35, 0, 0), (0.66, 1, 1), (0.89, 1, 1),
                           (1, 0.5, 0.5)),
                   'green': ((0.0, 0, 0), (0.025, 0, 0), (0.275, 1, 1), (0.54, 1, 1),
                             (0.81, 0, 0), (1, 0, 0)),
                   'blue': ((0.0, 0.5, 0.5), (0.11, 1, 1), (0.18, 1, 1), (0.4, 0, 0),
                            (1, 0, 0))}


class NoFileError(Exception):
    pass


#   plot using python basemap
def BmPlot(fileCoord, _sA, **kwargs):
    oroute = '../eps'
    try:
        prj = kwargs['projection']
    except KeyError:
        prj = 'hammer'
    mfroute = oroute
    otype = 'jpg'
    outname = os.path.join(mfroute, "SurfArea_%s" % (prj))
    sclField = sclS
    unit = "km^2.0"
    scl_ang = 180.0 / np.pi
    # read data
    data = np.genfromtxt(fileCoord, skip_header=1)
    ilon = data[:, 1] * scl_ang
    ilat = data[:, 0] * scl_ang
    ifield = _sA * sclField
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


def ReadParaFromFile():
    if os.path.isfile('inputFileSurf') is False:
        raise NoFileError('inputFile')
    _inD = read_input('inputFileSurf')
    return _inD


def main():
    inD = ReadParaFromFile()
    fileCoord = inD['fileCoord']
    routeA = inD['routeA']
    prefixA = inD['prefixA']
    SArea = SurfArea(inD)
    sA = SArea(routeA, prefixA, col=1)     # read surf_area file
    BmPlot(fileCoord, sA)


main()
