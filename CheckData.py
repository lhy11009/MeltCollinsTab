import os
import numpy as np
import matplotlib.pyplot as plt
import argparse


class NoFileError(Exception):
    pass


def DataCheck(_filename, _col, **kwargs):
    try:
        _header = kwargs['header']
    except KeyError:
        _header = 0
    if os.path.isfile(_filename) is False:
        raise NoFileError(_filename)
    _data = np.genfromtxt(_filename, skip_header=_header)
    if len(_data.shape) == 1:
        _data.resize(_data.shape[0], 1)
    _min = np.min(_data[:, _col])
    _max = np.max(_data[:, _col])
    pts = np.array(range(_data.shape[0]))
    fig, ax = plt.subplots()
    ax.plot(pts, _data[:, _col], 'r.')
    ax.set(xlabel='N', ylabel='Data column %d' % _col,
           title='Data Range %.4e %.4e' % (_min, _max))
    fig.tight_layout()
    fig.savefig("../eps/dataCheck.png")


def SameValue(_d0, _d1, lim):
    if abs(_d0) < 1e-8:
        divi = 1.0
    else:
        divi = abs(_d0)
    if abs(_d1 - _d0) / divi < lim:
        return True
    else:
        return False


def OutofRange(_d0):
    lat0 = _d0['lat']
    lon0 = _d0['lon']
    bool0 = lat0 < -np.pi / 2.0 or lat0 > np.pi / 2.0
    bool1 = lon0 < -np.pi or lon0 > np.pi
    if bool0 or bool1:
        return True
    else:
        return False


def CheckDup(_filename, **kwargs):
    try:
        _header = kwargs['header']
    except KeyError:
        _header = 0
    if os.access(_filename, os.R_OK) is False:
        raise NoFileError(_filename)
    _data = np.genfromtxt(_filename, skip_header=_header)
    dtype = [('index', int), ('lat', float), ('lon', float)]  # structed
    group = []
    for i in range(_data.shape[0]):
        this = (i, _data[i, 0], _data[i, 1])
        group.append(this)
    structData = np.array(group, dtype=dtype)
    print("CheckDup: Sort input")
    structData = np.sort(structData, order='lat', kind='mergesort')  # sort
    ofile = os.path.join("..", "duplicate_points")  # compare & output
    ofile1 = os.path.join("..", "out_range_points")
    print("CheckDup: Check dupicate and out of range")
    fo = open(ofile, "w")
    fo1 = open(ofile1, "w")
    for i in range(structData.size - 1):  # duplicate
        d0 = structData[i]
        lat0 = d0['lat']
        lon0 = d0['lon']
        for j in range(i + 1, structData.size):
            dThis = structData[j]
            latThis = dThis['lat']
            lonThis = dThis['lon']
            if SameValue(lat0, latThis, 1e-5) is False:
                break
            elif SameValue(lon0, lonThis, 1e-5) is True:
                fo.write("%d %d\n" % (d0['index'] + 2, dThis['index'] + 2))
    for i in range(structData.size - 1):  # out of range
        d0 = structData[i]
        if OutofRange(d0):
            fo1.write("%d\n" % (d0['index'] + 2))
    fo.close()
    fo1.close()


def main():
    parser = argparse.ArgumentParser(description='Do Data Check')
    parser.add_argument('-f', '--File', type=str,
                        default='MgH30v1e-2V5e20E100_hot/MgH30v1e-2V5e20E100_hot.MF_Total',
                        help="File name")
    parser.add_argument('-c', '--Column', type=int, default=0,
                        help="Column")
    parser.add_argument('-r', '--RowHeader', type=int, default=0, help="Header")
    parser.add_argument('-d', '--Duplicate', type=int, default=0, help="Whether check duplicate")
    arg = parser.parse_args()
    if arg.Duplicate == 0:
        DataCheck(arg.File, arg.Column, header=arg.RowHeader)
    else:
        CheckDup(arg.File, header=arg.RowHeader)


main()
