import os
import numpy as np
import matplotlib.pyplot as plt
import argparse


class NoFileError:
    pass


def DataCheck(_filename, _col, **kwargs):
    try:
        _header = kwargs['header']
    except KeyError:
        _header = 0
    if os.path.isfile(_filename) is False:
        raise NoFileError(_filename)
    _data = np.genfromtxt(_filename, skip_header=_header)
    _min = np.min(_data[:, _col])
    _max = np.max(_data[:, _col])
    pts = np.array(range(_data.shape[0]))
    fig, ax = plt.subplots()
    ax.plot(pts, _data[:, _col], 'r.')
    ax.set(xlabel='N', ylabel='Data column %d' % _col,
           title='Data Range %.4e %.4e' % (_min, _max))
    fig.tight_layout()
    fig.savefig("./eps/dataCheck.png")


def main():
    parser = argparse.ArgumentParser(description='Do Data Check')
    parser.add_argument('-f', '--File', type=str,
                        default='MgH30v1e-2V5e20E100_hot/MgH30v1e-2V5e20E100_hot.MF_Total',
                        help="File name")
    parser.add_argument('-c', '--Column', type=int, default=0,
                        help="Column")
    parser.add_argument('-r', '--RowHeader', type=int, default=0, help="Header")
    arg = parser.parse_args()
    DataCheck(arg.File, arg.Column, header=arg.RowHeader)


main()
