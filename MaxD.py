import numpy as np
import os
import argparse
from FGenLib import LaterBetter
# from FGenLib import SingularPoints


stepCol = 0


class NoFileError(Exception):
    pass


def MaxV(filename, col, header):
    lc = stepCol
    if os.access(filename, os.R_OK) is False:
        raise NoFileError(filename)
    data = np.genfromtxt(filename, skip_header=header)
    if lc >= 0:
        data = LaterBetter(data, lc)    # adjust restart
    # data = SingularPoints(data, col)  # delete singular points
    mN = np.argmax(data[:, col])    # max value
    mV = data[mN, col]
    print("Index: %.4e, Max Value: %.4e" % (mN, mV))
    np.savetxt(filename, data)    # output


def main():
    parser = argparse.ArgumentParser(description="Max value along first axits")
    parser.add_argument('-f', "--filename", type=str,
                        default="MgH30v1e-2V5e20E100_hot/MgH30v1e-2V5e20E100_hot.mf.dat",
                        help="File name")
    parser.add_argument('-c', '--column', type=int, default=0, help='Column')
    parser.add_argument('-r', '--row_header', type=int, default=0, help='Skip row header')
    parser.add_argument('-l', '--later_better', type=int, default=-1,
                        help='Column for LaterBetter')
    arg = parser.parse_args()
    MaxV(arg.filename, arg.column, arg.row_header)


main()
