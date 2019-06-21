import os
import numpy as np
from ReadFiles import read_input
from FGenLib import line_prepender
import argparse


class NoFileError(Exception):
    pass


class Cls():
    def __init__(self, _route, _caseName, _timeScale):
        self.route = _route
        self.caseName = _caseName
        self.timeScale = _timeScale


def ReadParaFromFile(inputFile):
    if os.path.isfile(inputFile) is False:
        raise NoFileError(inputFile)
    _inD = read_input(inputFile)
    return _inD


def SurfArea(inputFile, col):
    inD = ReadParaFromFile(inputFile)   # control parameters
    route = inD['route']
    oroute = inD['oroute']
    prefix = inD['prefix']
    oroute = os.path.join(oroute)
    if os.path.isdir(oroute) is False:
        os.mkdir(oroute)
    nprocx = int(inD['nprocx'])
    nprocz = int(inD['nprocz'])
    nproc = 12 * nprocx**2 * nprocz
    ifile = os.path.join(route, "%s.surf_area.%d" % (prefix, nprocz - 1))  # read in data
    print("SurfArea: %s" % ifile)
    if os.access(ifile, os.R_OK) is False:
        raise NoFileError(ifile)
    data = np.genfromtxt(ifile)
    for i in range(2 * nprocz - 1, nproc, nprocz):
        ifile = os.path.join(route, "%s.surf_area.%d" % (prefix, i))
        print("SurfArea: %s" % ifile)
        if os.access(ifile, os.R_OK) is False:
            raise NoFileError(ifile)
        data1 = np.genfromtxt(ifile)
        data = np.concatenate((data, data1))
    if len(data.shape) == 1:
        data.resize(data.shape[0], 1)   # fix 1-column bug
    ofile = os.path.join(oroute, "element_area.dat")
    fo = open(ofile, 'w')  # combine and output
    for j in range(data.shape[0]):
        fo.write("%.4e\n" % data[j, col])
    fo.close()
    sumV = data[:, col].sum()  # header, max and sum
    maxV = np.max(data[:, col])
    line = "%.4e %.4e\n" % (maxV, sumV)
    line_prepender(ofile, line)


def main():
    parser = argparse.ArgumentParser(description='Combine surf_area files')
    parser.add_argument('-f', '--File', type=str, default='inputFile1',
                        help="File name")
    parser.add_argument('-c', '--Col', type=int, default=0,
                        help="Colume of interest")
    arg = parser.parse_args()
    SurfArea(arg.File, arg.Col)


main()
