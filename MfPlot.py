import matplotlib.pyplot as plt
import numpy as np
import os
import argparse

cDict = {'step': 0, 'stime': 1, 'time': 2, 'melt': [0, 3, 4, 5, 6]}
color = ['k', 'r', 'b', 'g', 'c']
cname = ['Total', 'Lower Mantle', 'Upper Mantle', 'IBC', 'Crust']


def MfPlot(filename):
    oroute = '../eps'
    otype = '.eps'
    otype1 = '.eps'
    data = np.genfromtxt(filename)      # read in
    shape = data.shape
    stC = cDict['stime']
    stime = data[:, stC]
    tC = cDict['time']
    time = data[:, tC]
    melt = data[:, tC + 1: shape[1] + 1]
    mA = cDict['melt']
    wmelt = np.zeros(melt.shape)
    for i in range(1, shape[0]):
        wmelt[i, :] = wmelt[i - 1, :] + melt[i - 1, :] * stime[i - 1]     # total melt
    fig, ax = plt.subplots()    # plot
    ax = plt.subplot(2, 1, 1)
    for i in range(len(mA)):
        iC = mA[i]
        ax.plot(time / 1e6, melt[:, iC], '-'+color[i], label=cname[i])
    ax.set(xlabel='Time [Ma]', ylabel='Melting Rate [km/Ma]')
    ax = plt.subplot(2, 1, 2)
    for i in range(len(mA)):
        iC = mA[i]
        ax.plot(time / 1e6, wmelt[:, iC], '-'+color[i], label=cname[i])
    ax.set(xlabel='Time [Ma]', ylabel='Melt [km^3]')
    ax.legend()
    fig.tight_layout()
    ofile = os.path.join(oroute, 'mf')
    fig.savefig(ofile+otype)
    fig.savefig(ofile+otype1)


def main():
    parser = argparse.ArgumentParser(description='Plot mf file')
    parser.add_argument('-f', '--filename', type=str,
                        default='../MgH30v1e-2V5e20E100_hot/MgH30v1e-2V5e20E100_hot.mf.dat',
                        help='Filename')
    arg = parser.parse_args()
    MfPlot(arg.filename)


main()
