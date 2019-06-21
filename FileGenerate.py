import os
from FGenLib import CombineSurfFile
from FGenLib import SurfArea
from ReadFiles import read_input
from ReadFiles import get_steps
from ReadFiles import get_variable
from ReadFiles import read_time
import multiprocessing
from joblib import Parallel, delayed
import numpy as np

Ma = 1e6 * 365 * 24 * 3600
Ro = 1740
S = Ro**2.0


class NoFileError():
    pass


class Cls():
    def __init__(self, _route, _caseName, _timeScale):
        self.route = _route
        self.caseName = _caseName
        self.timeScale = _timeScale


def ReadParaFromFile():
    if os.path.isfile('inputFile') is False:
        raise NoFileError('inputFile')
    _inD = read_input('inputFile')
    return _inD


def StepTime(_oroute, _caseName, _stepTuple, _timeArray, _timeScale):
    oFile = os.path.join(_oroute, "StepTime")
    fo = open(oFile, 'w')
    for i in range(len(_stepTuple)):
        _step = _stepTuple[i]
        _time = _timeArray[_step] * _timeScale
        fo.write("%d %.4e\n" % (_step, _time))


def main():
    inD = ReadParaFromFile()
    route = inD['route']
    oroute = inD['oroute']
    prefix = inD['prefix']
    routeA = inD['routeA']
    prefixA = inD['prefixA']
    if os.path.isdir(oroute) is False:
        os.mkdir(oroute)
    nprocz = inD['nprocz']
    #   Case steps
    pp_file = os.path.join(route, "post_process")
    p_dict = read_input(pp_file)
    e_steps = get_variable(p_dict, 'episode_steps', 'int_list')
    step_tuple = get_steps(route, prefix, 100000)
    stepInEpisode = [step for step in step_tuple if step <= e_steps[-1] and step > 0]
    #   Combine MF files
    SArea = SurfArea(inD)
    sA = SArea(routeA, prefixA, col=1)     # read surf_area file
    sA *= S
    combineS = CombineSurfFile('MF', inD, surfArea=sA)
    # combineS(route, prefix, 1000, nprocz, "surf",   # test a single step
    #         obl=False, expandLon=True, oroute=oroute, area=1.0,
    #         eaMesh=True)
    num_cores = multiprocessing.cpu_count()
    print("start parallel: cores %d" % num_cores)
    Parallel(n_jobs=num_cores)(delayed(combineS)(route, prefix, step, nprocz, "surf",
                                                 obl=True, expandLon=True, oroute=oroute,
                                                 area=1.0, eaMesh=False)
                               for step in stepInEpisode)
    #   Total amount of melt
    pFile = os.path.join(route, "input_solidus_i", "in_65moon")
    pDict = read_input(pFile)
    Ro = get_variable(pDict, 'radius', 'float')
    Kappa = get_variable(pDict, 'thermdiff', 'float')
    timeScale = pow(Ro, 2.0) / Kappa
    _cls = Cls(oroute, prefix, timeScale)    # input class
    timeArray = read_time(route, prefix)    # time
    combineS.TotalAmount(_cls, stepInEpisode, timeArray, obl=True, area=1.0)
        # im = m.imshow(s, cmap=cm, vmin=0.0, vmax=mm[1])
    # combineS.EaGrid(oroute, prefix, nprocz, total=True)
    #   Time of Steps
    StepTime(oroute, prefix, step_tuple, timeArray, timeScale/Ma)


main()
