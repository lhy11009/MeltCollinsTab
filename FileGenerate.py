import os
from FGenLib import CombineSurfFile
from ReadFiles import read_input
from ReadFiles import get_steps
from ReadFiles import get_variable
from ReadFiles import read_time
import multiprocessing
from joblib import Parallel, delayed
import numpy as np

Ma = 1e6 * 365 * 24 * 3600


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
    _route = inD['route']
    _oroute = inD['oroute']
    _caseName = inD['caseName']
    _oroute = os.path.join(_oroute, _caseName)
    if os.path.isdir(_oroute) is False:
        os.mkdir(_oroute)
    nprocz = inD['nprocz']
    #   Case steps
    pp_file = os.path.join(_route, "post_process")
    p_dict = read_input(pp_file)
    e_steps = get_variable(p_dict, 'episode_steps', 'int_list')
    step_tuple = get_steps(_route, _caseName, 100000)
    stepInEpisode = [step for step in step_tuple if step <= e_steps[-1] and step > 0]
    #   Combine MF files
    combineS = CombineSurfFile('MF', inD)
    combineS(_route, _caseName, 1000, nprocz, "surf",
             obl=True, expandLon=False, oroute=_oroute, area=1.0)
    exit()
    num_cores = multiprocessing.cpu_count()
    print("start parallel: cores %d" % num_cores)
    Parallel(n_jobs=num_cores)(delayed(combineS)(_route, _caseName, step, nprocz, "surf",
                                                 obl=True, expandLon=False, oroute=_oroute,
                                                 area=1.0)
                               for step in stepInEpisode[1: len(stepInEpisode)+1])
    #   Total amount of melt
    pFile = os.path.join(_route, "input_solidus_i", "in_65moon")
    pDict = read_input(pFile)
    Ro = get_variable(pDict, 'radius', 'float')
    Kappa = get_variable(pDict, 'thermdiff', 'float')
    timeScale = pow(Ro, 2.0) / Kappa
    _cls = Cls(_oroute, _caseName, timeScale)    # input class
    timeArray = read_time(_route, _caseName)    # time
    combineS.TotalAmount(_cls, stepInEpisode, timeArray, obl=True, area=1.0)
    #   Time of Steps
    StepTime(_oroute, _caseName, step_tuple, timeArray, timeScale/Ma)


main()
