import os
from FGenLib import PlotSurf
from FGenLib import ReadParaFromFile
from FGenLib import MakeAni
from ReadFiles import get_steps
from ReadFiles import read_input
from ReadFiles import get_variable
from ReadFiles import read_time
import argparse
import multiprocessing
from joblib import Parallel, delayed


def MakeGif(filename, chemical, skipFig):
    inD = ReadParaFromFile(filename)    # init
    droute = inD['droute']
    prefix = inD['prefix']
    citcomsFile = inD['citcomsFile']
    pp_file = os.path.join(droute, "post_process")
    p_dict = read_input(pp_file)
    e_steps = get_variable(p_dict, 'episode_steps', 'int_list')
    stepTuple = get_steps(droute, prefix, 100000)
    stepTuple = [step for step in stepTuple if step <= e_steps[-1] and step > 0]
    PlotSurf0 = PlotSurf(inD, 'MF0', stepTuple)
    # plot figures
    if skipFig == 0:
        # maxV = PlotSurf0.GetUnvMax()     # get max value
        maxV = PlotSurf0.GetUnvMax1(1, header=1)     # get max value
        num_cores = multiprocessing.cpu_count()
        print("start parallel: cores %d" % num_cores)
        Parallel(n_jobs=num_cores)(delayed(PlotSurf0)(step, chemical=chemical, maxV=maxV, otype='jpg')
                                   for step in stepTuple)
    # make gif
    timeTuple = read_time(droute, prefix)
    citcomsD = read_input(citcomsFile)
    Ro = get_variable(citcomsD, 'radius', 'float')
    Kappa = get_variable(citcomsD, 'thermdiff', 'float')
    timeScale = pow(Ro, 2.0) / Kappa
    MakeAni0 = MakeAni(inD, timeScale, stepTuple, timeTuple)
    MakeAni0('MF_hammer', chemical=chemical, removeJpg=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--File', type=str, help='Input File', default='inputFileAni')
    parser.add_argument('-c', '--chemical', type=int, help='Chemical', default=0)
    parser.add_argument('-s', '--skipFig', type=int, help='Skip Figure', default=0)
    arg = parser.parse_args()
    MakeGif(arg.File, arg.chemical, arg.skipFig)


main()
