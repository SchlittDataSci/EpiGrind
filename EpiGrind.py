import numpy as np
import random as ra
import pandas as pd
import urllib
import zipfile
import sys
import os
import mpld3
import math
#from numba import jit

from copy import deepcopy
from itertools import chain
from multiprocessing import Pool, cpu_count
from geopandas.plotting import plot_dataframe
from os import listdir
from shapely.geometry import Point, shape
from shapely import speedups
from glob import glob
#from functools import partial

import matplotlib.pyplot as plt
import geopandas as gp

ra.seed(8675309)
crsESPG = 4326
badRun = 999999999999

speedups.enable()



### Import model if not ready

def getModel(ref):
    """Imports model from str if given str"""
    if type(ref) is str:
        model = __import__(ref)
        return model
    else:
        return ref

def runScriptedSim(script,modelRef='EGManoreEtAlWithTravel',duration='null'):
    """Runs a sim from str arguments"""
    scenario='null'
    compartments='null'
    exposures='null'
    
    model = getModel(modelRef)
    exec(script)
    
    try:
        if duration != 'null':
            scenario['days'] = duration
    except:
        pass
        
    return runSim(scenario,compartments,exposures,modelRef=modelRef)


### Model Execution

def runSim(scen,
           compartmentsIn,
           exposuresIn,
           modelRef='EGDefaultSEIR'):
    """Runs a sim"""
    model = getModel(modelRef)
    compartments = deepcopy(compartmentsIn)
    exposures = deepcopy(exposuresIn)
    
    try:
        shifts = scen['shifts']
    except:
        shifts = {}
    
    if type(shifts) is dict:
        useShifts = True

    compartments,scen,N = model.modelPrep(compartments,scen)

    epidemic = []
    
    for t in range(scen['days']):
        compartments,scen = model.modelExec(compartments,
                                            exposures,
                                            scen,
                                            shifts,
                                            useShifts,
                                            t)
        
        epidemic.append(deepcopy(compartments))
    
    return epidemic,N


    
### Base denaturing algorithm

def ezDenat(baseSim,template,expVars,
            valueFunc = "null",
            scoreFunc = "null",
            stagnant=500,
            verbose=True,
            recordEvery = 1,
            initMute = 3,
            select='null',
            shapes='null',
            showVars=False,
            externalCalib='null',
            calibCats='null',
            modelRef='null',
            duration='null',
            externalDetection=1.0,
            figureArgs={}):
    """Flexible stochastic search algorithim"""
    
    scores = []
    solutions = []
    plots = []

    bestScore = badRun
    gen = 0
    sinceImproved = 0
    improvements = 0
    bestMatch = deepcopy(expVars)
    testVars = deepcopy(expVars)
    surgeTime = stagnant
    quietTime = stagnant*5
    
    calibDetection = True
    
    if duration == 'null':
        duration = int(baseSim.split("'days':")[1].split(',')[0])
    if type(externalDetection) is tuple:
        calibDetection = True
        testDetection = deepcopy(externalDetection)
        externalDetectionOut = testDetection[0]
    else:
        calibDetection = False
        externalDetectionOut = externalDetection 
    
    targetCurves,N = runScriptedSim(baseSim.replace('$DETECTION',str(externalDetectionOut)),
                                    modelRef=modelRef,
                                    duration=duration)
                                    
    targetScores,targetPlots = extractOutputs(targetCurves,
                                              fx=valueFunc,
                                              n=N,
                                              calibCats=calibCats)
    
    if externalCalib != 'null':
        extCurves = pd.DataFrame.from_csv(externalCalib)
        if extCurves.columns.tolist() != 'all':
            extCurves = flattenCaseByGeo(extCurves,duration)
        targetScores = extCurves['all']
        targetScores = extCurves['all']*(1./externalDetectionOut)
    
    refineStage = False
    
    #if "#ENABLE CHLOROPLETH" in baseSim:
    #    plotCpleth(targetCurves,shapes,var='EPITOTAL',size=7,scheme='QUANTILES',title='Target results')
    
    while sinceImproved < stagnant or not refineStage:
        expCurves,N = runScriptedSim(prepScript(template,testVars).replace('$DETECTION',str(externalDetectionOut)),
                                                    modelRef=modelRef,
                                                    duration=duration)
                                                    
        expScores, expPlots = extractOutputs(expCurves,
                                                 fx=valueFunc,
                                                 n=N,
                                                 calibCats=calibCats)

        score = computeObjective(expScores,targetScores,fx=scoreFunc)
        
        if gen == surgeTime:
            if verbose:
                print "Ending surge mutation rate"
            initMute = 1
            sinceImproved = 0
        elif gen == quietTime:
            recordEvery = recordEvery*10
        
        if score < bestScore:
            scores.append([score,gen])
            bestMatch = deepcopy(testVars)
            sinceImproved = 0
            improvements += 1
            bestScore = score
            if calibDetection:
                externalDetection = testDetection
            
            if improvements % recordEvery == 0:
                solutions.append(deepcopy(expVars))
                plots.append(expScores)
                if verbose:
                    print "Gen: %s    Best Score: %s" % (gen,score)
                    if showVars:
                        print bestMatch
            
        gen += 1
        sinceImproved += 1
            
        if sinceImproved >= stagnant and not refineStage:
            if gen>surgeTime:
                if verbose:
                    print "Switching to refinement mode at gen: %s" % gen
                sinceImproved = 0
                refineStage = True
                select = '1'
                divisor = {True:5,False:1}
                initMute = 0.2
        
        testVars = shiftVars(bestMatch,select=select,mute=initMute)
        if calibDetection:
            testDetection = shiftVars({'temp':externalDetection},select='1',mute=initMute)['temp']
            targetScores = extCurves['all']*(1./testDetection[0])
    
    solutions.append(deepcopy(expVars))
    plots.append(expScores)
                                 
    print "Finishing sim at gen %s with best score of %s" % (gen,bestScore)
    if verbose:
        print "Best variable match found: %s" % bestMatch
    
    if "#ENABLE CHLOROPLETH" in baseSim:
        expCurves,N = runScriptedSim(prepScript(template,bestMatch),duration=duration)
        plotCpleth(expCurves,shapes,var='EPITOTAL',size=7,scheme='QUANTILES',title='Experimental results')
        
    plotExp(targetScores,plots,figureArgs=figureArgs)
    
    if calibDetection:
        print "Best detection rate found: %s" % externalDetection[0]
    
    return bestMatch,scores,solutions,plots,targetScores 
    

            
def tryAFewTest(baseScript,flexScript,expVars,
            tries=3,
            threads = 1,
            stagnant=50,
            verbose=False,
            **kwargs):
            
    """Runs multiple iteratinos of denaturing algorithm and returns the best"""
    
    p = Pool(threads)
    sols = []
    kwargs['verbose'] = False
    kwargs['pool'] = {}
    for i in range(tries):
        kwargs['iter'] = i
        p.apply_async(ezDenat,(baseScript,flexScript,expVars,),kwargs)
    p.close()
    p.join()
    #results = [sol.get() for sol in sols]
    return sols,False
        
    #data = triedToDf(tried,scores)
    #return data


def tryAFew(baseScript,flexScript,expVars,tries=3,verbose=True,**kwargs):
    """Runs multiple iteratinos of denaturing algorithm and returns the best"""
    tried = []
    scores = []
    for i in range(tries):
        if verbose:
            print "Attempt", i,
        bestN,scoresN,temp,plotsN,targetN = ezDenat(baseScript,flexScript,expVars,verbose=False,**kwargs)
        tried.append(bestN)
        scores.append(scoresN[-1])
        
    data = triedToDf(tried,scores)

    return data
    
    
    
def triedToDf(tried,scores):
    """Generates dataframe of calibration run final values"""
    keys = tried[0].keys()
    data = pd.DataFrame(columns=keys)
    for attempt in tried:
        dictOut = {key:value[0] for key,value in attempt.iteritems()}
        data = data.append(dictOut,ignore_index=True)
    data['score'] = [i[0]for i in scores]
    return data
    

      
    
def fillValsFromRow(expVars,data,row):
    varsOut = deepcopy(expVars)
    for key in expVars.keys():
        varsOut[key] = tuple([float(i) for i in[float(data[key][row])]+list(varsOut[key])[1:]])
    return varsOut
    
    
    
def sweepsToExpVars(data,priorVars,index='null'):
    """Generates set of experimental variables from tryAFew() result DF"""
    tRange = 2
    expVars =  dict()
    summaryStats = data.describe()
    if index == 'null':
        index = data['score'].idxmin()
    if priorVars != 'null':
        for key in list(data.columns):
            bestVar = data[key][index]
            if key != 'score':
                expVars[key] = (bestVar,
                           summaryStats[key]['std']/bestVar,
                           max(bestVar-summaryStats[key]['std']*2,0.,priorVars[key][2]),
                           min(bestVar+summaryStats[key]['std']*2,priorVars[key][3]))
    else:
        for key in list(data.columns):
            bestVar = data[key][index]
            if key != 'score':
                expVars[key] = (bestVar,0,0,0)
        
    try:
        del expVars['score']
    except:
        pass
    return expVars



def autoPilot(readyToRun,flexScript,expVars,
              tries=50,
              stagnant=200,
              layers=3,
              outFile='ODEAuto',
              **kwargs):
    if not os.path.exists('Output/'+outFile):
        os.makedirs('Output/'+outFile)
    for i in range(layers):
        print "\nStarting Layer %s" % (i+1)
        replicates = tryAFew(readyToRun,flexScript,expVars,
                               tries=tries,
                               stagnant=stagnant*(i+1),
                               **kwargs)
        replicates.to_csv('output/%s/%s%s.csv' % (outFile,outFile,i))
        expVars = sweepsToExpVars(replicates,expVars)
        plt.title("%s Layer %s Calibration Sweeps" % (outFile.title(),i+1))
        plt.show()
        plt.gca()
    print """AutoPilot Complete!"""
    return replicates
    
    
def multiWeekSweep(readyToRun,flexScript,sweeps,iterator,
                   outFile='multiWeekCalib',
                   **kwargs):
    for i in iterator:
        autoPilot(readyToRun,flexScript,sweeps,
             outFile=outFile+str(i),
             duration=i,
             **kwargs)


### Modular output and objective function derivation

def extractOutputs(curves,n=1.,fx="null",calibCats='null'):
    """Converts pop daily compartment percentages for plots and scoring"""
    catVars = {key:[] for key in curves[0].keys()}

    if fx == 'null':
        runFx = lambda x,n: x['newInf']*n
    else:
        runFx = eval('lambda x,n: '+fx)

    for epiDay in curves:
        for item in epiDay.items():
            result = runFx(item[1],n)
            catVars[item[0]].append(result)
    if calibCats == 'null':
        flatVars = list(chain(*[i[1] for i in sorted(catVars.items())]))
    else:
        flatVars = list(chain(*[i[1] for i in sorted(catVars.items()) if i[0] in calibCats]))
    return flatVars,catVars


def computeObjective(experimental,target,fx='null'):
    """Provides a score for denaturing algorithm"""
    
    if fx == 'null':
        runFx = lambda a,b: ((a-b)**2)/a
    else:
        runFx = eval('lambda a,b: '+fx)
    
    if type(target) is list:
        zipped = zip(target,experimental)
        results = [runFx(a,b) for (a,b) in zipped if a != 0]
        return np.mean(results)
    elif type(target) is pd.core.series.Series:
        target = target.to_dict()
    results = [runFx(a,experimental[i]) for (i,a) in target.items() if a != 0 and i < len(experimental)]
    return np.mean(results)
    

def flattenCaseByGeo(fileRef,duration):
    """Pulls a dataframe of cumulative cases by patch"""
    if type(fileRef) is str:
        caseData = pd.DataFrame.from_csv(glob(fileRef)[0])
    else:
        caseData = fileRef
    listed = [np.nan]*duration*len(caseData.columns)
    pos = 0
    for i in sorted(caseData.columns):
        for j in caseData.index:
            if not math.isnan(caseData[i][j]):
                listed[int(j)+pos] = caseData[i][j]
        pos += duration

    flattened = pd.DataFrame(columns=['all'])
    for i,j in enumerate(listed):
        if not math.isnan(j):
            flattened.set_value(i,'all',int(j))
    
    return flattened



def computeObjectiveOld(experimental,target,fx='null'):
    """Provides a score for denaturing algorithm"""
    zipped = zip(target,experimental)
    
    if fx == 'null':
        runFx = lambda a,b: ((a-b)**2)/a
    else:
        runFx = eval('lambda a,b: '+fx)

    try:
        results = [runFx(a,b) for (a,b) in zipped if a != 0]
        return np.mean(results)
    except:
        return badRun
    
    
    
### Base SEIR Model Functions

def setEpiCounts(epiRecord,N,stocks={'s','e','i','r'}):
    """Weights epi stocks by pop size"""
    for i,day in enumerate(epiRecord):
        epiRecord[i] = {item[0]:{stock:item[1][stock]*N for stock in stocks} for item in day.items()}
    return epiRecord


### Visualization Functions

def animateCpleth(curves,shapes,skip=10,title='',fx='null'):
    """Generates chloropleths by increments for patch model"""
    dirOut = 'plots/%s/' % title.replace(' ','')
    if not os.path.exists(dirOut):
        os.makedirs(dirOut)
    for t,day in enumerate(curves[::skip]+[curves[-1]]):
        figTitleOut = '%s%s%s.png' % (dirOut,title.replace(' ',''),str(t*skip).zfill(3))
        titleOut = 'Day %s %s' % (str(t*skip).zfill(4),title)
        plotCpleth([day],
                   shapes,
                   var='EPITOTAL',
                   size=7,
                   scheme=None,
                   title=titleOut,
                   figTitle=figTitleOut,
                   legend=True,
                   fx=fx)
    

def plotCpleth(data,
               shapes,
               var='EPITOTAL',
               size=7,
               scheme='QUANTILES',
               title='',
               legend=True,
               showAxes=False,
               figTitle='null',
               fx='null'):
    
    """Plots chloropleth"""
    names = shapes.index
    if fx == 'null':
        getInfected = lambda x: 1-(x['s']/(x['s']+x['e']+x['i']+x['r']))
    else:
        getInfected = eval('lambda x: '+fx)
    epiTotals = [sum([getInfected(t[name]) for t in data]) for name in names]
    shapes['EPITOTAL'] = epiTotals
    
    plt.figure(figsize=(size,size))
    shapes.plot(column=var,scheme=scheme, k=9, colormap='OrRd')
    
    if scheme is None:
        scheme = 'Default'
    if not showAxes:
        plt.axis('off')
    
    plotTitle = "%s %s by patch via %s" % (title, var.title(),scheme.lower())
    if figTitle == 'null':
        figTitle = plotTitle.title().replace(' ','')+'.png'
            
    plt.title(plotTitle)
    plt.savefig(figTitle)
    plt.show()
    
    
def plotExp(targetPlot,expPlots,figureArgs={}):
    """Line graph of experiment results with optimizing progress"""
    args = {'ylabel':'New Infections Daily',
             'xlabel':'Day',
             'title':'Calibration Scenario',
             'figsize':(10,7),
             'dpi':600,
             'ylim':'null'}
    if figureArgs != 'null':
        args.update(figureArgs)
    plt.figure(figsize=args['figsize'],
              dpi=args['dpi'])
    plt.plot(targetPlot,color='r',zorder=3)
    yMax = max(targetPlot)
    opacity = min(15./len(expPlots),0.3)
    for entry in expPlots[:-1]:
        plt.plot(entry,alpha=opacity,color='g',zorder=1)
    plt.plot(expPlots[-1],color='b',zorder=2)
    plt.ylabel(args['ylabel'])
    plt.xlabel(args['xlabel'])
    plt.title(args['title'])
    if args['ylim'] == 'null':
        plt.ylim((0,yMax*1.3))
    else:
        plt.ylim((0,args['ylim']))


### Denaturing Algorithm script parsers

def shiftVars(expVars,select='null',mute = 1):
    """Shifts experimental vars incrementally"""
    if expVars == {}:
        return expVars
    newVars = deepcopy(expVars)
    if select != 'null':
        iterator = ra.sample(expVars.items(),eval(select))
    else:
        iterator = ra.sample(expVars.items(),ra.randint(1,len(expVars)/2))
        
    for key,(base,change,lower,upper) in iterator:
        #if int(base) != base:
        if True:
            shifted = ra.uniform(max(base-change*base*mute,lower),min(base+change*base*mute,upper))
        else:
            M = max(change*mute,1)
            shifted = ra.randint(max(base-M,lower),min(base+M,upper))
        newVars[key] = (shifted,change,lower,upper)
    return newVars


def prepScript(template,expVars):
    """Prepares exp script based on exp vars"""
    varDict = {item[0]:str(item[1][0]) for item in expVars.items()}
    newScript = template % varDict
    return newScript


### 

