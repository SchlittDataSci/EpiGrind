#Model parameters and equations based upon the work of Manore et al.

#Manore, C. A., Hickmann, K. S., Xu, S., Wearing, H. J., & Hyman, J. M. (2014). Comparing dengue and chikungunya emergence and endemic transmission in A. aegypti and A. albopictus. Journal of theoretical biology, 356, 174-191.

#http://www.sciencedirect.com/science/article/pii/S002251931400263X

import numpy as np
import pandas as pd
from glob import glob
print "Importing Manore with travel model\n"

stocks = {'sH','eH','iH','rH','nV','sV','eV','iV'}
flowStocks = {'sH','eH','iH','rH','sV','eV','iV'}

###Requires prepped experiment with values:
## Stocks:
# sh, eH, iH, rH, nV, sV, eV, iV
## Pop specific quantitities:
# kV

#Daily model parameters from yearly model
#From Table 3
dengueModel = {'human':{'sigmaH':19.,
                        'nuH':1./5,
                        'gammaH':1./6,
                        'muH':1./(70*365.25),
                        'psiH':1./(70*365.25)},
               'Aedes aegypti':{'betaHV':0.33,
                               'betaVH':0.33,
                               'psiV':0.30,
                               'sigmaV':0.5,
                               'nuV':1./10,
                               'muV':1./14},
               'Aedes albopictus':{'betaHV':0.31,
                               'betaVH':0.31,
                               'psiV':0.24,
                               'sigmaV':0.26,
                               'nuV':1./10,
                               'muV':1./21}}

chikungunyaModel = {'human':{'sigmaH':19.,
                        'nuH':1./3,
                        'gammaH':1./6,
                        'muH':1./(70*365.25),
                        'psiH':1./(70*365.25)},
               'Aedes aegypti':{'betaHV':0.24,
                               'betaVH':0.24,
                               'psiV':0.30,
                               'sigmaV':0.5,
                               'nuV':1./11,
                               'muV':1./14},
               'Aedes albopictus':{'betaHV':0.33,
                               'betaVH':0.33,
                               'psiV':0.24,
                               'sigmaV':0.26,
                               'nuV':1./3.5,
                               'muV':1./21}}


models = {'dengue':dengueModel,
          'chikungunya':chikungunyaModel}
            
                    
                                    
def prepTravelMatrix(fileRef,transScale,pops):
    """Pivots travel list to from v to matrix"""
    if type(fileRef) is str:
        travelData = pd.DataFrame.from_csv(glob(fileRef)[0])
    else:
        travelData = fileRef
    #patches = set(travelData.index) | set(travelData['Destination'])
    #patches =  sorted(list(patches))
    patches = sorted(pops.keys())
    expRates = pd.DataFrame(0.,index=patches,columns=patches)
    for index,row in travelData.iterrows():
        value = row['Travel_Rate']
        expRates[index][row['Destination']] = value
    expRates = expRates*transScale
    for i in patches:
        if expRates[i][i] == 0:
            expRates[i][i] = pops[i]['nH']-sum(expRates[i])
            if expRates[i][i] < 0:
                print "ERROR, TOO HIGH OF TRAVEL WEIGHT, U LOSTS ALL THE PEOPLES", expRates[i][j]
    for i in patches:
        sumWeight = float(sum(expRates[i]))
        for j in patches:
            expRates[i][j] = (expRates[i][j]/sumWeight)*pops[i]['nH']
    return expRates


def prepPops(fileRef):
    """Assumes susceptibles are only the exposed pop, but keeps total pop size as 'totalPop' """
    popData = pd.DataFrame.from_csv(glob(fileRef)[0])
    pops = sorted(popData.index.tolist())
    popDict = {}
    for pop in pops:
        popDict[pop] = {'sH':1.*popData['Exposed_Pop'][pop],
                       'eH':0.,
                       'iH':0.,
                       'rH':0.,
                       'sV':2.*popData['Exposed_Pop'][pop],
                       'eV':0.,
                       'iV':0.,
                       'kV':2.*popData['Exposed_Pop'][pop],
                       'totalPop':int(1.*popData['Total_Pop'][pop])}
    
    return popDict    
    
def modelPrep(pops,scen,models=models,showTravel=False):
    """REQUIRED: Prepares model for execution at start of sim"""
    
    if 'modelParams' in scen.keys():
        scen.update(scen['modelParams'])
    else:
        if 'strModel' in scen.keys():
            models = eval(scen['strModel'])

        disease = scen['disease']
        species = scen['species']
        scen.update(models[disease]['human'])
        scen.update(models[disease][species])
    
    if type(pops) is str:
        pops = prepPops(pops)
        
    N = 0.0
    for key,item in pops.iteritems():
        if 'nH' not in item.keys():
            pops[key]['nH'] = item['sH']+item['eH']+item['iH']+item['rH']
        if 'nV' not in item.keys():
            pops[key]['nV'] = item['sV']+item['eV']+item['iV'] 
        N += pops[key]['nH']
    scen['N'] = N
        
    for key,item in pops.iteritems():
        for stock in stocks:
            pops[key][stock] = float(pops[key][stock])
        pops[key]['p'] = float(pops[key]['nH'])/N
    
    poplist = pops.keys()
    
    if 'transScale' in scen.keys():
        transScale = scen['transScale']
    else:
        transScale = 1.
    
    if 'travelMatrix' not in scen.keys():
        if 'travelFile' in scen.keys():
            scen['travelMatrix'] = pd.DataFrame.from_csv(scen['travelFile'])
            if len(scen['travelMatrix'].columns)>1:
                scen['travelMatrix'] = prepTravelMatrix(scen['travelMatrix'],transScale,pops)
        else:
            scen['travelMatrix'] = pd.DataFrame(0,index=poplist,columns=poplist)
            for pop in poplist:
                scen['travelMatrix'][pop][pop] = pops[pop]['nH']

    
    travel = scen['travelMatrix'] 
    if showTravel:
        print "Travel Matrix:"
        print travel
    
    scen['nIn'] = dict()
    scen['nOut'] = dict()
    scen['nEff'] = dict()
    scen['pHV'] = dict()
    scen['pVH'] = dict()
    scen['pops'] = pops
    
    compartmentalize = {'sigmaH','sigmaV','psiV','muV'}
    
    for pop in poplist:
        selfTravel = travel[pop][pop]
        popSize = travel[pop].sum()
        scen['nOut'][pop] = popSize - selfTravel
        scen['nIn'][pop] = travel.loc[pop].sum() - selfTravel
        scen['nEff'][pop] = selfTravel + scen['nIn'][pop]
        scen['pHV'][pop] = {i:float(travel[i][pop])/scen['nEff'][pop] for i in pops}
        scen['pVH'][pop] = {i:float(travel[pop][i])/popSize for i in pops}
        
        for variable in compartmentalize:
            if variable not in pops[pop].keys():
                pops[pop][variable] = scen[variable]
                

    return pops,scen,N


def modelExec(pops,exps,scen,shifts,useShifts,t):
    """REQUIRED: Runs one iteration of model"""
    flows = dict()
    for key,popIn in pops.iteritems():
        #From Section 2.1, 2.2
        
        lambdaHByPop = []
        lambdaVByPop = []
        
        iHWeighted = sum([pops[i]['iH']*scen['pHV'][key][i] for i in scen['pops']])
        
        for key2 in scen['pops']:
            popOut = pops[key2]
            nHIn = scen['nEff'][key]
            nHOut = scen['nEff'][key2]
            
            lambdaHNum = popOut['sigmaV']*popOut['sigmaH']*popOut['nV']
            lambdaHDenom = popOut['sigmaV']*popOut['nV']+popOut['sigmaH']*nHOut
            
            lambdaVNum = popIn['sigmaV']*popIn['sigmaH']*nHIn
            lambdaVDenom = popIn['sigmaV']*popIn['nV']+popIn['sigmaH']*nHIn
            
            try:
                lambdaHi = (lambdaHNum/lambdaHDenom)*scen['betaHV']*(popOut['iV']/popOut['nV'])
            except Exception as e:
                lambdaHi = 0.0

            try:
                lambdaVi = (lambdaVNum/lambdaVDenom)*scen['betaVH']*(iHWeighted/nHIn)
            except Exception as e:
                lambdaVi = 0.0
                
            lambdaHByPop.append(lambdaHi*scen['pVH'][key][key2])
            lambdaVByPop.append(lambdaVi*scen['pHV'][key][key2])
        
        lambdaH = sum(lambdaHByPop)
        lambdaV = sum(lambdaVByPop)
        
        popIn['rV'] = popIn['psiV'] - popIn['muV']
        hV = popIn['psiV'] - (popIn['rV']/popIn['kV']) * popIn['nV']
        
        flow = {'sH':scen['psiH']*popIn['nH'] - (lambdaH+scen['muH'])*popIn['sH'],
                'eH':lambdaH*popIn['sH'] - (scen['nuH']+scen['muH'])*popIn['eH'],
                'iH':scen['nuH']*popIn['eH'] - (scen['gammaH']+scen['muH'])*popIn['iH'],
                'rH':scen['gammaH']*popIn['iH'] - scen['muH']*popIn['rH'],
                'sV':hV*popIn['nV'] - (lambdaV+popIn['muV'])*popIn['sV'],
                'eV':lambdaV*popIn['sV']-(scen['nuV']+popIn['muV'])*popIn['eV'],
                'iV':scen['nuV']*popIn['eV']-popIn['muV']*popIn['iV']}
        
        flows[key] = flow
        
    applyFlows(pops,scen,flows,shifts,useShifts,t)
    return pops,scen

    

    
def applyFlows(pops,scen,flows,shifts,useShifts,t):
    """Applies flows and seeded infections"""
    if useShifts:
        try:
            actions = shifts[t]
            exec(actions)
            print "Applied scripted action at t=%s" % t
        except:
            pass
        try:
            exec(shifts['all'])
        except:
            pass
        
    for pop in pops.items():
        key = pop[0]
        item = pop[1]
        for stock in flowStocks:
            pops[key][stock] = item[stock]+flows[key][stock]
        pops[key]['newInf'] = -flows[key]['sH']
        pops[key]['nV'] = item['sV']+item['eV']+item['iV']

    
