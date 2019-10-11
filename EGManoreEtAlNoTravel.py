#Model parameters and equations based upon the work of Manore et al.

#Manore, C. A., Hickmann, K. S., Xu, S., Wearing, H. J., & Hyman, J. M. (2014). Comparing dengue and chikungunya emergence and endemic transmission in A. aegypti and A. albopictus. Journal of theoretical biology, 356, 174-191.

#http://www.sciencedirect.com/science/article/pii/S002251931400263X

import numpy as np
print "Importing Manore no travel model\n"

stocks = {'sH','eH','iH','rH','nV','sV','eV','iV'}
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
               'Aedes Aegypti':{'betaHV':0.33,
                               'betaVH':0.33,
                               'psiV':0.30,
                               'sigmaV':0.5,
                               'nuV':1./10,
                               'muV':1./14},
               'Aedes Albopictus':{'betaHV':0.31,
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
               'Aedes Aegypti':{'betaHV':0.24,
                               'betaVH':0.24,
                               'psiV':0.30,
                               'sigmaV':0.5,
                               'nuV':1./11,
                               'muV':1./14},
               'Aedes Albopictus':{'betaHV':0.33,
                               'betaVH':0.33,
                               'psiV':0.24,
                               'sigmaV':0.26,
                               'nuV':1./3.5,
                               'muV':1./21}}


models = {'dengue':dengueModel,
          'chikungunya':chikungunyaModel}
                    

    
def modelPrep(pops,scen,models=models):
    """REQUIRED: Prepares model for execution at start of sim"""
    
    if 'strModel' in scen.keys():
        models = eval(strModel)
        
    disease = scen['disease']
    species = scen['species']
    scen.update(models[disease]['human'])
    scen.update(models[disease][species])
    scen['rV'] = scen['psiV'] - scen['muV']
    scen['sigmaVsigmaH'] = scen['sigmaV']*scen['sigmaH']
    
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
            pops[key][stock] = float(pops[key][stock]):q
        pops[key]['p'] = float(pops[key]['nH'])/N
        
    return pops,scen,N


def modelExec(pops,exps,scen,seeds,useSeeds,t):
    """REQUIRED: Runs one iteration of model"""
    flows = dict()
    for key,pop in pops.iteritems():
        #From Section 2.1, 2.2
        
        lambdaHNum = scen['sigmaVsigmaH']*pop['nV']
        lambdaVNum = scen['sigmaVsigmaH']*pop['nH'] 
        lambdaDenom = scen['sigmaV']*pop['nV']+scen['sigmaH']*pop['nH']
        
        try:
            lambdaH = (lambdaHNum/lambdaDenom)*scen['betaHV']*(pop['iV']/pop['nV'])
        except:
            lambdaH = 0.0
        
        try:
            lambdaV = (lambdaVNum/lambdaDenom)*scen['betaVH']*(pop['iH']/pop['nH'])
        except:
            lambdaV = 0.0
            
        hVnV = scen['psiV'] - (scen['rV']/pop['kV']) * pop['nV']
        
        
        flow = {'sH':scen['psiH']*pop['nH'] - (lambdaH+scen['muH'])*pop['sH'],
                'eH':lambdaH*pop['sH'] - (scen['nuH']+scen['muH'])*pop['eH'],
                'iH':scen['nuH']*pop['eH'] - (scen['gammaH']+scen['muH'])*pop['iH'],
                'rH':scen['gammaH']*pop['iH'] - scen['muH']*pop['rH'],
                'nV':scen['rV'] * (1-pop['nV']/pop['kV'])*pop['nV'],
                'sV':hVnV*pop['nV'] - (lambdaV+scen['muV'])*pop['sV'],
                'eV':lambdaV*pop['sV']-(scen['nuV']+scen['muV'])*pop['eV'],
                'iV':scen['nuV']*pop['eV']-scen['muV']*pop['iV']}

        flows[key] = flow
        
    applyFlows(pops,flows,seeds,useSeeds,t,scen['N'])
    return pops,scen
    

    
def applyFlows(pops,flows,seeds,useSeeds,t,N):
    """Applies flows and seeded infections"""
    if useSeeds:
        try:
            (key,rate) = seeds[t]
            flows[key]['sH'] -= rate/N
            flows[key]['eH'] += rate/N
        except:
            pass
        
    for pop in pops.items():
        key = pop[0]
        item = pop[1]
        for stock in stocks:
            pops[key][stock] = item[stock]+flows[key][stock]
        pops[key]['newInf'] = -flows[key]['sH']

    
