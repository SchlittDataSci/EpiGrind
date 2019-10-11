import numpy as np

stocks = {'s','e','i','r'}

def modelPrep(pops,scen):
    """REQUIRED: Prepares model for execution at start of sim"""
    N = 0.0
    for key,item in pops.iteritems():
        if 'n' not in item.keys():
            pops[key]['n'] = sum(item.values())
        N += pops[key]['n']
        
    for key,item in pops.iteritems():
        for stock in stocks:
            pops[key][stock] = item[stock]/N
        pops[key]['p'] = float(pops[key]['n'])/N
        
        #Assume 72 years lifespan if no lifespan in data
        try:
            pops[key]['u'] = 1./item['tLife']
        except:
            pops[key]['u'] = 1./26280
            
        #Preps gamma and theta for SEIR model
        pops[key]['gamma'] =  prepRate(pops[key],scen,'tDwell')
        pops[key]['theta'] =  prepRate(pops[key],scen,'tInfect')
        scen['N'] = N
    return pops,scen,N


def modelExec(pops,exps,scen,seeds,useSeeds,t):
    """REQUIRED: Runs one iteration of model"""
    flows = dict()
    infectivity = scen['infectivity']
    for key,item in pops.iteritems():
        u = item['u'] 
        flow = {'s':(item['p']-item['s'])*u,
                'e':-item['e']*u,
                'i':-item['i']*u,
                'r':-item['r']*u}
            
        bySource = []

        for exp in exps[key]:
            (source,contactRate) = exp
            pInf = pops[source]['i']/pops[source]['p']
            force = contactRate*infectivity*pInf
            bySource.append(1-min(force,1))

        trueExp = item['s']*(1-np.prod(bySource))
        newInf = item['e']*item['gamma']
        newRec = item['i']*item['theta']
                
        flow['s']-=trueExp
        flow['e']+=trueExp-newInf
        flow['i']+=newInf-newRec
        flow['r']+=newRec
        flows[key] = flow
            
    applyFlows(pops,flows,seeds,useSeeds,t,scen['N'])
    return pops,scen
    

    
def applyFlows(pops,flows,seeds,useSeeds,t,N):
    """Applies flows and seeded infections"""
    if useSeeds:
        try:
            (key,rate) = seeds[t]
            flows[key]['s'] -= rate/N
            flows[key]['e'] += rate/N
        except:
            pass
        
    for pop in pops.items():
        key = pop[0]
        item = pop[1]
        for stock in stocks:
            pops[key][stock] = item[stock]+flows[key][stock]
        pops[key]['newInf'] = -flows[key]['s']
    
    
    
def prepRate(compartment,scenario,key):
    """Calculates gamma and theta using scenario value if individual values not available"""
    try:
        duration = compartment[key]
    except:
        duration = scenario[key]
    if duration == 0:
        return 1.
    elif duration == -1:
        return 0.
    return 1./duration
    
