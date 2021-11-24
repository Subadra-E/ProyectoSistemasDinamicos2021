# Código que genera los datos del modelo estocástico

import numpy as np
import matplotlib.pyplot as plt
import matplotlib

import pandas as pd

from tqdm import tqdm


#Agent based

_nAgents = 10000 #population
state = np.zeros(_nAgents)
data = {"state": state}
df = pd.DataFrame(data)
df.describe()

def infect(df, contacts, probability=1.0):
    unique, counts = np.unique(contacts, return_counts=True)
    
    roll = np.random.uniform(0,1,len(unique))

    state = df.loc[unique,"state"]
    probability = 1 - np.power(1-probability, counts)
    change = np.array(roll <= probability).astype(int)
    
    df.loc[unique,"state"] = state + np.maximum(change*(1-state),0)
    
    
def init(nAgents=1000, nPatientZero=1):
    state = np.zeros(nAgents)

    neighborhood = np.zeros(nAgents)
    data = {"state": state, "neighborhood": neighborhood}

    df = pd.DataFrame(data)
    patientZero = np.random.choice(df.index, nPatientZero, replace=False)
    infect(df, patientZero, probability=1.0)
    return df

def recover(df, probability=1.0):    
    roll = np.random.uniform(0,1,len(df[df["state"] == 1]))
    chance = np.array(roll <= probability).astype(int)
    
    df.loc[df["state"] == 1,"state"] = 1 + chance

def step(df):
    nInfected = np.sum(df["state"] == 1)
    contacts = np.random.choice(df.index, _randomContacts * nInfected, replace=True)

    infect(df, contacts, _chanceOfInfection)
    recover(df, _chanceOfRecovery)
    
def simulate(df, stats, nSteps=100, mode="random", nRandomContacts=0, plotLattice=False):
    for i in tqdm(range(nSteps)):        
        step(df)
            
        stats["nSusceptible"].append(np.sum(df["state"] == 0))    
        stats["nInfected"].append(np.sum(df["state"] == 1))
        stats["nRemoved"].append(np.sum(df["state"] == 2))
        
  
_randomContacts = 8 #number of contacts (susceptible)
_chanceOfInfection = 0.03 #beta
_daysCuredAfter = 10 
_chanceOfRecovery = 1./_daysCuredAfter #gamma

_nExperiments = 10 #number of times the experiment is repeated 
_nAgents = 10000 #population
_nSteps = 150

_nPatientZero = 5 #initial number of infected population

#days
x = np.linspace(0,_nSteps-1,_nSteps)
allStats = []

for iExp in range(_nExperiments):
    print("Starting Experiment:",iExp+1,"/",_nExperiments)
    st = {"nInfected": [], "nRemoved": [], "nSusceptible": []}

    df = init(_nAgents, _nPatientZero)

    simulate(df, stats=st, nSteps=_nSteps)
    
    allStats.append(st)

def calculateStats(allStats):
    medianStats = dict()
    lowerStats = dict()
    higherStats = dict()

    for key in allStats[0]:
        l = []
        for st in allStats:
            l.append(st[key])
        a = np.stack(l)
        medianStats[key] = np.median(a, axis=0)
        lowerStats[key] = np.quantile(a, 0.25, axis=0)
        higherStats[key] = np.quantile(a, 0.75, axis=0)
    
    return medianStats, lowerStats, higherStats

def plotSIR(x,mdianStats,lowerStats,higherStats,figName="tmp.png"):
    fig=plt.figure(figsize=(8, 6))
    plt.plot(x,medianStats["nSusceptible"], color = "green", label="Susceptible")
    plt.plot(x,medianStats["nInfected"], color="red", label="Infected")
    plt.plot(x,medianStats["nRemoved"], color="blue", label="Recovered")
    plt.fill_between(x, lowerStats["nSusceptible"], higherStats["nSusceptible"],
                     color='green', alpha=0.1)
    plt.fill_between(x, lowerStats["nInfected"], higherStats["nInfected"],
                     color='red', alpha=0.1)
    plt.fill_between(x, lowerStats["nRemoved"], higherStats["nRemoved"],
                     color='blue', alpha=0.1)


    plt.xlabel("Time steps [days]",size=15)
    plt.ylabel("Number of cases",size=15)

    lgd = plt.legend(bbox_to_anchor=(1.01,0.65), loc="center left",fontsize=15)
    plt.tight_layout()
    
    plt.savefig(figName, bbox_extra_artists=(lgd,), bbox_inches='tight')
   
    plt.show()

medianStats, lowerStats, higherStats = calculateStats(allStats)
 
 #In this section the data will be saved in files .txt
with open("MedianStats.txt","w") as f0: 
 for i in range(0, len(list(medianStats["nSusceptible"]))):
    #The data corresponding to MedianStats are stored for each of the cases in three columns 
    f0.write("{0}\t{1}\t{2}\t{3}\n".format(x[i],medianStats["nSusceptible"][i],medianStats["nInfected"][i], medianStats["nRemoved"][i]))
f0.close()

with open("LowerStats.txt","w") as f1: 
 for i in range(0, len(list(lowerStats["nSusceptible"]))):
    #The data corresponding to lowerStats are saved for each of the cases in three columns 
    f1.write("{0}\t{1}\t{2}\t{3}\n".format(x[i],lowerStats["nSusceptible"][i],lowerStats["nInfected"][i], lowerStats["nRemoved"][i]))
f1.close()

with open("HigherStats.txt","w") as f2: 
 for i in range(0, len(list(higherStats["nSusceptible"]))):
    #The data corresponding to HigherStats are stored for each of the cases in three columns 
    f2.write("{0}\t{1}\t{2}\t{3}\n".format(x[i],higherStats["nSusceptible"][i],higherStats["nInfected"][i], higherStats["nRemoved"][i]))
f2.close()

plotSIR(x,medianStats,lowerStats,higherStats,figName="SimulacionEst1.png")

