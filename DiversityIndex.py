# -*- coding: utf-8 -*-
"""
Created on Tue Feb  8 19:02:16 2022

@author: vmcx3
"""
import pandas as pd
import numpy as np
import math
from scipy.stats import norm
from matplotlib import pyplot as plt
import math
from datetime import datetime
import random
import copy

df = pd.read_csv('C:\\Users\\vmcx3\\Desktop\\Data sets\\Electricity data\\TOSG_cont\\2016.csv')

Meters = pd.unique(df["dataid"])
Meter_Readings = df["use"]*1000 
Meter_Times = df["localminute"]
date_format = "%Y-%m-%d %H:%M:%S"
Total_readings = len(df)

N = len(Meters)                                            # Number of Meters  
M=60;                                                      # Number of compromised meters
SW = 100                                                   # SW is the species width in watts.
R = 0                                                      # Initialising richness (Number of unique species) to 0.
q = 0.5

Meter_Readings[Meter_Readings > 3000] = 3000
Meter_Readings[Meter_Readings <= 0] = 50                   #adjusting min and max values

Maximum_Reading = max(Meter_Readings);                     # Retrieves maximum reading to find richness.
Minimum_Reading = min(Meter_Readings);                     # Retrieves minimum reading to find richness.
R = math.ceil((Maximum_Reading-Minimum_Reading)/SW);       # Calculates the value R (Number of unique species)

Frame_Size = 20                                             # Frame size is in number of days where we calculate Diversity Index
max_time = datetime.strptime(max(Meter_Times), date_format)
min_time = datetime.strptime(min(Meter_Times), date_format)
F = math.ceil((max_time-min_time).days/Frame_Size)      #total number of frames in data


Hist_D_Count = np.zeros((F,N,R))

for i in range(Total_readings):                           # Declaring matrix to store probability of R species of N meters
    current_time = datetime.strptime(Meter_Times[i], date_format)
    Frame_Number = math.ceil((current_time-min_time).days/Frame_Size) - 1
    Meter_Number = np.where(Meters == df["dataid"][i])[0][0]
    Species_Number = math.ceil(Meter_Readings[i]/SW) - 1
    Hist_D_Count[Frame_Number, Meter_Number, Species_Number] = Hist_D_Count[Frame_Number, Meter_Number, Species_Number] + 1
    
Hist_D_Prob = np.zeros((F,N,R))                          # Declaring matrix to store probability of R species of N meters

for i in range(F):
    for j in range(N):
        for k in range(R):
            Hist_D_Prob[i, j, k] = (Hist_D_Count[i, j, k] + 1)/(sum(Hist_D_Count[i, j, :])  + R)
            
            






Affected = random.sample(range(0, M), M)                           # Select M random numbers for selecting columns of attacked meters. Change M to N if you want true random instead of first M
Affected.sort()                                                    # Sorts the vector in previous step and loads to new vector affected_meters
Affected_Meters = Meters[Affected]

dmin=120;                                                           # Loading Delta min for attack
dmax=200;                                                           # Loading Delta max for attack


House_Partition = {}

for i in range(len(Meter_Readings)):
    if df['dataid'][i] in House_Partition:
        House_Partition[df['dataid'][i]].append(Meter_Readings[i])
    else:
        House_Partition[df['dataid'][i]] = [Meter_Readings[i]]
    
    
Meterwise_Unaffected = copy.deepcopy(House_Partition)


for i in House_Partition:
    if i in Affected_Meters:
        for j in range(len(House_Partition[i])):
            House_Partition[i][j] = House_Partition[i][j] + np.random.randint(dmin,dmax)
        
Meterwise_Affected = copy.deepcopy(House_Partition)   
Meter_Readings_Affected = copy.deepcopy(Meter_Readings)   

for i in range(len(Meter_Readings)):
    Meter_Readings_Affected[i] = House_Partition[df['dataid'][i]][0]
    House_Partition[df['dataid'][i]].pop(0)
   
Meter_Readings_Affected[Meter_Readings_Affected > 3000] = 3000
Meter_Readings_Affected[Meter_Readings_Affected <= 0] = 50                   #adjusting min and max values   


Current_D_Count = np.zeros((F,N,R))

for i in range(Total_readings):                           # Declaring matrix to store probability of R species of N meters
    current_time = datetime.strptime(Meter_Times[i], date_format)
    Frame_Number = math.ceil((current_time-min_time).days/Frame_Size) - 1
    Meter_Number = np.where(Meters == df["dataid"][i])[0][0]
    Species_Number = math.ceil(Meter_Readings_Affected[i]/SW) - 1
    Current_D_Count[Frame_Number, Meter_Number, Species_Number] = Current_D_Count[Frame_Number, Meter_Number, Species_Number] + 1
    
Current_D_Prob = np.zeros((F,N,R))                          # Declaring matrix to store probability of R species of N meters

for i in range(F):
    for j in range(N):
        for k in range(R):
            Current_D_Prob[i, j, k] = (Current_D_Count[i, j, k] + 1)/(sum(Current_D_Count[i, j, :])  + R)


Ab=0.29;
Bb=0.08;
Nu=0.05;                            # sigmoid parameters


Hist_Frame = 1
Current_Frame = 2

D_f = {}

for i in range(N):
    for j in range(R):
        temp = abs(Current_D_Prob[Current_Frame, i, j] - Hist_D_Prob[Hist_Frame, i, j])
        temp = 1/(math.pow(1+Ab*math.pow(math.e,(-1*Bb*100*temp)), 1/Nu))
        if Meters[i] in D_f:
            D_f[Meters[i]].append(temp)
        else:
            D_f[Meters[i]] = [temp]
            
E_D = {}

for i in D_f:
    for j in range(R): 
        if i in E_D:
            E_D[i].append(math.pow(D_f[i][j] * Hist_D_Prob[Hist_Frame, np.where(Meters == i)[0][0], j], 0.5))
        else:
            E_D[i] = [math.pow(D_f[i][j] * Hist_D_Prob[Hist_Frame, np.where(Meters == i)[0][0], j], 0.5)]
            
            
DI = {}

for i in E_D:
    for j in range(R): 
        if i in DI:
            DI[i].append(E_D[i][j] * (1 - Hist_D_Prob[Hist_Frame, np.where(Meters == i)[0][0], j]))
        else:
            DI[i] = [E_D[i][j] * (1 - Hist_D_Prob[Hist_Frame, np.where(Meters == i)[0][0], j])]
            

for i in DI:
    DI[i] = sum(DI[i])
            
               
xax = np.zeros(shape=N)

for i in range(N):
    xax[i] = i

plt.hold(True)
for i in range(N):
    if i<M-1:
        plt.plot(xax[i],DI[Meters[i]],'bo')
    elif i==M-1:
        plt.plot(xax[i],DI[Meters[i]],'bo', label = "Honest Smart Meter")
    elif i==N-1:
        plt.plot(xax[i],DI[Meters[i]],'r*', label = "Malicious Smart Meter")
    else:
        plt.plot(xax[i],DI[Meters[i]],'r*')
        
plt.axhline(y=0.52, color='g', linestyle='-')

plt.xlabel('Smart meter ID')
plt.ylabel('Trust Score')
plt.legend()
#plt.savefig('Attack_Evasion')


        
        
        
        
        
        
        
        
        
        
    


