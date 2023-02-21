#!/usr/bin/env python
import os, re
import commands
import math, time
import sys
import argparse
import subprocess

##########################################################
### laser runs

runRanges = [

#----th1 only---- high att ---
#    (4843,4847),
#    (4853,4857)
#    (4858,4862),
#    (4868,4872),
#    (4873,4877),
#    (4883,4887),
#    (4888,4892),
#    (4898,4902),
#    (4903,4907),
#    (4913,4917),
#    (4918,4922),
#    (4928,4932),
#    (4933,4937),
#    (4943,4947),
#    (4948,4952),
#    (4958,4962),
#    (4963,4967),
#    (4973,4977),
#    (4978,4982),
#    (4988,4992)

#----th1 only---- low att --- scaling 0
#    (5154,5158),
#    (5164,5168),
#    (5169,5173),
#    (5179,5183),
#    (5184,5188),
#    (5194,5198),
#    (5199,5203),
#    (5209,5213)

#----th1 only--low att---scaling 1---
#    (5224,5228),
#    (5239,5243),
#    (5254,5258),
#    (5269,5273),
#    (5284,5288)
#    (5299,5303),
#    (5314,5318)
#    (5505, 5505)


#    (5566, 5570),
#    (5576, 5580),
#    (5586, 5590),
#    (5596, 5600)
    
#    (6573,6577),
#    (6583,6587),
#    (6593,6597),
#    (6603,6607),
#    (6613,6617),
#    (6623,6627),
#    (6633,6637),
#    (6643,6647),
#    (6653,6657),
#    (6663,6667),
#    (6673,6677),
#    (6683,6687),
#    (6693,6697),
#    (6703,6707),

#    (6745,6745),
#    (6746,6746),
#    (6753,6753),
#    (6754,6754),
#    (6763,6763),
#    (6764,6764),
#    (6765,6765),
#    (6766,6766),
    (6767,6767),
]






runList=''
for runRange in runRanges:
    runList+=str(runRange[0])+'-'+str(runRange[1])+','
runList=runList[:-1]











#############
##step1, ith1
command = './submit_moduleCharacterization_TOFHIR2.py -l test -b /home/cmsdaq/Programs/Lab5015Analysis/ -e bin/moduleCharacterization_step1.exe -r '+runList+' -c cfg/moduleCharacterization.cfg  --vth vth1 -i /data/Lab5015Analysis/moduleCharacterization/TOFHIR2B/ -o /data/Lab5015Analysis/moduleCharacterization/TOFHIR2B/ --submit -p 30'
print command
os.system(command)

#############
##step1, ith2
#command = './submit_moduleCharacterization_TOFHIR2.py -l test -b /home/cmsdaq/Programs/Lab5015Analysis/ -e bin/moduleCharacterization_step1.exe -r '+runList+' -c cfg/moduleCharacterization.cfg  --vth vth2 -i /data/Lab5015Analysis/moduleCharacterization/TOFHIR2B/ -o /data/Lab5015Analysis/moduleCharacterization/TOFHIR2B/ --submit -p 30'
#print command
#os.system(command)


############
## hadd fils
for runRange in runRanges:
    runs = str(runRange[0])+'-'+str(runRange[1])
    command = 'hadd -f /data/Lab5015Analysis/moduleCharacterization/TOFHIR2B/moduleCharacterization_step1_run'+runs+'.root '
    for run in range(runRange[0],runRange[1]+1):
        command += '/data/Lab5015Analysis/moduleCharacterization/TOFHIR2B/moduleCharacterization_step1_run'+str(run)+'.root '
    print command
    os.system(command)


######
#step2
for runRange in runRanges:
    runs = str(runRange[0])+'-'+str(runRange[1])
    command = 'mkdir /var/www/html/TOFHIR2B/ModuleCharacterization/run'+runs
    os.system(command)
    
command = './submit_moduleCharacterization_TOFHIR2.py -l test -b /home/cmsdaq/Programs/Lab5015Analysis/ -e bin/moduleCharacterization_step2.exe -R '+runList+' -c cfg/moduleCharacterization.cfg  --vth vth2 -i /data/Lab5015Analysis/moduleCharacterization/TOFHIR2B/ -o /data/Lab5015Analysis/moduleCharacterization/TOFHIR2B/ --submit -p 30'
print command
os.system(command)


##############
#summary plots
for runRange in runRanges:
    runList = str(runRange[0])+'-'+str(runRange[1])
    
    command = 'python /home/cmsdaq/Programs/Lab5015Analysis/macros/moduleCharacterizationSummaryPlots.py -m2 -i run'+runList+' -o /var/www/html/TOFHIR2B/ModuleCharacterization/run'+runList+'/' 
    print command
    os.system(command)
