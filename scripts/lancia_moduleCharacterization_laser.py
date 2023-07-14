#!/usr/bin/env python
import os, re
import commands
import math, time
import sys
import argparse
import subprocess



multipleRuns = False
vth = 'vth1'

# ----  TOFHIR 2X
# dataFolder = '/data/Lab5015Analysis/moduleCharacterization/TOFHIR2X/'
# plotFolder = '/var/www/html/TOFHIR2X/ModuleCharacterization'

# ----  TOFHIR 2C
dataFolder = '/home/cmsdaq/Lab5015Analysis_new/904/Lab5015Analysis/plots/'      #    904 measurements
plotFolder = '/var/www/html/TOFHIR2C/904/ModuleCharacterization/'


runRanges = [
    #(it,it) for it in range(163,169,1)
    #(it,it) for it in range(175,235,1)
    #(it,it) for it in range(235,289,1)
    #(it,it) for it in range(289,343,1)
    #(it,it) for it in range(439,493,1)
    #(it,it) for it in range(554,626,1)
    #(it,it) for it in range(626,644,1)
    #    (it,it) for it in range(979,1094,1)
    (it,it) for it in range(1094,1100,1)
]




runList=''
for runRange in runRanges:
    if runRange[1] == runRange[0]:
        runList += str(runRange[0])+','
    else:
        runList+=str(runRange[0])+'-'+str(runRange[1])+','
runList = runList[:-1]
print(runList)




parser = argparse.ArgumentParser(description='This script splits submits moduleCharacterization jobs')
parser.add_argument("-s", "--submit", help="submit jobs", action='store_true')
args = parser.parse_args()




#########
### step1
#command = './submit_moduleCharacterization_TOFHIR2.py -l test -b /home/cmsdaq/Programs/Lab5015Analysis/ -e bin/moduleCharacterization_step1.exe -r '+runList+' -c cfg/moduleCharacterization.cfg  --vth '+vth+' -i '+dataFolder+' -o '+dataFolder+' --submit -p 4'
command = './submit_moduleCharacterization_TOFHIR2.py -l test -b /home/cmsdaq/Lab5015Analysis_new/904/Lab5015Analysis/ -e bin/moduleCharacterization_step1.exe -r '+runList+' -c cfg/moduleCharacterization_laser.cfg  --vth '+vth+' -i '+dataFolder+' -o '+dataFolder+' --submit -p 4'
print command
if args.submit:
    os.system(command)




#############
### hadd fils
if multipleRuns:
    for runRange in runRanges:
        runs = str(runRange[0])+'-'+str(runRange[1])
        command = 'hadd -f '+dataFolder+'/moduleCharacterization_step1_run'+runs+'.root '
        for run in range(runRange[0],runRange[1]+1):
            command += dataFolder+'/moduleCharacterization_step1_run'+str(run)+'.root '
        print command
        if args.submit:
            os.system(command)




#########
### step2
for runRange in runRanges:
    if runRange[1] == runRange[0]:
        runs = str(runRange[0])
    else:
        runs = str(runRange[0])+'-'+str(runRange[1])
    command = 'mkdir '+plotFolder+'/run'+runs
    if args.submit:
        os.system(command)

runCommand = '-r'
if multipleRuns: runCommand = '-R'
#command = './submit_moduleCharacterization_TOFHIR2.py -l test -b /home/cmsdaq/Programs/Lab5015Analysis/ -e bin/moduleCharacterization_step2.exe '+runCommand+' '+runList+' -c cfg/moduleCharacterization.cfg  --vth '+vth+' -i '+dataFolder+' -o '+dataFolder+' --submit -p 4'
command = './submit_moduleCharacterization_TOFHIR2.py -l test -b /home/cmsdaq/Lab5015Analysis_new/904/Lab5015Analysis/ -e bin/moduleCharacterization_step2.exe '+runCommand+' '+runList+' -c cfg/moduleCharacterization_laser.cfg  --vth '+vth+' -i '+dataFolder+' -o '+dataFolder+' --submit -p 4'
print command
if args.submit:
    os.system(command)




#################
### summary plots
for runRange in runRanges:
    if runRange[1] == runRange[0]:
        runs = str(runRange[0])
    else:
        runs = str(runRange[0])+'-'+str(runRange[1])
    #command = 'python /home/cmsdaq/Programs/Lab5015Analysis/macros/moduleCharacterizationSummaryPlots.py -m 2 -i run'+runs+' -o run'+runs
    command = 'python /home/cmsdaq/Lab5015Analysis_new/904/Lab5015Analysis/macros/moduleCharacterizationSummaryPlots.py -m 2 -i run'+runs+' -o run'+runs
    print command
    if args.submit:
        os.system(command)
