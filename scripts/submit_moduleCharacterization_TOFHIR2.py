#!/usr/bin/env python
import os, re
import commands
import math, time
import sys
import argparse
import subprocess

parser = argparse.ArgumentParser(description='This script splits moduleCharacterization_step* tasks in multiple parallel jobs')

parser.add_argument("-l",  "--label",        required=True,  type=str, help="job label")
parser.add_argument("-i",  "--inputFolder",  required=True,  type=str, help="input file folder")
parser.add_argument("-o",  "--outputFolder", required=True,  type=str, help="output file folder")
parser.add_argument("-b",  "--baseFolder",   required=True,  type=str, help="base folder")
parser.add_argument("-e",  "--exeName",      required=True,  type=str, help="absolute path of executable")
parser.add_argument("-r",  "--runs",         required=False,  type=str, help="comma-separated list of runs to be processed")
parser.add_argument("-R",  "--ranges",       required=False,  type=str, help="comma-separated list of ranges to be processed")
parser.add_argument("-c",  "--configFile",   required=True,  type=str, help="config file")
parser.add_argument("-t",  "--vth",          required=False, type=str, help="threshold type (vth1 or vth2)")
parser.add_argument("-s",  "--submit",                                help="submit jobs", action='store_true')
parser.add_argument("-v",  "--verbose",        		              help="increase output verbosity", action='store_true')
parser.add_argument("-p",  "--numParal",     required=True, type=int, help='number of parallel jobs, 0 for all')

args = parser.parse_args()

preruns = []

runs = []
if args.runs != None:
   comma_list = args.runs.split(',')
   for item in comma_list:
      hyphen_list = item.split('-')
      if len(hyphen_list) > 1:
         for i in range(int(hyphen_list[0]), int(hyphen_list[1])+1):
            runs.append(i)
      else:
         runs.append(int(hyphen_list[0]))
if args.ranges != None:
   comma_list = args.ranges.split(',')
   for item in comma_list:
      runs.append(item)
print runs



print 
print 'START'
print 


currDir = os.getcwd()

print

try:
   subprocess.check_output(['mkdir',args.baseFolder+'/scripts/jobs'])
except subprocess.CalledProcessError as e:
  print e.output
try:
   subprocess.check_output(['mkdir',args.baseFolder+'/scripts/jobs/'+args.label])
except subprocess.CalledProcessError as e:
   print e.output
      
   
parallelCommand = "parallel --bar --jobs "+str(args.numParal)+" --results "+args.baseFolder+"/scripts/jobs/"+args.label+"/ "+args.baseFolder+"/"+args.exeName+" ::: "


##### creates directory and file list for job #######
jobDir = args.baseFolder+'/scripts/jobs/'+args.label+'/'
os.system('mkdir '+jobDir)
os.chdir(jobDir)

for run in runs:
   ##### creates config file #######
   with open(args.baseFolder + '/'+args.configFile) as fi:
      contents = fi.read()
      configFileName = jobDir+"/config_run"+str(run)+".cfg"
      with open(configFileName, "w") as fo:
         fo.write(contents)
      command = 'sed -i \"s%^runs .*$%runs '+str(run)+'%\" '+configFileName
      os.system(command)
      command = 'sed -i \"s%^step1FileName .*$%step1FileName '+args.inputFolder+'/moduleCharacterization_step1_run'+str(run)+'.root%\" '+configFileName
      os.system(command)
      command = 'sed -i \"s%^outFileNameStep1 .*$%outFileNameStep1 '+args.outputFolder+'/moduleCharacterization_step1_run'+str(run)+'.root%\" '+configFileName
      os.system(command)
      command = 'sed -i \"s%^outFileNameStep2 .*$%outFileNameStep2 '+args.outputFolder+'/moduleCharacterization_step2_run'+str(run)+'.root%\" '+configFileName
      os.system(command)
      command = 'sed -i \"s%^outFileNameStep3 .*$%outFileNameStep3 '+args.outputFolder+'/moduleCharacterization_step3_run'+str(run)+'.root%\" '+configFileName
      os.system(command)
      if args.vth != None:
         command = 'sed -i \"s%^vth .*$%vth '+args.vth+'%\" '+configFileName
         os.system(command)
      # command = 'sed -i \"s%^plotDir .*$%plotDir '+'/var/www/html/TOFHIR2B/ModuleCharacterization/run'+str(run)+ '%\" '+configFileName
      command = 'sed -i \"s%^plotDir .*$%plotDir '+'/var/www/html/TOFHIR2X/ModuleCharacterization/run'+str(run)+ '%\" '+configFileName
      #command = 'sed -i \"s%^plotDir .*$%plotDir '+'/var/www/html/TOFHIR2A/MTDTB_CERN_Jul21/ModuleCharacterization/run'+str(run)+ '%\" '+configFileName
      os.system(command)
      
      parallelCommand += configFileName + " "
      
print(parallelCommand)
        

##### creates job file #######
with open(jobDir+'/jobs.sh', 'w') as fout:
   fout.write("#!/bin/sh\n")
   fout.write("echo\n")
   fout.write("echo 'START---------------'\n")
   fout.write("echo 'current dir: ' ${PWD}\n")
   fout.write("cd "+str(args.baseFolder)+"\n")
   fout.write("echo 'current dir: ' ${PWD}\n")
   fout.write("source scripts/setup.sh\n")
   fout.write(parallelCommand+"\n")
   fout.write("echo 'STOP---------------'\n")
   fout.write("echo\n")
   fout.write("echo\n")
   
if args.submit:
   os.system("source "+jobDir+"/jobs.sh")

print
print 'END'
print
