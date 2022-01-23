#!/usr/bin/env python
import os, re
import commands
import math, time
import sys
import argparse
import subprocess

parser = argparse.ArgumentParser(description='This script splits moduleCharacterization_step* tasks in multiple parallel jobs')

parser.add_argument("-l",  "--label",          required=True, type=str, help="job label")
parser.add_argument("-b",  "--baseFolder",     required=True, type=str, help="base folder")
parser.add_argument("-e",  "--exeName",        required=True, type=str, help="absolute path of executable")
parser.add_argument("-r",  "--runs",           required=True, type=str, help="comma-separated list of runs to be processed")
parser.add_argument("-c",  "--configFile",     required=True, type=str, help="config file")
parser.add_argument("-j",  "--jobs",           required=True, type=str, help="indicate the number of CPU to be used")
parser.add_argument("-s",  "--submit",                                  help="submit jobs", action='store_true')
parser.add_argument("-v",  "--verbose",                                 help="increase output verbosity", action='store_true')

args = parser.parse_args()


runs = []
comma_list = args.runs.split(',')
for item in comma_list:
   hyphen_list = item.split('-')
   if len(hyphen_list) > 1:
      for i in range(int(hyphen_list[0]), int(hyphen_list[1])+1):
         runs.append(i)
   else:
      runs.append(int(hyphen_list[0]))
print runs


print 
print 'START'
print 

currDir = os.getcwd()

print

try:
   subprocess.check_output(['mkdir','jobs'])
except subprocess.CalledProcessError as e:
   print e.output
try:
   subprocess.check_output(['mkdir','jobs/'+args.label])
except subprocess.CalledProcessError as e:
   print e.output


parallelCommand = "parallel --bar --jobs "+args.jobs+" --results ./jobs/"+args.label+"/ "+args.baseFolder+"/"+args.exeName+" ::: "  #limit to 2 cpu on lip laptop


##### creates directory and file list for job #######
jobDir = currDir+'/jobs/'+args.label+'/'
os.system('mkdir '+jobDir)
os.chdir(jobDir)

for intrun in runs:
   ##### creates config file #######
   run = '{0:04d}'.format(intrun)
   with open(args.baseFolder+'/'+args.configFile) as fi:
      contents = fi.read()
      configFileName = jobDir+"/config_run"+str(run)+".cfg"
      with open(configFileName, "w") as fo:
         fo.write(contents)
      command = 'sed -i \"s%^runs .*$%runs '+str(run)+'%\" '+configFileName
      os.system(command)
      command = 'sed -i \"s%^step1FileName .*$%step1FileName '+args.baseFolder+'/plots_fede/moduleCharacterization_step1_run'+str(run)+'.root%\" '+configFileName
      os.system(command)
      command = 'sed -i \"s%^outFileNameStep1 .*$%outFileNameStep1 '+args.baseFolder+'/plots_fede/moduleCharacterization_step1_run'+str(run)+'.root%\" '+configFileName
      os.system(command)
      command = 'sed -i \"s%^outFileNameStep2 .*$%outFileNameStep2 '+args.baseFolder+'/plots_fede/moduleCharacterization_step2_run'+str(run)+'.root%\" '+configFileName
      os.system(command)
      command = 'sed -i \"s%^plotDir .*$%plotDir /var/www/html/MTDST_CERN_Oct21/CCv2/ModuleCharacterization//run'+str(run)+'/%\" '+configFileName
      os.system(command)
      
   parallelCommand += configFileName + " "


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
os.system("chmod 755 "+jobDir+"/jobs.sh")
   
if args.submit:
   os.system("source "+jobDir+"/jobs.sh")

print
print 'END'
print
