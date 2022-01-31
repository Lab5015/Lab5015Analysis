#!/usr/bin/env python
import os, re
import commands
import math, time
import sys
import argparse
import subprocess
import collections

parser = argparse.ArgumentParser(description='This script runs the moduleCharacterization steps and the drawPulseShape analysis')

parser.add_argument("-r",         "--rangeFile", required=True, type=str, help="input file containing the run ranges to be processed")
parser.add_argument("-s",         "--submit",                             help="submit jobs", action='store_true')
args = parser.parse_args()



baseFolder    = "/home/data/mtd/RUptoro2/Lab5015Analysis_MTDST_CERN_Oct21"
jobsFolder    = baseFolder+"/jobs/fullAnalysis"

moduleCharCfg = baseFolder+"/cfg/moduleCharacterization_fede.cfg"
moduleCharExe1= baseFolder+"/bin/moduleCharacterization_step1.exe"
moduleCharExe2= baseFolder+"/bin/moduleCharacterization_step2.exe"

drawPulseCfg  = baseFolder+"/cfg/drawPulseShape_fede.cfg"
drawPulseExe  = baseFolder+"/bin/drawPulseShape_fede.exe"



parallelModChar1  = "parallel --bar --jobs 8 --results "+jobsFolder+"/ " + moduleCharExe1 + " ::: "  #limit to 8 cpu on ceacmsfw server
parallelModChar2  = "parallel --bar --jobs 8 --results "+jobsFolder+"/ " + moduleCharExe2 + " ::: "  #limit to 8 cpu on ceacmsfw server
parallelPulseShape = "parallel --bar --jobs 8 --results "+jobsFolder+"/ " + drawPulseExe + " ::: "   #limit to 8 cpu on ceacmsfw server



#parse input file
runs_dict = {}
ch1 = -1
ch2 = -1
pngName = "test.png"
cfgFile = open(baseFolder+'/'+args.rangeFile)
lines = cfgFile.readlines()

for line in lines:
   if line[0] == "#":
      continue

   line = line.strip()
   if line == "":
      continue

   line = line.split()

   if "ch1" in line:
      ch1 = line[1]
      continue
   if "ch2" in line:
      ch2 = line[1]
      continue
   if "pngName" in line:
      pngName = line[1]
      continue

   runs_dict[line[0]] = [line[1], float(line[2]), list(map(float, line[3].split(",")))]


for run_range, params in sorted(runs_dict.items()):
   os.system("mkdir -p "+jobsFolder)

   configFileName = ""

   #prepare moduleCharacterization config
   with open(moduleCharCfg) as fi:
      contents = fi.read()
      configFileName = jobsFolder+"/configModChar_run"+str(run_range)+".cfg"
      with open(configFileName, "w") as fo:
         fo.write(contents)
      command = 'sed -i \"s%^runs .*$%runs '+str(run_range)+'%\" '+configFileName
      os.system(command)
      command = 'sed -i \"s%^step1FileName .*$%step1FileName '+baseFolder+'/plots_fede/moduleCharacterization_step1_run'+str(run_range)+'.root%\" '+configFileName
      os.system(command)
      command = 'sed -i \"s%^outFileNameStep1 .*$%outFileNameStep1 '+baseFolder+'/plots_fede/moduleCharacterization_step1_run'+str(run_range)+'.root%\" '+configFileName
      os.system(command)
      command = 'sed -i \"s%^outFileNameStep2 .*$%outFileNameStep2 '+baseFolder+'/plots_fede/moduleCharacterization_step2_run'+str(run_range)+'.root%\" '+configFileName
      os.system(command)
      command = 'sed -i \"s%^plotDir .*$%plotDir /var/www/html/MTDST_CERN_Oct21/CCv2/ModuleCharacterization//run'+str(run_range)+'/%\" '+configFileName
      os.system(command)
      command = 'sed -i \"s%^channelMapping .*$%channelMapping '+ch1+' '+ch2+'%\" '+configFileName
      os.system(command)
      

   parallelModChar1 += configFileName + " " 
   parallelModChar2 += configFileName + " " 
   submitCommandModChar = moduleCharExe1 +" "+ configFileName +" 0; "+moduleCharExe2 +" "+ configFileName +" 0"


   #prepare drawPulseShape config
   with open(drawPulseCfg) as fi:
      contents = fi.read()
      configFileName = jobsFolder+"/configDrawPulse_run"+str(run_range)+".cfg"
      with open(configFileName, "w") as fo:
         fo.write(contents)
      command = 'sed -i \"s%^runs .*$%runs '+str(run_range)+'%\" '+configFileName
      os.system(command)
      command = 'sed -i \"s%^step1FileName .*$%step1FileName '+baseFolder+'/plots_fede/pulseShape_run'+str(run_range)+'.root%\" '+configFileName
      os.system(command)
      command = 'sed -i \"s%^outFileNameStep1 .*$%outFileNameStep1 '+baseFolder+'/plots_fede/pulseShape_run'+str(run_range)+'.root%\" '+configFileName
      os.system(command)
      command = 'sed -i \"s%^outFileNameStep2 .*$%outFileNameStep2 '+baseFolder+'/plots_fede/pulseShape_run'+str(run_range)+'.root%\" '+configFileName
      os.system(command)
      command = 'sed -i \"s%^plotDir .*$%plotDir /var/www/html/TOFHIR2X/MTDST_CERN_Oct21/pulseShapes/run'+str(run_range)+'/%\" '+configFileName
      os.system(command)
      command = 'sed -i \"s%^ch1 .*$%ch1 '+ch1+'%\" '+configFileName
      os.system(command)
      command = 'sed -i \"s%^ch2 .*$%ch2 '+ch2+'%\" '+configFileName
      os.system(command)

   parallelPulseShape += configFileName + " " 
   submitCommandPulseShape = drawPulseExe +" "+ configFileName +" 0"


##### creates job file #######
#jobs_modChar = "jobs_modChar_run" + run_range + ".sh"
jobs_modChar = "jobs_modChar.sh"
with open(jobsFolder+'/'+jobs_modChar, 'w') as fout:
   fout.write("#!/bin/sh\n")
   fout.write("echo\n")
   fout.write("echo 'START---------------'\n")
   fout.write("echo 'current dir: ' ${PWD}\n")
   fout.write("cd "+str(baseFolder)+"\n")
   fout.write("echo 'current dir: ' ${PWD}\n")
   fout.write("source scripts/setup.sh\n")
   #fout.write(submitCommandModChar+"\n")
   fout.write(parallelModChar1+"\n")
   fout.write(parallelModChar2+"\n")
   fout.write("echo 'STOP---------------'\n")
   fout.write("echo\n")
   fout.write("echo\n")
os.system("chmod 755 "+jobsFolder+"/"+jobs_modChar)

##### creates job file #######
#jobs_drawPS = "jobs_drawPS_run" + run_range + ".sh"
jobs_drawPS = "jobs_drawPS.sh"
with open(jobsFolder+'/'+jobs_drawPS, 'w') as fout:
   fout.write("#!/bin/sh\n")
   fout.write("echo\n")
   fout.write("echo 'START---------------'\n")
   fout.write("echo 'current dir: ' ${PWD}\n")
   fout.write("cd "+str(baseFolder)+"\n")
   fout.write("echo 'current dir: ' ${PWD}\n")
   fout.write("source scripts/setup.sh\n")
   #fout.write(submitCommandPulseShape+"\n")
   fout.write(parallelPulseShape+"\n")
   fout.write("echo 'STOP---------------'\n")
   fout.write("echo\n")
   fout.write("echo\n")
os.system("chmod 755 "+jobsFolder+"/"+jobs_drawPS)
      
if args.submit:
   print("++++ SUBMIT MODULE CHARACTERIZATION  ++++")
   os.system("source "+jobsFolder+"/"+jobs_modChar)
   print("++++ SUBMIT PULSE SHAPE ++++")
   os.system("source "+jobsFolder+"/"+jobs_drawPS)

   

#submit the drawLaser_vs_SR summary plot
if args.submit:
   print("++++ CREATE SUMMARY PLOT ++++")
   os.system("python macros/drawLaserCTR_vs_SR_fede.py --pngName "+pngName+" --rangeFile "+args.rangeFile)
