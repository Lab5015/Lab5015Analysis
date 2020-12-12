#!/usr/bin/env python
import os, re
import commands
import math, time
import sys
import argparse
import subprocess

parser = argparse.ArgumentParser(description='This script splits moduleCharacterization_step* tasks in multiple parallel jobs')

parser.add_argument("-l",  "--label",          required=True, type=str, help="job label")
parser.add_argument("-i",  "--inputFolder",     required=True, type=str, help="input file folder")
parser.add_argument("-o",  "--outputFolder",     required=True, type=str, help="output file folder")
parser.add_argument("-b",  "--baseFolder",     required=True, type=str, help="base folder")
parser.add_argument("-e",  "--exeName",        required=True, type=str, help="absolute path of executable")
parser.add_argument("-r",  "--runs",           required=True, type=str, help="comma-separated list of runs to be processed")
parser.add_argument("-c",  "--configFile",     required=True, type=str, help="config file")
parser.add_argument("-s",  "--submit",         													help="submit jobs", action='store_true')
parser.add_argument("-v",  "--verbose",        													help="increase output verbosity", action='store_true')
parser.add_argument("-p",  "--numParal", 			 required=True, type=int, help='number of parallel jobs, 0 for all')

args = parser.parse_args()

preruns = []

runs = []
comma_list = args.runs.split(',')
for item in comma_list:
	hyphen_list = item.split('-')
	if len(hyphen_list) > 1:
		for i in range(int(hyphen_list[0]), int(hyphen_list[1])+1):
			preruns.append(i)
	else:
		preruns.append(int(hyphen_list[0]))

if args.numParal == 0:
	listRuns= []
	for i in range (len(preruns)):
		listRuns.append(int(preruns[i]))
	runs.append(listRuns);

if args.numParal !=0:
	step = (int)(len(preruns)/args.numParal)
	
	for i in range (step): 
		listRuns= []
		for j in range ( args.numParal * i , args.numParal * (i+1)):
			listRuns.append(int(preruns[j]))
		runs.append(listRuns)
		del listRuns
	
	if len(preruns) % args.numParal != 0:
		listRuns= []
		for j in range ( step * args.numParal, len(preruns)):
			listRuns.append(int(preruns[j]))
		runs.append(listRuns)
		
			
print runs

print 
print 'START'
print 



for i in range (len(runs)):
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


	parallelCommand = "parallel --results ./jobs/"+args.label+"/ "+args.baseFolder+"/"+args.exeName+" ::: "


	##### creates directory and file list for job #######
	jobDir = args.baseFolder+'/jobs/'+args.label+'/'
	os.system('mkdir '+jobDir)
	os.chdir(jobDir)
	for run in runs[i]:
	##### creates config file #######
		with open(args.baseFolder + '/'+args.configFile) as fi:
			contents = fi.read()
			configFileName = jobDir+"/config_run"+str(run)+".cfg"
			with open(configFileName, "w") as fo:
				fo.write(contents)
			command = 'sed -i \"s%^runs .*$%runs '+str(run)+'%\" '+configFileName
			os.system(command)
			command = 'sed -i \"s%^step1FileName .*$%step1FileName '+args.inputFolder+'/plots/moduleCharacterization_step1_run'+str(run)+'.root%\" '+configFileName
			os.system(command)
			command = 'sed -i \"s%^outFileNameStep1 .*$%outFileNameStep1 '+args.outputFolder+'/plots/moduleCharacterization_step1_run'+str(run)+'.root%\" '+configFileName
			os.system(command)
			command = 'sed -i \"s%^outFileNameStep2 .*$%outFileNameStep2 '+args.outputFolder+'/plots/moduleCharacterization_step2_run'+str(run)+'.root%\" '+configFileName
			os.system(command)
			#command = 'sed -i \"s%^plotDir .*$%plotDir '+'/var/www/html/ModuleCharacterization/ALIO/Co60/FBK_thinQuartz/Data/'+str(run)+ '%\" '+configFileName
			#os.system(command)
      
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
	   os.system("chmod 755 job.sh")
   
	if args.submit:
  	 os.system("source "+jobDir+"/jobs.sh")
	print 'end job'
print
print 'END'
print
