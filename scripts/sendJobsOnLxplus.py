#!/usr/bin/env python
import os, re
import commands
import math, time
import sys
import argparse
import subprocess

from itertools import chain



def parse_range(rng):
   parts = rng.split('-')
   if 1 > len(parts) > 2:
      raise ValueError("Bad range: '%s'" % (rng,))
   parts = [int(i) for i in parts]
   start = parts[0]
   end = start if len(parts) == 1 else parts[1]
   if start > end:
      end, start = start, end
   return range(start, end + 1)



def parse_range_list(rngs):
   return sorted(set(chain(*[parse_range(rng) for rng in rngs.split(',')])))



# ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----

if __name__ == '__main__':
   
   parser = argparse.ArgumentParser(description='This script splits tasks in multiple jobs and sends them on htcondor')
   parser.add_argument("-l",  "--label",          required=True,         type=str,   help="job label")
   parser.add_argument("-b",  "--baseFolder",     required=True,         type=str,   help="base folder")
   parser.add_argument("-e",  "--exeName",        required=True,         type=str,   help="absolute path of executable")
   parser.add_argument("-o",  "--outputFileName", required=True,         type=str,   help="outputFileName")
   parser.add_argument("-O",  "--outputFolder",   required=True,         type=str,   help="outputFolder")
   parser.add_argument("-c",  "--configFile",     required=True,         type=str,   help="config file")
   parser.add_argument("-q",  "--queue",          default="microcentury",type=str,   help="htcondor queue to use")
   parser.add_argument("-s",  "--submit",                                            help="submit jobs", action='store_true')
   parser.add_argument("-v",  "--verbose",                                           help="increase output verbosity", action='store_true')
   parser.add_argument("-r",  "--runs",           required=True,          type=str,  help="run ranges")
   parser.add_argument("-n",  "--nRunsPerJob",    required=True,          type=int,  help="number of runs per job")
   args = parser.parse_args()
   
   
   print 
   print 'START'
   print 
   
   currDir = os.getcwd()
   
   if not os.path.exists('jobs/'+args.label):
      os.makedirs('jobs/'+args.label)
      
   ##### creates directory and file list for job #######
   jobDir = currDir+'/jobs/'+args.label+'/'
   os.chdir(jobDir)
   
      
   runList = parse_range_list(args.runs)
   nRuns = len(runList)
   print runList
   print "number of runs to process: " + str(nRuns)
   
   nJobs = int(math.ceil(1.*nRuns/args.nRunsPerJob))
   print "number of jobs: " + str(nJobs)
   
   for iJob in range(1,nJobs+1):
      
      print "\ndoing job "+str(iJob)+" / "+str(nJobs)
      
      runs = ""
      for jj in range(0,args.nRunsPerJob):
         if (jj+(iJob-1)*args.nRunsPerJob) < len(runList):
            runs += str(runList[jj+(iJob-1)*args.nRunsPerJob]) + ","
      runs = runs[:-1]
      
      ##### creates config file #######
      with open(args.baseFolder+'/'+args.configFile) as fi:
         contents = fi.read()
         configFileName = "config_"+str(iJob)+".cfg"
         with open(jobDir+"/"+configFileName, "w") as fo:
            fo.write(contents)
         command = 'sed -i \"s%^runs .*$%runs '+runs+'%\" '+jobDir+'/'+configFileName
         os.system(command)
         command = 'sed -i \"s%^outFileName .*$%outFileName '+args.outputFileName+'_'+str(iJob)+'%\" '+jobDir+'/'+configFileName
         os.system(command)
         
         ##### creates jobs #######
         with open('job_'+str(iJob)+'.sh', 'w') as fout:
            fout.write("#!/bin/sh\n")
            fout.write("echo\n")
            fout.write("echo 'START---------------'\n")
            fout.write("echo 'current dir: ' ${PWD}\n")
            fout.write("cd "+str(args.baseFolder)+"\n")
            fout.write("echo 'current dir: ' ${PWD}\n")
            fout.write("source scripts/setup.sh\n")
            fout.write("cd -\n")
            fout.write(args.baseFolder+"/"+args.exeName+" "+jobDir+"/"+configFileName+"\n")
            fout.write("mv "+args.outputFileName+"* "+args.outputFolder+"/\n")
            fout.write("echo 'STOP---------------'\n")
            fout.write("echo\n")
            fout.write("echo\n")
            fout.close()
            os.system("chmod 755 job_"+str(iJob)+".sh")
            
         #### creates htcondor job ######
         with open('htcondor_job_'+str(iJob)+'.sub', 'w') as fout:
            fout.write("executable            = %s/job_%s.sh\n"%(jobDir,str(iJob)))
            fout.write("arguments             = $(ClusterId) $(ProcId)\n")
            fout.write("output                = %s/job_%s.$(ClusterId).$(ProcId).out\n"%(jobDir,str(iJob)) )
            fout.write("error                 = %s/job_%s.$(ClusterId).$(ProcId).err\n"%(jobDir,str(iJob)) )
            fout.write("log                   = %s/job_%s.$(ClusterId).log\n"%(jobDir,str(iJob)) )
            fout.write("transfer_output_files = \"\"\n" )
            fout.write("+JobFlavour           = \"%s\"\n"%(args.queue) )
            fout.write("queue \n")
            fout.close()
         
         ###### sends bjobs ######
         if args.submit:
            command = "condor_submit "+jobDir+"/htcondor_job_"+str(iJob)+".sub"
            print command
            os.system(command)
            print "job " + str(args.label) + " submitted"
            
   os.chdir("../..")
            
   print
   print "your jobs:"
   os.system("condor_q")
   print
   print 'END'
   print
