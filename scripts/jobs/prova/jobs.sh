#!/bin/sh
echo
echo 'START---------------'
echo 'current dir: ' ${PWD}
cd /home/cmsdaq/Lab5015Analysis/
echo 'current dir: ' ${PWD}
source scripts/setup.sh
parallel --results ./jobs/prova/ /home/cmsdaq/Lab5015Analysis//bin/moduleCharacterization_step1.exe ::: /home/cmsdaq/Lab5015Analysis/scripts/jobs/prova//config_run2428.cfg /home/cmsdaq/Lab5015Analysis/scripts/jobs/prova//config_run2429.cfg /home/cmsdaq/Lab5015Analysis/scripts/jobs/prova//config_run2430.cfg /home/cmsdaq/Lab5015Analysis/scripts/jobs/prova//config_run2431.cfg /home/cmsdaq/Lab5015Analysis/scripts/jobs/prova//config_run2432.cfg /home/cmsdaq/Lab5015Analysis/scripts/jobs/prova//config_run2433.cfg /home/cmsdaq/Lab5015Analysis/scripts/jobs/prova//config_run2434.cfg /home/cmsdaq/Lab5015Analysis/scripts/jobs/prova//config_run2435.cfg /home/cmsdaq/Lab5015Analysis/scripts/jobs/prova//config_run2436.cfg /home/cmsdaq/Lab5015Analysis/scripts/jobs/prova//config_run2437.cfg /home/cmsdaq/Lab5015Analysis/scripts/jobs/prova//config_run2438.cfg /home/cmsdaq/Lab5015Analysis/scripts/jobs/prova//config_run2439.cfg /home/cmsdaq/Lab5015Analysis/scripts/jobs/prova//config_run2440.cfg /home/cmsdaq/Lab5015Analysis/scripts/jobs/prova//config_run2441.cfg /home/cmsdaq/Lab5015Analysis/scripts/jobs/prova//config_run2442.cfg /home/cmsdaq/Lab5015Analysis/scripts/jobs/prova//config_run2443.cfg /home/cmsdaq/Lab5015Analysis/scripts/jobs/prova//config_run2444.cfg /home/cmsdaq/Lab5015Analysis/scripts/jobs/prova//config_run2445.cfg /home/cmsdaq/Lab5015Analysis/scripts/jobs/prova//config_run2446.cfg /home/cmsdaq/Lab5015Analysis/scripts/jobs/prova//config_run2447.cfg /home/cmsdaq/Lab5015Analysis/scripts/jobs/prova//config_run2448.cfg /home/cmsdaq/Lab5015Analysis/scripts/jobs/prova//config_run2449.cfg /home/cmsdaq/Lab5015Analysis/scripts/jobs/prova//config_run2450.cfg /home/cmsdaq/Lab5015Analysis/scripts/jobs/prova//config_run2451.cfg /home/cmsdaq/Lab5015Analysis/scripts/jobs/prova//config_run2452.cfg /home/cmsdaq/Lab5015Analysis/scripts/jobs/prova//config_run2453.cfg 
echo 'STOP---------------'
echo
echo
