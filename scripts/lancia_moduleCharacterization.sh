##########################################################
### 22Na runs

###firstRun=6801
###lastRun=6801
###label='run6801'
###
####step1, ith1
###./submit_moduleCharacterization_TOFHIR2.py -l test -b /home/cmsdaq/Programs/Lab5015Analysis/ -e bin/moduleCharacterization_step1.exe -r $firstRun-$lastRun -c cfg/moduleCharacterization_22Na.cfg  --vth vth1 -i /data/Lab5015Analysis/moduleCharacterization/TOFHIR2X/ -o /data/Lab5015Analysis/moduleCharacterization/TOFHIR2X/ --submit -p 30
###
###./submit_moduleCharacterization_TOFHIR2.py -l test -b /home/cmsdaq/Programs/Lab5015Analysis/ -e bin/moduleCharacterization_step2.exe -r $firstRun-$lastRun  -c cfg/moduleCharacterization_22Na.cfg  --vth vth1 -i /data/Lab5015Analysis/moduleCharacterization/TOFHIR2X/ -o /data/Lab5015Analysis/moduleCharacterization/TOFHIR2X/ --submit -p 24
###
###plotFolder=/var/www/html/TOFHIR2X/ModuleCharacterization/$label/
###
###command="mkdir $plotFolder"
###eval $command
###
###command="python ../macros/moduleCharacterizationSummaryPlots.py -m2 -o $plotFolder -i "
###for run in `seq $firstRun $lastRun`; do
###    command+="run$run,"
###done;
###command2=${command::-1}
###eval $command2



##########################################################
### laser runs

#step1, ith1
#./submit_moduleCharacterization_TOFHIR2.py -l test -b /home/cmsdaq/Programs/Lab5015Analysis/ -e bin/moduleCharacterization_step1.exe -r 8276-8287 -c cfg/moduleCharacterization.cfg  --vth vth1 -i /data/Lab5015Analysis/moduleCharacterization/TOFHIR2X/ -o /data/Lab5015Analysis/moduleCharacterization/TOFHIR2X/ --submit -p 30
./submit_moduleCharacterization_TOFHIR2.py -l test -b /home/cmsdaq/Programs/Lab5015Analysis/ -e bin/moduleCharacterization_step1.exe -r 8292 -c cfg/moduleCharacterization.cfg  --vth vth1 -i /data/Lab5015Analysis/moduleCharacterization/TOFHIR2X/ -o /data/Lab5015Analysis/moduleCharacterization/TOFHIR2X/ --submit -p 30

#step1, ith2
#./submit_moduleCharacterization_TOFHIR2.py -l test -b /home/cmsdaq/Programs/Lab5015Analysis/ -e bin/moduleCharacterization_step1.exe -r 2468,2470,2472,2474,2476,2478,2480,2482 -c cfg/moduleCharacterization.cfg  --vth vth2 -i /data/Lab5015Analysis/moduleCharacterization/TOFHIR2X/ -o /data/Lab5015Analysis/moduleCharacterization/TOFHIR2X/ --submit -p 24

#step2
#./submit_moduleCharacterization_TOFHIR2.py -l test -b /home/cmsdaq/Programs/Lab5015Analysis/ -e bin/moduleCharacterization_step2.exe -r 8276-8287  -c cfg/moduleCharacterization.cfg  --vth vth1 -i /data/Lab5015Analysis/moduleCharacterization/TOFHIR2X/ -o /data/Lab5015Analysis/moduleCharacterization/TOFHIR2X/ --submit -p 24
./submit_moduleCharacterization_TOFHIR2.py -l test -b /home/cmsdaq/Programs/Lab5015Analysis/ -e bin/moduleCharacterization_step2.exe -r 8292  -c cfg/moduleCharacterization.cfg  --vth vth1 -i /data/Lab5015Analysis/moduleCharacterization/TOFHIR2X/ -o /data/Lab5015Analysis/moduleCharacterization/TOFHIR2X/ --submit -p 24
