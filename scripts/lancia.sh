#./submit_moduleCharacterization_TOFHIR2.py -l test -b /home/cmsdaq/Lab5015Analysis/ -e bin/moduleCharacterization_step1.exe -r 1073 -c cfg/moduleCharacterization.cfg  --vth vth1 -i /data/Lab5015Analysis/moduleCharacterization/ -o /data/Lab5015Analysis/moduleCharacterization/ --submit -p 6

#./submit_moduleCharacterization_TOFHIR2.py -l test -b /home/cmsdaq/Lab5015Analysis/ -e bin/moduleCharacterization_step2.exe -r 1072-1073 -c cfg/moduleCharacterization.cfg -i /data/Lab5015Analysis/moduleCharacterization/ -o /data/Lab5015Analysis/moduleCharacterization/ --submit -p 6

./submit_moduleCharacterization_TOFHIR2.py -l test -b /home/cmsdaq/Lab5015Analysis/ -e bin/moduleCharacterization_step3.exe -r 1072-1073 -c cfg/moduleCharacterization_step3.cfg -i /data/Lab5015Analysis/moduleCharacterization/ -o /data/Lab5015Analysis/moduleCharacterization/ --submit -p 6
