python scripts/submit_moduleCharacterization.py --label modChar1 -b /home/petsys/TOFHiR2X/sw_analysis/Lab5015Analysis -e bin/moduleCharacterization_step1.exe -r 98-105 -c cfg/moduleCharacterization_fede.cfg -j 2 --submit
python scripts/submit_moduleCharacterization.py --label modChar2 -b /home/petsys/TOFHiR2X/sw_analysis/Lab5015Analysis -e bin/moduleCharacterization_step2.exe -r 98-105 -c cfg/moduleCharacterization_fede.cfg -j 2 --submit

python scripts/submit_drawPulseShape.py --label pulseShape -b /home/petsys/TOFHiR2X/sw_analysis/Lab5015Analysis -e bin/drawPulseShape_fede.exe -r 98-105 -c cfg/drawPulseShape_fede.cfg -j 2 --submit
