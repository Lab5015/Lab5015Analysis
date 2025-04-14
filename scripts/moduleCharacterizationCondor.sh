#!/bin/bash
#!/bin/sh
echo
echo 'START---------------'
echo 'current dir: ' ${PWD}
cd /eos/home-f/ftonetto/Lab5015Analysis/
echo 'current dir: ' ${PWD}
source scripts/setup.sh
./bin/moduleCharacterization_step1.exe $1
./bin/moduleCharacterization_step2.exe $1
echo 'STOP---------------'
echo
echo
