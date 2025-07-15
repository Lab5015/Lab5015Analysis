#!/bin/bash
#!/bin/sh
echo
echo 'START---------------'
echo 'current dir: ' ${PWD}
cd /afs/cern.ch/user/f/ftonetto/francesco/Lab5015Analysis/
echo 'current dir: ' ${PWD}
source scripts/setup.sh
./bin/moduleCharacterization_step1.exe $1
#./bin/moduleCharacterization_step2ES.exe $1
echo 'STOP---------------'
echo
echo
