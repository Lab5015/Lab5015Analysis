#python sendJobsOnLxplus.py -l conf9.2_Caltech -b /afs/cern.ch/user/m/mtd/Lab5015Analysis/ -e bin/analyzeTOFPET2_new_ab.exe -o analyzeTOFPET1_conf9.2_CaltechMatrix_ab -O /eos/user/m/mtd/Lab5015Analysis/ -c cfg/analyzeTOFPET_TB_conf9.2_CaltechMatrix_ab.cfg -r 25000-25250,25270-25436 -n 10 -s 

#python sendJobsOnLxplus.py -l conf9.2_Milano -b /afs/cern.ch/user/m/mtd/Lab5015Analysis/ -e bin/analyzeTOFPET2_new_ab.exe -o analyzeTOFPET1_conf9.2_MilanoMatrix_ab -O /eos/user/m/mtd/Lab5015Analysis/ -c cfg/analyzeTOFPET_TB_conf9.2_MilanoMatrix_ab.cfg -r 25000-25250,25270-25436 -n 10 -s 

#python sendJobsOnLxplus.py -l conf9.2_Milano -b /afs/cern.ch/work/m/malberti/MTD/TBatFNALFeb2020/Lab5015Analysis/ -e bin/analyzeTOFPET2_step1.exe -o analyzeTOFPET_step1_conf9.2_MilanoMatrix_ab -O /afs/cern.ch/work/m/malberti/MTD/TBatFNALFeb2020/Lab5015Analysis/ -c cfg/analyzeTOFPET_TB_conf9.2_MilanoMatrix_ab.cfg -r 25000-25250,25270-25436 -n 10 -s 

python sendJobsOnLxplus.py -l conf4 -b /afs/cern.ch/work/m/malberti/MTD/TBatFNALFeb2020/Lab5015Analysis/ -e bin/analyzeTOFPET2_step1.exe -o analyzeTOFPET_step1_conf4 -O /afs/cern.ch/work/m/malberti/MTD/TBatFNALFeb2020/Lab5015Analysis/ -c cfg/analyzeTOFPET_TB_conf4_mm.cfg -r 21433-21875,21954-22008 -n 10 -s 

python sendJobsOnLxplus.py -l conf5 -b /afs/cern.ch/work/m/malberti/MTD/TBatFNALFeb2020/Lab5015Analysis/ -e bin/analyzeTOFPET2_step1.exe -o analyzeTOFPET_step1_conf5 -O /afs/cern.ch/work/m/malberti/MTD/TBatFNALFeb2020/Lab5015Analysis/ -c cfg/analyzeTOFPET_TB_conf5_mm.cfg -r 22086-22200 -n 10 -s 

python sendJobsOnLxplus.py -l conf7 -b /afs/cern.ch/work/m/malberti/MTD/TBatFNALFeb2020/Lab5015Analysis/ -e bin/analyzeTOFPET2_step1.exe -o analyzeTOFPET_step1_conf7 -O /afs/cern.ch/work/m/malberti/MTD/TBatFNALFeb2020/Lab5015Analysis/ -c cfg/analyzeTOFPET_TB_conf7_mm.cfg -r 23204-23585 -n 10 -s 

python sendJobsOnLxplus.py -l conf9.1_Milano -b /afs/cern.ch/work/m/malberti/MTD/TBatFNALFeb2020/Lab5015Analysis/ -e bin/analyzeTOFPET2_step1.exe -o analyzeTOFPET_step1_conf9.1_MilanoMatrix -O /afs/cern.ch/work/m/malberti/MTD/TBatFNALFeb2020/Lab5015Analysis/ -c cfg/analyzeTOFPET_TB_conf9.1_MilanoMatrix_mm.cfg -r 24926-25250 -n 10 -s 

python sendJobsOnLxplus.py -l conf9.1_Caltech -b /afs/cern.ch/work/m/malberti/MTD/TBatFNALFeb2020/Lab5015Analysis/ -e bin/analyzeTOFPET2_step1.exe -o analyzeTOFPET_step1_conf9.1_CaltechMatrix -O /afs/cern.ch/work/m/malberti/MTD/TBatFNALFeb2020/Lab5015Analysis/ -c cfg/analyzeTOFPET_TB_conf9.1_CaltechMatrix_mm.cfg -r 24926-25250 -n 10 -s 
