#!/usr/bin/env python
import os, re
import commands
import sys


runs = { 2289 : 90,
         2290 : 88,
         2291 : 89,
         2292 : 87,
         2293 : 86,
         2294 : 85,
         2295 : 83,
         2296 : 80,
         2297 : 70,
         2298 : 60,
         2299 : 50,
         2300 : 40,
         2301 : 30,
         2302 : 20,
         2303 : 10
}



for run, tune in runs.items():
    print run, tune
    command1 = 'cp  analyzeTOFPET2_wirelessBar_HDR2_UVlaser_TEMPLATE.cfg analyzeTOFPET2_wirelessBar_HDR2_UVlaser_tune%d.cfg'%(tune)
    command2 = 'sed -i -- "s/MYRUN/%d/g" analyzeTOFPET2_wirelessBar_HDR2_UVlaser_tune%d.cfg'%(run,tune)
    command3 = 'sed -i -- "s/MYTUNE/%d/g" analyzeTOFPET2_wirelessBar_HDR2_UVlaser_tune%d.cfg'%(tune,tune)
    print command1
    print command2
    print command3
    os.system(command1)
    os.system(command2)
    os.system(command3)

print 'Done!'
