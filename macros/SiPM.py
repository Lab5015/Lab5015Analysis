#! /usr/bin/env python
import math

def PDE(ov, sipm, irr='0'):
    k = 1.
    if (irr == '2E14' and 'HPK' in sipm): k = 0.85 # 15% PDE reduction for HPK SiPMs irradiated 2E14   
    if ('HPK' in sipm):
        return k * 1.0228 * 0.384 * ( 1. - math.exp(-1.*0.583*ov) ) # 1.0228 factor to account for LYSO emission spectrum
    if ('FBK' in sipm):
        return k * 0.8847*0.466 * ( 1. - math.exp(-1.*0.314*ov) ) # 0.8847 factor to account for LYSO emission spectrum

def Gain(ov, sipm, irr='0'):
    k = 1.
    #if (irr == '2E14'): k = 0.7 # gain reduction for 2E14 irradiated SiPMs 
    if ('HPK' in sipm):
        return k*(36890. + 97602.*ov) # HPK
    if ('FBK' in sipm):
        return k*(50739. + 95149.*ov) # FBK
    
    
def sigma_noise(sr):
    noise_single = math.sqrt( pow(420./sr,2) + 16.7*16.7 )
    return noise_single / math.sqrt(2)
                                                                                    
                            
