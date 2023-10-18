#! /usr/bin/env python
import math

def PDE(sipm_type, ov):

    if   (sipm_type == "HPK-MS"): 
        return 0.389 * ( 1. - math.exp(-1.*0.593*ov) )
    elif (sipm_type == "FBK-MS"): 
        return 0.419907 * ( 1. - math.exp(-1.*0.3046*ov) )
    elif (sipm_type == "FBK-W4S"): 
        return 0.411 * ( 1. - math.exp(-1.*0.191*ov) )/1.071
    elif (sipm_type == "FBK-W4C"): 
        return 0.490 * ( 1. - math.exp(-1.*0.225*ov) )/1.071

    elif (sipm_type == "HPK-PIT-C20-ES2"):
        return 0.576 * ( 1. - math.exp(-1.*0.625*ov) )
    elif (sipm_type == "HPK-PIT-C25-ES2"):
        return 0.638 * ( 1. - math.exp(-1.*0.651*ov) )
    elif (sipm_type == "HPK-PIT-C20-ES3"):
        return 0.568 * ( 1. - math.exp(-1.*0.588*ov) )
    elif (sipm_type == "HPK-PIT-C25-ES3"):
        return 0.638 * ( 1. - math.exp(-1.*0.589*ov) )
    elif (sipm_type == "HPK-PIT-C30-ES3"):
        return 0.653 * ( 1. - math.exp(-1.*0.728*ov) )
        
    else:
        print("error!!! PDE not specified!")
        #return 0.389 * ( 1. - math.exp(-1.*0.593*ov) )
        return 0


def Gain(sipm_type, ov):
    if   (sipm_type == "HPK-MS"): 
        return 97602.9*(ov+0.377962)
    elif (sipm_type == "FBK-MS"): 
        return 94954.6*(ov+0.512167)
    elif (sipm_type == "FBK-W4S"): 
        return 91541.7*(ov+0.408182)
    elif (sipm_type == "FBK-W4C"): 
        return 91541.7*(ov+0.408182)
    elif   (sipm_type == "HPK-PIT-C20-ES2"): 
        return 6.234E04 + 1.787E05*ov
    elif   (sipm_type == "HPK-PIT-C25-ES2"): 
        return 7.044E04 + 2.895E05*ov
    elif   (sipm_type == "HPK-PIT-C20-ES3"): 
        return 5.731E04 + 1.759E05*ov
    elif   (sipm_type == "HPK-PIT-C25-ES3"): 
        return 7.857E04 + 2.836E05*ov
    elif   (sipm_type == "HPK-PIT-C30-ES3"): 
        return 9.067E04 + 4.020E05*ov
    elif   (sipm_type == "HPK-PIT-C30-ES2"): #BOH 
        return 9.067E04 + 4.020E05*ov

    else:
        print("error!!! Gain not specified!")
        return 1
        
    #        return 36890.225 + 97602.904*ov
    

        
def ECF(sipm_type, ov):
    if   (sipm_type == "HPK-MS"): 
        return 1 + 0.000790089*ov + 0.00226734*ov*ov
    elif (sipm_type == "FBK-MS"): 
        return 1 + 0.00215668*ov +0.00303006*ov*ov
    elif (sipm_type == "FBK-W4S"): 
        return 1 +  0.00215668*ov +0.00303006*ov*ov
    elif (sipm_type == "FBK-W4C"): 
        return 1 +  0.00215668*ov +0.00303006*ov*ov
    elif   ("HPK-PIT" in sipm_type): 
        return 1 + 0.000790089*ov + 0.00226734*ov*ov # place-holder
        
    else:
        print("error!!! ECF not specified!")
        return 100 
        
        #return 1. -2.60076e-02*ov + 9.10258e-03*ov*ov
        #return 1.




def GetSaturationCorrection(Ncells, Edep, LY, PDE, LCE):
    Npe    = Edep * LY * PDE * LCE
    Nfired = Ncells * (1 - math.exp(-Npe/Ncells))
    k = Nfired/Npe
    return k

def fit_PDE(x, par):
        xx = x[0]
        return par[0] * ( PDE(sipm_type, xx-par[1]) )

def fit_PDE_ECF_Gain(x, par):
        xx = x[0]
        return par[0] * ( PDE(sipm_type, xx-par[1]) * ECF(sipm_type, xx-par[1]) * Gain(sipm_type, xx-par[1]) )

