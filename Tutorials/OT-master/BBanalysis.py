import numpy as np
import time
import sys
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.special import erfc
from scipy.special import erf
import matplotlib.cm as cm
import matplotlib.backends.backend_pdf as pltpdf
#import seaborn as sns
import csv
import math
#sys.path.append('myScripts/') #python 2.7 ?
#from mpa_configurations import * # python 2.7 ?
#from mpa_configurations import *
import os.path
import glob
import re
import ROOT


ROOT.TH1.AddDirectory(ROOT.kFALSE);

def LoadHistogram(mapsa,chip,v2=False):
    h = ROOT.TH1F()
    htemp = ROOT.TH1F()
    #hfile = ROOT.TFile.Open("BBnoisefiles/" + "BBstudy_"+mapsa+".root")
    hfile = ROOT.TFile.Open("BBstudy_"+mapsa+".root")
    clonename = "noise_"+chip
    if v2:
        htemp = hfile.Get("h_noisev2_"+chip)
        clonename = "noisev2_"+chip
    else:
        #print ("h_noise_"+chip)
        htemp = hfile.Get("h_noise_"+chip)
    #hfile.Close()
    h = htemp.Clone(clonename)
    #print(type(h))
    return h

def CloneHistogram(hist):
    return hist.Clone(hist.GetName()+"_clone")

#TMath.Gaus(x,mean,sigma)
#TMath.Landau(x,mpv,sigma)
#TMath.BreitWigner(x,mean,sigma)
###TMath.Voigt(x,sigma,lg)
def FitHistogram(hist,mapsa):
    name = hist.GetName()
    hist.SetLineColor(ROOT.kBlack)
    h = CloneHistogram(hist)
    h.SetBinContent(0,0)#clear fit failures
    h.SetBinContent(1,0)#clear fit failures
    for i in range(1,h.GetNbinsX()+1):
        if h.GetBinContent(i)>=10 and 0.25*h.GetBinContent(i)>h.GetBinContent(i+1) and 0.25*h.GetBinContent(i)>h.GetBinContent(i-1):
            h.SetBinContent(h.GetMaximumBin(),max(h.GetBinContent(i-1),h.GetBinContent(i+1)))
    while (h.GetNbinsX()>128 and h.GetBinContent(h.GetMaximumBin())<30):
        #print(h.GetNbinsX())
        h.Rebin(2)
    #if h.GetBinCenter(h.GetMaximumBin())>50:
    #    h.Rebin(10)
    #elif h.GetBinCenter(h.GetMaximumBin())>20:
    #    h.Rebin(3)
    #h.Rebin(3)
    #first removed spike, then rebin
    #h.SetBinContent(1,0)
    #h.SetBinContent(h.GetNbinsX(),0)    
    #print("Maximum: ",h.GetMaximum(),h.GetMaximumBin())
    #remove non-scurve-fit spike
    #if 0.25*h.GetBinContent(h.GetMaximumBin())>h.GetBinContent(h.GetMaximumBin()+1) and 0.25*h.GetBinContent(h.GetMaximumBin())>h.GetBinContent(h.GetMaximumBin()+1):
    #    h.SetBinContent(h.GetMaximumBin(),0)
        #print("Mod")
    #print("Maximum: ",h.GetMaximum(),h.GetMaximumBin())    
    #print("Maximum: ",h.GetBinContent(h.GetMaximumBin()),h.GetMaximumBin())
    meanbin = h.GetMaximumBin()
    widthbinl, widthbinh = 20, 20
    for i in reversed(xrange(h.GetMaximumBin())):
        if h.GetBinContent(i) <= 0.5 * h.GetMaximum():
            widthbinl = h.GetMaximumBin() - i
            break
    for i in xrange(h.GetMaximumBin()+1,h.GetNbinsX()+1):
        if h.GetBinContent(i) <= 0.5 * h.GetMaximum():
            widthbinh = i-h.GetMaximumBin()
            break
    meanbinx = h.GetBinCenter(meanbin)
    widthx = max(widthbinl,widthbinh) * h.GetBinWidth(meanbin)
    if widthx<0: print("WTF ",widthx,widthbin, h.GetMaximumBin())
    #funcG = ROOT.TF1("Gauss_"+name, "[0]*ROOT.TMath.Gaus(x,[1],[2])",meanbinx-3.*widthx,meanbinx+3.*widthx)
    funcG = ROOT.TF1(name+"_Gauss", "[0]*TMath::Gaus(x,[1],[2])",meanbinx-5.*widthx,meanbinx+5.*widthx)
    #funcG = ROOT.TF1(name+"_Gauss", "[0]*TMath::Gaus(x,[1],[2])",0.,255.)
    funcG.SetLineWidth(4)
    funcG.SetLineColor(ROOT.kRed)
    funcG.SetParameters(h.GetMaximum(),meanbinx,widthx)
    funcG.SetParLimits(0,h.GetMaximum()/3.,h.GetMaximum()*3.)
    funcG.SetParLimits(1,max(0,meanbinx-3.*widthx),min(255.,meanbinx+3.*widthx))
    funcG.SetParLimits(2,0.,5.*widthx)
    h.Fit(funcG,"RQ+")
    """
    #funcL = ROOT.TF1(name+"_Landau", "[0]*TMath::Landau(x,[1],[2])",meanbinx-5.*widthx,meanbinx+5.*widthx)
    funcL = ROOT.TF1(name+"_Landau", "[0]*TMath::Landau(x,[1],[2])",meanbinx-5.*widthx,255.)#makes only sense with a tail
    #funcL = ROOT.TF1(name+"_Landau", "[0]*TMath::Landau(x,[1],[2])",0.,255.)
    funcL.SetLineWidth(4)
    funcL.SetLineColor(ROOT.kMagenta)
    funcL.SetParameters(h.GetMaximum(),meanbinx,widthx)
    #funcL.SetParLimits(0,h.GetMaximum()/3.,h.GetMaximum()*3.)
    funcL.SetParLimits(1,max(0,meanbinx-3.*widthx),min(255.,meanbinx+3.*widthx))
    #funcL.SetParLimits(2,0.,5.*widthx)
    h.Fit(funcL,"RQ+")
    """
    funcBW = ROOT.TF1(name+"_BreitWigner", "[0]*TMath::BreitWigner(x,[1],[2])",meanbinx-5.*widthx,meanbinx+5.*widthx)
    #funcBW = ROOT.TF1(name+"_BreitWigner", "[0]*TMath::BreitWigner(x,[1],[2])",0.,255.)
    funcBW.SetLineWidth(4)
    funcBW.SetLineColor(ROOT.kGreen+2)
    funcBW.SetParameters(h.GetMaximum(),meanbinx,widthx)
    funcBW.SetParLimits(0,h.GetMaximum()/5.,h.GetMaximum()*50.)
    funcBW.SetParLimits(1,max(0,meanbinx-5.*widthx),min(255.,meanbinx+5.*widthx))
    funcBW.SetParLimits(2,0.,5.*math.sqrt(math.pi)*widthx)
    h.Fit(funcBW,"RQ+")
    
    badpixelsG, badpixelsL,badpixelsBW = 0,0,0
    #for i in range(0,hist.GetNbinsX()):
    for i in range(2,hist.GetNbinsX()):
        if hist.GetBinCenter(i)< funcG.GetParameter(1)-5.* abs(funcG.GetParameter(2)):
            badpixelsG  += hist.GetBinContent(i)
        #if hist.GetBinCenter(i)< funcL.GetParameter(1)-5.* abs(funcL.GetParameter(2)):
        #    badpixelsL  += hist.GetBinContent(i)
        if hist.GetBinCenter(i)<funcBW.GetParameter(1)-5.*abs(funcBW.GetParameter(2)/math.sqrt(math.pi)): #divide by sqrt(pi) to approximate Gaussian sigma
            badpixelsBW += hist.GetBinContent(i)
    
    filetest = ROOT.TFile("FitTest_"+mapsa+".root","update")
    filetest.cd()
    hist.Write(hist.GetName(),ROOT.TObject.kOverwrite)
    h.Write(h.GetName(),ROOT.TObject.kOverwrite)
    funcG.Write(funcG.GetName(),ROOT.TObject.kOverwrite)
    #funcL.Write(funcL.GetName(),ROOT.TObject.kOverwrite)
    funcBW.Write(funcBW.GetName(),ROOT.TObject.kOverwrite)
    filetest.Close()

    #return badpixelsG, badpixelsL,badpixelsBW, funcG.GetParameter(1),abs(funcG.GetParameter(2)), funcL.GetParameter(1),abs(funcL.GetParameter(2)), funcBW.GetParameter(1),abs(funcBW.GetParameter(2))
    return badpixelsG,badpixelsBW, funcG.GetParameter(1),abs(funcG.GetParameter(2)), funcBW.GetParameter(1),abs(funcBW.GetParameter(2))
    


def runAnalysis(v2=False):
    chiplist = ["Chip1","Chip2","Chip3","Chip4","Chip5","Chip6","Chip7","Chip8","Chip9","Chip10","Chip11","Chip12","Chip13","Chip14","Chip15","Chip16"]
    allmapsas = ["AEMTec1","AEMTec2","AEMTec3","AEMTec4","AEMTec5","AEMTec6","AEMTec7","AEMTec8","AEMTec9","AEMTec10","HPK2_1","HPK15_1","HPK15_2","HPK16_1","HPK17_1","HPK18_1","HPK19_1","HPK25_2","HPK26_1","HPK26_2"]
    badpixellistG, badpixellistL, badpixellistBW  = [], [], []
    meanlistG, meanlistL, meanlistBW  = [], [], []
    widthlistG, widthlistL, widthlistBW  = [], [], []

    for mapsa in allmapsas:
        tempG, tempL, tempBW = [], [], []
        meanG, meanL, meanBW = [], [], []
        widthG, widthL, widthBW = [], [], []
        for chip in chiplist:
            ##G,L,BW, mG,wG, mL,wL, mBW,wBW = FitHistogram(LoadHistogram(mapsa,chip),mapsa)
            ###G,L,BW, mG,wG, mL,wL, mBW,wBW = FitHistogram(LoadHistogram(mapsa,chip),mapsa, True)
            #G,BW, mG,wG, mBW,wBW = FitHistogram(LoadHistogram(mapsa,chip),mapsa)
            G,BW, mG,wG, mBW,wBW = FitHistogram(LoadHistogram(mapsa,chip, v2),mapsa)
            tempG.append(G)
            #tempL.append(L)
            tempBW.append(BW)
            meanG.append(mG)
            #meanL.append(mL)
            meanBW.append(mBW)
            widthG.append(wG)
            #widthL.append(wL)
            widthBW.append(wBW)
        badpixellistG.append(tempG)
        #badpixellistL.append(tempL)
        badpixellistBW.append(tempBW)
        meanlistG.append(meanG)
        #meanlistL.append(meanL)
        meanlistBW.append(meanBW)
        widthlistG.append(widthG)
        #widthlistL.append(widthL)
        widthlistBW.append(widthBW)

    bordervalues = []

    for i in range(len(allmapsas)):
        for j in range(len(chiplist)):
            #print("For MaPSA "+allmapsas[i]+" and "+chiplist[j]+" I find bad bumps: Gaussian = "+str(badpixellistG[i][j])+" Landau = "+str(badpixellistL[i][j]) + " BreitWigner = "+str(badpixellistBW[i][j]))
            #print("          : Gaussian = "+str(round(meanlistG[i][j],2))+"+/-"+str(round(widthlistG[i][j],2))+" Landau = "+str(round(meanlistL[i][j],2))+"+/-"+str(round(widthlistL[i][j],2)) + " BreitWigner = "+str(round(meanlistBW[i][j],2))+"+/-"+str(round(widthlistBW[i][j],2))) + " ("+str(round(widthlistBW[i][j]/math.sqrt(math.pi),2))+")"
            print("For MaPSA "+allmapsas[i]+" and "+chiplist[j]+" I find bad bumps: Gaussian = "+str(int(badpixellistG[i][j])) + ", BreitWigner = "+str(int(badpixellistBW[i][j])) + \
            "          - Gaussian = "+str(round(meanlistG[i][j],2))+"+/-"+str(round(widthlistG[i][j],2)) + ",   BreitWigner = "+str(round(meanlistBW[i][j],2))+"+/-"+str(round(widthlistBW[i][j],2))) + " ("+str(round(widthlistBW[i][j]/math.sqrt(math.pi),2))+")"
            bordervalues.append([allmapsas[i],chiplist[j],meanlistG[i][j],widthlistG[i][j]])

    csvfilename = "../GaussianFitValues.csv"
    if v2: csvfilename = "../GaussianFitValuesv2.csv"
    with open(csvfilename, "w") as f:
        writer = csv.writer(f)
        writer.writerows(bordervalues)

runAnalysis(True)
runAnalysis(False)
