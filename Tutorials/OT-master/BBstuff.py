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
from mpa_configurations import *
import os.path
import glob
import re
import ROOT



def loadValuesFromCSV(csvfilename):
    #print(csvfilename)
    #valuedict = dict()
    values = []
    with open(csvfilename, 'r') as f:
        reader = csv.reader(f, delimiter=',')
        for row in reader:
            if row[0] == '':
                continue
            pixedid = int(row[0])
            value = float(row[1])
            #valuedict[pixedid] = value
            values.append(value)
    #return valuedict
    return values

def AnalyzeBBallchips(modulename, show_plot=True, save_plot=True):
    moduleid = "mpa_test_"+modulename
    thepath = modulename+"/"
    if not os.path.isdir(thepath):
        print("The directory  "+thepath+"  does not exist - cannot plot module maps.")
        return
    chipnames  = ['','','','','','','','','','','','','','','','']
    for i in range(0,16):
        teststring = moduleid + "_Chip" + str(i+1)+"_"
        filelist = glob.glob(thepath+teststring+"*.csv")
        #filelist = glob.glob(thepath+teststring+"*_refitted.csv")
        maxstr = "" #string with highest number
        maxnbr = 0
        for f in filelist:
            numbers_from_string = filter(str.isdigit, f)
            if numbers_from_string > maxnbr:
                maxnbr = numbers_from_string
                maxstr = f
        #print(maxstr)

        reducedmaxstr = maxstr[maxstr.find('Chip'):-4]#ensure that I don't make a mistake given that AEMTec have one number, and HPK have 2 numbers
        #print(reducedmaxstr)
        numberlist = re.findall('\d+', reducedmaxstr )
        if len(numberlist)>=7:#1 chip + 6 date
            chipstring = "Chip"+numberlist[0]
            for j in range(1,7):
                chipstring += "_" + numberlist[j]
            #print(chipstring)
            chipnames[i] = chipstring
        else:
            print("Could not find any csv file for chip "+ str(i+1)+" for the module "+modulename)
            return
    #print(chipnames)
    for chipname in chipnames:
        #AnalyzeBBonechip(inpath=thepath,mapsaname=modulename,chip=chipname)
        AnalyzeBBonechipv2(inpath=thepath,mapsaname=modulename,chip=chipname)

def AnalyzeBBonechip(inpath,mapsaname,chip):
    moduleid = "mpa_test_"+mapsaname
    thepath = mapsaname+"/"
    if not os.path.isdir(thepath):
        print("The directory  "+thepath+"  does not exist - cannot plot module maps.")
        return
    #bbscurvecsv  = thepath + moduleid + "_" + chip + "_BumpBonding_SCurve_BadBump.csv"
    #bbscurvecsv2 = thepath + moduleid + "_" + chip + "_BumpBondingv2_SCurve_BadBump.csv"
    bbnoisecsv   = thepath + moduleid + "_" + chip + "_BumpBonding_Noise_BadBump_refitted.csv"
    bbnoisecsv2  = thepath + moduleid + "_" + chip + "_BumpBondingv2_Noise_BadBump_refitted.csv"

    
    hist = ROOT.TH1F("h_noise_"+chip, "", 2550, 0, 255)
    histv2 = ROOT.TH1F("h_noisev2_"+chip, "", 2550, 0, 255)
    if(os.path.isfile(bbnoisecsv)):
        bbnoise = loadValuesFromCSV(bbnoisecsv)
        for n in bbnoise:
            if n < 0: hist.Fill(0.0000001)
            elif n>255: hist.Fill(254.9999999)
            else: hist.Fill(n)
    if(os.path.isfile(bbnoisecsv2)):
        bbnoisev2 = loadValuesFromCSV(bbnoisecsv2)
        for n in bbnoisev2:
            if n < 0: histv2.Fill(0.0000001)
            elif n>255: histv2.Fill(254.9999999)
            else: histv2.Fill(n)
    outfile = ROOT.TFile(thepath + "BBstudy_"+mapsaname+".root","update")
    outfile.cd()
    hist.Write(hist.GetName(),ROOT.TH1F.kOverwrite);
    histv2.Write(histv2.GetName(),ROOT.TH1F.kOverwrite);
    outfile.Close()

    
def AnalyzeBBonechipv2(inpath,mapsaname,chip):
    moduleid = "mpa_test_"+mapsaname
    thepath = mapsaname+"/"
    if not os.path.isdir(thepath):
        print("The directory  "+thepath+"  does not exist - cannot plot module maps.")
        return
    #bbscurvecsv  = thepath + moduleid + "_" + chip + "_BumpBonding_SCurve_BadBump.csv"
    #bbscurvecsv2 = thepath + moduleid + "_" + chip + "_BumpBondingv2_SCurve_BadBump.csv"
    bbnoisecsv   = thepath + moduleid + "_" + chip + "_BumpBonding_Noise_BadBump_refitted.csv"
    bbnoisecsv2  = thepath + moduleid + "_" + chip + "_BumpBondingv2_Noise_BadBump_refitted.csv"

    chipshortname = chip[chip.find('Chip'):-20]
    
    hist = ROOT.TH1F("h_noise_"+chipshortname, "", 2560, 0, 256)
    histv2 = ROOT.TH1F("h_noisev2_"+chipshortname, "", 2560, 0, 256)
    if(os.path.isfile(bbnoisecsv)):
        bbnoise = loadValuesFromCSV(bbnoisecsv)
        for n in bbnoise:
            if n < 0: hist.Fill(0.0000001)
            elif n>255: hist.Fill(254.9999999)
            else: hist.Fill(n)
    if(os.path.isfile(bbnoisecsv2)):
        bbnoisev2 = loadValuesFromCSV(bbnoisecsv2)
        for n in bbnoisev2:
            if n < 0: histv2.Fill(0.0000001)
            elif n>255: histv2.Fill(254.9999999)
            else: histv2.Fill(n)
    outfile = ROOT.TFile("BBnoisefiles/" + "BBstudy_"+mapsaname+".root","update")
    outfile.cd()
    hist.Write(hist.GetName(),ROOT.TH1F.kOverwrite);
    histv2.Write(histv2.GetName(),ROOT.TH1F.kOverwrite);
    outfile.Close()


allmapsas = ["AEMTec1","AEMTec2","AEMTec3","AEMTec4","AEMTec5","AEMTec6","AEMTec7","AEMTec8","AEMTec9","AEMTec10","HPK2_1","HPK15_1","HPK15_2","HPK16_1","HPK17_1","HPK18_1","HPK19_1","HPK25_2","HPK26_1","HPK26_2"]
for mapsa in allmapsas:
    AnalyzeBBallchips(mapsa)
