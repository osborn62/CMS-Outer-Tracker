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
sys.path.append('MPA_Test/myScripts/') #python 2.7 ?
from mpa_configurations import * # python 2.7 ?
#from MPA_Test.myScripts.mpa_configurations import *
from ArrayToCSV import *
import os.path
import glob
import re
import ROOT

conf = mpa_configurations()

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

def loadSCurvesFromCSV(csvfilename):
    #print(csvfilename)
    #valuedict = dict()
    values = []
    with open(csvfilename, 'r') as f:
        reader = csv.reader(f, delimiter=',')
        for row in reader:
            if row[0] == '':
                continue
            pixedid = int(row[0])
            #valuedict[pixedid] = value
            #values.append(value)
            value = [float(i) for i in row]
            value.pop(0)
            values.append(value)
    #return valuedict
    return values


def AnalyzeBBonechip(mapsaname,chip,v2=False):
    moduleid = "mpa_test_"+mapsaname
    thepath = mapsaname+"/"
    if not os.path.isdir(thepath):
        print("The directory  "+thepath+"  does not exist - cannot plot module maps.")
        return
    bbnoisecsv   = thepath + moduleid + "_" + chip + "_BumpBonding_Noise_BadBump_refitted.csv"
    gaussianCSVfilename = "GaussianFitValues.csv"
    if v2:
        bbnoisecsv   = thepath + moduleid + "_" + chip + "_BumpBondingv2_Noise_BadBump_refitted.csv"
        gaussianCSVfilename = "GaussianFitValuesv2.csv"
        
    gaussianCSV   = open(gaussianCSVfilename)
    gaussiandata  = csv.reader(gaussianCSV)
    gaussianarray = list(gaussiandata)

    #print(gaussianarray)
    extractedvalue = 3.#current default is vcal<=3
    for entry in gaussianarray:
        if len(entry)!=4: print("WTF: ",entry)
        if entry[0]!=mapsaname or entry[1]!=chip:
            continue
        #print(entry)
        extractedvalue = float(entry[2])-5.*float(entry[3])
        break
    RedoBBCutoffOneChip(extractedvalue, mapsaname, chip,v2)

    
def RedoBBCutoffOneChip(extractedvalue,mapsa, chip,v2=False):
    moduleid = "mpa_test_"+mapsa
    thepath = mapsa+"/"
    teststring = moduleid + "_"+chip+"_"
    bbstring = "_BumpBonding_Noise_BadBump_refitted"
    if v2: bbstring = "_BumpBondingv2_Noise_BadBump_refitted"
    filelist = glob.glob(thepath+teststring+"*"+bbstring+".csv")
    maxstr = "" #string with highest number
    maxnbr = 0
    #print (thepath+teststring+"*"+bbstring+".csv")
    if len(filelist) == 0:
        print("No file found: "+thepath+teststring+"*"+bbstring+".csv")
        return False
    for f in filelist:
        numbers_from_string = filter(str.isdigit, f)
        if numbers_from_string > maxnbr:
            maxnbr = numbers_from_string
            maxstr = f
    #print(maxstr)
    #print(maxstr[0:-4])
    csvfilename = maxstr
    newcsvfilename = csvfilename
    newcsvfilename = newcsvfilename.replace("Noise_BadBump","BadBumpGaussian")
    #print(csvfilename,len(filelist))
    data_array = loadSCurvesFromCSV(csvfilename)
    #data_array.insert(0, ["",0])
    newdata_array = []
    for i in range(len(data_array)):
        #print(i,data_array[i])
        if data_array[i][0]=="":
            newdata_array.append(["",""])
        else:
            if extractedvalue<=0:
                newdata_array.append([i,-3])
            elif data_array[i][0]<0:
                newdata_array.append([i,-2])
            elif data_array[i][0]==0:
                newdata_array.append([i,-1])
            elif data_array[i][0]<=3:
                newdata_array.append([i,1])
            else:
                newdata_array.append([i,0])

    newdata_array.insert(0, ["",0])
    #print(newcsvfilename)
    #with open("out.csv", "w", newline="") as f:
    with open(newcsvfilename, "w") as f:
        writer = csv.writer(f)
        writer.writerows(newdata_array)
    
def RedoBBGaussianCutoffs():
    chiplist = ["Chip1","Chip2","Chip3","Chip4","Chip5","Chip6","Chip7","Chip8","Chip9","Chip10","Chip11","Chip12","Chip13","Chip14","Chip15","Chip16"]
    #allmapsas = ["AEMTec1","AEMTec2","AEMTec3","AEMTec4","AEMTec5","AEMTec6","AEMTec7","AEMTec8","AEMTec9","AEMTec10","HPK2_1","HPK15_1","HPK15_2","HPK16_1","HPK17_1","HPK18_1","HPK19_1","HPK25_2","HPK26_1","HPK26_2"]
    allmapsas = ["AEMTec1","AEMTec2","AEMTec3","AEMTec4","AEMTec5","AEMTec6","AEMTec7","AEMTec8","AEMTec9","AEMTec10","HPK2_1","HPK15_1","HPK15_2","HPK16_1","HPK17_1","HPK18_1","HPK19_1","HPK25_2","HPK26_1","HPK26_2"]
    for mapsa in allmapsas:
        for chip in chiplist:
            print("Redo Bump Bonding Gaussian for " + chip + " of MaPSA " + mapsa)
            AnalyzeBBonechip(mapsa,chip,False)
            print("Redo Bump Bonding v2 Gaussian for " + chip + " of MaPSA " + mapsa)
            AnalyzeBBonechip(mapsa,chip,True)

        
RedoBBGaussianCutoffs()

