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

class mpa_cal_utility():
    def __init__(self):
        self.conf = conf

    def errorf(self,x, *p):
        a, mu, sigma, offset = p
        return 0.5*a*(1.0+erf((x-mu)/sigma)) + offset #XXX added offset
    def line(self,x, *p):
        g, offset = p
        return  np.array(x) *g + offset
    def gauss(self,x, *p):
        A, mu, sigma = p
        return A*np.exp(-(x-mu)**2/(2.*sigma**2))
    def errorfc(self,x, *p):
        a, mu, sigma, offset = p
        return a*0.5*erfc((x-mu)/sigma) + offset #XXX added offset

##### plot_extract_scurve function take scurve data and extract threhsold and noise data. If plot = 1, it also plot scurves and histograms
    def plot_extract_scurve(self,row, col, s_type, scurve, n_pulse, nominal_DAC, start, stop, extract, plot, save_plot, israw=False, printout=True, filename="../Results_MPATesting/plot_extract_scurve"):
        if len(scurve.shape)!=2:
            print("Expect a 2D scurve array")
            return
        if israw:
            if scurve.shape[0]<1920:
                print("Data array is too short, expected 1920 pixels, found "+str(scurve.shape[0]))
                return
            ydim = scurve.shape[1]
            data = np.zeros(self.conf.npixsnom, ydim, dtype = np.int )
            for p in range(self.conf.npixsraw):
                pixel = self.conf.getnompix(p)
                for y in range(ydim):
                    data[pixel,y] = scurve[p,y]
            for r in range(len(row)):
                row[r] = self.conf.getnomrow(row[r])
            for c in range(len(col)):
                col[c] = self.conf.getnomcol(col[c])
        else:
            if scurve.shape[0]<1888:
                print("Data array is too short, expected 1888 pixels, found "+str(scurve.shape[0]))
                return
            if scurve.shape[0]==1888:
                data = scurve
            else:
                data = np.zeros(self.conf.npixsnom, ydim, dtype = np.int )#ensuring correct length of data array
                for p in range(self.conf.npixsnom):
                    pixel = self.conf.getnompix(p)
                    for y in range(ydim):
                        data[pixel,y] = scurve[p,y]            
            
        #th_array = np.zeros(2040, dtype = np.int )
        #noise_array = np.zeros(2040, dtype = np.float )
        th_array = np.zeros(self.conf.npixsnom, dtype = np.int )
        noise_array = np.zeros(self.conf.npixsnom, dtype = np.float )
        nonfittedpixel = []
        for r in row:
            for c in col:
                pixel = self.conf.pixelidnom(r,c)
                plt.plot(range(start,stop), scurve[pixel,0:(stop-start)],'-')
                # Noise and Spread study
                if extract:
                    try:
                        if s_type == "THR":
                            start_DAC = np.argmax(scurve[pixel,:])+10
                            #if pixel == 532: print "pixel 532",start_DAC
                            while start_DAC <  (stop-start):
                                #print start_DAC,scurve[pixel,start_DAC],n_pulse+math.sqrt(n_pulse)
                                if math.fabs(scurve[pixel,start_DAC] - n_pulse) < math.sqrt(n_pulse):
                                    break
                                start_DAC += 1
                            if start_DAC >= stop-start-1:
                                start_DAC = min(np.argmax(scurve[pixel,:])+10,stop-start-5)
                            offset = np.mean(scurve[pixel,-12:-2]) #XXX added offset
                            #print "Fit Range (",start_DAC," to ", (stop-start), ") with initial parameters ", [n_pulse, nominal_DAC, 2],"for pixel",pixel
                            par, cov = curve_fit(self.errorfc, range(start_DAC, (stop-start)), scurve[pixel,start_DAC + 1 :(stop-start) + 1], p0= [n_pulse, (start_DAC+(stop-start))/2, 2, offset]) #XXX added offset
                            #print "Paramters from fit ",par
                            #print "Covariance ", cov
                        elif s_type == "CAL":
                            offset = np.mean(scurve[pixel,1:11]) #XXX added offset
                            start_DAC = 0
                            #print "Fit Range (",start_DAC," to ", (stop-start), ") with initial parameters ", [n_pulse, nominal_DAC, 2]
                            par, cov = curve_fit(self.errorf,  range(start_DAC, (stop-start)), scurve[pixel,start_DAC + 1 :(stop-start) + 1], p0= [n_pulse, nominal_DAC, 2, offset]) #XXX added offset
                            #print "Paramters from fit ",par
                            #print "Covariance ", cov
                        if par[0] == n_pulse and par[1] == nominal_DAC and par[2] == 2:
                            plt.figure(10)
                            print start_DAC + 1,(stop-start) + 1
                            plt.plot(range(start_DAC + 1,(stop-start) + 1), scurve[pixel,start_DAC + 1 :(stop-start) + 1],'-')
                            xaxis = np.arange(start_DAC + 1,(stop-start) + 1,1)
                            fittedscurve = [self.errorfc(i, *[n_pulse, nominal_DAC, 2]) for i in xaxis]
                            plt.plot(range(start_DAC + 1,(stop-start) + 1), fittedscurve)
                            plt.show()
                        if not math.isnan(par[1]): th_array[pixel] = int(round(par[1]))
                        else: th_array[pixel] = -1
                        if not math.isnan(par[2]): noise_array[pixel] = par[2]
                        else: noise_array[pixel] = -1
                    except RuntimeError or TypeError:
                        #print "Fit Range (",start_DAC," to ", (stop-start), ") with initial parameters ", [n_pulse, nominal_DAC, 2]
                        ###for index in range( start_DAC + 1 ,(stop-start) + 1): print "At x = ", index + start," s curve ",  scurve[pixel,index]
                        #print "Fitting failed on column ", c , " row: " ,r
                        ###plt.figure(8)
                        ###for f in nonfittedpixel:
                        ###    plt.plot(range(start_DAC + 1,(stop-start) + 1), scurve[pixel,start_DAC + 1 :(stop-start) + 1],'-')    
                        ###    xaxis = np.arange(start_DAC + 1,(stop-start) + 1,1)
                        ###    fittedscurve = [self.errorfc(i, *[n_pulse, nominal_DAC, 2]) for i in xaxis]
                        ###    plt.plot(range(start_DAC + 1,(stop-start) + 1), fittedscurve)
                        ###    plt.show()
                        nonfittedpixel.append((r,c))
                        ###th_array[pixel] = nominal_DAC
        if extract:
            if printout: print "Non fitable pixels " + str(len(nonfittedpixel))+" (fitting failed)."
        #th_array = [a for a in th_array if a != 0]
        #noise_array = [b for b in noise_array if b != 0]
        
        #plt.show()

        if printout: print("Plotting")
        if s_type == "THR":
            plot_xlabel = "Threshold DAC value"
        if s_type == "CAL":
            plot_xlabel = "Calibration DAC value"
        plot_ylabel = "Counter Value"
        plt.title("SCurves")
        plt.xlabel(plot_xlabel)
        plt.ylabel(plot_ylabel)
        if save_plot:
            plt.savefig(filename+"_curves.png")
            plt.close()
        if extract:
            th_array_formean = [a for a in th_array if a != 0]
            noise_array_formean = [b for b in noise_array if b != 0]
            if printout: print("Figure 2")
            plt.figure(2)
            plt.close()
            if len(th_array) == 1:
                plt.plot(col, th_array, 'o')
            elif len(th_array) != 0:
                if printout: print("Plot Figure 2")
                plt.hist(th_array, bins=range(min(th_array), max(th_array) + 1, 1))
            plt.title("Threshold")
            plt.xlabel(plot_xlabel)
            plt.ylabel(plot_ylabel)
            if save_plot:
                plt.savefig(filename+"_threshold.png")
                plt.close()
            plt.figure(3)
            if printout: print("Figure 3")
            if len(noise_array) == 1:
                plt.plot(col, noise_array, 'o')
            elif len(noise_array) != 0:
                if printout: print("Plot Figure 3")
                plt.hist(noise_array, bins=np.arange(min(noise_array), max(noise_array) + 1 , 0.1))
            #print "Threshold Average: ", np.mean(th_array), " Spread SD: ", np.std(th_array)
            #print "Noise Average: ", np.mean(noise_array), " Spread SD: ", np.std(noise_array)
            plt.title("Noise")
            plt.xlabel(plot_xlabel)
            plt.ylabel(plot_ylabel)
            if save_plot:
                plt.savefig(filename+"_noise.png")
                plt.close()
            if printout: print "Threshold Average: ", np.mean(th_array_formean), " Spread SD: ", np.std(th_array_formean)
            if printout: print "Noise Average: ", np.mean(noise_array_formean), " Spread SD: ", np.std(noise_array_formean)
        if extract:
            return     th_array, noise_array, nonfittedpixel
        return nonfittedpixel
# Readout Counters current method
#uses RAW pixel id

cal = mpa_cal_utility()
    
def RefitSCurvesOneChip(mapsa, chip,v2=False):
    moduleid = "mpa_test_"+mapsa
    thepath = mapsa+"/"
    teststring = moduleid + "_"+chip+"_"
    bbstring = "_BumpBonding_SCurve_BadBump"
    if v2: bbstring = "_BumpBondingv2_SCurve_BadBump"
    filelist = glob.glob(thepath+teststring+"*"+bbstring+".csv")
    maxstr = "" #string with highest number
    maxnbr = 0
    if len(filelist) == 0:
        print("No S Curve found.")
        return False
    for f in filelist:
        numbers_from_string = filter(str.isdigit, f)
        if numbers_from_string > maxnbr:
            maxnbr = numbers_from_string
            maxstr = f
    #print(maxstr)
    #print(maxstr[0:-4])
    csvfilename = maxstr
    data_array = loadSCurvesFromCSV(csvfilename)
    extract_val = int(250. * 95./218.)
    #print(np.array(cal.conf.rowsnom))
    th_array, noise_array, notfittedpixels = cal.plot_extract_scurve(row = np.array(cal.conf.rowsnom), col = np.array(cal.conf.colsnom), s_type = "CAL", scurve = np.array(data_array) , n_pulse = 1000, nominal_DAC = extract_val, start = 0, stop = 256, extract = True, plot = False, save_plot = True, israw = False, printout = 0, filename = "./x")
    newcsvfilename = csvfilename[0:-4] + "_refitted.csv"
    newcsvfilename = newcsvfilename.replace("SCurve","Noise")
    CSV.ArrayToCSV (noise_array, newcsvfilename)
    return True

def RefitSCurves():
    chiplist = ["Chip1","Chip2","Chip3","Chip4","Chip5","Chip6","Chip7","Chip8","Chip9","Chip10","Chip11","Chip12","Chip13","Chip14","Chip15","Chip16"]
    #allmapsas = ["AEMTec1","AEMTec2","AEMTec3","AEMTec4","AEMTec5","AEMTec6","AEMTec7","AEMTec8","AEMTec9","AEMTec10","HPK2_1","HPK15_1","HPK15_2","HPK16_1","HPK17_1","HPK18_1","HPK19_1","HPK25_2","HPK26_1","HPK26_2"]
    allmapsas = ["AEMTec2","AEMTec3","AEMTec4","AEMTec5","AEMTec6","AEMTec7","AEMTec8","AEMTec9","AEMTec10","HPK2_1","HPK15_1","HPK15_2","HPK16_1","HPK17_1","HPK18_1","HPK19_1","HPK25_2","HPK26_1","HPK26_2"]
    for mapsa in allmapsas:
        for chip in chiplist:
            print("Refit Bump Bonding SCurve for " + chip + " of MaPSA " + mapsa)
            RefitSCurvesOneChip(mapsa,chip,False)
            print("Refit Bump Bonding v2 SCurve for " + chip + " of MaPSA " + mapsa)
            RefitSCurvesOneChip(mapsa,chip,True)

def RedoVCal3CutoffOneChip(mapsa, chip,v2=False):
    moduleid = "mpa_test_"+mapsa
    thepath = mapsa+"/"
    teststring = moduleid + "_"+chip+"_"
    bbstring = "_BumpBonding_Noise_BadBump_refitted"
    if v2: bbstring = "_BumpBondingv2_Noise_BadBump_refitted"
    filelist = glob.glob(thepath+teststring+"*"+bbstring+".csv")
    maxstr = "" #string with highest number
    maxnbr = 0
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
    newcsvfilename = newcsvfilename.replace("Noise_BadBump","BadBumpVCal3")
    #print(csvfilename,len(filelist))
    data_array = loadSCurvesFromCSV(csvfilename)
    #data_array.insert(0, ["",0])
    newdata_array = []
    for i in range(len(data_array)):
        #print(i,data_array[i])
        if data_array[i][0]=="":
            newdata_array.append(["",""])
        else:
            if data_array[i][0]<0:
                newdata_array.append([i,-2])
            elif data_array[i][0]==0:
                newdata_array.append([i,-1])
            elif data_array[i][0]<=3:
                newdata_array.append([i,1])
            else:
                newdata_array.append([i,0])

    newdata_array.insert(0, ["",0])
    #with open("out.csv", "w", newline="") as f:
    with open(newcsvfilename, "w") as f:
        writer = csv.writer(f)
        writer.writerows(newdata_array)
    
def RedoVCal3Cutoffs():
    chiplist = ["Chip1","Chip2","Chip3","Chip4","Chip5","Chip6","Chip7","Chip8","Chip9","Chip10","Chip11","Chip12","Chip13","Chip14","Chip15","Chip16"]
    #allmapsas = ["AEMTec1","AEMTec2","AEMTec3","AEMTec4","AEMTec5","AEMTec6","AEMTec7","AEMTec8","AEMTec9","AEMTec10","HPK2_1","HPK15_1","HPK15_2","HPK16_1","HPK17_1","HPK18_1","HPK19_1","HPK25_2","HPK26_1","HPK26_2"]
    allmapsas = ["AEMTec1","AEMTec2","AEMTec3","AEMTec4","AEMTec5","AEMTec6","AEMTec7","AEMTec8","AEMTec9","AEMTec10","HPK2_1","HPK15_1","HPK15_2","HPK16_1","HPK17_1","HPK18_1","HPK19_1","HPK25_2","HPK26_1","HPK26_2"]
    for mapsa in allmapsas:
        for chip in chiplist:
            print("Redo Bump Bonding VCal3 for " + chip + " of MaPSA " + mapsa)
            RedoVCal3CutoffOneChip(mapsa,chip,False)
            print("Redo Bump Bonding v2 VCal3 for " + chip + " of MaPSA " + mapsa)
            RedoVCal3CutoffOneChip(mapsa,chip,True)

        
#RefitSCurves()
RedoVCal3Cutoffs()

