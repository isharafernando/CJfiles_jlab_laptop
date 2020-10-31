# -*- coding: utf-8 -*-
"""
Created on Mon Jun 17 11:05:15 2019

@author: isharafernando
"""

#import numpy as np
import matplotlib     
matplotlib.rc('xtick', labelsize=7)     
matplotlib.rc('ytick', labelsize=7)
import glob



import matplotlib.pyplot as plt
#import numpy as np
import ioff1definitions
import chi2_tables

#FileNames=['./Data/CJ15','./Data/fit_ioff1_CJ15-2-new','./Data/fit_ioff1_CJ15-3-new','./Data/fit_ioff1_CJ15-2-noW','./Data/fit_ioff1_CJ15-3-noW','./Data/fit_ioff1_CJ15-2-prime','./Data/fit_ioff1_CJ15-3-prime','./Data/fit_ioff1_CJ15-2-noW-prime','./Data/fit_ioff1_CJ15-3-noW-prime'];
#filenames=['CJ15','fit_ioff1_CJ15-2-new','fit_ioff1_CJ15-3-new','fit_ioff1_CJ15-2-noW','fit_ioff1_CJ15-3-noW','fit_ioff1_CJ15-2-prime','fit_ioff1_CJ15-3-prime','fit_ioff1_CJ15-2-noW-prime','fit_ioff1_CJ15-3-noW-prime'];
#filenames=['CJ15','fit_ioff1_CJ15-2-new','fit_ioff1_CJ15-3-new','fit_ioff1_CJ15-2-noW'];
#FileNames=['./Data/'+i for i in filenames]
#ALabels=['CJ15 (ref)','ioff1_CJ15-2','ioff1_CJ15-3','ioff1_CJ15-2-noW','ioff1_CJ15-3-noW',"ioff1_CJ15'-2","ioff1_CJ15'-3","ioff1_CJ15'-2-noW","ioff1_CJ15'-3-noW"]
#ALabels=['CJ15 (ref)','ioff1_CJ15-2','ioff1_CJ15-3','ioff1_CJ15-2-noW']
#,'i9_CJ15',,'ioff9_CJ15','i1_CJ15-3-NoJL','ioff1_CJ15-3-NoJLabData',-3+JLabAC-Wasy, 'i1n_CJ15-JAC-3'
filenames=['CJ15', 'i1n_CJ15-3-NoJAC+BONUS'];
FileNames=['./Data/'+i for i in filenames]
ALabels=['CJ15 (ref)',"ioff1new_2+JLabAC"]
PlotOS=[1,1]

test_par_files = glob.glob('./Data/par_'+filenames[0]+'/*.par')
#test_par =[i+'.par' for i in FileNames]
#print test_par_files[0]
#print len(test_par_files)
#print ioff1definitions.ioffsforErr(test_par_files[0],'off2')
#print ioff1definitions.OSforSingleX(0.5,test_par_files[0],1)
#print ioff1definitions.OSforSingleX(0.5,test_par_files[1])
#print ioff1definitions.OSforSingleFile(test_par_files[0],1)
#print ioff1definitions.deltaOSPlus(filenames[0],1)
#print ioff1definitions.deltaOSMinus(filenames[1],2)

#print test_par
#print ioff1definitions.parline(test_par[0],'Covariance')
#print ioff1definitions.parline(test_par[0],'Covariance')

#print ioff1definitions.deltaOSPlusTest(filenames[1])



#FileNames=['./Data/CJ15','./Data/fit_ioff1_CJ15-2','./Data/fit_ioff1_CJ15-3','./Data/fit_ioff1_CJ15-2-noW','./Data/fit_ioff1_CJ15-3-noW','./Data/fit_ioff1_CJ15-2-prime','./Data/fit_ioff1_CJ15-3-prime','./Data/fit_ioff1_CJ15-2-noW-prime','./Data/fit_ioff1_CJ15-3-noW-prime','./Data/fit_ioff1_CJ15-3-new'];
#FileNames=['./Data/CJ15','./Data/fit_ioff1_CJ15-2','./Data/fit_ioff1_CJ15-3'];
#ALabels=['CJ15 (ref)','ioff1_CJ15-2','ioff1_CJ15-3']
#ALabels=['CJ15 (ref)','ioff1_CJ15-2','ioff1_CJ15-3','ioff1_CJ15-2-noW','ioff1_CJ15-3-noW',"ioff1_CJ15'-2","ioff1_CJ15'-3","ioff1_CJ15'-2-noW","ioff1_CJ15'-3-noW",'ioff1_CJ15-3-tol-0.1']
# line width
#LW=0.7

DataFileNames=[i+'.pdf' for i in FileNames]

OutFileNames=[i+'.out' for i in FileNames]

ParFileNames=[i+'.par' for i in FileNames]

#print ioff1definitions.parline(ParFileNames[0],'Covariance')
#print ioff1definitions.parcov(ParFileNames[0],'off2','off3')
print ioff1definitions.parcov(ParFileNames[1],'off1','off2')
print ioff1definitions.parcov(ParFileNames[1],'off1','off3')
#print ioff1definitions.parcov(ParFileNames[1],'off1','off4')
print ioff1definitions.parcov(ParFileNames[1],'off2','off3')
#print ioff1definitions.parcov(ParFileNames[1],'off2','off4')
#print ioff1definitions.parcov(ParFileNames[1],'off3','off4')
#print ioff1definitions.parcov(ParFileNames[1],'off1','off5')
#print ioff1definitions.parcov(ParFileNames[1],'off2','off5')
#print ioff1definitions.parcov(ParFileNames[1],'off3','off5')
#print ioff1definitions.parcov(ParFileNames[1],'off4','off5')
#print ioff1definitions.parcov(ParFileNames[1],'off1','off6')
#print ioff1definitions.parcov(ParFileNames[1],'off2','off6')
#print ioff1definitions.parcov(ParFileNames[1],'off3','off6')
#print ioff1definitions.parcov(ParFileNames[1],'off4','off6')
#print ioff1definitions.parcov(ParFileNames[1],'off5','off6')
#print ioff1definitions.parcov(ParFileNames[1],'off2','off3')


DatFiles=[]
lists=[]
for i in range(0,len(DataFileNames)):
    templistDat=[]
    DatFiles.append(open(DataFileNames[i]))
    for line in open(DataFileNames[i]):
        templistDat += [line.split()]
    lists.append(templistDat)    
    
OutFiles=[]
Outlists=[]
for i in range(0,len(OutFileNames)):
    templistOut=[]
    OutFiles.append(open(OutFileNames[i]))
    for line in open(OutFileNames[i]):
        templistOut += [line.split()]
    Outlists.append(templistOut)    

ParFiles=[]
Parlists=[]
for i in range(0,len(ParFileNames)):
    templistPar=[]
    OutFiles.append(open(ParFileNames[i]))
    for line in open(ParFileNames[i]):
        templistPar += [line.split()]
    Parlists.append(templistPar)    


#print ioff1definitions.OffShellPlot(PlotOS,Outlists,OutFileNames,Parlists,ParFileNames,filenames,-0.7,0.3,ALabels)
   

def duPlotsNew(TN=2,alignment=1.0):
    plt.subplot(2,2,1)
    ioff1definitions.IndividualPlot(lists,DataFileNames,TN,5,0.85,1.15,ALabels)
    plt.subplot(2,2,2)
    #-1.0,0.6, ,-0.7,0.3
    ioff1definitions.OffShellPlot(PlotOS,Outlists,OutFileNames,Parlists,ParFileNames,filenames,lists,DataFileNames,-0.7,0.3,ALabels)
    plt.subplot(2,2,3)
    ioff1definitions.SinglePlot(lists,DataFileNames,TN,5,0.0,1.0,ALabels)
    plt.subplot(2,2,4)
    ioff1definitions.Legendplot(Outlists,DataFileNames,OutFileNames,alignment,ALabels)
    plt.suptitle('Initial conditions: pure CJ15')
    
    
def ugPlotsNew(TN=2):
    plt.subplot(2,2,1)
    ioff1definitions.IndividualPlot(lists,DataFileNames,TN,9,0.98,1.02,ALabels)
    plt.subplot(2,2,2)
    ioff1definitions.IndividualPlot(lists,DataFileNames,TN,11,0.95,1.5,ALabels)
    plt.subplot(2,2,3)
    ioff1definitions.RelativeRatioPlot(lists,DataFileNames,TN,9,0.9,1.4,ALabels)
    plt.subplot(2,2,4)
    ioff1definitions.RelativeRatioPlot(lists,DataFileNames,TN,11,0.9,1.4,ALabels)
    plt.suptitle('Initial conditions: pure CJ15')        
  

def dNduPlots(TN=2,alignment=1.0):
    plt.subplot(2,2,1)
    ioff1definitions.IndividualPlot(lists,DataFileNames,TN,7,0.85,1.05,ALabels)
    plt.subplot(2,2,2)
    ioff1definitions.IndividualPlot(lists,DataFileNames,TN,5,0.85,1.05,ALabels)
    plt.subplot(2,2,3)
    #ioff1definitions.Legendplot(Outlists,DataFileNames,OutFileNames,alignment,ALabels) 
    #plt.suptitle('Initial conditions: pure CJ15')
    ioff1definitions.RelativeRatioPlot(lists,DataFileNames,TN,7,0.4,3.0,ALabels)
    plt.subplot(2,2,4)
    ioff1definitions.RelativeRatioPlot(lists,DataFileNames,TN,5,0.4,3.0,ALabels)
    plt.suptitle('Initial conditions: pure CJ15')        

f1=plt.figure(1)
f1.subplots_adjust(hspace=0.4, wspace=0.4)
duPlotsNew(2,1.05)

f2=plt.figure(2)
f2.subplots_adjust(hspace=0.4, wspace=0.4)
ugPlotsNew(2)    

f3=plt.figure(3)
f3.subplots_adjust(hspace=0.4, wspace=0.4)
dNduPlots(2,0.8)    

import matplotlib.backends.backend_pdf
pdf = matplotlib.backends.backend_pdf.PdfPages("./Resultsioff1/iOff1-Set1new.pdf")
for fig in xrange(1,4):
    pdf.savefig( fig )
pdf.close()
print chi2_tables.summary_chi2(FileNames)


#import fpdf
#
#PDF = fpdf.FPDF(format='letter')
#PDF.add_page()
#PDF.set_font("Arial", size=12)
#PDF.write(chi2_tables.summary_chi2(FileNames))
#f3=plt.figure(3)
#f3.subplots_adjust(hspace=0.4, wspace=0.4)
#print chi2_tables.summary_chi2(FileNames) 




      
