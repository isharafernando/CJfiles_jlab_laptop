# -*- coding: utf-8 -*-
"""
Created on Mon Jun 17 11:05:15 2019

@author: isharafernando
"""


import matplotlib     
matplotlib.rc('xtick', labelsize=7)     
matplotlib.rc('ytick', labelsize=7)

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.cm import get_cmap
import glob
#cmap = plt.cm.Spectral
#cmap = plt.cm.get_cmap("Paired,Accent,Set3")
#cmap = matplotlib.cm.get_cmap("brewer_Set3_12")
#cmap = plt.cm.cmap_d
cmap = plt.cm.get_cmap("Paired")

#print len(cmap)

LW=0.7



pdflabel=['$x$','$\kappa$','$\delta \kappa$','db/ub','$\delta db/ub$','$d/u$','$\delta d/u$','$xd$','$\delta xd$','$xu$','$\delta xu$','$xg$','$\delta xg$','$xub$','$\dela xub$','$xdb$','$\delta xdb$','$xs$','$\delta xs$','$xc$','$\delta xc$','$xb$','$\delta xb$']
tablelabel=['$Q^2$=1.69 GeV$^2$','$Q^2$=2.0 GeV$^2$','$Q^2$=10.0 GeV$^2$','$Q^2$=25.0 GeV$^2$','$Q^2$=64.0 GeV$^2$','$Q^2$=100.0 GeV$^2$','$Q^2$=1000.0 GeV$^2$','$Q^2$=8300.0 GeV$^2$','$x=0.1$','$x=0.3$','$x=0.5$','$x=0.7$','$x=0.85$']
xlabel=['$x$','$Q^2$']    

    
# This function automatically identifies the line numbers of "Q2"
def linenumbers(filename):
    temparray1=[]
    with open(filename) as myFile:   
        for num, line in enumerate(myFile, 0):
            if 'Q2' in line:
                temparray1.append(num)          
    return temparray1

def parline(filename,word):
    temparray2=[]
    with open(filename) as myFile:   
        for num, line in enumerate(myFile, 0):
            if word in line:
                temparray2.append(num)          
    return temparray2    

def TOTALchi2(Outlists,OutFileNames,filename):
    tempchi2=float(Outlists[OutFileNames.index(filename)][parline(filename,'TOTAL+norm')[0]][1])
    #tempchi2=float(Outlists[OutFileNames.index(filename)][parline(filename,'TOTAL+norm')[0]][-1])
    return tempchi2
# This function will return the column we are interested in
# In this function you have to give the start and end line numbers as an input    
def pdfSelect(listname,column,start,stop):
    TempColumn=[listname[i][column] for i in range(start,stop)]
    return TempColumn
        

# This function will provide the "start" and "end" points of the tables in the file
# There are 13 tables (8 tabes for fixed Q2 and 5 tables for fixed x values)
def LinNumTbles(lists,filename,Numbr):
    lins=linenumbers(filename)
    temp1=lins[:len(lins)-5]
    temp2=[lins[i]-2 for i in range(8,len(lins))]
    newlins= temp1 + temp2
    TblStrtRef=[newlins[i]+3 for i in range(0,len(newlins))]
    TblEndRefTEMP=[newlins[i+1]-3 for i in range(0,len(newlins)-1)]
    TblEndRef=TblEndRefTEMP+[len(lists[0])]
    temparray=[TblStrtRef[Numbr],TblEndRef[Numbr]]
    return temparray

def ioffs(Parlists,ParFileNames,filename,val):
    tempN=float(Parlists[ParFileNames.index(filename)][parline(filename,val)[0]][1])
    tempErrN=float(Parlists[ParFileNames.index(filename)][parline(filename,val)[0]][2])
    return tempN, tempErrN        
    

def parselection(filename):
    with open(filename) as myFile:   
        for num, line in enumerate(myFile, 0):
            if "off4" in line:           
                temppar=3
            else:
                temppar=2
    return temppar


def lineoff4(filename):
    off4exist=False
    temppar=2
    with open(filename) as myFile:   
        for num, line in enumerate(myFile, 0):
            if 'off4' in line:
                off4exist=True
                temppar=3
    with open(filename) as myFile:            
        if off4exist==True:
            for num, line in enumerate(myFile, 0):
                if 'off4          0.0000        0.0000' in line:
                    off4exist=False
                    temppar=2      
    return temppar
    

def Nval(filename,Outlists,OutFileNames):
    tempN=float(Outlists[OutFileNames.index(filename)][parline(filename,'N =')[0]][2])
    tempErrN=float(Outlists[OutFileNames.index(filename)][parline(filename,'N =')[0]][4])
    return tempN, tempErrN
    
def x0val(filename,Outlists,OutFileNames):
    tempx0=float(Outlists[OutFileNames.index(filename)][parline(filename,'x0 =')[0]][2])
    tempErrx0=float(Outlists[OutFileNames.index(filename)][parline(filename,'x0 =')[0]][4])
    return tempx0, tempErrx0

def x1val(filename,Outlists,OutFileNames):
    tempx1=float(Outlists[OutFileNames.index(filename)][parline(filename,'x1 =')[0]][2])
    return tempx1 
    
    
def OffShell(x,N,x0,x1):
    tempOffShell=N*(x-x0)*(x-x1)*(1+x0-x)
    return tempOffShell
    
def OffShelliOff12(x,N,x0,x1):
    #tempOffShell=N*(x-x0)*(x-x1)
    tempOffShell= N + x0*x + x1*(x**2)
    return tempOffShell    
    
def OffShelliOff13(x,N,x0,x1,x2):
    #tempOffShell=N*(x-x0)*(x-x1)*(x-x2)
    tempOffShell= N + x0*x + x1*(x**2) + x2*(x**3)
    return tempOffShell

def deltaOffShell(x,N,dN,x0,dx0,x1,cov):
    dfdN=(x-x0)*(x-x1)*(1+x0-x)
    dfdx0=N*(x-x1)*(2*(x-x0)-1)
    #cov= -0.067
    dfdx1=-N*(1+x0-x)
    dx1=0.005
    tempdOffSh2=(dfdN**2)*(dN**2)+(dfdx0**2)*(dx0**2)+(2*(dfdN)*(dfdx0)*cov*dN*dx0)+(dfdx1**2)*(dx1**2)
    tempdOffSh=np.sqrt(np.abs(tempdOffSh2))
    return tempdOffSh


def ioffsforErr(filename,val):
    templistpar=[]
    for line in open(filename):
        templistpar += [line.split()] 
    tempNN=templistpar[parline(filename,val)[0]][1]
    #For the uncertainty value one can obtain the [2] of the above line instead of [1]    
    return float(tempNN)

def parcov(filename,par1,par2):
    templistpc=[]
    StartLine=parline(filename,'Covariance matrix')[0]
    par1line=parline(filename,par1)[0]
    par2line=parline(filename,par2)[0]
    parLAMBDA=parline(filename,'LAMBDA')[0]
    for line in open(filename):
        templistpc += [line.split()] 
    count1=0
    count2=0
    for i in range(parLAMBDA,par1line):
        if templistpc[i][2]!='0.0000':
            count1=count1+1 
    for i in range(parLAMBDA,par2line):
        if templistpc[i][2]!='0.0000':
            count2=count2+1     
        #print templistpc[i][2]
    tt1=float(templistpc[StartLine+count1+1][count1])
    tt2=float(templistpc[StartLine+count1+1][count2])
    tt3=float(templistpc[StartLine+count2+1][count1])
    tt4=float(templistpc[StartLine+count2+1][count2])
    
#    uu1=count1+1
#    uu2=count2+1
#    uu3=templistpc[StartLine+count1+1]
#    uu4=templistpc[StartLine+count2+1]
#    uu5=templistpc[StartLine+count1+1][count2]
#    uu6=templistpc[StartLine+count2+1][count1]
    
    s1=np.sqrt(abs(tt1))
    s2=np.sqrt(abs(tt4))
    
    ttt1=tt1/(s1*s1)
    ttt2=tt2/(s1*s2)
    ttt3=tt3/(s1*s2)
    ttt4=tt4/(s2*s2)
    #tt5=templistpc[parLAMBDA-1+35],count1,count2,tt5[0],tt1,tt3
    #,float(tt2),float(tt3),float(tt4)
    #ttt1,ttt2,ttt3,ttt4
    return ttt1,ttt2,ttt3,ttt4

#def OSforSingleX(x,filename):
#    tempord=lineoff4(filename)
#    tN=ioffsforErr(filename,'off1')
#    tx0=ioffsforErr(filename,'off2')
#    tx1=ioffsforErr(filename,'off3')
#    if tempord==3:
#        tx2=ioffsforErr(filename,'off4')
#        tempOSX=OffShelliOff13(x,tN,tx0,tx1,tx2)
#    else:
#        tx2=0
#        tempOSX=OffShelliOff12(x,tN,tx0,tx1)
#    return tempOSX  

#For convensional offshell use flag 1, and use flag 2 for anything else
def OSforSingleX(x,filename,y):
    tempord=lineoff4(filename)
    tN=ioffsforErr(filename,'off1')
    tx0=ioffsforErr(filename,'off2')
    tx1=ioffsforErr(filename,'off3')
    if y==1:
        tempOSX=OffShell(x,tN,tx0,tx1)
    elif y==2 and tempord==3:
        tx2=ioffsforErr(filename,'off4')
        tempOSX=OffShelliOff13(x,tN,tx0,tx1,tx2)
    else:
        tx2=0
        #Now we need to check whether this is ioff9
        #Because in ioff9, x1 is free and use the same OS function as CJ15
        templist=[]
        for line in open(filename):
            templist += [line.split()]
        ioffline=parline(filename,'inuke')[0]
        ioff9check=templist[ioffline][1][2]   
        if ioff9check=='9':    
            tempOSX=OffShell(x,tN,tx0,tx1)
        else:
            tempOSX=OffShelliOff12(x,tN,tx0,tx1)
        #tempOSX=OffShelliOff12(x,tN,tx0,tx1)
    return tempOSX    

def OSforSingleFile(filename,y):
    x=np.linspace(0,1,100)
    tempvals=[OSforSingleX(i,filename,y) for i in x]
    return tempvals

def deltaOSPlus(filename,y):
    test_files = glob.glob('./Data/par_'+filename+'/*.par')
    #print test_files[0]
    nn=len(test_files)
    nnn=(nn-1)/2
    x=np.linspace(0,1,100)
    tempvvv=[]
    for i in x:
        tempsum=0
        for j in range(1,nnn+1):
            tempvvv01=OSforSingleX(i,test_files[2*j-1],y)-OSforSingleX(i,test_files[0],y)
            tempvvv02=OSforSingleX(i,test_files[2*j],y)-OSforSingleX(i,test_files[0],y)
            tempvvv1=max(tempvvv01,tempvvv02,0)
            tempvvv2=np.square(tempvvv1)            
            tempsum=tempsum+tempvvv2
        tempvvv.append(np.sqrt(tempsum))
    return tempvvv
    
def deltaOSMinus(filename,y):
    test_files = glob.glob('./Data/par_'+filename+'/*.par')
    #print test_files[0]
    nn=len(test_files)
    nnn=(nn-1)/2
    x=np.linspace(0,1,100)
    tempvvv=[]
    for i in x:
        tempsum=0
        for j in range(1,nnn+1):
            tempvvv01=OSforSingleX(i,test_files[0],y)-OSforSingleX(i,test_files[2*j-1],y)
            tempvvv02=OSforSingleX(i,test_files[0],y)-OSforSingleX(i,test_files[2*j],y)
            tempvvv1=max(tempvvv01,tempvvv02,0)
            tempvvv2=np.square(tempvvv1)            
            tempsum=tempsum+tempvvv2
        tempvvv.append(np.sqrt(tempsum))
    return tempvvv    
 
    

#def ioffsforErr(filename,val):
#    templistpar=[]
#    for line in open(filename):
#        templistpar += [line.split()] 
#    tempNN=templistpar[parline(filename,val)[0]][1]
#    #For the uncertainty value one can obtain the [2] of the above line instead of [1]    
#    return float(tempNN)
    
        
def PlotData(lists,listname,filename,TblNo,ColNo):
    tempx=pdfSelect(listname,0,LinNumTbles(lists,filename,TblNo)[0],LinNumTbles(lists,filename,TblNo)[1])
    tempy=pdfSelect(listname,ColNo,LinNumTbles(lists,filename,TblNo)[0],LinNumTbles(lists,filename,TblNo)[1])
    tempyErr=pdfSelect(listname,ColNo+1,LinNumTbles(lists,filename,TblNo)[0],LinNumTbles(lists,filename,TblNo)[1])
    tx=np.array([float(i) for i in tempx])
    ty=np.array([float(i) for i in tempy])
    tyErr=np.array([float(i) for i in tempyErr])
    return tx,ty,tyErr
    
    
#def PlotData(lists,listname,filename,TblNo,ColNo):
#    tempx=pdfSelect(listname,0,LinNumTbles(lists,filename,TblNo)[0],LinNumTbles(lists,filename,TblNo)[1])
#    tempy=pdfSelect(listname,ColNo,LinNumTbles(lists,filename,TblNo)[0],LinNumTbles(lists,filename,TblNo)[1])
#    tempyErr=pdfSelect(listname,ColNo+1,LinNumTbles(lists,filename,TblNo)[0],LinNumTbles(lists,filename,TblNo)[1])
#    tx=np.array([float(i) for i in tempx])
#    ty=np.array([float(i) for i in tempy])
#    tyErr=np.array([float(i) for i in tempyErr])
#    return tx,ty,tyErr


#def SinglePlot(DataFileName1,DataFileName2,TN,CN,XL,AL):
#    x0=PlotData(lists[DataFileName1],DataFileNames[DataFileName1],TN,CN)[0]
#    y1=PlotData(lists[DataFileName1],DataFileNames[DataFileName1],TN,CN)[1]
#    y1Err=PlotData(lists[DataFileName1],DataFileNames[DataFileName1],TN,CN)[2]
#    y2=PlotData(lists[DataFileName2],DataFileNames[DataFileName2],TN,CN)[1]
#    y2Err=PlotData(lists[DataFileName2],DataFileNames[DataFileName2],TN,CN)[2]
#    plt.plot(x0,y2,'r-',label=AL[DataFileName2])
#    plt.plot(x0,y1,'b-',label=AL[DataFileName1])
#    plt.legend(bbox_to_anchor=(0.5, 1), loc=1, borderaxespad=0.01, fontsize=5)
#    #plt.legend(bbox_to_anchor=(0.96, 1), loc=1, borderaxespad=0.01, prop={'size': 3.5}, fontsize=16)    
#    #plt.ylim([1-np.max(y1Err)*6,1+np.max(y1Err)*6])
#    plt.fill_between(x0,(y1-y1Err),(y1+y1Err),color='lightblue',alpha=.9)
#    plt.fill_between(x0,(y2-y2Err),(y2+y2Err),color='pink',alpha=.7)
#    plt.xlabel(xlabel[XL])
#    plt.ylabel('$\phi$')
#    #plt.text(0.85,1+(np.max(y1Err))*3,pdflabel[CN],fontsize='10')
#    #plt.text(0.05,1-(np.max(y1Err))*5,tablelabel[TN],fontsize='8')
#    plt.rc('font',size=8)
#    #plt.rc('legend', fontsize=14)

def SinglePlot(lists,DataFileNames,TN,CN,ymin,ymax,ALabels):
    x0=PlotData(lists,lists[0],DataFileNames[0],TN,CN)[0]
    y1=PlotData(lists,lists[0],DataFileNames[0],TN,CN)[1]
    y1Err=PlotData(lists,lists[0],DataFileNames[0],TN,CN)[2]
    line_colors = cmap(np.linspace(0,1,12))
    for i in range(1,len(DataFileNames)):
        temp1=PlotData(lists,lists[i],DataFileNames[i],TN,CN)[1]
        #tempratio=temp1/y1
        plt.plot(x0,temp1,'-',label=ALabels[i],linewidth=LW,c=line_colors[i-1])
    plt.plot(x0,y1,'--',color='deepskyblue',label=ALabels[0],linewidth=LW)
    plt.ylim([ymin,ymax])
    plt.fill_between(x0,(y1-y1Err),(y1+y1Err),color='white',alpha=0.7,facecolor='deepskyblue')
    plt.xlabel(xlabel[0],fontsize='10')
    plt.ylabel('$\phi$',fontsize='10')
    plt.text(0.05,ymin+(ymax-ymin)*0.85,pdflabel[CN],fontsize='10')
    plt.text(0.05,ymin+(ymax-ymin)*0.1,tablelabel[TN],fontsize='8')



#len(ALabels)
def IndividualPlot(lists,DataFileNames,TN,CN,ymin,ymax,ALabels):
    x0=PlotData(lists,lists[0],DataFileNames[0],TN,CN)[0]
    y1=PlotData(lists,lists[0],DataFileNames[0],TN,CN)[1]
    y1Err=PlotData(lists,lists[0],DataFileNames[0],TN,CN)[2]
    line_colors = cmap(np.linspace(0,1,12))
    for i in range(1,len(DataFileNames)):
        temp1=PlotData(lists,lists[i],DataFileNames[i],TN,CN)[1]
        tempratio=temp1/y1
        plt.plot(x0,tempratio,'-',label=ALabels[i],linewidth=LW,c=line_colors[i-1])
    plt.plot(x0,y1/y1,'--',color='deepskyblue',label=ALabels[0],linewidth=LW)
    plt.ylim([ymin,ymax])
    plt.fill_between(x0,(y1-y1Err)/y1,(y1+y1Err)/y1,color='white',alpha=0.45,facecolor='deepskyblue')
    plt.xlabel(xlabel[0],fontsize='10')
    plt.ylabel('$\phi/\phi_{ref}$',fontsize='10')
    plt.text(0.05,ymin+(ymax-ymin)*0.85,pdflabel[CN],fontsize='10')
    plt.text(0.05,ymin+(ymax-ymin)*0.1,tablelabel[TN],fontsize='8')


#def OffShellPlot(PlotNot,Outlists,OutFileNames,Parlists,ParFileNames,filenames,ymin,ymax,ALabels):
#    x=np.linspace(0,1,100)
#    N1=Nval(OutFileNames[0],Outlists,OutFileNames)[0]
#    x01=x0val(OutFileNames[0],Outlists,OutFileNames)[0]
#    x11=x1val(OutFileNames[0],Outlists,OutFileNames)
#    Nvalues=[]
#    x0values=[]
#    x1values=[]
#    x2values=[]
#    for i in range(1,len(ParFileNames)):
#        Nvalues.append(ioffs(Parlists,ParFileNames,ParFileNames[i],'off1')[0])
#        x0values.append(ioffs(Parlists,ParFileNames,ParFileNames[i],'off2')[0])
#        x1values.append(ioffs(Parlists,ParFileNames,ParFileNames[i],'off3')[0])
#        if lineoff4(ParFileNames[i])==3:
#            x2values.append(ioffs(Parlists,ParFileNames,ParFileNames[i],'off4')[0])
#        else:
#            x2values.append(0)
#    line_colors=cmap(np.linspace(0,1,12))
#    KPN=8.1
#    KPdN=0.3
#    KPdNsys=0.526
#    KPx0=0.448
#    KPdx0=0.005
#    KPdx0sys=0.0067
#    KPx1=0.050
#    cov=-0.067
#    covsys=0
#    y1=OffShell(x,N1,x01,x11)
#    dy1p=deltaOSPlus(filenames[0],1)
#    dy1m=deltaOSMinus(filenames[0],1)
#    for i in range(0,len(ParFileNames)-1):
#        if lineoff4(ParFileNames[i+1])==2:
#            templist=[]
#            for line in open(ParFileNames[i+1]):
#                templist += [line.split()]
#            ioff9=parline(ParFileNames[i+1],'10914')[0]
#            print ioff9
#            #print templist[58]    
#            if ioff9:
#                temp1=OffShell(x,Nvalues[i],x0values[i],x1values[i])
#            else:
#                temp1=OffShelliOff12(x,Nvalues[i],x0values[i],x1values[i])
#            #temp1=OffShelliOff12(x,Nvalues[i],x0values[i],x1values[i])    
#            #iOff12
#            #print lineoff4(ParFileNames[i+1])
#            temp1yp=deltaOSPlus(filenames[i+1],2)
#            temp1ym=deltaOSMinus(filenames[i+1],2)
#            plt.plot(x,temp1,'-',label=ALabels[i],linewidth=LW,c=line_colors[i])
#            if PlotNot[i+1]==1:                
#                plt.fill_between(x,(temp1-temp1ym),(temp1+temp1yp),color=line_colors[i],alpha=0.5,linewidth=0.01)
#        else:
#            temp2=OffShelliOff13(x,Nvalues[i],x0values[i],x1values[i],x2values[i])
#            #print lineoff4(ParFileNames[i+1])
#            #plt.plot(x,temp2,'-',label=ALabels[i],linewidth=LW,c=line_colors[i])
#            temp2yp=deltaOSPlus(filenames[i+1],2)
#            temp2ym=deltaOSMinus(filenames[i+1],2)
#            plt.plot(x,temp2,'-',label=ALabels[i],linewidth=LW,c=line_colors[i])
#            if PlotNot[i+1]==1:
#                plt.fill_between(x,(temp2-temp2ym),(temp2+temp2yp),color=line_colors[i],alpha=0.5,linewidth=0.01)
#    plt.plot(x,y1,'--',color='deepskyblue',label=ALabels[0],linewidth=LW)
#    if PlotNot[0]==1:
#        plt.fill_between(x,(y1-dy1m),(y1+dy1p),color='deepskyblue',alpha=0.3,linewidth=0.01)
#    yKP=OffShell(x,KPN,KPx0,KPx1)
#    yKPErrStat=deltaOffShell(x,KPN,KPdN,KPx0,KPdx0,KPx1,cov)
#    yKPErrSys=deltaOffShell(x,KPN,KPdNsys,KPx0,KPdx0sys,KPx1,covsys)
#    yKPErr=np.sqrt(yKPErrStat**2+yKPErrSys**2)
#    plt.plot(x,yKP,'-',color='pink',linewidth=LW,label='Kulagin-Petti')
#    plt.fill_between(x,(yKP-yKPErr),(yKP+yKPErr),color='pink',alpha=0.7,linewidth=0.01)
#    #plt.plot(x,y1,'-',color='royalblue',label=ALabels[0],linewidth=LW)
#    tempx=np.linspace(0,1,50)
#    tempy1=0.0+0*tempx
#    #plt.plot(tempx,tempy1,'.',color='black',linewidth=0.6)
#    plt.plot(tempx,tempy1,'.',color='black',ms=2)
#    plt.ylim([ymin,ymax])
#    plt.xlabel('$x$',fontsize='10')
#    plt.ylabel('$\delta f^N$',fontsize='10')

def OffShellPlot(PlotNot,Outlists,OutFileNames,Parlists,ParFileNames,filenames,ymin,ymax,ALabels):
    x=np.linspace(0,1,100)
    N1=Nval(OutFileNames[0],Outlists,OutFileNames)[0]
    x01=x0val(OutFileNames[0],Outlists,OutFileNames)[0]
    x11=x1val(OutFileNames[0],Outlists,OutFileNames)
    Nvalues=[]
    x0values=[]
    x1values=[]
    x2values=[]
    for i in range(1,len(ParFileNames)):
        Nvalues.append(ioffs(Parlists,ParFileNames,ParFileNames[i],'off1')[0])
        x0values.append(ioffs(Parlists,ParFileNames,ParFileNames[i],'off2')[0])
        x1values.append(ioffs(Parlists,ParFileNames,ParFileNames[i],'off3')[0])
        if lineoff4(ParFileNames[i])==3:
            x2values.append(ioffs(Parlists,ParFileNames,ParFileNames[i],'off4')[0])
        else:
            x2values.append(0)
    line_colors=cmap(np.linspace(0,1,12))
    KPN=8.1
    KPdN=0.3
    KPdNsys=0.526
    KPx0=0.448
    KPdx0=0.005
    KPdx0sys=0.0067
    KPx1=0.050
    cov=-0.067
    covsys=0
    y1=OffShell(x,N1,x01,x11)
    dy1p=deltaOSPlus(filenames[0],1)
    dy1m=deltaOSMinus(filenames[0],1)
    for i in range(0,len(ParFileNames)-1):
        if lineoff4(ParFileNames[i+1])==2:
            templist=[]
            for line in open(ParFileNames[i+1]):
                templist += [line.split()]
                #parline(ParFileNames[i],'10914')[0]
            ioffline=parline(ParFileNames[i+1],'inuke')[0]
            ioff9check=templist[ioffline][1][2]
            print ioff9check
#            #print templist[58]    
            if ioff9check=='9':
                temp1=OffShell(x,Nvalues[i],x0values[i],x1values[i])
                #print ioff9check
            else:
                temp1=OffShelliOff12(x,Nvalues[i],x0values[i],x1values[i])
                print ioff9check
            #temp1=OffShelliOff12(x,Nvalues[i],x0values[i],x1values[i])    
            #iOff12
            #print lineoff4(ParFileNames[i+1])
            temp1yp=deltaOSPlus(filenames[i+1],2)
            temp1ym=deltaOSMinus(filenames[i+1],2)
            plt.plot(x,temp1,'-',label=ALabels[i],linewidth=LW,c=line_colors[i])
            if PlotNot[i+1]==1:                
                plt.fill_between(x,(temp1-temp1ym),(temp1+temp1yp),color=line_colors[i],alpha=0.5,linewidth=0.01)
        else:
            temp2=OffShelliOff13(x,Nvalues[i],x0values[i],x1values[i],x2values[i])
            #print lineoff4(ParFileNames[i+1])
            #plt.plot(x,temp2,'-',label=ALabels[i],linewidth=LW,c=line_colors[i])
            temp2yp=deltaOSPlus(filenames[i+1],2)
            temp2ym=deltaOSMinus(filenames[i+1],2)
            plt.plot(x,temp2,'-',label=ALabels[i],linewidth=LW,c=line_colors[i])
            if PlotNot[i+1]==1:
                plt.fill_between(x,(temp2-temp2ym),(temp2+temp2yp),color=line_colors[i],alpha=0.5,linewidth=0.01)
    plt.plot(x,y1,'--',color='deepskyblue',label=ALabels[0],linewidth=LW)
    if PlotNot[0]==1:
        plt.fill_between(x,(y1-dy1m),(y1+dy1p),color='deepskyblue',alpha=0.3,linewidth=0.01)
    yKP=OffShell(x,KPN,KPx0,KPx1)
    yKPErrStat=deltaOffShell(x,KPN,KPdN,KPx0,KPdx0,KPx1,cov)
    yKPErrSys=deltaOffShell(x,KPN,KPdNsys,KPx0,KPdx0sys,KPx1,covsys)
    yKPErr=np.sqrt(yKPErrStat**2+yKPErrSys**2)
    plt.plot(x,yKP,'-',color='pink',linewidth=LW,label='Kulagin-Petti')
    plt.fill_between(x,(yKP-yKPErr),(yKP+yKPErr),color='pink',alpha=0.7,linewidth=0.01)
    #plt.plot(x,y1,'-',color='royalblue',label=ALabels[0],linewidth=LW)
    tempx=np.linspace(0,1,50)
    tempy1=0.0+0*tempx
    #plt.plot(tempx,tempy1,'.',color='black',linewidth=0.6)
    plt.plot(tempx,tempy1,'.',color='black',ms=2)
    plt.ylim([ymin,ymax])
    plt.xlabel('$x$',fontsize='10')
    plt.ylabel('$\delta f^N$',fontsize='10')


#pink

def RelativeRatioPlot(lists,DataFileNames,TN,CN,ymin,ymax,ALabels):
    x0=PlotData(lists,lists[0],DataFileNames[0],TN,CN)[0]
    y1=PlotData(lists,lists[0],DataFileNames[0],TN,CN)[1]
    y1Err=PlotData(lists,lists[0],DataFileNames[0],TN,CN)[2]
    RE1=y1Err/y1
    RRE1=RE1/RE1
    line_colors = cmap(np.linspace(0,1,12))
    for i in range(1,len(DataFileNames)):
        tempy=PlotData(lists,lists[i],DataFileNames[i],TN,CN)[1]
        tempyErr=PlotData(lists,lists[i],DataFileNames[i],TN,CN)[2]
        temprelErr=tempyErr/tempy
        tempRRE=temprelErr/RE1
        plt.plot(x0,tempRRE,'-',label=ALabels[i],linewidth=LW,c=line_colors[i-1])
    plt.plot(x0,RRE1,'--',color='deepskyblue',label=ALabels[0],linewidth=LW)
    plt.ylim([ymin,ymax])
    plt.xlabel(xlabel[0],fontsize='10')
    plt.ylabel('$(\delta\phi/\phi)/(\delta\phi/\phi)_{ref}$',fontsize='10')
    plt.text(0.05,ymin+(ymax-ymin)*0.85,pdflabel[CN],fontsize='10')
    plt.text(0.05,ymin+(ymax-ymin)*0.1,tablelabel[TN],fontsize='8')
    

def Legendplot(Outlists,DataFileNames,OutFileNames,alignment,ALabels):
    tempx=np.linspace(0,1,50)
    tempy1=2.1+0*tempx
    tempyKP=3.0+0*tempx
    plt.plot(tempx,tempyKP,'-',color='pink',label='Kulagin-Petti',linewidth=LW)
    plt.plot(tempx,tempy1,'-',color='deepskyblue',label=ALabels[0],linewidth=LW)
    line_colors = cmap(np.linspace(0,1,12))
    for i in range(1,len(DataFileNames)):
        tempy=2.1+0.1*i+0*tempx
        plt.plot(tempx,tempy,'-',c=line_colors[i-1],label=ALabels[i]+str(' (')+((str(TOTALchi2(Outlists,OutFileNames,OutFileNames[i]))))+str(') '),linewidth=LW)
    plt.legend(bbox_to_anchor=(alignment, 1), loc=1, borderaxespad=0.01, fontsize=7)
    plt.ylim([0,1])
    plt.axis('off')
       
# plt.plot(tempx,tempy,'-',c=line_colors[i],label=ALabels[i]+ str(':chi2=')+(str(TOTALchi2(Outlists,OutFileNames,OutFileNames[i]))),linewidth=LW)
#def duPlotsNew(TN,alignment):
#    plt.subplot(2,2,1)
#    IndividualPlot(TN,5,0.85,1.05)
#    plt.subplot(2,2,2)
#    OffShellPlot(-1.0,0.6)
#    plt.subplot(2,2,3)
#    RelativeRatioPlot(TN,5,0.4,2.5)
#    plt.subplot(2,2,4)
#    Legendplot(alignment)  
#    
#def ugPlotsNew(TN):
#    plt.subplot(2,2,1)
#    IndividualPlot(TN,9,0.98,1.02)
#    plt.subplot(2,2,2)
#    IndividualPlot(TN,11,0.95,1.5)
#    plt.subplot(2,2,3)
#    RelativeRatioPlot(TN,9,0.98,1.3)
#    plt.subplot(2,2,4)
#    RelativeRatioPlot(TN,11,0.9,1.4)        

#f1=plt.figure(1)
#f1.subplots_adjust(hspace=0.4, wspace=0.4)
#duPlotsNew(2,1.0)
#
#f2=plt.figure(2)
#f2.subplots_adjust(hspace=0.4, wspace=0.4)
#ugPlotsNew(2)    
#
#
#import matplotlib.backends.backend_pdf
#pdf = matplotlib.backends.backend_pdf.PdfPages("./Resultsioff1/iOff1-Set1.pdf")
#for fig in xrange(1,3): ## will open an empty extra figure :(
#    pdf.savefig( fig )
#pdf.close()

#This produce only the Off-Shell function
#f10=plt.figure(1)
#f10.subplots_adjust(hspace=0.4, wspace=0.4)
#OnlyOffShellPlots(1,2,3,4,5,6,7,2,ALabels,1.7)
#f10.savefig('./ResultsOffShell/OnlyOffShell_Fits-Alberto.pdf')


      
