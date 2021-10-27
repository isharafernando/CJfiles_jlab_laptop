# -*- coding: utf-8 -*-
"""
Created on Mon Jun 17 11:05:15 2019

@author: isharafernando
"""

import numpy as np
import matplotlib     
matplotlib.rc('xtick', labelsize=7)     
matplotlib.rc('ytick', labelsize=7)
import glob

centralvaluespath_Mul='./Data/Wreal_Mul'
centralvaluespath_Add='./Data/Wreal_Add'
centralvaluespath_Mul_pDn='./Data/Wreal_Mul_pDn'
centralvaluespath_Add_pDn='./Data/Wreal_Add_pDn'



#REMEMBER TO change HTpforSingleX and HTnforSingleX variable names for pDn case
#mainfitlabel='CJ15'
mainfitname_mul='Wreal_Mul'
mainfitname_add='Wreal_Add'
mainfitname_mul_pDn='Wreal_Mul_pDn'
mainfitname_add_pDn='Wreal_Add_pDn'


import matplotlib.pyplot as plt
#import numpy as np

#Number_out_files = glob.glob(datapath+'/*.out')
#Number_par_files = glob.glob(datapath+'/*.par')
#Number_pdf_files = glob.glob(datapath+'/*.pdf')

#Pdf_file_CJ15 =(centralvaluespath+'/CJ15.pdf')
#CJ15_outfile=(centralvaluespath+'/CJ15.out')

#CJ15HTA_pdf=('./CJ15HTAdata/CJ15-HTA-fit0.pdf')
#CJ15HTA_out=('./CJ15HTAdata/CJ15-HTA-fit0.out')
#CJ15HTA_par=('./CJ15HTAdata/CJ15-HTA-fit0.par')


CJ15ioffXMul_pdf=(centralvaluespath_Mul+'/'+mainfitname_mul+'.pdf')
CJ15ioffXMul_out=(centralvaluespath_Mul+'/'+mainfitname_mul+'.out')
CJ15ioffXMul_par=(centralvaluespath_Mul+'/'+mainfitname_mul+'.par')

CJ15ioffXAdd_pdf=(centralvaluespath_Add+'/'+mainfitname_add+'.pdf')
CJ15ioffXAdd_out=(centralvaluespath_Add+'/'+mainfitname_add+'.out')
CJ15ioffXAdd_par=(centralvaluespath_Add+'/'+mainfitname_add+'.par')


CJ15ioffXMul_pDn_pdf=(centralvaluespath_Mul_pDn+'/'+mainfitname_mul_pDn+'.pdf')
CJ15ioffXMul_pDn_out=(centralvaluespath_Mul_pDn+'/'+mainfitname_mul_pDn+'.out')
CJ15ioffXMul_pDn_par=(centralvaluespath_Mul_pDn+'/'+mainfitname_mul_pDn+'.par')

CJ15ioffXAdd_pDn_pdf=(centralvaluespath_Add_pDn+'/'+mainfitname_add_pDn+'.pdf')
CJ15ioffXAdd_pDn_out=(centralvaluespath_Add_pDn+'/'+mainfitname_add_pDn+'.out')
CJ15ioffXAdd_pDn_par=(centralvaluespath_Add_pDn+'/'+mainfitname_add_pDn+'.par')



#print len(Number_out_files)
#print Number_out_files[0]

##################################### extracting chi2 values #################
def parline(filename,word):
    temparray2=[]
    with open(filename) as myFile:   
        for num, line in enumerate(myFile, 0):
            if word in line:
                temparray2.append(num)          
    return temparray2

#print(parline(Number_out_files[1],'TOLERANCE IS'))   
# For initiall chi2
def InChi2parametervalue(filename,val):
    templistpar=[]
    for line in open(filename):
        templistpar += [line.split()] 
    tempNN=templistpar[parline(filename,val)[0]][5]
    #For the uncertainty value one can obtain the [2] of the above line instead of [1]    
    return float(tempNN)    

#print(Number_out_files[0])    
#print(InChi2parametervalue(Number_out_files[0],'INITIAL'))
    
def parametervalue(filename,val):
    templistpar=[]
    for line in open(filename):
        templistpar += [line.split()] 
    tempNN=templistpar[parline(filename,val)[0]][1]
    #For the uncertainty value one can obtain the [2] of the above line instead of [1]    
    return float(tempNN)
    
#print parametervalue(Number_out_files[0],'TOTAL')
    
def InitialChi2values(filenames_array):
    tempChi2=[]
    tempx=[]
    #ll=len(filenames_array)
    #tempChi2=InChi2parametervalue(filenames_array[0],'INITIAL')
    for i in range(0,len(filenames_array)):
        tempx.append(i)
        tempChi2.append(InChi2parametervalue(filenames_array[i],'INITIAL'))
    return tempx,tempChi2

#print InitialChi2values(Number_out_files)    

def Chi2values(filenames_array):
    tempChi2=[]
    tempx=[]
    #ll=len(filenames_array)
    for i in range(0,len(filenames_array)):
        tempx.append(i)
        tempChi2.append(parametervalue(filenames_array[i],'TOTAL+norm'))
    return tempx,tempChi2

#print Chi2values(Number_out_files)


def InitialChi2Plot(filenames_array):
    tempx=InitialChi2values(filenames_array)[0]
    tempy=InitialChi2values(filenames_array)[1]
    sigma2=np.var(tempy)
    sigma=np.sqrt(sigma2)
    mean=np.mean(tempy)
    tempxx=np.linspace(0,1,len(tempx))
    meanarray=mean + 0*tempxx
    #print(tempy)
    print(sigma)
    #print(meanarray)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.xlim(0,max(tempx))
    plt.plot(tempx,tempy,'g-')
    plt.plot(tempx,meanarray,'--',color='darkgreen')
    #plt.fill_between(tempx,mean-sigma,mean+sigma,color='green',alpha=0.1,linewidth=0.01)
    #plt.ylabel('Initial $\chi^2$', fontsize=12)
    #plt.xlabel(" '.out' file number ", fontsize=12)

def Chi2Plot(filenames_array):
    tempx=Chi2values(filenames_array)[0]
    tempy=Chi2values(filenames_array)[1]
    sigma2=np.var(tempy)
    sigma=np.sqrt(sigma2)
    mean=np.mean(tempy)
    tempxx=np.linspace(0,1,len(tempx))
    meanarray=mean + 0*tempxx
    #print(tempy)
    print(sigma)
    #print(meanarray)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.xlim(0,max(tempx))
    plt.plot(tempx,tempy,'b-')
    plt.plot(tempx,meanarray,'--',color='darkblue')
    plt.fill_between(tempx,mean-sigma,mean+sigma,color='blue',alpha=0.3,linewidth=0.01)
    #plt.ylabel('Total $\chi^2$', fontsize=12)
    #plt.xlabel(" '.out' file number ", fontsize=12)
    #plt.ylim(4685,4715)

################################################################


################## Ploting C(x) function ################################

def HTfunction(x,h1,h2,h3,h4):
    #tempfunctionHT=(h1*x**h2)*(1+x*h3+x*h4**2)
    tempfunctionHT=(h1*x**h2)*(1+x*h3)*((1-x)**h4)
    #print(h1)
    return tempfunctionHT


def HTpforSingleX(x,filename):
    h1=parametervalue(filename,'ht1')
    h2=parametervalue(filename,'ht2')
    h3=parametervalue(filename,'ht3')
    h4=parametervalue(filename,'ht4')
    HTpTemp=HTfunction(x,h1,h2,h3,h4)
    return HTpTemp

#print(HTpforSingleX(1,'./CJ15Calc/NewFits/CJ15fit/par_CJ15pDn/CJ15pDn_00.par'))

def HTnforSingleX(x,filename):
    h1n=parametervalue(filename,'ht1')
    h2n=parametervalue(filename,'ht2')
    h3n=parametervalue(filename,'ht3')
    h4n=parametervalue(filename,'ht4')
    HTnTemp=HTfunction(x,h1n,h2n,h3n,h4n)
    return HTnTemp

#print(HTnforSingleX(1,'./CJ15Calc/NewFits/CJ15fit/par_CJ15pDn/CJ15pDn_00.par'))

def deltaHTpPlus(filename,xp):
    test_files = glob.glob(centralvaluespath_Mul+'/par_'+filename+'/*.par')
    nn=len(test_files)
    nnn=(nn-1)/2
    x=xp
    tempvvv=[]
    for i in x:
        tempsum=0
        for j in range(1,nnn+1):
            tempvvv01=HTpforSingleX(i,test_files[2*j-1])-HTpforSingleX(i,test_files[0])
            tempvvv02=HTpforSingleX(i,test_files[2*j])-HTpforSingleX(i,test_files[0])
            tempvvv1=max(tempvvv01,tempvvv02,0)
            tempvvv2=np.square(tempvvv1)
            tempsum=tempsum+tempvvv2
        tempvvv.append(np.sqrt(tempsum))
    return tempvvv

def deltaHTnPlus(filename,xp):
    test_files = glob.glob(centralvaluespath_Mul+'/par_'+filename+'/*.par')
    #print test_files[0]
    nn=len(test_files)
    nnn=(nn-1)/2
    #x=np.linspace(0,1,100)
    x=xp
    tempvvv=[]
    for i in x:
        tempsum=0
        for j in range(1,nnn+1):
            tempvvv01=HTnforSingleX(i,test_files[2*j-1])-HTnforSingleX(i,test_files[0])
            tempvvv02=HTnforSingleX(i,test_files[2*j])-HTnforSingleX(i,test_files[0])
            tempvvv1=max(tempvvv01,tempvvv02,0)
            tempvvv2=np.square(tempvvv1)
            tempsum=tempsum+tempvvv2
        tempvvv.append(np.sqrt(tempsum))
    return tempvvv

def deltaHTnMinus(filename,xp):
    test_files = glob.glob(centralvaluespath_Mul+'/par_'+filename+'/*.par')
    #print test_files[0]
    nn=len(test_files)
    nnn=(nn-1)/2
    #x=np.linspace(0,1,100)
    x=xp
    tempvvv=[]
    for i in x:
        tempsum=0
        for j in range(1,nnn+1):
            tempvvv01=HTnforSingleX(i,test_files[0])-HTnforSingleX(i,test_files[2*j-1])
            tempvvv02=HTnforSingleX(i,test_files[0])-HTnforSingleX(i,test_files[2*j])
            tempvvv1=max(tempvvv01,tempvvv02,0)
            tempvvv2=np.square(tempvvv1)
            tempsum=tempsum+tempvvv2
        tempvvv.append(np.sqrt(tempsum))
    return tempvvv



def deltaHTpMinus(filename,xp):
    test_files = glob.glob(centralvaluespath_Mul+'/par_'+filename+'/*.par')
    nn=len(test_files)
    nnn=(nn-1)/2
    x=xp
    tempvvv=[]
    for i in x:
        tempsum=0
        for j in range(1,nnn+1):
            tempvvv01=HTpforSingleX(i,test_files[0])-HTpforSingleX(i,test_files[2*j-1])
            tempvvv02=HTpforSingleX(i,test_files[0])-HTpforSingleX(i,test_files[2*j])
            tempvvv1=max(tempvvv01,tempvvv02,0)
            tempvvv2=np.square(tempvvv1)
            tempsum=tempsum+tempvvv2
        tempvvv.append(np.sqrt(tempsum))
    return tempvvv

def HTpCentralVal(path,parfilename,xp):
    x=xp
    tempvvv=[]
    for i in x:
        tempvvv01=HTpforSingleX(i,path+'/par_'+parfilename+'/'+parfilename+'_00.par')
        tempvvv.append(tempvvv01)
    return tempvvv
    
def HTnCentralVal(path,parfilename,xp):
    x=xp
    tempvvv=[]
    for i in x:
        tempvvv01=HTnforSingleX(i,path+'/par_'+parfilename+'/'+parfilename+'_00.par')
        tempvvv.append(tempvvv01)
    return tempvvv
    
     

def CpBand(path,parfilename,listname,Q2vv,clr,lbl):
    xs=[]
    for i in range(len(listname)):
        if listname[i][1]==Q2vv:
            xs.append(float(listname[i][0]))
    HvalueUP=np.array(deltaHTpPlus(parfilename,xs))
    HvalueDOWN=np.array(deltaHTpMinus(parfilename,xs))
    Hcentral=np.array(HTpCentralVal(path,parfilename,xs))
    plt.fill_between(xs,(Hcentral-HvalueDOWN),(Hcentral+HvalueUP),color=clr,alpha=0.15,linewidth=0.01)
    plt.plot(xs,Hcentral,linewidth=0.5,color=clr,label=lbl)
    

#CpBand('./CJ15Calc/NewFits/CJ15fit','CJ15pDn',CJ15_HT_P_LT,'2.000','blue')

def CnBand(path,parfilename,listname,Q2vv,clr,lbl):
    xs=[]
    for i in range(len(listname)):
        if listname[i][1]==Q2vv:
            xs.append(float(listname[i][0]))
    HvalueUP=np.array(deltaHTnPlus(parfilename,xs))
    HvalueDOWN=np.array(deltaHTnMinus(parfilename,xs))
    Hcentral=np.array(HTnCentralVal(path,parfilename,xs))
    plt.fill_between(xs,(Hcentral-HvalueDOWN),(Hcentral+HvalueUP),color=clr,alpha=0.15,linewidth=0.01)
    plt.plot(xs,Hcentral,linewidth=0.5,color=clr,label=lbl)  
  
    
def CpforSingleFile(filename):
    x=np.linspace(0.1,1,100)
    tempvals=[HTpforSingleX(i,filename) for i in x]
    return tempvals
    
#print(parametervalue(Number_par_files[0],'ht1'))    
#print(HTpforSingleFile(Number_par_files[0]))

def CnforSingleFile(filename):
    x=np.linspace(0.1,1,100)
    tempvals=[HTnforSingleX(i,filename) for i in x]
    return tempvals

#print(HTnforSingleFile('./CJ15Calc/NewFits/CJ15fit/par_CJ15pDn/CJ15pDn_00.par'))

def CpSystmaticBand(filenames_array):
    x=np.linspace(0.1,1,100)
    #tempyUP=[]
    #tempyDOWN=[]
    plt.plot(x,CpforSingleFile(filenames_array[0]),'-r')
    for i in range(1,len(filenames_array)):
        plt.plot(x,CpforSingleFile(filenames_array[i]),'-r')
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    

def CnSystmaticBand(filenames_array):
    x=np.linspace(0.1,1,100)
    #tempyUP=[]
    #tempyDOWN=[]
    plt.plot(x,CnforSingleFile(filenames_array[0]),'-b')
    for i in range(1,len(filenames_array)):
        plt.plot(x,CnforSingleFile(filenames_array[i]),'-b')
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    


############################## H function plots ######################

def smallLIST(name):
    tempdat=open(name)
    templist=[]
    for line in tempdat:
        templist += [line.split()]
    return templist 

def datalines(filename,keyword):
    temparray1=[]
    #temparray2=[]
    with open(filename) as myFile:
        for num, line in enumerate(myFile, 0):
            if keyword in line:
                if keyword == smallLIST(filename)[num][-1]:
                    temparray1.append(smallLIST(filename)[num])
    #Templist=[LIST(filename)[i] for i in range(startline,finishline+1)]
    return temparray1


#F2_P_LT=datalines(centralvaluespath+'/'+mainfitlabel+'-P-LT.out','test_P_all')
#F2_N_LT=datalines(centralvaluespath+'/'+mainfitlabel+'-N-LT.out','test_N_all')


#HTA_F2_P_LT=datalines('./CJ15HTAdata/CJ15-HTA-P-LT.out','test_P_all')
#HTA_F2_N_LT=datalines('./CJ15HTAdata/CJ15-HTA-N-LT.out','test_N_all')

HTA_F2_P_LT=datalines(centralvaluespath_Mul+'/'+mainfitname_mul+'-F2s.out','test_P_all')
HTA_F2_N_LT=datalines(centralvaluespath_Mul+'/'+mainfitname_mul+'-F2s.out','test_N_all')

HTA_F2_DN=datalines(centralvaluespath_Add+'/'+mainfitname_add+'-F2s.out','test_D_N')
HTM_F2_DN=datalines(centralvaluespath_Mul+'/'+mainfitname_mul+'-F2s.out','test_D_N')


HTM_F2_P_LT=datalines(centralvaluespath_Add+'/'+mainfitname_add+'-F2s.out','test_P_all')
HTM_F2_N_LT=datalines(centralvaluespath_Add+'/'+mainfitname_add+'-F2s.out','test_N_all')


F2_P_mul_on=datalines(centralvaluespath_Mul+'/'+mainfitname_mul+'-F2s_onshell.out','test_P')
F2_N_mul_on=datalines(centralvaluespath_Mul+'/'+mainfitname_mul+'-F2s_onshell.out','test_N')
F2_D_mul_on=datalines(centralvaluespath_Mul+'/'+mainfitname_mul+'-F2s_onshell.out','test_D')
F2_DN_mul_on=datalines(centralvaluespath_Mul+'/'+mainfitname_mul+'-F2s_onshell.out','test_D_N')


F2_P_add_on=datalines(centralvaluespath_Add+'/'+mainfitname_add+'-F2s_onshell.out','test_P')
F2_N_add_on=datalines(centralvaluespath_Add+'/'+mainfitname_add+'-F2s_onshell.out','test_N')
F2_D_add_on=datalines(centralvaluespath_Add+'/'+mainfitname_add+'-F2s_onshell.out','test_D')
F2_DN_add_on=datalines(centralvaluespath_Add+'/'+mainfitname_add+'-F2s_onshell.out','test_D_N')


F2_P_mul_off=datalines(centralvaluespath_Mul+'/'+mainfitname_mul+'-F2s.out_NLO','test_P')
F2_N_mul_off=datalines(centralvaluespath_Mul+'/'+mainfitname_mul+'-F2s.out_NLO','test_N')
F2_D_mul_off=datalines(centralvaluespath_Mul+'/'+mainfitname_mul+'-F2s.out_NLO','test_D')
F2_DN_mul_off=datalines(centralvaluespath_Mul+'/'+mainfitname_mul+'-F2s.out_NLO','test_D_N')


F2_P_add_off=datalines(centralvaluespath_Add+'/'+mainfitname_add+'-F2s.out_NLO','test_P')
F2_N_add_off=datalines(centralvaluespath_Add+'/'+mainfitname_add+'-F2s.out_NLO','test_N')
F2_D_add_off=datalines(centralvaluespath_Add+'/'+mainfitname_add+'-F2s.out_NLO','test_D')
F2_DN_add_off=datalines(centralvaluespath_Add+'/'+mainfitname_add+'-F2s.out_NLO','test_D_N')

    
def HTpforSingleFile(filename,xp):
    x=xp
    tempvals=[HTpforSingleX(i,filename) for i in x]
    return tempvals


def HTnforSingleFile(filename,xp):
    x=xp
    tempvals=[HTnforSingleX(i,filename) for i in x]
    return tempvals    


def HpCurve(path,parfilename,listname,Q2vv,clr,lbl):
    xs=[]
    F2val=[]
    for i in range(len(listname)):
        if listname[i][1]==Q2vv:
            xs.append(float(listname[i][0]))
            F2val.append(float(listname[i][3]))
    Hcentral=np.array(HTpCentralVal(path,parfilename,xs))*np.array(F2val)
    plt.plot(xs,Hcentral,linewidth=0.5,color=clr,label=lbl)    

#HpBand('./CJ15Calc/NewFits/CJ15fit','CJ15pDn',CJ15_HT_P_LT,'5.000','red')


def HnCurve(path,parfilename,listname,Q2vv,clr,lbl):
    xs=[]
    F2val=[]
    for i in range(len(listname)):
        if listname[i][1]==Q2vv:
            xs.append(float(listname[i][0]))
            F2val.append(float(listname[i][3]))
    Hcentral=np.array(HTnCentralVal(path,parfilename,xs))*np.array(F2val)
    plt.plot(xs,Hcentral,linewidth=0.5,color=clr,label=lbl)

def HpAdditiveCurve(path,parfilename,listname,Q2vv,clr,lbl):
    xs=[]
    F2val=[]
    for i in range(len(listname)):
        if listname[i][1]==Q2vv:
            xs.append(float(listname[i][0]))
            F2val.append(float(listname[i][3]))
    Hcentral=np.array(HTpCentralVal(path,parfilename,xs))*1
    plt.plot(xs,Hcentral,color=clr,label=lbl,linewidth=1)    

#HpBand('./CJ15Calc/NewFits/CJ15fit','CJ15pDn',CJ15_HT_P_LT,'5.000','red')


def HnAdditiveCurve(path,parfilename,listname,Q2vv,clr,lbl):
    xs=[]
    F2val=[]
    for i in range(len(listname)):
        if listname[i][1]==Q2vv:
            xs.append(float(listname[i][0]))
            F2val.append(float(listname[i][3]))
    Hcentral=np.array(HTnCentralVal(path,parfilename,xs))*1
    plt.plot(xs,Hcentral,color=clr,label=lbl,linewidth=1)
    
    
def CpTildeCurve(path,parfilename,listname,Q2vv,clr,lbl):
    xs=[]
    F2val=[]
    for i in range(len(listname)):
        if listname[i][1]==Q2vv:
            xs.append(float(listname[i][0]))
            F2val.append(float(listname[i][3]))
    Hcentral=np.array(HTpCentralVal(path,parfilename,xs))/np.array(F2val)
    plt.plot(xs,Hcentral,color=clr,label=lbl,linewidth=1)    

#HpBand('./CJ15Calc/NewFits/CJ15fit','CJ15pDn',CJ15_HT_P_LT,'5.000','red')


def CnTildeCurve(path,parfilename,listname,Q2vv,clr,lbl):
    xs=[]
    F2val=[]
    for i in range(len(listname)):
        if listname[i][1]==Q2vv:
            xs.append(float(listname[i][0]))
            F2val.append(float(listname[i][3]))
    Hcentral=np.array(HTnCentralVal(path,parfilename,xs))/np.array(F2val)
    plt.plot(xs,Hcentral,color=clr,label=lbl,linewidth=1)    
        

def HpSystmaticBand(filenames_array,listname,Q2vv):
    #x=np.linspace(0.1,1,100)
    #tempyUP=[]
    #tempyDOWN=[]
    #plt.plot(x,HTpforSingleFile(filenames_array[0]),'-r',label='C_p(x)')
    xs=[]
    F2val=[]
    for i in range(len(listname)):
        if listname[i][1]==Q2vv:
            xs.append(float(listname[i][0]))
            F2val.append(float(listname[i][3]))
    plt.plot(xs,np.array(HTpforSingleFile(filenames_array[0],xs))*np.array(F2val),'-r')        
    for i in range(0,len(filenames_array)):
        plt.plot(xs,np.array(HTpforSingleFile(filenames_array[i],xs))*np.array(F2val),'-r')
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    

def HnSystmaticBand(filenames_array,listname,Q2vv):
    #x=np.linspace(0.1,1,100)
    #tempyUP=[]
    #tempyDOWN=[]
    #plt.plot(x,HTnforSingleFile(filenames_array[0]),'-b',label='C_n(x)')
    xs=[]
    F2val=[]
    for i in range(len(listname)):
        if listname[i][1]==Q2vv:
            xs.append(float(listname[i][0]))
            F2val.append(float(listname[i][3]))
    plt.plot(xs,np.array(HTnforSingleFile(filenames_array[0],xs))*np.array(F2val),'-b')                
    for i in range(0,len(filenames_array)):
        plt.plot(xs,np.array(HTnforSingleFile(filenames_array[i],xs))*np.array(F2val),'-b')
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)


def HnHpSystematicDifferenceBand(filenames_array,listnameN,listnameP,Q2vv):
    xs=[]
    F2Nval=[]
    #HnHpdiff=[]
    for i in range(len(listnameN)):
        if listnameN[i][1]==Q2vv:
            xs.append(float(listnameN[i][0]))
            F2Nval.append(float(listnameN[i][3]))
    F2Pval=[]
    for i in range(len(listnameP)):
        if listnameP[i][1]==Q2vv:
            #xs.append(float(listnameP[i][0]))
            F2Pval.append(float(listnameP[i][3]))
    Hnvalues0=np.array(HTnforSingleFile(filenames_array[0],xs))*np.array(F2Nval)        
    Hpvalues0=np.array(HTpforSingleFile(filenames_array[0],xs))*np.array(F2Pval)
    HnHpSysdiff0=Hnvalues0-Hpvalues0
    plt.plot(xs,HnHpSysdiff0,'-g')
    for i in range(0,len(filenames_array)):        
        Hnvalues=np.array(HTnforSingleFile(filenames_array[i],xs))*np.array(F2Nval)        
        Hpvalues=np.array(HTpforSingleFile(filenames_array[i],xs))*np.array(F2Pval)
        HnHpSysdiff=Hnvalues-Hpvalues
        plt.plot(xs,HnHpSysdiff,'-g')
        #print(HnHpdiff)
    #print(Hnvalues)    
    #HnHpSysdiff=HnHpdiff
    #for i in range(0,len(filenames_array)):
    #plt.plot(xs,HnHpSysdiff,'-g')    
    #plt.ylim(0,2.0)

def HpBand(path,parfilename,listname,Q2vv,clr,lbl):
    xs=[]
    F2val=[]
    for i in range(len(listname)):
        if listname[i][1]==Q2vv:
            xs.append(float(listname[i][0]))
            F2val.append(float(listname[i][3]))
    HvalueUP=np.array(deltaHTpPlus(parfilename,xs))*np.array(F2val)
    HvalueDOWN=np.array(deltaHTpMinus(parfilename,xs))*np.array(F2val)
    Hcentral=np.array(HTpCentralVal(path,parfilename,xs))*np.array(F2val)
    plt.fill_between(xs,(Hcentral-HvalueDOWN),(Hcentral+HvalueUP),color=clr,alpha=0.15,linewidth=0.01)
    plt.plot(xs,Hcentral,linewidth=0.5,color=clr,label=lbl)    

#HpBand('./CJ15Calc/NewFits/CJ15fit','CJ15pDn',CJ15_HT_P_LT,'5.000','red')


def HnBand(path,parfilename,listname,Q2vv,clr,lbl):
    xs=[]
    F2val=[]
    for i in range(len(listname)):
        if listname[i][1]==Q2vv:
            xs.append(float(listname[i][0]))
            F2val.append(float(listname[i][3]))
    HvalueUP=np.array(deltaHTnPlus(parfilename,xs))*np.array(F2val)
    HvalueDOWN=np.array(deltaHTnMinus(parfilename,xs))*np.array(F2val)
    Hcentral=np.array(HTnCentralVal(path,parfilename,xs))*np.array(F2val)
    plt.fill_between(xs,(Hcentral-HvalueDOWN),(Hcentral+HvalueUP),color=clr,alpha=0.15,linewidth=0.01)
    plt.plot(xs,Hcentral,linewidth=0.5,color=clr,label=lbl)
    
    
#HnBand('./CJ15Calc/NewFits/CJ15fit','CJ15pDn',CJ15_HT_P_LT,'5.000','blue')
def HnHpDifferenceBand(path,parfilename,listnameN,listnameP,Q2vv,clr,lbl):
    xs=[]
    F2nval=[]
    F2pval=[]
    for i in range(len(listnameP)):
        if listnameP[i][1]==Q2vv:
            xs.append(float(listnameP[i][0]))
            F2nval.append(float(listnameN[i][3]))
            F2pval.append(float(listnameP[i][3]))
    HpvalueUP=np.array(deltaHTpPlus(parfilename,xs))*np.array(F2pval)
    HpvalueDOWN=np.array(deltaHTpMinus(parfilename,xs))*np.array(F2pval)
    Hpcentral=np.array(HTpCentralVal(path,parfilename,xs))*np.array(F2pval)
    ###########
    HnvalueUP=np.array(deltaHTnPlus(parfilename,xs))*np.array(F2nval)
    HnvalueDOWN=np.array(deltaHTnMinus(parfilename,xs))*np.array(F2nval)
    Hncentral=np.array(HTnCentralVal(path,parfilename,xs))*np.array(F2nval)
    ###########
    HnHpRcentral=Hncentral-Hpcentral
    HnHpRUP=np.sqrt((HnvalueUP)**2+(HpvalueUP)**2)
    HnHpRDOWN=np.sqrt((HnvalueDOWN)**2+(HpvalueDOWN)**2)
    plt.fill_between(xs,(HnHpRcentral-HnHpRDOWN),(HnHpRcentral+HnHpRUP),color=clr,alpha=0.3,linewidth=0.01)
    plt.plot(xs,HnHpRcentral,linewidth=0.5,color=clr,label=lbl)
    #plt.ylim(0,2.0)


#d/u ratio
#LW=10
#DatFiles=[]
#lists=[]
#for i in range(0,len(Number_pdf_files)):
#    templistDat=[]
#    DatFiles.append(open(Number_pdf_files[i]))
#    for line in open(Number_pdf_files[i]):
#        templistDat += [line.split()]
#    lists.append(templistDat)

def LIST(filename):
    templistDat=[]
    lists=[]
    for line in open(filename):
        templistDat += [line.split()]
    lists.append(templistDat)
    return lists


#print(LIST(Number_pdf_files[0])[0][12])
#print([LIST(Number_pdf_files[0])[0][i][0] for i in range(6,12)]) 
    
    
def linenumbers(filename):
    temparray1=[]
    with open(filename) as myFile:   
        for num, line in enumerate(myFile, 0):
            if 'Q2' in line:
                temparray1.append(num)          
    return temparray1

def pdfSelect(filename,column,start,stop):
    TempColumn=[LIST(filename)[0][i][column] for i in range(start,stop)]
    return TempColumn

#print(pdfSelect(Number_pdf_files[0],0,10,12))
#print(len(LIST(Number_pdf_files[0])[0]))
    
def LinNumTbles(filename,Numbr):
    lins=linenumbers(filename)
    temp1=lins[:len(lins)-5]
    temp2=[lins[i]-2 for i in range(8,len(lins))]
    newlins= temp1 + temp2
    TblStrtRef=[newlins[i]+3 for i in range(0,len(newlins))]
    TblEndRefTEMP=[newlins[i+1]-3 for i in range(0,len(newlins)-1)]
    TblEndRef=TblEndRefTEMP+[len(LIST(CJ15ioffXAdd_pdf)[0])]
    temparray=[TblStrtRef[Numbr],TblEndRef[Numbr]]
    return temparray

#print(LinNumTbles(Number_pdf_files[0],0)[1])    
    
def PlotData(filename,TblNo,ColNo):
    tempx=pdfSelect(filename,0,LinNumTbles(filename,TblNo)[0],LinNumTbles(filename,TblNo)[1])
    tempy=pdfSelect(filename,ColNo,LinNumTbles(filename,TblNo)[0],LinNumTbles(filename,TblNo)[1])
    tempyErr=pdfSelect(filename,ColNo+1,LinNumTbles(filename,TblNo)[0],LinNumTbles(filename,TblNo)[1])
    tx=np.array([float(i) for i in tempx])
    ty=np.array([float(i) for i in tempy])
    tyErr=np.array([float(i) for i in tempyErr])
    return tx,ty,tyErr

#print(PlotData(Number_pdf_files[0],10,0))

def SinglePDFPlot(filename,TN,CN,clr):
    x0=PlotData(filename,TN,CN)[0]
    y1=PlotData(filename,TN,CN)[1]
    #y1Err=PlotData(filename,TN,CN)[2]
    #print(y1)
    plt.plot(x0,y1,'-',color=clr)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    
def SinglePDFPlotWithBand(filename,TN,CN,clr,lbl):
    x0=PlotData(filename,TN,CN)[0]
    y1=PlotData(filename,TN,CN)[1]
    y1Err=PlotData(filename,TN,CN)[2]
    #print(y1)
    plt.plot(x0,y1,'-',color=clr,label=lbl)
    plt.fill_between(x0,y1-y1Err,y1+y1Err,alpha=0.2,color='white',facecolor=clr)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)    

def AllPDFPlots(folder):
    SinglePDFPlotWithBand(folder[0],2,5,'blue','CJ15')
    for i in range(0,len(folder)):
        SinglePDFPlot(folder[i],2,5,'blue')



#############  Off-Shell #########################

def ioffs(filename,val):
    tempN=float(LIST(filename)[0][parline(filename,val)[0]][2])
    #tempErrN=float(LIST(filename)[0][parline(filename,val)[0]][4])
    return tempN        
#print(ioffs(Number_out_files[0],'N =')[0])

def ioffsErr(filename,val):
    tempN=float(LIST(filename)[0][parline(filename,val)[0]][4])
    #tempErrN=float(LIST(filename)[0][parline(filename,val)[0]][4])
    return tempN      
    
def OffShell(x,N,x0,x1):
    tempOffShell=N*(x-x0)*(x-x1)*(1+x0-x)
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

def OffShellSinglePlot(filename,clr,lbl=''):
    x=np.linspace(0,1,100)
    N1=ioffs(filename,'N =')
    x01=ioffs(filename,'x0 =')
    x11=ioffs(filename,'x1 =')
    y1=OffShell(x,N1,x01,x11)
    plt.plot(x,y1,'-',color=clr,label=lbl)

#OffShellSinglePlot(Number_out_files[0]) 
   
def OffShellSinglePlotBand(filename,clr,lbl):
    x=np.linspace(0,1,100)
    N1=ioffs(filename,'N =')
    N1Err=ioffsErr(filename,'N =')
    x01=ioffs(filename,'x0 =')
    x01Err=ioffsErr(filename,'x0 =')
    x11=ioffs(filename,'x1 =')
    y1=OffShell(x,N1,x01,x11)
    y1err=deltaOffShell(x,N1,N1Err,x01,x01Err,x11,0)
    plt.plot(x,y1,'-',color=clr,label=lbl,alpha=0.3)
    plt.fill_between(x,(y1-y1err),(y1+y1err),color=clr,alpha=0.3,linewidth=0.01)
    
#OffShellSinglePlotBand(CJ15_outfile,'deepskyblue','CJ15')   


def AllOSPlots(folder,clr,lbl):
    OffShellSinglePlot(folder[0],clr,lbl) 
    for i in range(1,len(folder)):
        OffShellSinglePlot(folder[i],clr)    


def OffShelliOff12(x,N,x0,x1):
    #tempOffShell=N*(x-x0)*(x-x1)
    tempOffShell= N + x0*x + x1*(x**2)
    return tempOffShell    
    
def OffShelliOff13(x,N,x0,x1,x2):
    #tempOffShell=N*(x-x0)*(x-x1)*(x-x2)
    tempOffShell= N + x0*x + x1*(x**2) + x2*(x**3)
    return tempOffShell
    
def OffShelliOff14(x,N,x0,x1,x2,x3):
    #tempOffShell=N*(x-x0)*(x-x1)*(x-x2)
    tempOffShell= N + x0*x + x1*(x**2) + x2*(x**3) + x3*(x**3)
    return tempOffShell

def OffShelliOff15(x,N,x0,x1,x2,x3,x4):
    #tempOffShell=N*(x-x0)*(x-x1)*(x-x2)
    tempOffShell= N + x0*x + x1*(x**2) + x2*(x**3) + x3*(x**3) + x4*(x**4)
    return tempOffShell        

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


def OSforSingleX(x,filename,ioffFlg):
    if ioffFlg=='8' or '9':
        Nval=float(parametervalue(filename,'off1'))
        x0val=float(parametervalue(filename,'off2'))
        x1val=float(parametervalue(filename,'off3'))
        tempOSX=OffShell(x,Nval,x0val,x1val)
    if ioffFlg=='1':
        Nval=float(parametervalue(filename,'off1'))
        x0val=float(parametervalue(filename,'off2'))
        x1val=float(parametervalue(filename,'off3'))
        Polydegree=parline(filename,'off4')
        if Polydegree==parline(filename,'off3')[0]+1:
            x2val=float(parametervalue(filename,'off4'))
            tempOSX=OffShelliOff13(x,Nval,x0val,x1val,x2val)
        else:
            tempOSX=OffShelliOff12(x,Nval,x0val,x1val)
    return tempOSX



def OSforSingleFile(filename,y):
    x=np.linspace(0,1,100)
    tempvals=[OSforSingleX(i,filename,y) for i in x]
    return tempvals

def deltaOSPlus(path,filename,y,xp):
    test_files = glob.glob(path+'/par_'+filename+'/*.par')
    #print test_files[0]
    nn=len(test_files)
    nnn=int((nn-1)/2)
    #x=np.linspace(0,1,100)
    x=xp
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

def deltaOSMinus(path,filename,y,xp):
    test_files = glob.glob(path+'/par_'+filename+'/*.par')
    #print test_files[0]
    nn=len(test_files)
    nnn=int((nn-1)/2)
    #x=np.linspace(0,1,100)
    x=xp
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




def OffShellPlot(path,parfilename,clr,lbl):
    central_par_file=(path+'/'+parfilename+'.par')
    par_files = glob.glob(path+'/par_'+parfilename+'/*.par')
    #print(par_files[0])
    x=np.linspace(0,1,100)
    ioffFlag=LIST(par_files[0])[0][parline(par_files[0],'inuke')[0]][1][2]
    #print(ioffFlag)
    if ioffFlag=='9':
        Nval=float(parametervalue(central_par_file,'off1'))
        x0val=float(parametervalue(central_par_file,'off2'))
        x1val=float(parametervalue(central_par_file,'off3'))
        temp1=OffShell(x,Nval,x0val,x1val)
        plt.plot(x,temp1,color=clr,label=lbl)
        #plt.ylim(-0.7,0.3)
        temp1yp=deltaOSPlus(path,parfilename,ioffFlag,x)
        temp1ym=deltaOSMinus(path,parfilename,ioffFlag,x)
        plt.fill_between(x,(temp1-temp1ym),(temp1+temp1yp),color=clr,alpha=0.2,linewidth=0.01)
    if ioffFlag=='1':
        Nval=float(parametervalue(central_par_file,'off1'))
        x0val=float(parametervalue(central_par_file,'off2'))
        x1val=float(parametervalue(central_par_file,'off3'))
        Polydegree=parline(central_par_file,'off4')
        if Polydegree==parline(central_par_file,'off3')[0]+1:
            x2val=float(parametervalue(central_par_file,'off4'))
            temp1=OffShelliOff13(x,Nval,x0val,x1val,x2val)
            temp1yp=deltaOSPlus(path,parfilename,ioffFlag,x)
            temp1ym=deltaOSMinus(path,parfilename,ioffFlag,x)
            plt.plot(x,temp1,color=clr,label=lbl)
            plt.fill_between(x,(temp1-temp1ym),(temp1+temp1yp),color=clr,alpha=0.2,linewidth=0.01)
        else:
            temp1=OffShelliOff12(x,Nval,x0val,x1val)
            temp1yp=deltaOSPlus(path,parfilename,ioffFlag,x)
            temp1ym=deltaOSMinus(path,parfilename,ioffFlag,x)
            plt.plot(x,temp1,color=clr,label=lbl)
            plt.fill_between(x,(temp1-temp1ym),(temp1+temp1yp),color=clr,alpha=0.2,linewidth=0.01)            
    



#OffShellPlot(centralvaluespath_Add,'CJ15pDn_i9_A')


# F2 plots
def columns(filename,column):
    templist=[]
    for line in open(filename):
        templist += [line.split()]
    TempColumn=[float(templist[i][column]) for i in range(0,len(templist))]
    return np.array(TempColumn)

def plotband(filename1,filename2,LB):
    xpt=columns(filename1,0)
    ylow=columns(filename1,1)
    yup=columns(filename2,1)
    #plt.plot(xpt,(ylow+yup)/2,color='pink',linewidth=0.5,alpha=0.9,label=LB)
    #plt.fill_between(xpt,ylow,yup,color='pink',alpha=0.7,linewidth=0.01)
    plt.plot(xpt,yup,'--',color='black',linewidth=2,label=LB)
    plt.plot(xpt,ylow,'--',color='black',linewidth=2,label='')


def F2Plot(listname,Q2vv,clr,lbl):
    xs=[]
    F2val=[]
    F2valErr=[]
    for i in range(len(listname)):
        if listname[i][1]==Q2vv:
            xs.append(float(listname[i][0]))
            F2val.append(float(listname[i][3]))
            F2valErr.append(float(listname[i][4]))
    F2central=np.array(F2val)
    F2Error=np.array(F2valErr)
    plt.fill_between(xs,(F2central-F2Error),(F2central+F2Error),color=clr,alpha=0.3,linewidth=0.01)
    plt.plot(xs,F2central,linewidth=1,color=clr,label=lbl)



def F2ComparisonPlot(listnameON,listnameOFF,Q2vv,clr,lbl):
    xs=[]
    F2val_ON=[]
    F2val_OFF=[]
    for i in range(len(listnameON)):
        if listnameON[i][1]==Q2vv:
            xs.append(float(listnameON[i][0]))
            F2val_ON.append(float(listnameON[i][3]))
            F2val_OFF.append(float(listnameOFF[i][3]))
    F2ONcentral=np.array(F2val_ON)
    F2OFFcentral=np.array(F2val_OFF)
    F2Ratio=(F2ONcentral-F2OFFcentral)/F2ONcentral
    plt.plot(xs,F2Ratio,linewidth=1,color=clr,label=lbl)



def F2ComparisonPlotFixedX(listnameON,listnameOFF,Xval,clr,lbl):
    Q2val=[]
    F2val_ON=[]
    F2val_OFF=[]
    for i in range(len(listnameON)):
        if listnameON[i][0]==Xval:
            Q2val.append(float(listnameON[i][1]))
            F2val_ON.append(float(listnameON[i][3]))
            F2val_OFF.append(float(listnameOFF[i][3]))
    F2ONcentral=np.array(F2val_ON)
    F2OFFcentral=np.array(F2val_OFF)
    F2Ratio=(F2ONcentral-F2OFFcentral)/F2ONcentral
    plt.plot(Q2val,F2Ratio,linewidth=1,color=clr,label=lbl)


#### New Edits
def F2NP_Plot_FixedQ2(listnameN,listnameP,Q2val,clr,lbl):
    Xval=[]
    F2val_N=[]
    F2val_P=[]
    F2NvalErr=[]
    F2PvalErr=[]
    for i in range(len(listnameN)):
        if listnameN[i][1]==Q2val:
            Xval.append(float(listnameN[i][0]))
            F2val_N.append(float(listnameN[i][3]))
            F2NvalErr.append(float(listnameN[i][4]))
            F2val_P.append(float(listnameP[i][3]))
            F2PvalErr.append(float(listnameP[i][4]))
    F2N=np.array(F2val_N)
    dF2N=np.array(F2NvalErr)
    F2P=np.array(F2val_P)
    dF2P=np.array(F2PvalErr)
    F2Ratio=F2N/F2P
    F2RatioError=np.sqrt((np.square(dF2N/F2N)+np.square(dF2P/F2P))*np.square(F2Ratio))
    plt.fill_between(Xval,(F2Ratio-F2RatioError),(F2Ratio+F2RatioError),color=clr,alpha=0.3,linewidth=0.01)
    plt.plot(Xval,F2Ratio,linewidth=1,color=clr,label=lbl)



def plotpoints(filename,clr,LB,shift):
    xpt=columns(filename,0)
    ypt=columns(filename,1)-shift
    plt.plot(xpt,ypt,'--',color=clr,label=LB,linewidth=2)


#f3=plt.figure(3)
#plt.xlim(0.5,1.0)
#plt.ylim(0.0,0.25)
#SinglePDFPlotWithBand(CJ15ioffXMul_pdf,2,5,'green','Multiplicative_HT($p = n$)')
#SinglePDFPlotWithBand(CJ15ioffXAdd_pdf,2,5,'orange','Additive_HT($p = n$)')
#plt.ylabel('$d/u$', fontsize=25)
#plt.xlabel('$x$', fontsize=25)
#plt.xticks(fontsize=16)
#plt.yticks(fontsize=16)
#plt.legend(loc=1, fontsize=16,frameon=False, handlelength=3)
#f3.savefig('./Results/yWyLyB/pEQn_PDFs.pdf')
#
#
#
#f4=plt.figure(4)
#CpBand(centralvaluespath_Mul,mainfitname_mul,HTM_F2_P_LT,'5.000','red','$C_{(p)}$')
#CnBand(centralvaluespath_Mul,mainfitname_mul,HTM_F2_N_LT,'5.000','blue','$C_{(n)}$')
#CpTildeCurve(centralvaluespath_Add,mainfitname_add,HTM_F2_P_LT,'5.000','orange','$\\widetilde{C}_{(p)}$(from Additive)')
#CnTildeCurve(centralvaluespath_Add,mainfitname_add,HTM_F2_N_LT,'5.000','green','$\\widetilde{C}_{(n)}$(from Additive)')
#plt.ylabel('$C(x)$', fontsize=25)
#plt.xlabel("$x$", fontsize=25)
#plt.ylim(-0.2,1.5)
#plt.xticks(fontsize=16)
#plt.yticks(fontsize=16)
#plt.legend(loc=2, fontsize=16,frameon=False, handlelength=3)
#f4.savefig('./Results/yWyLyB/pEQn_Cs.pdf')
#
#
#
#f6=plt.figure(6)
##plt.suptitle('$H_p(x)$, $H_n(x)$ and $\\widetilde{H}_{(p)}$, $\\widetilde{H}_{(n)}$ vs $x$', fontsize=16)
#HpBand(centralvaluespath_Mul,mainfitname_mul,HTM_F2_P_LT,'5.000','red','$\\widetilde{H}_{(p)}$')
#HnBand(centralvaluespath_Mul,mainfitname_mul,HTM_F2_N_LT,'5.000','blue','$\\widetilde{H}_{(n)}$')
#HpAdditiveCurve(centralvaluespath_Add,mainfitname_add,HTA_F2_P_LT,'5.000','orange','$H_{(p)}$(Additive)')
#HnAdditiveCurve(centralvaluespath_Add,mainfitname_add,HTA_F2_N_LT,'5.000','green','$H_{(n)}$(Additive)')
#plt.ylabel('$H(x)$', fontsize=20)
#plt.xlabel("$x$", fontsize=25)
#plt.ylim(-0.06,0.04)
#plt.xticks(fontsize=16)
#plt.yticks(fontsize=16)
#plt.legend(loc=2, fontsize=16,frameon=False, handlelength=3)
#f6.savefig('./Results/yWyLyB/pEQn_Hs.pdf')
#
#
#
#f9=plt.figure(9)
#plt.suptitle('$\delta f$ vs $x$', fontsize=16)
#OffShellPlot(centralvaluespath_Mul,mainfitname_mul,'green','Multiplicative_HT($p = n$)')
#OffShellPlot(centralvaluespath_Add,mainfitname_add,'orange','Additive_HT($p = n$)')
#plotpoints('./Data/AKP-OS-UpperBound.txt','black','AKP (2017)',0.0)
#plotpoints('./Data/AKP-OS-LowerBound.txt','black','',0.0)
#plt.legend(loc=3, fontsize=16,frameon=False, handlelength=3)
#plt.ylim(-0.7,0.3)
#plt.xticks(fontsize=16)
#plt.yticks(fontsize=16)
#plt.xlabel("$x$", fontsize=25)
#plt.ylabel("$\delta f$", fontsize=25)
#f9.savefig('./Results/yWyLyB/pEQn_OSs.pdf')
##
#
#
#f10=plt.figure(10)
##plt.suptitle('$F_2^D/F_2^N vs $x$', fontsize=16)
#F2Plot(HTM_F2_DN,'5.000','green','Multiplicative_HT($p = n$)')
#F2Plot(HTA_F2_DN,'5.000','orange','Additive_HT($p = n$)')
#plotband('./Data/AKP-F2DN-Upper.txt','./Data/AKP-F2DN-Lower.txt','AKP (2017)')
#plt.legend(loc=3, fontsize=16,frameon=False, handlelength=3)
#plt.ylim(0.9,1.07)
#plt.xlim(0.1,0.9)
#plt.xticks(fontsize=16)
#plt.yticks(fontsize=16)
#plt.xlabel("$x$", fontsize=25)
#plt.ylabel("$F_2^D/F_2^N$", fontsize=25)
#f10.savefig('./Results/yWyLyB/pEQn_F2DNs.pdf')
#
#
#
#
#f11=plt.figure(11)
##plt.suptitle('$F_2^D/F_2^N vs $x$', fontsize=16)
#F2Plot(F2_DN_mul_off,'5.000','green','Multiplicative_HT($p = n$)')
#F2Plot(F2_DN_add_off,'5.000','orange','Additive_HT($p = n$)')
#plotband('./Data/AKP-F2DN-Upper.txt','./Data/AKP-F2DN-Lower.txt','AKP (2017)')
#plt.legend(loc=3, fontsize=16,frameon=False, handlelength=3)
#plt.ylim(0.9,1.07)
#plt.xlim(0.1,0.9)
#plt.xticks(fontsize=16)
#plt.yticks(fontsize=16)
#plt.xlabel("$x$", fontsize=25)
#plt.ylabel("$F_2^D/F_2^N$", fontsize=25)
#f11.savefig('./Results/yWyLyB/pEQn_F2DNs_errorbands.pdf')


#F2ComparisonPlot(F2_D_mul_on,F2_D_mul_off,'2.000','green','Multiplicative_HT 2($p = n$)')
#F2ComparisonPlot(F2_D_add_on,F2_D_add_off,'2.000','orange','Additive_HT 2($p = n$)')
#F2ComparisonPlot(F2_D_mul_on,F2_D_mul_off,'5.000','blue','Multiplicative_HT 5($p = n$)')
#F2ComparisonPlot(F2_D_add_on,F2_D_add_off,'5.000','red','Additive_HT 5($p = n$)')



#f12=plt.figure(12, figsize=(10, 8))
#plt.suptitle('$F_2^D$((OnShell-OffShell)/OnShell) vs $x$ (Multiplicative HT $p = n$)', fontsize=16)
#F2ComparisonPlot(F2_D_mul_on,F2_D_mul_off,'2.000','yellow','$Q^2=2GeV^2$($p = n$)')
#F2ComparisonPlot(F2_D_mul_on,F2_D_mul_off,'5.000','orange','$Q^2=5GeV^2$($p = n$)')
#F2ComparisonPlot(F2_D_mul_on,F2_D_mul_off,'10.000','red','$Q^2=10GeV^2$($p = n$)')
#F2ComparisonPlot(F2_D_mul_on,F2_D_mul_off,'20.000','darkred','$Q^2=20GeV^2$($p = n$)')
#F2ComparisonPlot(F2_D_mul_on,F2_D_mul_off,'50.000','brown','$Q^2=50GeV^2$($p = n$)')
#F2ComparisonPlot(F2_D_mul_on,F2_D_mul_off,'100.000','black','$Q^2=100GeV^2$($p = n$)')
#plt.legend(loc=3, fontsize=16,frameon=False, handlelength=3)
#plt.ylim(-0.20,0.20)
##plt.xlim(0.1,0.9)
#plt.xticks(fontsize=16)
#plt.yticks(fontsize=16)
#plt.xlabel("$x$", fontsize=25)
#plt.ylabel("$F_2^D$((OnShell-OffShell)/OnShell)", fontsize=25)
#f12.savefig('./Results/yWyLyB/pEQn_F2Ds_Off_vs_On_Mul.pdf')
#
#
#
#f13=plt.figure(13, figsize=(10, 8))
#plt.suptitle('$F_2^D$((OnShell-OffShell)/OnShell) vs $x$ (Additive HT $p = n$)', fontsize=16)
#F2ComparisonPlot(F2_D_add_on,F2_D_add_off,'2.000','greenyellow','$Q^2=2GeV^2$($p = n$)')
#F2ComparisonPlot(F2_D_add_on,F2_D_add_off,'5.000','lawngreen','$Q^2=5GeV^2$($p = n$)')
#F2ComparisonPlot(F2_D_add_on,F2_D_add_off,'10.000','limegreen','$Q^2=10GeV^2$($p = n$)')
#F2ComparisonPlot(F2_D_add_on,F2_D_add_off,'20.000','mediumseagreen','$Q^2=20GeV^2$($p = n$)')
#F2ComparisonPlot(F2_D_add_on,F2_D_add_off,'50.000','green','$Q^2=50GeV^2$($p = n$)')
#F2ComparisonPlot(F2_D_add_on,F2_D_add_off,'100.000','darkgreen','$Q^2=100GeV^2$($p = n$)')
#plt.legend(loc=3, fontsize=16,frameon=False, handlelength=3)
#plt.ylim(-0.20,0.20)
##plt.xlim(0.1,0.9)
#plt.xticks(fontsize=16)
#plt.yticks(fontsize=16)
#plt.xlabel("$x$", fontsize=25)
#plt.ylabel("$F_2^D$((OnShell-OffShell)/OnShell)", fontsize=25)
#f13.savefig('./Results/yWyLyB/pEQn_F2Ds_Off_vs_On_Add.pdf')



#f14=plt.figure(14, figsize=(12, 10))
#plt.suptitle('$F_2^D$((OnShell-OffShell)/OnShell) vs $x$ (Multiplicative HT $p = n$ at fixed x)', fontsize=16)
#F2ComparisonPlotFixedX(F2_D_mul_on,F2_D_mul_off,'.2000E+00','yellow','$x=0.2$($p = n$)')
#F2ComparisonPlotFixedX(F2_D_mul_on,F2_D_mul_off,'.4000E+00','orange','$x=0.4$($p = n$)')
#F2ComparisonPlotFixedX(F2_D_mul_on,F2_D_mul_off,'.6000E+00','red','$x=0.6$($p = n$)')
#F2ComparisonPlotFixedX(F2_D_mul_on,F2_D_mul_off,'.7000E+00','darkred','$x=0.7$($p = n$)')
#F2ComparisonPlotFixedX(F2_D_mul_on,F2_D_mul_off,'.8000E+00','brown','$x=0.8$($p = n$)')
#plt.legend(loc=3, fontsize=16,frameon=False, handlelength=3)
#plt.ylim(-0.05,0.05)
##plt.xlim(0.1,0.9)
#plt.xticks(fontsize=16)
#plt.yticks(fontsize=16)
#plt.xlabel("$Q^2$", fontsize=25)
#plt.ylabel("$F_2^D$((OnShell-OffShell)/OnShell)", fontsize=25)
#f14.savefig('./Results/pEQn_F2Ds_Off_vs_On_Mul_fixd_x.pdf')
#
#
#
#f15=plt.figure(15, figsize=(12, 10))
#plt.suptitle('$F_2^D$((OnShell-OffShell)/OnShell) vs $x$ (Additive HT $p = n$ ar fixed x)', fontsize=16)
#F2ComparisonPlotFixedX(F2_D_add_on,F2_D_add_off,'.2000E+00','greenyellow','$x=0.2$($p = n$)')
#F2ComparisonPlotFixedX(F2_D_add_on,F2_D_add_off,'.4000E+00','lawngreen','$x=0.4$($p = n$)')
#F2ComparisonPlotFixedX(F2_D_add_on,F2_D_add_off,'.6000E+00','limegreen','$x=0.6$($p = n$)')
#F2ComparisonPlotFixedX(F2_D_add_on,F2_D_add_off,'.7000E+00','mediumseagreen','$x=0.7$($p = n$)')
#F2ComparisonPlotFixedX(F2_D_add_on,F2_D_add_off,'.8000E+00','green','$x=0.8$($p = n$)')
#plt.legend(loc=3, fontsize=16,frameon=False, handlelength=3)
#plt.ylim(-0.05,0.05)
##plt.xlim(0.1,0.9)
#plt.xticks(fontsize=16)
#plt.yticks(fontsize=16)
#plt.xlabel("$Q^2$", fontsize=25)
#plt.ylabel("$F_2^D$((OnShell-OffShell)/OnShell)", fontsize=25)
#f15.savefig('./Results/pEQn_F2Ds_Off_vs_On_Add_fixd_x.pdf')


f1=plt.figure(1)
plt.xlim(0.2,1.0)
plt.ylim(0.0,0.4)
SinglePDFPlotWithBand(CJ15ioffXMul_pdf,2,5,'green','Multiplicative_HT($p = n$)')
SinglePDFPlotWithBand(CJ15ioffXAdd_pdf,2,5,'orange','Additive_HT($p = n$)')
plt.ylabel('$d/u$', fontsize=25)
plt.xlabel('$x$', fontsize=25)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.legend(loc=1, fontsize=16,frameon=False, handlelength=3)
f1.savefig('./Results/pEQn_du.pdf')


f2=plt.figure(2)
plt.xlim(0.0,0.8)
plt.ylim(0.8,1.6)
SinglePDFPlotWithBand(CJ15ioffXMul_pdf,2,3,'green','Multiplicative_HT($p = n$)')
SinglePDFPlotWithBand(CJ15ioffXAdd_pdf,2,3,'orange','Additive_HT($p = n$)')
plt.ylabel('$\\bar{d}/\\bar{u}$', fontsize=25)
plt.xlabel('$x$', fontsize=25)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.legend(loc=1, fontsize=16,frameon=False, handlelength=3)
f2.savefig('./Results/pEQn_dbub.pdf')



f3=plt.figure(3)
F2NP_Plot_FixedQ2(F2_N_mul_off,F2_P_mul_off,'10.000','green','Multiplicative_HT($p = n$)')
F2NP_Plot_FixedQ2(F2_N_add_off,F2_P_add_off,'10.000','orange','Additive_HT($p = n$)')
plt.legend(loc=3, fontsize=16,frameon=False, handlelength=3)
plt.ylim(0.2,0.6)
plt.xlim(0.2,1.0)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.xlabel("$x$", fontsize=25)
plt.ylabel("$F_2^n/F_2^p$", fontsize=25)
f3.savefig('./Results/pEQn_F2PNs.pdf')


f4=plt.figure(4)
plt.suptitle('Off-shell vs $x$', fontsize=16)
OffShellPlot(centralvaluespath_Mul,mainfitname_mul,'green','Multiplicative_HT($p = n$)')
OffShellPlot(centralvaluespath_Add,mainfitname_add,'orange','Additive_HT($p = n$)')
#plotpoints('./Data/AKP-OS-UpperBound.txt','black','AKP (2017)',0.0)
#plotpoints('./Data/AKP-OS-LowerBound.txt','black','',0.0)
plt.legend(loc=3, fontsize=16,frameon=False, handlelength=3)
plt.ylim(-0.7,0.3)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.xlabel("$x$", fontsize=25)
plt.ylabel("Off-shell", fontsize=25)
f4.savefig('./Results/pEQn_OSs.pdf')





#
#import matplotlib.backends.backend_pdf
#pdf = matplotlib.backends.backend_pdf.PdfPages("./Results/CJ15_dch1.pdf")
#for fig in xrange(1,10):
#    pdf.savefig( fig )
#pdf.close()
#    
#      
