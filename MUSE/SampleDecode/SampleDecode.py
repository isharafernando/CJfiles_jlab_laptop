# -*- coding: utf-8 -*-
"""
Created on Fri Aug 30 11:27:08 2019

@author: isharafernando
"""

import numpy as np
#import sys
import matplotlib.pyplot as plt

def parline(filename,word):
    temparray2=[]
    with open(filename) as myFile:   
        for num, line in enumerate(myFile, 0):
            if word in line:
                temparray2.append(num)          
    return temparray2    


def hexconversion(filename):
    templist=[]
    tempconvrt=[]
    tempwords=[]
    count=0
    BankNameline=parline(filename,'Bank:')[0]
    for line in open(filename):
        templist += [line.split()]
    Endline=parline(filename,'1577')[0]
    for i in range(BankNameline+1,Endline+1):
        linelength=len(templist[i])
        #tempelementsinline=[]
        for j in range(1,linelength):
            count=count+1
            tempword=templist[i][j]
            #stringtempword=str(tempword)
            selectstrings=tempword[7:10]
            inttempword=selectstrings
            #tempelementsinline.append(int(inttempword,16))
            tempwords.append([count,int(inttempword,16)])
            #yyy=str(xxxx)[1:4] 
            #yyyy=int(yyy)
            #print yyy
            #tempelementsinline.append(bin(int(templist[i][j],16)))
            #tempelementsinline.append(bin(int(templist[i][j],16)))
        #tempconvrt.append([templist[i][0],tempelementsinline])
        #tempconvrt.append([tempelementsinline])    
    return tempwords   
        
#BankNameline, templist[BankNameline+1][1],templist[-1] 
#int('0x0007019e',16)
        
hexconresult=hexconversion('SampleEvent.txt')        
#print (len(hexconresult)-4-3)/131

dataarrays=[]
pointer=4
for i in range(0,(len(hexconresult)-4-3)/131 +1):
    dataarrayssub=[]
    for j in range(pointer+i*131,pointer+(i+1)*131):
        dataarrayssub.append(hexconresult[j][1])
    dataarrays.append(dataarrayssub)    
        #dataarrays.append(dataarrayssub)    
    #print(hexconresult[4+i][1])




ii=0
print len(dataarrays[ii])
print len(dataarrays[ii][1:129])
print dataarrays[ii][1:129]
testarray=dataarrays[ii][1:129]
#x = [i for i in range(128)]

iii=0
if ii==0 or ii==4 or ii==8:
    iii=1
elif ii==1 or ii==5 or ii==9:
    iii=2
elif ii==2 or ii==6 or ii==10:
    iii=3
elif ii==3 or ii==7 or ii==11:
    iii=4
else:
    iii=0    
    

sample_num=1
if 0<=ii<4:
    sample_num=1
elif 4<=ii<8:
    sample_num=2
else:
    sample_num=3
    
print sample_num
    
#plt.plot(x,testarray,'.')
x = np.arange(128)
#plt.hist(x,height=np.array(testarray))
#plt.hist(testarray,bins=x,normed=True)
#plt.hist(testarray,bins=x)
f1=plt.figure(1)
plt.plot(x,testarray,'.')
plt.xlim(-1,128)
plt.xlabel('APV channel number')
plt.ylabel('ADC value (last 3 digits of hex word)')
plt.title('Data_Block#='+str(ii)+',APV#='+str(iii)+', Sample#='+str(sample_num))
plt.text(10, 850, 'Preliminary', fontsize=50, color='gray', alpha=0.5)
         

f2=plt.figure(2)
plt.hist(testarray)
plt.xlabel('ADC value (last 3 digits of hex word)')
plt.ylabel('Count')
plt.title('Data_Block#='+str(ii)+',APV#='+str(iii)+', Sample#='+str(sample_num))
plt.text(600, 10, 'Preliminary', fontsize=50, color='gray', alpha=0.5)
        
    
import matplotlib.backends.backend_pdf
pdf = matplotlib.backends.backend_pdf.PdfPages('./Results/DataBlock_'+str(ii)+'.pdf')
for fig in xrange(1,3):
    pdf.savefig( fig )
pdf.close()



#fname=open('sample_decode_last3dgts_decimals.txt','w')
#sys.stdout=fname
#print hextodecimal('SampleEvent.txt')
#fname.close()
#sys.stdout=sys._stdout_
#with open('sample_decode_last3dgts_decimals.txt') as file:
#    data=file.read()
#    file.close()

        