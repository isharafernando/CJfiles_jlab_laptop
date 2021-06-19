#!/usr/bin/env python2.7 
#-------------------------------------------
# summary chi^2 tables for a series of fits
#-------------------------------------------
import sys, os, glob
import CJ

#----------------2-------------------
# shell.style aliases
from os import getcwd as pwd
from os import chdir as cd

def ls():
    
    for f in os.listdir(pwd()):
        print f.ljust(max)


def summary_chi2(files):
# ---- summary chi^2 (all data sets) -----

    for i in range(len(files)):
        files[i] = files[i].replace('.out','')
    # reorder the negative fmKP off-shell values 
    negs = [f for f in files if '_m' in f] 
    if len(negs)>0:
        start = files.index(negs[0])
        negs.sort(reverse=True)
    else:
        start = 0
    files[start:start+len(negs)]=negs

    test=CJ.out(files)

    headers = []
    for f in files:
        headers.append(f.replace('.out','').split('_')[-1])

    test.summarychi2(headers)


def random_chi2(cut):
# ---- summary chi^2 for random seed runs -----
#   cut = upper bound on total chi^2 to consider

  # gets all .out files in the current folder
    cwd = os.getcwd()
    raw = os.listdir(cwd)
    files = [f for f in raw if '.out' in f]
    files.sort()
    prunedfiles=[]
    headers = []
    for f in files:
        prunedfiles.append(f.replace('.out',''))
        headers.append(f.replace('.out','').split('_')[-1])

    test=CJ.out(prunedfiles)
    test.summarychi2_stats(cut,headers)


def pbp():
# ---- point-by-point chi^2 summary, selected data sets -----

    files_sel = ['off_AV18_ON','off_AV18_KP','off_AV18_m6','off_AV18_0','off_AV18_4']

    files = ['off_AV18_bns_ON','off_AV18_bns_KP','off_AV18_bns_m6','off_AV18_bns_0','off_AV18_bns_4']
    headers = ['ON','KP','-0.6%','0%','0.4%']
    test=CJ.out(files_sel)
    test.table_result('D0_Wasy','theory/data',headers)
    test.table_result('D0_Wasy','chi2',headers)    
    test.table_result('BNS_F2nd','theory/data',headers)
    test.table_result('BNS_F2nd','chi2',headers)
    test.table_result('e866pd06xf','theory/data',headers)
    test.table_result('e866pd06xf','chi2',headers)
    test.table_result('NmcRat','theory/data',headers)
    test.table_result('NmcRat','chi2',headers)
    test.table_result('SLAC_d','theory/data',headers)
    test.table_result('SLAC_d','chi2',headers)
    

def split_calcs():
# ---- splits results of calcs in blocks

    calc = calcfile('calc_DN_AV18.out')
    headers = ['on-shell','KP model offshell','fmKP offshell = -0.9%','fmKP offshell = -0.8%','fmKP offshell = -0.7%','fmKP offshell = -0.6%','fmKP offshell = -0.5%','fmKP offshell = -0.4%','fmKP offshell = -0.3%','fmKP offshell = -0.2%','fmKP offshell = -0.1%','fmKP offshell = 0%','fmKP offshell = 0.1%','fmKP offshell = 0.2%','fmKP offshell = 0.3%','fmKP offshell = 0.4%']
    calc.add_headers(headers)

#print
files=['./data/CJ15/CJ15.out','./data/CJ15pDn/CJ15pDn.out','./data/CJ15pDn_os0/CJ15pDn_noOS.out']
summary_chi2(files)

#-------------------------------------------------------
#  MAIN
#-------------------------------------------------------

#print

### Anything .out, alphabetic order

#files = sorted(glob.glob('*out'))[:]

### offshell stuff

#files = ['CJ15','fit_ioff9_KPe','fit_ioff9_CJ15','fit_ioff9_KPa','fit_ioff9_KPc','fit_ioff9_KPd','fit_ioff9_KPe','fit_ioff9_KPf','fit_ioff9_KPg','fit_ioff9_KPh','fit_ioff9_KPi','fit_ioff9_KPj']

#files = ['CJ15','fit_ioff1_CJ15-2','fit_ioff1_CJ15-3','fit_ioff1_KPe-3','fit_ioff1_KPe-3noW']

#files = ['CJ15','fit_ioff9_KPe','fit_ioff9_KPf','fit_ioff9_KPg','fit_CJ15-prime','fit_ioff9_KPe-prime','fit_ioff9_KPf-prime','fit_ioff9_KPg-prime']

#files = ['CJ15','fit_CJ15-prime','fit_ioff9_CJ15-noW','fit_ioff9_CJ15-noW-prime']

#files = ['CJ15','fit_ioff9_KPa','fit_ioff9_KPc','fit_ioff9_KPd','fit_ioff9_KPa-prime','fit_ioff9_KPc-prime','fit_ioff9_KPd-prime']

#files = ['CJ15','fit_ioff9_CJ15','fit_ioff9_KPa','fit_ioff9_KPc','fit_ioff9_KPd']

#summary_chi2(files)
#exit()
