#!/usr/bin/env python2.7 
import sys
import random
import numpy as np

#----------------------------------------------------
class dat(object):
  #---- reads in and manipualtes a list of .out files

  def __init__(self,files):

    self.datfiles=datfiles

  def get(self,fname):

    try:
      F=open(fname+'.out','r')
    except IOError:
      print "WARNING: get_summarychi2 cannot open'"+fname+".out'"
    else:
      L=F.readlines()
      F.close()
      L=[l.strip() for l in L]
    return L

  def rmline(self,readfile,substring):
  # removes (in place) a line from the current read in file 
  # if it contains the given substring.


    
    for file in N:
      new = []
      for l in L:
        if substring not in readfile:
          new.append(l)
  
#  def print(self)

    return new

#----------------------------------------------------
class out(object):
  #---- reads in and manipualtes a list of .out files

  def __init__(self,files):

    self.N=files

  def mean(self,list):

    mean = 0
    for x in list: mean += x
    mean /= len(list) 
    return mean

  def stdv(self,list):
    
    stdv = 0
    xbar = self.mean(list) 
    for x in list: 
      stdv += (x - xbar)**2
    stdv = (stdv/(len(list)-1))**0.5

    return stdv

  def summarychi2(self,headers):

    table = self.maketable_summarychi2(headers)
    self.printtable_summarychi2(table)

  def summarychi2_stats(self,cut,headers):

    table = self.maketable_summarychi2(headers)
    newl = list(table[0])
    newl.insert(16,"    mean    sdev")
    table[0] = ''.join(newl)

    # ... finds Nan and obviously false minima (defined by TOTAL chi2>cut)
    index = -1
    remove = []
    totals = table[-1].split()[1:]
    for totchi2 in totals:  
      index += 1
      if float(totchi2)>cut or totchi2=='NaN': remove.append(index)
    
    newtable = []
    for l in table:
      if 'data' not in l and l!='':
        if 'TOTAL+norm' in l:
          start=1
        else:
          start=2
        chi2 = [float(i) for i in l.split()[start:]]
        for i in sorted(remove,reverse=True): 
          chi2.pop(i) # removes obviously false minima
        add = [self.mean(chi2),self.stdv(chi2)]

        newl = list(l)
        newl.insert(16,"%8.1f"%add[0]+"%8.1f"%add[1])
        l = ''.join(newl)
      newtable.append(l)
      
    self.printtable_summarychi2(newtable)
    print "REMOVED: ", [headers[i] for i in remove] 

  def printtable_summarychi2(self,table):

    L = table
    for j in range(len(L)): 
      if j==len(L)-4:
	print '-'*len(L[0])
      print L[j]

  def maketable_summarychi2_old(self,headers):

    # local retrieve
    fits=self.N

    CHI2={}
    for f in fits:
      H,names,npts,chi2 = self.get_summarychi2(f)
      if H != '---': 
        NMS = names
        NPTS = npts
        CHI2[f]=chi2
      
    L=[]
    print NMS
    for i in range(len(NMS)):
      row=[]
      row.append(NMS[i])
      row.append(NPTS[i])
      for f in fits:
	row.append(CHI2[f][i])
      L.append(row)

    H=['data set','npts']+[name.rjust(8) for name in headers]
    L.insert(0,H)

    table = []
    for l in L: 
      row =l[0].ljust(12)+l[1].rjust(4)
      for i in range(2,len(l)): row+=l[i].rjust(8)
      table.append(row)

    table.insert(1,'')

    return table

  def maketable_summarychi2(self,headers):

    # local retrieve
    fits=self.N

    CHI2={}
    NPTS={}
    NMS = []
    for f in fits:
      CHI2[f]={}
      H,names,npts,chi2 = self.get_summarychi2(f)
      #print '* names: ',names
      if H != '---': 
        if set(names).issuperset(set(NMS)):
          NMS = names
        #else:
        #  NMS = list(set(NMS+names))
        idata=-1
        for dataset in names:
          idata+=1
          NPTS[dataset] = npts[idata]
          CHI2[f][dataset] = chi2[idata]
          
    L=[]
    print '* NMS:   ',NMS 
    fitpts={}
    for f in fits: fitpts[f]= -int(float(NPTS['TOTAL']))
    for dataset in NMS:
      row=[]
      row.append(dataset)
      try:
        row.append(NPTS[dataset])
      except:
        row.append('')     
      for f in fits:
        try:
	  row.append(CHI2[f][dataset])
          try:
            fitpts[f] += int(float(NPTS[dataset]))
          except:
            print 'EXCEPT -- file=',f,' set=',dataset
        except:
          row.append('')
      L.append(row)

    row=['no. pts','']
    row2=['(TOT+nrm)/pts','']
    for f in fits:
      row.append(str(fitpts[f]))
      row2.append('%6.2f '%(float(CHI2[f]['TOTAL+norm'])/fitpts[f]))
    L.append(row)
    L.append(row2)

    H=['data set','npts']+[name.rjust(8) for name in headers]
    L.insert(0,H)

    print np.asmatrix(L)

    table = []
    for l in L: 
      row =l[0].ljust(12)+l[1].rjust(4)
      for i in range(2,len(l)): row+=l[i].rjust(8)
      table.append(row)

    table.insert(1,'')

    return table


  def get_summarychi2(self,fname):
  #Extracts the chi^2 summary table for a single out file
    
    try:
      F=open(fname+'.out','r')
    except IOError:
      print "WARNING: get_summarychi2 cannot open'"+fname+".out'"
      H='---'
      return H,0,0,0
    else:
      L=F.readlines()
      F.close()
      L=[l.strip() for l in L]

      ini_norm = 0
      ini = 0
      fin = 0
      for i in range(len(L)): 
        if 'Summary normalization chi2 contributions' in L[i]: ini_norm=i+1 
        if 'Summary chi2 per data set' in L[i]: ini=i+1
        if 'Momentum Sum Rule Check' in L[i]: fin=i

      #print '*',fname
      #print '* ini_norm, ini, fin =', ini_norm, ini, fin

      if ini>0 and ini_norm>0:  # found summary chi^2 tables
        
        # get the normalization chi^2

        Ln = L[ini_norm:ini-1]  # normalization chi^2 lines
        Ln=[l for l in Ln if l!='']  #  removes blank lines

        Hn = Ln[0]
        Ln1 = Ln[1:]
        norm_names=[l.split()[0] for l in Ln1]
        norm_npts = ['' for l in Ln1]
        norm_chi2=[l.split()[3] for l in Ln1]

        norm_names.extend(['',''])
        norm_npts.extend(['',''])
        norm_chi2.extend(['',''])

        # now the data chi^2

        L=L[ini:fin]   #  data chi^2 lines
        L=[l for l in L if l!='']  # strips empty lines

        H=L[0]         # headers
        L1=L[1:-3]     # data lines 
        L2=L[-2:]      # total chi^2 lines

        names=[l.split()[0] for l in L1]   # names of data sets 
        npts=[l.split()[4] for l in L1]    # no. of points of each data set
        chi2=[l.split()[1] for l in L1]    # total chi^2 for that data set 

        L2=[l.split() for l in L2]
        names.extend([l[0] for l in L2])
        npts.extend([L2[0][7],''])
        chi2.extend([l[1] for l in L2])
      else:
        print "WARNING (get_summarychi2) '"+fname+".out' is incomplete"
        H = '---'
        names = []
        npts = []
        chi2 = []

      F.close()

      #norm=True
      norm=False
      if norm:
        names=norm_names
        npts = norm_npts
        chi2=norm_chi2

      return H,names,npts,chi2


  def table_result(self,dataset,what,headers):
    # what: 'chi2' 'data/theory' 

    N=self.N # file names local retriev

    rdflag = 1
    for fil in N:
      D, legend = self.get_results(fil,dataset)     
      if fil == N[0]: table = D['kin'] 
      if what=='chi2':
        width = 10
        table.extend([D[what]])
      elif what=='data/theory' :
        width = 18
        d = [float(D['data'][i]) for i in range(len(D['data']))]
        t = [float(D['theory'][i]) for i in range(len(D['theory']))]
        r = [float(D['error'][i]) for i in range(len(D['theory']))]
        dt = [["%4.3f"%(d[i]/t[i])+' +- '+ "%4.3f"%(r[i]/d[i]) for i in range(len(t))]]
        table.extend(dt)
      elif what=='theory/data':
        width = 10
        t = [float(D['theory'][i]) for i in range(len(D['theory']))]
        d = [float(D['data'][i]) for i in range(len(D['data']))]
        r = [float(D['error'][i]) for i in range(len(D['theory']))]
        td = [["%4.3f"%(t[i]/d[i]) for i in range(len(d))]]
        if rdflag: 
          dd = [["%4.3f"%1. for i in range(len(d))]]
          rd = [["%4.3f"%(r[i]/d[i]) for i in range(len(d))]]
          table.extend(dd)
          table.extend(rd)
          rdflag = 0
        table.extend(td)
      else:
        print 'ERROR (table_results): option not available: ',what
        exit()

    rows = range(len(table[0]))
    chi2table = []
    for i in rows: 
      chi2table.extend([[l[i] for l in table]])

    print
    print dataset+'  -  '+what
    print
    row=''
    for i in legend:
      row += i.center(10)
    if what=='theory/data': 
      row += 'data'.center(width)+' exp.err'.center(width)
    for i in headers: 
      row += i.center(width)
    print row
    for line in chi2table:
      chi2row = ''
      for num in line:
        if '+-' in num:
          size = width
        else:
          size = 10
        chi2row += num.rjust(size)
      print chi2row
        

  def get_results(self,fname,dataset):

    F=open(fname+'.out','r')
    outf = F.readlines()  # .out file read in 
    F.close()
    outf = [l.strip() for l in outf]  # strip trailing stuff

    # find data vs theory block for specified dataset
    block=[line for line in outf if line.find(dataset) > 10]
    # find its legend
    index = -2
    for line in outf:
      index += 1
      if line in block: break 
    legend = outf[index]
    ## find its statistics summary
    #summary = []
    #found = 0
    #count = 0
    #for line in outf:
    #  noset = 0
    #  if line.find(dataset) > 10:
    #    found=1
    #  else: 
    #    if line != '':
    #      noset=1
    #  print found, noset, line
    #  if found and noset: 
    #    summary.append(line)
    #    count += 1
    # if count==3: break

    # how many kinematic columns?
    index = -1
    for name in legend.split():
      index += 1
      if name == 'THEORY': nkin = index
      if name == 'DATA'  : ndat = index
      if name == 'ERROR' : nerr = index
      if name in ['chi2','chi^2']: nchi2 = index
    # extracts those columns into rows, store them in a dictionary
    results = {}  # teh results dictionary
    cols = range(nkin)
    kin = []
    for i in cols: 
      kin.extend([[row.split()[i] for row in block]])
    results['kin'] = kin
    results['theory'] = [row.split()[nkin] for row in block]
    results['data']   = [row.split()[ndat] for row in block]
    results['error']  = [row.split()[nerr] for row in block]
    results['chi2']   = [row.split()[nchi2] for row in block]

    # returns the dicitonary and the legend for the kinematic variables
    return results, legend.split()[:nkin]



#-------------------------------------------------------------
class dat(object):
  # reads in and manipulates the input .dat file
 
  def __init__(self,datfile):

    F=open(datfile+'.dat','r')
    datf = F.readlines()  # .out file read in 
    F.close()
    self.datf = [l.rstrip() for l in datf]  # strip trailing \n


  def output(self,fname=False):
    # saves to file "fname" (prints to screen if fname omitted)

    datf = self.datf
    if fname:
      outf = open(fname+'.dat','wb')
      for l in datf: outf.write(l+'\n')
    else:
      for l in datf: print l

    
  def new_nuke(self,nuke=False):
    # In-place substitution of the inuke code 
    # 
    # INPUT: [nuke] = (s) the new nuclear label
    #                 [no action if absent)
    #
    # RETURNS: the inuke label before substition.

    datf=self.datf # local retrieve - but only of the refernece!!
                   # self.datf will be modified as well !!

    flag = False
    for i in range(len(datf)):
      if 'END' in datf[i]: 
        flag=True
        continue
      if 'dummy' in datf[i]: 
        flag=False
      if flag:
        old = datf[i].split(',')[4]
        if float(old) > 0 :
          if nuke:
            datf[i] = datf[i].replace(old.strip(),nuke)
    
    return old
      

  def randomize_seed(self,sigma):
    # Randomizes the initial parameters with a flat
    #  probability distribution within +-sigma * uncertainty 
    # (with uncertainty) and central values taken from the dat file
    #
    # INPUT:  sigma = (f) no. of std deviations
    # 
    # RETURNS: line-by-line list with randomized .dat file

    datf=self.datf # local retrieve
    
    rand = []
    mix = 1
    hdr = 1
    for l in datf:
      if hdr: 
        rand.append(l)
        hdr = 0
        continue
      if 'END' in l: 
        mix = 0 
      if mix: # randomization
        s = l.split()
        for i in range(1,3): s[i] = float(s[i])
        if s[2]!=0:  # randomizes the seed
          delta = sigma * s[2]
          s[1] = s[1] + random.uniform(-delta,delta)
        newl = ' '+s[0].ljust(10) + "%13.4e"%s[1] + "%13.4e"%s[2]
        rand.append(newl)
      else: # copies all lines after the parameters
        rand.append(l)
    
    return rand

#-------------------------------------------------------------
class calcfile(object):
  # reads in and manipulates the output of 
  # a calculation made with fitpack
 
  def __init__(self,calcfile):

    F=open(calcfile,'r')
    outf = F.readlines()  # .out file read in 
    F.close()
    self.outf = [l.rstrip() for l in outf]  # strip trailing \n

  def add_headers(self,headers):

    outf=self.outf # local retrieve

    start = 0
    ihdr = 0
    for l in outf:
      if 'THEORY' in l: 
        start=1
        legend = l[:]
        continue
      if 'CHI SQUARE=' in l: break
      if start:
        if '.1000E-01' in l:
          print 
          print headers[ihdr]
          print 
          print legend
          ihdr += 1
        print l

