#!/usr/bin/env python
import sys,os
import numpy as np
from tools import checkdir
import subprocess
import time
##########################################

class FARM:

  def __init__(self,conf):
    self.conf=conf
    checkdir('out')
    checkdir('err')
  
  def gen_xml(self,tags):
    L=[]
    L.append('<Request>')
    L.append('<Project name="cteqX"/>')
    #L.append('<Track name="reconstruction"/>')
    L.append('<Track name="debug"/>')
    L.append('<OS name="centos7"/>')
    L.append('<Name name="[[jobname]]"/>')
    L.append('<TimeLimit time ="[[time]]" unit="hours"/>')
    L.append('<Job>')
    L.append('<Command>')
    L.append('<![CDATA[')
    L.append('setenv cteqx [[path2fitpack]]')
    L.append('[[exe]]  [[input]]')
    L.append(']]>')
    L.append('</Command>')
    L.append('<Memory space="1050" unit="MB"/>')
    L.append('<Input src="[[path2exe]]" dest="[[exe]]"/>')
    L.append('<Input src="[[path2input]]/[[input]]" dest="[[input]]"/>')
    L.append('<Output src="[[out]]" dest="[[path2out]]/[[out]]"/>')
    L.append('<Output src="[[pdf]]" dest="[[path2pdf]]/[[pdf]]"/>')
    L.append('<Output src="[[par]]" dest="[[path2par]]/[[par]]"/>')
    L.append('<Stdout dest="[[root]]/out/[[jobname]].out"/>')
    L.append('<Stderr dest="[[root]]/err/[[jobname]].err"/>')
    L.append('</Job>')
    L.append('</Request>')
  
    for k in tags.keys():
      for i in range(len(L)):
        L[i]=L[i].replace(k,tags[k])
    L=[l+'\n' for l in L]
    F=open('job.xml','w')
    F.writelines(L)
    F.close()
  
  def send_job(self,jname):    
    conf=self.conf

    tags={}
    tags['[[time]]']='6'
    tags['[[root]]']=conf['cwd']
    tags['[[path2fitpack]]']=conf['path2fitpack']
    tags['[[exe]]']=conf['exe']
    tags['[[input]]']='%s.dat'%jname
    tags['[[out]]']='%s.out'%jname
    tags['[[pdf]]']='%s.pdf'%jname
    tags['[[par]]']='%s.par'%jname
    tags['[[path2exe]]']='%s/%s'%(conf['cwd'],conf['exe'])
    tags['[[path2input]]']='%s'%conf['cwd']
    tags['[[path2out]]']='%s'%conf['cwd']
    tags['[[path2pdf]]']='%s'%conf['cwd']
    tags['[[path2par]]']='%s'%conf['cwd']
    tags['[[jobname]]']=jname
    self.gen_xml(tags)
 
    cmd=['jsub','-xml','job.xml']
    delay=5
    while 1:
      out=subprocess.check_output(cmd)
      if ('error' in out)==False:
        print 'fit %s has been submitted'%jname
        break
      else:
        print 'error occured, waiting %d s for resubmission'%delay
        time.sleep(delay)


if __name__=='__main__':

  cwd=os.getcwd()
  conf={'cwd':cwd}
  conf['exe']='pro'
  conf['path2fitpack']='%s/../CJ-code/fitpack'%cwd

  farm=FARM(conf)
  farm.send_job('CJ15')  # <-- this is the inputfile *.dat    





