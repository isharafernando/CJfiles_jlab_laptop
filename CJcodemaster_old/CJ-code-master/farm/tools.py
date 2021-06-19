#!/usr/bin/env python
import sys,os
import numpy as np
import time
import fnmatch
import cPickle 
from operator import mul

def logo():
  print """
  =======================================
    ____    _ ____   ____    _    _   _ 
   / ___|  | / ___| / ___|  / \  | \ | |
  | |   _  | \___ \| |     / _ \ |  \| |
  | |__| |_| |___) | |___ / ___ \| |\  |
   \____\___/|____/ \____/_/   \_\_| \_|
  
  authors: 
    N. Sato
  =======================================
  """

def timer(f,*args):
  t0=time.time()
  f(*args)
  t1=time.time()
  print 't = %5.0e '%(t1-t0),f.func_name
  return t1-t0

def checkdir(path):
  if not os.path.exists(path): 
    os.makedirs(path)

def tex(x):
  return r'$\mathrm{'+x+'}$'

def save(data,name):  
  f=open(name,"w")
  cPickle.dump(data, f)
  f.close()

def load(name):  
  f=open(name,"r")
  data=cPickle.load(f)
  f.close()
  return data

def isnumeric(value):
  try:
    int(value)
    return True
  except:
    return False

  return r'$\mathrm{'+x+'}$'

class BAR(object):

  def __init__(self,msg,size):
    self.msg=msg
    self.size=size
    self.cnt=0

  def next(self):
    sys.stdout.write('\r')
    percentage=int(self.cnt/float(self.size)*100)
    sys.stdout.write('%s [%d%%]' % (self.msg,percentage))
    sys.stdout.flush()
    self.cnt+=1

  def finish(self):
    self.next()
    print 
