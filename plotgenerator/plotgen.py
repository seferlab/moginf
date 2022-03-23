import networkx as nx
import numpy as np
import scipy as sp
import math
import random
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm

plotmarkers=["s","o","+","*",".","D","x","s","^","v","s","o","+","*",".","D","x","s","^","v","s","o","+","*",".","D","x","s","^","v"]
plotcolors=["r","b","g","c","k","m","y","r","b","g","c","k","m","y","r","b","g","c","k","m","y"]
plotfolder="newplots"

filepath="scoreout_partial_overtime-f1div10-noise1.0-sample1.txt"
xtitle="Trace count"
ytitle="Average F0.1"
plt.clf()
plt.xlabel(xtitle)
plt.ylabel(ytitle)
file=open(filepath,"r")
index=0
topx=-1.0
for line in file:
    line=line.rstrip()
    datainfo,xstr,ystr=line.split("\t")
    xparts=xstr.split(",")
    xvals=[]
    for elem in xparts:
        xvals.append(int(float(elem)*500))
    if max(xvals)>=topx:
       topx=max(xvals) 
    yparts=ystr.split(",")
    yvals=[]
    for elem in yparts:
        yvals.append(float(elem))
    print len(xvals)
    print len(yvals)
    print datainfo
    plt.plot(xvals,yvals,marker=plotmarkers[index],linestyle='None',color=plotcolors[index],label=datainfo)
    index+=1
file.close()
plt.xlim([0,topx+2])
plt.ylim([0,1.0]) 
plotfilepath="contact-overtime-f0.1-noise1.0-sample1.0.png"
plt.legend(loc=4)   
plt.savefig(plotfilepath)

filepath="scoreout_partial_distinfo_fixedtrace-f1div10-1.0-1-0.2.txt"
xtitle="Exponential distribution mean"
ytitle="Average F0.1"
plt.clf()
plt.xlabel(xtitle)
plt.ylabel(ytitle)
file=open(filepath,"r")
index=0
topx=-1.0
bottomx=10.0
for line in file:
    line=line.rstrip()
    datainfo,xstr,ystr=line.split("\t")
    xparts=xstr.split(",")
    xvals=[]
    for elem in xparts:
        xvals.append(float(elem))
    if max(xvals)>=topx:
       topx=max(xvals)
    if min(xvals)<=bottomx:
       bottomx=min(xvals)    
    yparts=ystr.split(",")
    yvals=[]
    for elem in yparts:
        yvals.append(float(elem))
    print len(xvals)
    print len(yvals)
    print datainfo
    plt.plot(xvals,yvals,marker=plotmarkers[index],linestyle='None',color=plotcolors[index],label=datainfo)
    index+=1
file.close()
plt.xlim([bottomx-0.05,topx+0.05])
plt.ylim([0,1.0]) 
plotfilepath="email-distinfo-f0.1-noise0.8-sample1.0.png"
plt.legend(loc=4)   
plt.savefig(plotfilepath)

filepath="scoreout_partial_sampling_fixedtrace-f1div10-1.0-0.2.txt"
xtitle="Sampling intervals"
ytitle="Average F0.1"
plt.clf()
plt.xlabel(xtitle)
plt.ylabel(ytitle)
file=open(filepath,"r")
index=0
topx=-1.0
for line in file:
    line=line.rstrip()
    datainfo,xstr,ystr=line.split("\t")
    xparts=xstr.split(",")
    xvals=[]
    for elem in xparts:
        xvals.append(float(elem))
    if max(xvals)>=topx:
       topx=max(xvals) 
    yparts=ystr.split(",")
    yvals=[]
    for elem in yparts:
        yvals.append(float(elem))
    print len(xvals)
    print len(yvals)
    print datainfo
    plt.plot(xvals,yvals,marker=plotmarkers[index],linestyle='None',color=plotcolors[index],label=datainfo)
    index+=1
file.close()
plt.xlim([0,topx+0.05])
plt.ylim([0,1.0])
xticks = [1,2,3,4,5,6,8,12]
plt.xticks(xticks)
plotfilepath="contact-sampling-f0.1-noise1.0.png"
plt.legend(loc=1)   
plt.savefig(plotfilepath)

filepath="scoreout_partial_noise_fixedtrace-f1div10-1-0.2.txt"
xtitle="Partial observability(Uncertainity) degree"
ytitle="Average F0.1"
plt.clf()
plt.xlabel(xtitle)
plt.ylabel(ytitle)
file=open(filepath,"r")
index=0
topx=-1.0
for line in file:
    line=line.rstrip()
    datainfo,xstr,ystr=line.split("\t")
    xparts=xstr.split(",")
    xvals=[]
    for elem in xparts:
        xvals.append(1.0-float(elem))
    if max(xvals)>=topx:
       topx=max(xvals) 
    yparts=ystr.split(",")
    yvals=[]
    for elem in yparts:
        yvals.append(float(elem))
    print len(xvals)
    print len(yvals)
    print datainfo
    plt.plot(xvals,yvals,marker=plotmarkers[index],linestyle='None',color=plotcolors[index],label=datainfo)
    index+=1
file.close()
plt.xlim([0-0.05,topx+0.05])
plt.ylim([0,1.0]) 
plotfilepath="email-partial-f0.1-sample1.0.png"
plt.legend(loc=1)   
plt.savefig(plotfilepath)

filepath="scoreout_parameterauccurve-roc-noise1.0-sample1-frac0.6.txt"
xtitle="FPR"
ytitle="Sensitivity"
plt.clf()
plt.xlabel(xtitle)
plt.ylabel(ytitle)
file=open(filepath,"r")
index=0
topx=-1.0
for line in file:
    line=line.rstrip()
    datainfo,xstr,ystr=line.split("\t")
    xparts=xstr.split(",")
    xvals=[]
    for elem in xparts:
        xvals.append(float(elem))
    if max(xvals)>=topx:
       topx=max(xvals) 
    yparts=ystr.split(",")
    yvals=[]
    for elem in yparts:
        yvals.append(float(elem))
    print len(xvals)
    print len(yvals)
    print datainfo
    plt.plot(xvals,yvals,marker=plotmarkers[index],linestyle='None',color=plotcolors[index],label=datainfo)
    index+=1
file.close()
plt.xlim([0,1.0])
plt.ylim([0,1.0]) 
plotfilepath="parameterauccurve-roc-noise1.0-sample1-frac0.6.png"
plt.legend(loc=4)   
plt.savefig(plotfilepath)

filepath="scoreout_partial_noise_fixedtrace-f1-1-0.4.txt"
ytitle="Average F1"
xtitle="Partial observability(Uncertainity) degree"
plt.clf()
plt.xlabel(xtitle)
plt.ylabel(ytitle)
file=open(filepath,"r")
index=0
topx=-1.0
for line in file:
    line=line.rstrip()
    datainfo,xstr,ystr=line.split("\t")
    xparts=xstr.split(",")
    xvals=[]
    for elem in xparts:
        xvals.append(1.0-float(elem))
    if max(xvals)>=topx:
       topx=max(xvals) 
    yparts=ystr.split(",")
    yvals=[]
    for elem in yparts:
        yvals.append(float(elem))
    print len(xvals)
    print len(yvals)
    print datainfo
    plt.plot(xvals,yvals,marker=plotmarkers[index],linestyle='None',color=plotcolors[index],label=datainfo)
    index+=1
file.close()
plt.xlim([0-0.05,topx+0.05])
plt.ylim([0,1.0]) 
plotfilepath="ff-noise-f1-sample1.0-frac0.4.png"
plt.legend(loc=1)   
plt.savefig(plotfilepath)

filepath="scoreout_realdynamic_overtime-f1-noise1.0-sample1.txt"
xtitle="Trace count"
ytitle="Average F1"
plt.clf()
plt.xlabel(xtitle)
plt.ylabel(ytitle)
file=open(filepath,"r")
index=0
topx=-1.0
for line in file:
    line=line.rstrip()
    datainfo,xstr,ystr=line.split("\t")
    xparts=xstr.split(",")
    xvals=[]
    for elem in xparts:
        xvals.append(int(float(elem)*500))
    if max(xvals)>=topx:
       topx=max(xvals) 
    yparts=ystr.split(",")
    yvals=[]
    for elem in yparts:
        yvals.append(float(elem))
    print len(xvals)
    print len(yvals)
    print datainfo
    plt.plot(xvals,yvals,marker=plotmarkers[index],linestyle='None',color=plotcolors[index],label=datainfo)
    index+=1
file.close()
plt.xlim([0,topx+2])
plt.ylim([0,1.0]) 
plotfilepath="dynamiccontact-overtime-f1-noise1.0-sample1.0.png"
plt.legend(loc=4)   
plt.savefig(plotfilepath)

filepath="scoreout_syndynamic_overtime-f1-noise1.0-sample1.txt"
xtitle="Trace count"
ytitle="Average F1"
plt.clf()
plt.xlabel(xtitle)
plt.ylabel(ytitle)
file=open(filepath,"r")
index=0
topx=-1.0
for line in file:
    line=line.rstrip()
    datainfo,xstr,ystr=line.split("\t")
    xparts=xstr.split(",")
    xvals=[]
    for elem in xparts:
        xvals.append(int(float(elem)*500))
    if max(xvals)>=topx:
       topx=max(xvals) 
    yparts=ystr.split(",")
    yvals=[]
    for elem in yparts:
        yvals.append(float(elem))
    print len(xvals)
    print len(yvals)
    print datainfo
    plt.plot(xvals,yvals,marker=plotmarkers[index],linestyle='None',color=plotcolors[index],label=datainfo)
    index+=1
file.close()
plt.xlim([0,topx+2])
plt.ylim([0,1.0]) 
plotfilepath="dynamicsyn-overtime-f1-noise1.0-sample1.0.png"
plt.legend(loc=4)   
plt.savefig(plotfilepath)

