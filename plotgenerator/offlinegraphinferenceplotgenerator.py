import networkx as nx
import numpy as np
import scipy as sp
import random
import os
import math
import sys
import myutilities as myutil
import operator
import cPickle
import pickle
import string
import itertools
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm

def estimatenewrocscore(xlist,ylist,scoretype):
    crocdir="croc"
    codefilename="my{0}.py".format(scoretype)
    codepath="{0}/{1}".format(crocdir,codefilename)   
    tempfilename=''.join(random.choice(string.ascii_uppercase) for x in range(30))
    file=open(tempfilename,"w")
    for index in range(0,len(xlist)):
        file.write("{0}\t{1}\n".format(xlist[index],ylist[index]))
    file.close()
    outfilename="out_{0}".format(tempfilename)
    code="python {0} < {1} > {2}".format(codepath,tempfilename,outfilename)
    os.system(code)
    myscore=0.0
    file=open(outfilename,"r")
    for line in file:
        line=line.rstrip()
        myscore=float(line)
    file.close()
    code="rm -rf {0}".format(tempfilename)
    os.system(code)
    code="rm -rf {0}".format(outfilename)
    os.system(code)
    return myscore

def estimateaucscore(xlist,ylist,scoretype):
    if scoretype=="roc":
       return estimaterocscore(xlist,ylist)
    elif scoretype=="pr":
       return estimateprscore(xlist,ylist)
       #prscore=estimateprscore(xlist,ylist)
       #print estimatearea(xlist,ylist)
       #print prscore
       #print xlist
       #print ylist
       #exit(1)
    elif scoretype=="powerroc":
       return estimatenewrocscore(xlist,ylist,"power")
    elif scoretype=="bedroc": #bedroc is exponential transformation!
       return estimatenewrocscore(xlist,ylist,"exponential")
    elif scoretype=="logroc":
       return estimatenewrocscore(xlist,ylist,"logarithm")
    else:
       print "this type of {0} auc score is unknonw!!".format(scoretype)
       exit(1)
    
def estimaterocscore(xlist,ylist):
    #myscore=javaauc(tprs,fprs,poscount,negcount):
    initial=[0.0,0.0]
    last=[1.0,1.0]
    myscore=auc(ylist,xlist,initial,last)
    return myscore

def estimateprscore(xlist,ylist):
    initial=[1.0,1.0]
    last=[1.0,0.0]
    myscore=auc(ylist,xlist,initial,last)
    return myscore

#auc java code, due to poscount, negcount not using right now!!
def javaauc(tprs,fprs,poscount,negcount):
    chars=string.ascii_uppercase + string.digits
    size=12
    filename=''.join(random.choice(chars) for x in range(size))
    file=open(filename,"w")
    for index in range(0,len(tprs)):
        file.write("{0}\t{1}\n".format(fprs[index],tprs[index]))
    file.close()

    outfilename=''.join(random.choice(chars) for x in range(size))
    code="java -jar auc.jar {0} {1} {2} {3} > {4}".format(filename,"ROC",poscount,negcount,outfilename)
    os.system(code)
    
    file=open(outfilename,"r")
    for line in file:
        line=line.rstrip()
        if line.startswith("Area Under the Curve for ROC"):
           rocscore=float(line.split("is")[1])
    file.close()    
 
    code="rm -rf {0}".format(filename)
    os.system(code)
    code="rm -rf {0}".format(outfilename)
    os.system(code)
    code="rm -rf {0}.pr".format(filename)
    os.system(code)
    code="rm -rf {0}.spr".format(filename)
    os.system(code)
    code="rm -rf {0}.roc".format(filename)
    os.system(code)
    return rocscore

# there is also mirrored pr curve area of which is 1.0-auc of pr curve. so nothing extra for it has been implemented.
def auc(ylist, xlist,initial,last):
    xlist=list(xlist)
    ylist=list(ylist)
    if last[0]!=xlist[-1] or last[1]!=ylist[-1]:
       ylist.append(last[1])
       xlist.append(last[0])  
    if initial[0]!=xlist[0] or initial[1]!=ylist[0]:     
       ylist.insert(0,initial[1])
       xlist.insert(0,initial[0])
    #important change in the code!! not using numpy trapz anymore 
    #return float(np.trapz(ylist, x=xlist)) 
    return estimatearea(xlist,ylist,addinitial=False)

#first x coordinates in increasing order
#always add point 0,0!! so can not directly use for roc and pr score estimation
#addinitial is very imoportant for pr and roc scores
def estimatearea(xlist,ylist,addinitial=True):
    paramareascore=0.0
    if addinitial:
       tempxlist=[0.0]
       tempylist=[0.0]
       tempxlist.extend(xlist)
       tempylist.extend(ylist)
    else:
       tempxlist=list(xlist)
       tempylist=list(ylist) 
    sortedxlist=sorted(tempxlist)
    for index in range(1,len(sortedxlist)):
        prex=sortedxlist[index-1]
        curx=sortedxlist[index]
        prexindex=tempxlist.index(prex)
        curxindex=tempxlist.index(curx)
        prey=tempylist[prexindex]
        cury=tempylist[curxindex]
        partarea=float(curx-prex)*(cury+prey)/2.0
        paramareascore+=partarea
    assert paramareascore>=0
    return paramareascore

#weighted early retrieval score
def estimateearlyscore(fraclist,ylist):
    assert sorted(fraclist)==fraclist
    divsum=0.0
    if weight=="1-square":
       myscore=0.0
       for index in range(0,fraclist):
           myscore+=(((1.0-fraclist[index])**2)*ylist[index])
           divsum+=((1.0-fraclist[index])**2)
    elif weight=="1-normal":
       myscore=0.0
       for index in range(0,fraclist):
           myscore+=((1.0-fraclist[index])*ylist[index])
           divsum+=(1.0-fraclist[index])
    else:
       print "weight {0} is unknonw!!".format(weight)
       exit(1)
    return float(myscore)/divsum    
    
def plotlabelcreate(inferalgo,sentparameterstr,paramtype="-1"):
    if paramtype=="-1":
       if inferalgo in ["multitree","netrate","netinf","connie"]:
          datalabel="{0}".format(inferalgo)
       elif inferalgo in ["absel1","lsel1"]:
          datalabel="cefer_{0}".format(inferalgo)
       elif inferalgo in ["lsel1fused","absel1fused"]: #our algos dynamic
          datalabel="cefer_{0}_fused".format(inferalgo)  
       else:
          print "algo {0} is unknown for label creation!!".format(inferalgo)
          exit(1)
    elif paramtype=="algoparam":
       parameterstr,extend=sentparameterstr.split("+")
       parts=parameterstr.replace("result_","").split("_")
       if inferalgo in ["netinf","multitree"]:
          if parts[-1].find("categorical")!=-1:
             foundone="categorical?" 
          elif parts[-1].find("error")!=-1: 
             foundone="error?" 
          else:
             print "not known rounding method!!"
             exit(1) 
          algoparam=float(parts[-1].replace(foundone,""))
       elif inferalgo in ["netrate"]:
          algoparam="1" #there won't be param for netrate
       elif inferalgo in ["lsel1","absel1"]: #our algos
          algoparam=float(parts[-1].split("?")[0])
       elif inferalgo in ["lsel1fused","absel1fused"]: #our algos dynamic
          algoparam=float(parts[-1].split("?")[0]) 
       else:
          print "this algo is unknown {0}".format(sentalgo)
          exit(1)
       datalabel="-1"
       if inferalgo in ["multitree","netinf"]:
          datalabel="{0} k={1}".format(inferalgo,algoparam)
       elif inferalgo in ["netrate"]:
          datalabel="{0}".format(inferalgo)
       elif inferalgo in ["connie"]:
          datalabel="{0} lambda={1}".format(inferalgo,algoparam)   
       elif inferalgo in ["absel1","lsel1"]:
          datalabel="cefer_{0} lambda={1}".format(inferalgo,algoparam)
       elif inferalgo in ["lsel1fused","absel1fused"]: #our algos dynamic
          datalabel="cefer_{0} lambda={1}".format(inferalgo,algoparam) 
       else:
          print "algo {0} is unknown for label creation!!".format(inferalgo)
          exit(1)
    else:
       print "ERROR: paramtype {0} is unknonw!!".format(paramtype)
       exit(1)
    return datalabel
                

def paramchangeplotgenerate(sentscores,plottype,plottitle,plotfilepath):
    ykey=plottype
    if plottype in ["sen","fpr","recall","precision","acc","f1","f2","mcc","spec","f1overn","f1div100","f1div10","f1div2","f10","f20"]:
       maxormin="max" 
    else:
       print "this plottype {0} is unknown for max or min !!".format(plottype)
       exit(1) 
    myxlabel="Lambda"
    myxlabel2="Edge number"
    myylabel=plottype.capitalize()
    
    #preprocessing part
    plotxlist={}
    plotylist={}
    labelinfo={}
    for sentalgo in sentscores.keys():
        bestavgscore=-1.0
        bestylist=[]
        bestxlist=[]
        datalabel="-1"
        for sentparameterstr in sentscores[sentalgo].keys():
            #parameterstr,extend=sentparameterstr.split("+")
            tempylist=[] #y
            tempxlist=[] #x
            
            if sentalgo in ["lsel1","absel1"]:
               threslist=sorted(sentscores[sentalgo][sentparameterstr].keys(),reverse=True)
            elif sentalgo in ["multitree","netinf"]:
               threslist=sorted(sentscores[sentalgo][sentparameterstr].keys())
            elif sentalgo in ["netrate"]:
               threslist=["1"]
            else:
               print "this algo {0} is unknown".format(sentalgo)
               exit(1)
            for thres in threslist:
                if sentalgo in ["netinf","multitree"]:
                   if thres not in requirededgethres:
                      continue 
                elif sentalgo in ["lsel1","absel1"]:
                   if thres not in requiredlambdathres:
                      continue 
                elif sentalgo in ["netrate"]:
                   if thres not in "1":
                      continue 
                else:
                   print "sentalgo is unknown{0}".format(sentalgo)
                   exit(1)
                tempylist.append(sentscores[sentalgo][sentparameterstr][thres][ykey])
            avgscore=0.0
            for elem in tempylist:
                avgscore+=elem
            avgscore/=float(len(tempylist))
            if avgscore>=bestavgscore:
               bestavgscore=avgscore
               bestylist=list(tempylist)
               bestxlist=list(threslist)
               datalabel=plotlabelcreate(sentalgo,sentparameterstr)
        
        assert datalabel!="-1"
        plotxlist[sentalgo]=list(bestxlist)
        plotylist[sentalgo]=list(bestylist)
        labelinfo[sentalgo]=datalabel

    if plottype in ["bedroc","proc","roc(exp)"]: #we don't plot for this scoretypes!!
       return
   
    if len(labelinfo.keys())==0:
       print "ERROR: No data for plot having title:{0} !!".format(plottitle)
       exit(1)

    #output scores
    parts=plotfilepath.split("/")
    scoreoutfilepath="/".join(parts[0:-1])+"/scoreout_{0}".format(parts[-1].replace(".png",".txt"))
    file=open(scoreoutfilepath,"w")
    for sentalgo in plotscores.keys():
        file.write("{0}\t{1}\n".format(sentalgo,plotscores[sentalgo]))
    file.close()
    
    #plotting part
    plt.clf()
    plt.xlim([0,1.01])
    plt.ylim([0,1.01]) 
    plt.xlabel(myxlabel)
    plt.ylabel(myylabel)
    plt.title(plottitle)
    for index in range(0,len(plotxlist.keys())):
        datainfo=plotxlist.keys()[index] #datainfo=runalgo
        xlist=plotxlist[datainfo]
        ylist=plotylist[datainfo]
        plt.plot(xlist,ylist,marker=plotmarkers[index],linestyle='None',color=plotcolors[index],label=labelinfo[datainfo])
    plt.legend(loc=4)   
    plt.savefig(plotfilepath)


#auc curve generator for various parameters of the algorithm(edgecount or lambda)
#in case bedroc,proc and roc(exp), it only prints scores
#can not do bedroc right now since it is always increasing!!
def parameterauccurveplotgenerate(sentscores,plottype,plottitle,plotfilepath):
    maxormin="max"
    if plottype in ["roc","bedroc","proc","roc(exp)"]:
       xkey="fpr"
       ykey="sen"
       myxlabel="FPR"
       myylabel="Sensitivity"
    elif plottype in ["pr"]:
       xkey="recall"
       ykey="precision"
       myxlabel="Recall"
       myylabel="Precision"
    else:
      print "this plottype is unknown {0}".format(plottype)
      exit(1)
      
    print "Writing {0} auc PARAMETER SCORES for {1}".format(plottype,plottitle)
    #preprocessing part
    plotxlist={}
    plotylist={}
    plotscores={}
    labelinfo={}
    for sentalgo in sentscores.keys():
        bestaucscore=-1.0
        bestylist=[]
        bestxlist=[]
        datalabel="-1"
        #there is only single param for other algos since fraction is fixed
        #print "keysinfo: {0}".format(sentscores[sentalgo].keys())
        #if sentalgo in ["multitree","netinf","netrate","connie"]:
        #   assert len(sentscores[sentalgo].keys())==1
        for sentparameterstr in sentscores[sentalgo].keys():
            #parameterstr,extend=sentparameterstr.split("+")
            tempylist=[] #y
            tempxlist=[] #x
            
            if sentalgo in ["lsel1","absel1"]:
               threslist=sorted(sentscores[sentalgo][sentparameterstr].keys(),reverse=True)
            elif sentalgo in ["multitree","netinf"]:
               threslist=sorted(sentscores[sentalgo][sentparameterstr].keys())
            elif sentalgo in ["netrate"]:
               threslist=["1"]
            else:
               print "this algo {0} is unknown".format(sentalgo)
               exit(1)
            for thres in threslist:
                if sentalgo in ["netinf","multitree"]:
                   if thres not in requirededgethres:
                      continue 
                elif sentalgo in ["lsel1","absel1"]:
                   if thres not in requiredlambdathres:
                      continue 
                elif sentalgo in ["netrate"]:
                   if thres not in "1":
                      continue 
                else:
                   print "sentalgo is unknown{0}".format(sentalgo)
                   exit(1)
                tempylist.append(sentscores[sentalgo][sentparameterstr][thres][ykey])
                tempxlist.append(sentscores[sentalgo][sentparameterstr][thres][xkey])
            print "inside info:"    
            print sentalgo
            print sentparameterstr
            print threslist
            print tempylist
            print tempxlist
            paramaucscore=estimateaucscore(tempxlist,tempylist,plottype)
            if paramaucscore>=bestaucscore:
               bestaucscore=paramaucscore
               bestylist=list(tempylist)
               bestxlist=list(tempxlist)
               datalabel=plotlabelcreate(sentalgo,sentparameterstr)
            
        assert datalabel!="-1"
        plotxlist[sentalgo]=list(bestxlist)
        plotylist[sentalgo]=list(bestylist)
        plotscores[sentalgo]=bestaucscore
        labelinfo[sentalgo]=datalabel
        print "{0} {1}: {2}".format(plottype,sentalgo,bestaucscore)

    if plottype in ["bedroc","proc","roc(exp)"]: #we don't plot for this scoretypes!!
       return
   
    if len(labelinfo.keys())==0:
       print "ERROR: No data for plot having title:{0} !!".format(plottitle)
       exit(1)

    #output scores
    parts=plotfilepath.split("/")
    scoreoutfilepath="/".join(parts[0:-1])+"/scoreout_{0}".format(parts[-1].replace(".png",".txt"))
    file=open(scoreoutfilepath,"w")
    for datainfo in plotxlist.keys():
        xlist=plotxlist[datainfo]
        xliststr=str(xlist[0])
        for elem in xlist[1:]:
            xliststr+=",{0}".format(elem)
        ylist=plotylist[datainfo]
        yliststr=str(ylist[0])
        for elem in ylist[1:]:
            yliststr+=",{0}".format(elem)
        file.write("{0}\t{1}\t{2}\n".format(datainfo,xliststr,yliststr))
    file.close()
    
    #plotting part
    plt.clf()
    plt.xlim([0,1.01])
    plt.ylim([0,1.01]) 
    plt.xlabel(myxlabel)
    plt.ylabel(myylabel)
    plt.title(plottitle)
    for index in range(0,len(plotxlist.keys())):
        datainfo=plotxlist.keys()[index] #datainfo=runalgo
        xlist=plotxlist[datainfo]
        ylist=plotylist[datainfo]
        plt.plot(xlist,ylist,marker=plotmarkers[index],linestyle='None',color=plotcolors[index],label=labelinfo[datainfo])
    plt.legend(loc=4)   
    plt.savefig(plotfilepath)


#Plot generator over time(fraction)
#also plots roc curve and pr curve for changing fraction values
def overtimeplotgenerate(sentscores,plottype,plottitle,plotfilepath,xlabel,algoparaminfo,xaxis):
    ykey=plottype
    if plottype in ["sen","fpr","recall","precision","acc","f1","f2","mcc","spec","f1overn","f1div100","f1div10","f1div2","f10","f20","roc","pr"]:
       maxormin="max" 
    else:
       print "this plottype {0} is unknown for max or min !!".format(plottype)
       exit(1) 
    myxlabel=xlabel
    myylabel=plottype.capitalize()

    if plottype in ["roc","pr"]:
       if plottype in ["roc"]:
          subykey="sen"
          subxkey="fpr"
       elif plottype in ["pr"]:
          subykey="precision"
          subxkey="recall"
       sentscores=sentscores2paramwise(sentscores)
       plotxlist={}
       plotylist={}
       labelinfo={}
       for runalgo in sentscores.keys():
           print runalgo
           maxvalue=-1.0
           maxlabel="-1.0"
           bestxlist=[]
           bestylist=[]
           for sentparamstr in sentscores[runalgo].keys():
               print "here"
               overtimexlist=[]
               overtimeylist=[]
               for fraction in sorted(sentscores[runalgo][sentparamstr].keys()):
                   overtimexlist.append(fraction)
                   thresylist=[] #either sen or precision
                   thresxlist=[] #either fpr or recall
                   if runalgo in ["lsel1","absel1"]:
                      threslist=sorted(sentscores[runalgo][sentparamstr][fraction].keys(),reverse=True)
                   elif runalgo in ["lsel1fused","absel1fused"]: #our algos dynamic
                      threslist=sorted(sentscores[runalgo][sentparamstr][fraction].keys(),reverse=True) 
                   elif runalgo in ["multitree","netinf"]:
                      threslist=sorted(sentscores[runalgo][sentparamstr][fraction].keys())
                   elif runalgo in ["netrate"]:
                      threslist=["1"]
                   else:
                      print "this algo {0} is unknown".format(runalgo)
                      exit(1)
                   for thres in threslist:
                       if runalgo in ["netinf","multitree"]:
                          if thres not in requirededgethres:
                             continue 
                       elif runalgo in ["lsel1","absel1"]:
                          if thres not in requiredlambdathres:
                             continue 
                       elif runalgo in ["netrate"]:
                          if thres not in "1":
                             continue
                       elif runalgo in ["lsel1fused","absel1fused"]: #our algos dynamic
                          if thres not in requiredlambdathres:
                             continue 
                       else:
                          print "sentalgo is unknown{0}".format(runalgo)
                          exit(1)
                       thresylist.append(sentscores[runalgo][sentparamstr][fraction][thres][subykey])
                       thresxlist.append(sentscores[runalgo][sentparamstr][fraction][thres][subxkey])
                   tempscore=estimateaucscore(thresxlist,thresylist,plottype)
                   overtimeylist.append(tempscore)
               areavalue=estimatearea(overtimexlist,overtimeylist)
               if areavalue>=maxvalue:
                  maxvalue=areavalue
                  maxlabel=plotlabelcreate(runalgo,sentparamstr)
                  bestxlist=list(overtimexlist)
                  bestylist=list(overtimeylist)
           assert maxvalue!=-1.0
           plotxlist[runalgo]=bestxlist
           plotylist[runalgo]=bestylist
           labelinfo[runalgo]=maxlabel

       #plotting part
       parts=plotfilepath.split("/")
       partialplotfolder="/".join(parts[0:-1])
       partialplotfilepath="{0}/area_{1}_{2}".format(partialplotfolder,plottype,parts[-1])
       plt.clf()
       plt.xlim([0,max(requiredfractions)+0.05])
       plt.ylim([0,1.0]) 
       plt.xlabel(myxlabel)
       plt.ylabel(myylabel)
       plt.title(plottitle)
       for index in range(0,len(plotxlist.keys())):
           datainfo=plotxlist.keys()[index] 
           xlist=plotxlist[datainfo]
           ylist=plotylist[datainfo]
           plt.plot(xlist,ylist,marker=plotmarkers[index],linestyle='None',color=plotcolors[index],label=labelinfo[datainfo])
       plt.legend(loc=4)   
       plt.savefig(partialplotfilepath)

       #output scores
       parts=plotfilepath.split("/")
       scoreoutfilepath="/".join(parts[0:-1])+"/scoreout_area_{0}".format(parts[-1].replace(".png",".txt"))
       file=open(scoreoutfilepath,"w")
       for datainfo in plotxlist.keys():
           xlist=plotxlist[datainfo]
           xliststr=str(xlist[0])
           for elem in xlist[1:]:
               xliststr+=",{0}".format(elem)
           ylist=plotylist[datainfo]
           yliststr=str(ylist[0])
           for elem in ylist[1:]:
               yliststr+=",{0}".format(elem)
           file.write("{0}\t{1}\t{2}\n".format(datainfo,xliststr,yliststr))
       file.close()
       
       return
   
    if algoparaminfo!="-1":
       plotxlist={}
       plotylist={}
       labelinfo={} 
       for runalgo in algoparaminfo.keys():
           for algoparam in algoparaminfo[runalgo]:
               print "top info"
               print runalgo
               print algoparam
               print type(algoparam)
               foundparamstr="--"
               if not sentscores.has_key(runalgo):
                  continue
               #print "starting {0} {1}".format(runalgo,algoparam)
               for sentparameterstr in sentscores[runalgo].keys():
                   parameterstr,extend=sentparameterstr.split("+")
                   parts=parameterstr.replace("result_","").split("_")
                   print "mypart: {0}".format(parts)
                   if runalgo in ["netinf","multitree"]:
                      if parts[-1].find("categorical")!=-1:
                         foundone="categorical?" 
                      elif parts[-1].find("error")!=-1: 
                         foundone="error?" 
                      else:
                         print "not known rounding method!!"
                         exit(1)
                      tempparam=float(parts[-1].replace(foundone,""))
                   elif runalgo in ["netrate"]:
                      tempparam="1" #there won't be param for netrate
                   elif runalgo in ["lsel1","absel1"]: #our algos
                      tempparam=float(parts[-1].split("?")[0])
                   elif runalgo in ["lsel1fused","absel1fused"]: #our algos dynamic
                      tempparam=float(parts[-1].split("?")[0])
                      print "temp {0}".format(tempparam)
                   else:
                      print "this algo is unknoawn {0}".format(runalgo)
                      exit(1)
                   #result_uniform_undirected_0.6-randomnoise_epsilon-0.9_0.001??2.0?epsilon0.01?flowdifference_3
                   print "{0} -> {1}".format(runalgo,tempparam)
                   print "tempparam: {0}".format(tempparam)
                   print type(tempparam)
                   if tempparam==algoparam:
                      print "here breaking!!"
                      foundparamstr=sentparameterstr
                      break
               #print "lastinfo:"    
               #print runalgo
               #print algoparam
               assert foundparamstr!="--"
               xlist=[]
               ylist=[]
               for fraction in xaxis:    
                   if not sentscores[runalgo][foundparamstr].has_key(fraction):
                      continue
                   ylist.append(sentscores[runalgo][foundparamstr][fraction][plottype])
                   xlist.append(float(fraction))
               datalabel=plotlabelcreate(runalgo,foundparamstr,"algoparam")
               plotxlist[datalabel]=xlist
               plotylist[datalabel]=ylist
               labelinfo[datalabel]=datalabel

       #plotting part
       parts=plotfilepath.split("/")
       partialplotfolder="/".join(parts[0:-1])
       partialplotfilepath="{0}/partial_{1}".format(partialplotfolder,parts[-1])
       plt.clf()
       plt.xlim([0,max(xaxis)+0.05])
       plt.ylim([0,1.0]) 
       plt.xlabel(myxlabel)
       plt.ylabel(myylabel)
       plt.title(plottitle)
       for index in range(0,len(plotxlist.keys())):
           datainfo=plotxlist.keys()[index] 
           xlist=plotxlist[datainfo]
           ylist=plotylist[datainfo]
           plt.plot(xlist,ylist,marker=plotmarkers[index],linestyle='None',color=plotcolors[index],label=labelinfo[datainfo])
       plt.legend(loc=4)   
       plt.savefig(partialplotfilepath)

       #output scores
       parts=plotfilepath.split("/")
       scoreoutfilepath="/".join(parts[0:-1])+"/scoreout_partial_{0}".format(parts[-1].replace(".png",".txt"))
       file=open(scoreoutfilepath,"w")
       for datainfo in plotxlist.keys():
           xlist=plotxlist[datainfo]
           xliststr=str(xlist[0])
           for elem in xlist[1:]:
               xliststr+=",{0}".format(elem)
           ylist=plotylist[datainfo]
           yliststr=str(ylist[0])
           for elem in ylist[1:]:
               yliststr+=",{0}".format(elem)
           file.write("{0}\t{1}\t{2}\n".format(datainfo,xliststr,yliststr))
       file.close()
       return 
    
    #preprocessing part
    bestplotxlist={}
    bestplotylist={}
    worstplotxlist={}
    worstplotylist={}
    bestlabelinfo={}
    worstlabelinfo={}
    plotmeanlist={}
    plotstddevlist={}
    for runalgo in sentscores.keys():
        if maxormin=="max":
           bestscore=-1.0
           worstscore=100.0
        elif maxormin=="min":
           bestscore=100.0
           worstscore=-1.0
        bestylist=[]
        bestxlist=[]
        allylist=[] #all y list scores
        worstxlist=[]
        worstylist=[]
        #bestsen=[]
        #bestfpr=[]
        #bestaucscore=-1.0
        bestdatalabel="-1"
        worstdatalabel="-1"
        for sentparameterstr in sentscores[runalgo].keys():
            parameterstr,extend=sentparameterstr.split("+")
            #tempsenlist=[]
            #tempfprlist=[]
            #for fraction in sorted(sentscores[runalgo][sentparameterstr].keys()):
            #    if fraction not in requiredfractions:
            #       continue 
            #    tempsenlist.append(sentscores[runalgo][sentparameterstr][fraction]["sen"])
            #    tempfprlist.append(sentscores[runalgo][sentparameterstr][fraction]["fpr"])
            #paramaucscore=estimateaucscore(tempfprlist,tempsenlist,"roc")
            #if paramaucscore>=bestaucscore:
            #   bestaucscore=paramaucscore
            
            ylist=[]
            xlist=[]
            #for fraction in sorted(sentscores[runalgo][sentparameterstr].keys()):
            for fraction in xaxis:    
                #if fraction not in requiredfractions:
                #   continue
                if not sentscores[runalgo][sentparameterstr].has_key(fraction):
                   continue
                ylist.append(sentscores[runalgo][sentparameterstr][fraction][plottype])
                xlist.append(float(fraction))
            allylist.append(ylist)    
            paramareascore=estimatearea(xlist,ylist)
            #if plottype=="precision":
            #   print "area {0}".format(paramareascore)
            #   print "param: {0}".format(sentparameterstr)
            if maxormin=="max":
               if paramareascore>=bestscore:
                  bestscore=paramareascore
                  bestylist=list(ylist)
                  bestxlist=list(xlist)
                  bestdatalabel=plotlabelcreate(runalgo,sentparameterstr)
               if paramareascore<=worstscore:
                  worstscore=paramareascore
                  worstylist=list(ylist)
                  worstxlist=list(xlist)
                  worstdatalabel=plotlabelcreate(runalgo,sentparameterstr)   
            elif maxormin=="min":
               if paramareascore<=bestscore:
                  bestscore=paramareascore
                  bestylist=list(ylist)
                  bestxlist=list(xlist)
                  bestdatalabel=plotlabelcreate(runalgo,sentparameterstr)
               if paramareascore>=worstscore:
                  worstscore=paramareascore
                  worstylist=list(ylist)
                  worstxlist=list(xlist)
                  worstdatalabel=plotlabelcreate(runalgo,sentparameterstr)     

        meanlist=[]
        stddevlist=[]
        for fracindex in range(0,len(xaxis)):   
            locallist=[]
            for elem in allylist:
                if fracindex>=len(elem):
                   locallist.append(0.0)
                else:   
                   locallist.append(elem[fracindex])
            mean=np.mean(locallist)
            stddev=np.std(locallist)
            meanlist.append(mean)
            stddevlist.append(stddev)
            
        assert bestdatalabel!="-1"
        assert worstdatalabel!="-1"
        bestplotxlist[runalgo]=list(bestxlist)
        bestplotylist[runalgo]=list(bestylist)
        worstplotxlist[runalgo]=list(worstxlist)
        worstplotylist[runalgo]=list(worstylist)
        plotmeanlist[runalgo]=list(meanlist)
        plotstddevlist[runalgo]=list(stddevlist)
        bestlabelinfo[runalgo]=bestdatalabel
        worstlabelinfo[runalgo]=worstdatalabel

        #print "runalgo auc score {0}: {1}".format(runalgo,bestaucscore)
             
    if len(bestlabelinfo.keys())==0 or len(worstlabelinfo.keys())==0:
       print "No data for plot for title:{0} having path {1} !!".format(plottitle,plotfilepath)
       return

    #best performing plotting part
    parts=plotfilepath.split("/")
    bestplotfolder="/".join(parts[0:-1])
    bestplotfilepath="{0}/best_{1}".format(bestplotfolder,parts[-1])
    plt.clf()
    plt.xlim([0,max(xaxis)+0.05])
    plt.ylim([0,1.0]) 
    plt.xlabel(myxlabel)
    plt.ylabel(myylabel)
    plt.title(plottitle)
    for index in range(0,len(bestplotxlist.keys())):
        datainfo=bestplotxlist.keys()[index] 
        xlist=bestplotxlist[datainfo]
        ylist=bestplotylist[datainfo]
        #plt.plot(xlist,ylist,marker=plotmarkers[index],linestyle='--',color=plotcolors[index],label=labelinfo[index])
        plt.plot(xlist,ylist,marker=plotmarkers[index],linestyle='None',color=plotcolors[index],label=bestlabelinfo[datainfo])
    plt.legend(loc=4)   
    plt.savefig(bestplotfilepath)
    #output scores best performing
    parts=plotfilepath.split("/")
    scoreoutfilepath="/".join(parts[0:-1])+"/scoreout_best_{0}".format(parts[-1].replace(".png",".txt"))
    file=open(scoreoutfilepath,"w")
    for datainfo in bestplotxlist.keys():
        xlist=bestplotxlist[datainfo]
        xliststr=str(xlist[0])
        for elem in xlist[1:]:
            xliststr+=",{0}".format(elem)
        ylist=bestplotylist[datainfo]
        yliststr=str(ylist[0])
        for elem in ylist[1:]:
            yliststr+=",{0}".format(elem)
        file.write("{0}\t{1}\t{2}\n".format(datainfo,xliststr,yliststr))
    file.close()
       
    #worst performing plotting part
    parts=plotfilepath.split("/")
    worstplotfolder="/".join(parts[0:-1])
    worstplotfilepath="{0}/worst_{1}".format(worstplotfolder,parts[-1])
    plt.clf()
    plt.xlim([0,max(xaxis)+0.05])
    plt.ylim([0,1.0]) 
    plt.xlabel(myxlabel)
    plt.ylabel(myylabel)
    plt.title(plottitle)
    for index in range(0,len(worstplotxlist.keys())):
        datainfo=worstplotxlist.keys()[index] 
        xlist=worstplotxlist[datainfo]
        ylist=worstplotylist[datainfo]
        #plt.plot(xlist,ylist,marker=plotmarkers[index],linestyle='--',color=plotcolors[index],label=labelinfo[index])
        plt.plot(xlist,ylist,marker=plotmarkers[index],linestyle='None',color=plotcolors[index],label=worstlabelinfo[datainfo])
    plt.legend(loc=4)   
    plt.savefig(worstplotfilepath)
    #output scores worst performing
    parts=plotfilepath.split("/")
    scoreoutfilepath="/".join(parts[0:-1])+"/scoreout_worst_{0}".format(parts[-1].replace(".png",".txt"))
    file=open(scoreoutfilepath,"w")
    for datainfo in worstplotxlist.keys():
        xlist=worstplotxlist[datainfo]
        xliststr=str(xlist[0])
        for elem in xlist[1:]:
            xliststr+=",{0}".format(elem)
        ylist=worstplotylist[datainfo]
        yliststr=str(ylist[0])
        for elem in ylist[1:]:
            yliststr+=",{0}".format(elem)
        file.write("{0}\t{1}\t{2}\n".format(datainfo,xliststr,yliststr))
    file.close()

    if False:
     #average performing plotting part
     parts=plotfilepath.split("/")
     meanplotfolder="/".join(parts[0:-1])
     meanplotfilepath="{0}/mean_{1}".format(meanplotfolder,parts[-1])
     plt.clf()
     if type(xaxis)!=str:
        plt.xlim([0,max(xaxis)+0.05])
     plt.ylim([0,1.0]) 
     plt.xlabel(myxlabel)
     plt.ylabel(myylabel)
     plt.title(plottitle)
     for index in range(0,len(bestplotxlist.keys())):
        datainfo=bestplotxlist.keys()[index] 
        xlist=bestplotxlist[datainfo]
        meanylist=plotmeanlist[datainfo]
        #plt.plot(xlist,ylist,marker=plotmarkers[index],linestyle='--',color=plotcolors[index],label=labelinfo[index])
        plt.plot(xlist,meanylist,marker=plotmarkers[index],linestyle='None',color=plotcolors[index],label=bestlabelinfo[datainfo])
     plt.legend(loc=4)   
     plt.savefig(meanplotfilepath)
     #output scores mean performing
     parts=plotfilepath.split("/")
     scoreoutfilepath="/".join(parts[0:-1])+"/scoreout_mean_{0}".format(parts[-1].replace(".png",".txt"))
     file=open(scoreoutfilepath,"w")
     for datainfo in bestplotxlist.keys():
        xlist=bestplotxlist[datainfo]
        xliststr=str(xlist[0])
        for elem in xlist[1:]:
            xliststr+=",{0}".format(elem)
        ylist=plotmeanlist[datainfo]
        yliststr=str(ylist[0])
        for elem in ylist[1:]:
            yliststr+=",{0}".format(elem)
        file.write("{0}\t{1}\t{2}\n".format(datainfo,xliststr,yliststr))
     file.close()
    
     #error plotting part
     parts=plotfilepath.split("/")
     errorplotfolder="/".join(parts[0:-1])
     errorplotfilepath="{0}/error_{1}".format(errorplotfolder,parts[-1])
     plt.clf()
     if type(xaxis)!=str:
        plt.xlim([0,max(xaxis)+0.05])
     plt.ylim([0,1.0]) 
     plt.xlabel(myxlabel)
     plt.ylabel(myylabel)
     plt.title(plottitle)
     for index in range(0,len(bestplotxlist.keys())):
        datainfo=bestplotxlist.keys()[index] #datainfo=runalgo
        xlist=bestplotxlist[datainfo]
        meanylist=plotmeanlist[datainfo]
        stddevlist=plotstddevlist[datainfo]
        plt.errorbar(xlist, meanylist,yerr=stddevlist, fmt=plotmarkers[index], color=plotcolors[index],ecolor=plotcolors[index],label=bestlabelinfo[datainfo])
        #plt.errorbar(xlist, meanylist,fmt='ro', ecolor='g',yerr=stddevlist, marker=plotmarkers[index],linestyle='--',color=plotcolors[index],label=labelinfo[datainfo]) 
     plt.legend(loc=4)   
     plt.savefig(errorplotfilepath)
     #output scores errorbar performing
     parts=plotfilepath.split("/")
     scoreoutfilepath="/".join(parts[0:-1])+"/scoreout_error_{0}".format(parts[-1].replace(".png",".txt"))
     file=open(scoreoutfilepath,"w")
     for datainfo in bestplotxlist.keys():
        xlist=bestplotxlist[datainfo]
        xliststr=str(xlist[0])
        for elem in xlist[1:]:
            xliststr+=",{0}".format(elem)
        ylist=plotmeanlist[datainfo]
        yliststr=str(ylist[0])
        for elem in ylist[1:]:
            yliststr+=",{0}".format(elem)
        zlist=plotstddevlist[datainfo]
        zliststr=str(zlist[0])
        for elem in zlist[1:]:
            zliststr+=",{0}".format(elem)    
        file.write("{0}\t{1}\t{2}\t{3}\n".format(datainfo,xliststr,yliststr,zliststr))
     file.close()
       
#does everything except early one!!
#we fix number of traces!!    
#def noiseplotgenerate(sentscores,plottype,plottitle,plotfilepath,xlabel,fixedfraction):
def noiseplotgenerate(sentscores,plottype,plottitle,plotfilepath,xlabel):    
    #f1,f2 kind of scores will be average scores over number of traces(time)
    assert plottype in ["f1","f2","f1overn","f1div100","f1div10","f1div2","f10","f20","roc-trace","roc-param","pr-trace","pr-param","bep-trace","bep-param"] 
    
    firstpart,detail=plottype.split("-")
    maxormin="max" 
    myxlabel=xlabel
    myylabel=plottype.capitalize()
    
    
    ykey=plottype
    maxormin="max" 
    myxlabel=xlabel
    myylabel=plottype.capitalize()
    
    #preprocessing part
    plotxlist={}
    plotylist={}
    labelinfo={}
    for sentalgo in sentscores.keys():
        datalabel="-1"
        bestylist=[]
        bestxlist=[]
        if maxormin=="max": #can be either area or average f
           bestareascore=-1.0
        elif maxormin=="min":
           bestareascore=100.0
        for sentparameterstr in sentscores[sentalgo].keys():
            ylist=[]
            xlist=sorted(sentscores[sentalgo][sentparameterstr].keys())
            for noiselevel in sorted(sentscores[sentalgo][sentparameterstr].keys()):
                threslist=sorted(sentscores[sentalgo][sentparameterstr][noiselevel].keys())
                tempxlist=[]
                tempylist=[]
                for thres in threslist:
                    tempxlist.append(float(sentscores[sentalgo][sentparameterstr][noiselevel][thres]["fpr"]))
                    tempylist.append(float(sentscores[sentalgo][sentparameterstr][noiselevel][thres]["sen"]))
                paramaucscore=estimateaucscore(tempxlist,tempylist,"roc")
                ylist.append(paramaucscore)
            paramareascore=estimatearea(xlist,ylist)
            if paramareascore>=bestareascore:
               bestareascore=paramareascore
               bestxlist=list(xlist)
               bestylist=list(ylist)
               datalabel=plotlabelcreate(sentalgo,sentparameterstr)
             
        assert datalabel!="-1"
        plotxlist[sentalgo]=list(bestxlist)
        plotylist[sentalgo]=list(bestylist)
        labelinfo[sentalgo]=datalabel

    if len(labelinfo.keys())==0:
       print "No data for plot algo:{0} for title:{1} !!".format(runalgo,plottitle)
       return

    xmax=-1.0
    xmin=14000.0
    seenxs=set()
    for sentalgo in plotxlist.keys():
        for elem in plotxlist[sentalgo]:
           if elem<xmin:
              xmin=elem
           if elem>xmax:
              xmax=elem   
           seenxs.add(elem)
           
    #plotting part
    plt.clf()
    #plt.set_xticks(seenxs)
    plt.xlim([-0.5+xmin,xmax+0.5])
    plt.ylim([0,1.0]) 
    plt.xlabel(myxlabel)
    plt.ylabel(myylabel)
    plt.title(plottitle)
    for index in range(0,len(plotxlist.keys())):
        datainfo=plotxlist.keys()[index] #datainfo=runalgo
        xlist=plotxlist[datainfo]
        ylist=plotylist[datainfo]
        #plt.plot(xlist,ylist,marker=plotmarkers[index],linestyle='--',color=plotcolors[index],label=labelinfo[index])
        plt.plot(xlist,ylist,marker=plotmarkers[index],linestyle='None',color=plotcolors[index],label=labelinfo[datainfo])
    plt.legend(loc=1)   
    plt.savefig(plotfilepath)

#sampling interval plot generator(right now assumes plottype is always "roc-trace")
def samplingplotgenerate(sentscores,plottype,plottitle,plotfilepath,xlabel):    
    ykey=plottype
    maxormin="max" 
    myxlabel=xlabel
    myylabel=plottype.capitalize()
    
    #preprocessing part
    plotxlist={}
    plotylist={}
    labelinfo={}
    for sentalgo in sentscores.keys():
        datalabel="-1"
        bestylist=[]
        bestxlist=[]
        if maxormin=="max": #can be either area or average f
           bestareascore=-1.0
        elif maxormin=="min":
           bestareascore=100.0
        for sentparameterstr in sentscores[sentalgo].keys():
            ylist=[]
            xlist=sorted(sentscores[sentalgo][sentparameterstr].keys())
            for samplelevel in sorted(sentscores[sentalgo][sentparameterstr].keys()):
                threslist=sorted(sentscores[sentalgo][sentparameterstr][samplelevel].keys())
                tempxlist=[]
                tempylist=[]
                for thres in threslist:
                    tempxlist.append(float(sentscores[sentalgo][sentparameterstr][samplelevel][thres]["fpr"]))
                    tempylist.append(float(sentscores[sentalgo][sentparameterstr][samplelevel][thres]["sen"]))
                paramaucscore=estimateaucscore(tempxlist,tempylist,"roc")
                ylist.append(paramaucscore)
            paramareascore=estimatearea(xlist,ylist)
            if paramareascore>=bestareascore:
               bestareascore=paramareascore
               bestxlist=list(xlist)
               bestylist=list(ylist)
               datalabel=plotlabelcreate(sentalgo,sentparameterstr)
             
        assert datalabel!="-1"
        plotxlist[sentalgo]=list(bestxlist)
        plotylist[sentalgo]=list(bestylist)
        labelinfo[sentalgo]=datalabel

    if len(labelinfo.keys())==0:
       print "No data for plot algo:{0} for title:{1} !!".format(runalgo,plottitle)
       return

    xmax=-1.0
    xmin=12000.0
    seenxs=set()
    for sentalgo in plotxlist.keys():
        for elem in plotxlist[sentalgo]:
           if elem<xmin:
              xmin=elem
           if elem>xmax:
              xmax=elem   
           seenxs.add(elem)
           
    #plotting part
    plt.clf()
    #plt.set_xticks(seenxs)
    plt.xlim([-0.5+xmin,xmax+0.5])
    plt.ylim([0,1.0]) 
    plt.xlabel(myxlabel)
    plt.ylabel(myylabel)
    plt.title(plottitle)
    for index in range(0,len(plotxlist.keys())):
        datainfo=plotxlist.keys()[index] #datainfo=runalgo
        xlist=plotxlist[datainfo]
        ylist=plotylist[datainfo]
        #plt.plot(xlist,ylist,marker=plotmarkers[index],linestyle='--',color=plotcolors[index],label=labelinfo[index])
        plt.plot(xlist,ylist,marker=plotmarkers[index],linestyle='None',color=plotcolors[index],label=labelinfo[datainfo])
    plt.legend(loc=1)   
    plt.savefig(plotfilepath)


#prob dist plot generator, variable trace
def probdist_variabletraceplotgenerate(sentscores2,plottype,plottitle,plotfilepath,xlabel):
    ykey=plottype
    assert plottype in ["roc-trace","pr-trace","bep-trace"]
    firstpart,detail=plottype.split("-")
    maxormin="max" 
    myxlabel=xlabel
    myylabel=plottype.capitalize()
    
    #preprocessing part
    plotxlist={}
    plotylist={}
    plotmeanlist={}
    plotstddevlist={}
    labelinfo={}
    for sentalgo in sentscores2.keys():
        datalabel="-1"
        bestylist=[]
        bestxlist=[]
        if maxormin=="max": #can be either area or average f
           bestareascore=-1.0
        elif maxormin=="min":
           bestareascore=100.0
        for sentparameterstr in sentscores2[sentalgo].keys():
            ylist=[]
            xlist=sorted(sentscores2[sentalgo][sentparameterstr].keys())
            for plotdistinfo in sorted(sentscores2[sentalgo][sentparameterstr].keys()):
                if sentalgo in ["lsel1","absel1"]:
                   threslist=sorted(sentscores2[sentalgo][sentparameterstr][plotdistinfo].keys(),reverse=True)
                elif sentalgo in ["multitree","netinf"]:
                   threslist=sorted(sentscores2[sentalgo][sentparameterstr][plotdistinfo].keys())
                elif sentalgo in ["netrate"]:
                   pass
                else:
                   print "this algo {0} is unknown".format(sentalgo)
                   exit(1)
                if firstpart in ["f1","f2","f1overn","f1div100","f1div10","f1div2","f10","f20"]:
                   if detail=="mean":
                      avgscore=0.0   
                      for thres in threslist: #!!!!! WE MUST MAKE SURE EVERYONE HAS USED THE SAME AMOUNT OF TRACE!!
                         avgscore+=float(sentscores2[sentalgo][sentparameterstr][plotdistinfo][thres][ykey])
                      avgscore/=float(len(threslist))
                      sentscore=avgscore
                   elif detail=="error":   
                      tempscores=[]
                      for thres in threslist: #!!!!! WE MUST MAKE SURE EVERYONE HAS USED THE SAME AMOUNT OF TRACE!!
                         tempscore=float(sentscores2[sentalgo][sentparameterstr][plotdistinfo][thres][ykey])
                         tempscores.append(tempscore)
                      meanscore=np.mean(tempscores)
                      stddevscore=np.std(tempscores)
                      sentscore=(meanscore,stddevscore)
                   elif detail=="median":
                      tempscores=[]
                      for thres in threslist: #!!!!! WE MUST MAKE SURE EVERYONE HAS USED THE SAME AMOUNT OF TRACE!!
                         tempscore=float(sentscores2[sentalgo][sentparameterstr][plotdistinfo][thres][ykey])
                         tempscores.append(tempscore) 
                      medianscore=np.median(tempscores)
                      sentscore=medianscore
                   else:
                      print "unknown detail {0}".format(detail)
                      exit(1)
                ylist.append(sentscore)
    
            paramareascore=estimatearea(xlist,ylist)
            if paramareascore>=bestareascore:
               bestareascore=paramaucscore
               bestxlist=list(xlist)
               bestylist=list(ylist)
               datalabel=plotlabelcreate(sentalgo,sentparameterstr)   
            
        assert datalabel!="-1"
        plotxlist[sentalgo]=list(bestxlist)
        plotylist[sentalgo]=list(bestylist)
        labelinfo[sentalgo]=datalabel

    if len(labelinfo.keys())==0:
       print "No data for plot algo:{0} for title:{1} !!".format(runalgo,plottitle)
       return

    if False:
     xmax=-1.0
     xmin=12000.0
     seenxs=set()
     for sentalgo in plotxlist.keys():
        for elem in plotxlist[sentalgo]:
           if elem<xmin:
              xmin=elem
           if elem>xmax:
              xmax=elem   
           seenxs.add(elem)
           
     #plotting part
     plt.clf()
     #plt.set_xticks(seenxs)
     plt.xlim([-0.5+xmin,xmax+0.5])
    
    #plotting part
    plt.clf()
    #plt.xlim([0,0.9])
    plt.ylim([0,1.0]) 
    plt.xlabel(myxlabel)
    plt.ylabel(myylabel)
    plt.title(plottitle)
    for index in range(0,len(plotxlist.keys())):
        datainfo=plotxlist.keys()[index] #datainfo=runalgo
        xlist=plotxlist[datainfo]
        ylist=plotylist[datainfo]
        #plt.plot(xlist,ylist,marker=plotmarkers[index],linestyle='--',color=plotcolors[index],label=labelinfo[index])
        plt.plot(xlist,ylist,marker=plotmarkers[index],linestyle='None',color=plotcolors[index],label=labelinfo[datainfo])
    plt.legend(loc=1)   
    plt.savefig(plotfilepath)


#prob dist plot generator, fixed trace
def probdist_fixedtraceplotgenerate(sentscores2,plottype,plottitle,plotfilepath,xlabel):
    ykey=plottype
    assert plottype in ["f1","f1overn","f1div100","f1div10","f1div2","f20","f50","f100","roc-param","pr-param","bep-param"]
    firstpart,detail=plottype.split("-")
    maxormin="max" 
    myxlabel=xlabel
    myylabel=plottype.capitalize()
    
    #preprocessing part
    plotxlist={}
    plotylist={}
    plotmeanlist={}
    plotstddevlist={}
    labelinfo={}
    for sentalgo in sentscores2.keys():
        datalabel="-1"
        bestylist=[]
        bestxlist=[]
        if maxormin=="max": #can be either area or average f
           bestareascore=-1.0
        elif maxormin=="min":
           bestareascore=100.0
        for sentparameterstr in sentscores2[sentalgo].keys():
            ylist=[]
            xlist=sorted(sentscores2[sentalgo][sentparameterstr].keys())
            for plotdistinfo in sorted(sentscores2[sentalgo][sentparameterstr].keys()):
                if sentalgo in ["lsel1","absel1"]:
                   threslist=sorted(sentscores2[sentalgo][sentparameterstr][plotdistinfo].keys(),reverse=True)
                elif sentalgo in ["multitree","netinf"]:
                   threslist=sorted(sentscores2[sentalgo][sentparameterstr][plotdistinfo].keys())
                elif sentalgo in ["netrate"]:
                   pass
                else:
                   print "this algo {0} is unknown".format(sentalgo)
                   exit(1)
                if firstpart in ["f1","f2","f1overn","f1div100","f1div10","f1div2","f10","f20"]:
                   if detail=="mean":
                      avgscore=0.0   
                      for thres in threslist: #!!!!! WE MUST MAKE SURE EVERYONE HAS USED THE SAME AMOUNT OF TRACE!!
                         avgscore+=float(sentscores2[sentalgo][sentparameterstr][plotdistinfo][thres][ykey])
                      avgscore/=float(len(threslist))
                      sentscore=avgscore
                   elif detail=="error":   
                      tempscores=[]
                      for thres in threslist: #!!!!! WE MUST MAKE SURE EVERYONE HAS USED THE SAME AMOUNT OF TRACE!!
                         tempscore=float(sentscores2[sentalgo][sentparameterstr][plotdistinfo][thres][ykey])
                         tempscores.append(tempscore)
                      meanscore=np.mean(tempscores)
                      stddevscore=np.std(tempscores)
                      sentscore=(meanscore,stddevscore)
                   elif detail=="median":
                      tempscores=[]
                      for thres in threslist: #!!!!! WE MUST MAKE SURE EVERYONE HAS USED THE SAME AMOUNT OF TRACE!!
                         tempscore=float(sentscores2[sentalgo][sentparameterstr][plotdistinfo][thres][ykey])
                         tempscores.append(tempscore) 
                      medianscore=np.median(tempscores)
                      sentscore=medianscore
                   else:
                      print "unknown detail {0}".format(detail)
                      exit(1)
                ylist.append(sentscore)
    
            paramareascore=estimatearea(xlist,ylist)
            if paramareascore>=bestareascore:
               bestareascore=paramaucscore
               bestxlist=list(xlist)
               bestylist=list(ylist)
               datalabel=plotlabelcreate(sentalgo,sentparameterstr)   
            
        assert datalabel!="-1"
        plotxlist[sentalgo]=list(bestxlist)
        plotylist[sentalgo]=list(bestylist)
        labelinfo[sentalgo]=datalabel

    if len(labelinfo.keys())==0:
       print "No data for plot algo:{0} for title:{1} !!".format(runalgo,plottitle)
       return

    if False:
     xmax=-1.0
     xmin=12000.0
     seenxs=set()
     for sentalgo in plotxlist.keys():
        for elem in plotxlist[sentalgo]:
           if elem<xmin:
              xmin=elem
           if elem>xmax:
              xmax=elem   
           seenxs.add(elem)
           
     #plotting part
     plt.clf()
     #plt.set_xticks(seenxs)
     plt.xlim([-0.5+xmin,xmax+0.5])
    
    #plotting part
    plt.clf()
    #plt.xlim([0,0.9])
    plt.ylim([0,1.0]) 
    plt.xlabel(myxlabel)
    plt.ylabel(myylabel)
    plt.title(plottitle)
    for index in range(0,len(plotxlist.keys())):
        datainfo=plotxlist.keys()[index] #datainfo=runalgo
        xlist=plotxlist[datainfo]
        ylist=plotylist[datainfo]
        #plt.plot(xlist,ylist,marker=plotmarkers[index],linestyle='--',color=plotcolors[index],label=labelinfo[index])
        plt.plot(xlist,ylist,marker=plotmarkers[index],linestyle='None',color=plotcolors[index],label=labelinfo[datainfo])
    plt.legend(loc=1)   
    plt.savefig(plotfilepath)


#not used
def varparam_fixedtraceplotgenerate_notused(sentscores,plottype,plottitle,plotfilepath,xlabel):
    ykey=plottype
    assert plottype in ["f1","f1overn","f1div100","f1div10","f1div2","f20","f50","f100"]
    maxormin="max" 
    myxlabel=xlabel
    myylabel=plottype.capitalize()
    
    #preprocessing part
    plotxlist={}
    plotylist={}
    plotmeanlist={}
    plotstddevlist={}
    labelinfo={}
    for sentalgo in sentscores.keys():
        datalabel="-1"
        bestylist=[]
        bestxlist=[]
        bestylist=[]
        bestxlist=[]
        bestylist=[]
        bestxlist=[]
        if maxormin=="max": #can be either area or average f
           bestareascore=-1.0
        elif maxormin=="min":
           bestareascore=100.0
        for sentparameterstr in sentscores[sentalgo].keys():
            ylist=[]
            xlist=sorted(sentscores[sentalgo][sentparameterstr].keys())
            for varparam in xlist:
                ylist.append(sentscores[sentalgo][sentparameterstr][varparam])
            #max score

            #mean score    
            if detail=="mean":
                avgscore=0.0   
                for thres in threslist:
                    avgscore+=float(sentscores[sentalgo][sentparameterstr][varparam][thres][ykey])
                avgscore/=float(len(threslist))
                sentscore=avgscore
            elif detail=="error":   
                tempscores=[]
                for thres in threslist:
                    tempscore=float(sentscores[sentalgo][sentparameterstr][varparam][thres][ykey])
                    tempscores.append(tempscore)
                meanscore=np.mean(tempscores)
                stddevscore=np.std(tempscores)
                sentscore=(meanscore,stddevscore)
            elif detail=="median":
                tempscores=[]
                for thres in threslist:
                    tempscore=float(sentscores[sentalgo][sentparameterstr][varparam][thres][ykey])
                    tempscores.append(tempscore) 
                medianscore=np.median(tempscores)
                sentscore=medianscore
            else:
                print "unknown detail {0}".format(detail)
                exit(1)
            ylist.append(sentscore)
    
            paramareascore=estimatearea(xlist,ylist)
            if paramareascore>=bestareascore:
               bestareascore=paramaucscore
               bestxlist=list(xlist)
               bestylist=list(ylist)
               datalabel=plotlabelcreate(sentalgo,sentparameterstr)
               
                
            for varparam in sorted(sentscores[sentalgo][sentparameterstr].keys()):
                if sentalgo in ["lsel1","absel1"]:
                   threslist=sorted(sentscores[sentalgo][sentparameterstr][varparam].keys(),reverse=True)
                elif sentalgo in ["multitree","netinf"]:
                   threslist=sorted(sentscores[sentalgo][sentparameterstr][varparam].keys())
                elif sentalgo in ["netrate"]:
                   pass
                else:
                   print "this algo {0} is unknown".format(sentalgo)
                   exit(1)
                #if firstpart in ["f1","f2","f1overn","f1div100","f1div10","f1div2","f10","f20"]:
                   if detail=="mean":
                      avgscore=0.0   
                      for thres in threslist:
                         avgscore+=float(sentscores[sentalgo][sentparameterstr][varparam][thres][ykey])
                      avgscore/=float(len(threslist))
                      sentscore=avgscore
                   elif detail=="error":   
                      tempscores=[]
                      for thres in threslist:
                         tempscore=float(sentscores[sentalgo][sentparameterstr][varparam][thres][ykey])
                         tempscores.append(tempscore)
                      meanscore=np.mean(tempscores)
                      stddevscore=np.std(tempscores)
                      sentscore=(meanscore,stddevscore)
                   elif detail=="median":
                      tempscores=[]
                      for thres in threslist:
                         tempscore=float(sentscores[sentalgo][sentparameterstr][varparam][thres][ykey])
                         tempscores.append(tempscore) 
                      medianscore=np.median(tempscores)
                      sentscore=medianscore
                   else:
                      print "unknown detail {0}".format(detail)
                      exit(1)
                ylist.append(sentscore)
    
            paramareascore=estimatearea(xlist,ylist)
            if paramareascore>=bestareascore:
               bestareascore=paramaucscore
               bestxlist=list(xlist)
               bestylist=list(ylist)
               datalabel=plotlabelcreate(sentalgo,sentparameterstr)   
            
        assert datalabel!="-1"
        plotxlist[sentalgo]=list(bestxlist)
        plotylist[sentalgo]=list(bestylist)
        labelinfo[sentalgo]=datalabel

    if len(labelinfo.keys())==0:
       print "No data for plot algo:{0} for title:{1} !!".format(runalgo,plottitle)
       return

    if False:
     xmax=-1.0
     xmin=12000.0
     seenxs=set()
     for sentalgo in plotxlist.keys():
        for elem in plotxlist[sentalgo]:
           if elem<xmin:
              xmin=elem
           if elem>xmax:
              xmax=elem   
           seenxs.add(elem)
           
     #plotting part
     plt.clf()
     #plt.set_xticks(seenxs)
     plt.xlim([-0.5+xmin,xmax+0.5])
    
    #plotting part
    plt.clf()
    #plt.xlim([0,0.9])
    plt.ylim([0,1.0]) 
    plt.xlabel(myxlabel)
    plt.ylabel(myylabel)
    plt.title(plottitle)
    for index in range(0,len(plotxlist.keys())):
        datainfo=plotxlist.keys()[index] #datainfo=runalgo
        xlist=plotxlist[datainfo]
        ylist=plotylist[datainfo]
        #plt.plot(xlist,ylist,marker=plotmarkers[index],linestyle='--',color=plotcolors[index],label=labelinfo[index])
        plt.plot(xlist,ylist,marker=plotmarkers[index],linestyle='None',color=plotcolors[index],label=labelinfo[datainfo])
    plt.legend(loc=1)   
    plt.savefig(plotfilepath)


def plotheatmap(datahash,xtitle,ytitle,outfilepath):
    xlabels=sorted(datahash.keys())
    ylabels=sorted(datahash[xlabels[0]].keys())
    A=np.zeros((len(xlabels),len(ylabels)),dtype=np.float)
    for index1 in range(0,len(xlabels)):
        xlabel=xlabels[index1]
        for index2 in range(0,len(ylabels)):
            ylabel=ylabels[index2]
            A[index1,index2]=datahash[xlabel][ylabel]
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    # use dir(matplotlib.cm) to get a list of the installed colormaps
    # the "_r" means "reversed" and accounts for why zero values are plotted as white
    cmap = cm.get_cmap('gray_r', 10)
    im=ax1.imshow(A, interpolation="none", cmap=cmap, aspect='auto')
    fig.colorbar(im)
    plt.xlabel(xtitle)
    plt.ylabel(ytitle)
    plotxlabels=[]
    xlocs=range(0,len(xlabels))
    ylocs=range(0,len(ylabels))
    for elem in xlabels:
        plotxlabels.append(elem)
    plotylabels=[]
    for elem in ylabels:
        plotylabels.append(elem)
    plt.xticks(ylocs,plotylabels)
    plt.yticks(xlocs,plotxlabels)    
    #ax1.set_xticklabels(ylocs,plotylabels)
    #ax1.set_yticklabels(xlocs,plotxlabels)
    plt.savefig(outfilepath)



def varparam_fixedtraceplotgenerate(sentscores,plottype,plottitle,plotfilepath,xlabel,algoparaminfo):
    ykey=plottype
    if plottype in ["sen","fpr","recall","precision","acc","f1","f2","mcc","spec","f1overn","f1div100","f1div10","f1div2","f10","f20","roc","pr"]:
       maxormin="max" 
    else:
       print "this plottype {0} is unknown for max or min !!".format(plottype)
       exit(1) 
    myxlabel=xlabel
    myylabel=plottype.capitalize()

    if plottype in ["roc","pr"]:
       if plottype in ["roc"]:
          subykey="sen"
          subxkey="fpr"
       elif plottype in ["pr"]:
          subykey="precision"
          subxkey="recall" 
       sentscores=sentscores2paramwise(sentscores)
       plotxlist={}
       plotylist={}
       labelinfo={}
       for runalgo in sentscores.keys():
           print runalgo
           maxvalue=-1.0
           maxlabel="-1.0"
           bestxlist=[]
           bestylist=[]
           for sentparamstr in sentscores[runalgo].keys():
               print "here"
               overtimexlist=[]
               overtimeylist=[]
               for fraction in sorted(sentscores[runalgo][sentparamstr].keys()):
                   overtimexlist.append(fraction)
                   thresylist=[] #either sen or precision
                   thresxlist=[] #either fpr or recall
                   if runalgo in ["lsel1","absel1"]:
                      threslist=sorted(sentscores[runalgo][sentparamstr][fraction].keys(),reverse=True)
                   elif runalgo in ["multitree","netinf"]:
                      threslist=sorted(sentscores[runalgo][sentparamstr][fraction].keys())
                   elif runalgo in ["netrate"]:
                      threslist=["1"]
                   else:
                      print "this algo {0} is unknown".format(runalgo)
                      exit(1)
                   for thres in threslist:
                       if runalgo in ["netinf","multitree"]:
                          if thres not in requirededgethres:
                             continue 
                       elif runalgo in ["lsel1","absel1"]:
                          if thres not in requiredlambdathres:
                             continue 
                       elif runalgo in ["netrate"]:
                          if thres not in "1":
                             continue 
                       else:
                          print "sentalgo is unknown{0}".format(runalgo)
                          exit(1)
                       thresylist.append(sentscores[runalgo][sentparamstr][fraction][thres][subykey])
                       thresxlist.append(sentscores[runalgo][sentparamstr][fraction][thres][subxkey])
                   tempscore=estimateaucscore(thresxlist,thresylist,plottype)
                   overtimeylist.append(tempscore)
               areavalue=estimatearea(overtimexlist,overtimeylist)
               if areavalue>=maxvalue:
                  maxvalue=areavalue
                  maxlabel=plotlabelcreate(runalgo,sentparamstr)
                  bestxlist=list(overtimexlist)
                  bestylist=list(overtimeylist)
           assert maxvalue!=-1.0
           plotxlist[runalgo]=bestxlist
           plotylist[runalgo]=bestylist
           labelinfo[runalgo]=maxlabel

       #plotting part
       parts=plotfilepath.split("/")
       partialplotfolder="/".join(parts[0:-1])
       partialplotfilepath="{0}/area_{1}_{2}".format(partialplotfolder,plottype,parts[-1])
       plt.clf()
       plt.xlim([0,max(requiredfractions)+0.05])
       plt.ylim([0,1.0]) 
       plt.xlabel(myxlabel)
       plt.ylabel(myylabel)
       plt.title(plottitle)
       for index in range(0,len(plotxlist.keys())):
           datainfo=plotxlist.keys()[index] 
           xlist=plotxlist[datainfo]
           ylist=plotylist[datainfo]
           plt.plot(xlist,ylist,marker=plotmarkers[index],linestyle='None',color=plotcolors[index],label=labelinfo[datainfo])
       plt.legend(loc=4)   
       plt.savefig(partialplotfilepath)

       #output scores
       parts=plotfilepath.split("/")
       scoreoutfilepath="/".join(parts[0:-1])+"/scoreout_area_{0}".format(parts[-1].replace(".png",".txt"))
       file=open(scoreoutfilepath,"w")
       for datainfo in plotxlist.keys():
           xlist=plotxlist[datainfo]
           xliststr=str(xlist[0])
           for elem in xlist[1:]:
               xliststr+=",{0}".format(elem)
           ylist=plotylist[datainfo]
           yliststr=str(ylist[0])
           for elem in ylist[1:]:
               yliststr+=",{0}".format(elem)
           file.write("{0}\t{1}\t{2}\n".format(datainfo,xliststr,yliststr))
       file.close()
       
       return
   
    if algoparaminfo!="-1":
       plotxlist={}
       plotylist={}
       labelinfo={} 
       for runalgo in algoparaminfo.keys():
           for algoparam in algoparaminfo[runalgo]:
               foundparamstr="--"
               if not sentscores.has_key(runalgo):
                  continue
               #print "starting {0} {1}".format(runalgo,algoparam)
               for sentparameterstr in sentscores[runalgo].keys():
                   parameterstr,extend=sentparameterstr.split("+")
                   parts=parameterstr.replace("result_","").split("_")
                   if runalgo in ["netinf","multitree"]:
                      tempparam=float(parts[-1].replace("categorical?",""))
                   elif runalgo in ["netrate"]:
                      tempparam="1" #there won't be param for netrate
                   elif runalgo in ["lsel1","absel1"]: #our algos
                      tempparam=float(parts[-1].split("?")[0])
                   elif runalgo in ["lsel1fused","absel1fused"]: #our algos dynamic
                      tempparam=float(parts[-1].split("?")[0])
                   else:
                      print "this algo is unknown {0}".format(sentalgo)
                      exit(1)
                   #result_uniform_undirected_0.6-randomnoise_epsilon-0.9_0.001??2.0?epsilon0.01?flowdifference_3
                   #print "{0} -> {1}".format(runalgo,tempparam)
                   if tempparam==algoparam:
                      foundparamstr=sentparameterstr
                      break
               #print "lastinfo:"    
               #print runalgo
               #print algoparam
               assert foundparamstr!="--"
               xlist=[]
               ylist=[]
               for fraction in sorted(sentscores[runalgo][foundparamstr].keys()):
                   if fraction not in requiredfractions:
                      continue 
                   ylist.append(sentscores[runalgo][sentparameterstr][fraction][plottype])
                   xlist.append(float(fraction))
               datalabel=plotlabelcreate(runalgo,sentparameterstr,"algoparam")
               plotxlist[datalabel]=xlist
               plotylist[datalabel]=ylist
               labelinfo[datalabel]=datalabel

       #plotting part
       parts=plotfilepath.split("/")
       partialplotfolder="/".join(parts[0:-1])
       partialplotfilepath="{0}/partial_{1}".format(partialplotfolder,parts[-1])
       plt.clf()
       plt.xlim([0,max(requiredfractions)+0.05])
       plt.ylim([0,1.0]) 
       plt.xlabel(myxlabel)
       plt.ylabel(myylabel)
       plt.title(plottitle)
       for index in range(0,len(plotxlist.keys())):
           datainfo=plotxlist.keys()[index] 
           xlist=plotxlist[datainfo]
           ylist=plotylist[datainfo]
           plt.plot(xlist,ylist,marker=plotmarkers[index],linestyle='None',color=plotcolors[index],label=labelinfo[datainfo])
       plt.legend(loc=4)   
       plt.savefig(partialplotfilepath)

       #output scores
       parts=plotfilepath.split("/")
       scoreoutfilepath="/".join(parts[0:-1])+"/scoreout_partial_{0}".format(parts[-1].replace(".png",".txt"))
       file=open(scoreoutfilepath,"w")
       for datainfo in plotxlist.keys():
           xlist=plotxlist[datainfo]
           xliststr=str(xlist[0])
           for elem in xlist[1:]:
               xliststr+=",{0}".format(elem)
           ylist=plotylist[datainfo]
           yliststr=str(ylist[0])
           for elem in ylist[1:]:
               yliststr+=",{0}".format(elem)
           file.write("{0}\t{1}\t{2}\n".format(datainfo,xliststr,yliststr))
       file.close()
       
       return 
    
    #preprocessing part
    bestplotxlist={}
    bestplotylist={}
    worstplotxlist={}
    worstplotylist={}
    bestlabelinfo={}
    worstlabelinfo={}
    plotmeanlist={}
    plotstddevlist={}
    for runalgo in sentscores.keys():
        if maxormin=="max":
           bestscore=-1.0
           worstscore=100.0
        elif maxormin=="min":
           bestscore=100.0
           worstscore=-1.0
        bestylist=[]
        bestxlist=[]
        allylist=[] #all y list scores
        worstxlist=[]
        worstylist=[]
        #bestsen=[]
        #bestfpr=[]
        #bestaucscore=-1.0
        bestdatalabel="-1"
        worstdatalabel="-1"
        for sentparameterstr in sentscores[runalgo].keys():
            parameterstr,extend=sentparameterstr.split("+")
            #tempsenlist=[]
            #tempfprlist=[]
            #for fraction in sorted(sentscores[runalgo][sentparameterstr].keys()):
            #    if fraction not in requiredfractions:
            #       continue 
            #    tempsenlist.append(sentscores[runalgo][sentparameterstr][fraction]["sen"])
            #    tempfprlist.append(sentscores[runalgo][sentparameterstr][fraction]["fpr"])
            #paramaucscore=estimateaucscore(tempfprlist,tempsenlist,"roc")
            #if paramaucscore>=bestaucscore:
            #   bestaucscore=paramaucscore
            
            ylist=[]
            xlist=[]
            for fraction in sorted(sentscores[runalgo][sentparameterstr].keys()):
                if fraction not in requiredfractions:
                   continue 
                ylist.append(sentscores[runalgo][sentparameterstr][fraction][plottype])
                xlist.append(float(fraction))
            allylist.append(ylist)    
            paramareascore=estimatearea(xlist,ylist)
            #if plottype=="precision":
            #   print "area {0}".format(paramareascore)
            #   print "param: {0}".format(sentparameterstr)
            if maxormin=="max":
               if paramareascore>=bestscore:
                  bestscore=paramareascore
                  bestylist=list(ylist)
                  bestxlist=list(xlist)
                  bestdatalabel=plotlabelcreate(runalgo,sentparameterstr)
               if paramareascore<=worstscore:
                  worstscore=paramareascore
                  worstylist=list(ylist)
                  worstxlist=list(xlist)
                  worstdatalabel=plotlabelcreate(runalgo,sentparameterstr)   
            elif maxormin=="min":
               if paramareascore<=bestscore:
                  bestscore=paramareascore
                  bestylist=list(ylist)
                  bestxlist=list(xlist)
                  bestdatalabel=plotlabelcreate(runalgo,sentparameterstr)
               if paramareascore>=worstscore:
                  worstscore=paramareascore
                  worstylist=list(ylist)
                  worstxlist=list(xlist)
                  worstdatalabel=plotlabelcreate(runalgo,sentparameterstr)     

        meanlist=[]
        stddevlist=[]
        for fracindex in range(0,len(requiredfractions)):
            locallist=[]
            for elem in allylist:
                locallist.append(elem[fracindex])
            mean=np.mean(locallist)
            stddev=np.std(locallist)
            meanlist.append(mean)
            stddevlist.append(stddev)
            
        assert bestdatalabel!="-1"
        assert worstdatalabel!="-1"
        bestplotxlist[runalgo]=list(bestxlist)
        bestplotylist[runalgo]=list(bestylist)
        worstplotxlist[runalgo]=list(worstxlist)
        worstplotylist[runalgo]=list(worstylist)
        plotmeanlist[runalgo]=list(meanlist)
        plotstddevlist[runalgo]=list(stddevlist)
        bestlabelinfo[runalgo]=bestdatalabel
        worstlabelinfo[runalgo]=worstdatalabel

        #print "runalgo auc score {0}: {1}".format(runalgo,bestaucscore)
             
    if len(bestlabelinfo.keys())==0 or len(worstlabelinfo.keys())==0:
       print "No data for plot algo:{0} for title:{1} !!".format(runalgo,plottitle)
       return

    #best performing plotting part
    parts=plotfilepath.split("/")
    bestplotfolder="/".join(parts[0:-1])
    bestplotfilepath="{0}/best_{1}".format(bestplotfolder,parts[-1])
    plt.clf()
    plt.xlim([0,max(requiredfractions)+0.05])
    plt.ylim([0,1.0]) 
    plt.xlabel(myxlabel)
    plt.ylabel(myylabel)
    plt.title(plottitle)
    for index in range(0,len(bestplotxlist.keys())):
        datainfo=bestplotxlist.keys()[index] 
        xlist=bestplotxlist[datainfo]
        ylist=bestplotylist[datainfo]
        #plt.plot(xlist,ylist,marker=plotmarkers[index],linestyle='--',color=plotcolors[index],label=labelinfo[index])
        plt.plot(xlist,ylist,marker=plotmarkers[index],linestyle='None',color=plotcolors[index],label=bestlabelinfo[datainfo])
    plt.legend(loc=4)   
    plt.savefig(bestplotfilepath)
    #output scores best performing
    parts=plotfilepath.split("/")
    scoreoutfilepath="/".join(parts[0:-1])+"/scoreout_best_{0}".format(parts[-1].replace(".png",".txt"))
    file=open(scoreoutfilepath,"w")
    for datainfo in bestplotxlist.keys():
        xlist=bestplotxlist[datainfo]
        xliststr=str(xlist[0])
        for elem in xlist[1:]:
            xliststr+=",{0}".format(elem)
        ylist=bestplotylist[datainfo]
        yliststr=str(ylist[0])
        for elem in ylist[1:]:
            yliststr+=",{0}".format(elem)
        file.write("{0}\t{1}\t{2}\n".format(datainfo,xliststr,yliststr))
    file.close()
       
    #worst performing plotting part
    parts=plotfilepath.split("/")
    worstplotfolder="/".join(parts[0:-1])
    worstplotfilepath="{0}/worst_{1}".format(worstplotfolder,parts[-1])
    plt.clf()
    plt.xlim([0,max(requiredfractions)+0.05])
    plt.ylim([0,1.0]) 
    plt.xlabel(myxlabel)
    plt.ylabel(myylabel)
    plt.title(plottitle)
    for index in range(0,len(worstplotxlist.keys())):
        datainfo=worstplotxlist.keys()[index] 
        xlist=worstplotxlist[datainfo]
        ylist=worstplotylist[datainfo]
        #plt.plot(xlist,ylist,marker=plotmarkers[index],linestyle='--',color=plotcolors[index],label=labelinfo[index])
        plt.plot(xlist,ylist,marker=plotmarkers[index],linestyle='None',color=plotcolors[index],label=worstlabelinfo[datainfo])
    plt.legend(loc=4)   
    plt.savefig(worstplotfilepath)
    #output scores worst performing
    parts=plotfilepath.split("/")
    scoreoutfilepath="/".join(parts[0:-1])+"/scoreout_worst_{0}".format(parts[-1].replace(".png",".txt"))
    file=open(scoreoutfilepath,"w")
    for datainfo in worstplotxlist.keys():
        xlist=worstplotxlist[datainfo]
        xliststr=str(xlist[0])
        for elem in xlist[1:]:
            xliststr+=",{0}".format(elem)
        ylist=worstplotylist[datainfo]
        yliststr=str(ylist[0])
        for elem in ylist[1:]:
            yliststr+=",{0}".format(elem)
        file.write("{0}\t{1}\t{2}\n".format(datainfo,xliststr,yliststr))
    file.close()
    
    #average performing plotting part
    parts=plotfilepath.split("/")
    meanplotfolder="/".join(parts[0:-1])
    meanplotfilepath="{0}/mean_{1}".format(meanplotfolder,parts[-1])
    plt.clf()
    plt.xlim([0,max(requiredfractions)+0.05])
    plt.ylim([0,1.0]) 
    plt.xlabel(myxlabel)
    plt.ylabel(myylabel)
    plt.title(plottitle)
    for index in range(0,len(bestplotxlist.keys())):
        datainfo=bestplotxlist.keys()[index] 
        xlist=bestplotxlist[datainfo]
        meanylist=plotmeanlist[datainfo]
        #plt.plot(xlist,ylist,marker=plotmarkers[index],linestyle='--',color=plotcolors[index],label=labelinfo[index])
        plt.plot(xlist,meanylist,marker=plotmarkers[index],linestyle='None',color=plotcolors[index],label=bestlabelinfo[datainfo])
    plt.legend(loc=4)   
    plt.savefig(meanplotfilepath)
    #output scores mean performing
    parts=plotfilepath.split("/")
    scoreoutfilepath="/".join(parts[0:-1])+"/scoreout_mean_{0}".format(parts[-1].replace(".png",".txt"))
    file=open(scoreoutfilepath,"w")
    for datainfo in bestplotxlist.keys():
        xlist=bestplotxlist[datainfo]
        xliststr=str(xlist[0])
        for elem in xlist[1:]:
            xliststr+=",{0}".format(elem)
        ylist=plotmeanlist[datainfo]
        yliststr=str(ylist[0])
        for elem in ylist[1:]:
            yliststr+=",{0}".format(elem)
        file.write("{0}\t{1}\t{2}\n".format(datainfo,xliststr,yliststr))
    file.close()
    
    #error plotting part
    parts=plotfilepath.split("/")
    errorplotfolder="/".join(parts[0:-1])
    errorplotfilepath="{0}/error_{1}".format(errorplotfolder,parts[-1])
    plt.clf()
    plt.xlim([0,max(requiredfractions)+0.05])
    plt.ylim([0,1.0]) 
    plt.xlabel(myxlabel)
    plt.ylabel(myylabel)
    plt.title(plottitle)
    for index in range(0,len(bestplotxlist.keys())):
        datainfo=bestplotxlist.keys()[index] #datainfo=runalgo
        xlist=bestplotxlist[datainfo]
        meanylist=plotmeanlist[datainfo]
        stddevlist=plotstddevlist[datainfo]
        plt.errorbar(xlist, meanylist,yerr=stddevlist, fmt=plotmarkers[index], color=plotcolors[index],ecolor=plotcolors[index],label=bestlabelinfo[datainfo])
        #plt.errorbar(xlist, meanylist,fmt='ro', ecolor='g',yerr=stddevlist, marker=plotmarkers[index],linestyle='--',color=plotcolors[index],label=labelinfo[datainfo]) 
    plt.legend(loc=4)   
    plt.savefig(errorplotfilepath)
    #output scores errorbar performing
    parts=plotfilepath.split("/")
    scoreoutfilepath="/".join(parts[0:-1])+"/scoreout_error_{0}".format(parts[-1].replace(".png",".txt"))
    file=open(scoreoutfilepath,"w")
    for datainfo in bestplotxlist.keys():
        xlist=bestplotxlist[datainfo]
        xliststr=str(xlist[0])
        for elem in xlist[1:]:
            xliststr+=",{0}".format(elem)
        ylist=plotmeanlist[datainfo]
        yliststr=str(ylist[0])
        for elem in ylist[1:]:
            yliststr+=",{0}".format(elem)
        zlist=plotstddevlist[datainfo]
        zliststr=str(zlist[0])
        for elem in zlist[1:]:
            zliststr+=",{0}".format(elem)    
        file.write("{0}\t{1}\t{2}\t{3}\n".format(datainfo,xliststr,yliststr,zliststr))
    file.close()
    
    
    
    ykey=plottype
    assert plottype in ["f1","f1overn","f1div100","f1div10","f1div2","f20","f50","f100"]
    maxormin="max" 
    myxlabel=xlabel
    myylabel=plottype.capitalize()
    
    #preprocessing part
    plotxlist={}
    plotylist={}
    plotmeanlist={}
    plotstddevlist={}
    labelinfo={}
    for sentalgo in sentscores.keys():
        datalabel="-1"
        bestylist=[]
        bestxlist=[]
        bestylist=[]
        bestxlist=[]
        bestylist=[]
        bestxlist=[]
        if maxormin=="max": #can be either area or average f
           bestareascore=-1.0
        elif maxormin=="min":
           bestareascore=100.0
        for sentparameterstr in sentscores[sentalgo].keys():
            ylist=[]
            xlist=sorted(sentscores[sentalgo][sentparameterstr].keys())
            for varparam in xlist:
                ylist.append(sentscores[sentalgo][sentparameterstr][varparam])
            #max score

            #mean score    
            if detail=="mean":
                avgscore=0.0   
                for thres in threslist:
                    avgscore+=float(sentscores[sentalgo][sentparameterstr][varparam][thres][ykey])
                avgscore/=float(len(threslist))
                sentscore=avgscore
            elif detail=="error":   
                tempscores=[]
                for thres in threslist:
                    tempscore=float(sentscores[sentalgo][sentparameterstr][varparam][thres][ykey])
                    tempscores.append(tempscore)
                meanscore=np.mean(tempscores)
                stddevscore=np.std(tempscores)
                sentscore=(meanscore,stddevscore)
            elif detail=="median":
                tempscores=[]
                for thres in threslist:
                    tempscore=float(sentscores[sentalgo][sentparameterstr][varparam][thres][ykey])
                    tempscores.append(tempscore) 
                medianscore=np.median(tempscores)
                sentscore=medianscore
            else:
                print "unknown detail {0}".format(detail)
                exit(1)
            ylist.append(sentscore)
    
            paramareascore=estimatearea(xlist,ylist)
            if paramareascore>=bestareascore:
               bestareascore=paramaucscore
               bestxlist=list(xlist)
               bestylist=list(ylist)
               datalabel=plotlabelcreate(sentalgo,sentparameterstr)
               
                
            for varparam in sorted(sentscores[sentalgo][sentparameterstr].keys()):
                if sentalgo in ["lsel1","absel1"]:
                   threslist=sorted(sentscores[sentalgo][sentparameterstr][varparam].keys(),reverse=True)
                elif sentalgo in ["multitree","netinf"]:
                   threslist=sorted(sentscores[sentalgo][sentparameterstr][varparam].keys())
                elif sentalgo in ["netrate"]:
                   pass
                else:
                   print "this algo {0} is unknown".format(sentalgo)
                   exit(1)
                #if firstpart in ["f1","f2","f1overn","f1div100","f1div10","f1div2","f10","f20"]:
                   if detail=="mean":
                      avgscore=0.0   
                      for thres in threslist:
                         avgscore+=float(sentscores[sentalgo][sentparameterstr][varparam][thres][ykey])
                      avgscore/=float(len(threslist))
                      sentscore=avgscore
                   elif detail=="error":   
                      tempscores=[]
                      for thres in threslist:
                         tempscore=float(sentscores[sentalgo][sentparameterstr][varparam][thres][ykey])
                         tempscores.append(tempscore)
                      meanscore=np.mean(tempscores)
                      stddevscore=np.std(tempscores)
                      sentscore=(meanscore,stddevscore)
                   elif detail=="median":
                      tempscores=[]
                      for thres in threslist:
                         tempscore=float(sentscores[sentalgo][sentparameterstr][varparam][thres][ykey])
                         tempscores.append(tempscore) 
                      medianscore=np.median(tempscores)
                      sentscore=medianscore
                   else:
                      print "unknown detail {0}".format(detail)
                      exit(1)
                ylist.append(sentscore)
    
            paramareascore=estimatearea(xlist,ylist)
            if paramareascore>=bestareascore:
               bestareascore=paramaucscore
               bestxlist=list(xlist)
               bestylist=list(ylist)
               datalabel=plotlabelcreate(sentalgo,sentparameterstr)   
            
        assert datalabel!="-1"
        plotxlist[sentalgo]=list(bestxlist)
        plotylist[sentalgo]=list(bestylist)
        labelinfo[sentalgo]=datalabel

    if len(labelinfo.keys())==0:
       print "No data for plot algo:{0} for title:{1} !!".format(runalgo,plottitle)
       return

    if False:
     xmax=-1.0
     xmin=12000.0
     seenxs=set()
     for sentalgo in plotxlist.keys():
        for elem in plotxlist[sentalgo]:
           if elem<xmin:
              xmin=elem
           if elem>xmax:
              xmax=elem   
           seenxs.add(elem)
           
     #plotting part
     plt.clf()
     #plt.set_xticks(seenxs)
     plt.xlim([-0.5+xmin,xmax+0.5])
    
    #plotting part
    plt.clf()
    #plt.xlim([0,0.9])
    plt.ylim([0,1.0]) 
    plt.xlabel(myxlabel)
    plt.ylabel(myylabel)
    plt.title(plottitle)
    for index in range(0,len(plotxlist.keys())):
        datainfo=plotxlist.keys()[index] #datainfo=runalgo
        xlist=plotxlist[datainfo]
        ylist=plotylist[datainfo]
        #plt.plot(xlist,ylist,marker=plotmarkers[index],linestyle='--',color=plotcolors[index],label=labelinfo[index])
        plt.plot(xlist,ylist,marker=plotmarkers[index],linestyle='None',color=plotcolors[index],label=labelinfo[datainfo])
    plt.legend(loc=1)   
    plt.savefig(plotfilepath)




def returnbunchplots(key,sentparams):
    requiredplots=[]
    #plotinfo5=[key,"acc",sentparams]
    #requiredplots.append(plotinfo5)
    plotinfo6=[key,"f1",sentparams]
    requiredplots.append(plotinfo6)
    plotinfo6=[key,"f1div10",sentparams]
    requiredplots.append(plotinfo6) 
    plotinfo6=[key,"f1div2",sentparams]
    requiredplots.append(plotinfo6)
    #plotinfo6=[key,"f10",sentparams]
    #requiredplots.append(plotinfo6)
    #plotinfo6=[key,"f20",sentparams]
    #requiredplots.append(plotinfo6)
    ##plotinfo7=[key,"f1overn",sentparams]
    ##requiredplots.append(plotinfo7)
    #plotinfo8=[key,"spec",sentparams]
    #requiredplots.append(plotinfo8)
    #plotinfo9=[key,"sen",sentparams]
    #requiredplots.append(plotinfo9)
    #plotinfo10=[key,"precision",sentparams]
    #requiredplots.append(plotinfo10)
    #plotinfo11=[key,"recall",sentparams]
    #requiredplots.append(plotinfo11)
    return requiredplots


#only infertype will have a main impact on the plots generated!!!
def requiredplotreturner(realdata,graphevolution,spreadmodel,infertype):
    realstaticgraphnames=["enron-email-clustered-261-2444-7-2","stanford-contact-static-207-2962-300-1"]
    realdynamicgraphnames=["stanford-contact-dynamic-t24s2"]
    syngraphwisegraphnames=["datan250deg10/rds/13.sparse6","datan250deg10/lpa/13.sparse6","datan250deg10/undirected_ff/13.sparse6","datan250deg10/dmc/13.sparse6"]
    syndynamicgraphnames=["datan32s2t20/fusedlasso/intervaldistuniform_intervaldistparam0.02_startgraphalgords_startgraphalgoparam0.1/20","datan32s2t20/fusedlasso/intervaldistuniform_intervaldistparam0.1_startgraphalgords_startgraphalgoparam0.1/20","datan64s2t20/fusedlasso/intervaldistuniform_intervaldistparam0.02_startgraphalgords_startgraphalgoparam0.1/20","datan64s2t20/fusedlasso/intervaldistuniform_intervaldistparam0.1_startgraphalgords_startgraphalgoparam0.1/20"]
    requiredplots=[]
    if infertype=="edge":
       if (realdata=="real" and graphevolution=="static") or (realdata=="synthetic" and graphevolution=="graphwise"):
          if (realdata=="real" and graphevolution=="static"):
             graphnames=realstaticgraphnames
          elif(realdata=="synthetic" and graphevolution=="graphwise"):
             graphnames=syngraphwisegraphnames
          runalgos=["multitree","netinf","netrate"] #"lsel1mul","absel1mul" "lsel1sum","absel1sum",
          if spreadmodel in ["si"]:
             for graphname in graphnames:
                 #overtime plot
                 if graphname in ["stanford-contact-static-207-2962-300-1"]:   
                     distinfo="spreadprob_0.5-s2i_weibull_3.0_2.5"
                 elif graphname in ["enron-email-clustered-261-2444-7-2"]:
                     distinfo="spreadprob_0.5-s2i_expo_2.0"
                 else:
                     distinfo="spreadprob_0.5-s2i_expo_2.0"
                      
                 sampleinterval=1
                 noiselevel=1.0
                 algoparaminfo={}
                 algoparaminfo["lsel1"]=[0.01,1.0]
                 algoparaminfo["absel1"]=[0.01,0.05]
                 algoparaminfo["multitree"]=[1.2]
                 algoparaminfo["netinf"]=[1.2]
                 algoparaminfo["netrate"]=["1"]
                 sentparams=["si",distinfo,sampleinterval,runalgos,graphname,realdata,graphevolution,noiselevel,algoparaminfo]
                 bunchplot=returnbunchplots("overtime",sentparams)
                 requiredplots.extend(bunchplot)

                 sampleinterval=1
                 noiselevel=1.0
                 sentparams=["si",distinfo,sampleinterval,runalgos,graphname,realdata,graphevolution,noiselevel,"-1"]
                 bunchplot=returnbunchplots("overtime",sentparams)
                 requiredplots.extend(bunchplot)

                 if graphname in ["stanford-contact-static-207-2962-300-1"]:   
                    distinfo="spreadprob_0.1-s2i_weibull_3.0_2.5"
                 elif graphname in ["enron-email-clustered-261-2444-7-2"]:
                    distinfo="spreadprob_0.1-s2i_expo_2.0"
                 else:
                     distinfo="spreadprob_0.1-s2i_expo_2.0"
                     
                 sampleinterval=1
                 noiselevel=1.0
                 algoparaminfo={}
                 algoparaminfo["lsel1"]=[0.01,1.0]
                 algoparaminfo["absel1"]=[0.01,0.05]
                 algoparaminfo["multitree"]=[1.2]
                 algoparaminfo["netinf"]=[1.2]
                 algoparaminfo["netrate"]=["1"]
                 sentparams=["si",distinfo,sampleinterval,runalgos,graphname,realdata,graphevolution,noiselevel,algoparaminfo]
                 bunchplot=returnbunchplots("overtime",sentparams)
                 requiredplots.extend(bunchplot)

                 sampleinterval=1
                 noiselevel=1.0
                 sentparams=["si",distinfo,sampleinterval,runalgos,graphname,realdata,graphevolution,noiselevel,"-1"]
                 bunchplot=returnbunchplots("overtime",sentparams)
                 requiredplots.extend(bunchplot)

                 #parameter auc curve
                 if graphname in ["stanford-contact-static-207-2962-300-1"]:
                    distinfo="spreadprob_0.5-s2i_weibull_3.0_2.5"
                 else:    
                    distinfo="spreadprob_0.5-s2i_expo_2.0" 
                 sampleinterval=1
                 noiselevel=1.0
                 fixedfraction=0.15
                 sentparams=["si",distinfo,sampleinterval,runalgos,graphname,realdata,graphevolution,noiselevel,fixedfraction]
                 myplotinfo=["parameterauccurve","roc",sentparams]
                 requiredplots.append(myplotinfo)
                 #myplotinfo=["parameterauccurve","bedroc",sentparams]
                 #requiredplots.append(myplotinfo)
                 myplotinfo=["parameterauccurve","pr",sentparams]
                 requiredplots.append(myplotinfo)
          
                 sampleinterval=1
                 noiselevel=1.0
                 fixedfraction=0.1
                 sentparams=["si",distinfo,sampleinterval,runalgos,graphname,realdata,graphevolution,noiselevel,fixedfraction]
                 myplotinfo=["parameterauccurve","roc",sentparams]
                 requiredplots.append(myplotinfo)
                 #myplotinfo=["parameterauccurve","bedroc",sentparams]
                 #requiredplots.append(myplotinfo)
                 myplotinfo=["parameterauccurve","pr",sentparams]
                 requiredplots.append(myplotinfo)
          
                 sampleinterval=1
                 noiselevel=1.0
                 fixedfraction=0.07
                 sentparams=["si",distinfo,sampleinterval,runalgos,graphname,realdata,graphevolution,noiselevel,fixedfraction]
                 myplotinfo=["parameterauccurve","roc",sentparams]
                 requiredplots.append(myplotinfo)
                 #myplotinfo=["parameterauccurve","bedroc",sentparams]
                 #requiredplots.append(myplotinfo)
                 myplotinfo=["parameterauccurve","pr",sentparams]
                 requiredplots.append(myplotinfo)
          
                 sampleinterval=1
                 noiselevel=1.0
                 fixedfraction=0.05
                 sentparams=["si",distinfo,sampleinterval,runalgos,graphname,realdata,graphevolution,noiselevel,fixedfraction]
                 myplotinfo=["parameterauccurve","roc",sentparams]
                 requiredplots.append(myplotinfo)
                 #myplotinfo=["parameterauccurve","bedroc",sentparams]
                 #requiredplots.append(myplotinfo)
                 myplotinfo=["parameterauccurve","pr",sentparams]
                 requiredplots.append(myplotinfo)
          
                 sampleinterval=1
                 noiselevel=1.0
                 fixedfraction=0.2
                 sentparams=["si",distinfo,sampleinterval,runalgos,graphname,realdata,graphevolution,noiselevel,fixedfraction]
                 myplotinfo=["parameterauccurve","roc",sentparams]
                 requiredplots.append(myplotinfo)
                 #myplotinfo=["parameterauccurve","bedroc",sentparams]
                 #requiredplots.append(myplotinfo)
                 myplotinfo=["parameterauccurve","pr",sentparams]
                 requiredplots.append(myplotinfo)

                 #parameter auc curve
                 if graphname in ["stanford-contact-static-207-2962-300-1"]:
                    distinfo="spreadprob_0.1-s2i_weibull_3.0_2.5"
                 else:    
                    distinfo="spreadprob_0.1-s2i_expo_2.0"

                 sampleinterval=1
                 noiselevel=1.0
                 fixedfraction=0.8
                 sentparams=["si",distinfo,sampleinterval,runalgos,graphname,realdata,graphevolution,noiselevel,fixedfraction]
                 myplotinfo=["parameterauccurve","roc",sentparams]
                 requiredplots.append(myplotinfo)
                 #myplotinfo=["parameterauccurve","bedroc",sentparams]
                 #requiredplots.append(myplotinfo)
                 myplotinfo=["parameterauccurve","pr",sentparams]
                 requiredplots.append(myplotinfo)

                 sampleinterval=1
                 noiselevel=1.0
                 fixedfraction=0.6
                 sentparams=["si",distinfo,sampleinterval,runalgos,graphname,realdata,graphevolution,noiselevel,fixedfraction]
                 myplotinfo=["parameterauccurve","roc",sentparams]
                 requiredplots.append(myplotinfo)
                 #myplotinfo=["parameterauccurve","bedroc",sentparams]
                 #requiredplots.append(myplotinfo)
                 myplotinfo=["parameterauccurve","pr",sentparams]
                 requiredplots.append(myplotinfo)

                 sampleinterval=1
                 noiselevel=1.0
                 fixedfraction=0.5
                 sentparams=["si",distinfo,sampleinterval,runalgos,graphname,realdata,graphevolution,noiselevel,fixedfraction]
                 myplotinfo=["parameterauccurve","roc",sentparams]
                 requiredplots.append(myplotinfo)
                 #myplotinfo=["parameterauccurve","bedroc",sentparams]
                 #requiredplots.append(myplotinfo)
                 myplotinfo=["parameterauccurve","pr",sentparams]
                 requiredplots.append(myplotinfo)
                 
                 sampleinterval=1
                 noiselevel=1.0
                 fixedfraction=0.15
                 sentparams=["si",distinfo,sampleinterval,runalgos,graphname,realdata,graphevolution,noiselevel,fixedfraction]
                 myplotinfo=["parameterauccurve","roc",sentparams]
                 requiredplots.append(myplotinfo)
                 #myplotinfo=["parameterauccurve","bedroc",sentparams]
                 #requiredplots.append(myplotinfo)
                 myplotinfo=["parameterauccurve","pr",sentparams]
                 requiredplots.append(myplotinfo)

                 sampleinterval=1
                 noiselevel=1.0
                 fixedfraction=0.1
                 sentparams=["si",distinfo,sampleinterval,runalgos,graphname,realdata,graphevolution,noiselevel,fixedfraction]
                 myplotinfo=["parameterauccurve","roc",sentparams]
                 requiredplots.append(myplotinfo)
                 #myplotinfo=["parameterauccurve","bedroc",sentparams]
                 #requiredplots.append(myplotinfo)
                 myplotinfo=["parameterauccurve","pr",sentparams]
                 requiredplots.append(myplotinfo)
          
                 sampleinterval=1
                 noiselevel=1.0
                 fixedfraction=0.07
                 sentparams=["si",distinfo,sampleinterval,runalgos,graphname,realdata,graphevolution,noiselevel,fixedfraction]
                 myplotinfo=["parameterauccurve","roc",sentparams]
                 requiredplots.append(myplotinfo)
                 #myplotinfo=["parameterauccurve","bedroc",sentparams]
                 #requiredplots.append(myplotinfo)
                 myplotinfo=["parameterauccurve","pr",sentparams]
                 requiredplots.append(myplotinfo)
          
                 sampleinterval=1
                 noiselevel=1.0
                 fixedfraction=0.05
                 sentparams=["si",distinfo,sampleinterval,runalgos,graphname,realdata,graphevolution,noiselevel,fixedfraction]
                 myplotinfo=["parameterauccurve","roc",sentparams]
                 requiredplots.append(myplotinfo)
                 #myplotinfo=["parameterauccurve","bedroc",sentparams]
                 #requiredplots.append(myplotinfo)
                 myplotinfo=["parameterauccurve","pr",sentparams]
                 requiredplots.append(myplotinfo)
          
                 sampleinterval=1
                 noiselevel=1.0
                 fixedfraction=0.2
                 sentparams=["si",distinfo,sampleinterval,runalgos,graphname,realdata,graphevolution,noiselevel,fixedfraction]
                 myplotinfo=["parameterauccurve","roc",sentparams]
                 requiredplots.append(myplotinfo)
                 #myplotinfo=["parameterauccurve","bedroc",sentparams]
                 #requiredplots.append(myplotinfo)
                 myplotinfo=["parameterauccurve","pr",sentparams]
                 requiredplots.append(myplotinfo)
                 
                 #sampling plots,fixed trace count
                 if graphname in ["stanford-contact-static-207-2962-300-1"]:
                    distinfo="spreadprob_0.5-s2i_weibull_3.0_2.5"
                 else:    
                    distinfo="spreadprob_0.5-s2i_expo_2.0"
                 sampleintervals=[1,2,3,4,5,6,8,12]
                 noiselevel=1.0
                 fixedfraction=0.1
                 algoparaminfo={}
                 algoparaminfo["lsel1"]=[0.01,1.0]
                 algoparaminfo["absel1"]=[0.01,0.05]
                 algoparaminfo["multitree"]=[1.2]
                 algoparaminfo["netinf"]=[1.2]
                 algoparaminfo["netrate"]=["1"]
                 sentparams=["si",distinfo,sampleintervals,runalgos,graphname,realdata,graphevolution,noiselevel,fixedfraction,algoparaminfo]
                 bunchplot=returnbunchplots("sampling_fixedtrace",sentparams)
                 requiredplots.extend(bunchplot)

                 noiselevel=1.0
                 fixedfraction=0.1
                 sentparams=["si",distinfo,sampleintervals,runalgos,graphname,realdata,graphevolution,noiselevel,fixedfraction,"-1"]
                 bunchplot=returnbunchplots("sampling_fixedtrace",sentparams)
                 requiredplots.extend(bunchplot)
                 #pr and roc can not be run for specific parameters since the score is estimated over parameters
                 #plotinfo11=["sampling_fixedtrace","pr",sentparams] 
                 #requiredplots.append(plotinfo11)
                 #plotinfo11=["sampling_fixedtrace","roc",sentparams]
                 #requiredplots.append(plotinfo11)

                 noiselevel=1.0
                 fixedfraction=0.15
                 sentparams=["si",distinfo,sampleintervals,runalgos,graphname,realdata,graphevolution,noiselevel,fixedfraction,algoparaminfo]
                 bunchplot=returnbunchplots("sampling_fixedtrace",sentparams)
                 requiredplots.extend(bunchplot)

                 noiselevel=1.0
                 fixedfraction=0.15
                 sentparams=["si",distinfo,sampleintervals,runalgos,graphname,realdata,graphevolution,noiselevel,fixedfraction,"-1"]
                 bunchplot=returnbunchplots("sampling_fixedtrace",sentparams)
                 requiredplots.extend(bunchplot)
                 #pr and roc can not be run for specific parameters since the score is estimated over parameters
                 #plotinfo11=["sampling_fixedtrace","pr",sentparams] 
                 #requiredplots.append(plotinfo11)
                 #plotinfo11=["sampling_fixedtrace","roc",sentparams]
                 #requiredplots.append(plotinfo11)

                 noiselevel=1.0
                 fixedfraction=0.2
                 sentparams=["si",distinfo,sampleintervals,runalgos,graphname,realdata,graphevolution,noiselevel,fixedfraction,algoparaminfo]
                 bunchplot=returnbunchplots("sampling_fixedtrace",sentparams)
                 requiredplots.extend(bunchplot)

                 noiselevel=1.0
                 fixedfraction=0.2
                 sentparams=["si",distinfo,sampleintervals,runalgos,graphname,realdata,graphevolution,noiselevel,fixedfraction,"-1"]
                 bunchplot=returnbunchplots("sampling_fixedtrace",sentparams)
                 requiredplots.extend(bunchplot)
                 #pr and roc can not be run for specific parameters since the score is estimated over parameters
                 #plotinfo11=["sampling_fixedtrace","pr",sentparams] 
                 #requiredplots.append(plotinfo11)
                 #plotinfo11=["sampling_fixedtrace","roc",sentparams]
                 #requiredplots.append(plotinfo11)

                 noiselevel=1.0
                 fixedfraction=0.07
                 sentparams=["si",distinfo,sampleintervals,runalgos,graphname,realdata,graphevolution,noiselevel,fixedfraction,algoparaminfo]
                 bunchplot=returnbunchplots("sampling_fixedtrace",sentparams)
                 requiredplots.extend(bunchplot)

                 noiselevel=1.0
                 fixedfraction=0.07
                 sentparams=["si",distinfo,sampleintervals,runalgos,graphname,realdata,graphevolution,noiselevel,fixedfraction,"-1"]
                 bunchplot=returnbunchplots("sampling_fixedtrace",sentparams)
                 requiredplots.extend(bunchplot)
                 #pr and roc can not be run for specific parameters since the score is estimated over parameters
                 #plotinfo11=["sampling_fixedtrace","pr",sentparams] 
                 #requiredplots.append(plotinfo11)
                 #plotinfo11=["sampling_fixedtrace","roc",sentparams]
                 #requiredplots.append(plotinfo11)

                 noiselevel=1.0
                 fixedfraction=0.05
                 sentparams=["si",distinfo,sampleintervals,runalgos,graphname,realdata,graphevolution,noiselevel,fixedfraction,algoparaminfo]
                 bunchplot=returnbunchplots("sampling_fixedtrace",sentparams)
                 requiredplots.extend(bunchplot)

                 noiselevel=1.0
                 fixedfraction=0.05
                 sentparams=["si",distinfo,sampleintervals,runalgos,graphname,realdata,graphevolution,noiselevel,fixedfraction,"-1"]
                 bunchplot=returnbunchplots("sampling_fixedtrace",sentparams)
                 requiredplots.extend(bunchplot)
                 #pr and roc can not be run for specific parameters since the score is estimated over parameters
                 #plotinfo11=["sampling_fixedtrace","pr",sentparams] 
                 #requiredplots.append(plotinfo11)
                 #plotinfo11=["sampling_fixedtrace","roc",sentparams]
                 #requiredplots.append(plotinfo11)

                 if graphname in ["stanford-contact-static-207-2962-300-1"]:
                    distinfo="spreadprob_0.1-s2i_weibull_3.0_2.5"
                 else:
                    distinfo="spreadprob_0.1-s2i_expo_2.0"
                 noiselevel=1.0
                 fixedfraction=0.1
                 algoparaminfo={}
                 algoparaminfo["lsel1"]=[0.01,1.0]
                 algoparaminfo["absel1"]=[0.01,0.05]
                 algoparaminfo["multitree"]=[1.2]
                 algoparaminfo["netinf"]=[1.2]
                 algoparaminfo["netrate"]=["1"]
                 sentparams=["si",distinfo,sampleintervals,runalgos,graphname,realdata,graphevolution,noiselevel,fixedfraction,algoparaminfo]
                 bunchplot=returnbunchplots("sampling_fixedtrace",sentparams)
                 requiredplots.extend(bunchplot)

                 noiselevel=1.0
                 fixedfraction=0.1
                 sentparams=["si",distinfo,sampleintervals,runalgos,graphname,realdata,graphevolution,noiselevel,fixedfraction,"-1"]
                 bunchplot=returnbunchplots("sampling_fixedtrace",sentparams)
                 requiredplots.extend(bunchplot)
                 #pr and roc can not be run for specific parameters since the score is estimated over parameters
                 #plotinfo11=["sampling_fixedtrace","pr",sentparams] 
                 #requiredplots.append(plotinfo11)
                 #plotinfo11=["sampling_fixedtrace","roc",sentparams]
                 #requiredplots.append(plotinfo11)

                 
                 noiselevel=1.0
                 fixedfraction=0.15
                 sentparams=["si",distinfo,sampleintervals,runalgos,graphname,realdata,graphevolution,noiselevel,fixedfraction,algoparaminfo]
                 bunchplot=returnbunchplots("sampling_fixedtrace",sentparams)
                 requiredplots.extend(bunchplot)

                 noiselevel=1.0
                 fixedfraction=0.15
                 sentparams=["si",distinfo,sampleintervals,runalgos,graphname,realdata,graphevolution,noiselevel,fixedfraction,"-1"]
                 bunchplot=returnbunchplots("sampling_fixedtrace",sentparams)
                 requiredplots.extend(bunchplot)
                 #pr and roc can not be run for specific parameters since the score is estimated over parameters
                 #plotinfo11=["sampling_fixedtrace","pr",sentparams] 
                 #requiredplots.append(plotinfo11)
                 #plotinfo11=["sampling_fixedtrace","roc",sentparams]
                 #requiredplots.append(plotinfo11)

                 
                 noiselevel=1.0
                 fixedfraction=0.2
                 sentparams=["si",distinfo,sampleintervals,runalgos,graphname,realdata,graphevolution,noiselevel,fixedfraction,algoparaminfo]
                 bunchplot=returnbunchplots("sampling_fixedtrace",sentparams)
                 requiredplots.extend(bunchplot)

                 noiselevel=1.0
                 fixedfraction=0.2
                 sentparams=["si",distinfo,sampleintervals,runalgos,graphname,realdata,graphevolution,noiselevel,fixedfraction,"-1"]
                 bunchplot=returnbunchplots("sampling_fixedtrace",sentparams)
                 requiredplots.extend(bunchplot)
                 #pr and roc can not be run for specific parameters since the score is estimated over parameters
                 #plotinfo11=["sampling_fixedtrace","pr",sentparams] 
                 #requiredplots.append(plotinfo11)
                 #plotinfo11=["sampling_fixedtrace","roc",sentparams]
                 #requiredplots.append(plotinfo11)

                 
                 noiselevel=1.0
                 fixedfraction=0.07
                 sentparams=["si",distinfo,sampleintervals,runalgos,graphname,realdata,graphevolution,noiselevel,fixedfraction,algoparaminfo]
                 bunchplot=returnbunchplots("sampling_fixedtrace",sentparams)
                 requiredplots.extend(bunchplot)

                 noiselevel=1.0
                 fixedfraction=0.07
                 sentparams=["si",distinfo,sampleintervals,runalgos,graphname,realdata,graphevolution,noiselevel,fixedfraction,"-1"]
                 bunchplot=returnbunchplots("sampling_fixedtrace",sentparams)
                 requiredplots.extend(bunchplot)
                 #pr and roc can not be run for specific parameters since the score is estimated over parameters
                 #plotinfo11=["sampling_fixedtrace","pr",sentparams] 
                 #requiredplots.append(plotinfo11)
                 #plotinfo11=["sampling_fixedtrace","roc",sentparams]
                 #requiredplots.append(plotinfo11)

                 noiselevel=1.0
                 fixedfraction=0.05
                 sentparams=["si",distinfo,sampleintervals,runalgos,graphname,realdata,graphevolution,noiselevel,fixedfraction,algoparaminfo]
                 bunchplot=returnbunchplots("sampling_fixedtrace",sentparams)
                 requiredplots.extend(bunchplot)

                 noiselevel=1.0
                 fixedfraction=0.05
                 sentparams=["si",distinfo,sampleintervals,runalgos,graphname,realdata,graphevolution,noiselevel,fixedfraction,"-1"]
                 bunchplot=returnbunchplots("sampling_fixedtrace",sentparams)
                 requiredplots.extend(bunchplot)

                 #dist info plot
                 distinfos=["spreadprob_0.1-s2i_expo_2.0","spreadprob_0.1-s2i_expo_3.0","spreadprob_0.1-s2i_expo_4.0","spreadprob_0.1-s2i_expo_6.0","spreadprob_0.1-s2i_expo_8.0"]
                 sampleinterval=1
                 noiselevel=1.0
                 fixedfraction=0.1
                 algoparaminfo={}
                 algoparaminfo["lsel1"]=[0.01,1.0]
                 algoparaminfo["absel1"]=[0.01,0.05]
                 algoparaminfo["multitree"]=[1.2]
                 algoparaminfo["netinf"]=[1.2]
                 algoparaminfo["netrate"]=["1"]
                 sentparams=["si",distinfos,sampleinterval,runalgos,graphname,realdata,graphevolution,noiselevel,fixedfraction,algoparaminfo]
                 bunchplot=returnbunchplots("distinfo_fixedtrace",sentparams)
                 requiredplots.extend(bunchplot)

                 fixedfraction=0.1
                 sentparams=["si",distinfos,sampleinterval,runalgos,graphname,realdata,graphevolution,noiselevel,fixedfraction,"-1"]
                 bunchplot=returnbunchplots("distinfo_fixedtrace",sentparams)
                 requiredplots.extend(bunchplot)
                 #pr and roc can not be run for specific parameters since the score is estimated over parameters
                 #plotinfo11=["distinfo_fixedtrace","pr",sentparams] 
                 #requiredplots.append(plotinfo11)
                 #plotinfo11=["distinfo_fixedtrace","roc",sentparams]
                 #requiredplots.append(plotinfo11)

                 fixedfraction=0.15
                 sentparams=["si",distinfos,sampleinterval,runalgos,graphname,realdata,graphevolution,noiselevel,fixedfraction,algoparaminfo]
                 bunchplot=returnbunchplots("distinfo_fixedtrace",sentparams)
                 requiredplots.extend(bunchplot)

                 fixedfraction=0.15
                 sentparams=["si",distinfos,sampleinterval,runalgos,graphname,realdata,graphevolution,noiselevel,fixedfraction,"-1"]
                 bunchplot=returnbunchplots("distinfo_fixedtrace",sentparams)
                 requiredplots.extend(bunchplot)
                 #pr and roc can not be run for specific parameters since the score is estimated over parameters
                 #plotinfo11=["distinfo_fixedtrace","pr",sentparams] 
                 #requiredplots.append(plotinfo11)
                 #plotinfo11=["distinfo_fixedtrace","roc",sentparams]
                 #requiredplots.append(plotinfo11)

                 fixedfraction=0.2
                 sentparams=["si",distinfos,sampleinterval,runalgos,graphname,realdata,graphevolution,noiselevel,fixedfraction,algoparaminfo]
                 bunchplot=returnbunchplots("distinfo_fixedtrace",sentparams)
                 requiredplots.extend(bunchplot)

                 fixedfraction=0.2
                 sentparams=["si",distinfos,sampleinterval,runalgos,graphname,realdata,graphevolution,noiselevel,fixedfraction,"-1"]
                 bunchplot=returnbunchplots("distinfo_fixedtrace",sentparams)
                 requiredplots.extend(bunchplot)
                 #pr and roc can not be run for specific parameters since the score is estimated over parameters
                 #plotinfo11=["distinfo_fixedtrace","pr",sentparams] 
                 #requiredplots.append(plotinfo11)
                 #plotinfo11=["distinfo_fixedtrace","roc",sentparams]
                 #requiredplots.append(plotinfo11)

                 fixedfraction=0.07
                 sentparams=["si",distinfos,sampleinterval,runalgos,graphname,realdata,graphevolution,noiselevel,fixedfraction,algoparaminfo]
                 bunchplot=returnbunchplots("distinfo_fixedtrace",sentparams)
                 requiredplots.extend(bunchplot)

                 fixedfraction=0.07
                 sentparams=["si",distinfos,sampleinterval,runalgos,graphname,realdata,graphevolution,noiselevel,fixedfraction,"-1"]
                 bunchplot=returnbunchplots("distinfo_fixedtrace",sentparams)
                 requiredplots.extend(bunchplot)
                 #pr and roc can not be run for specific parameters since the score is estimated over parameters
                 #plotinfo11=["distinfo_fixedtrace","pr",sentparams] 
                 #requiredplots.append(plotinfo11)
                 #plotinfo11=["distinfo_fixedtrace","roc",sentparams]
                 #requiredplots.append(plotinfo11)

                 #noise plot generator fixed trace
                 if graphname in ["stanford-contact-static-207-2962-300-1"]:
                    distinfo="spreadprob_0.1-s2i_weibull_3.0_2.5"
                 else:
                    distinfo="spreadprob_0.1-s2i_expo_2.0"
                 sampleinterval=1
                 noiselevels=[1.0,0.9,0.8,0.7,0.6,0.5,0.4,0.3]
                 fixedfraction=0.1
                 algoparaminfo={}
                 algoparaminfo["lsel1"]=[0.01,1.0]
                 algoparaminfo["absel1"]=[0.01,0.05]
                 algoparaminfo["multitree"]=[1.2]
                 algoparaminfo["netinf"]=[1.2]
                 algoparaminfo["netrate"]=["1"]
                 sentparams=["si",distinfo,sampleinterval,runalgos,graphname,realdata,graphevolution,noiselevels,fixedfraction,algoparaminfo]
                 bunchplot=returnbunchplots("noise_fixedtrace",sentparams)
                 requiredplots.extend(bunchplot)

                 fixedfraction=0.1
                 sentparams=["si",distinfo,sampleinterval,runalgos,graphname,realdata,graphevolution,noiselevels,fixedfraction,"-1"]
                 bunchplot=returnbunchplots("noise_fixedtrace",sentparams)
                 requiredplots.extend(bunchplot)
                 #pr and roc can not be run for specific parameters since the score is estimated over parameters
                 #plotinfo11=["noise_fixedtrace","pr",sentparams] 
                 #requiredplots.append(plotinfo11)
                 #plotinfo11=["noise_fixedtrace","roc",sentparams]
                 #requiredplots.append(plotinfo11)

                 fixedfraction=0.15
                 sentparams=["si",distinfo,sampleinterval,runalgos,graphname,realdata,graphevolution,noiselevels,fixedfraction,algoparaminfo]
                 bunchplot=returnbunchplots("noise_fixedtrace",sentparams)
                 requiredplots.extend(bunchplot)

                 fixedfraction=0.15
                 sentparams=["si",distinfo,sampleinterval,runalgos,graphname,realdata,graphevolution,noiselevels,fixedfraction,"-1"]
                 bunchplot=returnbunchplots("noise_fixedtrace",sentparams)
                 requiredplots.extend(bunchplot)
                 #pr and roc can not be run for specific parameters since the score is estimated over parameters
                 #plotinfo11=["noise_fixedtrace","pr",sentparams] 
                 #requiredplots.append(plotinfo11)
                 #plotinfo11=["noise_fixedtrace","roc",sentparams]
                 #requiredplots.append(plotinfo11)

                 fixedfraction=0.2
                 sentparams=["si",distinfo,sampleinterval,runalgos,graphname,realdata,graphevolution,noiselevels,fixedfraction,algoparaminfo]
                 bunchplot=returnbunchplots("noise_fixedtrace",sentparams)
                 requiredplots.extend(bunchplot)

                 fixedfraction=0.2
                 sentparams=["si",distinfo,sampleinterval,runalgos,graphname,realdata,graphevolution,noiselevels,fixedfraction,"-1"]
                 bunchplot=returnbunchplots("noise_fixedtrace",sentparams)
                 requiredplots.extend(bunchplot)
                 #pr and roc can not be run for specific parameters since the score is estimated over parameters
                 #plotinfo11=["noise_fixedtrace","pr",sentparams] 
                 #requiredplots.append(plotinfo11)
                 #plotinfo11=["noise_fixedtrace","roc",sentparams]
                 #requiredplots.append(plotinfo11)

                 fixedfraction=0.07
                 sentparams=["si",distinfo,sampleinterval,runalgos,graphname,realdata,graphevolution,noiselevels,fixedfraction,algoparaminfo]
                 bunchplot=returnbunchplots("noise_fixedtrace",sentparams)
                 requiredplots.extend(bunchplot)

                 fixedfraction=0.07
                 sentparams=["si",distinfo,sampleinterval,runalgos,graphname,realdata,graphevolution,noiselevels,fixedfraction,"-1"]
                 bunchplot=returnbunchplots("noise_fixedtrace",sentparams)
                 requiredplots.extend(bunchplot)
                 #pr and roc can not be run for specific parameters since the score is estimated over parameters
                 #plotinfo11=["noise_fixedtrace","pr",sentparams] 
                 #requiredplots.append(plotinfo11)
                 #plotinfo11=["noise_fixedtrace","roc",sentparams]
                 #requiredplots.append(plotinfo11)
       elif (realdata=="real" and graphevolution=="dynamic") or (realdata=="synthetic" and graphevolution=="dynamic"):
          if (realdata=="real" and graphevolution=="dynamic"):
             graphnames=realdynamicgraphnames
          if(realdata=="synthetic" and graphevolution=="dynamic"):
             graphnames=syndynamicgraphnames
          runalgos=["lsel1sumfused","absel1sumfused"] #"lsel1mulfused","absel1mulfused"
          if spreadmodel in ["si"]:
             for graphname in graphnames:
                 if graphname.find("stanford")!=-1 or graphname.find("Stanford")!=-1:
                    maindistinfo=["spreadprob_0.5-s2i_weibull_3.0_2.5","spreadprob_0.1-s2i_weibull_3.0_2.5"]
                 else:
                    maindistinfo=["spreadprob_0.5-s2i_expo_2.0","spreadprob_0.1-s2i_expo_2.0"] 
                 #overtime plot
                 distinfo=maindistinfo[0]
                 sampleinterval=1
                 noiselevel=1.0
                 algoparaminfo={}
                 algoparaminfo["lsel1fused"]=[0.01,0.005,0.05]
                 algoparaminfo["absel1fused"]=[0.01,0.05]
                 sentparams=["si",distinfo,sampleinterval,runalgos,graphname,realdata,graphevolution,noiselevel,algoparaminfo]
                 bunchplot=returnbunchplots("overtime",sentparams)
                 requiredplots.extend(bunchplot)

                 sampleinterval=1
                 noiselevel=1.0
                 sentparams=["si",distinfo,sampleinterval,runalgos,graphname,realdata,graphevolution,noiselevel,"-1"]
                 bunchplot=returnbunchplots("overtime",sentparams)
                 requiredplots.extend(bunchplot)
                 #pr and roc can not be run for specific parameters since the score is estimated over parameters
                 #plotinfo11=["overtime","pr",sentparams] 
                 #requiredplots.append(plotinfo11)
                 #plotinfo11=["overtime","roc",sentparams]
                 #requiredplots.append(plotinfo11)

                 distinfo=maindistinfo[1]
                 sampleinterval=1
                 noiselevel=1.0
                 algoparaminfo={}
                 algoparaminfo["lsel1fused"]=[0.01,0.005,0.05]
                 algoparaminfo["absel1fused"]=[0.01,0.05]
                 sentparams=["si",distinfo,sampleinterval,runalgos,graphname,realdata,graphevolution,noiselevel,algoparaminfo]
                 bunchplot=returnbunchplots("overtime",sentparams)
                 requiredplots.extend(bunchplot)

                 sampleinterval=1
                 noiselevel=1.0
                 sentparams=["si",distinfo,sampleinterval,runalgos,graphname,realdata,graphevolution,noiselevel,"-1"]
                 bunchplot=returnbunchplots("overtime",sentparams)
                 requiredplots.extend(bunchplot)
                 #pr and roc can not be run for specific parameters since the score is estimated over parameters
                 #plotinfo11=["overtime","pr",sentparams] 
                 #requiredplots.append(plotinfo11)
                 #plotinfo11=["overtime","roc",sentparams]
                 #requiredplots.append(plotinfo11)

                 if False:
                  #parameter auc curve 
                  distinfo=maindistinfo[0]
                  sampleinterval=1
                  noiselevel=1.0
                  fixedfraction=0.15
                  sentparams=["si",distinfo,sampleinterval,runalgos,graphname,realdata,graphevolution,noiselevel,fixedfraction]
                  myplotinfo=["parameterauccurve","roc",sentparams]
                  requiredplots.append(myplotinfo)
                  #myplotinfo=["parameterauccurve","bedroc",sentparams]
                  #requiredplots.append(myplotinfo)
                  myplotinfo=["parameterauccurve","pr",sentparams]
                  requiredplots.append(myplotinfo)
          
                  distinfo=maindistinfo[1]
                  sampleinterval=1
                  noiselevel=1.0
                  fixedfraction=0.15
                  sentparams=["si",distinfo,sampleinterval,runalgos,graphname,realdata,graphevolution,noiselevel,fixedfraction]
                  myplotinfo=["parameterauccurve","roc",sentparams]
                  requiredplots.append(myplotinfo)
                  #myplotinfo=["parameterauccurve","bedroc",sentparams]
                  #requiredplots.append(myplotinfo)
                  myplotinfo=["parameterauccurve","pr",sentparams]
                  requiredplots.append(myplotinfo)
                 
                  #parameter auc curve 
                  distinfo=maindistinfo[0]
                  sampleinterval=1
                  noiselevel=1.0
                  fixedfraction=0.1
                  sentparams=["si",distinfo,sampleinterval,runalgos,graphname,realdata,graphevolution,noiselevel,fixedfraction]
                  myplotinfo=["parameterauccurve","roc",sentparams]
                  requiredplots.append(myplotinfo)
                  #myplotinfo=["parameterauccurve","bedroc",sentparams]
                  #requiredplots.append(myplotinfo)
                  myplotinfo=["parameterauccurve","pr",sentparams]
                  requiredplots.append(myplotinfo)
          
                  distinfo=maindistinfo[1]
                  sampleinterval=1
                  noiselevel=1.0
                  fixedfraction=0.1
                  sentparams=["si",distinfo,sampleinterval,runalgos,graphname,realdata,graphevolution,noiselevel,fixedfraction]
                  myplotinfo=["parameterauccurve","roc",sentparams]
                  requiredplots.append(myplotinfo)
                  #myplotinfo=["parameterauccurve","bedroc",sentparams]
                  #requiredplots.append(myplotinfo)
                  myplotinfo=["parameterauccurve","pr",sentparams]
                  requiredplots.append(myplotinfo)

                  #parameter auc curve 
                  distinfo=maindistinfo[0]
                  sampleinterval=1
                  noiselevel=1.0
                  fixedfraction=0.07
                  sentparams=["si",distinfo,sampleinterval,runalgos,graphname,realdata,graphevolution,noiselevel,fixedfraction]
                  myplotinfo=["parameterauccurve","roc",sentparams]
                  requiredplots.append(myplotinfo)
                  #myplotinfo=["parameterauccurve","bedroc",sentparams]
                  #requiredplots.append(myplotinfo)
                  myplotinfo=["parameterauccurve","pr",sentparams]
                  requiredplots.append(myplotinfo)
          
                  distinfo=maindistinfo[1]
                  sampleinterval=1
                  noiselevel=1.0
                  fixedfraction=0.07
                  sentparams=["si",distinfo,sampleinterval,runalgos,graphname,realdata,graphevolution,noiselevel,fixedfraction]
                  myplotinfo=["parameterauccurve","roc",sentparams]
                  requiredplots.append(myplotinfo)
                  #myplotinfo=["parameterauccurve","bedroc",sentparams]
                  #requiredplots.append(myplotinfo)
                  myplotinfo=["parameterauccurve","pr",sentparams]
                  requiredplots.append(myplotinfo)

                  #parameter auc curve 
                  distinfo=maindistinfo[0]
                  sampleinterval=1
                  noiselevel=1.0
                  fixedfraction=0.05
                  sentparams=["si",distinfo,sampleinterval,runalgos,graphname,realdata,graphevolution,noiselevel,fixedfraction]
                  myplotinfo=["parameterauccurve","roc",sentparams]
                  requiredplots.append(myplotinfo)
                  #myplotinfo=["parameterauccurve","bedroc",sentparams]
                  #requiredplots.append(myplotinfo)
                  myplotinfo=["parameterauccurve","pr",sentparams]
                  requiredplots.append(myplotinfo)
          
                  distinfo=maindistinfo[1]
                  sampleinterval=1
                  noiselevel=1.0
                  fixedfraction=0.05
                  sentparams=["si",distinfo,sampleinterval,runalgos,graphname,realdata,graphevolution,noiselevel,fixedfraction]
                  myplotinfo=["parameterauccurve","roc",sentparams]
                  requiredplots.append(myplotinfo)
                  #myplotinfo=["parameterauccurve","bedroc",sentparams]
                  #requiredplots.append(myplotinfo)
                  myplotinfo=["parameterauccurve","pr",sentparams]
                  requiredplots.append(myplotinfo)
             
                  #parameter auc curve 
                  distinfo=maindistinfo[0]
                  sampleinterval=1
                  noiselevel=1.0
                  fixedfraction=0.2
                  sentparams=["si",distinfo,sampleinterval,runalgos,graphname,realdata,graphevolution,noiselevel,fixedfraction]
                  myplotinfo=["parameterauccurve","roc",sentparams]
                  requiredplots.append(myplotinfo)
                  #myplotinfo=["parameterauccurve","bedroc",sentparams]
                  #requiredplots.append(myplotinfo)
                  myplotinfo=["parameterauccurve","pr",sentparams]
                  requiredplots.append(myplotinfo)
          
                  distinfo=maindistinfo[1]
                  sampleinterval=1
                  noiselevel=1.0
                  fixedfraction=0.2
                  sentparams=["si",distinfo,sampleinterval,runalgos,graphname,realdata,graphevolution,noiselevel,fixedfraction]
                  myplotinfo=["parameterauccurve","roc",sentparams]
                  requiredplots.append(myplotinfo)
                  #myplotinfo=["parameterauccurve","bedroc",sentparams]
                  #requiredplots.append(myplotinfo)
                  myplotinfo=["parameterauccurve","pr",sentparams]
                  requiredplots.append(myplotinfo)

                  #dist info plot
                  #right now we don't do for dynamic graphs
                 
                 #noise plot generator fixed trace
                 distinfo=maindistinfo[1]   
                 sampleinterval=1
                 noiselevels=[1.0,0.9,0.8,0.7,0.6,0.5,0.4,0.3]
                 fixedfraction=0.1
                 algoparaminfo={}
                 algoparaminfo["lsel1fused"]=[0.01,0.005,0.05]
                 algoparaminfo["absel1fused"]=[0.01,0.05]
                 sentparams=["si",distinfo,sampleinterval,runalgos,graphname,realdata,graphevolution,noiselevels,fixedfraction,algoparaminfo]
                 bunchplot=returnbunchplots("noise_fixedtrace",sentparams)
                 requiredplots.extend(bunchplot)

                 fixedfraction=0.1
                 algoparaminfo="-1"
                 sentparams=["si",distinfo,sampleinterval,runalgos,graphname,realdata,graphevolution,noiselevels,fixedfraction,"-1"]
                 bunchplot=returnbunchplots("noise_fixedtrace",sentparams)
                 requiredplots.extend(bunchplot)
                 #pr and roc can not be run for specific parameters since the score is estimated over parameters
                 plotinfo11=["noise_fixedtrace","pr",sentparams] 
                 requiredplots.append(plotinfo11)
                 plotinfo11=["noise_fixedtrace","roc",sentparams]
                 requiredplots.append(plotinfo11)

                 fixedfraction=0.15
                 sentparams=["si",distinfo,sampleinterval,runalgos,graphname,realdata,graphevolution,noiselevels,fixedfraction,algoparaminfo]
                 bunchplot=returnbunchplots("noise_fixedtrace",sentparams)
                 requiredplots.extend(bunchplot)

                 fixedfraction=0.15
                 algoparaminfo="-1"
                 sentparams=["si",distinfo,sampleinterval,runalgos,graphname,realdata,graphevolution,noiselevels,fixedfraction,algoparaminfo]
                 bunchplot=returnbunchplots("noise_fixedtrace",sentparams)
                 requiredplots.extend(bunchplot)
                 #pr and roc can not be run for specific parameters since the score is estimated over parameters
                 plotinfo11=["noise_fixedtrace","pr",sentparams] 
                 requiredplots.append(plotinfo11)
                 plotinfo11=["noise_fixedtrace","roc",sentparams]
                 requiredplots.append(plotinfo11)

                 fixedfraction=0.2
                 sentparams=["si",distinfo,sampleinterval,runalgos,graphname,realdata,graphevolution,noiselevels,fixedfraction,algoparaminfo]
                 bunchplot=returnbunchplots("noise_fixedtrace",sentparams)
                 requiredplots.extend(bunchplot)

                 fixedfraction=0.2
                 algoparaminfo="-1"
                 sentparams=["si",distinfo,sampleinterval,runalgos,graphname,realdata,graphevolution,noiselevels,fixedfraction,algoparaminfo]
                 bunchplot=returnbunchplots("noise_fixedtrace",sentparams)
                 requiredplots.extend(bunchplot)
                 #pr and roc can not be run for specific parameters since the score is estimated over parameters
                 plotinfo11=["noise_fixedtrace","pr",sentparams] 
                 requiredplots.append(plotinfo11)
                 plotinfo11=["noise_fixedtrace","roc",sentparams]
                 requiredplots.append(plotinfo11)

                 fixedfraction=0.07
                 sentparams=["si",distinfo,sampleinterval,runalgos,graphname,realdata,graphevolution,noiselevels,fixedfraction,algoparaminfo]
                 bunchplot=returnbunchplots("noise_fixedtrace",sentparams)
                 requiredplots.extend(bunchplot)

                 fixedfraction=0.07
                 algoparaminfo="-1"
                 sentparams=["si",distinfo,sampleinterval,runalgos,graphname,realdata,graphevolution,noiselevels,fixedfraction,algoparaminfo]
                 bunchplot=returnbunchplots("noise_fixedtrace",sentparams)
                 requiredplots.extend(bunchplot)
                 #pr and roc can not be run for specific parameters since the score is estimated over parameters
                 plotinfo11=["noise_fixedtrace","pr",sentparams] 
                 requiredplots.append(plotinfo11)
                 plotinfo11=["noise_fixedtrace","roc",sentparams]
                 requiredplots.append(plotinfo11)
                 
       return requiredplots
    else:
       print "This infertype is unknown {0}".format(infertype) 
       exit(1)
    return requiredplots


def graphplotnameassigner(graphname):        
    if graphname.find("lpa")!=-1:
       globals()["graph2plotname"][graphname]="Lpa" 
    elif graphname.find("dmc")!=-1:
       globals()["graph2plotname"][graphname]="Dmc"
    elif graphname.find("rds")!=-1:
       globals()["graph2plotname"][graphname]="Rds"
    elif graphname.find("smw")!=-1:
       globals()["graph2plotname"][graphname]="Smw"
    elif graphname.find("undirected_ff")!=-1:
       globals()["graph2plotname"][graphname]="Forest Fire"       
    elif graphname.find("enron")!=-1 or graphname.find("Enron")!=-1:
       globals()["graph2plotname"][graphname]="Enron Email"
    elif graphname.find("stanford")!=-1 or graphname.find("Stanford")!=-1:
       globals()["graph2plotname"][graphname]="Human Contact"
    else:
       print "this graphanem i unnknonw!! {0}".format(graphname)
       exit(1)

def returncounts(sen,fpr,acc,totalcount):
    A=np.zeros((4,4),dtype=np.float)
    A[0,0]=1.0-sen
    A[0,2]=-1.0*sen 
    A[1,1]=1.0-fpr
    A[1,3]=-1.0*fpr
    A[2,0]=1.0-acc
    A[2,1]=-1.0*acc
    A[2,2]=-1.0*acc
    A[2,3]=1.0-acc
    A[3,0]=1.0
    A[3,1]=1.0
    A[3,2]=1.0
    A[3,3]=1.0
    right=[0.0,0.0,0.0,float(totalcount)]
    assert np.linalg.matrix_rank(A)==4
    ret=np.linalg.solve(A,right)
    return ret #tp,fp,fn,tn

def assignscores(tp,fp,fn,tn):                                    
    sen=float(tp)/(tp+fn)
    fpr=float(fp)/(fp+tn)
    recall=sen
    if (tp+fp)!=0:
       precision=float(tp)/(tp+fp)
    else:
       precision=0.0
    acc=float(tp+tn)/(tp+tn+fn+fp)
    spec=1.0-fpr #true negative rate

    #there is no true negative rate on f formula, but there is on mcc
    if precision+recall!=0.0:
       f1score=float(2.0*precision*recall)/(precision+recall)
       f2score=float(5.0*precision*recall)/((4.0*precision)+recall)
       #avgnodenum=0.0
       #for time in G.keys():
       #    avgnodenum+=G[time].number_of_nodes()
       #avgnodenum/=float(len(G.keys()))    
       #beta=1.0/float(avgnodenum)
       #betasq=beta**2
       #f1overnscore=float((1.0+betasq)*precision*recall)/((betasq*precision)+recall)
    else:
       f1score=0.0
       f2score=0.0
       #f1overnscore=0.0

    if ((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))!=0:
       mcc=float((tp*tn)-(fp*fn))/((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))  
    else:
       mcc=0.0
       
    scores={}
    scores["sen"]=sen
    scores["fpr"]=fpr
    scores["recall"]=recall
    scores["precision"]=precision
    scores["acc"]=acc
    scores["spec"]=spec
    scores["f1"]=f1score
    scores["f2"]=f2score
    #scores["f1overn"]=f1overnscore
    scores["mcc"]=mcc
    return scores

def estimateaverage(allscores):       
       allscores2={}
       newkey2files={}
       if realdata=="synthetic" and graphevolution=="dynamic":
          pass
       if realdata=="synthetic" and graphevolution=="graphwise":
          pass 
       else:
          print "Average is not defined for synthetic graphs {0} {1}".format(realdata,graphevolution)
          exit(1)
           
       for newkey in newkey2files.keys():
           tempscores[newkey]={}
           for filepathkey in newkey2files.keys():
               for traceparam in allscores[filepathkey].keys():
                   if not tempscores[newkey].has_key(traceparam):
                       tempscores[newkey][traceparam]={}
                   for inferalgo in allscores[filepathkey][traceparam].keys():
                       if not tempscores[newkey][traceparam].has_key(inferalgo):
                          tempscores[newkey][traceparam][inferalgo]={}
                       for sampleinterval in allscores[filepathkey][traceparam][inferalgo].keys():
                           if not tempscores[newkey][traceparam][inferalgo].has_key(sampleinterval):
                              tempscores[newkey][traceparam][inferalgo][sampleinterval]={}
                           for fraction in allscores[filepathkey][traceparam][inferalgo][sampleinterval].keys():
                               if not tempscores[newkey][traceparam][inferalgo][sampleinterval].has_key(fraction):
                                  tempscores[newkey][traceparam][inferalgo][sampleinterval][fraction]={}
                               for noiseratio in allscores[filepathkey][traceparam][inferalgo][sampleinterval][fraction].keys():
                                   if not tempscores[newkey][traceparam][inferalgo][sampleinterval][fraction].has_key(noiseratio):
                                      tempscores[newkey][traceparam][inferalgo][sampleinterval][fraction][noiseratio]={}
                                   for paramstr in allscores[filepathkey][traceparam][inferalgo][sampleinterval][fraction][noiseratio].keys():
                                       if not tempscores[newkey][traceparam][inferalgo][sampleinterval][fraction][noiseratio].has_key(paramstr):
                                          tempscores[filepathkey][traceparam][inferalgo][sampleinterval][fraction][noiseratio][paramstr]={}
                                       tempscores.append(dict(allscores[filepathkey][traceparam][inferalgo][sampleinterval][fraction][noiseratio][paramstr]))
           allscores={}
           for filepathkey in tempscores.keys():
               for traceparam in tempscores[filepathkey].keys():
                   for inferalgo in tempscores[filepathkey][traceparam].keys():
                       for sampleinterval in tempscores[filepathkey][traceparam][inferalgo].keys():
                           for fraction in tempscores[filepathkey][traceparam][inferalgo][sampleinterval].keys():
                               for noiseratio in tempscores[filepathkey][traceparam][inferalgo][sampleinterval][fraction].keys():
                                   for paramstr in tempscores[filepathkey][traceparam][inferalgo][sampleinterval][fraction][noiseratio].keys():
                                       results=tempscores[filepathkey][traceparam][inferalgo][sampleinterval][fraction][noiseratio][paramstr]
                                       avgscore={}
                                       for keyname in results[0].keys():
                                           avgscore[keyname]=0.0
                                       incount=0    
                                       for result in results:
                                           if result=="-1":
                                              continue
                                           else:
                                              incount+=1 
                                              for keyname in result.keys():
                                                  avgscore[keyname]+=result[keyname]
                                       if incount!=0:
                                          for keyname in avgscore.keys():
                                              avgscore[keyname]/=float(incount)
                                       else:
                                          avgscore="-1"
                                       allscores[filepathkey][traceparam][inferalgo][sampleinterval][fraction][noiseratio][paramstr]=avgscore
       return allscores2

def estimatefractionwisescores(spreamodel,distinfo,sampleinterval,runalgos,graphname,realdata,graphevolution,noiselevel,allscores): 
           sentscores={}
           filepathkey=graphname
           for runalgo in runalgos:
                  if runalgo.find("sum")!=-1:
                     extend="sum"
                     sentalgo=runalgo.replace("sum","")
                  elif runalgo.find("multi")==-1 and runalgo.find("mul")!=-1:
                     extend="mul"
                     sentalgo=runalgo.replace("mul","")
                  else: 
                     extend=""
                     sentalgo=runalgo
                  if not sentscores.has_key(sentalgo):
                     sentscores[sentalgo]={}

                  print "top info"
                  print runalgo
                  print sampleinterval
                  print sentalgo
                  if not allscores[filepathkey][distinfo].has_key(runalgo):
                     continue 
                  print "fractions {0}".format(allscores[filepathkey][distinfo][runalgo][sampleinterval].keys())
                  for fraction in allscores[filepathkey][distinfo][runalgo][sampleinterval].keys():
                      if fraction not in requiredfractions:  ##??
                         continue 
                      print "infom"
                      print fraction
                      print filepathkey
                      print runalgo
                      print sampleinterval
                      print distinfo
                      print noiselevel
                      print allscores[filepathkey][distinfo][runalgo][sampleinterval][fraction].keys()
                      if not allscores[filepathkey][distinfo][runalgo][sampleinterval][fraction].has_key(noiselevel):
                         continue
                      for paramstr in allscores[filepathkey][distinfo][runalgo][sampleinterval][fraction][noiselevel].keys():
                          sentparamstr="{0}+{1}".format(paramstr,extend)   
                          if not sentscores[sentalgo].has_key(sentparamstr):
                             sentscores[sentalgo][sentparamstr]={}
                          sentscores[sentalgo][sentparamstr][fraction]=allscores[filepathkey][distinfo][runalgo][sampleinterval][fraction][noiselevel][paramstr]
           #eliminate results that do not have results for all required fractions
           for sentalgo in sentscores.keys():
               myflag=False
               for sentparamstr in sentscores[sentalgo].keys():
                   tempfractions=set(sentscores[sentalgo][sentparamstr].keys())
                   removeitems=set()
                   for frac in tempfractions:
                       if sentscores[sentalgo][sentparamstr][frac]=="-1":
                          removeitems.add(frac)
                   tempfractions-=removeitems 
                   if len(set(requiredfractions).difference(tempfractions))!=0:
                      del sentscores[sentalgo][sentparamstr]
           for sentalgo in sentscores.keys():
               for sentparamstr in sentscores[sentalgo].keys():
                   if len(sentscores[sentalgo][sentparamstr])==0:
                      del sentscores[sentalgo][sentparamstr] 
               if len(sentscores[sentalgo])==0:
                  del sentscores[sentalgo]    
           return sentscores


def sentscores2paramwise(sentscores):
    sentscores2={}
    for sentalgo in sentscores.keys():
        assert not sentscores2.has_key(sentalgo)
        sentscores2[sentalgo]={}
        for parameterstr in sentscores[sentalgo].keys():
            parts=parameterstr.replace("result_","").split("_")
            print parts
            if sentalgo in ["netinf","multitree"]:
               if parts[-1].find("categorical")!=-1:
                  foundone="categorical?" 
               elif parts[-1].find("error")!=-1: 
                  foundone="error?" 
               else:
                  print "not known rounding method!!"
                  exit(1)
               curthres=float(parts[-1].replace(foundone,"").replace("+",""))
               reststr=parts[-1].replace(str(curthres),"")
            elif sentalgo in ["netrate"]:
               curthres="1" #there won't be param for netrate
               reststr=parts[-1]
            elif sentalgo in ["lsel1","absel1"]: #our algos
               curthres=float(parts[-1].split("?")[0])
               reststr="?".join(parts[-1].split("?")[1:])
            elif sentalgo in ["lsel1fused","absel1fused"]: #our algos dynamic
               curthres=float(parts[-1].split("?")[0])
               reststr="?".join(parts[-1].split("?")[1:])   
            else:
               print "this algo is unknown {0}".format(sentalgo)
               exit(1)
            replacedstr="result_"+"_".join(parts[0:-1])+"_"+reststr #this is very important
            #sentparamstr="{0}+{1}".format(replacedstr,extend)
            sentparamstr=replacedstr
                      
            if not sentscores2[sentalgo].has_key(sentparamstr):
               sentscores2[sentalgo][sentparamstr]={}
            for varparam in sentscores[sentalgo][parameterstr].keys(): #varparam can be fraction,plotdistinfo etc
                if not sentscores2[sentalgo][sentparamstr].has_key(varparam):
                   sentscores2[sentalgo][sentparamstr][varparam]={}   
                assert not sentscores2[sentalgo][sentparamstr][varparam].has_key(curthres)  
                sentscores2[sentalgo][sentparamstr][varparam][curthres]=sentscores[sentalgo][parameterstr][varparam]
                #sentscores[sentalgo][sentparamstr][plotdistinfo]=allscores[filepathkey][distinfo][runalgo][sampleinterval][fixedfraction][noiselevel][paramstr]
    return sentscores2        

          
def estimateparamwisescores(spreamodel,distinfo,sampleinterval,runalgos,graphname,realdata,graphevolution,noiselevel,fixedfraction,allscores):      
    sentscores={}
    filepathkey=graphname
    seenthres=set()
    for runalgo in runalgos:
        if runalgo.find("sum")!=-1:
           extend="sum"
           sentalgo=runalgo.replace("sum","")
        elif runalgo.find("multi")==-1 and runalgo.find("mul")!=-1:
           extend="mul"
           sentalgo=runalgo.replace("mul","")
        else: 
           extend=""
           sentalgo=runalgo

        #print "inside info:"             
        #print filepathkey
        #print runalgo
        #print sampleinterval
        #print fixedfraction
        if not allscores[filepathkey][distinfo].has_key(runalgo):
           continue 
        print "alfractions are {0}".format( allscores[filepathkey][distinfo][runalgo][sampleinterval].keys())
        if not allscores[filepathkey][distinfo][runalgo][sampleinterval].has_key(fixedfraction):
           continue
        if not allscores[filepathkey][distinfo][runalgo][sampleinterval][fixedfraction].has_key(noiselevel):
           continue 
        if not sentscores.has_key(sentalgo):
           sentscores[sentalgo]={}
                  
        for paramstr in allscores[filepathkey][distinfo][runalgo][sampleinterval][fixedfraction][noiselevel].keys():
            parts=paramstr.replace("result_","").split("_")
            print parts
            if sentalgo in ["netinf","multitree"]:
               if parts[-1].find("categorical")!=-1:
                  foundone="categorical?" 
               elif parts[-1].find("error")!=-1: 
                  foundone="error?" 
               else:
                  print "not known rounding method!!"
                  exit(1)
               curthres=float(parts[-1].replace(foundone,""))   
               if curthres not in requirededgethres:
                  continue 
               reststr=parts[-1].replace(str(curthres),"")
            elif sentalgo in ["netrate"]:
               curthres="1" #there won't be param for netrate
               reststr=parts[-1]
            elif sentalgo in ["lsel1","absel1"]: #our algos
               curthres=float(parts[-1].split("?")[0])
               if curthres not in requiredlambdathres:
                  continue 
               reststr="?".join(parts[-1].split("?")[1:])
            elif sentalgo in ["lsel1fused","absel1fused"]: #our algos
               curthres=float(parts[-1].split("?")[0])
               if curthres not in requiredlambdathres:
                  continue 
               reststr="?".join(parts[-1].split("?")[1:])   
            else:
               print "this algo is unknown {0}".format(sentalgo)
               exit(1)
            replacedstr="result_"+"_".join(parts[0:-1])+"_"+reststr #this is very important
            sentparamstr="{0}+{1}".format(replacedstr,extend)
                      
            if not sentscores[sentalgo].has_key(sentparamstr):
               sentscores[sentalgo][sentparamstr]={}
            print sentalgo
            print replacedstr
            print curthres
            assert not sentscores[sentalgo][sentparamstr].has_key(curthres)  
            sentscores[sentalgo][sentparamstr][curthres]=allscores[filepathkey][distinfo][runalgo][sampleinterval][fixedfraction][noiselevel][paramstr] 
    #eliminate results that do not have results for all required fractions
    for sentalgo in sentscores.keys():
        for sentparamstr in sentscores[sentalgo].keys():
            allthres=set(sentscores[sentalgo][sentparamstr].keys())
            removeitems=set()
            for tempthres in allthres:
                if sentscores[sentalgo][sentparamstr][tempthres]=="-1":
                   removeitems.add(tempthres)
            allthres-=removeitems
            for removethres in removeitems: # important!!
                del sentscores[sentalgo][sentparamstr][removethres]
            #if sentalgo in ["netinf","multitree"]:
            #   if len(set(requirededgethres).difference(allthres))!=0:
            #      del sentscores[sentalgo][sentparamstr]
            #elif sentalgo in ["lsel1","absel1"]: #our algo
            #   if len(set(requiredlambdathres).difference(allthres))!=0:
            #      del sentscores[sentalgo][sentparamstr]
            #elif sentalgo in ["netrate"]:
            #   assert len(allthres)<=1
            #   if len(allthres)==0:
            #      del sentscores[sentalgo][sentparamstr] 
            #else:
            #   print "sentalgo is unknown{0}".format(sentalgo)
            #   exit(1)
    for sentalgo in sentscores.keys():
        for sentparamstr in sentscores[sentalgo].keys():
            if len(sentscores[sentalgo][sentparamstr])==0:
               del sentscores[sentalgo][sentparamstr] 
        if len(sentscores[sentalgo])==0:
           del sentscores[sentalgo]
    return sentscores


def returnscores(infofilepaths,curdistinfos,cursampleintervals,currunalgos,curnoiselevels,curgraphnames,curfractions,algoparaminfo="-1"):   
    allscores={}
    for filepathkey in infofilepaths.keys():
        if filepathkey not in curgraphnames:
           continue 
        if filepathkey.find("offlineresults")!=-1:
           continue
        #graphplotnameassigner(filepathkey)
        [resultfilepath,graphfilepath]=infofilepaths[filepathkey]
        allscores[filepathkey]={}
        for traceparam in myutil.listdirectories(resultfilepath):
            print traceparam
            if traceparam not in curdistinfos:
               continue
            tracedirname="{0}/{1}".format(resultfilepath,traceparam)
            allscores[filepathkey][traceparam]={} 
            for inferalgo in myutil.listdirectories(tracedirname):
                if inferalgo not in currunalgos:
                   continue 
                algopath="{0}/{1}".format(tracedirname,inferalgo)
                allscores[filepathkey][traceparam][inferalgo]={}
                for sampleinterval in myutil.listdirectories(algopath):
                    sampleinterval=int(sampleinterval)
                    if sampleinterval not in cursampleintervals:
                       continue
                    samplepath="{0}/{1}".format(algopath,sampleinterval)
                    allscores[filepathkey][traceparam][inferalgo][sampleinterval]={}
                    for fraction in myutil.listdirectories(samplepath):
                        fraction=float(fraction)
                        if fraction not in curfractions:
                           continue 
                        fractionpath="{0}/{1}".format(samplepath,fraction)
                        allscores[filepathkey][traceparam][inferalgo][sampleinterval][fraction]={}
                        #find all possible noise ratios
                        noiseratios=set()
                        for resultfilename in myutil.listfiles(fractionpath):
                            parts=resultfilename.split("_")
                            noiseratio=float(parts[3].replace("-randomnoise",""))
                            if noiseratio not in curnoiselevels:
                               continue 
                            noiseratios.add(noiseratio)
                            assert len(parts)==7
                        for noiseratio in noiseratios:
                            allscores[filepathkey][traceparam][inferalgo][sampleinterval][fraction][noiseratio]={}
                        #average score for each parameter
                        #when there is not enough repeat, score will be "-1"    
                        for noiseratio in allscores[filepathkey][traceparam][inferalgo][sampleinterval][fraction].keys():
                            uniqueones={}
                            for resultfilename in myutil.listfiles(fractionpath):
                                parts=resultfilename.split("_")
                                assert len(parts)==7
                                tempnoiseratio=float(parts[3].replace("-randomnoise",""))
                                if noiseratio!=tempnoiseratio:
                                   continue
                                if algoparaminfo=="-1":
                                   pass
                                elif not algoparaminfo.has_key(inferalgo):
                                   pass
                                else:
                                   if inferalgo in ["netinf","multitree"]:
                                      if parts[-2].find("categorical")!=-1:
                                         foundone="categorical?" 
                                      elif parts[-2].find("error")!=-1: 
                                         foundone="error?" 
                                      else:
                                         print "not known rounding method!!"
                                         exit(1)
                                      curthres=float(parts[-2].replace(foundone,""))
                                   elif inferalgo in ["netrate"]:
                                      curthres="1" #there won't be param for netrate
                                   elif inferalgo in ["lsel1","absel1"]: #our algos
                                      curthres=float(parts[-2].split("?")[0])
                                   elif inferalgo in ["lsel1fused","absel1fused"]: #our algos
                                      curthres=float(parts[-2].split("?")[0])
                                   else:
                                      print "this algo is unknown {0}".format(inferalgo)
                                      exit(1)
                                   if curthres not in algoparaminfo[inferalgo]:
                                      continue
                                #print "startinfo"
                                #print fractionpath
                                #print resultfilename
                                #print filepathkey
                                #print noiseratio
                                #print fraction
                                #print sampleinterval
                                #print inferalgo
                                #print traceparam
                                #print parts[-1]
                                #print parts

                                #print fractionpath
                                #print resultfilename

                                repeatindex=int(parts[-1])
                                paramstr="_".join(parts[0:-1])
                                if not uniqueones.has_key(paramstr): 
                                   uniqueones[paramstr]=set()
                                uniqueones[paramstr].add(repeatindex)
                            #print uniqueones.keys()
                            avgscores={}
                            for paramstr in uniqueones.keys():
                                if len(uniqueones[paramstr])<minresultcount:
                                   allscores[filepathkey][traceparam][inferalgo][sampleinterval][fraction][noiseratio][paramstr]="-1"
                                   continue
                                tp=0.0
                                fp=0.0
                                fn=0.0
                                tn=0.0
                                shufflelist=list(uniqueones[paramstr])
                                random.shuffle(shufflelist)
                                if useallresults:
                                   usefulindices=list(shufflelist)
                                else:    
                                   usefulindices=shufflelist[0:minresultcount]
                                #initial=True
                                #realusefulindices=set()
                                totalcount=1000.0
                                for usefulindex in usefulindices:
                                    outfilename="{0}_{1}".format(paramstr,usefulindex)
                                    resultpath="{0}/{1}".format(fractionpath,outfilename)
                                    print "here index {0}".format(usefulindex)
                                    print resultpath
                                    infile=open(resultpath,"rb")
                                    tempscores=cPickle.load(infile)
                                    infile.close()
                                    #if initial:
                                    #   for score in tempscores.keys():
                                    #       avgscores[score]=0.0
                                    #   initial=False
                                    #idenfity special no prediction made case
                                    sen=tempscores["sen"]
                                    fpr=tempscores["fpr"]
                                    acc=tempscores["acc"]
                                    localtp,localfp,localfn,localtn=returncounts(sen,fpr,acc,totalcount)
                                    tp+=localtp
                                    fp+=localfp
                                    fn+=localfn
                                    tn+=localtn
                                    #if tempscores["precision"]==0.0 and tempscores["recall"]==0.0 and tempscores["spec"]!=0.0 and tempscores["acc"]!=0.0 and tempscores["f1"]==0.0 and tempscores["f1"]==0.0 and tempscores["mcc"]==0.0:
                                    #   continue
                                    #realusefulindices.add(usefulindex) 
                                    #for score in tempscores.keys():
                                    #    avgscores[score]+=tempscores[score]
                                assert round(tp+fp+fn+tn,4)==round(len(usefulindices)*totalcount,4)
                                avgscores=dict(assignscores(tp,fp,fn,tn))
                                if avgscores["precision"]==0.0 and avgscores["recall"]==0.0 and avgscores["spec"]!=0.0 and avgscores["acc"]!=0.0 and avgscores["f1"]==0.0 and avgscores["f1"]==0.0 and avgscores["mcc"]==0.0:
                                   allscores[filepathkey][traceparam][inferalgo][sampleinterval][fraction][noiseratio][paramstr]="-1" 
                                   continue
                                
                                #if len(realusefulindices)==0:
                                #   allscores[filepathkey][traceparam][inferalgo][sampleinterval][fraction][noiseratio][paramstr]="-1"
                                #   continue
                                #for score in avgscores.keys():
                                #    avgscores[score]/=float(len(realusefulindices))
                                
                                #also estimate other f scores
                                tempprecision=avgscores["precision"]
                                temprecall=avgscores["recall"]
                                if tempprecision==0.0 and temprecall==0.0:
                                   avgscores["f1div100"]=0.0
                                   avgscores["f1div10"]=0.0
                                   avgscores["f1div2"]=0.0
                                   avgscores["f10"]=0.0
                                   avgscores["f20"]=0.0
                                else:   
                                   coef=0.01
                                   f1div100=((1.0+(coef**2))*tempprecision*temprecall)/float(((coef**2)*tempprecision)+temprecall)
                                   avgscores["f1div100"]=f1div100
                                   coef=0.1
                                   f1div10=((1.0+(coef**2))*tempprecision*temprecall)/float(((coef**2)*tempprecision)+temprecall)
                                   avgscores["f1div10"]=f1div10
                                   coef=0.5
                                   f1div2=((1.0+(coef**2))*tempprecision*temprecall)/float(((coef**2)*tempprecision)+temprecall)
                                   avgscores["f1div2"]=f1div2
                                   coef=10.0
                                   f10=((1.0+(coef**2))*tempprecision*temprecall)/float(((coef**2)*tempprecision)+temprecall)
                                   avgscores["f10"]=f10
                                   coef=20.0
                                   f20=((1.0+(coef**2))*tempprecision*temprecall)/float(((coef**2)*tempprecision)+temprecall)
                                   avgscores["f20"]=f20
                                allscores[filepathkey][traceparam][inferalgo][sampleinterval][fraction][noiseratio][paramstr]=dict(avgscores)
    return allscores


#graph plot name assigner  
graph2plotname={}
plotfolderstatic="offlineplots"
resultfolderstatic="offlineresults"
tracefolderstatic="traces"
#datatype=[("real","dynamic"),("real","static"),("synthetic","graphwise"),("synthetic,""paramwise",),("synthetic","dynamic")]
datatype=("synthetic","graphwise")
spreadmodel="si"
infertype="edge" #difprobpartial","difprobunknown","spreadprobunknown,"macropartial","lag"
probtype="discrete"
#plotmarkers=["s","o","+","*",".","D","x","s","^","v"]
#plotcolors=["r","b","g","c","k","m","y"]
plotmarkers=["s","o","+","*",".","D","x","s","^","v","s","o","+","*",".","D","x","s","^","v","s","o","+","*",".","D","x","s","^","v"]
plotcolors=["r","b","g","c","k","m","y","r","b","g","c","k","m","y","r","b","g","c","k","m","y"]
#requiredfractions=[0.01,0.02,0.03,0.04,0.05,0.07,0.1,0.15,0.2] #for dynamic graphs
#requiredfractions=[0.01,0.02,0.03,0.04,0.05,0.07,0.1,0.15,0.2]
requiredfractions=[0.01,0.02,0.03,0.04,0.05,0.07,0.1,0.15,0.2] #0.25,0.3,0.35,0.4,0.5,0.6,0.8 
requirededgethres=[0.2,0.3,0.6,0.8,1.2,3.0,5.0,10.0,20.0,30.0,40.0]
#[0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0,1.2,1.4,1.7,2.0,2.5,3.0,3.5,4.0,4.5,10.0,20.0,40.0]
requiredlambdathres=[0.01,0.05,0.1,0.5,1.0,1.5,5.0,10.0,15.0,30.0,50.0,100.0]
#requiredlambdathres=[1.0]
minresultcount=3
useallresults=False #if false use onlt minresultcount number of results
takeaveragesyn=False #estimate average for synthetic graphs
nodecount=500

if __name__ == "__main__":
    realdata=datatype[0]
    graphevolution=datatype[1]
    if realdata=="real" and graphevolution=="static":
       graphfolderprefix="../staticrealdata"
       #tracefolderprefix="../tracegenerator"
       resultfolderprefix="../runner"
    elif realdata=="real" and graphevolution=="dynamic":
       graphfolderprefix="../dynamicrealdata"
       #tracefolderprefix="../tracegenerator"
       resultfolderprefix="../runner"
    elif realdata=="synthetic" and graphevolution=="graphwise":
       graphfolderprefix="../syndatagen"
       #tracefolderprefix="../tracegenerator"
       resultfolderprefix="../runner"
    elif realdata=="synthetic" and graphevolution=="paramwise":
       graphfolderprefix="../syndatagen"
       #tracefolderprefix="../tracegenerator"
       resultfolderprefix="../runner"
    elif realdata=="synthetic" and graphevolution=="dynamic":
       graphfolderprefix="../dynamicsyndatagen"
       #tracefolderprefix="../tracegenerator"
       resultfolderprefix="../runner"

    if realdata=="synthetic":
       graphfolder="{0}_{1}_{2}".format(realdata,graphevolution,"graphs")
       graphfolder="{0}/{1}".format(graphfolderprefix,graphfolder)
    elif realdata=="real":
       graphfolder="{0}".format(graphfolderprefix)
       
    #tracefolder="{0}_{1}_{2}_{3}_{4}".format(tracefolderstatic,realdata,graphevolution,spreadmodel,infertype)
    #tracefolder="{0}/{1}".format(tracefolderprefix,tracefolder)
    resultfolder="{0}_{1}_{2}_{3}_{4}_{5}".format(resultfolderstatic,realdata,graphevolution,spreadmodel,infertype,probtype)
    resultfolder="{0}/{1}".format(resultfolderprefix,resultfolder)
    plotfolder="{0}_{1}_{2}_{3}_{4}_{5}".format(plotfolderstatic,realdata,graphevolution,spreadmodel,infertype,probtype)
    if not os.path.exists(plotfolder):
       os.makedirs(plotfolder)
       
    infofilepaths={} #no need to store G, can retrieve from filepath info
    if realdata=="real" and graphevolution=="static":   
       filenames=myutil.listdirectories(resultfolder)
       for filename in filenames:
           resultfilepath="{0}/{1}".format(resultfolder,filename)
           graphfilepath="{0}/{1}".format(graphfolder,filename)
           filepathkey=filename.replace(".gml","")
           infofilepaths[filepathkey]=[resultfilepath,graphfilepath]
    elif realdata=="real" and graphevolution=="dynamic":
       graphdirnames=myutil.listdirectories(resultfolder)
       for graphdirname in graphdirnames:
           resultfilepath="{0}/{1}".format(resultfolder,graphdirname)
           graphfilepath="{0}/{1}".format(graphfolder,graphdirname)
           filepathkey=graphdirname.replace(".gml","")
           infofilepaths[filepathkey]=[resultfilepath,graphfilepath]
    elif realdata=="synthetic" and graphevolution=="graphwise":
       datas=myutil.listdirectories(resultfolder)
       for datainfo in datas:
           datadirname="{0}/{1}".format(resultfolder,datainfo)
           graphalgos=myutil.listdirectories(datadirname)
           for graphalgo in graphalgos:
               graphdirname="{0}/{1}".format(datadirname,graphalgo)
               filenames=myutil.listdirectories(graphdirname) 
               for filename in filenames:
                   resultfilepath="{0}/{1}/{2}/{3}".format(resultfolder,datainfo,graphalgo,filename)
                   graphfilepath="{0}/{1}/{2}/{3}".format(graphfolder,datainfo,graphalgo,filename)
                   filepathkey="{0}/{1}/{2}".format(datainfo,graphalgo,filename.replace(".pkl",""))
                   infofilepaths[filepathkey]=[resultfilepath,graphfilepath]
    elif realdata=="synthetic" and graphevolution=="dynamic":
       datas=myutil.listdirectories(resultfolder)
       for datainfo in datas:
           datadirname="{0}/{1}".format(resultfolder,datainfo)
           graphalgos=myutil.listdirectories(datadirname)
           for graphalgo in graphalgos:
               graphdirname="{0}/{1}".format(datadirname,graphalgo)
               paramfilenames=myutil.listdirectories(graphdirname)
               for paraminfo in paramfilenames:
                   paramdirname="{0}/{1}".format(graphdirname,paraminfo)
                   filenames=myutil.listdirectories(paramdirname)  
                   for filename in filenames:
                       resultfilepath="{0}/{1}/{2}/{3}/{4}".format(resultfolder,datainfo,graphalgo,paraminfo,filename)
                       graphfilepath="{0}/{1}/{2}/{3}/{4}".format(graphfolder,datainfo,graphalgo,paraminfo,filename)
                       filepathkey="{0}/{1}/{2}/{3}".format(datainfo,graphalgo,paraminfo,filename.replace(".pkl",""))
                       infofilepaths[filepathkey]=[resultfilepath,graphfilepath]
                   
    print infofilepaths.keys()[0:10]
    print len(infofilepaths.keys())
    print infofilepaths.values()[0:10]
    print "started reading"

    #depending on evolution model and data being real or synthetic, we take average here!!
    if takeaveragesyn:
       allscores=dict(estimateaverage(allscores))
    
    #specific params plot generator      
    requiredplots=requiredplotreturner(realdata,graphevolution,spreadmodel,infertype)
    for plotinfo in requiredplots:
        print "starting:"
        print plotinfo
        mainplot,plottype,plotparams=plotinfo
        plotmethodname="{0}plotgenerate".format(mainplot)
        if mainplot not in ["overtime"]: #overtime","parameterauccurve","sampling_fixedtrace","sampling_fixedtrace","distinfo_fixedtrace","noise_fixedtrace"
           continue
        if mainplot=="paramchange": #not being used right now!
           [spreamodel,distinfo,sampleinterval,runalgos,graphname,realdata,graphevolution,noiselevel,fixedfraction]=plotparams
           allscores=returnscores(infofilepaths,[distinfo],[sampleinterval],runalgos,[noiselevel],[graphname],[fixedfraction])
           sentscores=estimateparamwisescores(spreamodel,distinfo,sampleinterval,runalgos,graphname,realdata,graphevolution,noiselevel,fixedfraction,allscores)
           print "paramchange: info pre plot"
           for sentalgo in sentscores.keys():
               print sentalgo
               for sentparamstr in sentscores[sentalgo].keys():
                   print sentparamstr
                   for thres in sentscores[sentalgo][sentparamstr].keys():
                       print thres
           realplotfolder="{0}/{1}/{2}/{3}/{4}".format(plotfolder,graphname,mainplot,plottype,distinfo)
           if not os.path.exists(realplotfolder):
              os.makedirs(realplotfolder)
           plotfilename="{0}-{1}-noise{2}-sample{3}-frac{4}.png".format(mainplot,plottype,noiselevel,sampleinterval,fixedfraction)
           plotfilepath="{0}/{1}".format(realplotfolder,plotfilename)
           plottitle=""
           param=[sentscores,plottype,plottitle,plotfilepath]
           globals()[plotmethodname](*param)
        elif mainplot=="parameterauccurve": #this will be for fixed fraction!!
           [spreamodel,distinfo,sampleinterval,runalgos,graphname,realdata,graphevolution,noiselevel,fixedfraction]=plotparams
           allscores=returnscores(infofilepaths,[distinfo],[sampleinterval],runalgos,[noiselevel],[graphname],[fixedfraction])
           sentscores=estimateparamwisescores(spreamodel,distinfo,sampleinterval,runalgos,graphname,realdata,graphevolution,noiselevel,fixedfraction,allscores)
           print "info pre plot"
           for sentalgo in sentscores.keys():
               print sentalgo
               for sentparamstr in sentscores[sentalgo].keys():
                   print sentparamstr
                   for thres in sentscores[sentalgo][sentparamstr].keys():
                       print thres
           print "Creating AUC plots for algos: {0}".format(sentscores.keys())      
           realplotfolder="{0}/{1}/{2}/{3}/{4}".format(plotfolder,graphname,mainplot,plottype,distinfo)
           if not os.path.exists(realplotfolder):
              os.makedirs(realplotfolder)
           plotfilename="{0}-{1}-noise{2}-sample{3}-frac{4}.png".format(mainplot,plottype,noiselevel,sampleinterval,fixedfraction)
           plotfilepath="{0}/{1}".format(realplotfolder,plotfilename)
           #graphinfo=graphname
           #plottitle="{0} plot for {1} interval:{2}, noise level {3} fraction:{4}".format(plottype,graphinfo,sampleinterval,noiselevel,fixedfraction)
           plottitle=""
           param=[sentscores,plottype,plottitle,plotfilepath]
           globals()[plotmethodname](*param)
        elif mainplot=="overtime":
           print plotparams
           [spreamodel,distinfo,sampleinterval,runalgos,graphname,realdata,graphevolution,noiselevel,algoparaminfo]=plotparams
           allscores=returnscores(infofilepaths,[distinfo],[sampleinterval],runalgos,[noiselevel],[graphname],requiredfractions,algoparaminfo)
           sentscores=estimatefractionwisescores(spreamodel,distinfo,sampleinterval,runalgos,graphname,realdata,graphevolution,noiselevel,allscores)
           #maybe eliminate some here!!
           print "Creating plots for algos: {0}".format(sentscores.keys())      
           realplotfolder="{0}/{1}/{2}/{3}/{4}".format(plotfolder,graphname,mainplot,plottype,distinfo)
           if not os.path.exists(realplotfolder):
              os.makedirs(realplotfolder)
           plotfilename="{0}-{1}-noise{2}-sample{3}.png".format(mainplot,plottype,noiselevel,sampleinterval)
           plotfilepath="{0}/{1}".format(realplotfolder,plotfilename)
           #graphinfo=graphname
           #plottitle="{0} plot for {1} interval:{2}, noise level {3}".format(plottype,graphinfo,sampleinterval,noiselevel)
           plottitle=""
           xlabel="Trace count"
           xaxis=requiredfractions
           param=[sentscores,plottype,plottitle,plotfilepath,xlabel,algoparaminfo,xaxis]
           globals()[plotmethodname](*param)
        elif mainplot in ["noise_fixedtrace","sampling_fixedtrace","distinfo_fixedtrace"]:
           #assert plottype in ["f1","f1overn","f1div100","f1div10","f1div2","f20","f50","f100","roc","pr","bep"]
           if mainplot=="noise_fixedtrace":
              print "here noise plot" 
              [spreadmodel,distinfo,sampleinterval,runalgos,graphname,realdata,graphevolution,noiselevels,fixedfraction,algoparaminfo]=plotparams
              allscores=returnscores(infofilepaths,[distinfo],[sampleinterval],runalgos,noiselevels,[graphname],[fixedfraction],algoparaminfo)
              xaxis=noiselevels
           elif mainplot=="sampling_fixedtrace":
              [spreadmodel,distinfo,sampleintervals,runalgos,graphname,realdata,graphevolution,noiselevel,fixedfraction,algoparaminfo]=plotparams
              allscores=returnscores(infofilepaths,[distinfo],sampleintervals,runalgos,[noiselevel],[graphname],[fixedfraction],algoparaminfo)
              xaxis=sampleintervals
           elif mainplot=="distinfo_fixedtrace":
              print "inside probdist" 
              [spreadmodel,distinfos,sampleinterval,runalgos,graphname,realdata,graphevolution,noiselevel,fixedfraction,algoparaminfo]=plotparams
              allscores=returnscores(infofilepaths,distinfos,[sampleinterval],runalgos,[noiselevel],[graphname],[fixedfraction],algoparaminfo)
              #xaxis=distinfos
              xaxis=[]
              #validity check
              splitted=distinfos[0].split("-")
              spreadprob=float(splitted[0].replace("spreadprob_",""))
              s2idist=splitted[1].replace("s2i_","").split("_")[0]
              s2idistparam=float(splitted[1].replace("s2i_","").split("_")[1])
              xaxis.append(s2idistparam)
              distparammap={}
              distparammap[distinfos[0]]=s2idistparam
              for distinfo in distinfos[1:]:
                  tempsplitted=distinfo.split("-")
                  tempspreadprob=float(tempsplitted[0].replace("spreadprob_",""))
                  temps2idist=tempsplitted[1].replace("s2i_","").split("_")[0]
                  temps2idistparam=float(tempsplitted[1].replace("s2i_","").split("_")[1])
                  xaxis.append(temps2idistparam)
                  distparammap[distinfo]=temps2idistparam
                  assert tempspreadprob==spreadprob and temps2idist==s2idist
           sentscores={}
           filepathkey=graphname
           for runalgo in runalgos:
               if runalgo.find("sum")!=-1:
                  extend="sum"
                  sentalgo=runalgo.replace("sum","")
               elif runalgo.find("multi")==-1 and runalgo.find("mul")!=-1:
                  extend="mul"
                  sentalgo=runalgo.replace("mul","")
               else: 
                  extend=""
                  sentalgo=runalgo
               if not sentscores.has_key(sentalgo):
                  sentscores[sentalgo]={}

               if mainplot=="noise_fixedtrace":
                  if not allscores[filepathkey].has_key(distinfo):
                     continue
                  if not allscores[filepathkey][distinfo][runalgo][sampleinterval].has_key(fixedfraction):
                     continue
                  for noiselevel in noiselevels:
                      if not allscores[filepathkey][distinfo][runalgo][sampleinterval][fixedfraction].has_key(noiselevel):
                         continue 
                      for paramstr in allscores[filepathkey][distinfo][runalgo][sampleinterval][fixedfraction][noiselevel].keys():
                          print paramstr
                          tempparts=paramstr.split("_")
                          temp3=tempparts[3].replace("-randomnoise","")
                          delnoise=float(temp3)
                          print delnoise
                          tempparts[3]=tempparts[3].replace(str(delnoise),"")
                          curparamstr="_".join(tempparts)
                          sentparamstr="{0}+{1}".format(curparamstr,extend)
                          if not sentscores[sentalgo].has_key(sentparamstr):
                             sentscores[sentalgo][sentparamstr]={}
                          assert not sentscores[sentalgo][sentparamstr].has_key(noiselevel)
                          sentscores[sentalgo][sentparamstr][noiselevel]=allscores[filepathkey][distinfo][runalgo][sampleinterval][fixedfraction][noiselevel][paramstr]
               elif mainplot=="sampling_fixedtrace":
                  if not allscores[filepathkey].has_key(distinfo):
                     continue
                  for sampleinterval in sampleintervals:
                      if not allscores[filepathkey][distinfo][runalgo].has_key(sampleinterval):
                         continue
                      print runalgo
                      print filepathkey
                      print sampleinterval
                      print fixedfraction
                      if not allscores[filepathkey][distinfo][runalgo][sampleinterval].has_key(fixedfraction):
                         continue
                      if not allscores[filepathkey][distinfo][runalgo][sampleinterval][fixedfraction].has_key(noiselevel):
                         continue
                      for paramstr in allscores[filepathkey][distinfo][runalgo][sampleinterval][fixedfraction][noiselevel].keys():
                          sentparamstr="{0}+{1}".format(paramstr,extend)
                          if not sentscores[sentalgo].has_key(sentparamstr):
                             sentscores[sentalgo][sentparamstr]={}
                          assert not sentscores[sentalgo][sentparamstr].has_key(sampleinterval)
                          sentscores[sentalgo][sentparamstr][sampleinterval]=allscores[filepathkey][distinfo][runalgo][sampleinterval][fixedfraction][noiselevel][paramstr]
               elif mainplot=="distinfo_fixedtrace":
                  for distinfo in distinfos:
                      if not allscores[filepathkey].has_key(distinfo):
                         continue
                      if not allscores[filepathkey][distinfo][runalgo][sampleinterval].has_key(fixedfraction):
                         continue
                      if not allscores[filepathkey][distinfo][runalgo][sampleinterval][fixedfraction].has_key(noiselevel):
                         continue 
                      plotdistinfo=distparammap[distinfo]
                      for paramstr in allscores[filepathkey][distinfo][runalgo][sampleinterval][fixedfraction][noiselevel].keys():
                          sentparamstr="{0}+{1}".format(paramstr,extend)   
                          if not sentscores[sentalgo].has_key(sentparamstr):
                             sentscores[sentalgo][sentparamstr]={}
                          assert not sentscores[sentalgo][sentparamstr].has_key(plotdistinfo)   
                          sentscores[sentalgo][sentparamstr][plotdistinfo]=allscores[filepathkey][distinfo][runalgo][sampleinterval][fixedfraction][noiselevel][paramstr]
           #eliminate some results (we will only check fraction validity)
           #varparam can be either probdist,sampleinterval,noiselevel
           for sentalgo in sentscores.keys():
               for sentparamstr in sentscores[sentalgo].keys():
                   for varparam in sentscores[sentalgo][sentparamstr].keys():
                       if sentscores[sentalgo][sentparamstr][varparam]=="-1":
                          del sentscores[sentalgo][sentparamstr][varparam]
           for sentalgo in sentscores.keys():
               for sentparamstr in sentscores[sentalgo].keys():
                   if len(sentscores[sentalgo][sentparamstr])==0:
                      del sentscores[sentalgo][sentparamstr]
               if len(sentscores[sentalgo])==0:
                  del sentscores[sentalgo]
           print "Creating plots for algos: {0}".format(sentscores.keys())
           if mainplot=="distinfo_fixedtrace":
              realplotfolder="{0}/{1}/{2}/{3}/{4}-{5}".format(plotfolder,graphname,mainplot,plottype,spreadprob,s2idist)
              plotfilename="{0}-{1}-{2}-{3}-{4}.png".format(mainplot,plottype,noiselevel,sampleinterval,fixedfraction)
              #plottitle="{0} plot for {1} interval:{2}, noise level {3}, fixedfraction {4}".format(plottype,graphinfo,sampleinterval,noiselevel,fixedfraction)
              xlabel="{0} distribution parameters".format(s2idist)
           elif mainplot=="sampling_fixedtrace":
              realplotfolder="{0}/{1}/{2}/{3}/{4}".format(plotfolder,graphname,mainplot,plottype,distinfo)
              plotfilename="{0}-{1}-{2}-{3}.png".format(mainplot,plottype,noiselevel,fixedfraction)
              #plottitle="{0} plot for {1}, noise level {2}, fixedfraction:{3}".format(plottype,graphinfo,noiselevel,fixedfraction)
              xlabel="Sampling intervals"
           elif mainplot=="noise_fixedtrace":
              realplotfolder="{0}/{1}/{2}/{3}/{4}".format(plotfolder,graphname,mainplot,plottype,distinfo)
              plotfilename="{0}-{1}-{2}-{3}.png".format(mainplot,plottype,sampleinterval,fixedfraction)
              #plottitle="{0} plot for {1} interval:{2}, fixedfraction:{3}".format(plottype,graphinfo,sampleinterval,fixedfraction)
              xlabel="Partial observability(Uncertainity) degree"
           plottitle=""
           if not os.path.exists(realplotfolder):
              os.makedirs(realplotfolder)   
           plotfilepath="{0}/{1}".format(realplotfolder,plotfilename)
           graphinfo=graphname
           #graphinfo=graphinfo.replace("126","512")
           #probdistplotgenerate(sentscores2,plottype,plottitle,plotfilepath,xlabel)
           overtimeplotgenerate(sentscores,plottype,plottitle,plotfilepath,xlabel,algoparaminfo,xaxis)
           #param=[sentscores,plottype,plottitle,plotfilepath,xlabel,algoparaminfo]
           #globals()[plotmethodname](*param)

        elif mainplot in ["noise_variabletrace","sampling_variabletrace","probdist_variabletrace"]: #variable trace can only happen for early retrieval right now
           assert plottype in [(("early","1-square"),"f1"),(("early","1-square"),"f1overn"),(("early","1-square"),"f1div100"),(("early","1-square"),"f1div10"),(("early","1-square"),"f1div2"),(("early","1-square"),"f20"),(("early","1-square"),"f50"),(("early","1-square"),"f100"),(("early","1-square"),"roc"),(("early","1-square"),"pr"),(("early","1-square"),"bep"),(("early","1-normal"),"f1"),(("early","1-normal"),"f1overn"),(("early","1-normal"),"f1div100"),(("early","1-normal"),"f1div10"),(("early","1-normal"),"f1div2"),(("early","1-normal"),"f20"),(("early","1-normal"),"f50"),(("early","1-normal"),"f100"),(("early","1-normal"),"roc"),(("early","1-normal"),"pr"),(("early","1-normal"),"bep")]
           if mainplot=="noisy_variabletrace":
              [spreadmodel,distinfo,sampleinterval,runalgos,graphname,realdata,graphevolution,noiselevels]=plotparams
           elif mainplot=="sampling_variabletrace":
              [spreadmodel,distinfo,sampleintervals,runalgos,graphname,realdata,graphevolution,noiselevel]=plotparams
           elif mainplot=="probdist_variabletrace":
              [spreadmodel,distinfos,sampleinterval,runalgos,graphname,realdata,graphevolution,noiselevel]=plotparams
              #validity check
              splitted=distinfos[0].split("-")
              spreadprob=float(splitted[0].replace("spreadprob_",""))
              s2idist=splitted[1].replace("s2i_","").split("_")[0]
              s2idistparam=float(splitted[1].replace("s2i_","").split("_")[1])
              distparammap={}
              distparammap[distinfos[0]]=s2idistparam
              for distinfo in distinfos[1:]:
                  tempsplitted=distinfo.split("-")
                  tempspreadprob=float(tempsplitted[0].replace("spreadprob_",""))
                  temps2idist=tempsplitted[1].replace("s2i_","").split("_")[0]
                  temps2idistparam=float(tempsplitted[1].replace("s2i_","").split("_")[1])
                  distparammap[distinfo]=temps2idistparam
                  assert tempspreadprob==spreadprob and temps2idist==s2idist
           earlytype=plottype[0]
           scoretype=plottype[1]    
           sentscores={}
           filepathkey=graphname
           for runalgo in runalgos:
               if runalgo.find("sum")!=-1:
                  extend="sum"
                  sentalgo=runalgo.replace("sum","")
               elif runalgo.find("multi")==-1 and runalgo.find("mul")!=-1:
                  extend="mul"
                  sentalgo=runalgo.replace("mul","")
               else: 
                  extend=""
                  sentalgo=runalgo
               if not sentscores.has_key(sentalgo):
                  sentscores[sentalgo]={}

               if mainplot=="noisy_variabletrace":
                  if not allscores[filepathkey].has_key(distinfo):
                    continue
                  for noiselevel in noiselevels:
                      for fraction in allscores[filepathkey][distinfo][runalgo][sampleinterval].keys():
                          if fraction not in requiredfractions:
                             continue
                          if not allscores[filepathkey][distinfo][runalgo][sampleinterval][fraction].has_key(noiselevel):
                             continue
                          for paramstr in allscores[filepathkey][distinfo][runalgo][sampleinterval][fraction][noiselevel].keys():
                              sentparamstr="{0}+{1}".format(paramstr,extend)
                              if not sentscores[sentalgo].has_key(sentparamstr):
                                 sentscores[sentalgo][sentparamstr]={}
                              if not sentscores[sentalgo][sentparamstr].has_key(noiselevel):
                                 sentscores[sentalgo][sentparamstr][noiselevel]={}
                              assert not sentscores[sentalgo][sentparamstr][noiselevel].has_key(fraction)
                              sentscores[sentalgo][sentparamstr][noiselevel][fraction]=allscores[filepathkey][distinfo][runalgo][sampleinterval][fraction][noiselevel][paramstr]
               elif mainplot=="sampling_variabletrace":
                  if not allscores[filepathkey].has_key(distinfo):
                     continue
                  for sampleinterval in sampleintervals:
                      assert allscores[filepathkey][distinfo][runalgo].has_key(sampleinterval)
                      for fraction in allscores[filepathkey][distinfo][runalgo][sampleinterval].keys():
                          if fraction not in requiredfractions:
                             continue
                          if not allscores[filepathkey][distinfo][runalgo][sampleinterval][fraction].has_key(noiselevel):
                             continue
                          for paramstr in allscores[filepathkey][distinfo][runalgo][sampleinterval][fraction][noiselevel].keys():
                              sentparamstr="{0}+{1}".format(paramstr,extend)
                              if not sentscores[sentalgo].has_key(sentparamstr):
                                 sentscores[sentalgo][sentparamstr]={}
                              if not sentscores[sentalgo][sentparamstr].has_key(sampleinterval):
                                 sentscores[sentalgo][sentparamstr][sampleinterval]={}
                              assert not sentscores[sentalgo][sentparamstr][sampleinterval].has_key(fraction)
                              sentscores[sentalgo][sentparamstr][sampleinterval][fraction]=allscores[filepathkey][distinfo][runalgo][sampleinterval][fraction][noiselevel][paramstr] 
               elif mainplot=="distprob_variabletrace":
                  for distinfo in distinfos:
                      if not allscores[filepathkey].has_key(distinfo):
                         continue 
                      plotdistinfo=distparammap[distinfo]
                      for fraction in allscores[filepathkey][distinfo][runalgo][sampleinterval].keys():
                          if fraction not in requiredfractions:
                             continue
                          if not allscores[filepathkey][distinfo][runalgo][sampleinterval][fraction].has_key(noiselevel):
                             continue
                          for paramstr in allscores[filepathkey][distinfo][runalgo][sampleinterval][fraction][noiselevel].keys():
                              sentparamstr="{0}+{1}".format(paramstr,extend)   
                              if not sentscores[sentalgo].has_key(sentparamstr):
                                 sentscores[sentalgo][sentparamstr]={}
                              if not sentscores[sentalgo][sentparamstr].has_key(plotdistinfo):
                                 sentscores[sentalgo][sentparamstr][plotdistinfo]={}
                              assert not sentscores[sentalgo][sentparamstr][plotdistinfo].has_key(fraction)   
                              sentscores[sentalgo][sentparamstr][plotdistinfo][fraction]=allscores[filepathkey][distinfo][runalgo][sampleinterval][fraction][noiselevel][paramstr]
           #eliminate some results (we will only check fraction validity)
           #eliminate algorithms that does not have result for all distinfos given
           for sentalgo in sentscores.keys():
               for sentparamstr in sentscores[sentalgo].keys():
                   for varparam in sentscores[sentalgo][sentparamstr].keys(): 
                       tempfractions=set(sentscores[sentalgo][sentparamstr][varparam].keys())
                       removeitems=set()
                       for frac in tempfractions:
                           if sentscores[sentalgo][sentparamstr][varparam][frac]=="-1":
                              removeitems.add(frac)
                       tempfractions-=removeitems 
                       if len(set(requiredfractions).difference(tempfractions))!=0:
                          del sentscores[sentalgo][sentparamstr][varparam]
           for sentalgo in sentscores.keys():
               for sentparamstr in sentscores[sentalgo].keys(): 
                   for varparam in sentscores[sentalgo][sentparamstr].keys():
                       if len(sentscores[sentalgo][sentparamstr][varparam])==0:
                          del sentscores[sentalgo][sentparamstr][varparam] 
                   if len(sentscores[sentalgo][sentparamstr])==0:
                      del sentscores[sentalgo][sentparamstr]
               if len(sentscores[sentalgo])==0:
                  del sentscores[sentalgo]

           print "Creating plots for algos: {0}".format(sentscores.keys())  
           if mainplot=="probdist_variabletrace":
              realplotfolder="{0}/{1}/{2}/{3}/{4}-{5}".format(plotfolder,graphname,mainplot,plottype,spreadprob,s2idist)
              plotfilename="{0}-{1}-{2}-{3}.png".format(mainplot,plottype,noiselevel,sampleinterval)
              plottitle="{0} plot for {1} interval:{2}, noise level {3}".format(plottype,graphinfo,sampleinterval,noiselevel)
              xlabel="{0} distribution parameters".format(s2idist)
           elif mainplot=="sampling_variabletrace":
              realplotfolder="{0}/{1}/{2}/{3}/{4}".format(plotfolder,graphname,mainplot,plottype,distinfo)
              plotfilename="{0}-{1}-{2}.png".format(mainplot,plottype,noiselevel)
              plottitle="{0} plot for {1}, noise level {2}".format(plottype,graphinfo,noiselevel)
              xlabel="Sampling intervals"
           elif mainplot=="noisy_variabletrace":
              realplotfolder="{0}/{1}/{2}/{3}/{4}".format(plotfolder,graphname,mainplot,plottype,distinfo)
              plotfilename="{0}-{1}-{2}.png".format(mainplot,plottype,sampleinterval)
              plottitle="{0} plot for {1} interval:{2}".format(plottype,graphinfo,sampleinterval)
              xlabel="Noise levels"

           plotfilepath="{0}/{1}".format(realplotfolder,plotfilename)
           graphinfo=graphname
           #graphinfo=graphinfo.replace("126","512")
           #noiseplotgenerate(sentscores,plottype,plottitle,plotfilepath,xlabel)
           param=[sentscores,plottype,plottitle,plotfilepath,xlabel]
           globals()[plotmethodname](*param)

        else:
            print "plottype {0} is unknown".format(mainplot)
            exit(1)
                  
