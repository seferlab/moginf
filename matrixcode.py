


def trace2Matrix(algoinfo,algopar,traces,infernodes,smodel,dists,sideparams,iprobmethod):
    """Generates constraint matrix out of trace data
    Args:
       algo: infer method
       algopar: method parameters
       alltraces: all trace data
       smodel: diffusion model
       infernodes: nodes which will be inferred(only important when run in parallel mode)
       dists: distributions
       sideparams: sideparams for CEFER execution(EPSILON,LOGBASE)
       iprobmethod: infection probability assignment method
    Returns:
       consmat: Trace Constraint matrix
       consbmat: Trace cons. right hand side
       covermat: Covering Constraint matrix
       coverbmat: Covering right hand side
       node2loc: maps node 2 their locations(indexing purpose) 
    """

    EPSILON,LOGBASE = sideparams
    node2loc = {infernodes[index]: index for index in xrange(len(infernodes))}
    rownum = 10000000
    colnum = len(infernodes) * len(infernodes)
    consmat = scipy.sparse.csc_matrix((rownum,colnum),dtype = np.float)
    consbmat = scipy.sparse.csc_matrix((rownum,1),dtype = np.float)
    covermat = scipy.sparse.csc_matrix((rownum,colnum),dtype = np.float)
    consnum = 0
    covernum = 0
    consrow, conscol, consdata, coverrow, covercol, coverdata = [],[],[],[],[],[] 
    for trace in traces:
        itimes = set([trace[node][Trace.INFECTED] for node in trace.keys() if trace[node].has_key(Trace.INFECTED)])
        maxitime = max(itimes)
        if smodel in ["si","sir"]:       
           key1 = Trace.INFECTED
        elif smodel == "seir":
           key1 = Trace.EXPOSED
        for node in infernodes:
            if trace.has_key(node):
               assert trace[node].has_key(key1) 
               poseffectors=[]
               for sendernode in trace.keys():
                   if sendernode == node or not trace[sendernode].has_key(INFECTED):
                      continue
                   if trace[sendernode][INFECTED] < trace[node][key1]:
                      if smodel in ["seir","sir"]:
                         if trace[node][key1] >= trace[sendernode][RECOVERED]:
                            timedif = trace[sendernode][RECOVERED] - trace[sendernode][INFECTED]
                         else:
                            timedif = trace[node][key1] - trace[sendernode][INFECTED] 
                      else:
                         timedif = trace[node][key1] - trace[sendernode][INFECTED] 
                      poseffectors.append(sendernode)
                      approx = 0.001
                      pretimedif = timedif-approx
                      diststr = "s2{0}".format(key1)
                      val = 1.0 - tracegen.retCdf(dists[diststr][0],[dists[diststr][1]],timedif,dists["sprob"])
                      preval = 1.0 - tracegen.retCdf(dists[diststr][0],[dists[diststr][1]],pretimedif,dists["sprob"])
                      precoef = math.log(preval,LOGBASE)
                      coef = math.log(val,LOGBASE)
                      coef -= precoef
                      coef /= approx
                      assert coef <= 0.00000000001
                      colpos = (len(infernodes) * node2loc[sendernode]) + node2loc[node]
                      #consmat[consnum,colpos] = coef
                      consrow.append(consnum)
                      conscol.append(colpos)
                      consdata.append(coef)
                   consbmat[consnum,0] = math.log(EPSILON,LOGBASE)
                   consnum += 1
               if len(poseffectors) > 0 and "cover" in algoinfo:
                  for node in poseffectors:
                      colpos = (len(infernodes) * node2loc[sendernode]) + node2loc[node]
                      #covermat[covernum,colpos] = 1.0
                      coverrow.append(covernum)
                      covercol.append(colpos)
                      coverdata.append(1.0)
                  covernum += 1
            else:
               for sendernode in trace.keys():
                   if sendernode == node or not trace[sendernode].has_key(INFECTED):
                      continue
                   timedif = maxitime - trace[sendernode][INFECTED]
                   diststr = "s2{0}".format(key1)
                   val = 1.0-tracegen.retCdf(dists[diststr][0],[dists[diststr][1]],timedif,dists["sprob"])
                   coef = math.log(val,LOGBASE)
                   colpos = (len(infernodes) * node2loc[sendernode]) + node2loc[node]
                   #consmat[consnum,colpos] = coef
                   consrow.append(consnum)
                   conscol.append(colpos)
                   consdata.append(coef)
               consbmat[consnum,0] = 0.0
               consnum += 1
    consmat=scipy.sparse.csc_matrix((consdata,(consrow,conscol)), shape=(consnum,colnum))
    covermat=scipy.sparse.csc_matrix((coverdata,(coverrow,covercol)), shape=(covernum,colnum))
    #consmat = sp.delete(consmat,np.s_[consnum:],0)
    #covermat = sp.delete(covermat,np.s_[consnum:],0)
    return (consmat,consbmat,covermat,node2loc)       
    

def matrix2Lp(consmat,consbmat,covermat,node2loc,algoinfo,algopar):
    """creates lp format program out of given constraint matrix, objective matrix and node mappings
    Args:
       consmat: Trace Constraint matrix
       consbmat: Trace cons. right hand side
       covermat: Covering Constraint matrix
       node2loc: maps node 2 their locations(indexing purpose) 
       algoinfo: algo parameter info
       algopar: algo parameters
    Returns:
       consstr: trace constraint string
       objstr: objective function constraint string
       boundstr: boundary string
    """
    (errortype,sparsetype,fusedtype,cover) = algoinfo
    assert errortype in ["abse","lse"]
    loc2node = dict([[v,k] for k,v in node2loc.items()])

    #trace constraints
    consrownum,conscolnum = sp.shape(consmat)
    consbrownum,temp = sp.shape(consbmat)
    assert consrownum == consbrownum
    constr = ""
    for rowindex in xrange(consrownum):
        linestr = ""
        for colindex in xrange(conscolnum):
            node1 = loc2node[colindex/len(node2loc.keys())]
            node2 = loc2node[colindex%len(node2loc.keys())]
            varname = "x{0}?{1}".format(node1,node2)
            linestr += " {0} {1}".format(consmat[rowindex,colindex],varname)
        plusvar = "p?{0}".format(rowindex)
        minusvar = "m?{0}".format(rowindex)
        linestr += " {0} - {1} = {2} \n ".format(plusvar,minusvar,consbmat[rowindex,0])    
        consstr += linestr

    #set covering constraints
    coverrownum,covercolnum = sp.shape(covermat)
    for rowindex in xrange(coverrownum):
        linestr = ""
        for colindex in xrange(covercolnum):
            node1 = loc2node[colindex / len(node2loc.keys())]
            node2 = loc2node[colindex % len(node2loc.keys())]
            varname = "x{0}?{1}".format(node1,node2)
            linestr += " {0} {1}".format(covermat[rowindex,colindex],varname)
        linestr += " >= 1.0 \n "    
        consstr += linestr

    #objective function string     
    objstr = "Minimize\n obj: "
    if errortype == "abse":
       for rowindex in xrange(consrownum):
           coef = 1.0
           varname = "p?{0}".format(rowindex+1)
           objstr += " + {0} {1}".format(coef,varname)
           varname = "m?{0}".format(rowindex+1)
           objstr += " + {0} {1}".format(coef,varname)
    elif errortype == "lse":
       objstr += " [ " 
       for rowindex in xrange(consrownum):
           coef = 2.0
           varname = "p?{0}".format(rowindex+1)
           objstr += " + {0} {1} * {2}".format(coef,varname,varname)
           varname = "m?{0}".format(rowindex+1)
           objstr += " + {0} {1} * {2}".format(coef,varname,varname)
       objstr += " ] / 2"  
         
    if sparsetype in ["l1","l1l2"]:
       for node1 in node2loc.keys():
           for node2 in node2loc.keys():
               if node1 == node2:
                  continue     
               varname = "x{0}?{1}".format(node1,node2)
               objstr += " + {0} {1}".format(algopar["lambda1"],varname)
    if sparsetype in ["l2","l1l2"]:
       objstr += " [ "  
       for node1 in node2loc.keys():
           for node2 in node2loc.keys():
               if node1 == node2:
                  continue     
               varname = "x{0}?{1}".format(node1,node2)
               objstr += " + {0} {1} * {1}".format(2.0*algopar["lambda2"],varname)
       objstr += " ] / 2"
    
    #boundary constraints
    boundstr = "Bounds\n"
    for node1 in node2loc.keys():
        for node2 in node2loc.keys():
            if node1 == node2:
               continue 
            varname = "x{0}?{1}".format(node1,node2)
            boundstr += "0 <= {0} <= 1\n".format(varname)         
    for rowindex in xrange(consrownum):
        varname = "p?{0}".format(rowindex+1)
        boundstr += "{0} >= 0 \n".format(varname)
        varname = "m?{0}".format(rowindex+1)
        boundstr += "{0} >= 0 \n".format(varname)
    
    return (consstr,objstr,boundstr)
