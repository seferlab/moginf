# MOGINF
Modular Graph Inference from SEIR Dynamics

Installation:
pip install -r requirements.txt


Graph and Community Datasets:
1- Processed human contact interaction data is on file "human_contact_network.edgelist", whereas
its prior community structure is on file "human_contact_network.community" which comes from spatial
closeness of human while tracking them. Raw datasets are in "rawcontactdata" directory.

2- LFR Benchmarks: We have generated 20 synthetic benchmarks under "lfr_benchmark"
folder. Each has graph as edgelist format, community as .community
file. Community file is tab separated list of nodes for each community.



In order to run MOGINF:
1- Make sure CPLEX is installed in your system(If CPLEX is not available but you have another solver for LP format, contact us or modify code accordingly)
2- Diffusion data is available in plain text format(You can create
data from our diffusion data generator for various models and
undersampled/noisy cases)

A- Trace Generation:

a- Trace Info:

MOGINF assumes traces are stored in plain text format. When data is
perfect, it is simply,

node1 state time 
node2 state time 
...

When data is noisy, the data is stored in the following format,

node time state probability
...

b- Generating Traces:

In order to generate traces, run as following

python tracegen.py trace.config

trace.config file has the following parameters with default values in
paranthesis

ongraph: whether trace information is provided as graph attributes (False)
traceoutfolder: trace output folder(cefertraces)
graphfilename: graph file (graph.gml)
tracecount: number of traces(10)
samplerate: sampling rate(0)
startnodecount: number of start nodes(1)
noise: noise rate(0.0)
smodel: spreading model(si)
evol: graph evolution(static)
dist: continous(cont) or discrete distribution(dis)
s2i: transition distribution for susceptible to infected state(expo 1.0)
i2r: transition distribution for infected to recovered state(expo 2.0)
e2i: transition distribution for exposed to infected state(expo 2.0)
i2s: transition distribution for infected to susceptible state(expo 2.0)
s2e: transition distribution for susceptible to exposed state(expo 1.0)
sprob: spreading probability(0.1)


B- Graph Inference

MOGINF reads traces and uses either default probability estimation for
second step(default) or Kernel method(Kernel). When data is perfect,
it infers graph with only first step. It also generates png file of
graph returned. For more details about those methods, look at our paper.

Run as following:
python moginfcplex.py moginf.config

moginf.config is configuration file. It includes the following parameters with default values in paranthesis

tracefolder: folder including traces
resultfilename: file for saving output 
tracecount: number of traces to be used(50)
smodel: spreading model (si)
samplerate: sampling rate(0)
evol: type of graph evolution(static)
sprob: spreading probability(0.1)
cover: covering constraints exist or not(cover)
extra:  Algorithm for the second step (Kernel)
errortype: type of error to be optimized(abse)
s2i: transition distribution for susceptible to infected state(expo 1.0)
i2r: transition distribution for infected to recovered state(expo 2.0)
e2i: transition distribution for exposed to infected state(expo 2.0)
i2s: transition distribution for infected to susceptible state(expo 2.0)
s2e: transition distribution for susceptible to exposed state(expo 1.0)
runmode: run mode(parallel or serial)(serial)
parallelcount: number of parallel partitions(200)
graphfilename: original graphfile to estimate score(graph.gml)
printscore: score to be output, can be also be None(f01)

C- Further info:
- MOGINF returns graphs in edgelist format for static case
- Tracefiles can be generated by tracegen.py code associated with this package

If you use this code, please cite our following paper:



