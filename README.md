# MOGINF
Modular Graph Inference from SEIR Dynamics

Installation:
pip install -r requirements.txt


Run MOGINF:



Graph and Community Datasets:
1- Processed human contact interaction data is on file "human_contact_network.edgelist", whereas
its prior community structure is on file "human_contact_network.community" which comes from spatial
closeness of human while tracking them. Raw datasets are in "rawcontactdata" directory.

2- LFR Benchmarks: We have generated 20 synthetic benchmarks under "lfr_benchmark"
folder. Each has graph as edgelist format, community as .community
file. Community file is tab separated list of nodes for each community.


Diffusion Trace Datasets:
We can generate synthetic diffusion traces over graph by:

