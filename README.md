# BiHiDef
Bipartite adaptation of HiDeF using BRIM. The goal of this program is to adapt the HiDef (cite paper/package) to work with bipartite networks and BRIM instead of unipartite networks and Leiden/Louvain. This package is dependent on the HiDef package (it uses its classes and the weaver module) and the netZooPy package for the pycondor implementation of BRIM.

HiDef creates a proximity graph from a range of resolutions and runs Leiden on each of these resolution. Communities in different proximal resolutions are merged if they have a Jaccard score over 0.75 to obtain persistent communities across resolutions. The persistent communities are finally organised into a hierarchy through containment indices. Therefore HiDef obtains hierarchical overlapping communities obtained from different scales.

BiHiDef uses the same approach to bipartite networks. Running BRIM on a bipartite network gives a community structure on each of the nodesets. On a given resolution these community structures are highly dependent on each other. We apply BRIM along the range of resolutions, but the community structures on each nodeset are kept track of sepparately. The merging of communities and hierarchy are constructed independently for each nodeset. Therefore BiHiDef obtains two hierarchical overlapping community structures on a bipartite network. 

## Install instructions
Warning: The package requires a version of NetZooPy not yet available. Hopefully soon it will be merged.
```
git clone https://github.com/genisott/BiHiDef.git
cd bihidef
pip install .
```
## Usage
```
import bihidef
bihidef("edgelist.txt")
```
The file should be a comma separated edgelist without header. Two or three columns. If no weights are given they will be initialized as 1.
