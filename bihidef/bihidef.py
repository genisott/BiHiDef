import numpy as np
import networkx as nx
from netZooPy import condor
import pandas as pd
import networkx as nx
import igraph as ig
import os
import pickle
import copy
from scipy.sparse import *
import scipy as sp
import multiprocessing as mp


from hidef.hidef_finder import ClusterGraph,update_resolution_graph
from hidef.hidef_finder import collapse_cluster_graph,consensus
from hidef.hidef_finder import output_all
from hidef.utils import jaccard_matrix
from hidef import weaver


def create_resolution_graph(minres=0.001,maxres=10,density=0.1,neighbors=10,min_diff_resolution = 0.001):
    
    resolution_graph = nx.Graph()

    stack_res_range = []
    stack_res_range.append((minres, maxres))   

    _ = update_resolution_graph(resolution_graph, minres, density, neighbors)
    _ = update_resolution_graph(resolution_graph, maxres, density, neighbors)

    all_resolutions = [minres, maxres]

    while stack_res_range:
            current_range = stack_res_range.pop(0)
            resname1, resname2 = '{:.4f}'.format(current_range[0]), '{:.4f}'.format(current_range[1])
            # LOGGER.debug('Current resolution range:{} {}'.format(resname1, resname2))

            if round(current_range[1] - current_range[0], 4) <= min_diff_resolution:
                # #LOGGER.debug('Reaching the minimum difference between resolutions')
                continue
            if resolution_graph.nodes[resname1]['padded'] and resolution_graph.nodes[resname2]['padded']:
                continue

            # sample new resolutions and generate more partitions
            new_resolution = np.round(np.sqrt(current_range[1] * current_range[0]), 4)

            stack_res_range.append((current_range[0], new_resolution))
            stack_res_range.append((new_resolution, current_range[1]))

            all_resolutions.append(new_resolution)

            _ = update_resolution_graph(resolution_graph, new_resolution, density, neighbors)
    resolution_graphR = copy.deepcopy(resolution_graph)
    
    return resolution_graph,resolution_graphR,all_resolutions


def run_alg(condor_object,resolution):
    condor_object.initial_community(resolution=resolution)
    condor_object.brim(resolution=resolution)
    clT = sorted(condor_object.tar_memb["community"].unique())
    clR = sorted(condor_object.reg_memb["community"].unique())
    T = sp.sparse.coo_matrix(np.matrix([np.array(condor_object.tar_memb["community"]==i).astype(int) for i in clT])).tocsr()
    R = sp.sparse.coo_matrix(np.matrix([np.array(condor_object.reg_memb["community"]==i).astype(int) for i in clR])).tocsr()
    print("Resolution: "+str(resolution)+" NComs: "+str(len(condor_object.tar_memb["community"].unique()))+" Modularity: "+str(condor_object.modularity))
    return T,R



def run(filename,jaccard,resolution_graph,resolution_graphR,all_resolutions,processes=10):

    minres = all_resolutions.pop(0)
    # This is done only once.
    network = pd.read_csv(filename,header=None)
    condor_object = condor.condor_object(dataframe=network,silent=True)
    
    condor_object.initial_community(resolution=minres)
    condor_object.brim(resolution=minres)
    maxc = max(max(condor_object.tar_memb["community"]),max(condor_object.reg_memb["community"]))
    B,_,_,_,gn,rg = condor_object.matrices(maxc+6,1)
    A,_,_,_,_,_ = condor_object.matrices(maxc+6,0)

    clT = sorted(condor_object.tar_memb["community"].unique())
    clR = sorted(condor_object.reg_memb["community"].unique())
    T = sp.sparse.coo_matrix(np.matrix([np.array(condor_object.tar_memb["community"]==clT[i]).astype(int) for i in clT])).tocsr()
    R = sp.sparse.coo_matrix(np.matrix([np.array(condor_object.reg_memb["community"]==clR[i]).astype(int) for i in clR])).tocsr()

    cluT = ClusterGraph()
    cluR = ClusterGraph()
    cluT.graph['sim_threshold'] = jaccard
    cluT.graph['num_leaves'] = len(gn)

    cluR.graph['sim_threshold'] = jaccard
    cluR.graph['num_leaves'] = len(rg)

    _arg_tuples = [(condor_object, res) for res in all_resolutions]

    #_arg_tuples = [(condor_object, res, maxc+3) for res in all_resolutions]

    with mp.Pool(processes=processes) as pool:
            results = pool.starmap(run_alg, _arg_tuples)  # results contains "partition" class
    all_resolutions = [minres]+all_resolutions
    results = [[sp.sparse.coo_matrix(T).tocsr(),sp.sparse.coo_matrix(R).tocsr()]]+results

    for i in range(len(all_resolutions)):
            nodename = '{:.4f}'.format(all_resolutions[i])
            resolution_graph.nodes[nodename]['matrix'] = results[i][0]
            resolution_graphR.nodes[nodename]['matrix'] = results[i][1]
            cluT.add_clusters(resolution_graph, all_resolutions[i])
            cluR.add_clusters(resolution_graphR, all_resolutions[i])
    
    return cluT,cluR,gn,rg,A,B


def weave_and_out(T,R,gn,rg,oR,oT,A):
    
    cluR_collapsed = [x[0] for x in R]
    len_componentR = [x[1] for x in R]
    cluR_collapsed.insert(0, np.ones(len(cluR_collapsed[0]), ))
    len_componentR.insert(0, 0)

    cluT_collapsed = [x[0] for x in T]
    len_componentT = [x[1] for x in T]
    cluT_collapsed.insert(0, np.ones(len(cluT_collapsed[0]), ))
    len_componentT.insert(0, 0)
    
    wvR = weaver.Weaver()
    wvT = weaver.Weaver()
    
    T = wvT.weave(cluT_collapsed, terminals=list(gn.keys()), boolean=True, levels=False, merge=True, cutoff=0.75)
    output_all(wvT, list(gn.keys()), oT, persistence=len_componentT)
    
    R = wvR.weave(cluR_collapsed, terminals=list(rg.keys()), boolean=True, levels=False, merge=True, cutoff=0.75)
    output_all(wvR, list(rg.keys()), oR, persistence=len_componentR)
    
    
    fileR = open("cocluster_"+oR+"_Reg.txt","w")
    fileT = open("cocluster_"+oT+"_Tar.txt","w")

    fileR.write("ClusterR\tClusterT\tLevel\n")
    fileT.write("ClusterT\tClusterR\tLevel\n")

    for k in range(1,10):
        ccR,ccT = co_cluster_k(wvR,wvT,k,A)
        if (ccR,ccT) == (0,0): break

        fileR.writelines(["Cluster"+str(k)+"-"+str(i)+"\t"+"Cluster"+str(k)+"-"+str(ccR[i])+"\t"+str(k)+"\n" for i in range(0,len(ccR))])


        fileT.writelines(["Cluster"+str(k)+"-"+str(i)+"\t"+"Cluster"+str(k)+"-"+str(ccT[i])+"\t"+str(k)+"\n" for i in range(0,len(ccT))])

    fileR.close()
    fileT.close()



def wv_clust(wv): 
    wv_clusts = []
    for v, vdata in wv.hier.nodes(data=True):
        if not isinstance(v, tuple): #Skip the leaves.
            continue
        ind = vdata['index']
        name = 'Cluster{}-{}'.format(str(v[0]), str(v[1]))
        if isinstance(ind, int):
            wv_clusts.append([name, wv._assignment[ind]])
        else:
            wv_clusts.append([name, wv._assignment[ind[0]]])
    wv_clusts = sorted(wv_clusts, key=lambda x: np.sum(x[1]), reverse=True)
    return wv_clusts

def level_k_coms(wv_clusts,k):
    return [wv_clusts[i][1] for i in range(0,len(wv_clusts)) if int(wv_clusts[i][0][7])==k]




def co_cluster_k(wv1,wv2,k,matrix):
    cl1 = level_k_coms(wv_clust(wv1),k)
    cl2 = level_k_coms(wv_clust(wv2),k)
    
    if cl1 == [] or cl2 == []: return (0,0)
    
    cc1 = list()
    for i in range(0,len(cl1)):
        shw_i = [matrix[:,cl1[i]].transpose()[:,cl2[j]].sum() for j in range(0,len(cl2))]
        cc1.append(shw_i.index(max(shw_i)))
        
    cc2 = list()
    for i in range(0,len(cl2)):
        shw_i = [matrix[:,cl1[j]].transpose()[:,cl2[i]].sum() for j in range(0,len(cl1))]
        cc2.append(shw_i.index(max(shw_i)))
    return cc1,cc2




def bihidef(filename,
            jaccard=0.75,
            minres=0.001,
            maxres=10,
            density=0.1,
            processes=10,
            neighbors=10,
            min_diff_resolution=0.001,
            k=2,p=50,
            oR="pvr",oT="pvg"):
    
    resolution_graph,resolution_graphR,all_resolutions = create_resolution_graph(minres=minres,maxres=maxres,density=density,neighbors=neighbors,min_diff_resolution = min_diff_resolution)
    print("Computing community structure for "+str(len(all_resolutions))+" resolution points")
    
    cluT,cluR,gn,rg,A,B = run(filename=filename,jaccard=jaccard,resolution_graph=resolution_graph,resolution_graphR=resolution_graphR,all_resolutions=all_resolutions,processes=processes)
    
    consensusR = consensus(cluR,k=k,p=p)
    consensusT = consensus(cluT,k=k,p=p)
    
    gn = {k[4:]:gn[k] for k in gn.keys()}
    rg = {k[4:]:rg[k] for k in rg.keys()}
    
    weave_and_out(consensusT,consensusR,gn,rg,oR,oT,A)   




