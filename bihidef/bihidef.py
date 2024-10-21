import numpy as np
import networkx as nx
from netZooPy import condor
import pandas as pd
import copy
from scipy.sparse import *
import scipy as sp
import multiprocessing as mp


from hidef.hidef_finder import ClusterGraph,update_resolution_graph
from hidef.hidef_finder import collapse_cluster_graph,consensus
from hidef.hidef_finder import output_all
from hidef.utils import jaccard_matrix
from hidef import weaver


def create_resolution_graph(minres=0.001, maxres=10, density=0.1, neighbors=10, min_diff_resolution=0.001):
    """
    Creates two resolution graphs for a range of resolutions. 
    These graphs are used to identify community structures at multiple resolutions 
    and are attached to each node set in a bipartite graph.

    Parameters:
    - minres (float): The minimum resolution to initiate the process (larger communities).
    - maxres (float): The maximum resolution to stop the process (smaller communities).
    - density (float): The density threshold that defines connections between resolutions.
    - neighbors (int): The number of neighboring resolutions considered for each update.
    - min_diff_resolution (float): Minimum allowed difference between resolutions; stops splitting.

    Returns:
    - resolution_graph (networkx.Graph): The generated resolution graph for community detection.
    - resolution_graphR (networkx.Graph): A deep copy of the resolution graph for separate use.
    - all_resolutions (list): A list containing all resolutions processed in the graph.
    """

    # Initialize an empty graph to store all resolution nodes
    resolution_graph = nx.Graph()

    # Initialize the list with a single range of resolutions to explore (from minres to maxres)
    stack_res_range = []
    stack_res_range.append((minres, maxres))   

    # Add the initial boundary resolutions to the graph
    _ = update_resolution_graph(resolution_graph, minres, density, neighbors)
    _ = update_resolution_graph(resolution_graph, maxres, density, neighbors)

    # Track all resolutions that have been added to the graph
    all_resolutions = [minres, maxres]

    # Continue until all ranges in the stack have been processed
    while stack_res_range:
        # Retrieve and remove the first range of resolutions from the stack
        current_range = stack_res_range.pop(0)
        resname1, resname2 = '{:.4f}'.format(current_range[0]), '{:.4f}'.format(current_range[1])

        # If the difference between resolutions is too small, stop further splitting
        if round(current_range[1] - current_range[0], 4) <= min_diff_resolution:
            continue

        # Check if the nodes for these resolutions are already padded; if so, skip
        if resolution_graph.nodes[resname1].get('padded', False) and resolution_graph.nodes[resname2].get('padded', False):
            continue

        # Calculate a new intermediate resolution between the current range endpoints
        new_resolution = np.round(np.sqrt(current_range[1] * current_range[0]), 4)

        # Add the new resolution range back onto the stack for further exploration
        stack_res_range.append((current_range[0], new_resolution))
        stack_res_range.append((new_resolution, current_range[1]))

        # Append the new resolution to the list and update the graph with it
        all_resolutions.append(new_resolution)
        _ = update_resolution_graph(resolution_graph, new_resolution, density, neighbors)

    # Create a deep copy of the original resolution graph to maintain two independent versions
    resolution_graphR = copy.deepcopy(resolution_graph)
    
    # Return both graphs and the list of all resolutions
    return resolution_graph, resolution_graphR, all_resolutions



def run_alg(condor_object, resolution):
    """
    Executes the community detection algorithm on the provided condor_object using a specified resolution. 
    It initializes communities and calculates membership matrices for target and regulator members.

    Parameters:
    - condor_object: An object containing the data and methods required for community detection.
    - resolution (float): The resolution parameter that influences community structure detection.

    Returns:
    - T (scipy.sparse.csr_matrix): A sparse matrix indicating membership of target nodes in detected communities.
    - R (scipy.sparse.csr_matrix): A sparse matrix indicating membership of regulator nodes in detected communities.
    """

    # Initialize the community detection process for the specified resolution
    condor_object.initial_community(resolution=resolution)

    # Apply the BRIM algorithm on the condor_object for the specified resolution
    condor_object.brim(resolution=resolution)

    # Extract unique community identifiers for target and regulator members and sort them
    clT = sorted(condor_object.tar_memb["community"].unique())
    clR = sorted(condor_object.reg_memb["community"].unique())

    # Create a sparse matrix for target communities (T), where each column represents a community
    T = sp.sparse.coo_matrix(np.matrix([np.array(condor_object.tar_memb["community"] == i).astype(int) for i in clT])).tocsr()

    # Create a sparse matrix for regulator communities (R), where each column represents a community
    R = sp.sparse.coo_matrix(np.matrix([np.array(condor_object.reg_memb["community"] == i).astype(int) for i in clR])).tocsr()

    # Print the resolution used, number of unique communities found, and the modularity of the condor_object
    print("Resolution: " + str(resolution) + " NComs: " + str(len(condor_object.tar_memb["community"].unique())) + " Modularity: " + str(condor_object.modularity))

    # Return the sparse matrices for target and regulator communities
    return T, R




def run(filename, jaccard, resolution_graph, resolution_graphR, all_resolutions, processes=10):
    """
    Executes the community detection and similarity graph construction process for a network defined in the given CSV file.
    
    Parameters:
    - filename (str): The path to the CSV file containing network data.
    - jaccard (float): The Jaccard similarity threshold used for cluster graph creation.
    - resolution_graph (networkx.Graph): The graph to store resolution-based information for target communities.
    - resolution_graphR (networkx.Graph): The graph to store resolution-based information for regulator communities.
    - all_resolutions (list of float): A list of resolution values to be used in community detection.
    - processes (int): The number of parallel processes to use for running the algorithm (default is 10).
    
    Returns:
    - cluT (ClusterGraph): A graph representing clusters based on target communities.
    - cluR (ClusterGraph): A graph representing clusters based on regulator communities.
    - gn (array-like): Data related to target nodes.
    - rg (array-like): Data related to regulator nodes.
    - A (scipy.sparse.csr_matrix): The adjacency matrix for the network.
    - B (scipy.sparse.csr_matrix): The bipartite graph matrix for the network.
    """

    # Extract the minimum resolution from the list of all resolutions
    minres = all_resolutions.pop(0)

    # Load the network data from the specified CSV file
    network = pd.read_csv(filename, header=None)

    # Create a condor_object from the loaded network data
    condor_object = condor.condor_object(dataframe=network, silent=True)

    # Initialize the community detection process with the minimum resolution
    condor_object.initial_community(resolution=minres)

    # Apply the BRIM algorithm to detect communities using the minimum resolution
    condor_object.brim(resolution=minres)

    # Determine the maximum community identifier for target and regulator memberships
    maxc = max(max(condor_object.tar_memb["community"]), max(condor_object.reg_memb["community"]))

    # Generate matrices for community detection; B is for bipartite graph, A is for adjacency
    # i have no idea why + 6
    B, _, _, _, gn, rg = condor_object.matrices(maxc + 6, 1)
    A, _, _, _, _, _ = condor_object.matrices(maxc + 6, 0)

    # Extract unique community identifiers for target and regulator members and sort them
    clT = sorted(condor_object.tar_memb["community"].unique())
    clR = sorted(condor_object.reg_memb["community"].unique())

    # Create sparse matrices representing the membership of target and regulator communities
    T = sp.sparse.coo_matrix(np.matrix([np.array(condor_object.tar_memb["community"] == clT[i]).astype(int) for i in clT])).tocsr()
    R = sp.sparse.coo_matrix(np.matrix([np.array(condor_object.reg_memb["community"] == clR[i]).astype(int) for i in clR])).tocsr()

    # Initialize ClusterGraph objects for target and regulator communities
    cluT = ClusterGraph()
    cluR = ClusterGraph()

    # Set graph attributes for the clustering process
    cluT.graph['sim_threshold'] = jaccard
    cluT.graph['num_leaves'] = len(gn)

    cluR.graph['sim_threshold'] = jaccard
    cluR.graph['num_leaves'] = len(rg)

    # Prepare arguments for parallel processing with the community detection algorithm
    _arg_tuples = [(condor_object, res) for res in all_resolutions]

    # Use multiprocessing to execute the run_alg function in parallel
    with mp.Pool(processes=processes) as pool:
        results = pool.starmap(run_alg, _arg_tuples)  # results contains "partition" class

    # Reconstruct the resolutions list by adding the minimum resolution back
    all_resolutions = [minres] + all_resolutions
    
    # Include the initial matrices for the first resolution in results
    results = [[sp.sparse.coo_matrix(T).tocsr(), sp.sparse.coo_matrix(R).tocsr()]] + results

    # Update the resolution graphs with matrices and add clusters for each resolution
    for i in range(len(all_resolutions)):
        nodename = '{:.4f}'.format(all_resolutions[i])
        resolution_graph.nodes[nodename]['matrix'] = results[i][0]
        resolution_graphR.nodes[nodename]['matrix'] = results[i][1]
        cluT.add_clusters(resolution_graph, all_resolutions[i])
        cluR.add_clusters(resolution_graphR, all_resolutions[i])

    # Return the cluster graphs and matrices for further analysis
    return cluT, cluR, gn, rg, A, B



def weave_and_out(T, R, gn, rg, oR, oT, A):
    """
    Performs weaving operations on the given community matrices T and R, and outputs the co-clustering results 
    for both target and regulator communities.
    
    Parameters:
    - T (list of tuples): List of tuples representing the target community structure, 
                          where each tuple contains the cluster representation and its length.
    - R (list of tuples): List of tuples representing the regulator community structure, 
                          where each tuple contains the cluster representation and its length.
    - gn (dict): A dictionary mapping node identifiers to target node properties.
    - rg (dict): A dictionary mapping node identifiers to regulator node properties.
    - oR (str): Output file name suffix for regulator clusters.
    - oT (str): Output file name suffix for target clusters.
    - A (scipy.sparse.csr_matrix): The adjacency matrix of the network.
    
    This function writes the co-clustered information to text files for both target and regulator clusters.
    """

    # Extract collapsed cluster representations and lengths for regulator communities
    cluR_collapsed = [x[0] for x in R]  # Get the first element (clusters) from R
    len_componentR = [x[1] for x in R]  # Get the second element (lengths) from R
    cluR_collapsed.insert(0, np.ones(len(cluR_collapsed[0]), ))  # Insert a vector of ones at the start
    len_componentR.insert(0, 0)  # Insert a zero length at the start

    # Extract collapsed cluster representations and lengths for target communities
    cluT_collapsed = [x[0] for x in T]  # Get the first element (clusters) from T
    len_componentT = [x[1] for x in T]  # Get the second element (lengths) from T
    cluT_collapsed.insert(0, np.ones(len(cluT_collapsed[0]), ))  # Insert a vector of ones at the start
    len_componentT.insert(0, 0)  # Insert a zero length at the start
    
    # Initialize Weaver objects for weaving operations
    wvR = weaver.Weaver()  # Weaver for regulator communities
    wvT = weaver.Weaver()  # Weaver for target communities
    
    # Weave the target communities using the specified parameters
    T = wvT.weave(cluT_collapsed, terminals=list(gn.keys()), boolean=True, levels=False, merge=True, cutoff=0.75)
    # Output results for the target weaving operation
    output_all(wvT, list(gn.keys()), oT, persistence=len_componentT)

    # Weave the regulator communities using the specified parameters
    R = wvR.weave(cluR_collapsed, terminals=list(rg.keys()), boolean=True, levels=False, merge=True, cutoff=0.75)
    # Output results for the regulator weaving operation
    output_all(wvR, list(rg.keys()), oR, persistence=len_componentR)

    # Open files for writing the co-clustered results for regulator and target clusters
    fileR = open("cocluster_" + oR + "_Reg.txt", "w")
    fileT = open("cocluster_" + oT + "_Tar.txt", "w")

    # Write header lines for the output files
    fileR.write("ClusterR\tClusterT\tLevel\n")
    fileT.write("ClusterT\tClusterR\tLevel\n")

    # Iterate through cluster levels to extract co-clustered data
    for k in range(1, 10):
        ccR, ccT = co_cluster_k(wvR, wvT, k, A)  # Get co-clustering results for level k
        if (ccR, ccT) == (0, 0):  # Break the loop if no further clusters are found
            break

        # Write co-clustering results for regulator clusters to fileR
        fileR.writelines(["Cluster" + str(k) + "-" + str(i) + "\t" + "Cluster" + str(k) + "-" + str(ccR[i]) + "\t" + str(k) + "\n" 
                          for i in range(0, len(ccR))])

        # Write co-clustering results for target clusters to fileT
        fileT.writelines(["Cluster" + str(k) + "-" + str(i) + "\t" + "Cluster" + str(k) + "-" + str(ccT[i]) + "\t" + str(k) + "\n" 
                          for i in range(0, len(ccT))])

    # Close the output files after writing all data
    fileR.close()
    fileT.close()



def wv_clust(wv):
    """
    Extracts and organizes cluster assignments from the given Weaver object.

    Parameters:
    - wv (Weaver): An instance of the Weaver class that contains hierarchical cluster data.

    Returns:
    - List[List[str, ndarray]]: A sorted list of clusters, each represented by a name and their corresponding assignments.
    """
    
    wv_clusts = []  # Initialize an empty list to hold cluster names and their assignments
    
    # Iterate through the nodes in the Weaver object's hierarchical structure
    for v, vdata in wv.hier.nodes(data=True):
        if not isinstance(v, tuple):  # Skip nodes that are not tuples (i.e., leaves)
            continue
        
        ind = vdata['index']  # Retrieve the index of the current node
        name = 'Cluster{}-{}'.format(str(v[0]), str(v[1]))  # Create a cluster name based on the node tuple

        # Check if the index is a single integer or a tuple
        if isinstance(ind, int):
            # Append the cluster name and its assignment to the list
            wv_clusts.append([name, wv._assignment[ind]])
        else:
            # Append the cluster name and the assignment corresponding to the first index of the tuple
            wv_clusts.append([name, wv._assignment[ind[0]]])

    # Sort the clusters based on the sum of their assignments in descending order
    wv_clusts = sorted(wv_clusts, key=lambda x: np.sum(x[1]), reverse=True)
    
    return wv_clusts  # Return the sorted list of clusters


def level_k_coms(wv_clusts, k):
    """
    Retrieves community assignments corresponding to clusters at a specified level.

    Parameters:
    - wv_clusts (List[List[str, ndarray]]): A list of clusters, where each cluster is represented by a name and its assignments.
    - k (int): The level of clusters to filter by, where the level is indicated in the cluster name.

    Returns:
    - List[ndarray]: A list of community assignments for clusters at the specified level.
    """

    # Extract community assignments from clusters that match the specified level k
    return [
        wv_clusts[i][1]  # Get the community assignments
        for i in range(0, len(wv_clusts))  # Iterate through all clusters
        if int(wv_clusts[i][0][7]) == k  # Check if the level (extracted from the name) matches k
    ]





def co_cluster_k(wv1, wv2, k, matrix):
    """
    Identifies the co-clustering relationships between two hierarchical clusterings at a specified level.

    Parameters:
    - wv1: The first Weaver object containing the first set of clusters.
    - wv2: The second Weaver object containing the second set of clusters.
    - k (int): The level of clusters to consider for co-clustering.
    - matrix (ndarray): A matrix representing the pairwise relationships between nodes.

    Returns:
    - Tuple[List[int], List[int]]: Two lists indicating the indices of the most similar clusters in
      `wv1` and `wv2` for the specified level `k`. Returns (0, 0) if either cluster list is empty.
    """

    # Retrieve community assignments for clusters at level k from both Weaver objects
    cl1 = level_k_coms(wv_clust(wv1), k)
    cl2 = level_k_coms(wv_clust(wv2), k)

    # Return (0, 0) if either list of clusters is empty
    if cl1 == [] or cl2 == []:
        return (0, 0)

    # List to store the index of the most similar cluster in cl2 for each cluster in cl1
    cc1 = list()
    for i in range(len(cl1)):
        # Calculate the similarity scores for each cluster in cl2 relative to the current cluster in cl1
        shw_i = [matrix[:, cl1[i]].transpose()[:, cl2[j]].sum() for j in range(len(cl2))]
        # Append the index of the most similar cluster in cl2 to cc1
        cc1.append(shw_i.index(max(shw_i)))

    # List to store the index of the most similar cluster in cl1 for each cluster in cl2
    cc2 = list()
    for i in range(len(cl2)):
        # Calculate the similarity scores for each cluster in cl1 relative to the current cluster in cl2
        shw_i = [matrix[:, cl1[j]].transpose()[:, cl2[i]].sum() for j in range(len(cl1))]
        # Append the index of the most similar cluster in cl1 to cc2
        cc2.append(shw_i.index(max(shw_i)))

    return cc1, cc2





def bihidef(
    filename,
    jaccard=0.75,
    minres=0.001,
    maxres=10,
    density=0.1,
    processes=10,
    neighbors=10,
    min_diff_resolution=0.001,
    k=2,
    p=50,
    oR="pvr",
    oT="pvg"
):
    """
    Computes community structures from a graph defined in an edgelist format and analyzes communities at multiple resolutions.

    Parameters:
    -----------
    filename : str
        Path to the edgelist file containing node connections (comma-separated). The file should not have a header.
        An optional third column may contain weights for the edges.

    jaccard : float, optional
        The Jaccard score threshold for considering two nodes as similar. Default is 0.75.

    minres : float, optional
        The minimum resolution value, which typically results in one large community. Default is 0.001.

    maxres : float, optional
        The maximum resolution value, allowing for the detection of many smaller communities. Default is 10.

    density : float, optional
        The density threshold to determine if two resolutions are considered proximal. Default is 0.1.

    processes : int, optional
        The number of threads to use for multiprocessing during community detection. Default is 10.

    neighbors : int, optional
        The number of closely related resolutions to include in the resolution graph. Default is 10.

    min_diff_resolution : float, optional
        The minimum difference between resolutions to consider when generating the resolution graph. Default is 0.001.

    k : int, optional
        The number of clusters to consider in the consensus step. Default is 2.

    p : int, optional
        A parameter related to the condensation step (check the hidef package for details). Default is 50.

    oR : str, optional
        Output prefix for the regulator graph results. Default is "pvr".

    oT : str, optional
        Output prefix for the target graph results. Default is "pvg".
    """
    
    # Create resolution graphs and resolution range.
    resolution_graph, resolution_graphR, all_resolutions = create_resolution_graph(
        minres=minres,
        maxres=maxres,
        density=density,
        neighbors=neighbors,
        min_diff_resolution=min_diff_resolution
    )
    print("Computing community structure for " + str(len(all_resolutions)) + " resolution points")
    
    # Create the cluster objects. Run condor on each resolution.
    cluT, cluR, gn, rg, A, B = run(
        filename=filename,
        jaccard=jaccard,
        resolution_graph=resolution_graph,
        resolution_graphR=resolution_graphR,
        all_resolutions=all_resolutions,
        processes=processes
    )
    
    # Run the persistence of community across resolution step.
    consensusR = consensus(cluR, k=k, p=p)
    consensusT = consensus(cluT, k=k, p=p)
    
    # Recover the original names of the nodes (remove tar_, reg_ added by condor).
    gn = {k[4:]: gn[k] for k in gn.keys()}
    rg = {k[4:]: rg[k] for k in rg.keys()}
    
    # Weaver (build the hierarchy) and output the community structure.
    weave_and_out(consensusT, consensusR, gn, rg, oR, oT, A)





