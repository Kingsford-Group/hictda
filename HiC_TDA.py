import numpy as np
import re
import os
import gudhi as gd
import time
#import dionysus as d
#from scipy.spatial.distance import squareform
import sys
import csv
import argparse
import networkx as nx

def read_raw_HiC_data(file):
    resolution=re.split('[_.]',os.path.basename(file).strip())[1]
    if(resolution[-2:]=='kb'):
        resolution=int(resolution[:-2])*1000
    elif(resolution[-2:]=='mb'):
        resolution=int(resolution[:-2])*1000000
    print(resolution)
    all_data_list=[]
    with open(file,"r")as f:
        for line in f:
            line = line.strip()
            line = re.sub(r'\s+',' ', line)
            line = line.split(' ')
            #print(line)
            all_data_list.append([float(line[i]) for i in range(len(line))])
    raw_data_matrix=np.array(all_data_list)
    dim=int(max(raw_data_matrix[:,0].max(),raw_data_matrix[:,1].max())/resolution+1)
    #print(dim)
    data_frequency_matrix=np.zeros((dim,dim))
    for x,y,freq in all_data_list:
        data_frequency_matrix[int(x/resolution),int(y/resolution)]=freq
        data_frequency_matrix[int(y/resolution),int(x/resolution)]=freq
    #print(data_frequency_matrix)
    return (resolution,data_frequency_matrix)

def read_raw_HiC_data_no_split(file,reso):
    resolution=int(reso)
    print('resolution: ',resolution)
    all_data_list=[]
    with open(file,"r")as f:
        for line in f:
            line = line.strip()
            line = re.sub(r'\s+',' ', line)
            line = line.split(' ')
            #print(line)
            all_data_list.append([float(line[i]) for i in range(len(line))])
    raw_data_matrix=np.array(all_data_list)
    temp_min=min(raw_data_matrix[:,0].min(),raw_data_matrix[:,1].min())
    temp_max=max(raw_data_matrix[:,0].max(),raw_data_matrix[:,1].max())
    dim=int((temp_max-temp_min)/resolution+1)
    #print('Hi-C matrix size =',dim)
    data_frequency_matrix=np.zeros((dim,dim))
    for x,y,freq in all_data_list:
        data_frequency_matrix[int((x-temp_min)/resolution),int((y-temp_min)/resolution)]=freq
        data_frequency_matrix[int((y-temp_min)/resolution),int((x-temp_min)/resolution)]=freq
    #print(data_frequency_matrix)
    return (resolution,data_frequency_matrix,temp_min,temp_max)

def split_TAD(freq_matrix,TAD_result_file,resol):
    TAD_matrix_list=[]
    with open(TAD_result_file,"r")as f:
        for line in f:
            line = line.strip()
            line = re.sub(r'\s+',' ', line)
            line = line.split(' ')
            #print([line[1],line[2]])
            index_x=int(int(line[1])/resol)
            index_y=int(int(line[2])/resol)+1
            TAD_matrix_list.append(freq_matrix[index_x:index_y,index_x:index_y])
            #print(freq_matrix[index_x:index_y,index_x:index_y])
            print(freq_matrix[index_x:index_y,index_x:index_y].shape)
    return TAD_matrix_list

def matrix_normalize(TAD_matrix_all):
    distance_matrix_all=[]
    for TAD_matrix in TAD_matrix_all:
        TAD_matrix=np.log(TAD_matrix+1)
        max_num = TAD_matrix.max()
        #print(max_num)
        TAD_matrix = TAD_matrix/(1.01*max_num)
        #print(TAD_matrix.shape[0])
        for i in range(TAD_matrix.shape[0]):
            TAD_matrix[i,i]=1.0
        #print(TAD_matrix)
        distance_matrix=1-TAD_matrix
        distance_matrix_all.append(distance_matrix)
        #print(distance_matrix)
    return distance_matrix_all

def TDA_func(distance_matrix, persfilename,skelfilename):
    print('distance matrix size =',distance_matrix.shape)
    rips_complex = gd.RipsComplex(distance_matrix=distance_matrix,max_edge_length=1.1)
    simplex_tree = rips_complex.create_simplex_tree(max_dimension=2)
    print('done creating simplex tree')
    #filtration = simplex_tree.get_filtration()
    skeleton = simplex_tree.get_skeleton(2)
    writeSkeletonToFile(skeleton,skelfilename)
    diag = simplex_tree.persistence()
    pairs = simplex_tree.persistence_pairs()
    fullpersinfo = []
    for pair in pairs:
        #print(pair)
        btime = simplex_tree.filtration(pair[0])
        dtime = simplex_tree.filtration(pair[1])
        #print(btime,dtime)
        try:
            diag.index((0,(btime,dtime)))
            htype = 0
            fullpersinfo.append([htype, btime, dtime, pair])
        except:
            diag.index((1,(btime,dtime)))
            htype = 1
            fullpersinfo.append([htype, btime, dtime, pair])
            #print('couldnt find persistence pair matching birth/death times', pair, (btime,dtime))
    writePersistencePairsToFile(fullpersinfo,persfilename)
    #simplex_tree.write_persistence_diagram(persfilename)
    #print('wrote persistence diagram to',persfilename)
    #print(diag)
    #data_0_dim=np.array([list(diag[i][1]) for i in range(len(diag)) if diag[i][0]==0])
    #data_1_dim=np.array([list(diag[i][1]) for i in range(len(diag)) if diag[i][0]==1])
    return (skeleton,diag,fullpersinfo)

def writeSkeletonToFile(skeleton, filename):
    with open(filename,'w') as f:
        fwriter = csv.writer(f, delimiter='\t')
        for simplex in skeleton:
            if len(simplex) > 1:
                fwriter.writerow(simplex)
    print('wrote simplex skeleton to',filename)
    return

def writePersistencePairsToFile(perspairs, filename):
    with open(filename,'w') as f:
        fwriter = csv.writer(f, delimiter='\t')
        for pers in perspairs:
            fwriter.writerow(pers)
    print('wrote persistence pairs to',filename)
    return

def randomlyPermuteDistMat(distance_matrix,flag=''):
    distmatsize = distance_matrix.shape
    n = distmatsize[0]
    permmat = np.zeros(distmatsize)
    if flag == 'edge':
        # randomly permute by row
        randperm = np.random.permutation(n)
        for rownum in range(n):
            permmat[randperm[rownum],:] = distance_matrix[rownum,randperm]
    elif flag == 'rand':
        # randomly permute all dist values (in upper section, to preserve symmetry)
        permidx = np.triu_indices(n,1)
        alldistvals = distance_matrix[permidx]
        permvals = np.random.permutation(alldistvals)
        permmat[permidx] = permvals
        permidx_lower = np.tril_indices(n,-1)
        permmat[permidx_lower] = permmat.T[permidx_lower]
    elif flag == 'dist':
        # matrix is purely distance dependent (same averages along non-main diagonals as original) + noise
        for diagnum in range(n-1):
            diagvals = np.diag(distance_matrix,diagnum+1)
            avgval = np.mean(diagvals)
            std = np.std(diagvals)
            newdiag = np.random.normal(avgval,std,len(diagvals))
            diagmat = np.diag(newdiag,diagnum+1)
            permmat += diagmat
            permmat += diagmat.T
    return permmat

def writeDistMatToFile(distmat,filename):
    with open(filename,'w') as f:
        fwriter = csv.writer(f,delimiter='\t')
        for row in distmat:
            fwriter.writerow(row)
    print('wrote distance matrix to',filename)
    return

def writeLoopToFile(path_data,filename):
    with open(filename,'w')as f:
        fwriter = csv.writer(f,delimiter='\t')
        for path in path_data:
            fwriter.writerow([path[0],"other edges with the same length:"+str(path[1])])
            fwriter.writerow(path[2])
    print('wrote the information of loops to ', filename)
    return


def generate_1_dim_simp_list_from_dist_mat(dist_matrix):
    mat_dim = dist_matrix.shape[0]
    dim_1_simp_list = [[i,j,dist_matrix[i,j]] for i in range(0,mat_dim-1) for j in range(i+1,mat_dim)]
    dim_1_simp_list.sort(key=lambda x: [x[2],max(x[0],x[1]),min(x[0],x[1])])
    return dim_1_simp_list

def get_one_dim_persis(fullpersis):
    persis_1_dim = []
    persis_1_dim_dict = {}
    for dim,birth,death,simp_pair in fullpersis:
        if(dim==1):
            persis_1_dim.append([birth,death,simp_pair[0]])
            if(tuple(simp_pair[0]) not in persis_1_dim_dict.keys()):
                persis_1_dim_dict[tuple(simp_pair[0])] = [[birth,death,simp_pair[1]]]
            else:
                persis_1_dim_dict[tuple(simp_pair[0])].append([birth,death,simp_pair[1]])
    persis_1_dim.sort(key=lambda x: [x[0]-x[1]])
    return persis_1_dim,persis_1_dim_dict

def check_same_length_edge(index,length,dim_1_simp_list):
    edge_num = 1
    temp_index = index-1
    all_simp = len(dim_1_simp_list)
    while True:
        if(dim_1_simp_list[temp_index][2]-length==0):
            edge_num+=1
            temp_index-=1
            if(temp_index<0):
                break
        else:
            break
    temp_index = index+1
    while True:
        if(dim_1_simp_list[temp_index][2]-length==0):
            edge_num+=1
            temp_index+=1
            if(temp_index>=all_simp):
                break
        else:
            break
    return edge_num

def get_loop(reso,minimal_bin,dist_mat,dim_1_simp_list,persis_1_dim_list,output_filename):
    loop_num = len(persis_1_dim_list)
    path_info = []
    for i in range(loop_num):
        essen_edge = persis_1_dim_list[i][2]
        essen_index = dim_1_simp_list.index([min(essen_edge[0],essen_edge[1]),max(essen_edge[0],essen_edge[1]),dist_mat[essen_edge[0],essen_edge[1]]])
        #check if there are other edges having the same length as the essential edge
        multi_edge = check_same_length_edge(essen_index,dist_mat[essen_edge[0],essen_edge[1]],dim_1_simp_list)
        if(multi_edge>1):
            print("essential edge: bin: %d bin: %d length: %f"%(int(minimal_bin+reso*essen_edge[0]),int(minimal_bin+reso*essen_edge[1]),dist_mat[essen_edge[0],essen_edge[1]]))
            print("another %d edges have the same length as the essential edge..."%(multi_edge-1))
        #get shortest path
        G=nx.Graph()
        G.add_weighted_edges_from(dim_1_simp_list[:essen_index])
        path_info.append([[int(minimal_bin+reso*i) for i in essen_edge],multi_edge-1,[int(minimal_bin+reso*i) for i in nx.shortest_path(G, source=essen_edge[0], target=essen_edge[1],weight="weight")]])
    writeLoopToFile(path_info,output_filename)
    return



def main(input_file,output_path,output_name,resol):
    t0 = time.time()
    resolut,freq_mat,min_bin,_=read_raw_HiC_data_no_split(input_file,resol)
    print("Generate distance matrix...")
    distance_matrix=matrix_normalize([freq_mat])[0]
    writeDistMatToFile(distance_matrix,output_path+output_name+'_distmat.txt')
    print("Calculate persistent homology...")
    simplex_skel,persisdiag,full_persis_info=TDA_func(distance_matrix,output_path+output_name+'_persisdiagram.txt',output_path+output_name+'_skeleton.txt')
    print("Loop trace back...")
    persis_dim_1_list,_ = get_one_dim_persis(full_persis_info)
    dim_1_simp = generate_1_dim_simp_list_from_dist_mat(distance_matrix)
    get_loop(resolut,min_bin,distance_matrix,dim_1_simp,persis_dim_1_list,output_path+output_name+'_loop_information.txt')
    print('total time (s):',time.time()-t0)
    return

if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument('-i', type=str,help= "input file of HiC contact matrix")
    parser.add_argument('-o',type=str, help='the name of output files')
    parser.add_argument('-p',type=str,help="the path of output files")
    parser.add_argument('-r',type=str,help="resolution of HiC input file")

    args=parser.parse_args()
    main(args.i,args.p,args.o,args.r)

