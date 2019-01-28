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
    return (resolution,data_frequency_matrix)

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
    return (skeleton,diag)

def writeSkeletonToFile(skeleton, filename):
    with open(filename,'w') as f:
        fwriter = csv.writer(f, delimiter='\t')
        for simplex in skeleton:
            if len(simplex) > 1:
                fwriter.writerow(simplex)
    print('wrote simplex skeleton to',filename)

def writePersistencePairsToFile(perspairs, filename):
    with open(filename,'w') as f:
        fwriter = csv.writer(f, delimiter='\t')
        for pers in perspairs:
            fwriter.writerow(pers)
    print('wrote persistence pairs to',filename)

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

def main(input_file,output_path,output_name,resol):
    resolut,freq_mat=read_raw_HiC_data_no_split(input_file,resol)
    print("generate distance matrix...")
    distance_matrix=matrix_normalize([freq_mat])[0]
    writeDistMatToFile(distance_matrix,output_path+output_name+'_distmat.txt')
    t0=time.time()
    simplex_skel,persisdiag=TDA_func(distance_matrix,output_path+output_name+'_persisdiagram.txt',output_path+output_name+'_skeleton.txt')
    print('total time (s):',time.time()-t0)

if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument('-i', type=str,help= "input file of HiC contact matrix")
    parser.add_argument('-o',type=str, help='the name of output files')
    parser.add_argument('-p',type=str,help="the path of output files")
    parser.add_argument('-r',type=str,help="resolution of HiC input file")

    args=parser.parse_args()
    main(args.i,args.p,args.o,args.r)

