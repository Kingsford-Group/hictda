import numpy as np
import re
import gudhi as gd
import networkx as nx
import csv
import time
import sys
import ast
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import argparse
import scipy.stats
import glob
#import seaborn as sns

# read in persistence file: return persistence diagram as just [[b,d],[b,d],...] for use in bottleneckdist
def readPersisFile(filename):
    persisdiag = []
    fullpersisinfo = []
    with open(filename,'rt') as f:
        freader = csv.reader(f,delimiter='\t')
        for line in freader:
            persisdiag.append((int(line[0]),(float(line[1]),float(line[2]))))
            pair = ast.literal_eval(line[3])
            #pair = [x.strip() for x in pair]
            fullpersisinfo.append([int(line[0]), float(line[1]), float(line[2]), pair])
    return persisdiag,fullpersisinfo

def readDistMat(filename):
    distmat = []
    with open(filename,'rt') as f:
        freader = csv.reader(f,delimiter='\t')
        for line in freader:
            row = [float(x) for x in line]
            distmat.append(row)
    return np.array(distmat)

def computeWassersteinDist(persis1,persis2,p=1):
    # compute the p-Wasserstein distance for 2 persistence diagrams by finding maximal matching, then matching all un-matched vertices to the closest point on the diagonal, then adding up all L-infinity distances between matched vertices
    #first turn persistence diagrams into a bipartite graph, nodes are (birth,death) tuples of only H1 (loop) structures
    bg = nx.Graph()
    persis1nodes = [(pers[1][0],pers[1][1]) for pers in persis1 if pers[0] == 1]
    persis2nodes = [(pers[1][0],pers[1][1]) for pers in persis2 if pers[0] == 1]
    bg.add_nodes_from(persis1nodes, bipartite=0)
    bg.add_nodes_from(persis2nodes, bipartite=1)
    # add all edges w/ edge weight (l inf norm)
    edgelist = []
    for node1 in persis1nodes:
        for node2 in persis2nodes:
            edgelist.append((node1,node2,np.linalg.norm(np.array(node1)-np.array(node2),ord=np.inf)))
    bg.add_weighted_edges_from(edgelist)
    if not nx.is_connected(bg):
        print('something went wrong, graph is not connected')
        sys.exit()
    maxmatch = nx.bipartite.maximum_matching(bg)
    #maxmatch = nx.hopcroft_karp_matching(bg)
    #print('done with max matching')
    #add in matching to diagonal if node wasn't matched
    for idx,node in enumerate(persis1nodes+persis2nodes):
        if node not in maxmatch:
            diagmatch = (0.5*(node[0]+node[1]),0.5*(node[0]+node[1]))
            maxmatch[node] = diagmatch
            maxmatch[diagmatch] = node
            if idx < len(persis1nodes):
                bipart = 1
            else:
                bipart = 0
            bg.add_node(diagmatch,bipartite=bipart)
            bg.add_edge(node,diagmatch,weight=np.linalg.norm(np.array(node)-np.array(diagmatch),ord=np.inf))
    #now add up all distances from maxmatch and divide by 2 (matching dictionary includes each edge twice)
    dist = 0
    for node1,node2 in maxmatch.iteritems():
        dist += (bg.edge[node1][node2]['weight'])**p
    dist = (dist/2)/(len(maxmatch)/2)
    return dist

# compute bottleneck distance between 2 persistence diagrams
def computeBottleneckDist(persis1, persis2, e=0):
    # bottleneck distance takes in list of [birth,death] lists, doesn't include dimension
    persis1list = []
    persis2list = []
    for pair in persis1:
        persis1list.append([pair[1][0],pair[1][1]])
    for pair in persis2:
        persis2list.append([pair[1][0],pair[1][1]])
    dist = gd.bottleneck_distance(persis1list,persis2list,e)
    return dist

def plotBarcode(persis,filename):
    colors=['#fc8d59', '#91bfdb']
    maxy = len(persis)
    # make sure persis is ordered first by all H0 (ordered by death time), then all H1
    persissorted = sorted(persis,key=lambda x: x[0])
    for idx,bar in enumerate(persissorted):
        if bar[1][1] == np.inf:
            plt.plot([0,1],[maxy-idx,maxy-idx],c=colors[bar[0]])
        else:
            plt.plot(bar[1],[maxy-idx,maxy-idx],c=colors[bar[0]])
    plt.axis([0,1,0,idx])
    plt.yticks(fontsize=16)
    plt.xticks(fontsize=12)
    plt.xlabel('Radius', fontsize=12)
    plt.legend(['H0','H1'],frameon=False)
    plt.savefig(filename,bbox_inches='tight')
    plt.close()
    print('saved barcode plot to',filename)

def plotPersisDiagram(persis,filename):
    colors=['#fc8d59', '#91bfdb']
    for pt in persis:
        if pt[1][1] == np.inf:
            plt.scatter(0,1,c=colors[pt[0]],s=5)
        else:
            plt.scatter(pt[1][0],pt[1][1],c=colors[pt[0]],s=5)
    plt.plot([0,1],[0,1],c='k')
    plt.axis([-0.05,1.05,-0.05,1.05])
    plt.legend(['H0','H1'],frameon=False)
    plt.xlabel('birth radius')
    plt.ylabel('death radius')
    plt.savefig(filename,bbox_inches='tight')
    plt.close()
    print('saved persistence diagram to',filename)

def plotHeatMap(avgdists,celltypelabels,filename):
    fig,ax1 = plt.subplots(figsize=(8,8))
    #sns.heatmap(avgbndists,cmap='YlGnBu',xticklabels=celltypelabels,yticklabels=celltypelabels,square=True)
    plt.imshow(avgdists,cmap='YlGnBu')
    plt.colorbar()
    ax1.grid(False)
    plt.xticks(np.arange(len(celltypelabels)), celltypelabels,rotation=90)
    plt.yticks(np.arange(len(celltypelabels)), celltypelabels)
    plt.tick_params(top=False,right=False,left=False,bottom=False)
    plt.tight_layout()
    plt.savefig(filename,bbox_inches='tight')
    print('heatmap of average distances saved as',filename)

def writeDistMattoFile(dists, labels, filename):
    with open(filename,'w') as f:
        fwriter = csv.writer(f,delimiter='\t')
        for chrnum in range(dists.shape[-1]):
            fwriter.writerow(['chr'+str(chrnum+1)])
            fwriter.writerow(['']+labels)
            for rownum in range(dists.shape[0]):
                fwriter.writerow([labels[rownum]] + list(dists[rownum,:,chrnum]))
    print('wrote distance matrix to',filename)

def readDistsFromFile(filename,celltypelabels,bndists=[]):
    if len(bndists) == 0:
        bndists = np.zeros((len(celltypelabels),len(celltypelabels),22))
    bndict = {}
    with open(filename,'rt') as f:
        freader = csv.reader(f,delimiter='\t')
        for row in freader:
            if len(row) == 1 or row[0][:3] == 'chr':
                chrnum = int(row[0][3:])
            elif row[0] == '':
                collabels = row[1:]
            else:
                rowlabel = row[0]
                for idx,distval in enumerate(row[1:]):
                    col = collabels[idx]
                    if (rowlabel,col,chrnum) not in bndict and (col,rowlabel,chrnum) not in bndict:
                        bndict[(rowlabel,col,chrnum)] = float(distval)
    for key,dist in bndict.items():
        if key[0] in celltypelabels and key[1] in celltypelabels:
            idx1 = celltypelabels.index(key[0])
            idx2 = celltypelabels.index(key[1])
            bndists[idx1,idx2,key[2]-1] = dist
            bndists[idx2,idx1,key[2]-1] = dist
    return bndists

def plotDistsByChr(dists, metric, filename):
    xvals = []
    yvals = []
    colorpattern = ['#bdc9e1','#74a9cf','#2b8cbe','#045a8d']
    for x in range(dists.shape[0]-1):
        xlist = []
        ylist = []
        for chrnum in range(1,23):
            xlist.extend([chrnum])
            if 'H1' not in filename:
                ylist.extend([dists[x,x+1,chrnum-1]])
            else:
                ylist.extend([dists[0,x+1,chrnum-1]])
        xvals.append(xlist)
        yvals.append(ylist)
    for idx,xlist in enumerate(xvals):
        plt.scatter(xlist,yvals[idx],c=colorpattern[idx])
    if 'H1' in filename:
        legendlist = ('ESC to ME', 'ESC to MS', 'ESC to NP', 'ESC to TB')
    elif 'WTC' in filename:
        legendlist = ('PSC to MES','MES to CP','CP to CM')
    else:
        legendlist = ('ESC to MES','MES to CP','CP to CM','CM to FH')
    plt.legend(legendlist[:idx+1],frameon=False)
    ax = plt.gca()
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    if metric == 'Bottleneck':
        ax.set_ylim([0,0.65])
    plt.xticks(np.arange(1,23))
    plt.xlabel('Chromosome')
    plt.ylabel(metric+' distance')
    plt.tight_layout()
    plt.savefig(filename,bbox_inches='tight')
    plt.close()
    print('saved distances by chromosome figure to',filename)
    #return [xvals,yvals]


def main(fileseeds,celltypelabels,bnfilename,wassfilename,outputloc):

    if len(celltypelabels) == 0:
        celltypelabels = [idx for idx,_ in enumerate(fileseeds)]
    elif len(celltypelabels) != len(fileseeds):
        print('please make sure there is exactly one label per filename')

    for chrnum in range(1,23):
        persis = {}
        for idx,fileseed in enumerate(fileseeds):
            persisdiag,_ = readPersisFile(fileseed+str(chrnum)+'_persisdiagram.txt')
            persis[celltypelabels[idx],idx] = persisdiag
            barcodefilename = outputloc+celltypelabels[idx]+'_chr'+str(chrnum)+'_barcodeplot.png'
            plotBarcode(persisdiag,barcodefilename)
            persisfilename = outputloc+celltypelabels[idx]+'_chr'+str(chrnum)+'_persistenceplot.png'
            plotPersisDiagram(persisdiag,persisfilename)

    if len(bnfilename) == 0:
        bndists = np.zeros((len(celltypelabels),len(celltypelabels),22))
        for chrnum in range(1,23):
            for key1,persis1 in persis.items():
                for key2,persis2 in persis.items():
                    idx1 = key1[1]
                    idx2 = key2[1]
                    if bndists[idx1,idx2,chrnum-1] != 0 or key1 == key2: continue
                    dist = computeBottleneckDist(persis1,persis2)
                    bndists[idx1,idx2,chrnum-1] = dist
                    bndists[idx2,idx1,chrnum-1] = dist
            writeDistMattoFile(bndists,celltypelabels,outputloc+'allsamples_bottleneckdists.txt')
    else:
        bndists = []
        for bnfile in bnfilename:
            bndists = readDistsFromFile(bnfile,celltypelabels,bndists)            
    if len(wassfilename) == 0:
        wassdists = np.zeros((len(celltypelabels),len(celltypelabels),22))
        for key1,persis1 in persis.items():
            for key2,persis2 in persis.items():
                idx1 = key1[1]
                idx2 = key2[1]
                if wassdists[idx1,idx2,chrnum-1] != 0 or key1 == key2: continue
                dist = computeWassersteinDist(persis1,persis2)
                wassdists[idx1,idx2,chrnum-1] = dist
                wassdists[idx2,idx1,chrnum-1] = dist
        writeDistMattoFile(wassdists,celltypelabels,outputloc+'allsamples_wassersteindists.txt')
    else:
        wassdists = []
        for wassfile in wassfilename:
            wassdists = readDistsFromFile(wassfile,celltypelabels,wassdists)
    
    # need to break up by cell type for these
    celltypes = [ct.split('_')[0] for ct in celltypelabels]
    startidx = 0
    for idx,celltype in enumerate(celltypes[1:]):
        if celltype != celltypes[idx]:
            plotDistsByChr(bndists[startidx:idx+1,startidx:idx+1,:], 'Bottleneck', outputloc+celltypes[idx]+'_bottleneck_bychr.png')
            plotDistsByChr(wassdists[startidx:idx+1,startidx:idx+1,:], 'Wasserstein', outputloc+celltypes[idx]+'_wasserstein_bychr.png')
            startidx = idx+1
    plotDistsByChr(bndists[startidx:,startidx:,:], 'Bottleneck', outputloc+celltype+'_bottleneck_bychr.png')
    plotDistsByChr(wassdists[startidx:,startidx:,:], 'Wasserstein', outputloc+celltype+'_wasserstein_bychr.png')

    avgbndists = bndists.mean(2)
    plotHeatMap(avgbndists,celltypelabels,outputloc+'allsamples_avgbottleneck_heatmap.pdf')
    avgwassdists = wassdists.mean(2)
    plotHeatMap(avgwassdists,celltypelabels,outputloc+'allsamples_avgwasserstein_heatmap.pdf')


if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument('-i', nargs='+',help='Filenames to analyze. Will plot barcodes/persistence diagrams for all, and compute all pairwise bottleneck/wasserstein distances.')
    parser.add_argument('-l', nargs='+', default=[], help='Names of cell types, in same order as filenames')
    parser.add_argument('-b', default='', nargs='+', help='If bottleneck distances have been precomputed, input filename')
    parser.add_argument('-w', default='', nargs='+', help='If wasserstein distances have been precomputed, input filename')
    parser.add_argument('-o',type=str, help='Location of output files')

    args=parser.parse_args()
    main(args.i,args.l,args.b,args.w,args.o)
