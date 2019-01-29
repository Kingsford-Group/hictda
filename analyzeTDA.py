import numpy as np
import re
#import gudhi as gd
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

# read in simplex file
def readSimplexFile(filename):
    simplexskel = []
    with open(filename,'rt') as f:
        freader = csv.reader(f,delimiter='\t')
        for line in freader:
            nodes = ast.literal_eval(line[0])
            #nodes = [x.strip() for x in nodes]
            simplexskel.append((nodes,float(line[1])))
    return simplexskel

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
    #print('done adding diagonal points')
    dist = 0
    for node1,node2 in maxmatch.iteritems():
        dist += (bg.edge[node1][node2]['weight'])**p
    dist = (dist/2)/(len(maxmatch)/2)
    return dist

# compute bottleneck distance between 2 persistence diagrams
#def computeBottleneckDist(persis1, persis2, e=0):
#    # bottleneck distance takes in list of [birth,death] lists, doesn't include dimension
#    persis1list = []
#    persis2list = []
#    for pair in persis1:
#        persis1list.append([pair[1][0],pair[1][1]])
#    for pair in persis2:
#        persis2list.append([pair[1][0],pair[1][1]])
#    dist = gd.bottleneck_distance(persis1list,persis2list,e)
#    return dist

# plot barcode/persistence diagrams
def plotBarcode(persis,filename):
    barcodefig = gd.plot_persistence_barcode(persistence=persis)
    plt.savefig(filename,bbox_inches='tight')
    plt.close()
    print('saved barcode plot to',filename)

def plotBarcodeCustom(persis,filename,persisdict={}):
    colors=['#7a5195']
    nullcolors=['#7a5195','#ffa600','#ef5675','#003f5c'] # edge, random, dist, original
    #persis = [x for x in persis if x[0] == 0]
    maxy = len(persis)
    if len(persisdict) == 0:
        # make sure persis is ordered first by all H0 (ordered by death time), then all H1
        persissorted = sorted(persis,key=lambda x: x[0])
        for idx,bar in enumerate(persissorted):
            if bar[1][1] == np.inf:
                plt.plot([0,1],[maxy-idx,maxy-idx],c=nullcolors[2])
            else:
                plt.plot(bar[1],[maxy-idx,maxy-idx],c=nullcolors[2])
    else:
        maxsize = len(persisdict[(0,np.inf)])
        for idx,bar in enumerate(persis):
            if bar[0] == 1: continue
            structsize = len(persisdict[bar[1]])
            #if structsize < 2: continue
            if bar[1][1] == np.inf:
                plt.plot([0,1],[idx,idx],c=colors[0],alpha=(structsize)*10/maxsize)
                #plt.plot([0,1],[maxy-idx,maxy-idx],c=colors[0],alpha=structsize/maxsize)
            else:
                plt.plot(bar[1],[idx,idx],c=colors[0],alpha=(structsize)*10/maxsize)
                #plt.plot(bar[1],[maxy-idx,maxy-idx],c=colors[0],alpha=structsize/maxsize) 
    plt.axis([0,1,0,idx])
    plt.yticks(fontsize=16)
    plt.xticks(fontsize=12)
    #plt.xlabel('Radius', fontsize=12)
    #plt.legend(['H0','H1'],frameon=False)
    plt.savefig(filename,bbox_inches='tight')
    plt.close()
    print('saved barcode plot to',filename)

def plotPersisDiagram(persis,filename):
    persisdiag = gd.plot_persistence_diagram(persistence=persis,legend=True)
    plt.savefig(filename,bbox_inches='tight')
    plt.close()
    print('saved persistence diagram to',filename)

def plotPersisDiagramCustom(persis,filename):
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
    pass

def plotHeatMap(avgbndists,celltypelabels,filename):
    fig,ax1 = plt.subplots(figsize=(8,8))
    #sns.heatmap(avgbndists,cmap='YlGnBu',xticklabels=celltypelabels,yticklabels=celltypelabels,square=True)
    plt.imshow(avgbndists,cmap='YlGnBu')
    plt.colorbar()
    ax1.grid(False)
    plt.xticks(np.arange(len(celltypelabels)), celltypelabels,rotation=90)
    plt.yticks(np.arange(len(celltypelabels)), celltypelabels)
    #ax1.set_xticklabels(celltypelabels)
    #ax1.set_yticklabels(celltypelabels)
    plt.tick_params(top=False,right=False,left=False,bottom=False)
    plt.tight_layout()
    plt.savefig(filename,bbox_inches='tight')
    print('heatmap of average bottleneck distances saved as',filename)

def writeBottleneckDisttoFile(bndists, labels, filename):
    with open(filename,'w') as f:
        fwriter = csv.writer(f,delimiter='\t')
        for chrnum in range(bndists.shape[-1]):
            fwriter.writerow(['chr'+str(chrnum+1)])
            fwriter.writerow(['']+labels)
            for rownum in range(bndists.shape[0]):
                fwriter.writerow([labels[rownum]] + list(bndists[rownum,:,chrnum]))
    print('wrote bottleneck distances to',filename)

def readBottleneckDistsFromFile(filename,celltypelabels,bndists=[]):
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

def plotBottleneckByChr(bndists, filename):
    xvals = []
    yvals = []
    colorpattern = ['#bdc9e1','#74a9cf','#2b8cbe','#045a8d']
    for x in range(bndists.shape[0]-1):
        xlist = []
        ylist = []
        for chrnum in range(1,23):
            xlist.extend([chrnum])
            if 'H1' not in filename:
                ylist.extend([bndists[x,x+1,chrnum-1]])
            else:
                ylist.extend([bndists[0,x+1,chrnum-1]])
        xvals.append(xlist)
        yvals.append(ylist)
    #yvalsbychr = list(map(list,zip(*yvals)))
    #for chrnum,dists in enumerate(yvalsbychr):
    #    if np.max(dists) < 0.2:
    #        print('all dist values low on chr',chrnum+1)
    #    elif np.min(dists) > 0.2:
    #        print('all dist values high on chr',chrnum+1)
    #sys.exit()
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
    ax.set_ylim([0,0.65])
    plt.xticks(np.arange(1,23))
    plt.xlabel('Chromosome')
    plt.ylabel('Bottleneck distance')
    #plt.ylabel('Wasserstein distance')
    plt.tight_layout()
    plt.savefig(filename,bbox_inches='tight')
    plt.close()
    print('saved bottleneck distance by chromosome figure to',filename)
    return [xvals,yvals]

def readGeneFile(filename):
    genelocs = []
    with open(filename,'rt') as f:
        freader = csv.reader(f,delimiter='\t')
        next(freader)
        for line in freader:
            try:
                genelocs.append([line[0], int(line[1])])
            except:
                continue
    return genelocs

def plotGenesVBottleneck(bnbychr,genelocs,filename):
    chrs = [x[1] for x in genelocs]
    genedistdata = {}
    for chrnum in range(1,23):
        bndists = [y[chrnum-1] for y in bnbychr[1]]
        if np.max(bndists) < 0.3:
            print(chrnum)
        genedistdata[chrnum] = [chrs.count(chrnum),np.max(bndists)]
    xvals, yvals = [],[]
    for _,data in genedistdata.items():
        xvals.extend([data[0]])
        yvals.extend([data[1]])
    plt.scatter(xvals,yvals)
    plt.xlabel('Number of GO genes on a chromosome')
    plt.ylabel('Max bottleneck distance of chromosome')
    plt.tight_layout()
    plt.savefig(filename,bbox_inches='tight')
    print('saved gene vs bottleneck figure to',filename)

def findH0structures(fullpersisinfo, distmat):
    h0dict = {}
    for persis in fullpersisinfo:
        if persis[0] != 0: continue
        nodeschecked = np.zeros(distmat.shape[0])
        deathtime = persis[2]
        birthstruct = persis[3][0]
        h0nodes = np.array(birthstruct)
        while any(nodeschecked[h0nodes] == 0):
            nodestocheck = np.where(nodeschecked[h0nodes] == 0)[0]
            nodestocheck = [h0nodes[x] for x in nodestocheck]
            for node in nodestocheck:
                connectednodes = np.where(distmat[node,:] < deathtime)[0]
                connectednodes = [int(x) for x in connectednodes if x not in h0nodes]
                if len(connectednodes) > 0:
                    h0nodes = np.append(h0nodes,connectednodes)
                nodeschecked[node] = 1
        # need to check each connected node for its connections too (check to make sure not visiting same nodes more than once)
        h0dict[(persis[1],deathtime)] = h0nodes
    h0list = []
    for key,h0struct in h0dict.items():
        h0list.append([key,h0struct])
    sortedh0 = sorted(h0list,key=lambda x: x[0][1])
    return h0dict,sortedh0

def plotH0(sortedh0,chrnum,filename):
    centromerelocs = [[121535434, 124535434],[92326171, 95326171],[90504854,93504854],[49660117, 52660117],[46405641, 49405641],[58830166, 61830166],[58054331, 61054331],[43838887, 46838887],[47367679, 50367679],[39254935,42254935],[51644205,54644205],[34856694,37856694],[16000000,19000000],[16000000,19000000],[17000000,20000000],[35335801,38335801],[22263006,25263006],[15460898,18460898],[24681782,27681782],[26369569,29369569],[11288129,14288129],[13000000,16000000]] #hg19 centromere locations from http://genome.ucsc.edu/cgi-bin/hgTables
    centromerelocs = [[np.floor(x[0]/100000.0), np.ceil(x[1]/100000.0)] for x in centromerelocs]
    #persis = [x for x in persis if x[0] == 0]
    for h0 in sortedh0:
        if len(h0[1]) < 5: continue
        plt.plot([min(h0[1]),max(h0[1])],[h0[0][1],h0[0][1]],c='#7a5195')
        plt.tight_layout()
    plt.plot([centromerelocs[chrnum-1][0], centromerelocs[chrnum-1][0]], [0, 1],color='red')
    plt.plot([centromerelocs[chrnum-1][1], centromerelocs[chrnum-1][1]], [0, 1],color='red')
    plt.axis([0,h0[1][-1],0,1])
    plt.xlabel('Chromosome bin (100kb)')
    plt.ylabel('Death radius')
    plt.savefig(filename,bbox_inches='tight')
    plt.close()
    print('saved h0 figure to',filename)

def countH1Structures(persis):
    numh1 = 0
    for p in persis:
        if p[0] == 1:
            numh1 += 1
    return numh1

def plotH1Counts(h1counts, filename):
    colorpattern = ['lightcoral','red','firebrick','darkred','black']
    for idx,ylist in enumerate(list(h1counts)):
        plt.scatter(range(1,23),ylist,c=colorpattern[idx])
    if 'RUES' in filename:
        legendlist = ('ESC', 'MES','CP','CM','FH')
    elif 'WTC' in filename:
        legendlist = ('PSC', 'MES','CP','CM')
    else:
        legendlist = ('ESC', 'ME', 'MS', 'NP', 'TB')
    plt.legend(legendlist[:idx+1])
    plt.xticks(np.arange(1,23))
    plt.xlabel('Chromosome')
    plt.ylabel('Number of loops')
    plt.tight_layout()
    plt.savefig(filename,bbox_inches='tight')
    plt.close()
    print('saved H1 count by chromosome figure to',filename)

def main(fileseeds,celltypelabels,bnfilename,genefilename,outputloc):
    if len(fileseeds) == 1:
        for chrnum in range(22,23):
            if len(glob.glob(fileseeds[0]+str(chrnum)+'_persisdiagram.txt')) == 1:
                persis1,fullpersisinfo = readPersisFile(fileseeds[0]+str(chrnum)+'_persisdiagram.txt')
                #distmat = readDistMat(fileseeds[0]+str(chrnum)+'_distmat.txt')
                persis1 = [x[1] for x in persis1 if x[0] == 1]
            else:
                persis1,fullpersisinfo = readPersisFile(fileseeds[0])
                #totlen = len(persis1)
                persis1 = [x for x in persis1 if x[0] == 1]
                #print('number of h0 =',len(persis1))
                #print('number of h1 =',totlen - len(persis1))
            #h0dict,sortedh0 = findH0structures(fullpersisinfo, distmat)
            #plotH0(sortedh0,chrnum, outputloc+str(chrnum)+'_H0vis_bydeathradius.png')
            #continue
            barcodefilename = outputloc+'_H1only_barcodeplot.png'
            plotBarcodeCustom(persis1,barcodefilename)#,h0dict)
            #persisfilename = outputloc+str(chrnum)+'_persistenceplot.png'
            #plotPersisDiagramCustom(persis1,persisfilename)
    elif len(fileseeds) > 1:
        if len(celltypelabels) == 0:
            celltypelabels = [idx for idx,_ in enumerate(fileseeds)]
        elif len(celltypelabels) != len(fileseeds):
            print('please make sure there is exactly one label per filename')
            sys.exit()
        if len(bnfilename) == 0:
            bndists = np.zeros((len(celltypelabels),len(celltypelabels),22))
            h1counts = np.zeros((len(fileseeds),22))
            for chrnum in range(1,23):
                persis = {}
                for idx,fileseed in enumerate(fileseeds):
                    persisdiag,_ = readPersisFile(fileseed+str(chrnum)+'_persisdiagram.txt')
                    persis[celltypelabels[idx],idx] = persisdiag
                    #h1counts[idx,chrnum-1] = countH1Structures(persisdiag)
                if True:
                    #compute pairwise bottleneck distances for each chromosome
                    for key1,persis1 in persis.items():
                        for key2,persis2 in persis.items():
                            idx1 = key1[1]
                            idx2 = key2[1]
                            #if abs(idx1-idx2) > 1: continue
                            if bndists[idx1,idx2,chrnum-1] != 0 or key1 == key2: continue
                            dist = computeWassersteinDist(persis1,persis2)
                            #dist = computeBottleneckDist(persis1,persis2)
                            bndists[idx1,idx2,chrnum-1] = dist
                            bndists[idx2,idx1,chrnum-1] = dist
                            #print(key1[0],key2[0],dist)
                print('done with chr'+str(chrnum))
            #plotH1Counts(h1counts,outputloc+'_H1counts.png')
            writeBottleneckDisttoFile(bndists,celltypelabels,outputloc+'_wassersteindists.txt')
        else:
            bndists = []
            for bnfile in bnfilename:
                bndists = readBottleneckDistsFromFile(bnfile,celltypelabels,bndists)
            #if any(np.sum(np.sum(bndists,axis=2),axis=1) == 0):
            # still gotta compute some distances
            #idxtocompute = [i for i,x in enumerate(np.sum(np.sum(bndists,axis=2),axis=1)) if x == 0]
            for chrnum in range(22,0,-1):
                print('chrnum',chrnum)
                missingvals = np.where(bndists[:,:,chrnum-1] == 0)
                if len(missingvals[0]) > len(celltypelabels):            
                    idxtocompute = [(missingvals[0][idx], missingvals[1][idx]) for idx in range(len(missingvals[0])) if missingvals[0][idx] != missingvals[1][idx]]
                    persis = {}
                    for idx,fileseed in enumerate(fileseeds):
                        persisdiag,_ = readPersisFile(fileseed+str(chrnum)+'_persisdiagram.txt')
                        persis[celltypelabels[idx],idx] = persisdiag
                    for pair in idxtocompute:
                        idx1 = pair[0]
                        persis1 = persis[celltypelabels[idx1],idx1]
                        #for key2,persis2 in persis.items():
                        idx2 = pair[1]
                        persis2 = persis[celltypelabels[idx2],idx2]
                        if bndists[idx1,idx2,chrnum-1] != 0 or idx1 == idx2: continue
                        #dist = computeBottleneckDist(persis1,persis2)
                        dist = computeWassersteinDist(persis1,persis2)
                        bndists[idx1,idx2,chrnum-1] = dist
                        bndists[idx2,idx1,chrnum-1] = dist
                        print('done with',celltypelabels[idx1],'vs',celltypelabels[idx2],'on chr',str(chrnum))
                        writeBottleneckDisttoFile(bndists,celltypelabels,outputloc+'_wassersteindists.txt')
                    #print('done with chr'+str(chrnum))
                #writeBottleneckDisttoFile(bndists,celltypelabels,outputloc+'_wassersteindists.txt')
        #only plot bottleneck dists by chromosome if only 1 cell type is input, otherwise you want the full heatmap
        bnbychr = plotBottleneckByChr(bndists, outputloc+'_bottleneck_bychr.png')
        sys.exit()
        if len(genefilename) > 0:
            genelocs = readGeneFile(genefilename)
            plotGenesVBottleneck(bnbychr,genelocs,outputloc+'_gogenesVbottleneck.png')
        avgbndists = bndists.mean(2)
        plotHeatMap(avgbndists,celltypelabels,outputloc+'_avgwasserstein_heatmap.pdf')
    #print(bottleneckdists)
    #print('average bottleneck dist =',np.mean(bottleneckdists))

if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument('-i', nargs='+',help='Filenames to analyze. If just 1, plots and saves persistence and barcode diagrams. If >=2, computes bottleneck distances')
    parser.add_argument('-l', nargs='+', default=[], help='Names of cell types, in same order as filenames')
    parser.add_argument('-b', default='', nargs='+', help='If bottleneck distances have been precomputed, input filename')
    parser.add_argument('-gf',default='',help='Gene file if plotting bottleneck dists vs genes')
    parser.add_argument('-o',type=str, help='Location of output files')

    args=parser.parse_args()
    main(args.i,args.l,args.b,args.gf,args.o)
