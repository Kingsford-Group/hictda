import csv
import glob
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

#need to read in persistence diagrams
def readPersisFile(filename):
    persisdiag = []
    with open(filename,'rt') as f:
        freader = csv.reader(f,delimiter='\t')
        for line in freader:
            if int(line[0]) != 1: continue
            persisdiag.append((float(line[1]),float(line[2])))
    return persisdiag

#plot birth times
def plotBirthTimes(origpersis,randperms,edgeperms,distperms,filename):

    origtimes,randtimes,edgetimes,disttimes = [],[],[],[]
    for chrnum in range(13,23):
        origtimes.extend([x[0] for celltype in origpersis[chrnum] for x in celltype[1]])
        if chrnum in randperms:
            randtimes.extend([x[0] for perm in randperms[chrnum] for x in perm[1]])
        else:
            print 'no random permutation data for chr',chrnum
        if chrnum in edgeperms:
            edgetimes.extend([x[0] for perm in edgeperms[chrnum] for x in perm[1]])
        else:
            print 'no edge permutation data for chr',chrnum
        if chrnum in distperms:
            disttimes.extend([x[0] for perm in distperms[chrnum] for x in perm[1]])
        else:
            print 'no distance permutation data for chr',chrnum
    print len(origtimes),len(randtimes),len(edgetimes),len(disttimes)
    fig,ax = plt.subplots()

    sns.distplot(randtimes, label="Random permutation", ax=ax, kde=True,kde_kws={'linewidth': 2}, norm_hist=True,color='#ffa600',hist_kws={"edgecolor": "none","alpha":0.1})#"histtype": "step","alpha":1,"linewidth":2})
    sns.distplot(edgetimes, label="Edge permutation", ax=ax, kde=True,kde_kws={'linewidth': 2},norm_hist=True,color='#7a5195',hist_kws={"edgecolor": "none","alpha":0.1})#"histtype": "step","alpha":1,"linewidth":2})#
    sns.distplot(disttimes, label="Linear dependence", ax=ax, kde=True,kde_kws={'linewidth': 2},norm_hist=True,color='#ef5675',hist_kws={"edgecolor": "none","alpha":0.1})#"histtype": "step","alpha":1,"linewidth":2})#
    sns.distplot(origtimes, label="Hi-C", ax=ax, kde=True,kde_kws={'linewidth': 2}, norm_hist=True,color='#003f5c',hist_kws={"edgecolor": "none","alpha":0.1})#"histtype": "step","alpha":1,"linewidth":2})#
    ax.set_xlim([0,1])
    ax.set_yticklabels([])
    plt.xlabel('Loop birth time')
    ax.legend(frameon=False,loc='upper left')
    plt.savefig(filename,bbox_inches='tight')
    print 'saved birth times histogram to',filename

#plot lifespans
def plotLifeSpans(origpersis,randperms,edgeperms,distperms,filename):
    origtimes,randtimes,edgetimes,disttimes = [],[],[],[]
    for chrnum in range(13,23):
        origtimes.extend([x[1]-x[0] for celltype in origpersis[chrnum] for x in celltype[1]])
        if chrnum in randperms:
            randtimes.extend([x[1]-x[0] for perm in randperms[chrnum] for x in perm[1]])
        else:
            print 'no random permutation data for chr',chrnum
        if chrnum in edgeperms:
            edgetimes.extend([x[1]-x[0] for perm in edgeperms[chrnum] for x in perm[1]])
        else:
            print 'no edge permutation data for chr',chrnum
        if chrnum in distperms:
            disttimes.extend([x[1]-x[0] for perm in distperms[chrnum] for x in perm[1]])
        else:
            print 'no distance permutation data for chr',chrnum
    print len(origtimes),len(randtimes),len(edgetimes),len(disttimes)
    fig,ax = plt.subplots()

    sns.distplot(randtimes, label="Random permutation", ax=ax, kde=True,kde_kws={'linewidth': 2},norm_hist=True,color = '#ffa600',hist_kws={"edgecolor": "none","alpha":0.1})#"histtype": "step","alpha":1,"linewidth":2})#
    sns.distplot(edgetimes, label="Edge permutation", ax=ax, kde=True,kde_kws={'linewidth': 2},norm_hist=True,color='#7a5195',hist_kws={"edgecolor": "none","alpha":0.1})#{"histtype": "step","alpha":1,"linewidth":2})#
    sns.distplot(disttimes, label="Linear dependence", ax=ax, kde=True,kde_kws={'linewidth': 2},norm_hist=True,color='#ef5675',hist_kws={"edgecolor": "none","alpha":0.1})#{"histtype": "step","alpha":1,"linewidth":2})#
    sns.distplot(origtimes, label="Hi-C", ax=ax, kde=True,kde_kws={'linewidth': 2},norm_hist=True,color='#003f5c',hist_kws={"edgecolor": "none","alpha":0.1})#{"histtype": "step","alpha":1,"linewidth":2})#
    ax.set_xlim([0,0.6])
    ax.set_yticklabels([])
    plt.xlabel('Loop life span')
    ax.legend(frameon=False,loc='upper right')
    plt.savefig(filename,bbox_inches='tight')
    print 'saved life spans histogram to',filename

#plot numbers of loops (ratio w/ original)
def plotLoopNumbers(origpersis,randperms,edgeperms,distperms,filename):
    randtimes,edgetimes,disttimes = [],[],[]
    for chrnum in range(13,23):
        origloopnums = {}
        for celltypes in origpersis[chrnum]:
            origloopnums[celltypes[0]] = len(celltypes[1])
        #origtimes.extend([x[1]-x[0] for celltype in origpersis[chrnum] for x in celltype])
        if chrnum in randperms:
            for perm in randperms[chrnum]:
                celltype = perm[0]
                randtimes.extend([float(len(perm[1]))/origloopnums[celltype]])
        else:
            print 'no random permutation data for chr',chrnum
        if chrnum in edgeperms:
            for perm in edgeperms[chrnum]:
                celltype = perm[0]
                edgetimes.extend([float(len(perm[1]))/origloopnums[celltype]])
            #edgetimes.extend([x[1]-x[0] for perm in edgeperms[chrnum] for x in perm])
        else:
            print 'no edge permutation data for chr',chrnum
        if chrnum in distperms:
            for perm in distperms[chrnum]:
                celltype = perm[0]
                disttimes.extend([float(len(perm[1]))/origloopnums[celltype]])
            #disttimes.extend([x[1]-x[0] for perm in distperms[chrnum] for x in perm])
        else:
            print 'no distance permutation data for chr',chrnum
    print len(randtimes),len(edgetimes),len(disttimes)
    fig,ax = plt.subplots()

    sns.distplot(randtimes, label="Random permutation", ax=ax, kde=False,norm_hist=True,color='#ffa600',hist_kws={"edgecolor": "none","alpha":0.75})
    sns.distplot(edgetimes, label="Edge permutation", ax=ax, kde=False,norm_hist=True,color='#7a5195',hist_kws={"edgecolor": "none","alpha":0.75})
    sns.distplot(disttimes, label="Linear dependence", ax=ax, kde=False,norm_hist=True,color='#ef5675',hist_kws={"edgecolor": "none","alpha":0.75})
    #sns.distplot(origtimes, label="Hi-C", ax=ax, kde=False,norm_hist=True)
    #ax.set_xlim([0,1])
    ax.set_yticklabels([])
    plt.xticks(np.arange(0,np.ceil(max(edgetimes)),step=1))
    plt.xlabel('Ratio of loop count in model to corresponding Hi-C')
    ax.legend(frameon=False,loc='upper right')
    plt.savefig(filename,bbox_inches='tight')
    print 'saved loop counts histogram to',filename


fileseed = sys.argv[1]
outputloc = sys.argv[2]
origpersis = {}
distperms = {}
randperms = {}
edgeperms = {}
for chrnum in range(13,23):
    origfilenames = glob.glob(fileseed+str(chrnum)+'_persisdiagram.txt')
    for origfile in origfilenames:
        celltype = origfile.split('_')[1]
        if chrnum in origpersis:
            origpersis[chrnum].append((celltype,readPersisFile(origfile)))
        else:
            origpersis[chrnum] = [(celltype,readPersisFile(origfile))]
    dpfilenames = glob.glob(fileseed+str(chrnum)+'_distperm*_persisdiagram.txt')
    #print dpfilenames
    if len(dpfilenames) > 0:
        for dpfile in dpfilenames:
            celltype = dpfile.split('_')[1]
            if chrnum in distperms:
                distperms[chrnum].append((celltype,readPersisFile(dpfile)))
            else:
                distperms[chrnum] = [(celltype,readPersisFile(dpfile))]
    rpfilenames = glob.glob(fileseed+str(chrnum)+'_randperm*_persisdiagram.txt')
    #print rpfilenames
    if len(rpfilenames) > 0:
        for rpfile in rpfilenames:
            celltype = rpfile.split('_')[1]
            if chrnum in randperms:
                randperms[chrnum].append((celltype,readPersisFile(rpfile)))
            else:
                randperms[chrnum] = [(celltype,readPersisFile(rpfile))]
    epfilenames = glob.glob(fileseed+str(chrnum)+'_edgeperm*_persisdiagram.txt')
    #print epfilenames
    if len(epfilenames) > 0:
        for epfile in epfilenames:
            celltype = epfile.split('_')[1]
            if chrnum in edgeperms:
                edgeperms[chrnum].append((celltype,readPersisFile(epfile)))
            else:
                edgeperms[chrnum] = [(celltype,readPersisFile(epfile))]
plotBirthTimes(origpersis,randperms,edgeperms,distperms,outputloc+'_birthtimes.pdf')
plotLifeSpans(origpersis,randperms,edgeperms,distperms, outputloc+'_lifespans.pdf')
plotLoopNumbers(origpersis,randperms,edgeperms,distperms,outputloc+'_loopcounts.pdf')
