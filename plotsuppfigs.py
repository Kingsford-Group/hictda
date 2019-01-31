import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import csv
import numpy as np
import glob

def readPersisFile(filename,dimwanted):
    persisdiag = []
    fullpersisinfo = []
    with open(filename,'rt') as f:
        freader = csv.reader(f,delimiter='\t')
        for line in freader:
            if int(line[0]) not in dimwanted: continue
            persisdiag.append((int(line[0]),(float(line[1]),float(line[2]))))
    return persisdiag

def plotAllNullBarcode(persis_esc,persis_mes,persis_cp,persis_cm,persis_fh,filename):
    if 'edge' in filename:
        ncol = '#7a5195'
    elif 'rand' in filename:
        ncol = '#ffa600'
    elif 'dist' in filename:
        ncol = '#ef5675'
    else:
        ncol = '#003f5c'
    fig,axs = plt.subplots(10,10,figsize=(30,30))
    for chrnum in range(13,23):
        pltrow = chrnum-13
        pltcol = 0
        idxmax = max([max([len(persis[0]), len(persis[1])]) for persis in [persis_esc[chrnum],persis_mes[chrnum],persis_cp[chrnum],persis_cm[chrnum],persis_fh[chrnum]]])
        for persis in persis_esc[chrnum][:2]:
            maxy = len(persis)
            # make sure persis is ordered first by all H0 (ordered by death time), then all H1
            persissorted = sorted(persis,key=lambda x: x[0])
            for idx,bar in enumerate(persissorted):
                if bar[1][1] == np.inf:
                    axs[pltrow,pltcol].plot([0,1],[maxy-idx,maxy-idx],c=ncol)
                else:
                    axs[pltrow,pltcol].plot(bar[1],[maxy-idx,maxy-idx],c=ncol)
            if pltrow < 9: axs[pltrow,pltcol].set_xticklabels([])
            if pltcol > 0: axs[pltrow,pltcol].set_yticklabels([])
            pltcol+=1
            plt.axis([0,1,0,idxmax])
            plt.yticks(fontsize=6)
            plt.xticks(fontsize=6)
        for persis in persis_mes[chrnum][:2]:
            maxy = len(persis)
            # make sure persis is ordered first by all H0 (ordered by death time), then all H1
            persissorted = sorted(persis,key=lambda x: x[0])
            for idx,bar in enumerate(persissorted):
                if bar[1][1] == np.inf:
                    axs[pltrow,pltcol].plot([0,1],[maxy-idx,maxy-idx],c=ncol)
                else:
                    axs[pltrow,pltcol].plot(bar[1],[maxy-idx,maxy-idx],c=ncol)
            if pltrow < 9: axs[pltrow,pltcol].set_xticklabels([])
            if pltcol > 0: axs[pltrow,pltcol].set_yticklabels([])
            pltcol+=1
            plt.axis([0,1,0,idxmax])
            plt.yticks(fontsize=6)
            plt.xticks(fontsize=6)
        for persis in persis_cp[chrnum][:2]:
            maxy = len(persis)
            # make sure persis is ordered first by all H0 (ordered by death time), then all H1
            persissorted = sorted(persis,key=lambda x: x[0])
            for idx,bar in enumerate(persissorted):
                if bar[1][1] == np.inf:
                    axs[pltrow,pltcol].plot([0,1],[maxy-idx,maxy-idx],c=ncol)
                else:
                    axs[pltrow,pltcol].plot(bar[1],[maxy-idx,maxy-idx],c=ncol)
            plt.axis([0,1,0,idxmax])
            if pltrow < 9: axs[pltrow,pltcol].set_xticklabels([])
            if pltcol > 0: axs[pltrow,pltcol].set_yticklabels([])
            plt.yticks(fontsize=6)
            plt.xticks(fontsize=6)
            pltcol+=1
        for persis in persis_cm[chrnum][:2]:
            maxy = len(persis)
            # make sure persis is ordered first by all H0 (ordered by death time), then all H1
            persissorted = sorted(persis,key=lambda x: x[0])
            for idx,bar in enumerate(persissorted):
                if bar[1][1] == np.inf:
                    axs[pltrow,pltcol].plot([0,1],[maxy-idx,maxy-idx],c=ncol)
                else:
                    axs[pltrow,pltcol].plot(bar[1],[maxy-idx,maxy-idx],c=ncol)
            plt.axis([0,1,0,idxmax])
            if pltrow < 9: axs[pltrow,pltcol].set_xticklabels([])
            if pltcol > 0: axs[pltrow,pltcol].set_yticklabels([])
            plt.yticks(fontsize=6)
            plt.xticks(fontsize=6)
            pltcol+=1
        for persis in persis_fh[chrnum][:2]:
            maxy = len(persis)
            # make sure persis is ordered first by all H0 (ordered by death time), then all H1
            persissorted = sorted(persis,key=lambda x: x[0])
            for idx,bar in enumerate(persissorted):
                if bar[1][1] == np.inf:
                    axs[pltrow,pltcol].plot([0,1],[maxy-idx,maxy-idx],c=ncol)
                else:
                    axs[pltrow,pltcol].plot(bar[1],[maxy-idx,maxy-idx],c=ncol)
            plt.axis([0,1,0,idxmax])
            if pltrow < 9: axs[pltrow,pltcol].set_xticklabels([])
            if pltcol > 0: axs[pltrow,pltcol].set_yticklabels([])
            plt.yticks(fontsize=6)
            plt.xticks(fontsize=6)
            pltcol+=1
        print('done plotting figs for chr',chrnum)
    #plt.yticks(fontsize=12)
    #plt.xticks(fontsize=12)
    #plt.xlabel('Radius', fontsize=12)
    #plt.legend(['H0','H1'],frameon=False)
    plt.savefig(filename,bbox_inches='tight')
    plt.close()
    print('saved barcode plot to',filename)

def plotAllHiCBarcodes(persisdict,filename):
    colors=['#fc8d59', '#91bfdb']
    fig,axs = plt.subplots(6,4,figsize=(15,20))
    figidx = 0
    for chrnum in range(1,23):
        pltcol = (chrnum-1) % 4
        pltrow = int(np.floor(chrnum-1)/4)
        persis = persisdict[chrnum]
        maxy = len(persis)
        # make sure persis is ordered first by all H0 (ordered by death time), then all H1
        persissorted = sorted(persis,key=lambda x: x[0])
        for idx,bar in enumerate(persissorted):
            if bar[1][1] == np.inf:
                axs[pltrow,pltcol].plot([0,1],[maxy-idx,maxy-idx],c=colors[bar[0]])
            else:
                axs[pltrow,pltcol].plot(bar[1],[maxy-idx,maxy-idx],c=colors[bar[0]])
        plt.axis([0,1,0,idx])
        if pltrow < 4 or (pltrow==4 and pltcol < 2): axs[pltrow,pltcol].set_xticklabels([])
        axs[pltrow,pltcol].set_title('chr'+str(chrnum))
    plt.savefig(filename,bbox_inches='tight')
    plt.close()
    print('saved barcode plot to',filename)


# need dictionary persis_esc[chrnum] = [perm0,perm1,perm2,...]
permtype = 'edge'
persis_esc,persis_mes,persis_cp,persis_cm,persis_fh = {},{},{},{},{}
for chrnum in range(13,23):
    escfiles = glob.glob('TDAresults/RUES2_ESC_100kb_chr'+str(chrnum)+'_'+permtype+'perm*_persisdiagram.txt')
    mesfiles = glob.glob('TDAresults/RUES2_MES_100kb_chr'+str(chrnum)+'_'+permtype+'perm*_persisdiagram.txt')
    cpfiles = glob.glob('TDAresults/RUES2_CP_100kb_chr'+str(chrnum)+'_'+permtype+'perm*_persisdiagram.txt')
    cmfiles = glob.glob('TDAresults/RUES2_CM_100kb_chr'+str(chrnum)+'_'+permtype+'perm*_persisdiagram.txt')
    fhfiles = glob.glob('TDAresults/RUES2_FetalHeart_100kb_chr'+str(chrnum)+'_'+permtype+'perm*_persisdiagram.txt')
    persis_esc[chrnum] = [readPersisFile(escfiles[0],[1])]
    persis_esc[chrnum].append(readPersisFile(escfiles[1],[1]))
    persis_mes[chrnum] = [readPersisFile(mesfiles[0],[1])]
    persis_mes[chrnum].append(readPersisFile(mesfiles[1],[1]))
    persis_cp[chrnum] = [readPersisFile(cpfiles[0],[1])]
    persis_cp[chrnum].append(readPersisFile(cpfiles[1],[1]))
    persis_cm[chrnum] = [readPersisFile(cmfiles[0],[1])]
    persis_cm[chrnum].append(readPersisFile(cmfiles[1],[1]))
    persis_fh[chrnum] = [readPersisFile(fhfiles[0],[1])]
    persis_fh[chrnum].append(readPersisFile(fhfiles[1],[1]))
plotAllNullBarcode(persis_esc,persis_mes,persis_cp,persis_cm,persis_fh,'TDAresults/analysis/'+permtype+'perm_allbarcodes.pdf')
celltypes = ['RUES2_ESC','RUES2_MES','RUES2_CP','RUES2_CM','RUES2_FetalHeart','WTC11_PSC','WTC11_MES','WTC11_CP','WTC11_CM','H1_ESC','H1_ME','H1_MS','H1_NP','H1_TB']
#celltype = 'RUES2_ESC'
for celltype in celltypes:
    persisdict = {}
    for chrnum in range(1,23):
        filename = glob.glob('TDAresults/'+celltype+'_100kb_chr'+str(chrnum)+'_persisdiagram.txt')
        persisdict[chrnum] = readPersisFile(filename[0],[0,1])
plotAllHiCBarcodes(persisdict,'TDAresults/analysis/'+celltype+'_allchr_barcodes.pdf')
