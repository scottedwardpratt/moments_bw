import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import os
from pylab import *
from matplotlib import ticker
from matplotlib.ticker import ScalarFormatter
sformatter=ScalarFormatter(useOffset=True,useMathText=True)
sformatter.set_scientific(True)
sformatter.set_powerlimits((-2,3))

#plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

font = {'family' : 'serif',
        'weight' : 'normal',
        'size'   : 14}
plt.rc('font', **font)
plt.rc('text', usetex=False)
plt.figure(figsize=(6,5))
fig = plt.figure(1)
ax = fig.add_axes([0.14,0.125,0.82,0.85])


nroots=7
nruns=10
tags=["eta0.3","eta0.4","eta0.5","eta1.0"]
colors=['r','g','b','k']
roots=[7.7,11.5,19.6,27.0,39.0,62.4,200.0]
n=0
stderrK=[]
stderrS=[]
stderrKq=[]
stderrSq=[]

for tag in tags:
    stderrK.append([])
    stderrS.append([])
    Ssigma_avg=[]
    Ksigma2_avg=[]
    stderrKq.append([])
    stderrSq.append([])
    Ssigmaq_avg=[]
    Ksigma2q_avg=[]
    for i in range(nroots):
        sumS=0
        sumK=0
        sumSq=0
        sumKq=0
        Omega=[]
        pbar=[]
        sigma2p=[]
        Ssigmap=[]
        Ksigma2p=[]
        Skellamp=[]
        qbar=[]
        sigma2q=[]
        Ssigmaq=[]
        Ksigma2q=[]
        Skellamq=[]
        Ssigma_avg.append(0)
        Ksigma2_avg.append(0)
        Ssigmaq_avg.append(0)
        Ksigma2q_avg.append(0)
        file="../momentsdata/fixedff"+str(i)+tag+".dat";
        mydata = np.loadtxt(file,skiprows=1,unpack=True)
        for run in range(nruns):
            Omega.append(mydata[0][run])

            pbar.append(mydata[5][run])
            sigma2p.append(mydata[6][run])
            Ssigmap.append(mydata[7][run])
            Ksigma2p.append(mydata[8][run])

            qbar.append(mydata[9][run])
            sigma2q.append(mydata[10][run])
            Ssigmaq.append(mydata[11][run])
            Ksigma2q.append(mydata[12][run])

            Skellamp.append(sigma2p[run]/(pbar[run]*Omega[run]))
            Skellamq.append(sigma2q[run]/(qbar[run]*Omega[run]))

            Ssigma_avg[i]+=Ssigmap[run]*Skellamp[run]
            Ksigma2_avg[i]+=Ksigma2p[run]

            Ssigmaq_avg[i]+=Ssigmaq[run]*Skellamq[run]
            Ksigma2q_avg[i]+=Ksigma2q[run]

        Ssigma_avg[i]*=1/nruns
        Ksigma2_avg[i]*=1/nruns
        Ssigmaq_avg[i]*=1/nruns
        Ksigma2q_avg[i]*=1/nruns

        for run in range(nruns):
            sumS+=(Ssigmap[run]*Skellamp[run]-Ssigma_avg[i])**2
            sumK+=(Ksigma2p[run]-Ksigma2_avg[i])**2
            sumSq+=(Ssigmaq[run]*Skellamq[run]-Ssigmaq_avg[i])**2
            sumKq+=(Ksigma2q[run]-Ksigma2q_avg[i])**2

        stderrS[n].append((1/nroots)*np.sqrt(sumS))
        stderrK[n].append((1/nroots)*np.sqrt(sumK))
        stderrSq[n].append((1/nroots)*np.sqrt(sumSq))
        stderrKq[n].append((1/nroots)*np.sqrt(sumKq))

        #print("Ksigma2 error =",stderrK,", Ssigma error =",stderrS)
    plt.errorbar(roots,Ssigma_avg,stderrS[n],linestyle='-',linewidth=2,color=colors[n],markersize=8, marker='s', markerfacecolor=None, markeredgecolor=None,label=tag+': $C_3/C_1$')
    plt.errorbar(roots,Ksigma2_avg,stderrK[n],linestyle='--',linewidth=2,color=colors[n],markersize=10, marker='^', markerfacecolor=None, markeredgecolor=None,label=tag+': $C_4/C_2$')

    #plt.errorbar(roots,Ssigmap*Skellamp,stderrS[n],linestyle='-',linewidth=2,color=colors[n],markersize=8, marker='s', markerfacecolor=None, markeredgecolor=None,label=tag+': $C_3/C_1$')
    #plt.errorbar(roots,Ksigma2p,stderrK[n],linestyle='--',linewidth=2,color=colors[n],markersize=10, marker='^', markerfacecolor=None, markeredgecolor=None,label=tag+': $C_4/C_2$')

    n+=1

ax.tick_params(axis='both', which='major', labelsize=14)
ax.set_xticks(np.arange(0,250,50), minor=False)
ax.set_xticklabels(np.arange(0,250,25), minor=False, family='serif')
ax.set_xticks(np.arange(0,250,50), minor=True)
ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%d'))
ax.xaxis.set_major_formatter(sformatter)
plt.xlim(0,210)

ax.set_yticks(np.arange(-1,1.5,0.5), minor=False)
ax.set_yticklabels(np.arange(-1,1.5,0.5), minor=False, family='serif')
ax.set_yticks(np.arange(-1,1.5,0.05), minor=True)
plt.ylim(0.0,1.05)
ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%1f'))
ax.yaxis.set_major_formatter(sformatter)

ax.legend(loc=(0.52,0.1));

plt.xlabel('$\sqrt{s}_{NN}$ (GeV)',fontsize=18 , weight='normal')
plt.ylabel('$S\sigma$,  $\kappa\sigma^2 (Q)$', fontsize=24, weight='normal')
plt.savefig('moments_bw_vsroots_fixedff.pdf',format='pdf')
os.system('xdg-open moments_bw_vsroots_fixedff.pdf')



#plt.show()
quit()
