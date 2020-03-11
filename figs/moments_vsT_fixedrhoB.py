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
rhoB=0.0373777 #[0.143511,0.0942362,0.0525105,0.0373777,0.0255775,0.0159063,0.00497411]
T=[143.45,149.53,153.4,154.45,155.08,155.44,155.66]
colors=['r','g','b','k']
n=0
stderrK=[]
stderrS=[]

for tag in tags:
    stderrK.append([])
    stderrS.append([])
    Ssigma_avg=[]
    Ksigma2_avg=[]
    for i in range(nroots):
        sumS=0
        sumK=0
        Omega=[]
        pbar=[]
        sigma2p=[]
        Ssigmap=[]
        Ksigma2p=[]
        Skellamp=[]
        Ssigma_avg.append(0)
        Ksigma2_avg.append(0)

        file="../moments_rhoB/T"+str(i)+tag+".dat";
        mydata = np.loadtxt(file,skiprows=1,unpack=True)
        for run in range(nruns):
            Omega.append(mydata[0][run])

            pbar.append(mydata[5][run])
            sigma2p.append(mydata[6][run])
            Ssigmap.append(mydata[7][run])
            Ksigma2p.append(mydata[8][run])

            Skellamp.append(sigma2p[run]/(pbar[run]*Omega[run]))

            Ssigma_avg[i]+=Ssigmap[run]*Skellamp[run]
            Ksigma2_avg[i]+=Ksigma2p[run]

        Ssigma_avg[i]*=1/nruns
        Ksigma2_avg[i]*=1/nruns

        for run in range(nruns):
            sumS+=(Ssigmap[run]*Skellamp[run]-Ssigma_avg[i])**2
            sumK+=(Ksigma2p[run]-Ksigma2_avg[i])**2

        stderrS[n].append((1/nroots)*np.sqrt(sumS))
        stderrK[n].append((1/nroots)*np.sqrt(sumK))

        #print("Ksigma2 error =",stderrK,", Ssigma error =",stderrS)

    plt.errorbar(T,Ssigma_avg,stderrS[n],linestyle='-',linewidth=2,color=colors[n],markersize=8, marker='s', markerfacecolor=None, markeredgecolor=None,label=tag+': $C_3/C_1$')
    plt.errorbar(T,Ksigma2_avg,stderrK[n],linestyle='--',linewidth=2,color=colors[n],markersize=10, marker='^', markerfacecolor=None, markeredgecolor=None,label=tag+': $C_4/C_2$')

    n+=1

ax.tick_params(axis='both', which='major', labelsize=14)
ax.set_xticks(np.arange(140,160,5), minor=False)
ax.set_xticklabels(np.arange(140,160,5), minor=False, family='serif')
ax.set_xticks(np.arange(140,160,5), minor=True)
ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%d'))
ax.xaxis.set_major_formatter(sformatter)
plt.xlim(140,160)

ax.set_yticks(np.arange(-1,1.5,0.5), minor=False)
ax.set_yticklabels(np.arange(-1,1.5,0.5), minor=False, family='serif')
ax.set_yticks(np.arange(-1,1.5,0.05), minor=True)
plt.ylim(0.0,1.05)
ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%1f'))
ax.yaxis.set_major_formatter(sformatter)

ax.legend(loc=(0.1,0.05));

plt.xlabel('T (MeV)',fontsize=18 , weight='normal')
plt.ylabel('$S\sigma$,  $\kappa\sigma^2$', fontsize=24, weight='normal')
plt.savefig('moments_bw_vsT_fixedrhoB.pdf',format='pdf')
os.system('xdg-open moments_bw_vsT_fixedrhoB.pdf')



#plt.show()
quit()
