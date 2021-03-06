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

nruns=10
nroots=1
tags=["rhoB0.02","rhoB0.04","rhoB0.06","rhoB0.08","rhoB0.1","rhoB0.12","rhoB0.14","rhoB0.16","rhoB0.18","rhoB0.2"];
figname='../figs/mpi_vs_rhoB_BQvar_bose.pdf'
acceptance=[0.01,.2,.4,.5,.6,.8,1.0]
rhoB=[0.02,0.04,.06,.08,.1,.12,.14,.16,.18,.2]

"""
["Acc0.2","Acc0.4","Acc0.5","Acc0.6","Acc0.8","Acc0","Acc1"]
["eta0.2","eta0.6","eta1.0"]
["Omega50","Omega100","Omega250"]
["rhoB0.05","rhoB0.1","rhoB0.15","rhoB0.2"]
["T110","T125","T140","T155"]
"""

colors=['r','orange','y','g','b','purple','k']
roots=[27.0] #[7.7,11.5,19.6,27.0,39.0,62.4,200.0]
stderrK=[]
stderrS=[]
stderrKq=[]
stderrSq=[]
stderrKb=[]
stderrSb=[]
stderrKpi=[]
stderrSpi=[]

for n in range(nroots):
    stderrK.append([])
    stderrS.append([])
    Ssigma_avg=[]
    Ksigma2_avg=[]
    stderrKq.append([])
    stderrSq.append([])
    Ssigmaq_avg=[]
    Ksigma2q_avg=[]
    stderrKb.append([])
    stderrSb.append([])
    Ssigmab_avg=[]
    Ksigma2b_avg=[]
    stderrKpi.append([])
    stderrSpi.append([])
    Ssigmapi_avg=[]
    Ksigma2pi_avg=[]
    i=0
    for tag in tags:
        sumS=0
        sumK=0
        sumSq=0
        sumKq=0
        sumSb=0
        sumKb=0
        sumSpi=0
        sumKpi=0
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
        bbar=[]
        sigma2b=[]
        Ssigmab=[]
        Ksigma2b=[]
        Skellamb=[]
        pibar=[]
        sigma2pi=[]
        Ssigmapi=[]
        Ksigma2pi=[]
        Skellampi=[]
        Ssigma_avg.append(0)
        Ksigma2_avg.append(0)
        Ssigmaq_avg.append(0)
        Ksigma2q_avg.append(0)
        Ssigmab_avg.append(0)
        Ksigma2b_avg.append(0)
        Ssigmapi_avg.append(0)
        Ksigma2pi_avg.append(0)
        file="../data/BQvar"+tag+".dat";
        mydata = np.loadtxt(file,skiprows=1,unpack=True)
        l=len(mydata[0])
        for run in range(nruns):
            Omega.append(mydata[0][l-nruns+run])

            pbar.append(mydata[5][l-nruns+run])
            sigma2p.append(mydata[6][l-nruns+run])
            Ssigmap.append(mydata[7][l-nruns+run])
            Ksigma2p.append(mydata[8][l-nruns+run])

            qbar.append(mydata[9][l-nruns+run])
            sigma2q.append(mydata[10][l-nruns+run])
            Ssigmaq.append(mydata[11][l-nruns+run])
            Ksigma2q.append(mydata[12][l-nruns+run])

            bbar.append(mydata[13][l-nruns+run])
            sigma2b.append(mydata[14][l-nruns+run])
            Ssigmab.append(mydata[15][l-nruns+run])
            Ksigma2b.append(mydata[16][l-nruns+run])

            pibar.append(mydata[21][l-nruns+run])
            sigma2pi.append(mydata[22][l-nruns+run])
            Ssigmapi.append(mydata[23][l-nruns+run])
            Ksigma2pi.append(mydata[24][l-nruns+run])

            if qbar[run]!=0:
                #Skellamp.append(sigma2p[run]/(pbar[run]*Omega[run]))
                Skellamq.append(sigma2q[run]/(qbar[run]*Omega[run]))

                #Ssigma_avg[i]+=Ssigmap[run]
                #Ksigma2_avg[i]+=Ksigma2p[run]

                Ssigmaq_avg[i]+=Ssigmaq[run]
                Ksigma2q_avg[i]+=Ksigma2q[run]

                Ssigmab_avg[i]+=Ssigmab[run]
                Ksigma2b_avg[i]+=Ksigma2b[run]

                Ssigmapi_avg[i]+=Ssigmapi[run]
                Ksigma2pi_avg[i]+=Ksigma2pi[run]

            else:
                print(tag,i,run,qbar[run])
                exit(1)

        Ssigma_avg[i]*=1/nruns
        Ksigma2_avg[i]*=1/nruns
        Ssigmaq_avg[i]*=1/nruns
        Ksigma2q_avg[i]*=1/nruns
        Ssigmab_avg[i]*=1/nruns
        Ksigma2b_avg[i]*=1/nruns
        Ssigmapi_avg[i]*=1/nruns
        Ksigma2pi_avg[i]*=1/nruns

        for run in range(nruns):
            sumS+=(Ssigmap[run]-Ssigma_avg[i])**2
            sumK+=(Ksigma2p[run]-Ksigma2_avg[i])**2
            sumSq+=(Ssigmaq[run]-Ssigmaq_avg[i])**2
            sumKq+=(Ksigma2q[run]-Ksigma2q_avg[i])**2
            sumSb+=(Ssigmab[run]-Ssigmab_avg[i])**2
            sumKb+=(Ksigma2b[run]-Ksigma2b_avg[i])**2
            sumSpi+=(Ssigmapi[run]-Ssigmapi_avg[i])**2
            sumKpi+=(Ksigma2pi[run]-Ksigma2pi_avg[i])**2

        stderrS[n].append((1/nruns)*np.sqrt(sumS))
        stderrK[n].append((1/nruns)*np.sqrt(sumK))
        stderrSq[n].append((1/nruns)*np.sqrt(sumSq))
        stderrKq[n].append((1/nruns)*np.sqrt(sumKq))
        stderrSb[n].append((1/nruns)*np.sqrt(sumSb))
        stderrKb[n].append((1/nruns)*np.sqrt(sumKb))
        stderrSpi[n].append((1/nruns)*np.sqrt(sumSpi))
        stderrKpi[n].append((1/nruns)*np.sqrt(sumKpi))
        i+=1

        #print("Ksigma2 error =",stderrK,", Ssigma error =",stderrS)
    #plt.errorbar(rhoB,Ssigmaq_avg,stderrSq[n],linestyle='-',linewidth=2,color='r',markersize=8, marker='^', markerfacecolor=None, markeredgecolor=None,label='$S\sigma$ (Q)')
    #plt.errorbar(rhoB,Ksigma2q_avg,stderrKq[n],linestyle='--',linewidth=2,color='r',markersize=10, marker='s', markerfacecolor=None, markeredgecolor=None,label='$\kappa\sigma^2$ (Q)')

    #plt.errorbar(rhoB,Ssigmab_avg,stderrSb[n],linestyle='-',linewidth=2,color='b',markersize=8, marker='^', markerfacecolor=None, markeredgecolor=None,label='$S\sigma$ (B)')
    #plt.errorbar(rhoB,Ksigma2b_avg,stderrKb[n],linestyle='--',linewidth=2,color='b',markersize=10, marker='s', markerfacecolor=None, markeredgecolor=None,label='$\kappa\sigma^2$ (B)')

    plt.errorbar(rhoB,Ssigmapi_avg,stderrSpi[n],linestyle='-',linewidth=2,color='purple',markersize=8, marker='^', markerfacecolor=None, markeredgecolor=None,label='$S\sigma$ (B)')
    plt.errorbar(rhoB,Ksigma2pi_avg,stderrKpi[n],linestyle='--',linewidth=2,color='purple',markersize=10, marker='s', markerfacecolor=None, markeredgecolor=None,label='$\kappa\sigma^2$ (B)')

stderrK=[]
stderrS=[]
stderrKq=[]
stderrSq=[]
stderrKb=[]
stderrSb=[]
stderrKpi=[]
stderrSpi=[]

for n in range(nroots):
    stderrK.append([])
    stderrS.append([])
    Ssigma_avg=[]
    Ksigma2_avg=[]
    stderrKq.append([])
    stderrSq.append([])
    Ssigmaq_avg=[]
    Ksigma2q_avg=[]
    stderrKb.append([])
    stderrSb.append([])
    Ssigmab_avg=[]
    Ksigma2b_avg=[]
    stderrKpi.append([])
    stderrSpi.append([])
    Ssigmapi_avg=[]
    Ksigma2pi_avg=[]
    i=0
    for tag in tags:
        sumS=0
        sumK=0
        sumSq=0
        sumKq=0
        sumSb=0
        sumKb=0
        sumSpi=0
        sumKpi=0
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
        bbar=[]
        sigma2b=[]
        Ssigmab=[]
        Ksigma2b=[]
        Skellamb=[]
        pibar=[]
        sigma2pi=[]
        Ssigmapi=[]
        Ksigma2pi=[]
        Skellampi=[]
        Ssigma_avg.append(0)
        Ksigma2_avg.append(0)
        Ssigmaq_avg.append(0)
        Ksigma2q_avg.append(0)
        Ssigmab_avg.append(0)
        Ksigma2b_avg.append(0)
        Ssigmapi_avg.append(0)
        Ksigma2pi_avg.append(0)
        file="../data/boseBQvar"+tag+".dat";
        mydata = np.loadtxt(file,skiprows=1,unpack=True)
        l=len(mydata[0])
        for run in range(nruns):
            Omega.append(mydata[0][l-nruns+run])

            pbar.append(mydata[5][l-nruns+run])
            sigma2p.append(mydata[6][l-nruns+run])
            Ssigmap.append(mydata[7][l-nruns+run])
            Ksigma2p.append(mydata[8][l-nruns+run])

            qbar.append(mydata[9][l-nruns+run])
            sigma2q.append(mydata[10][l-nruns+run])
            Ssigmaq.append(mydata[11][l-nruns+run])
            Ksigma2q.append(mydata[12][l-nruns+run])

            bbar.append(mydata[13][l-nruns+run])
            sigma2b.append(mydata[14][l-nruns+run])
            Ssigmab.append(mydata[15][l-nruns+run])
            Ksigma2b.append(mydata[16][l-nruns+run])

            pibar.append(mydata[21][l-nruns+run])
            sigma2pi.append(mydata[22][l-nruns+run])
            Ssigmapi.append(mydata[23][l-nruns+run])
            Ksigma2pi.append(mydata[24][l-nruns+run])

            if qbar[run]!=0:
                #Skellamp.append(sigma2p[run]/(pbar[run]*Omega[run]))
                Skellamq.append(sigma2q[run]/(qbar[run]*Omega[run]))

                #Ssigma_avg[i]+=Ssigmap[run]
                #Ksigma2_avg[i]+=Ksigma2p[run]

                Ssigmaq_avg[i]+=Ssigmaq[run]
                Ksigma2q_avg[i]+=Ksigma2q[run]

                Ssigmab_avg[i]+=Ssigmab[run]
                Ksigma2b_avg[i]+=Ksigma2b[run]

                Ssigmapi_avg[i]+=Ssigmapi[run]
                Ksigma2pi_avg[i]+=Ksigma2pi[run]

            else:
                print(tag,i,run,qbar[run])
                exit(1)

        Ssigma_avg[i]*=1/nruns
        Ksigma2_avg[i]*=1/nruns
        Ssigmaq_avg[i]*=1/nruns
        Ksigma2q_avg[i]*=1/nruns
        Ssigmab_avg[i]*=1/nruns
        Ksigma2b_avg[i]*=1/nruns
        Ssigmapi_avg[i]*=1/nruns
        Ksigma2pi_avg[i]*=1/nruns

        for run in range(nruns):
            sumS+=(Ssigmap[run]-Ssigma_avg[i])**2
            sumK+=(Ksigma2p[run]-Ksigma2_avg[i])**2
            sumSq+=(Ssigmaq[run]-Ssigmaq_avg[i])**2
            sumKq+=(Ksigma2q[run]-Ksigma2q_avg[i])**2
            sumSb+=(Ssigmab[run]-Ssigmab_avg[i])**2
            sumKb+=(Ksigma2b[run]-Ksigma2b_avg[i])**2
            sumSpi+=(Ssigmapi[run]-Ssigmapi_avg[i])**2
            sumKpi+=(Ksigma2pi[run]-Ksigma2pi_avg[i])**2

        stderrS[n].append((1/nruns)*np.sqrt(sumS))
        stderrK[n].append((1/nruns)*np.sqrt(sumK))
        stderrSq[n].append((1/nruns)*np.sqrt(sumSq))
        stderrKq[n].append((1/nruns)*np.sqrt(sumKq))
        stderrSb[n].append((1/nruns)*np.sqrt(sumSb))
        stderrKb[n].append((1/nruns)*np.sqrt(sumKb))
        stderrSpi[n].append((1/nruns)*np.sqrt(sumSpi))
        stderrKpi[n].append((1/nruns)*np.sqrt(sumKpi))
        i+=1

        #print("Ksigma2 error =",stderrK,", Ssigma error =",stderrS)
    #plt.errorbar(rhoB,Ssigmaq_avg,stderrSq[n],linestyle='-',linewidth=2,color='orange',markersize=8, marker='^', markerfacecolor=None, markeredgecolor=None,label='$S\sigma$ (Q)')
    #plt.errorbar(rhoB,Ksigma2q_avg,stderrKq[n],linestyle='--',linewidth=2,color='orange',markersize=10, marker='s', markerfacecolor=None, markeredgecolor=None,label='$\kappa\sigma^2$ (Q)')

    #plt.errorbar(rhoB,Ssigmab_avg,stderrSb[n],linestyle='-',linewidth=2,color='c',markersize=8, marker='^', markerfacecolor=None, markeredgecolor=None,label='$S\sigma$ (B)')
    #plt.errorbar(rhoB,Ksigma2b_avg,stderrKb[n],linestyle='--',linewidth=2,color='c',markersize=10, marker='s', markerfacecolor=None, markeredgecolor=None,label='$\kappa\sigma^2$ (B)')

    plt.errorbar(rhoB,Ssigmapi_avg,stderrSpi[n],linestyle='-',linewidth=2,color='plum',markersize=8, marker='^', markerfacecolor=None, markeredgecolor=None,label='$S\sigma$ (B)')
    plt.errorbar(rhoB,Ksigma2pi_avg,stderrKpi[n],linestyle='--',linewidth=2,color='plum',markersize=10, marker='s', markerfacecolor=None, markeredgecolor=None,label='$\kappa\sigma^2$ (B)')

ax.tick_params(axis='both', which='major', labelsize=14)
ax.set_xticks(np.arange(0,.21,.04), minor=False)
ax.set_xticklabels(np.arange(0,.21,.04), minor=False, family='serif')
ax.set_xticks(np.arange(0,.21,.02), minor=True)
ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%d'))
ax.xaxis.set_major_formatter(sformatter)
plt.xlim(0,.21)

ax.set_yticks(np.arange(-1,1.5,0.5), minor=False)
ax.set_yticklabels(np.arange(-1,1.5,0.5), minor=False, family='serif')
ax.set_yticks(np.arange(-1,1.5,0.1), minor=True)
plt.ylim(-.2,1.05)
ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%1f'))
ax.yaxis.set_major_formatter(sformatter)

#ax.legend(loc=(0.1,0.1));

plt.xlabel('$\\rho_B$ $(fm^{-3})$',fontsize=18 , weight='normal')
plt.ylabel('$S\sigma$,  $\kappa\sigma^2$', fontsize=24, weight='normal')
plt.savefig(figname,format='pdf',bbox_inches = "tight")
os.system('xdg-open '+figname)



#plt.show()
quit()
