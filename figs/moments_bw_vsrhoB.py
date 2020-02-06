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

mydata = np.loadtxt('moments_bw_Omega50.dat',skiprows=1,unpack=True)
Omega=mydata[0]
rhoB=mydata[3]
Pbar=mydata[5]
sigma2P=mydata[6]
SsigmaP=mydata[7]
Ksigma2P=mydata[8]
Qbar=mydata[9]
sigma2Q=mydata[10]
SsigmaQ=mydata[11]
Ksigma2Q=mydata[10]
#SkellamP=mydata[19]/Pbar
#SkellamQ=mydata[20]/Qbar
SkellamP=sigma2P/(Pbar*Omega)
#SkellamQ=sigma2Q/(Qbar*Omega)

plt.plot(rhoB,Ksigma2P,linestyle='-',linewidth=2,color='r',markersize=5, marker='o', markerfacecolor=None, markeredgecolor=None)
#plt.plot(rhoB,Ksigma2Q,linestyle='-',linewidth=2,color='r',markersize=5, marker='s', markerfacecolor=None, markeredgecolor=None)

plt.plot(rhoB,SsigmaP*SkellamP,linestyle='-',linewidth=2,color='b',markersize=5, marker='o', markerfacecolor=None, markeredgecolor=None)
#plt.plot(rhoB,SsigmaQ*SkellamQ,linestyle='-',linewidth=2,color='b',markersize=5, marker='s', markerfacecolor=None, markeredgecolor=None)

#plt.plot(rhoB,SkellamP,linestyle='-',linewidth=2,color='g',markersize=5, marker='o', markerfacecolor=None, markeredgecolor=None)
#plt.plot(rhoB,SkellamQ,linestyle='-',linewidth=2,color='g',markersize=5, marker='s', markerfacecolor=None, markeredgecolor=None)

ax.tick_params(axis='both', which='major', labelsize=14)
ax.set_xticks(np.arange(0,0.2,0.05), minor=False)
ax.set_xticklabels(np.arange(0,0.2,0.05), minor=False, family='serif')
ax.set_xticks(np.arange(0,0.2,0.01), minor=True)
ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%d'))
ax.xaxis.set_major_formatter(sformatter)
plt.xlim(0,0.15)

ax.set_yticks(np.arange(-1,1.5,0.5), minor=False)
ax.set_yticklabels(np.arange(-1,1.5,0.5), minor=False, family='serif')
ax.set_yticks(np.arange(-1,1.5,0.05), minor=True)
plt.ylim(-0.05,1.05)
ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%1f'))
ax.yaxis.set_major_formatter(sformatter)

plt.xlabel('$\\rho_B$ (fm$^{3}$)',fontsize=18 , weight='normal')
plt.ylabel('$S\sigma$,  $\kappa\sigma^2$', fontsize=24, weight='normal')

plt.savefig('moments_bw_vsrho.pdf',format='pdf')
os.system('xdg-open moments_bw_vsrho.pdf')
#plt.show()
quit()
