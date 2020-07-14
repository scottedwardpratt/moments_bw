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
plt.figure(figsize=(6,9))
fig = plt.figure(1)

Omega=100
Mydata1=np.loadtxt('../results/sigmaeta0.1_Omega100.dat',skiprows=1,unpack=True)
roots1=Mydata1[1]
pbar1=Mydata1[5]
sigma2p1=Mydata1[6]
Ssigmap1=Mydata1[7]
Ksigma2p1=Mydata1[8]
qbar1=Mydata1[9]
sigma2q1=Mydata1[10]
Ssigmaq1=Mydata1[11]
Ksigma2q1=Mydata1[12]
Skellamp1=sigma2p1/(pbar1*Omega)
Skellamq1=sigma2q1/(qbar1*Omega)

Mydata3=np.loadtxt('../results/sigmaeta0.3_Omega100.dat',skiprows=1,unpack=True)
roots3=Mydata3[1]
pbar3=Mydata3[5]
sigma2p3=Mydata3[6]
Ssigmap3=Mydata3[7]
Ksigma2p3=Mydata3[8]
qbar3=Mydata3[9]
sigma2q3=Mydata3[10]
Ssigmaq3=Mydata3[11]
Ksigma2q3=Mydata3[12]
Skellamp3=sigma2p3/(pbar3*Omega)
Skellamq3=sigma2q3/(qbar3*Omega)

Mydata5=np.loadtxt('../results/sigmaeta0.5_Omega100.dat',skiprows=1,unpack=True)
roots5=Mydata5[1]
pbar5=Mydata5[5]
sigma2p5=Mydata5[6]
Ssigmap5=Mydata5[7]
Ksigma2p5=Mydata5[8]
qbar5=Mydata5[9]
sigma2q5=Mydata5[10]
Ssigmaq5=Mydata5[11]
Ksigma2q5=Mydata5[12]
Skellamp5=sigma2p5/(pbar5*Omega)
Skellamq5=sigma2q5/(qbar5*Omega)

Mydata7=np.loadtxt('../results/sigmaeta0.3_Omega100.dat',skiprows=1,unpack=True)
roots7=Mydata7[1]
pbar7=Mydata7[5]
sigma2p7=Mydata7[6]
Ssigmap7=Mydata7[7]
Ksigma2p7=Mydata7[8]
qbar7=Mydata7[9]
sigma2q7=Mydata7[10]
Ssigmaq7=Mydata7[11]
Ksigma2q7=Mydata7[12]
Skellamp7=sigma2p7/(pbar7*Omega)
Skellamq7=sigma2q7/(qbar7*Omega)

#################################################################
######## LOWER PANEL protons
ax = fig.add_axes([0.15,0.12,0.84,0.43])

#plt.plot(roots,Ssigmap*Skellamp,linestyle='-',linewidth=2,color='r',markersize=8, marker='s', markerfacecolor=None, markeredgecolor=None,label='ETA=0.3: $C_3/C_1$')
plt.plot(roots1,Ssigmap1,linestyle='--',linewidth=2,color='r',markersize=10, marker='s', markerfacecolor='r', markeredgecolor='r',label='$\sigma_\eta=0.1$')
plt.plot(roots3,Ssigmap3,linestyle='--',linewidth=2,color='k',markersize=10, marker='s', markerfacecolor='k', markeredgecolor='k',label='$\sigma_\eta=0.3$')
plt.plot(roots5,Ssigmap5,linestyle='--',linewidth=2,color='g',markersize=10, marker='s', markerfacecolor='g', markeredgecolor='g',label='$\sigma_\eta=0.5$')
plt.plot(roots7,Ssigmap7,linestyle='--',linewidth=2,color='b',markersize=10, marker='s', markerfacecolor='b', markeredgecolor='b',label='$\sigma_\eta=0.7$')

ax.tick_params(axis='both', which='major', labelsize=14)
ax.set_xticks(np.arange(0,250,50), minor=False)
ax.set_xticklabels(np.arange(0,250,25), minor=False, family='serif')
ax.set_xticks(np.arange(0,250,50), minor=True)
ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%d'))
ax.xaxis.set_major_formatter(sformatter)
plt.xlim(0,210)

ax.set_yticks(np.arange(-1,2.5,0.2), minor=False)
ax.set_yticklabels(np.arange(-1,2.5,0.2), minor=False, family='serif')
ax.set_yticks(np.arange(-1,2.5,0.05), minor=True)
plt.ylim(0.0,1.05)
ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%1f'))
ax.yaxis.set_major_formatter(sformatter)

ax.legend(loc=(0.6,0.38),fontsize=18);

plt.xlabel('$\sqrt{s}_{NN}$ (GeV)',fontsize=22 , weight='normal')
plt.ylabel('$S\sigma=C_3/C_1$', fontsize=22, weight='normal')
text(200,0.95,'net protons',fontsize=22,ha='right')

######## Upper Panel charge
ax = fig.add_axes([0.15,0.55,0.84,0.43])

#plt.plot(roots,Ssigmap*Skellamp,linestyle='-',linewidth=2,color='r',markersize=8, marker='s', markerfacecolor=None, markeredgecolor=None,label='ETA=0.3: $C_3/C_1$')
plt.plot(roots1,Ssigmaq1,linestyle='--',linewidth=2,color='r',markersize=10, marker='s', markerfacecolor='r', markeredgecolor='r',label='$\sigma_\eta=0.1$')
plt.plot(roots3,Ssigmaq3,linestyle='--',linewidth=2,color='k',markersize=10, marker='s', markerfacecolor='k', markeredgecolor='k',label='$\sigma_\eta=0.3$')
plt.plot(roots5,Ssigmaq5,linestyle='--',linewidth=2,color='g',markersize=10, marker='s', markerfacecolor='g', markeredgecolor='g',label='$\sigma_\eta=0.5$')
plt.plot(roots7,Ssigmaq7,linestyle='--',linewidth=2,color='b',markersize=10, marker='s', markerfacecolor='b', markeredgecolor='b',label='$\sigma_\eta=0.7$')

ax.tick_params(axis='both', which='major', labelsize=14)
ax.set_xticks(np.arange(0,250,50), minor=False)
ax.set_xticklabels([])
ax.set_xticks(np.arange(0,250,50), minor=True)
#ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%d'))
#ax.xaxis.set_major_formatter(sformatter)
plt.xlim(0,210)

ax.set_yticks(np.arange(-1,2.5,0.1), minor=False)
ax.set_yticklabels(np.arange(-1,2.5,0.1), minor=False, family='serif')
ax.set_yticks(np.arange(-1,2.5,0.01), minor=True)
plt.ylim(0.0,0.25)
ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%1f'))
ax.yaxis.set_major_formatter(sformatter)

ax.legend(loc=(0.6,0.38),fontsize=18);

#plt.xlabel('$\sqrt{s}_{NN}$ (GeV)',fontsize=18 , weight='normal')
plt.ylabel('$S\sigma=C_3/C_1$', fontsize=22, weight='normal')
text(200,0.23,'net charge',fontsize=22,ha='right')

#########################################
plt.savefig('bw_skewness_sigmaeta.pdf',format='pdf')
os.system('open -a Preview bw_skewness_sigmaeta.pdf')

#plt.show()
quit()
