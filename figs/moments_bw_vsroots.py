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


mydata = np.loadtxt('../moments_eta0.3.dat',skiprows=1,unpack=True)
Omega=mydata[0]
roots=mydata[1]
Temp=mydata[2]
rhoB=mydata[3]
rhoQ=mydata[4]

pbar=mydata[5]
sigma2p=mydata[6]
Ssigmap=mydata[7]
Ksigma2p=mydata[8]

qbar=mydata[9]
sigma2q=mydata[10]
Ssigmaq=mydata[11]
Ksigma2q=mydata[12]

kbar=mydata[13]
sigma2k=mydata[14]
Ssigmak=mydata[15]
Ksigma2k=mydata[16]

pibar=mydata[17]
sigma2pi=mydata[18]
Ssigmapi=mydata[19]
Ksigma2pi=mydata[20]

topp=mydata[21]
totk=mydata[22]
toppi=mydata[23]


Skellamp=sigma2p/(pbar*Omega)
Skellamq=sigma2q/(qbar*Omega)

plt.plot(roots,Ssigmap*Skellamp,linestyle='-',linewidth=2,color='r',markersize=8, marker='s', markerfacecolor=None, markeredgecolor=None,label='ETA=0.3: $C_3/C_1$')
plt.plot(roots,Ksigma2p,linestyle='--',linewidth=2,color='r',markersize=10, marker='^', markerfacecolor=None, markeredgecolor=None,label='ETA=0.3: $C_4/C_2$')

mydata = np.loadtxt('../moments_eta0.4.dat',skiprows=1,unpack=True)
Omega=mydata[0]
roots=mydata[1]
Temp=mydata[2]
rhoB=mydata[3]
rhoQ=mydata[4]

pbar=mydata[5]
sigma2p=mydata[6]
Ssigmap=mydata[7]
Ksigma2p=mydata[8]

qbar=mydata[9]
sigma2q=mydata[10]
Ssigmaq=mydata[11]
Ksigma2q=mydata[12]

kbar=mydata[13]
sigma2k=mydata[14]
Ssigmak=mydata[15]
Ksigma2k=mydata[16]

pibar=mydata[17]
sigma2pi=mydata[18]
Ssigmapi=mydata[19]
Ksigma2pi=mydata[20]

topp=mydata[21]
totk=mydata[22]
toppi=mydata[23]


Skellamp=sigma2p/(pbar*Omega)
Skellamq=sigma2q/(qbar*Omega)

plt.plot(roots,Ssigmap*Skellamp,linestyle='-',linewidth=2,color='b',markersize=8, marker='s', markerfacecolor=None, markeredgecolor=None,label='ETA=0.4: $C_3/C_1$')
plt.plot(roots,Ksigma2p,linestyle='--',linewidth=2,color='b',markersize=10, marker='^', markerfacecolor=None, markeredgecolor=None,label='ETA=0.4: $C_4/C_2$')

mydata = np.loadtxt('../moments_eta0.5.dat',skiprows=1,unpack=True)
Omega=mydata[0]
roots=mydata[1]
Temp=mydata[2]
rhoB=mydata[3]
rhoQ=mydata[4]

pbar=mydata[5]
sigma2p=mydata[6]
Ssigmap=mydata[7]
Ksigma2p=mydata[8]

qbar=mydata[9]
sigma2q=mydata[10]
Ssigmaq=mydata[11]
Ksigma2q=mydata[12]

kbar=mydata[13]
sigma2k=mydata[14]
Ssigmak=mydata[15]
Ksigma2k=mydata[16]

pibar=mydata[17]
sigma2pi=mydata[18]
Ssigmapi=mydata[19]
Ksigma2pi=mydata[20]

topp=mydata[21]
totk=mydata[22]
toppi=mydata[23]


Skellamp=sigma2p/(pbar*Omega)
Skellamq=sigma2q/(qbar*Omega)

plt.plot(roots,Ssigmap*Skellamp,linestyle='-',linewidth=2,color='g',markersize=8, marker='s', markerfacecolor=None, markeredgecolor=None,label='ETA=0.5: $C_3/C_1$')
plt.plot(roots,Ksigma2p,linestyle='--',linewidth=2,color='g',markersize=10, marker='^', markerfacecolor=None, markeredgecolor=None,label='ETA=0.5: $C_4/C_2$')

mydata = np.loadtxt('../moments_eta1.0.dat',skiprows=1,unpack=True)
Omega=mydata[0]
roots=mydata[1]
Temp=mydata[2]
rhoB=mydata[3]
rhoQ=mydata[4]

pbar=mydata[5]
sigma2p=mydata[6]
Ssigmap=mydata[7]
Ksigma2p=mydata[8]

qbar=mydata[9]
sigma2q=mydata[10]
Ssigmaq=mydata[11]
Ksigma2q=mydata[12]

kbar=mydata[13]
sigma2k=mydata[14]
Ssigmak=mydata[15]
Ksigma2k=mydata[16]

pibar=mydata[17]
sigma2pi=mydata[18]
Ssigmapi=mydata[19]
Ksigma2pi=mydata[20]

topp=mydata[21]
totk=mydata[22]
toppi=mydata[23]


Skellamp=sigma2p/(pbar*Omega)
Skellamq=sigma2q/(qbar*Omega)

plt.plot(roots,Ssigmap*Skellamp,linestyle='-',linewidth=2,color='k',markersize=8, marker='s', markerfacecolor=None, markeredgecolor=None,label='ETA=1.0: $C_3/C_1$')
plt.plot(roots,Ksigma2p,linestyle='--',linewidth=2,color='k',markersize=10, marker='^', markerfacecolor=None, markeredgecolor=None,label='ETA=1.0: $C_4/C_2$')


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
plt.ylabel('$S\sigma$,  $\kappa\sigma^2$', fontsize=24, weight='normal')
plt.savefig('moments_bw_vsroots_with1.0.pdf',format='pdf')
os.system('xdg-open moments_bw_vsroots_with1.0.pdf')



#plt.show()
quit()
