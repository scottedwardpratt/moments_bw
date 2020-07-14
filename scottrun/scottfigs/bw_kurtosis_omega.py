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

Omega=25
Mydata25=np.loadtxt('../results/sigmaeta0.3_Omega25.dat',skiprows=1,unpack=True)
roots25=Mydata25[1]
pbar25=Mydata25[5]
sigma2p25=Mydata25[6]
Ssigmap25=Mydata25[7]
Ksigma2p25=Mydata25[8]
qbar25=Mydata25[9]
sigma2q25=Mydata25[10]
Ssigmaq25=Mydata25[11]
Ksigma2q25=Mydata25[12]
Skellamp25=sigma2p25/(pbar25*Omega)
Skellamq25=sigma2q25/(qbar25*Omega)

Omega=50
Mydata50=np.loadtxt('../results/sigmaeta0.3_Omega50.dat',skiprows=1,unpack=True)
roots50=Mydata50[1]
pbar50=Mydata50[5]
sigma2p50=Mydata50[6]
Ssigmap50=Mydata50[7]
Ksigma2p50=Mydata50[8]
qbar50=Mydata50[9]
sigma2q50=Mydata50[10]
Ssigmaq50=Mydata50[11]
Ksigma2q50=Mydata50[12]
Skellamp50=sigma2p50/(pbar50*Omega)
Skellamq50=sigma2q50/(qbar50*Omega)

Omega=100
Mydata100=np.loadtxt('../results/sigmaeta0.3_Omega100.dat',skiprows=1,unpack=True)
roots100=Mydata100[1]
pbar100=Mydata100[5]
sigma2p100=Mydata100[6]
Ssigmap100=Mydata100[7]
Ksigma2p100=Mydata100[8]
qbar100=Mydata100[9]
sigma2q100=Mydata100[10]
Ssigmaq100=Mydata100[11]
Ksigma2q100=Mydata100[12]
Skellamp100=sigma2p100/(pbar100*Omega)
Skellamq100=sigma2q100/(qbar100*Omega)

Omega=200
Mydata200=np.loadtxt('../results/sigmaeta0.3_Omega200.dat',skiprows=1,unpack=True)
roots200=Mydata200[1]
pbar200=Mydata200[5]
sigma2p200=Mydata200[6]
Ssigmap200=Mydata200[7]
Ksigma2p200=Mydata200[8]
qbar200=Mydata200[9]
sigma2q200=Mydata200[10]
Ssigmaq200=Mydata200[11]
Ksigma2q200=Mydata200[12]
Skellamp200=sigma2p200/(pbar200*Omega)
Skellamq200=sigma2q200/(qbar200*Omega)

Omega=400
Mydata400=np.loadtxt('../results/sigmaeta0.3_Omega400.dat',skiprows=1,unpack=True)
roots400=Mydata400[1]
pbar400=Mydata400[5]
sigma2p400=Mydata400[6]
Ssigmap400=Mydata400[7]
Ksigma2p400=Mydata400[8]
qbar400=Mydata400[9]
sigma2q400=Mydata400[10]
Ssigmaq400=Mydata400[11]
Ksigma2q400=Mydata400[12]
Skellamp400=sigma2p400/(pbar400*Omega)
Skellamq400=sigma2q400/(qbar400*Omega)

#################################################################
######## LOWER PANEL protons
ax = fig.add_axes([0.15,0.12,0.84,0.43])

#plt.plot(roots,Ssigmap*Skellamp,linestyle='-',linewidth=2,color='r',markersize=8, marker='s', markerfacecolor=None, markeredgecolor=None,label='ETA=0.3: $C_3/C_1$')
plt.plot(roots25,Ksigma2p25,linestyle='--',linewidth=2,color='r',markersize=10, marker='s', markerfacecolor='r', markeredgecolor='r',label='$\Omega=25$')
plt.plot(roots50,Ksigma2p50,linestyle='--',linewidth=2,color='k',markersize=10, marker='s', markerfacecolor='k', markeredgecolor='k',label='$\Omega=50$')
plt.plot(roots100,Ksigma2p100,linestyle='--',linewidth=2,color='g',markersize=10, marker='s', markerfacecolor='g', markeredgecolor='g',label='$\Omega=100$')
plt.plot(roots200,Ksigma2p200,linestyle='--',linewidth=2,color='b',markersize=10, marker='s', markerfacecolor='b', markeredgecolor='b',label='$\Omega=200$')
plt.plot(roots400,Ksigma2p400,linestyle='--',linewidth=2,color='cyan',markersize=10, marker='s', markerfacecolor='cyan', markeredgecolor='cyan',label='$\Omega=400$')

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
plt.ylim(0.25,1.05)
ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%1f'))
ax.yaxis.set_major_formatter(sformatter)

ax.legend(loc=(0.6,0.01),fontsize=18);

plt.xlabel('$\sqrt{s}_{NN}$ (GeV)',fontsize=18 , weight='normal')
plt.ylabel('$K\sigma^2=C_4/C_2$', fontsize=24, weight='normal')
text(120,0.45,'net protons',fontsize=22,ha='right')

######## Upper Panel charge
ax = fig.add_axes([0.15,0.55,0.84,0.43])

#plt.plot(roots,Ssigmap*Skellamp,linestyle='-',linewidth=2,color='r',markersize=8, marker='s', markerfacecolor=None, markeredgecolor=None,label='ETA=0.3: $C_3/C_1$')
plt.plot(roots25,Ksigma2q25,linestyle='--',linewidth=2,color='r',markersize=10, marker='s', markerfacecolor='r', markeredgecolor='r',label='$\Omega=25$')
plt.plot(roots50,Ksigma2q50,linestyle='--',linewidth=2,color='k',markersize=10, marker='s', markerfacecolor='k', markeredgecolor='k',label='$\Omega=50$')
plt.plot(roots100,Ksigma2q100,linestyle='--',linewidth=2,color='g',markersize=10, marker='s', markerfacecolor='g', markeredgecolor='g',label='$\Omega=100$')
plt.plot(roots200,Ksigma2q200,linestyle='--',linewidth=2,color='b',markersize=10, marker='s', markerfacecolor='b', markeredgecolor='b',label='$\Omega=200$')
plt.plot(roots400,Ksigma2q400,linestyle='--',linewidth=2,color='cyan',markersize=10, marker='s', markerfacecolor='cyan', markeredgecolor='cyan',label='$\Omega=400$')

ax.tick_params(axis='both', which='major', labelsize=14)
ax.set_xticks(np.arange(0,250,50), minor=False)
ax.set_xticklabels([])
ax.set_xticks(np.arange(0,250,50), minor=True)
#ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%d'))
#ax.xaxis.set_major_formatter(sformatter)
plt.xlim(0,210)

ax.set_yticks(np.arange(-1,2.5,0.2), minor=False)
ax.set_yticklabels(np.arange(-1,2.5,0.2), minor=False, family='serif')
ax.set_yticks(np.arange(-1,2.5,0.05), minor=True)
plt.ylim(0.25,1.05)
ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%1f'))
ax.yaxis.set_major_formatter(sformatter)

ax.legend(loc=(0.6,0.01),fontsize=18);

#plt.xlabel('$\sqrt{s}_{NN}$ (GeV)',fontsize=18 , weight='normal')
plt.ylabel('$K\sigma^2=C_4/C_2$', fontsize=24, weight='normal')
text(120,0.45,'net charge',fontsize=22,ha='right')

#########################################
plt.savefig('bw_kurtosis_omega.pdf',format='pdf')
os.system('open -a Preview bw_kurtosis_omega.pdf')



#plt.show()
quit()
