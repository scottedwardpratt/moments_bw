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
plt.figure(figsize=(6,12))
fig = plt.figure(1)

Omega=25
Mydata25=np.loadtxt('data/sigmaeta00.3_Omega25.dat',skiprows=1,unpack=True)
roots25=Mydata25[1]
pbar25=Mydata25[5]
sigma2p25=Mydata25[6]
Ssigmap25=Mydata25[7]
Ksigma2p25=Mydata25[8]
qbar25=Mydata25[9]
sigma2q25=Mydata25[10]
Ssigmaq25=Mydata25[11]
Ksigma2q25=Mydata25[12]
kbar25=Mydata25[17]
sigma2k25=Mydata25[18]
Ssigmak25=Mydata25[19]
Ksigma2k25=Mydata25[20]
pibar25=Mydata25[21]
sigma2pi25=Mydata25[22]
Ssigmapi25=Mydata25[23]
Ksigma2pi25=Mydata25[24]
Skellampi25=sigma2pi25/(pibar25*Omega)
Skellamk25=sigma2k25/(kbar25*Omega)
Skellamp25=sigma2p25/(pbar25*Omega)
Skellamq25=sigma2q25/(qbar25*Omega)


Omega=50
Mydata50=np.loadtxt('data/sigmaeta00.3_Omega50.dat',skiprows=1,unpack=True)
roots50=Mydata50[1]
pbar50=Mydata50[5]
sigma2p50=Mydata50[6]
Ssigmap50=Mydata50[7]
Ksigma2p50=Mydata50[8]
qbar50=Mydata50[9]
sigma2q50=Mydata50[10]
Ssigmaq50=Mydata50[11]
Ksigma2q50=Mydata50[12]
kbar50=Mydata50[17]
sigma2k50=Mydata50[18]
Ssigmak50=Mydata50[19]
Ksigma2k50=Mydata50[20]
pibar50=Mydata50[21]
sigma2pi50=Mydata50[22]
Ssigmapi50=Mydata50[23]
Ksigma2pi50=Mydata50[24]
Skellampi50=sigma2pi50/(pibar50*Omega)
Skellamk50=sigma2k50/(kbar50*Omega)
Skellamp50=sigma2p50/(pbar50*Omega)
Skellamq50=sigma2q50/(qbar50*Omega)

Omega=100
Mydata100=np.loadtxt('data/sigmaeta00.3_Omega100.dat',skiprows=1,unpack=True)
roots100=Mydata100[1]
pbar100=Mydata100[5]
sigma2p100=Mydata100[6]
Ssigmap100=Mydata100[7]
Ksigma2p100=Mydata100[8]
qbar100=Mydata100[9]
sigma2q100=Mydata100[10]
Ssigmaq100=Mydata100[11]
Ksigma2q100=Mydata100[12]
kbar100=Mydata100[17]
sigma2k100=Mydata100[18]
Ssigmak100=Mydata100[19]
Ksigma2k100=Mydata100[20]
pibar100=Mydata100[21]
sigma2pi100=Mydata100[22]
Ssigmapi100=Mydata100[23]
Ksigma2pi100=Mydata100[24]
Skellampi100=sigma2pi100/(pibar100*Omega)
Skellamk100=sigma2k100/(kbar100*Omega)
Skellamp100=sigma2p100/(pbar100*Omega)
Skellamq100=sigma2q100/(qbar100*Omega)

Omega=200
Mydata200=np.loadtxt('data/sigmaeta00.3_Omega200.dat',skiprows=1,unpack=True)
roots200=Mydata200[1]
pbar200=Mydata200[5]
sigma2p200=Mydata200[6]
Ssigmap200=Mydata200[7]
Ksigma2p200=Mydata200[8]
qbar200=Mydata200[9]
sigma2q200=Mydata200[10]
Ssigmaq200=Mydata200[11]
Ksigma2q200=Mydata200[12]
kbar200=Mydata200[17]
sigma2k200=Mydata200[18]
Ssigmak200=Mydata200[19]
Ksigma2k200=Mydata200[20]
pibar200=Mydata200[21]
sigma2pi200=Mydata200[22]
Ssigmapi200=Mydata200[23]
Ksigma2pi200=Mydata200[24]
Skellampi200=sigma2pi200/(pibar200*Omega)
Skellamk200=sigma2k200/(kbar200*Omega)
Skellamp200=sigma2p200/(pbar200*Omega)
Skellamq200=sigma2q200/(qbar200*Omega)

Omega=400
Mydata400=np.loadtxt('data/sigmaeta00.3_Omega400.dat',skiprows=1,unpack=True)
roots400=Mydata400[1]
pbar400=Mydata400[5]
sigma2p400=Mydata400[6]
Ssigmap400=Mydata400[7]
Ksigma2p400=Mydata400[8]
qbar400=Mydata400[9]
sigma2q400=Mydata400[10]
Ssigmaq400=Mydata400[11]
Ksigma2q400=Mydata400[12]
kbar400=Mydata400[17]
sigma2k400=Mydata400[18]
Ssigmak400=Mydata400[19]
Ksigma2k400=Mydata400[20]
pibar400=Mydata400[21]
sigma2pi400=Mydata400[22]
Ssigmapi400=Mydata400[23]
Ksigma2pi400=Mydata400[24]
Skellampi400=sigma2pi400/(pibar400*Omega)
Skellamk400=sigma2k400/(kbar400*Omega)
Skellamp400=sigma2p400/(pbar400*Omega)
Skellamq400=sigma2q400/(qbar400*Omega)

#stardata = np.loadtxt('data/starmoments_netp.txt',skiprows=1,unpack=True)
#roots_star=stardata[0]
#Ksigma2=stardata[4]
#Kerror=sqrt(stardata[5]*stardata[5]+stardata[6]*stardata[6])

stardata = np.loadtxt('data/star_netprotons_c1c2c3c4.txt',unpack=False,skiprows=1)
roots_star=stardata[0]

c1_protons_star=stardata[1]
stat1=stardata[2]
sys1=stardata[3]
c1_error_star=sqrt(stat1*stat1+sys1*sys1)
c1_relerror_star=c1_error_star/c1_protons_star

c2_protons_star=stardata[4]
stat2=stardata[5]
sys2=stardata[6]
c2_error_star=sqrt(stat2*stat2+sys2*sys2)
c2_relerror_star=c2_error_star/c2_protons_star

c3_protons_star=stardata[7]
stat3=stardata[8]
sys3=stardata[9]
c3_error_star=sqrt(stat3*stat3+sys3*sys3)
c3_relerror_star=c3_error_star/c3_protons_star

c4_protons_star=stardata[10]
stat4=stardata[11]
sys4=stardata[12]
c4_error_star=sqrt(stat4*stat4+sys4*sys4)
c4_relerror_star=c4_error_star/c4_protons_star

c2c1_protons_star=c1_protons_star/c2_protons_star
c2c1_protons_error_star=c2c1_protons_star*sqrt(c1_relerror_star*c1_relerror_star+c2_relerror_star*c2_relerror_star)

c3c1_protons_star=c3_protons_star/c1_protons_star
c3c1_protons_error_star=c3c1_protons_star*sqrt(c1_relerror_star*c1_relerror_star+c3_relerror_star*c3_relerror_star)

c4c2_protons_star=c4_protons_star/c2_protons_star
c4c2_protons_error_star=c4c2_protons_star*sqrt(c2_relerror_star*c2_relerror_star+c4_relerror_star*c4_relerror_star)

c3c2_protons_star=c3_protons_star/c2_protons_star
c3c2_protons_error_star=c3c2_protons_star*sqrt(c2_relerror_star*c2_relerror_star+c3_relerror_star*c3_relerror_star)

#################################################################
######## LOWER PANEL protons
ax = fig.add_axes([0.17,0.06,0.82,0.31])

#plt.plot(roots,Ssigmap*Skellamp,linestyle='-',linewidth=2,color='r',markersize=8, marker='s', markerfacecolor=None, markeredgecolor=None,label='ETA=0.3: $C_3/C_1$')
plt.plot(roots25,Ksigma2p25,linestyle='--',linewidth=2,color='r',markersize=10, marker='s', markerfacecolor='r', markeredgecolor='r',label='$\Omega=25$')
plt.plot(roots50,Ksigma2p50,linestyle='-',linewidth=2,color='k',markersize=10, marker='s', markerfacecolor='k', markeredgecolor='k',label='$\Omega=50$')
plt.plot(roots100,Ksigma2p100,linestyle=(0, (3, 10, 1, 10)),linewidth=2,color='g',markersize=10, marker='s', markerfacecolor='g', markeredgecolor='g',label='$\Omega=100$')
plt.plot(roots200,Ksigma2p200,linestyle=':',linewidth=2,color='b',markersize=10, marker='s', markerfacecolor='b', markeredgecolor='b',label='$\Omega=200$')
plt.plot(roots400,Ksigma2p400,linestyle=(0, (3, 1, 1, 1, 1, 1)),linewidth=2,color='cyan',markersize=10, marker='s', markerfacecolor='cyan', markeredgecolor='cyan',label='$\Omega=400$')
plt.errorbar(roots_star,c4c2_protons_star,c4c2_protons_error_star,linestyle=' ',linewidth=2,color='purple',markersize=11, marker='*', markerfacecolor=None, markeredgecolor=None,label='STAR\n(preliminary)')

ax.tick_params(axis='both', which='major', labelsize=16)
ax.set_xticks(np.arange(0,250,50), minor=False)
ax.set_xticklabels(np.arange(0,250,50), minor=False, family='serif',fontsize="14")
ax.set_xticks(np.arange(0,250,25), minor=True)
ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%d'))
ax.xaxis.set_major_formatter(sformatter)
plt.xlim(0,210)

ax.set_yticks(np.arange(-1,4,0.5), minor=False)
ax.set_yticklabels(np.arange(-1,4,0.5), minor=False, family='serif',fontsize="14")
ax.set_yticks(np.arange(-1,4,0.25), minor=True)
plt.ylim(0,3.15)
ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%1f'))
ax.yaxis.set_major_formatter(sformatter)

ax.legend(loc=(0.4,0.3),fontsize=18);

plt.xlabel('$\sqrt{s}_{NN}$ (GeV)',fontsize=24 , weight='normal')
#plt.ylabel('$K\sigma^2=C_4/C_2$', fontsize=24, weight='normal')
text(195,0.1,'(i) net protons',fontsize=22,ha='right')

######## Middle Panel kaons
ax = fig.add_axes([0.17,0.37,0.82,0.31])

stardata = np.loadtxt('data/starmoments_netk.txt',skiprows=1,unpack=True)
roots_star=stardata[0]
Ksigma2=stardata[4]
Kerror=sqrt(stardata[5]*stardata[5]+stardata[6]*stardata[6])

#plt.plot(roots,Ssigmap*Skellamp,linestyle='-',linewidth=2,color='r',markersize=8, marker='s', markerfacecolor=None, markeredgecolor=None,label='ETA=0.3: $C_3/C_1$')
plt.plot(roots25,Ksigma2k25,linestyle='--',linewidth=2,color='r',markersize=10, marker='s', markerfacecolor='r', markeredgecolor='r',label='$\Omega=25$')
plt.plot(roots50,Ksigma2k50,linestyle='-',linewidth=2,color='k',markersize=10, marker='s', markerfacecolor='k', markeredgecolor='k',label='$\Omega=50$')
plt.plot(roots100,Ksigma2k100,linestyle=(0, (3, 10, 1, 10)),linewidth=2,color='g',markersize=10, marker='s', markerfacecolor='g', markeredgecolor='g',label='$\Omega=100$')
plt.plot(roots200,Ksigma2k200,linestyle=':',linewidth=2,color='b',markersize=10, marker='s', markerfacecolor='b', markeredgecolor='b',label='$\Omega=200$')
plt.plot(roots400,Ksigma2k400,linestyle=(0, (3, 1, 1, 1, 1, 1)),linewidth=2,color='cyan',markersize=10, marker='s', markerfacecolor='cyan', markeredgecolor='cyan',label='$\Omega=400$')
plt.errorbar(roots_star,Ksigma2,Kerror,linestyle=' ',linewidth=2,color='purple',markersize=11, marker='*', markerfacecolor=None, markeredgecolor=None)

ax.tick_params(axis='both', which='major', labelsize=14)
ax.set_xticks(np.arange(0,250,50), minor=False)
ax.set_xticklabels([])
ax.set_xticks(np.arange(0,250,50), minor=True)
#ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%d'))
#ax.xaxis.set_major_formatter(sformatter)
plt.xlim(0,210)

ax.set_yticks(np.arange(-2,2.5,0.5), minor=False)
ax.set_yticklabels(np.arange(-2,2.5,0.5), minor=False, family='serif',fontsize="14")
ax.set_yticks(np.arange(-2,2.5,0.25), minor=True)
plt.ylim(-1.5,1.7)
ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%1f'))
ax.yaxis.set_major_formatter(sformatter)

#ax.legend(loc=(0.6,0.01),fontsize=18);

#plt.xlabel('$\sqrt{s}_{NN}$ (GeV)',fontsize=18 , weight='normal')
plt.ylabel('$C_4/C_2$', fontsize=24, weight='normal')
text(195,-1.3,'(h) net kaons',fontsize=22,ha='right')

######## Upper Panel charge
ax = fig.add_axes([0.17,0.68,0.82,0.31])

stardata = np.loadtxt('data/starmoments_netq.txt',skiprows=1,unpack=True)
roots_star=stardata[0]
Ksigma2=stardata[4]
Kerror=sqrt(stardata[5]*stardata[5]+stardata[6]*stardata[6])

#plt.plot(roots,Ssigmap*Skellamp,linestyle='-',linewidth=2,color='r',markersize=8, marker='s', markerfacecolor=None, markeredgecolor=None,label='ETA=0.3: $C_3/C_1$')
plt.plot(roots25,Ksigma2q25,linestyle='--',linewidth=2,color='r',markersize=10, marker='s', markerfacecolor='r', markeredgecolor='r')
plt.plot(roots50,Ksigma2q50,linestyle='-',linewidth=2,color='k',markersize=10, marker='s', markerfacecolor='k', markeredgecolor='k')
plt.plot(roots100,Ksigma2q100,linestyle=(0, (3, 10, 1, 10)),linewidth=2,color='g',markersize=10, marker='s', markerfacecolor='g', markeredgecolor='g')
plt.plot(roots200,Ksigma2q200,linestyle=':',linewidth=2,color='b',markersize=10, marker='s', markerfacecolor='b', markeredgecolor='b')
plt.plot(roots400,Ksigma2q400,linestyle=(0, (3, 1, 1, 1, 1, 1)),linewidth=2,color='cyan',markersize=10, marker='s', markerfacecolor='cyan', markeredgecolor='cyan')
plt.errorbar(roots_star,Ksigma2,Kerror,linestyle=' ',linewidth=2,color='purple',markersize=11, marker='*', markerfacecolor=None, markeredgecolor=None,label='STAR\n(preliminary)')

ax.tick_params(axis='both', which='major', labelsize=14)
ax.set_xticks(np.arange(0,250,50), minor=False)
ax.set_xticklabels([])
ax.set_xticks(np.arange(0,250,50), minor=True)
#ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%d'))
#ax.xaxis.set_major_formatter(sformatter)
plt.xlim(0,210)

ax.set_yticks(np.arange(-14,14,2), minor=False)
ax.set_yticklabels(np.arange(-14,14,2), minor=False, family='serif', fontsize="14")
ax.set_yticks(np.arange(-14,14,1), minor=True)
plt.ylim(-5,5)
ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%1f'))
ax.yaxis.set_major_formatter(sformatter)

#plt.xlabel('$\sqrt{s}_{NN}$ (GeV)',fontsize=18 , weight='normal')
#plt.ylabel('$K\sigma^2=C_4/C_2$', fontsize=24, weight='normal')
text(195,-4.7,'(g) net charge',fontsize=22,ha='right')

#########################################
plt.savefig('bw_kurtosis_omega.pdf',format='pdf')
#os.system('xdg-open bw_kurtosis_omega.pdf')
os.system('open -a Preview bw_kurtosis_omega.pdf')



#plt.show()
quit()
