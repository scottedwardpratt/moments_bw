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
c2c1k25=sigma2k25/(kbar25*Omega)
c2c1p25=sigma2p25/(pbar25*Omega)
c2c1q25=sigma2q25/(qbar25*Omega)
c3c1p25=Ssigmap25*sigma2p25/(Omega*pbar25)
c3c1k25=Ssigmak25*sigma2k25/(Omega*kbar25)
c3c1q25=Ssigmaq25*sigma2q25/(Omega*qbar25)


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
c2c1k50=sigma2k50/(kbar50*Omega)
c2c1p50=sigma2p50/(pbar50*Omega)
c2c1q50=sigma2q50/(qbar50*Omega)
c3c1p50=Ssigmap50*sigma2p50/(Omega*pbar50)
c3c1k50=Ssigmak50*sigma2k50/(Omega*kbar50)
c3c1q50=Ssigmaq50*sigma2q50/(Omega*qbar50)

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
c2c1k100=sigma2k100/(kbar100*Omega)
c2c1p100=sigma2p100/(pbar100*Omega)
c2c1q100=sigma2q100/(qbar100*Omega)
c3c1p100=Ssigmap100*sigma2p100/(Omega*pbar100)
c3c1k100=Ssigmak100*sigma2k100/(Omega*kbar100)
c3c1q100=Ssigmaq100*sigma2q100/(Omega*qbar100)

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
c2c1k200=sigma2k200/(kbar200*Omega)
c2c1p200=sigma2p200/(pbar200*Omega)
c2c1q200=sigma2q200/(qbar200*Omega)
c3c1p200=Ssigmap200*sigma2p200/(Omega*pbar200)
c3c1k200=Ssigmak200*sigma2k200/(Omega*kbar200)
c3c1q200=Ssigmaq200*sigma2q200/(Omega*qbar200)

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
c2c1k400=sigma2k400/(kbar400*Omega)
c2c1p400=sigma2p400/(pbar400*Omega)
c2c1q400=sigma2q400/(qbar400*Omega)
c3c1p400=Ssigmap400*sigma2p400/(Omega*pbar400)
c3c1k400=Ssigmak400*sigma2k400/(Omega*kbar400)
c3c1q400=Ssigmaq400*sigma2q400/(Omega*qbar400)


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

#plt.plot(roots,Ssigmap*c2c1p,linestyle='-',linewidth=2,color='r',markersize=8, marker='s', markerfacecolor=None, markeredgecolor=None,label='ETA=0.3: $C_3/C_1$')
plt.plot(roots25,c3c1p25,linestyle='--',linewidth=2,color='r',markersize=10, marker='s', markerfacecolor='r', markeredgecolor='r',label='$\Omega=25$')
plt.plot(roots50,c3c1p50,linestyle='-',linewidth=2,color='k',markersize=10, marker='s', markerfacecolor='k', markeredgecolor='k',label='$\Omega=50$')
plt.plot(roots100,c3c1p100,linestyle=(0, (3, 10, 1, 10)),linewidth=2,color='g',markersize=10, marker='s', markerfacecolor='g', markeredgecolor='g',label='$\Omega=100$')
plt.plot(roots200,c3c1p200,linestyle=':',linewidth=2,color='b',markersize=10, marker='s', markerfacecolor='b', markeredgecolor='b',label='$\Omega=200$')
plt.plot(roots400,c3c1p400,linestyle=(0, (3, 1, 1, 1, 1, 1)),linewidth=2,color='cyan',markersize=10, marker='s', markerfacecolor='cyan', markeredgecolor='cyan',label='$\Omega=400$')
plt.errorbar(roots_star,c3c1_protons_star,c3c1_protons_error_star,linestyle=' ',linewidth=2,color='purple',markersize=11, marker='*', markerfacecolor=None, markeredgecolor=None)

ax.tick_params(axis='both', which='major', labelsize=16)
ax.set_xticks(np.arange(0,250,50), minor=False)
ax.set_xticklabels(np.arange(0,250,50), minor=False, family='serif',fontsize="14")
ax.set_xticks(np.arange(0,250,25), minor=True)
ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%d'))
ax.xaxis.set_major_formatter(sformatter)
plt.xlim(0,210)

ax.set_yticks(np.arange(-1,2.5,0.2), minor=False)
ax.set_yticklabels(np.arange(-1,2.5,0.2), minor=False, family='serif',fontsize="14")
ax.set_yticks(np.arange(-1,2.5,0.05), minor=True)
plt.ylim(0.4,1.1)
ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%1f'))
ax.yaxis.set_major_formatter(sformatter)

#ax.legend(loc=(0.6,0.4),fontsize=18);

plt.xlabel('$\sqrt{s}_{NN}$ (GeV)',fontsize=24 , weight='normal')
text(100,0.43,'(f) net protons',fontsize=22,ha='left')

######## Middle Panel kaons
ax = fig.add_axes([0.17,0.37,0.82,0.31])

stardata = np.loadtxt('data/star_netkaons_c3c1.txt',skiprows=1,unpack=True)
roots_star=stardata[0]
c12_kaons_star=stardata[1]
c32_kaons_star=stardata[3]
c3c1_kaons_star=c32_kaons_star/c12_kaons_star
error12=stardata[2]/c12_kaons_star
error32=stardata[4]/c32_kaons_star
Serror=sqrt(error12*error12+error32*error32)

print(c32_kaons_star)

#plt.plot(roots,Ssigmap*c2c1p,linestyle='-',linewidth=2,color='r',markersize=8, marker='s', markerfacecolor=None, markeredgecolor=None,label='ETA=0.3: $C_3/C_1$')
plt.plot(roots25,c3c1k25,linestyle='--',linewidth=2,color='r',markersize=10, marker='s', markerfacecolor='r', markeredgecolor='r',label='$\Omega=25$')
plt.plot(roots50,c3c1k50,linestyle='-',linewidth=2,color='k',markersize=10, marker='s', markerfacecolor='k', markeredgecolor='k',label='$\Omega=50$')
plt.plot(roots100,c3c1k100,linestyle=(0, (3, 10, 1, 10)),linewidth=2,color='g',markersize=10, marker='s', markerfacecolor='g', markeredgecolor='g',label='$\Omega=100$')
plt.plot(roots200,c3c1k200,linestyle=':',linewidth=2,color='b',markersize=10, marker='s', markerfacecolor='b', markeredgecolor='b',label='$\Omega=200$')
plt.plot(roots400,c3c1k400,linestyle=(0, (3, 1, 1, 1, 1, 1)),linewidth=2,color='cyan',markersize=10, marker='s', markerfacecolor='cyan', markeredgecolor='cyan',label='$\Omega=400$')
plt.errorbar(roots_star,c3c1_kaons_star,Serror,linestyle=' ',linewidth=2,color='purple',markersize=11, marker='*', markerfacecolor=None, markeredgecolor=None)

ax.tick_params(axis='both', which='major', labelsize=14)
ax.set_xticks(np.arange(0,250,50), minor=False)
ax.set_xticklabels([])
ax.set_xticks(np.arange(0,250,50), minor=True)
#ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%d'))
#ax.xaxis.set_major_formatter(sformatter)
plt.xlim(0,210)

ax.set_yticks(np.arange(-1,2.5,0.5), minor=False)
ax.set_yticklabels(np.arange(-1,2.5,0.5), minor=False, family='serif',fontsize="14")
ax.set_yticks(np.arange(-1,2.5,0.1), minor=True)
plt.ylim(0,1.8)
ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%1f'))
ax.yaxis.set_major_formatter(sformatter)


#plt.xlabel('$\sqrt{s}_{NN}$ (GeV)',fontsize=18 , weight='normal')
plt.ylabel('$C_3/C_1$', fontsize=24, weight='normal')
text(100,0.1,'(e) net kaons',fontsize=22,ha='left')
######## Upper Panel charge
ax = fig.add_axes([0.17,0.68,0.82,0.31])

stardata = np.loadtxt('data/star_netcharge_c3c1.txt',skiprows=1,unpack=True)
roots_star=stardata[0]
c21_charge_star=stardata[1]
c32_charge_star=stardata[4]
c3c1_charge_star=c32_charge_star*c21_charge_star
staterror1=stardata[2]/c21_charge_star
syserror1=stardata[3]/c21_charge_star
staterror3=stardata[5]/c32_charge_star
syserror3=stardata[6]/c32_charge_star
Serror=sqrt(syserror3*syserror3+staterror3*staterror3+syserror1*syserror1+staterror1*staterror1)

#plt.plot(roots,Ssigmap*c2c1p,linestyle='-',linewidth=2,color='r',markersize=8, marker='s', markerfacecolor=None, markeredgecolor=None,label='ETA=0.3: $C_3/C_1$')
plt.plot(roots25,c3c1q25,linestyle='--',linewidth=2,color='r',markersize=10, marker='s', markerfacecolor='r', markeredgecolor='r',label='$\Omega=25$')
plt.plot(roots50,c3c1q50,linestyle='-',linewidth=2,color='k',markersize=10, marker='s', markerfacecolor='k', markeredgecolor='k',label='$\Omega=50$')
plt.plot(roots100,c3c1q100,linestyle=(0, (3, 10, 1, 10)),linewidth=2,color='g',markersize=10, marker='s', markerfacecolor='g', markeredgecolor='g',label='$\Omega=100$')
plt.plot(roots200,c3c1q200,linestyle=':',linewidth=2,color='b',markersize=10, marker='s', markerfacecolor='b', markeredgecolor='b',label='$\Omega=200$')
plt.plot(roots400,c3c1q400,linestyle=(0, (3, 1, 1, 1, 1, 1)),linewidth=2,color='cyan',markersize=10, marker='s', markerfacecolor='cyan', markeredgecolor='cyan',label='$\Omega=400$')
plt.errorbar(roots_star,c3c1_charge_star,Serror,linestyle=' ',linewidth=2,color='purple',markersize=11, marker='*', markerfacecolor=None, markeredgecolor=None,label='STAR\n(preliminary)')

ax.tick_params(axis='both', which='major', labelsize=14)
ax.set_xticks(np.arange(0,250,50), minor=False)
ax.set_xticklabels([])
ax.set_xticks(np.arange(0,250,50), minor=True)
#ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%d'))
#ax.xaxis.set_major_formatter(sformatter)
plt.xlim(0,210)

ax.set_yticks(np.arange(-1,6.5,1.0), minor=False)
ax.set_yticklabels(np.arange(-1,6.5,1.0), minor=False, family='serif',fontsize="14")
ax.set_yticks(np.arange(-1,6.5,0.5), minor=True)
plt.ylim(-0.5,4.5)
ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%1f'))
ax.yaxis.set_major_formatter(sformatter)

#plt.xlabel('$\sqrt{s}_{NN}$ (GeV)',fontsize=18 , weight='normal')
text(100,-0.2,'(d) net charge',fontsize=22,ha='left')

ax.legend(loc=(0.4,0.285),fontsize=18);

#########################################
plt.savefig('bw_c3c1_omega.pdf',format='pdf')
os.system('open -a Preview bw_c3c1_omega.pdf')
#os.system('open -a Preview bw_c3c1_omega.pdf')



#plt.show()
quit()
