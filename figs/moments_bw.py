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

#	Omega,roots,T,rhoB,rhoQ,pbar/Omega,sigma2p,Ssigmap,Ksigma2p,qbar/Omega,sigma2q,Ssigmaq,Ksigma2q,
#	kbar/Omega,sigma2k,Ssigmak,Ksigma2k,pibar/Omega,sigma2pi,Ssigmapi,Ksigma2pi,totp/Omega,totq/Omega,totk/Omega,totpi/Omega

fluc='charge'
#fluc='proton'

mydata = np.loadtxt('moments_bw_Omega50.dat',skiprows=1,unpack=True)
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

if fluc == 'charge':
  plt.plot(roots,Ssigmaq*Skellamq,linestyle='-',linewidth=2,color='r',markersize=8, marker='s', markerfacecolor=None, markeredgecolor=None,label='BW: $C_3/C_1$')
  plt.plot(roots,Ksigma2q,linestyle='--',linewidth=2,color='b',markersize=10, marker='^', markerfacecolor=None, markeredgecolor=None,label='BW: $C_4/C_2$')

if fluc == 'proton':
  plt.plot(roots,Ssigmap*Skellamp,linestyle='-',linewidth=2,color='r',markersize=8, marker='s', markerfacecolor=None, markeredgecolor=None,label='BW: $C_3/C_1$')
  plt.plot(roots,Ksigma2p,linestyle='--',linewidth=2,color='b',markersize=10, marker='^', markerfacecolor=None, markeredgecolor=None,label='BW: $C_4/C_2$')

if fluc == 'proton':
  stardata = np.loadtxt('starmoments_netp.txt',skiprows=1,unpack=True)
  roots_star=stardata[0]
  C1=stardata[4]
  E1=sqrt(stardata[5]*stardata[5]+stardata[6]*stardata[6])
  C2=stardata[7]
  E2=sqrt(stardata[8]*stardata[8]+stardata[9]*stardata[9])
  C3=stardata[10]
  E3=sqrt(stardata[11]*stardata[11]+stardata[12]*stardata[12])
  C4=stardata[13]
  E4=sqrt(stardata[14]*stardata[14]+stardata[15]*stardata[15])
  C4overC2=C4/C2
  E42=C4overC2*(sqrt(E4*E4/(C4*C4) +E2*E2/(C2*C2)))
  C3overC1=C3/C1
  E31=C3overC1*(sqrt(E3*E3/(C3*C3) +E1*E1/(C1*C1)))

if fluc == 'charge':
  stardata = np.loadtxt('starmoments_netq.txt',skiprows=1,unpack=True)
  roots_star=stardata[0]
  C1=stardata[1]
  E1=sqrt(stardata[2]*stardata[2]+stardata[3]*stardata[3])
  C2=stardata[4]
  E2=sqrt(stardata[5]*stardata[5]+stardata[6]*stardata[6])
  C3=stardata[7]*pow(C2,1.5)
  E3=sqrt(stardata[8]*stardata[8]+stardata[8]*stardata[8])
  C4=stardata[9]*C2*C2
  E4=sqrt(stardata[10]*stardata[10]+stardata[11]*stardata[11])
  C4overC2=C4/C2
  E42=C4overC2*(sqrt(E4*E4/(C4*C4) +E2*E2/(C2*C2)))
  C3overC1=C3/C1
  E31=C3overC1*(sqrt(E3*E3/(C3*C3) +E1*E1/(C1*C1)))


plt.errorbar(roots_star,C3overC1,E31,linestyle='--',linewidth=2,color='r',markersize=14, marker='*', markerfacecolor=None, markeredgecolor=None,label='STAR: $C_3/C_1$')
plt.errorbar(roots_star,C4overC2,E42,linestyle='-',linewidth=2,color='b',markersize=14, marker='*', markerfacecolor=None, markeredgecolor=None,label='STAR: $C_4/C_2$')

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

ax.legend(loc=(0.52,0.48));

plt.xlabel('$\sqrt{s}_{NN}$ (GeV)',fontsize=18 , weight='normal')
if fluc == 'proton':
  plt.ylabel('Net Proton Fluctuations', fontsize=24, weight='normal')
  plt.savefig('moments_bw_proton.pdf',format='pdf')
  os.system('xdg-open moments_bw_proton.pdf')
if fluc == 'charge':
  plt.ylabel('Net Charge Fluctuations', fontsize=24, weight='normal')
  plt.savefig('moments_bw_charge.pdf',format='pdf')
  os.system('xdg-open moments_bw_charge.pdf')



#plt.show()
quit()
