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
ax = fig.add_axes([0.17,0.12,0.8,0.8])

mydata = np.loadtxt('../momentsdata/Omega50_rhoB8_rhoQ4_rhoS0.dat',skiprows=0,unpack=True)
acceptance=100*mydata[0]
bbar_b0=mydata[1]
sigma2b_b0=mydata[2]
Ssigmab_b0=mydata[3]
Ksigma2b_b0=mydata[4]
qbar_b0=mydata[5]
sigma2q_b0=mydata[6]
Ssigmaq_b0=mydata[7]
Ksigma2q_b0=mydata[8]
plt.plot(acceptance,Ksigma2q_b0,linestyle='-',linewidth=2,color='r',markersize=3, marker='o', markerfacecolor=None, markeredgecolor=None)
plt.plot(acceptance,Ssigmaq_b0,linestyle='-',linewidth=2,color='b',markersize=3, marker='o', markerfacecolor=None, markeredgecolor=None)
mydata = np.loadtxt('../momentsdata/Omega100_rhoB8_rhoQ4_rhoS0.dat',skiprows=0,unpack=True)
acceptance=100*mydata[0]
bbar_b0=mydata[1]
sigma2b_b0=mydata[2]
Ssigmab_b0=mydata[3]
Ksigma2b_b0=mydata[4]
qbar_b0=mydata[5]
sigma2q_b0=mydata[6]
Ssigmaq_b0=mydata[7]
Ksigma2q_b0=mydata[8]
plt.plot(acceptance,Ksigma2q_b0,linestyle='-',linewidth=2,color='r',markersize=3, marker='s', markerfacecolor=None, markeredgecolor=None)
plt.plot(acceptance,Ssigmaq_b0,linestyle='-',linewidth=2,color='b',markersize=3, marker='s', markerfacecolor=None, markeredgecolor=None)
mydata = np.loadtxt('../momentsdata/Omega200_rhoB8_rhoQ4_rhoS0.dat',skiprows=0,unpack=True)
acceptance=100*mydata[0]
bbar_b0=mydata[1]
sigma2b_b0=mydata[2]
Ssigmab_b0=mydata[3]
Ksigma2b_b0=mydata[4]
qbar_b0=mydata[5]
sigma2q_b0=mydata[6]
Ssigmaq_b0=mydata[7]
Ksigma2q_b0=mydata[8]
plt.plot(acceptance,Ksigma2q_b0,linestyle='-',linewidth=2,color='r',markersize=3, marker='^', markerfacecolor=None, markeredgecolor=None)
plt.plot(acceptance,Ssigmaq_b0,linestyle='-',linewidth=2,color='b',markersize=3, marker='', markerfacecolor=None, markeredgecolor=None)
mydata = np.loadtxt('../momentsdata/Omega400_rhoB8_rhoQ4_rhoS0.dat',skiprows=0,unpack=True)
acceptance=100*mydata[0]
bbar_b0=mydata[1]
sigma2b_b0=mydata[2]
Ssigmab_b0=mydata[3]
Ksigma2b_b0=mydata[4]
qbar_b0=mydata[5]
sigma2q_b0=mydata[6]
Ssigmaq_b0=mydata[7]
Ksigma2q_b0=mydata[8]
plt.plot(acceptance,Ksigma2q_b0,linestyle='-',linewidth=2,color='r',markersize=3, marker='*', markerfacecolor=None, markeredgecolor=None)
plt.plot(acceptance,Ssigmaq_b0,linestyle='-',linewidth=2,color='b',markersize=3, marker='*', markerfacecolor=None, markeredgecolor=None)


ax.tick_params(axis='both', which='major', labelsize=14)
ax.set_xticks(np.arange(0,110,20), minor=False)
ax.set_xticklabels(np.arange(0,110,20), minor=False, family='serif')
ax.set_xticks(np.arange(0,110,10), minor=True)
ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%d'))
ax.xaxis.set_major_formatter(sformatter)
plt.xlim(0,100)

ax.set_yticks(np.arange(-1,1.5,0.5), minor=False)
ax.set_yticklabels(np.arange(-1,1.5,0.5), minor=False, family='serif')
ax.set_yticks(np.arange(-1,1.5,0.05), minor=True)
plt.ylim(-0.2,1.1)
ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%1f'))
ax.yaxis.set_major_formatter(sformatter)

plt.xlabel('acceptance  (percent)',fontsize=18 , weight='normal')
plt.ylabel('$K\sigma^2$ (MeV)', fontsize=18, weight='normal')

plt.savefig('moments_vs_acceptance.pdf',format='pdf')
os.system('open -a Preview moments_vs_Omega.pdf')
#plt.show()
quit()
