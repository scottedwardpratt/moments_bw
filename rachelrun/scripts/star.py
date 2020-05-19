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

figname="../figs/stardata.pdf"
stardata = np.loadtxt('../data/starmoments_netq.txt',skiprows=1,unpack=True)

roots_star=stardata[0]

Ssigma=stardata[13]
Serror=sqrt(stardata[14]*stardata[14]+stardata[15]*stardata[15])
Ksigma2=stardata[16]
Kerror=sqrt(stardata[17]*stardata[17]+stardata[18]*stardata[18])

"""
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
"""

plt.errorbar(roots_star,Ssigma,Serror,linestyle='-',linewidth=2,color='r',markersize=11, marker='^', markerfacecolor=None, markeredgecolor=None)
plt.errorbar(roots_star,Ksigma2,Kerror,linestyle='--',linewidth=2,color='r',markersize=11, marker='s', markerfacecolor=None, markeredgecolor=None)

ax.tick_params(axis='both', which='major', labelsize=14)
ax.set_xticks(np.arange(0,250,50), minor=False)
ax.set_xticklabels(np.arange(0,250,25), minor=False, family='serif')
ax.set_xticks(np.arange(0,250,50), minor=True)
ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%d'))
ax.xaxis.set_major_formatter(sformatter)
plt.xlim(0,210)

ax.set_yticks(np.arange(-1,2.5,0.5), minor=False)
ax.set_yticklabels(np.arange(-1,2.5,0.5), minor=False, family='serif')
ax.set_yticks(np.arange(-1,2.5,0.05), minor=True)
plt.ylim(-0.5,2.55)
ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%1f'))
ax.yaxis.set_major_formatter(sformatter)

plt.xlabel('$\sqrt{s}$ (GeV)',fontsize=18 , weight='normal')
plt.ylabel('$S\sigma$,  $\kappa\sigma^2$', fontsize=24, weight='normal')
plt.savefig(figname,format='pdf',bbox_inches = "tight")
os.system('xdg-open '+figname)

quit()
