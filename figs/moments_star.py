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

mydata = np.loadtxt('cumulants_vs_E.txt',skiprows=1,unpack=True)
roots=mydata[0]
C1=mydata[4]
E1=sqrt(mydata[5]*mydata[5]+mydata[6]*mydata[6])
C2=mydata[7]
E2=sqrt(mydata[8]*mydata[8]+mydata[9]*mydata[9])
C3=mydata[10]
E3=sqrt(mydata[11]*mydata[11]+mydata[12]*mydata[12])
C4=mydata[13]
E4=sqrt(mydata[14]*mydata[14]+mydata[15]*mydata[15])
C4overC2=C4/C2
E42=C4overC2*(sqrt(E4*E4/(C4*C4) +E2*E2/(C2*C2)))
C3overC1=C3/C1
E31=C3overC1*(sqrt(E3*E3/(C3*C3) +E1*E1/(C1*C1)))

#plt.errorbar(roots,C4overC2,E42,linestyle='-',linewidth=2,color='r',markersize=5, marker='o', markerfacecolor=None, markeredgecolor=None)
#plt.errorbar(roots,C3overC1,E31,linestyle='-',linewidth=2,color='b',markersize=5, marker='o', markerfacecolor=None, markeredgecolor=None)

a=0.166
b=0.139
c=0.053
d=1.308
e=0.273
rhoN=0.08436
muB=d/(1.0+e*roots)
T=a-b*muB*muB-c*muB*muB*muB*muB
rhoB=rhoN*sinh(muB/T)
for i in range (0,7):
  print(roots[i],1000*T[i],1000*muB[i],rhoB[i])

plt.errorbar(rhoB,C4overC2,E42,linestyle='-',linewidth=2,color='r',markersize=5, marker='o', markerfacecolor=None, markeredgecolor=None)
plt.errorbar(rhoB,C3overC1,E31,linestyle='-',linewidth=2,color='b',markersize=5, marker='o', markerfacecolor=None, markeredgecolor=None)

#plt.xscale('log')
#ax.tick_params(axis='both', which='major', labelsize=14)
#ax.set_xticks([1.0,2.0,5.0,10.0,20.0,50.0,100.0,200.0], minor=False)
#ax.set_xticks([1.0,2.0,5.0,10.0,20.0,50.0,100.0,200.0], minor=True)
#ax.set_xticklabels([1.0,2.0,5.0,10.0,20.0,50.0,100.0,200.0], minor=False, family='serif')
#ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%d'))

ax.tick_params(axis='both', which='major', labelsize=14)
ax.set_xticks(np.arange(0.0,0.2,0.05), minor=False)
ax.set_xticks(np.arange(0.0,0.2,0.01), minor=True)
ax.set_xticklabels(np.arange(0.0,0.2,0.05), minor=False,family='serif')
ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%4.2f'))

ax.xaxis.set_major_formatter(sformatter)
#plt.xlim(7,210)
plt.xlim(0,1)


ax.set_yticks(np.arange(-1,1.5,0.5), minor=False)
ax.set_yticklabels(np.arange(-1,1.5,0.5), minor=False, family='serif')
ax.set_yticks(np.arange(-1,1.5,0.05), minor=True)
plt.ylim(-0.05,1.05)
ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%1f'))
ax.yaxis.set_major_formatter(sformatter)

plt.xlabel('$\sqrt{s}_{NN}$ (GeV)',fontsize=18 , weight='normal')
plt.ylabel('$C_4/C_2$,  $C_3/C_1$', fontsize=24, weight='normal')

plt.savefig('moments_star.pdf',format='pdf')
os.system('xdg-open moments_star.pdf')
#plt.show()
quit()
