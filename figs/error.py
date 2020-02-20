import numpy as np

nroots=7
nruns=10

for i in range(nroots):
    sumS=0
    sumK=0
    Omega=[]
    pbar=[]
    sigma2p=[]
    Ssigmap=[]
    Ksigma2p=[]
    Skellamp=[]
    Ssigma_avg=0
    Ksigma2_avg=0
    file="../errbars/moments_roots"+str(i)+".dat";
    mydata = np.loadtxt(file,skiprows=1,unpack=True)
    for run in range(nruns):
        Omega.append(mydata[0][run])

        pbar.append(mydata[5][run])
        sigma2p.append(mydata[6][run])
        Ssigmap.append(mydata[7][run])
        Ksigma2p.append(mydata[8][run])

        Skellamp.append(sigma2p[run]/(pbar[run]*Omega[run]))

        Ssigma_avg+=Ssigmap[run]*Skellamp[run]
        Ksigma2_avg+=Ksigma2p[run]

    Ssigma_avg*=1/nroots
    Ksigma2_avg*=1/nroots

    for i in range(nroots):
        sumS+=(Ssigmap[i]*Skellamp[i]-Ssigma_avg)**2
        sumK+=(Ksigma2p[i]-Ksigma2_avg)**2

    stderrS=(1/nroots)*np.sqrt(sumS)
    stderrK=(1/nroots)*np.sqrt(sumK)

    print("Ksigma2 error =",stderrK,", Ssigma error =",stderrS)
