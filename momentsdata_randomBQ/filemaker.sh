#!/bin/sh
a=0.25
Omega=200.0
/bin/rm -f moments_vs_rhoB.dat
echo \# acceptance=${a}, Omega=${Omega} > moments_vs_rhoB.dat

rhoB=0
rhoQ=0
filename=Acceptance${a}_rhoB${rhoB}_rhoQ${rhoQ}.dat
echo filename=${filename}
grep ${Omega} ${filename} >> moments_vs_rhoB.dat

rhoB=0.01
rhoQ=0.005
filename=Acceptance${a}_rhoB${rhoB}_rhoQ${rhoQ}.dat
echo filename=${filename}
grep ${Omega} ${filename} >> moments_vs_rhoB.dat

rhoB=0.02
rhoQ=0.01
filename=Acceptance${a}_rhoB${rhoB}_rhoQ${rhoQ}.dat
echo filename=${filename}
grep ${Omega} ${filename} >> moments_vs_rhoB.dat

rhoB=0.03
rhoQ=0.015
filename=Acceptance${a}_rhoB${rhoB}_rhoQ${rhoQ}.dat
echo filename=${filename}
grep ${Omega} ${filename} >> moments_vs_rhoB.dat

rhoB=0.04
rhoQ=0.02
filename=Acceptance${a}_rhoB${rhoB}_rhoQ${rhoQ}.dat
echo filename=${filename}
grep ${Omega} ${filename} >> moments_vs_rhoB.dat

rhoB=0.08
rhoQ=0.04
filename=Acceptance${a}_rhoB${rhoB}_rhoQ${rhoQ}.dat
echo filename=${filename}
grep ${Omega} ${filename} >> moments_vs_rhoB.dat

rhoB=0.12
rhoQ=0.06
filename=Acceptance${a}_rhoB${rhoB}_rhoQ${rhoQ}.dat
echo filename=${filename}
grep ${Omega} ${filename} >> moments_vs_rhoB.dat