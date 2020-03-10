
cd ../../..

./bw figs/paper/pars/vsT/T110.dat > figs/paper/crap/crap1.dat &
./bw figs/paper/pars/vsT/T125.dat > figs/paper/crap/crap2.dat &
./bw figs/paper/pars/vsT/T140.dat > figs/paper/crap/crap3.dat &
./bw figs/paper/pars/vsT/T155.dat > figs/paper/crap/crap4.dat &

./bw figs/paper/pars/vsrhoB/rhoB0.1.dat > figs/paper/crap/crap5.dat &
./bw figs/paper/pars/vsrhoB/rhoB0.2.dat > figs/paper/crap/crap6.dat &
./bw figs/paper/pars/vsrhoB/rhoB0.15.dat > figs/paper/crap/crap8.dat &
./bw figs/paper/pars/vsrhoB/rhoB0.dat > figs/paper/crap/crap9.dat &

./bw figs/paper/pars/vsOmega/Omega50.dat > figs/paper/crap/crap10.dat &
./bw figs/paper/pars/vsOmega/Omega250.dat > figs/paper/crap/crap12.dat &

./bw figs/paper/pars/vsEta/eta0.2.dat > figs/paper/crap/crap13.dat &
./bw figs/paper/pars/vsEta/eta1.0.dat > figs/paper/crap/crap13.dat &

./bw figs/paper/pars/vsAcc/Acc0.2.dat > figs/paper/crap/crap14.dat &
./bw figs/paper/pars/vsAcc/Acc0.5.dat > figs/paper/crap/crap15.dat &
./bw figs/paper/pars/vsAcc/Acc0.6.dat > figs/paper/crap/crap16.dat &
./bw figs/paper/pars/vsAcc/Acc0.8.dat > figs/paper/crap/crap17.dat &
./bw figs/paper/pars/vsAcc/Acc0.dat > figs/paper/crap/crap18.dat &
./bw figs/paper/pars/vsAcc/Acc1.dat > figs/paper/crap/crap19.dat &

cd figs/paper/data

for i in 0 1 2 3 4 5 6
  cp roots${i}_T140.dat roots${i}_Omega100.dat
  cp roots${i}_T140.dat roots${i}_rhoB0.05.dat
  cp roots${i}_T140.dat roots${i}_eta0.6.dat
  cp roots${i}_T140.dat roots${i}_Acc0.4.dat
