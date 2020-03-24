cd ../altdata

for i in 0 1 2 3 4 5 6
do
  cp roots${i}_T140.dat roots${i}_Omega100.dat
  cp roots${i}_T140.dat roots${i}_rhoB0.05.dat
  cp roots${i}_T140.dat roots${i}_eta0.6.dat
  cp roots${i}_T140.dat roots${i}_Acc0.4.dat
done
