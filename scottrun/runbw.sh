#!/bin/bash
cd ${HOME}/git/moments_bw/scottrun
for ((iroots=0;iroots<3;iroots++));
do
	echo ${iroots} | bw default > logfiles/log${iroots}.txt &
done