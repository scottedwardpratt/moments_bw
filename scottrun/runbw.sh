#!/bin/bash
cd /home/scott/git/moments_bw/scottrun
for ((iroots=4;iroots<7;iroots++));
do
	echo ${iroots} | bw default > logfiles/log${iroots}.txt &
done