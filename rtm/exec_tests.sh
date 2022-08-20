#!/bin/bash

NUMBER_OF_THREADS=24

#modelo=p
modelo=g

log=log_small_3.0

if [ -e "$log" ]; then
	rm $log
fi

#USANDO O TIME
#for (( i = 1 ; i <= $NUMBER_OF_THREADS ; i++ ))
#do
#	echo "Running RTM migration with $i/$NUMBER_OF_THREADS threads..."
#	rtime="$( TIMEFORMAT='%3R';time ( ./main -$modelo $i ) 2>&1 1>/dev/null )"
#	echo "Total Execution Time: $rtime"
#	echo ""
#	echo "$i"",""$rtime" >> $log
#done

#USANDO O STDERR
for (( i = 1 ; i <= $NUMBER_OF_THREADS ; i++ ))
do
	echo "|---------------------------------------------------------------------------|"
	echo "Running RTM migration with $i/$NUMBER_OF_THREADS threads..."
	rtime=" $( ./main -$modelo $i 2>&1 1>/dev/null )"
	echo "Total Execution Time: $rtime"
	echo "$i"",""$rtime" >> $log
done
