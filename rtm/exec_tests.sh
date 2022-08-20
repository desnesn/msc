#!/bin/bash

NUMBER_OF_THREADS=24

#for (( i = 1 ; i <= $NUMBER_OF_THREADS ; i++ ))
for (( i = NUMBER_OF_THREADS ; i >= 1 ; i-- ))
do
	echo "Running RTM migration with $i/$NUMBER_OF_THREADS threads..."
	rtime="$( TIMEFORMAT='%3R';time ( ./main $i ) 2>&1 1>/dev/null )"
	echo "Total Execution Time: $rtime"
	echo ""
	echo "$i"",""$rtime" >> log_big
done

