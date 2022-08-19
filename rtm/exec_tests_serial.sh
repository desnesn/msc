#!/bin/bash

freq=12.0

max_number_of_threads=24

number_of_shots=1

vel[1]="velocity_smooth_05_05_downsampled[151,461,0.0,0.0,20.0,20.0].bin"
vel[2]="velocity_smooth_05_05_downsampled[188,576,0.0,0.0,16.0,16.0].bin"
vel[3]="velocity_smooth_05_05_downsampled[251,767,0.0,0.0,12.0,12.0].bin"
vel[4]="velocity_smooth_05_05_downsampled[376,1151,0.0,0.0,8.0,8.0].bin"
vel[5]="velocity_smooth_05_05_downsampled[751,2301,0.0,0.0,4.0,4.0].bin"

##############################################################################

log=final_execution_serial

if [ -e "$log" ]; then
	rm $log
fi

echo "RTM Serial Execution"

for (( exec_time = 1 ; $exec_time <= 10 ; exec_time++ ))
do
	for (( v = 1 ; v <= 5 ; v++ ))
	do
		vel_file=${vel[$v]}
		
		echo "|-------------------------------------------------------------------------|"
		echo "[$exec_time] serial execution with $vel_file"
		rtime=" $( ./main $vel_file $freq $number_of_shots 1 2>&1 1>/dev/null )"
		echo "Total Execution Time: $rtime"
		echo "$exec_time","$vel_file","$rtime" >> $log
		
		rm -f rtm_migrated_b1.0.su
	done
done

##############################################################################

