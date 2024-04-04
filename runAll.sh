#!/bin/bash
#72
for i in {0..7}
	do 
		sbatch slurm.batch $i
		sleep 5 
	done

