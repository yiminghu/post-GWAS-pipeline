#!/bin/bash
mv ../*part* .
for iter in {1..10}
do
	sqCreateScript -n 256 -c 1 -m 10G -w 5:00:00 -q scavenge -N gnova_ns.task_part_part${iter} --logdir=gnova_ns.task_part_part${iter}.log gnova_ns.task_part${iter} > gnova_ns.task_part${iter}_slurm.sh
	sbatch < gnova_ns.task_part${iter}_slurm.sh
done

iter=11
sqCreateScript -n 67 -c 1 -m 10G -w 5:00:00 -q scavenge -N gnova_ns.task_part_part${iter} --logdir=gnova_ns.task_part_part${iter}.log gnova_ns.task_part${iter} > gnova_ns.task_part${iter}_slurm.sh
sbatch < gnova_ns.task_part${iter}_slurm.sh
