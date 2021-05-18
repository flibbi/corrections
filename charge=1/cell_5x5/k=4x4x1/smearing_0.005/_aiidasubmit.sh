#!/bin/bash 
#SBATCH --no-requeue
#SBATCH --job-name="aiida-150850"
#SBATCH --get-user-env
#SBATCH --output=_scheduler-stdout.txt
#SBATCH --error=_scheduler-stderr.txt
#SBATCH --nodes=3
#SBATCH --ntasks-per-node=28
#SBATCH --time=20:00:00


'srun' '-n' '84' '/home/libbi/qe-6.1.0/bin/pw.x' '-npool' '2' '-in' 'aiida.in'  > 'aiida.out' 
