sh
#!/bin/bash
#SBATCH --job-name=Tests/512x128-MPI-4               # Job name
#SBATCH --mail-type=END,FAIL                 # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=maa621@york.ac.uk          # Where to send mail
#SBATCH --ntasks=4                          # Run sixteen tasks...
#SBATCH --cpus-per-task=1                    # ...with one core each
#SBATCH --mem-per-cpu=4gb                  # Memory per processor
#SBATCH --time=00:10:00                      # Time limit hrs:min:sec
#SBATCH --output=Tests/512x128-MPI-4.log          # Standard output and error log
         
echo "Running Tests/512x128-MPI-4 on ${SLURM_NTASKS} CPU cores"
                
mpirun -n ${SLURM_NTASKS} ./MPI/'VS_P_Working'/vortex-shedding/vortex -x 512 -y 128 -f 1000 -o Tests/512x128-MPI-4
            