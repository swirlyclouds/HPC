import argparse
import os, sys
 
class MPI_Test:
    PATH = "./MPI/'VS_P_Working'/vortex-shedding"
    def __init__(self, imax, jmax, n) -> None:
        self.imax = imax
        self.jmax = jmax
        self.n = n
        self.output_filename = f"{self.imax}x{self.jmax}-MPI-{self.n}"
        self.command = f"mpirun -n {self.n} {self.PATH}/vortex -x {self.imax} -y {self.jmax} -f 1000 -o Tests/{self.output_filename}"

    def run(self):
        print(f"running test -> {self.output_filename}")
        os.system(self.command)

    def generate_slurm(self):
        param = "${SLURM_NTASKS}"
        header = f"""sh
#!/bin/bash
#SBATCH --job-name={self.output_filename}               # Job name
#SBATCH --mail-type=END,FAIL                 # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=maa621@york.ac.uk          # Where to send mail
#SBATCH --ntasks={self.n}                          # Run sixteen tasks...
#SBATCH --cpus-per-task=1                    # ...with one core each
#SBATCH --mem-per-cpu=4gb                  # Memory per processor
#SBATCH --time=00:10:00                      # Time limit hrs:min:sec
#SBATCH --output={self.output_filename}.log          # Standard output and error log
         
echo "Running {self.output_filename} on {param} CPU cores"
                
mpirun -n {param} {self.PATH}/vortex -x {self.imax} -y {self.jmax} -f 1000 -o {self.output_filename}
            """
        file = open(f"Jobs/{self.output_filename}.job", 'w+')
        file.write(header)
        file.close()

class OMP_Test:
    PATH = "./OMP/vortex-shedding(1)/vortex-shedding"
    def __init__(self, imax, jmax, n) -> None:
        self.imax = imax
        self.jmax = jmax
        self.n = n
        self.output_filename = f"{self.imax}x{self.jmax}-OMP-{self.n}"
        self.command = f"export OMP_NUM_TASKS={self.n} {self.PATH}/vortex -x {self.imax} -y {self.jmax} -f 1000 -o Tests/{self.output_filename}"
    
    def generate_slurm(self):
        param = "${SLURM_CPUS_PER_TASK}"
        header = f"""sh
#!/bin/bash
#SBATCH --job-name={self.output_filename}               # Job name
#SBATCH --mail-type=END,FAIL                 # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=maa621@york.ac.uk          # Where to send mail
#SBATCH --cpus-per-task={self.n}                    # ...with one core each
#SBATCH --mem-per-cpu=4gb                  # Memory per processor
#SBATCH --time=00:10:00                      # Time limit hrs:min:sec
#SBATCH --output={self.output_filename}.log          # Standard output and error log
         
echo "Running {self.output_filename} on {param} CPU cores"
export OMP_NUM_TASKS={param}                
{self.PATH}/vortex -x {self.imax} -y {self.jmax} -f 1000 -o Tests/{self.output_filename}
            """
        file = open(f"Jobs/{self.output_filename}.job", 'w+')
        file.write(header)
        file.close()

class CUDA_Test:
    PATH = "./CUDA/vortex-shedding"
    def __init__(self, imax, jmax) -> None:
        self.imax = imax
        self.jmax = jmax
        self.output_filename = f"{self.imax}x{self.jmax}-CUDA"
        self.command = f"{self.PATH}/vortex -x {self.imax} -y {self.jmax} -f 1000 -o Tests/{self.output_filename}"
    
    def generate_slurm(self):
        out = f"""#!/bin/bash
#SBATCH --job-name={self.output_filename}                    # Job name
#SBATCH --mail-type=END,FAIL                   # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=maa621@york.ac.uk          # Where to send mail
#SBATCH --ntasks=1                             # Run a single task...
#SBATCH --cpus-per-task=1                      # ...with a single CPU
#SBATCH --mem=4gb                            # Job memory request
#SBATCH --time=00:10:00                        # Time limit hrs:min:sec
#SBATCH --output={self.output_filename}.log               # Standard output and error log
#SBATCH --partition=gpu                        # Select the GPU nodes...
#SBATCH --gres=gpu:1                           # ...and a single GPU
  
module load system/CUDA/11.7.0
 
echo `date`: executing gpu_test on host $HOSTNAME with $SLURM_CPUS_ON_NODE cpu cores
echo
cudaDevs=$(echo $CUDA_VISIBLE_DEVICES | sed -e 's/,/ /g')
echo I can see GPU devices $CUDA_VISIBLE_DEVICES
echo
 
{self.PATH}/vortex -x {self.imax} -y {self.jmax} -f 1000 -o Tests/{self.output_filename}"""
        file = open(f"Jobs/{self.output_filename}.job", 'w+')
        file.write(out)
        file.close()
        
def main():
    # run tests
    # different sizes :
    """
        Original (128 x 512)
        Small: 64x256
        Medium: 256 x 1024
        Large: 512 x 2048
        Collosall: 1024 x 4096
    """

    n_vals = range(1,4)
    domain_configurations = [
        (512, 128),
        (256, 64),
        (1024,256),
        (2048,512),
        (4096,1024)
    ]
    d = domain_configurations[0]
    n_ = 4
   # MPI_Test(d[0],d[1],n_).generate_slurm()
    OMP_Test(d[0],d[1],n_).generate_slurm()
   # CUDA_Test(d[0],d[1],n_).generate_slurm()
    for n in n_vals:
        for domain in domain_configurations:
            pass

main()
