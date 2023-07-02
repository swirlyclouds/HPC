import argparse
import os, sys
 
class MPI_Test:
    PATH = "./MPI/'VS_P_Working'/vortex-shedding"
    def __init__(self, imax, jmax, n) -> None:
        self.imax = imax
        self.jmax = jmax
        self.n = n
        self.output_filename = f"Tests/{self.imax}x{self.jmax}-MPI-{self.n}"
        self.command = f"mpirun -n {self.n} {self.PATH}/vortex -x {self.imax} -y {self.jmax} -f 1000 -o {self.output_filename}"

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
        file = open(f"{self.output_filename}.job", 'w+')
        file.write(header)
        file.close()

class OMP_Test:
    PATH = "./OMP/vortex-shedding(1)/vortex-shedding"
    def __init__(self, imax, jmax, n) -> None:
        self.imax = imax
        self.jmax = jmax
        self.n = n
        self.output_filename = f"Tests/{self.imax}x{self.jmax}-OMP-{self.n}"
        self.command = f"export OMP_NUM_TASKS={self.n} {self.n} {self.PATH}/vortex -x {self.imax} -y {self.jmax} -f 1000 -o {self.output_filename}"
    
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
OMP_NUM_TASKS={param}                
{self.PATH}/vortex -x {self.imax} -y {self.jmax} -f 1000 -o {self.output_filename}
            """
        file = open(f"{self.output_filename}.job", 'w+')
        file.write(header)
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
    MPI_Test(d[0],d[1],n_).generate_slurm()
    for n in n_vals:
        for domain in domain_configurations:
            pass

main()
