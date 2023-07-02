import os, sys
PATH = "./MPI/'VS_P_Working'/vortex-shedding"
 
class MPI_Test:
    def __init__(self, imax, jmax, n) -> None:
        self.imax = imax
        self.jmax = jmax
        self.n = n
        self.output_filename = f"Tests/{self.imax}x{self.jmax}-MPI-{self.n}"
        self.command = f"mpirun -n {self.n} {PATH}/vortex -x {self.imax} -y {self.jmax} -f 1000 -o {self.output_filename}"

    def run(self):
        print(f"running test -> {self.output_filename}")
        os.system(self.command)

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
    for n in n_vals:
        for domain in domain_configurations:
            MPI_Test(domain[0],domain[1],n).run()

main()
