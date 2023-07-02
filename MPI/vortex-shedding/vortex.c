#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <errno.h>
#include <fcntl.h>
#include <math.h>
#include <mpi.h>

#include "data.h"
#include "vtk.h"
#include "setup.h"
#include "boundary.h"
#include "args.h"

/**
 * @brief Computation of tentative velocity field (f, g)
 * 
 */
void compute_tentative_velocity(int rank, int size) {
    for (int i = 1; i < imax; i++) {
        for (int j = 1; j < jmax+1; j++) {
            /* only if both adjacent cells are fluid cells */
            if ((flag[i][j] & C_F) && (flag[i+1][j] & C_F)) {
                double du2dx = ((u[i][j] + u[i+1][j]) * (u[i][j] + u[i+1][j]) +
                                y * fabs(u[i][j] + u[i+1][j]) * (u[i][j] - u[i+1][j]) -
                                (u[i-1][j] + u[i][j]) * (u[i-1][j] + u[i][j]) -
                                y * fabs(u[i-1][j] + u[i][j]) * (u[i-1][j]-u[i][j]))
                                / (4.0 * delx);
                double duvdy = ((v[i][j] + v[i+1][j]) * (u[i][j] + u[i][j+1]) +
                                y * fabs(v[i][j] + v[i+1][j]) * (u[i][j] - u[i][j+1]) -
                                (v[i][j-1] + v[i+1][j-1]) * (u[i][j-1] + u[i][j]) -
                                y * fabs(v[i][j-1] + v[i+1][j-1]) * (u[i][j-1] - u[i][j]))
                                / (4.0 * dely);
                double laplu = (u[i+1][j] - 2.0 * u[i][j] + u[i-1][j]) / delx / delx +
                                (u[i][j+1] - 2.0 * u[i][j] + u[i][j-1]) / dely / dely;
   
                f[i][j] = u[i][j] + del_t * (laplu / Re - du2dx - duvdy);
            } else {
                f[i][j] = u[i][j];
            }
        }
    }
    for (int i = 1; i < imax+1; i++) {
        for (int j = 1; j < jmax; j++) {
            /* only if both adjacent cells are fluid cells */
            if ((flag[i][j] & C_F) && (flag[i][j+1] & C_F)) {
                double duvdx = ((u[i][j] + u[i][j+1]) * (v[i][j] + v[i+1][j]) +
                                y * fabs(u[i][j] + u[i][j+1]) * (v[i][j] - v[i+1][j]) -
                                (u[i-1][j] + u[i-1][j+1]) * (v[i-1][j] + v[i][j]) -
                                y * fabs(u[i-1][j] + u[i-1][j+1]) * (v[i-1][j]-v[i][j]))
                                / (4.0 * delx);
                double dv2dy = ((v[i][j] + v[i][j+1]) * (v[i][j] + v[i][j+1]) +
                                y * fabs(v[i][j] + v[i][j+1]) * (v[i][j] - v[i][j+1]) -
                                (v[i][j-1] + v[i][j]) * (v[i][j-1] + v[i][j]) -
                                y * fabs(v[i][j-1] + v[i][j]) * (v[i][j-1] - v[i][j]))
                                / (4.0 * dely);
                double laplv = (v[i+1][j] - 2.0 * v[i][j] + v[i-1][j]) / delx / delx +
                                (v[i][j+1] - 2.0 * v[i][j] + v[i][j-1]) / dely / dely;

                g[i][j] = v[i][j] + del_t * (laplv / Re - duvdx - dv2dy);
            } else {
                g[i][j] = v[i][j];
            }
        }
    }

    /* f & g at external boundaries */
    for (int j = 1; j < jmax+1; j++) {
        f[0][j]    = u[0][j];
        f[imax][j] = u[imax][j];
    }
    for (int i = 1; i < imax+1; i++) {
        g[i][0]    = v[i][0];
        g[i][jmax] = v[i][jmax];
    }
}


/**
 * @brief Calculate the right hand side of the pressure equation 
 * 
 */
void compute_rhs(int rank, int size) {
    for (int i = 1; i < imax+1; i++) {
        for (int j = 1;j < jmax+1; j++) {
            if (flag[i][j] & C_F) {
                /* only for fluid and non-surface cells */
                rhs[i][j] = ((f[i][j] - f[i-1][j]) / delx + 
                             (g[i][j] - g[i][j-1]) / dely)
                             / del_t;
            }
        }
    }
}

double get_p0(int rank, int size){
    double p0 = 0.0;
    /* Calculate sum of squares */
    for (int i = 1; i < imax+1; i++) {
        for (int j = 1; j < jmax+1; j++) {
            if (flag[i][j] & C_F) { p0 += p[i][j] * p[i][j]; }
        }
    }
   
    p0 = sqrt(p0 / fluid_cells); 
    if (p0 < 0.0001) { p0 = 1.0; }

    return p0;
}

double calc_res(int rank, int size, double p0){

    double rdx2 = 1.0 / (delx * delx);
    double rdy2 = 1.0 / (dely * dely);
    double beta_2 = -omega / (2.0 * (rdx2 + rdy2));

   
    int number_of_rows_to_send = p_size_x/size;
    int rem = p_size_x % size;
    int sum = 0;

    int* sendcounts = malloc(sizeof(int)*size*2);
    int* displs = malloc(sizeof(int)*size*2);
    double* rec_buff= malloc(sizeof(double)*p_size_y*p_size_x*2);
    
    /* Calculate Displacements and Send Counts */
    for (int i = 0; i < size; i++) {
        sendcounts[i] = number_of_rows_to_send * p_size_y;
        if (rem > 0) {
            sendcounts[i]+= p_size_y;
            rem--;
        }

        displs[i] = sum;
        sum += sendcounts[i];
    }

    /* Print Scatter  */
    if(rank == 0){
        for (int i = 0; i < size; i++) {
        //    printf("sendcounts[%d] = %d\tdispls[%d] = %d\n", i, sendcounts[i], i, displs[i]);
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);

    /* Scatter */
    MPI_Barrier(MPI_COMM_WORLD);
   // printf("SCATTER\n");
    if(size > 1)
        MPI_Scatterv(&(p[0][0]), sendcounts, displs, MPI_DOUBLE, &(rec_buff[p_size_y]), p_size_x*p_size_y, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    int size_n = sendcounts[rank];
    if(rank > 0){
        double * row = calloc(p_size_y,sizeof(double));
        // printf("sending top ghost rows\n");
        // printf("%d %d\n", rank, (displs[rank] - p_size_y));
        MPI_Send(&(p[0][displs[rank]]),p_size_y,MPI_DOUBLE, rank-1, 02, MPI_COMM_WORLD);
        MPI_Recv(&(rec_buff[displs[rank] + p_size_y]),p_size_y,MPI_DOUBLE, rank-1, 03, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      //  printf("%d received %f\n", rank, row[0]);
        
    }
    if(rank < (size - 1)){
        //printf("sending bottom ghost rows\n");
    //    printf("VAL %d %d %f\n", rank, (displs[rank+1]), data[0][displs[rank+1]]);
        MPI_Recv(&(rec_buff[sendcounts[rank]+p_size_y]),p_size_y,MPI_DOUBLE, rank+1, 02, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Send(&(p[0][displs[rank+1]-p_size_y]),p_size_y,MPI_DOUBLE, rank+1, 03, MPI_COMM_WORLD);
    }
    
    int start = 1; // 1st row (0 + 1)
    int end = sendcounts[rank] / p_size_y + start; 

    /* Send Halo*/
  //  int size_n = sendcounts[rank];
    if(rank == 0)
        start = 1;
    if(rank == size - 1){
      //  end -= 1;
    }
    MPI_Barrier(MPI_COMM_WORLD);
    // printf("p0 %f\n", p0);

    /* Red/Black SOR-iteration */
    int iter;
    double res = 0.0;
    double tmp_res = 0.0;
    for (iter = 0; iter < itermax; iter++) {
        for (int rb = 0; rb < 2; rb++) {

            if(rank > 0){
            double * row = calloc(p_size_y,sizeof(double));
            // printf("sending top ghost rows\n");
            // printf("%d %d\n", rank, (displs[rank] - p_size_y));
            MPI_Send(&(rec_buff[displs[rank]]),p_size_y,MPI_DOUBLE, rank-1, 02, MPI_COMM_WORLD);
            MPI_Recv(&(rec_buff[displs[rank] + p_size_y]),p_size_y,MPI_DOUBLE, rank-1, 03, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            //  printf("%d received %f\n", rank, row[0]);
                
            }
            if(rank < (size - 1)){
                //printf("sending bottom ghost rows\n");
            //    printf("VAL %d %d %f\n", rank, (displs[rank+1]), data[0][displs[rank+1]]);
                MPI_Recv(&(rec_buff[sendcounts[rank]+p_size_y]),p_size_y,MPI_DOUBLE, rank+1, 02, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
               // MPI_Send(&(p[0][displs[rank+1]-p_size_y]),p_size_y,MPI_DOUBLE, rank+1, 03, MPI_COMM_WORLD);
               MPI_Send(&(rec_buff[sendcounts[rank]]),p_size_y,MPI_DOUBLE, rank+1, 03, MPI_COMM_WORLD);
            }

            for (int a = start; a < end; a++) {
                for (int j = 1; j < p_size_y-1; j++) {
                    int i = (a-1) + displs[rank]/p_size_y;
                    if ((i + j) % 2 != rb) { continue; }
                    if (flag[i][j] == (C_F | B_NSEW)) {
                        /* five point star for interior fluid cells */
                        rec_buff[a*p_size_y + j] = (1.0 - omega) * rec_buff[a*p_size_y + j]- 
                              beta_2 * ((rec_buff[(a+1)*p_size_y + j] + rec_buff[(a-1)*p_size_y + j] ) *rdx2
                                  + (rec_buff[a*p_size_y + j+1] + rec_buff[a*p_size_y + j-1]) * rdy2
                                  - rhs[i][j]);
                    } else if (flag[i][j] & C_F) { 
                        /* modified star near boundary */

                        double eps_E = ((flag[i+1][j] & C_F) ? 1.0 : 0.0);
                        double eps_W = ((flag[i-1][j] & C_F) ? 1.0 : 0.0);
                        double eps_N = ((flag[i][j+1] & C_F) ? 1.0 : 0.0);
                        double eps_S = ((flag[i][j-1] & C_F) ? 1.0 : 0.0);

                        double beta_mod = -omega / ((eps_E + eps_W) * rdx2 + (eps_N + eps_S) * rdy2);
                        rec_buff[a*p_size_y + j] = (1.0 - omega) * rec_buff[a*p_size_y + j] -
                            beta_mod * ((eps_E * rec_buff[(a+1)*p_size_y + j] + eps_W * rec_buff[(a-1)*p_size_y + j]) * rdx2
                                + (eps_N * rec_buff[a*p_size_y + j+1] + eps_S * rec_buff[a*p_size_y + j-1]) * rdy2
                                - rhs[i][j]);
                    }
                }
            }
        }    
        
        
        /* computation of residual */
        for (int a = start; a < end; a++) {
            for (int j = 1; j < jmax+1; j++) {
                int i = (a-1) + displs[rank]/p_size_y;
                if (flag[i][j] & C_F) {
                    double eps_E = ((flag[i+1][j] & C_F) ? 1.0 : 0.0);
                    double eps_W = ((flag[i-1][j] & C_F) ? 1.0 : 0.0);
                    double eps_N = ((flag[i][j+1] & C_F) ? 1.0 : 0.0);
                    double eps_S = ((flag[i][j-1] & C_F) ? 1.0 : 0.0);

                    /* only fluid cells */
                    double add = (eps_E * (rec_buff[(a+1)*p_size_y + j] - rec_buff[a*p_size_y + j]) - 
                        eps_W * (rec_buff[a*p_size_y + j] - rec_buff[(a-1)*p_size_y + j])) * rdx2  +
                        (eps_N * (rec_buff[a*p_size_y + j+1] - rec_buff[a*p_size_y + j]) -
                        eps_S * (rec_buff[a*p_size_y + j] - rec_buff[a*p_size_y + j-1])) * rdy2  -  rhs[i][j];
                    res += add * add;
                }
            }
        }

        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Allreduce(&res,&tmp_res,1,MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);
        res = sqrt(tmp_res / fluid_cells) / p0;
        tmp_res = 0;
        /* convergence? */
        if (res < eps) break;
    }

      MPI_Gatherv(&(rec_buff[p_size_y]), sendcounts[rank], MPI_DOUBLE, &(p[0][0]), sendcounts, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    

    return res;
}

/**
 * @brief Red/Black SOR to solve the poisson equation.
 * 
 * @return Calculated residual of the computation
 * 
 */
double poisson(int rank, int size) {
    double p0  = 0.0;
    /* Get p0 value */
    if(rank == 0)
        p0= get_p0(rank, size);
    MPI_Bcast(&p0,1,MPI_DOUBLE, 0, MPI_COMM_WORLD);
    /* Red/Black SOR-iteration */
    
    return calc_res(rank, size, p0);
}


/**
 * @brief Update the velocity values based on the tentative
 * velocity values and the new pressure matrix
 */
void update_velocity(int rank, int size) {   
    for (int i = 1; i < imax-2; i++) {
        for (int j = 1; j < jmax-1; j++) {
            /* only if both adjacent cells are fluid cells */
            if ((flag[i][j] & C_F) && (flag[i+1][j] & C_F)) {
                u[i][j] = f[i][j] - (p[i+1][j] - p[i][j]) * del_t / delx;
            }
        }
    }
    
    for (int i = 1; i < imax-1; i++) {
        for (int j = 1; j < jmax-2; j++) {
            /* only if both adjacent cells are fluid cells */
            if ((flag[i][j] & C_F) && (flag[i][j+1] & C_F)) {
                v[i][j] = g[i][j] - (p[i][j+1] - p[i][j]) * del_t / dely;
            }
        }
    }
}


/**
 * @brief Set the timestep size so that we satisfy the Courant-Friedrichs-Lewy
 * conditions. Otherwise the simulation becomes unstable.
 */
void set_timestep_interval(int rank, int size) {
    /* del_t satisfying CFL conditions */
    if (tau >= 1.0e-10) { /* else no time stepsize control */
        double umax = 1.0e-10;
        double vmax = 1.0e-10; 
        
        for (int i = 0; i < imax+2; i++) {
            for (int j = 1; j < jmax+2; j++) {
                umax = fmax(fabs(u[i][j]), umax);
            }
        }

        for (int i = 1; i < imax+2; i++) {
            for (int j = 0; j < jmax+2; j++) {
                vmax = fmax(fabs(v[i][j]), vmax);
            }
        }

        double deltu = delx / umax;
        double deltv = dely / vmax; 
        double deltRe = 1.0 / (1.0 / (delx * delx) + 1 / (dely * dely)) * Re / 2.0;

        if (deltu < deltv) {
            del_t = fmin(deltu, deltRe);
        } else {
            del_t = fmin(deltv, deltRe);
        }
        del_t = tau * del_t; /* multiply by safety factor */
    }
}

/**
 * @brief The main routine that sets up the problem and executes the solving routines routines
 * 
 * @param argc The number of arguments passed to the program
 * @param argv An array of the arguments passed to the program
 * @return int The return value of the application
 */
int main(int argc, char *argv[]) {
    int rank, size;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    set_defaults(); 
    parse_args(argc, argv);
    setup();

    if (verbose) print_opts();

    allocate_arrays();
    problem_set_up();

    double res;

    /* Main loop */
    int iters = 0;
    double t;
    for (t = 0.0; t < t_end; t += del_t, iters++) {
        if (!fixed_dt)
            set_timestep_interval(rank, size);

        compute_tentative_velocity(rank, size);

        compute_rhs(rank, size);

        res = poisson(rank, size);

        update_velocity(rank, size);

        apply_boundary_conditions(rank, size);

        if ((iters % output_freq == 0)) {
            printf("[%d] Step %8d, Time: %14.8e (del_t: %14.8e), Residual: %14.8e\n" ,rank, iters, t+del_t, del_t, res);
 
            if ((!no_output) && (enable_checkpoints))
                write_checkpoint(iters, t+del_t);
        }
    } /* End of main loop */

    printf("Step %8d, Time: %14.8e, Residual: %14.8e\n", iters, t, res);
    printf("Simulation complete.\n");

    if (!no_output)
        write_result(iters, t);

    free_arrays();

    return 0;
}

