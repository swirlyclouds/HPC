
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include <getopt.h>
#include <errno.h>
#include <fcntl.h>
#include <math.h>
#include <time.h>

double **alloc_2d_array(int m, int n) {
  	double **x;
  	int i;

  	x = (double **)malloc(m*sizeof(double *));
  	x[0] = (double *)calloc(m*n,sizeof(double));
  	for ( i = 1; i < m; i++ )
    	x[i] = &x[0][i*n];
	return x;
}

#define p_size_y 5
#define p_size_x 5
void test_A(int rank, int size){
    int sum = 0;
    int rem = (p_size_y*p_size_y)%size;
    int *sendcounts = malloc(sizeof(int)*size);
    int *displs = malloc(sizeof(int)*size);
    double* rec_buff= malloc(sizeof(double)*p_size_y*p_size_y);
    double** data = alloc_2d_array(p_size_y,p_size_y);
    
    // calculate send counts and displacements
    for (int i = 0; i < size; i++) {
        sendcounts[i] = (p_size_y*p_size_y)/size;
        if (rem > 0) {
            sendcounts[i]++;
            rem--;
        }

        displs[i] = sum;
        sum += sendcounts[i];
    }
    
    if(rank == 0){
        for (int i = 0; i < size; i++) {
            printf("sendcounts[%d] = %d\tdispls[%d] = %d\n", i, sendcounts[i], i, displs[i]);
        }
        
    }
    for(int i = 0; i < p_size_y; i++){
        for(int j = 0; j < p_size_y; j++){
            data[i][j] = (i+ 1) * p_size_y + j;
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Scatterv(&(data[0][0]), sendcounts, displs, MPI_DOUBLE, &(rec_buff[0]), 100, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    printf("%d: ", rank);
    for (int i = 0; i < sendcounts[rank]; i++) {
        printf("%f\t", rec_buff[i]);
    }
    printf("\n");

    if(rank == 0){
        printf("sending ghost rows\n");
        for(int d = 1; d < size; d++){
            printf("%d %d\n", d, (displs[d] - 1));
            MPI_Send(&(data[0][0]) + sizeof(double)*(displs[d] - 1),sendcounts[d],MPI_DOUBLE, d, 02, MPI_COMM_WORLD);
        }
    }
    else{

    }

    MPI_Finalize();

    free(sendcounts);
    free(displs);
    free(rec_buff);
    free(data[0]);
    free(data);

}

int main(int argc, char *argv[]) {
    int rank, size;
    int sum = 0;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Datatype newtype;
    double** data = alloc_2d_array(p_size_x,p_size_y);
    double** data2 = alloc_2d_array(p_size_x,p_size_y);
    int number_of_rows_to_send = p_size_x/size;
    int rem = (p_size_x)%size;
    int *sendcounts = malloc(sizeof(int)*size*2);
    int *displs = malloc(sizeof(int)*size*2);
    double* rec_buff= malloc(sizeof(double)*p_size_y*p_size_x*2);
    
    // calculate send counts and displacements
    for (int i = 0; i < size; i++) {
        sendcounts[i] = number_of_rows_to_send * p_size_y;
        if (rem > 0) {
            sendcounts[i]+= p_size_y;
            rem--;
        }

        displs[i] = sum;
        sum += sendcounts[i];
    }
    if(rank == 0){
        for (int i = 0; i < size; i++) {
            printf("sendcounts[%d] = %d\tdispls[%d] = %d\n", i, sendcounts[i], i, displs[i]);
        }
        
    for(int i = 0; i < p_size_x; i++){
        for(int j = 0; j < p_size_y; j++){
            data[i][j] = (i) * p_size_y + j + 1;
            printf("%f, ", data[i][j]);
        }
        printf("\n");
    }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    printf("SCATTER\n");
    if(size > 1)
        MPI_Scatterv(&(data[0][0]), sendcounts, displs, MPI_DOUBLE, &(rec_buff[p_size_y]), p_size_x*p_size_y, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    else{
        memcpy(&(rec_buff[p_size_y]),&(data[0][0]),p_size_y*p_size_x*sizeof(double));    
    }
    int size_n = sendcounts[rank];
    if(rank > 0){
        double * row = calloc(p_size_y,sizeof(double));
        printf("sending top ghost rows\n");
        printf("%d %d\n", rank, (displs[rank] - p_size_y));
        MPI_Send(&(data[0][displs[rank]]),p_size_y,MPI_DOUBLE, rank-1, 02, MPI_COMM_WORLD);
        MPI_Recv(&(rec_buff[displs[rank] - p_size_y]),p_size_y,MPI_DOUBLE, rank-1, 03, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      //  printf("%d received %f\n", rank, row[0]);
        
    }
    if(rank < (size - 1)){
        printf("sending bottom ghost rows\n");
    //    printf("VAL %d %d %f\n", rank, (displs[rank+1]), data[0][displs[rank+1]]);
        MPI_Recv(&(rec_buff[sendcounts[rank]+p_size_y]),p_size_y,MPI_DOUBLE, rank+1, 02, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Send(&(data[0][displs[rank+1]-p_size_y]),p_size_y,MPI_DOUBLE, rank+1, 03, MPI_COMM_WORLD);
    }
    //printf("OUTPUT\n");
    //printf("%d: ", rank);
    
    MPI_Barrier(MPI_COMM_WORLD);
    for (int i = 1; i < sendcounts[rank]/p_size_y + 1; i++) {
        for (int j = 1; j < p_size_y-1; j++){
          //  printf("%f (%f, %f) (%f, %f)", rec_buff[i*p_size_y + j], rec_buff[(i-1)*p_size_y +j],rec_buff[(i+1)*p_size_y +j],rec_buff[i*p_size_y + j+1],rec_buff[i*p_size_y + j-1]);
            rec_buff[i*p_size_y + j] += rec_buff[(i-1)*p_size_y + j] + rec_buff[(i+1)*p_size_y + j] + rec_buff[(i)*p_size_y + j + 1] + rec_buff[(i)*p_size_y + j - 1];
        }
    }
    //printf("\n");
    MPI_Barrier(MPI_COMM_WORLD);
  //  printf("GATHER\n");
   // if(size > 1)
    MPI_Gatherv(&(rec_buff[p_size_y]), sendcounts[rank], MPI_DOUBLE, &(data2[0][0]), sendcounts, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    MPI_Barrier(MPI_COMM_WORLD);
    if(rank == 0){
        printf("\n\nOUTPUT: \n");
        
        for(int i = 0; i < p_size_x; i++){
            for(int j = 0; j < p_size_y; j++){
                printf("%f, ", data2[i][j]);
            }
            printf("\n");
        }
    }
    // printf("GATHER\n");
    // if(rank == 1&& size > 1){
    //     printf("\n\nOUTPUT: \n");
    //     for(int i = 0; i < p_size_x; i++){
    //         for(int j = 0; j < p_size_y; j++){
    //             printf("%f, ", data[i][j]);
    //         }
    //     printf("\n");
    //     }
    // }
    printf("GATHER\n");
    if(rank == 1 && size > 1){
    // printf("out");
    // for(int i = 0; i < sendcounts[rank] +2*p_size_y; i++){
    //     printf("%d.) %f, \n", i, rec_buff[i]);
    //     }
     }
    MPI_Finalize();

    free(sendcounts);
    free(displs);
    free(rec_buff);
    //free(data[0]);
   // free(data);
    return 0;
}
