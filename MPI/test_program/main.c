#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

double **alloc_2d_array(int m, int n) {
  	double **x;
  	int i;

  	x = (double **)malloc(m*sizeof(double *));
  	x[0] = (double *)calloc(m*n,sizeof(double));
  	for ( i = 1; i < m; i++ )
    	x[i] = &x[0][i*n];
	return x;
}

void addArrayInsideVals(int rank, int size){
    int p_size_x = 5;
    int p_size_y = 5;
    double** p = alloc_2d_array(5,5);

    for(int i = 0; i < p_size_x * p_size_y; i++)
        p[0][i] = (i + 1);


    int sum = 0;
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
    }
    for(int i = 0; i < 1; i++){
    MPI_Barrier(MPI_COMM_WORLD);
    //printf("SCATTER\n");
    if(size > 1)
        MPI_Scatterv(&(p[0][0]), sendcounts, displs, MPI_DOUBLE, &(rec_buff[p_size_y]), p_size_x*p_size_y, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    int size_n = sendcounts[rank];
    if(rank > 0){
        double * row = calloc(p_size_y,sizeof(double));
        //printf("sending top ghost rows\n");
        //printf("%d %d\n", rank, (displs[rank] - p_size_y));
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
    MPI_Barrier(MPI_COMM_WORLD);
    
         
    
    int start = displs[rank];
    int end =  sendcounts[rank];
    /* Send Halo*/
   // int size_n = sendcounts[rank];
    start = 1; // 1st row (0 + 1)
    end = sendcounts[rank] / p_size_y + start; // last row (secon
    // for(int i = start; i < end; i++){
    //     printf("%d| ", rank);
    //     for(int j = 1; j < p_size_y - 1; j++){
    //         printf("%f, ", rec_buff[i*p_size_y + j]);
    //     }
    //     printf("\n");
    // }

    if(rank == 0){ start += 1; }
    if(rank == size - 1) end -= 1;
    for(int i = start; i < end; i++){
        for(int j = 1; j < p_size_y - 1; j++){
            int a = (i-1) + displs[rank]/p_size_y;
          //  rec_buff[i*p_size_y + j] += p[0][(i-2) * p_size_y + j + displs[rank]] + p[0][(i) * p_size_y + j + displs[rank]] + p[0][(i-1) * p_size_y + j+1 + displs[rank]] + p[0][(i-1) * p_size_y + j - 1 + displs[rank]];            

           // rec_buff[i*p_size_y + j] += p[0][a - p_size_y + j] + p[0][a + p_size_y + j] + p[0][a + j + 1] + p[0][a + j - 1];   
            rec_buff[i*p_size_y + j] += p[a-1][j] + p[a+1][j] + p[a][j + 1] + p[a][j - 1];   
        }
    }
    
    MPI_Barrier(MPI_COMM_WORLD);
    start = 1; // 1st row (0 + 1)
    end = sendcounts[rank] / p_size_y + start;
    printf("--\n");
    for(int i = start; i < end; i++){
        printf("%d| ", rank);
        for(int j = 0; j < p_size_y; j++){
            printf("%f, ", rec_buff[i*p_size_y + j]);
        }
        printf("\n");
    }

    MPI_Gatherv(&(rec_buff[p_size_y]),sendcounts[rank],MPI_DOUBLE,p[0],sendcounts,displs,MPI_DOUBLE,0,MPI_COMM_WORLD);
    if(rank == 0){
        for(int i  = 0; i < 5; i++){
            for(int j = 0; j < 5; j++){
                printf("%f, " , p[i][j]);
            }
            printf("\n");
        }
    }
    }
}

/* Experiment 2: SendRecv */

void sendRecvGhost(int rank, int size){
    int p_size_x = 514;
    int p_size_y = 128;
    
    /* Create Array p */
    double** p = alloc_2d_array(p_size_x, p_size_y);
    
    /* Initialise P values */
    for(int i = 0; i < p_size_x * p_size_y; i++)
        p[0][i] = (i + 1);

    /* Determine how to scatter among processes */

    int sum = 0;
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
            if(rank == 0)  printf("%d : %d | %d\n", i, displs[i], sendcounts[i]);
        }
    if(rank == 1){
        for(int i = 0; i < p_size_x * p_size_y; i++){
         //   printf("%f, ", p[0][i]);
        }
        printf("\n");
    }
    /* Scatter among processes */

        // OR DON'T :P

    /* Trade Ghost Cells */


    int peer = (rank == 0) ? 1 : 0;

    if(rank == 0){
        // increment all values by 2
        for(int i = 0; i < 25; i++){
            p[0][i] += 0;
        }
    }

    MPI_Status status;
    //swap only the first row
    // MPI_Sendrecv(&(p[0][displs[rank]]), p_size_y, MPI_DOUBLE, peer, 1, // Send first row of rank to peer
    
    //             &(p[0][displs[peer]]), p_size_y, MPI_DOUBLE, peer, 1, MPI_COMM_WORLD,&status); // Place for from peer in first row of rank
    //     printf("MPI process %d received value %f from rank %d, with tag %d and error code %d.\n", 

    //            rank,

    //            p[0][displs[rank]],

    //            status.MPI_SOURCE,

    //            status.MPI_TAG,

    //            status.MPI_ERROR);

    /* Perform Operations */
    // if(rank == 1){
    //     for(int i = 0; i < p_size_x * p_size_y; i++){
    //         printf("%f, ", p[0][i]);
    //     }
    //     printf("\n");
    // }
    int start = displs[rank];
    int end = displs[rank] + sendcounts[rank];
    end /= p_size_y;
    start /= p_size_y;
    if(rank == 0){
        start += 1;
    }
    if(rank == size - 1){
        end -= 1;
    }
    for(int rb = 0; rb < 2; rb++){
        /* Replace Top Halo */
        //  MPI_Sendrecv(&(p[0][end]), p_size_y, MPI_DOUBLE, peer, 1, // Send first row of rank to peer
        //              &(p[0][start - 1]), p_size_y, MPI_DOUBLE, peer, 1, MPI_COMM_WORLD,&status); // Place for from peer in first row of rank
        
        for(int i = start; i < end; i++){
            for(int j = 1; j < p_size_y-1; j++){
                if ((i + j) % 2 != rb) { continue; }
                else{
                // if(rank == 1){
                //     printf("%d %d %d : %f = %f + %f %f %f %f\n",i, j, rb, (p[i][j] + p[i+1][j] + p[i-1][j] + p[i][j-1] + p[i][j+1]), p[i][j], p[i+1][j], p[i-1][j], p[i][j-1],p[i][j+1]);
                // }
                p[i][j] += (p[i+1][j] + p[i-1][j] + p[i][j-1] + p[i][j+1]);
                }
             
            }
        }
        if(size > 1){
        if(rank == 0){
        MPI_Send(&(p[end-1][0]), p_size_y, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD);
        MPI_Recv(&p[end][0], p_size_y, MPI_DOUBLE, 1, 3, MPI_COMM_WORLD, &status);
        }
        }
        if(rank == 1)
        {
            MPI_Send(&p[start][0], p_size_y, MPI_DOUBLE, 0, 3, MPI_COMM_WORLD);
            MPI_Recv(&(p[start - 1][0]), p_size_y, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
            if(size == 3){
                printf("sending...\n");
                MPI_Send(&(p[end-1][0]), p_size_y, MPI_DOUBLE, 2, 1, MPI_COMM_WORLD);
                MPI_Recv(&p[end][0], p_size_y, MPI_DOUBLE, 2, 2, MPI_COMM_WORLD, &status);
            }
        }
        if(rank == 2){
            MPI_Recv(&(p[start - 1][0]), p_size_y, MPI_DOUBLE, 1,1, MPI_COMM_WORLD, &status);
            MPI_Send(&p[start][0], p_size_y, MPI_DOUBLE, 1, 2, MPI_COMM_WORLD);
        }
        
        
    }



    MPI_Barrier(MPI_COMM_WORLD);

    /* Return Output*/
    for(int r = 0 ; r < size; r++){
    
    MPI_Barrier(MPI_COMM_WORLD);
    // if(rank == r){
    //     printf("r: %d\n", r);
    //     for(int i = 0; i < 5; i++){
    //         for(int j = 0; j < 5; j++){
    //             printf("%f, ", p[i][j]);
    //         }
    //         printf("\n");
    //     }
    // }
    MPI_Barrier(MPI_COMM_WORLD);


    }
    /* Gather scattered data */
    MPI_Gatherv(&p[0][displs[rank]],sendcounts[rank],MPI_DOUBLE, &(p[0][0]),sendcounts,displs,MPI_DOUBLE, 0, MPI_COMM_WORLD);

    MPI_Barrier(MPI_COMM_WORLD);
    
    MPI_Barrier(MPI_COMM_WORLD);
    if(rank == 0){
        printf("r: %d\n", 0);
        for(int i = 0; i < 5; i++){
            for(int j = 0; j < 5; j++){
                printf("%f, ", p[i][j]);
            }
            printf("\n");
        }
    }

}


void sendrecv_testA(int my_rank, int size){

    int buffer_send = (my_rank == 0) ? 12345 : 67890;

    int buffer_recv;

    int tag_send = 0;

    int tag_recv = tag_send;

    int peer = (my_rank == 0) ? 1 : 0;

 

    // Issue the send + receive at the same time

    printf("MPI process %d sends value %d to MPI process %d.\n", my_rank, buffer_send, peer);

    MPI_Sendrecv(&buffer_send, 1, MPI_INT, peer, tag_send,

                 &buffer_recv, 1, MPI_INT, peer, tag_recv, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    printf("MPI process %d received value %d from MPI process %d.\n", my_rank, buffer_recv, peer);
}

int main(int argc, char **argv)
{
    // Initialize the MPI environment
    MPI_Init(&argc, &argv);
  
    // Get the number of processes ssociated with the communicator
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // Get the rank of the calling process
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    // Get the name of the processor
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int name_len;
    MPI_Get_processor_name(processor_name, &name_len);

    printf("Hello world from process %s with rank %d out of %d processors\n", processor_name, world_rank, world_size);
    //addArrayInsideVals(world_rank, world_size);
    sendRecvGhost(world_rank, world_size);
    // Finalize: Any resources allocated for MPI can be freed
    MPI_Finalize();
}