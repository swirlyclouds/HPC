#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

void perform_work(int thread_id) {
	printf("Thread %d has started\n", thread_id);
}

int main(int argc, char *argv[]) {
	int my_val = 0;

	// create threads
	#pragma omp parallel
	{
		int thread_id = omp_get_thread_num();
		perform_work(thread_id);
		my_val += thread_id;
	}

	printf("The final value of my_val is: %d\n", my_val);

	exit(0);
}
