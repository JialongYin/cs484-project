#include <cmath>
#include <algorithm>
#include <iostream>
#include <utility>
#include <queue>
#include <vector>

#include <cassert>
#include <cstring>

#include "basic_defs.h"
#include "databasics.hpp"
#include "solution.hpp"



void rebalance(const dist_sort_t *data, const dist_sort_size_t myDataCount, dist_sort_t **rebalancedData, dist_sort_size_t *rCount) {
/*
	See the header file ('solution.hpp') for Doxygen docstrings explaining this function and its parameters.
*/
		int rank, nprocs;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

		dist_sort_size_t global_N;
		MPI_Allreduce(&myDataCount, &global_N, 1, MPI_TYPE_DIST_SORT_SIZE_T, MPI_SUM, MPI_COMM_WORLD);
		if (rank < nprocs-1) {
				*rCount = ceil(float(global_N)/float(nprocs));
		} else {
				*rCount = global_N - ceil(float(global_N)/float(nprocs)) * (nprocs-1);
		}
		*rebalancedData = (dist_sort_t*)malloc((*rCount)*sizeof(dist_sort_t));


		dist_sort_size_t global_count = 0;
    MPI_Exscan(&myDataCount, &global_count, 1, MPI_TYPE_DIST_SORT_SIZE_T, MPI_SUM, MPI_COMM_WORLD);
    MPI_Win win;
    MPI_Win_create(*rebalancedData, (*rCount) * sizeof(dist_sort_t), sizeof(dist_sort_t), MPI_INFO_NULL, MPI_COMM_WORLD, &win);
    MPI_Win_fence(MPI_MODE_NOPRECEDE, win); //fence - there are no epochs before this
    dist_sort_size_t i = 0;
		dist_sort_size_t max_size = ceil((float)global_N/(float)nprocs);
    while (i < myDataCount)
    {
        int target_rank = int((float)(global_count+i) / (float)max_size);
        dist_sort_size_t displacement = (global_count+i) % max_size;
        if (target_rank != rank) {
            MPI_Put(&(data[i]), 1, MPI_TYPE_DIST_SORT_T, target_rank, displacement, 1, MPI_TYPE_DIST_SORT_T, win);
        } else {
            (*rebalancedData)[displacement] = data[i];
        }
        ++i;
    }
    MPI_Win_fence(0, win);
    MPI_Win_fence(MPI_MODE_NOSUCCEED, win);
}





void findSplitters(const dist_sort_t *data, const dist_sort_size_t data_size, dist_sort_t *splitters, dist_sort_size_t *counts, int numSplitters) {
/*
	See the header file ('solution.hpp') for Doxygen docstrings explaining this function and its parameters.
*/
		int rank, nprocs;
		MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		dist_sort_size_t global_N;
		MPI_Allreduce(&data_size, &global_N, 1, MPI_TYPE_DIST_SORT_SIZE_T, MPI_SUM, MPI_COMM_WORLD);

		dist_sort_size_t *splitter_index;
		if (rank == 0) {
			splitter_index = (dist_sort_size_t*)malloc(numSplitters*sizeof(dist_sort_size_t));
			dist_sort_size_t interval = ceil(float(data_size)/float(numSplitters))
			for (dist_sort_size_t i = 0; i < data_size; ++i) {
					if ((i+1) % interval == 0 || i == data_size-1) {
							splitters[i/interval] = data[i];
							splitter_index[i/interval] = i;
					}
			}
		}

		while (true) {
				MPI_Bcast(splitters, numSplitters, MPI_TYPE_DIST_SORT_T, 0, MPI_COMM_WORLD);

				memset(counts, 0, numSplitters);
				dist_sort_size_t j = 0;
				for (dist_sort_size_t i = 0; i < data_size; ++i) {
						if (data[i] <= splitters[j]) {
								counts[j] += 1;
						} else {
								j += 1;
						}
				}

				dist_sort_size_t *counts_buffer = NULL;
				if (rank == 0) {
						counts_buffer = (dist_sort_size_t*)malloc(nProcs*numSplitters*sizeof(dist_sort_size_t));
				}
				MPI_Gather(counts, numSplitters, MPI_TYPE_DIST_SORT_SIZE_T, counts_buffer, numSplitters, MPI_TYPE_DIST_SORT_SIZE_T, 0, MPI_COMM_WORLD);

				if (rank == 0) {
						memset(counts, 0, numSplitters);
						for (dist_sort_size_t i = 0; i < nProcs*numSplitters; ++i) {
								counts[i%numSplitters] += counts_buffer[i];
						}
						dist_sort_size_t prefix_counts = 0;
						done = true;
						for (dist_sort_size_t i = 0; i < numSplitters; ++i) {
								prefix_counts += counts[i];
								if (i < ceil((float)(prefix_counts)/(float)global_N*numSplitters)) {
										--(splitter_index[i]);
										done = false;
								} else if (i > ceil((float)(prefix_counts)/(float)global_N*numSplitters)) {
										++(splitter_index[i]);
										done = false;
								}
						}
						if (done) {
								free(counts_buffer);
								break;
						}
						for (dist_sort_size_t i = 0; i < numSplitters; ++i) {
									dist_sort_size_t index = splitter_index[i];
									splitters[i] = data[index];
						}
				}
		}
}

void moveData(const dist_sort_t *const sendData, const dist_sort_size_t sDataCount,
		dist_sort_t **recvData, dist_sort_size_t *rDataCount,
		const dist_sort_t *const splitters, const dist_sort_t *const counts, int numSplitters) {
/*
	See the header file ('solution.hpp') for Doxygen docstrings explaining this function and its parameters.
*/


}

void sort(dist_sort_t *data, dist_sort_size_t size) {
	// You are welcome to use this sort function.
	// If you wish to implement your own, you can do that too.
	// Don't use bubblesort.
	std::sort(data,data+size);
}
