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
		std::cout << "nprocs:" << nprocs << ";rank:" << rank << std::endl << std::flush;
		dist_sort_size_t global_N;
		MPI_Allreduce(&myDataCount, &global_N, 1, MPI_TYPE_DIST_SORT_SIZE_T, MPI_SUM, MPI_COMM_WORLD);
		std::cout << "global_N:" << global_N << ";myDataCount:" << myDataCount << ";rank:" << rank << std::endl << std::flush;
		if (rank < nprocs-1) {
				*rCount = ceil(global_N/nprocs);
		} else {
				*rCount = global_N - ceil(global_N/nprocs) * (nprocs-1);
		}
		*rebalancedData = (dist_sort_t*)malloc((*rCount)*sizeof(dist_sort_t));

		dist_sort_size_t global_count = 0;
    MPI_Exscan(&myDataCount, &global_count, 1, MPI_TYPE_DIST_SORT_SIZE_T, MPI_SUM, MPI_COMM_WORLD);
		std::cout << "global_count:" << global_count << ";rank:" << rank << std::endl << std:flush;
    MPI_Win win;
    MPI_Win_create(*rebalancedData, (*rCount) * sizeof(dist_sort_t), sizeof(dist_sort_t), MPI_INFO_NULL, MPI_COMM_WORLD, &win);
    MPI_Win_fence(MPI_MODE_NOPRECEDE, win); //fence - there are no epochs before this
    int i = 0;
		dist_sort_size_t max_size = ceil(global_N/nprocs);
    while (i < myDataCount)
    {
        int target_rank = int((global_count+i) / max_size);
				if (target_rank == 4) {
						std::cout << global_count+i << "," << max_size << std::endl << std::flush;
				}
				// std::cout << target_rank << "," << global_count << "," << i << "," << max_size <<std::endl << std::flush;
        int displacement = (global_count+i) % max_size;
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
