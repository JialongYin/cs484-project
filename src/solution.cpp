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
		// uint64_t DEBUG = 10000000000000000;
		int rank, nprocs;
		MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);

		// Find global number of keys in each process
		dist_sort_size_t global_N;
		MPI_Allreduce(&data_size, &global_N, 1, MPI_TYPE_DIST_SORT_SIZE_T, MPI_SUM, MPI_COMM_WORLD);
		// Find global maximum in rank0 process
		dist_sort_t max_value = 0;
		for (dist_sort_size_t i = 0; i < data_size; ++i) {
				if (data[i] > max_value) {
						max_value = data[i];
				}
		}
		dist_sort_t global_max;
		MPI_Reduce(&max_value, &global_max, 1, MPI_TYPE_DIST_SORT_T, MPI_MAX, 0, MPI_COMM_WORLD);

		// Initialize splitters and assign global_max as last spliiter
		if (rank == 0) {
			dist_sort_size_t interval = ceil(float(data_size)/float(numSplitters));
			for (dist_sort_size_t i = interval-1; i < data_size; i += interval) {
					splitters[i/interval] = data[i];
			}
			splitters[numSplitters-1] = global_max;
		}
		// std::cerr << "pass here 1" << std::endl;

		// for (int i = 0; i < data_size; ++i) {
		// 		std::cerr << "data" << i << ":" << data[i]/DEBUG << ":rank:" << rank << std::endl;
		// }
		// int debug = 0;

		// Initialize upper/lowwer bound for each splitters
		dist_sort_t *lowerBound, *upperBound;
		if (rank == 0) {
				lowerBound = (dist_sort_t*)malloc(numSplitters*sizeof(dist_sort_t));
				upperBound = (dist_sort_t*)malloc(numSplitters*sizeof(dist_sort_t));
				for (dist_sort_size_t i = 0; i < numSplitters; ++i) {
						lowerBound[i] = 0;
						upperBound[i] = global_max;
				}
		}


		// Implement binary search to update them in each iterator
		while (true) {
				// std::cerr << "pass here 1" << std::endl;
				// if (rank == 0) {
				// 		for (int i = 0; i < numSplitters; ++i) {
				// 				std::cerr << "splitters" << i << ":" << splitters[i]/DEBUG << ":rank:" << rank << std::endl;
				// 		}
				// }

				std::cerr << "pass here 1.1:" << rank << std::endl;
				MPI_Bcast(splitters, numSplitters, MPI_TYPE_DIST_SORT_T, 0, MPI_COMM_WORLD);

				// std::cerr << "pass here 1.1" << std::endl;

				memset(counts, 0, numSplitters*sizeof(dist_sort_size_t));
				dist_sort_size_t j = 0;
				for (dist_sort_size_t i = 0; i < data_size; ++i) {
						while (data[i] > splitters[j]) {
								j += 1;
						}
						counts[j] += 1;
				}

				// std::cerr << "pass here 2" << std::endl;

				dist_sort_size_t *counts_buffer = NULL;
				if (rank == 0) {
						counts_buffer = (dist_sort_size_t*)malloc(nprocs*numSplitters*sizeof(dist_sort_size_t));
				}
				// std::cerr << "pass here 2.1" << std::endl;
				MPI_Gather(counts, numSplitters, MPI_TYPE_DIST_SORT_SIZE_T, counts_buffer, numSplitters, MPI_TYPE_DIST_SORT_SIZE_T, 0, MPI_COMM_WORLD);
				// std::cerr << "pass here 3" << std::endl;

				if (rank == 0) {
						memset(counts, 0, numSplitters*sizeof(dist_sort_size_t));
						// for (dist_sort_size_t i = 0; i < nprocs*numSplitters; ++i) {
						// 		counts[i%numSplitters] += counts_buffer[i];
						// 		std::cerr << "counts_buffer" << i << ":" << counts_buffer[i] << ":rank:" << rank << std::endl;
						// }

						// for (int i = 0; i < numSplitters; ++i) {
						// 		std::cerr << "counts" << i << ":" << counts[i] << ":rank:" << rank << std::endl;
						// }

						dist_sort_size_t prefix_counts[numSplitters];
						prefix_counts[0] = counts[0];
						for (dist_sort_size_t i = 1; i < numSplitters; ++i) {
								prefix_counts[i] += counts[i];
						}
						bool done = true;
						dist_sort_t new_splitters[numSplitters];
						// std::cerr << "pass here 3.1" << std::endl;
						for (dist_sort_size_t i = 0; i < numSplitters-1; ++i) {
								// std::cerr << "pass here 3.1" << std::endl;
								if (i+1 < ceil((float)(prefix_counts[i])/(float)global_N*numSplitters)) {
										// if (splitter_index[i] > 0) {
										// 	--(splitter_index[i]);
										// }
										// std::cerr << "pass here 3.1.1:" << i << ":" << (float)(prefix_counts[i])/(float)global_N*numSplitters << ":" << (float)(prefix_counts[i]) << ":" << (float)global_N << ":" << numSplitters << std::endl;
										dist_sort_size_t k = i;
										while (k > 0 && i < ceil((float)(prefix_counts[k])/(float)global_N*numSplitters)) {
												--k;
										}
										if (k != i)
												lowerBound[i] = std::max(lowerBound[i], splitters[k]);
										upperBound[i] = std::min(upperBound[i], splitters[i]);
										// std::cerr << "pass here 3.1.1.4:" << i << std::endl;
										done = false;
								} else if (i+1 > ceil((float)(prefix_counts[i])/(float)global_N*numSplitters)) {
										// if (splitter_index[i] < data_size-1) {
										// 	++(splitter_index[i]);
										// }
										// std::cerr << "pass here 3.1.2" << std::endl;
										dist_sort_size_t k = i;
										while (k < numSplitters-1 && i > ceil((float)(prefix_counts[k])/(float)global_N*numSplitters)) {
												++k;
										}
										lowerBound[i] = std::max(lowerBound[i], splitters[i]);
										if (k != i)
												upperBound[i] = std::min(upperBound[i], splitters[k]);
										done = false;
								} else {
										// std::cerr << "pass here 3.1.3" << std::endl;
										lowerBound[i] = splitters[i];
										upperBound[i] = splitters[i];
								}
								new_splitters[i] = (lowerBound[i] + upperBound[i]) / 2;
						}

						// for (int i = 0; i < numSplitters-1; ++i) {
						// 		std::cerr << "lowerBound" << i << ":" << lowerBound[i]/DEBUG << ":rank:" << rank << std::endl;
						// 		std::cerr << "upperBound" << i << ":" << upperBound[i]/DEBUG << ":rank:" << rank << std::endl;
						// }
						// std::cerr << "pass here 3.2" << std::endl;
						if (done) {
								// std::cerr << "pass here 3.1.4" << std::endl;
								free(counts_buffer);
								free(lowerBound);
								free(upperBound);
								break;
						}
						// std::cerr << "pass here 3.3" << std::endl;
						for (dist_sort_size_t i = 0; i < numSplitters-1; ++i) {
									splitters[i] = new_splitters[i];
						}
				}

				// debug++;
				// if (debug == 5) break;
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
