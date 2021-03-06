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


		dist_sort_size_t global_counts = 0;
    MPI_Exscan(&myDataCount, &global_counts, 1, MPI_TYPE_DIST_SORT_SIZE_T, MPI_SUM, MPI_COMM_WORLD);
    MPI_Win win;
    MPI_Win_create(*rebalancedData, (*rCount) * sizeof(dist_sort_t), sizeof(dist_sort_t), MPI_INFO_NULL, MPI_COMM_WORLD, &win);
    MPI_Win_fence(MPI_MODE_NOPRECEDE, win); //fence - there are no epochs before this
    dist_sort_size_t i = 0;
		dist_sort_size_t max_size = ceil((float)global_N/(float)nprocs);
		// TODO: add #pragma parallel for
    while (i < myDataCount)
    {
        int target_rank = int((float)(global_counts+i) / (float)max_size);
        dist_sort_size_t displacement = (global_counts+i) % max_size;
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
		uint64_t DEBUG = 10000000000000000;
		int rank, nprocs;
		MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);

		// Find global number of keys in each process
		dist_sort_size_t global_N;
		MPI_Allreduce(&data_size, &global_N, 1, MPI_TYPE_DIST_SORT_SIZE_T, MPI_SUM, MPI_COMM_WORLD);
		// Find global maximum in rank0 process
		// TODO: add #pragma parallel for
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
			dist_sort_size_t interval = floor(float(data_size)/float(numSplitters));
			// TODO: add #pragma parallel for
			for (dist_sort_size_t i = interval-1; i < data_size; i += interval) {
					splitters[i/interval] = data[i];
			}
			splitters[numSplitters-1] = global_max;
		}

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
				// if (rank == 0) {
				// 		for (int i = 0; i < numSplitters; ++i) {
				// 				std::cerr << "splitters" << i << ":" << splitters[i]/DEBUG << ":rank:" << rank << std::endl;
				// 		}
				// }
				MPI_Bcast(splitters, numSplitters, MPI_TYPE_DIST_SORT_T, 0, MPI_COMM_WORLD);
				memset(counts, 0, numSplitters*sizeof(dist_sort_size_t));
				dist_sort_size_t j = 0;
				// TODO: add #pragma parallel for
				for (dist_sort_size_t i = 0; i < data_size; ++i) {
						while (data[i] > splitters[j]) {
								j += 1;
						}
						counts[j] += 1;
				}

				dist_sort_size_t *counts_buffer = NULL;
				if (rank == 0) {
						counts_buffer = (dist_sort_size_t*)malloc(nprocs*numSplitters*sizeof(dist_sort_size_t));
				}
				MPI_Gather(counts, numSplitters, MPI_TYPE_DIST_SORT_SIZE_T, counts_buffer, numSplitters, MPI_TYPE_DIST_SORT_SIZE_T, 0, MPI_COMM_WORLD);

				bool done;
				if (rank == 0) {
						memset(counts, 0, numSplitters*sizeof(dist_sort_size_t));
						for (dist_sort_size_t i = 0; i < nprocs*numSplitters; ++i) {
								counts[i%numSplitters] += counts_buffer[i];
						}

						// for (int i = 0; i < numSplitters; ++i) {
						// 		std::cerr << "counts" << i << ":" << counts[i] << ":rank:" << rank << std::endl;
						// }

						dist_sort_size_t prefix_counts[numSplitters];
						prefix_counts[0] = counts[0];
						for (dist_sort_size_t i = 1; i < numSplitters; ++i) {
								prefix_counts[i] = prefix_counts[i-1] + counts[i];
						}
						dist_sort_t new_splitters[numSplitters];
						done = true;
						for (dist_sort_size_t i = 0; i < numSplitters-1; ++i) {
								if ((i+1) * global_N + numSplitters <= prefix_counts[i] * numSplitters) { // (i+1) * global_N < prefix_counts[i] * numSplitters
										// std::cerr << "pass here 3.1.1:" << (i+1) << ":" << global_N << ":" << prefix_counts[i] << ":" << numSplitters << std::endl;
										// std::cerr << "pass here 3.1.1" << std::endl;
										dist_sort_size_t k = i;
										while (k > 0 && (i+1) * global_N + numSplitters <= prefix_counts[k] * numSplitters) {
												--k;
										}
										if (k != i && (i+1) * global_N + numSplitters > prefix_counts[k] * numSplitters)
												lowerBound[i] = std::max(lowerBound[i], splitters[k]);
										upperBound[i] = std::min(upperBound[i], splitters[i]);
										done = false;
								} else if ((i+1) * global_N >= prefix_counts[i] * numSplitters + numSplitters) { // (i+1) * global_N > prefix_counts[i] * numSplitters
										// std::cerr << "pass here 3.1.2:" << (i+1) << ":" << global_N << ":" << prefix_counts[i] << ":" << numSplitters << std::endl;
										// std::cerr << "pass here 3.1.2" << std::endl;
										dist_sort_size_t k = i;
										while (k < numSplitters-1 && (i+1) * global_N >= prefix_counts[k] * numSplitters + numSplitters) {
												++k;
										}
										lowerBound[i] = std::max(lowerBound[i], splitters[i]);
										if (k != i && (i+1) * global_N < prefix_counts[k] * numSplitters + numSplitters)
												upperBound[i] = std::min(upperBound[i], splitters[k]);
										done = false;
								} else {
										// std::cerr << "pass here 3.1.3:" << (i+1) << ":" << global_N << ":" << prefix_counts[i] << ":" << numSplitters << std::endl;
										// std::cerr << "pass here 3.1.3" << std::endl;
										lowerBound[i] = splitters[i];
										upperBound[i] = splitters[i];
								}
								new_splitters[i] = lowerBound[i] + ((upperBound[i] - lowerBound[i]) / 2);
								// std::cerr << "lowerBound" << i << ":" << lowerBound[i]/DEBUG << ":rank:" << rank << std::endl;
								// std::cerr << "upperBound" << i << ":" << upperBound[i]/DEBUG << ":rank:" << rank << std::endl;
						}
						for (dist_sort_size_t i = 0; i < numSplitters-1; ++i) {
									splitters[i] = new_splitters[i];
						}
				}
				MPI_Bcast(&done, 1, MPI_C_BOOL, 0, MPI_COMM_WORLD);
				if (done) {
						// std::cerr << "pass here 3.1.5:" << rank << std::endl;
						if (rank == 0) {
							free(counts_buffer);
							free(lowerBound);
							free(upperBound);
						}
						break;
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
		uint64_t DEBUG = 10000000000000000;
		int rank, nprocs;
		MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);

		// if (rank == 0) {
		// 		for (int i = 0; i < numSplitters; ++i) {
		// 				std::cerr << "splitters" << i << ":" << splitters[i]/DEBUG << ":rank:" << rank << std::endl;
		// 		}
		// 		for (int i = 0; i < sDataCount; ++i) {
		// 				std::cerr << "counts" << i << ":" << counts[i] << ":rank:" << rank << std::endl;
		// 		}
		// }


		// std::cerr << "sDataCount:" << sDataCount << ":rank:" << rank << std::endl;
		// for (dist_sort_size_t i = 0; i < sDataCount; ++i) {
		// 		std::cerr << "data" << i << ":" << sendData[i]/DEBUG << ":rank:" << rank << std::endl;
		// }

		dist_sort_t *splittersBuffer = (dist_sort_t*)malloc(numSplitters*sizeof(dist_sort_t));
		dist_sort_t *countsBuffer = (dist_sort_t*)malloc(numSplitters*sizeof(dist_sort_t));
		if (rank == 0) {
				for (int i = 0; i < numSplitters; ++i) {
						splittersBuffer[i] = splitters[i];
						countsBuffer[i] = counts[i];
				}
		}
		MPI_Bcast(splittersBuffer, numSplitters, MPI_TYPE_DIST_SORT_T, 0, MPI_COMM_WORLD);
		MPI_Bcast(countsBuffer, numSplitters, MPI_TYPE_DIST_SORT_SIZE_T, 0, MPI_COMM_WORLD);
		// std::cerr << "pass here 3.1:" << rank << std::endl;

		*recvData = (dist_sort_t*)malloc(countsBuffer[rank]*sizeof(dist_sort_t));
		*rDataCount = countsBuffer[rank];
		// for (int i = 0; i < numSplitters; ++i) {
		// 		std::cerr << "countsBuffer" << i << ":" << countsBuffer[i] << ":rank:" << rank << std::endl;
		// }

		dist_sort_size_t local_counts[numSplitters] = {0};
		dist_sort_size_t j = 0;
		// TODO: add #pragma parallel for
		for (dist_sort_size_t i = 0; i < sDataCount; ++i) {
				while (sendData[i] > splittersBuffer[j]) {
						++j;
				}
				local_counts[j] += 1;
		}
		dist_sort_size_t global_counts[numSplitters] = {0};
		for (dist_sort_size_t i = 0; i < numSplitters; ++i) {
				MPI_Exscan(&local_counts[i], &global_counts[i], 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
		}
		// std::cerr << "pass here 3.2:" << rank << std::endl;

		MPI_Win win;
		MPI_Win_create(*recvData, countsBuffer[rank]*sizeof(dist_sort_t), sizeof(dist_sort_t), MPI_INFO_NULL, MPI_COMM_WORLD, &win);
		MPI_Win_fence(MPI_MODE_NOPRECEDE, win);
		dist_sort_size_t i = 0;
		dist_sort_size_t send_counts[nprocs] = {0};
		j = 0;
		// std::cerr << "pass here 3.3:" << rank << std::endl;

		// for (int i = 0; i < numSplitters; ++i) {
		// 		std::cerr << "local_counts" << i << ":" << local_counts[i] << ":rank:" << rank << std::endl;
		// }

		dist_sort_size_t pre_local_counts = local_counts[0];
		// TODO: add #pragma parallel for
		while (i < sDataCount)
		{
				// std::cerr << "i" << i << ";local_counts" << j << ":" << local_counts[j] << ":rank:" << rank << std::endl;
				if (i >= pre_local_counts) {
						++j;
						pre_local_counts += local_counts[j];
				}
				int target_rank = j;
				dist_sort_size_t displacement = global_counts[target_rank] + send_counts[target_rank];
				if (target_rank != rank) {
						// std::cerr << "target_rank:" << target_rank << ";rank:" << rank << std::endl;
						MPI_Put(&(sendData[i]), 1, MPI_TYPE_DIST_SORT_T, target_rank, displacement, 1, MPI_TYPE_DIST_SORT_T, win);
						// std::cerr << "attention:" << rank << std::endl;
				} else {
						(*recvData)[displacement] = sendData[i];
				}
				send_counts[target_rank] += 1;
				++i;
				// std::cerr << "pass here 3.4:" << rank << std::endl;
		}
		MPI_Win_fence(0, win);
		MPI_Win_fence(MPI_MODE_NOSUCCEED, win);
		free(countsBuffer);
		free(splittersBuffer);
		// std::cerr << "pass here 3.5:" << rank << std::endl;

		// for (dist_sort_size_t i = 0; i < *rDataCount; ++i) {
		// 		std::cerr << "data" << i << ":" << (*recvData)[i]/DEBUG << ":rank:" << rank << std::endl;
		// }
		// std::cerr << "rDataCount:" << *rDataCount << ":rank:" << rank << std::endl;
}

void sort(dist_sort_t *data, dist_sort_size_t size) {
	// You are welcome to use this sort function.
	// If you wish to implement your own, you can do that too.
	// Don't use bubblesort.
	std::sort(data,data+size);
}
