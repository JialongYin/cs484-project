void rebalance(const dist_sort_t *data, const dist_sort_size_t myDataCount, dist_sort_t **rebalancedData, dist_sort_size_t *rCount) {
/*
	See the header file ('solution.hpp') for Doxygen docstrings explaining this function and its parameters.
*/
		int rank, nprocs;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

		dist_sort_size_t global_N;
		MPI_Allreduce(&myDataCount, &global_N, 1, MPI_TYPE_DIST_SORT_SIZE_T, MPI_SUM, MPI_COMM_WORLD);
		*rCount = ceil(global_N / nprocs);
		*rebalancedData = (dist_sort_t*)malloc((*rCount)*sizeof(dist_sort_t));

		dist_sort_size_t global_count = 0;
    MPI_Exscan(&myDataCount, &global_count, 1, MPI_TYPE_DIST_SORT_SIZE_T, MPI_SUM, MPI_COMM_WORLD);

    MPI_Win win;
    *myResult = (item *)malloc((*nOut) * sizeof(item));
    MPI_Win_create(*myResult, (*nOut) * sizeof(item), sizeof(item), MPI_INFO_NULL, MPI_COMM_WORLD, &win);
    MPI_Win_fence(MPI_MODE_NOPRECEDE, win); //fence - there are no epochs before this
    int i = 0;
    int send_count[nprocs] = {0};
    while (i < N)
    {
        int target_rank = myItems[i].key;
        int displacement = global_count[target_rank] + send_count[target_rank];
        if (target_rank != rank) {
            MPI_Put(&(myItems[i]), 2, MPI_INT, target_rank, displacement, 2, MPI_INT, win);
        } else {
            (*myResult)[displacement].key = myItems[i].key;
            (*myResult)[displacement].value = myItems[i].value;
        }
        send_count[target_rank] += 1;
        ++i;

        // std::cout << myItems[i].key << "," << myItems[i].value << std::endl << std::flush;
        // std::cout << "rank:" << rank << ";item:" << myItems[i].key << "," << myItems[i].value << ";target:" << target_rank << ";displacement:" << displacement << std::endl
        //           << std::flush;
    }
    MPI_Win_fence(0, win);
    MPI_Win_fence(MPI_MODE_NOSUCCEED, win);

    // for (int i = 0; i < (*nOut); ++i) {
    //     std::cout << rank << ":" << (*myResult)[i].key << std::endl;
    // }

}
