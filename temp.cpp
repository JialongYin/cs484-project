void my_sort(int N, item *myItems, int *nOut, item **myResult)
{
    int rank, nprocs;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int local_count[nprocs] = {0};
    for (int i = 0; i < N; ++i)
    {
        int index = myItems[i].key;
        ++local_count[index];
    }
    for (int i = 0; i < nprocs; ++i) {
        MPI_Reduce(&local_count[i], nOut, 1, MPI_INT, MPI_SUM, i, MPI_COMM_WORLD);
    }
    // std::cout << (*nOut) << std::endl;
    int global_count[nprocs] = {0};
    for (int i = 0; i < nprocs; ++i)
    {
        MPI_Exscan(&local_count[i], &global_count[i], 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        // std::cout << global_count[i] << std::endl;
    }
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
