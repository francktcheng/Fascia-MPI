
#define CHUNK_SIZE 134217728

/*
MPI_Alltoallv_chunk(float* buffer, int* sendcounts, int* recvcounts)

  buf_v, comm->sendcounts, 
                  comm->sdispls, MPI_UINT64_T, 
                  comm->recvbuf_vert+sum_recv, comm->recvcounts, 
                  comm->rdispls, MPI_UINT64_T, MPI_COMM_WORLD)
{
  for (int32_t i = 0; i < nprocs; ++i)
    comm->recvcounts_temp[i] = 0;
  for (int32_t i = 0; i < nprocs; ++i)
    comm->sdispls_temp[i] -= comm->sendcounts_temp[i];

  MPI_Alltoall(comm->sendcounts_temp, 1, MPI_UINT64_T, 
               comm->recvcounts_temp, 1, MPI_UINT64_T, MPI_COMM_WORLD);

  comm->total_recv = 0;
  for (int i = 0; i < nprocs; ++i)
    comm->total_recv += comm->recvcounts_temp[i];

  comm->recvbuf_vert = (uint64_t*)malloc(comm->total_recv*sizeof(uint64_t));
  comm->recvbuf_data = (int32_t*)malloc(comm->total_recv*sizeof(uint32_t));
  if (comm->recvbuf_vert == NULL || comm->sendbuf_vert == NULL)
    throw_err("exchange_vert_data() unable to allocate comm buffers", procid);


  comm->global_queue_size = 0;
  uint64_t task_queue_size = comm->total_send;
  MPI_Allreduce(&task_queue_size, &comm->global_queue_size, 1, 
                MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);
  
  uint64_t num_comms = comm->global_queue_size / (uint64_t)MAX_SEND_SIZE + 1;
  uint64_t sum_recv = 0;
  uint64_t sum_send = 0;
  for (uint64_t c = 0; c < num_comms; ++c)
  {
    for (int32_t i = 0; i < nprocs; ++i)
    {
      uint64_t send_begin = (comm->sendcounts_temp[i] * c) / num_comms;
      uint64_t send_end = (comm->sendcounts_temp[i] * (c + 1)) / num_comms;
      if (c == (num_comms-1))
        send_end = comm->sendcounts_temp[i];
      comm->sendcounts[i] = (int32_t)(send_end - send_begin);
      assert(comm->sendcounts[i] >= 0);
    }

    MPI_Alltoall(comm->sendcounts, 1, MPI_INT32_T, 
                 comm->recvcounts, 1, MPI_INT32_T, MPI_COMM_WORLD);

    comm->sdispls[0] = 0;
    comm->sdispls_cpy[0] = 0;
    comm->rdispls[0] = 0;
    for (int32_t i = 1; i < nprocs; ++i)
    {
      comm->sdispls[i] = comm->sdispls[i-1] + comm->sendcounts[i-1];
      comm->rdispls[i] = comm->rdispls[i-1] + comm->recvcounts[i-1];
      comm->sdispls_cpy[i] = comm->sdispls[i];
    }

    int32_t cur_send = comm->sdispls[nprocs-1] + comm->sendcounts[nprocs-1];
    int32_t cur_recv = comm->rdispls[nprocs-1] + comm->recvcounts[nprocs-1];
    uint64_t* buf_v = (uint64_t*)malloc((uint64_t)(cur_send)*sizeof(uint64_t));
    int32_t* buf_d = (int32_t*)malloc((int32_t)(cur_send)*sizeof(int32_t));
    if (buf_v == NULL || buf_d == NULL)
      throw_err("exchange_verts(), unable to allocate comm buffers", procid);

    for (int32_t i = 0; i < nprocs; ++i)
    {
      uint64_t send_begin = (comm->sendcounts_temp[i] * c) / num_comms;
      uint64_t send_end = (comm->sendcounts_temp[i] * (c + 1)) / num_comms;
      if (c == (num_comms-1))
        send_end = comm->sendcounts_temp[i];

      for (uint64_t j = send_begin; j < send_end; ++j)
      {
        uint64_t vert = comm->sendbuf_vert[comm->sdispls_temp[i]+j];
        int32_t data = comm->sendbuf_data[comm->sdispls_temp[i]+j];
        buf_v[comm->sdispls_cpy[i]] = vert;
        buf_d[comm->sdispls_cpy[i]++] = data;
      }
    }

    MPI_Alltoallv(buf_v, comm->sendcounts, 
                  comm->sdispls, MPI_UINT64_T, 
                  comm->recvbuf_vert+sum_recv, comm->recvcounts, 
                  comm->rdispls, MPI_UINT64_T, MPI_COMM_WORLD);
    MPI_Alltoallv(buf_d, comm->sendcounts, 
                  comm->sdispls, MPI_INT32_T, 
                  comm->recvbuf_data+sum_recv, comm->recvcounts, 
                  comm->rdispls, MPI_INT32_T, MPI_COMM_WORLD);
    free(buf_v);
    free(buf_d);
    sum_recv += cur_recv;
    sum_send += cur_send;
  }
  free(comm->sendbuf_data);
  free(comm->sendbuf_vert);

  assert(sum_recv == comm->total_recv);
  assert(sum_send == comm->total_send);

  comm->global_queue_size = 0;
  task_queue_size = comm->total_recv + q->next_size;
  MPI_Allreduce(&task_queue_size, &comm->global_queue_size, 1, 
                MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);  

  q->next_size = 0;
  q->send_size = 0;
}
*/

void MPI_Bcast_chunk(int* arr, unsigned long length, int root, int rank)
{
#if DEBUG
  printf("%d MPI_B %d %lu\n", rank, root, length);
#endif

  if (length < CHUNK_SIZE)
  {
#if DEBUG
    printf("%d MPI_B no chunk %d %lu\n", rank, root, length);
#endif
    MPI_Bcast(arr, length, MPI_INT, root, MPI_COMM_WORLD);
  }
  else
  {
    unsigned long num_bcasts = length / CHUNK_SIZE;
    unsigned long cur_off = 0;
#if DEBUG
    printf("%d chunk numb %d %lu\n", rank, root, num_bcasts);
#endif
    for (int i = 0; i < num_bcasts; ++i)
    {      
#if DEBUG
      printf("%d doing b %d %lu\n", rank, root, cur_off);
#endif      
      MPI_Bcast(&arr[cur_off], CHUNK_SIZE, MPI_INT, root, MPI_COMM_WORLD);
      cur_off += CHUNK_SIZE;
    }
    unsigned long final_size = length - cur_off;
    assert(final_size > 0);
#if DEBUG
    printf("%d final b %d %lu %lu\n", rank, root, cur_off, final_size);
#endif  
    MPI_Bcast(&arr[cur_off], final_size, MPI_INT, root, MPI_COMM_WORLD);
  }
#if DEBUG
  printf("%d done\n", rank);
#endif
}

void MPI_Send_chunk(float* arr, unsigned long length, int to, int rank)
{
#if DEBUG
  printf("%d MPI_S f %d %lu\n", rank, to, length);
#endif
  if (length < CHUNK_SIZE) 
  {
#if DEBUG
    printf("%d MPI_S no chunk %d %lu\n", rank, to, length);
#endif
    MPI_Send(arr, length, MPI_FLOAT, to, 0, MPI_COMM_WORLD);
  }
  else
  {
    unsigned long num_bcasts = length / CHUNK_SIZE;
    unsigned long cur_off = 0;    
    for (int i = 0; i < num_bcasts; ++i)
    {            
#if DEBUG
      printf("%d MPI_S chunk %d %d %lu\n", rank, to, i, cur_off);
#endif
      MPI_Send(&arr[cur_off], CHUNK_SIZE, MPI_FLOAT, to, 0, MPI_COMM_WORLD);
      cur_off += CHUNK_SIZE;
    }
    unsigned long final_size = length - cur_off;
#if DEBUG
    printf("%d MPI_S chunk final %d %lu %lu\n", rank, to, final_size, cur_off);
#endif    
    assert(final_size > 0);
    MPI_Send(&arr[cur_off], final_size, MPI_FLOAT, to, 0, MPI_COMM_WORLD);
  }
#if DEBUG
  printf("%d done\n", rank);
#endif
}

void MPI_Send_chunk(short* arr, unsigned long length, int to, int rank)
{
#if DEBUG
  printf("%d MPI_S s %d %lu\n", rank, to, length);
#endif
  if (length < CHUNK_SIZE) 
  {
#if DEBUG
    printf("%d MPI_S no chunk %d %lu\n", rank, to, length);
#endif
    MPI_Send(arr, length, MPI_SHORT, to, 0, MPI_COMM_WORLD);
  }
  else
  {
    unsigned long num_bcasts = length / CHUNK_SIZE;
    unsigned long cur_off = 0;    
    for (int i = 0; i < num_bcasts; ++i)
    {            
#if DEBUG
      printf("%d MPI_S chunk %d %d %lu\n", rank, to, i, cur_off);
#endif
      MPI_Send(&arr[cur_off], CHUNK_SIZE, MPI_SHORT, to, 0, MPI_COMM_WORLD);
      cur_off += CHUNK_SIZE;
    }
    unsigned long final_size = length - cur_off;
#if DEBUG
    printf("%d MPI_S chunk final %d %lu %lu\n", rank, to, final_size, cur_off);
#endif    
    assert(final_size > 0);
    MPI_Send(&arr[cur_off], final_size, MPI_SHORT, to, 0, MPI_COMM_WORLD);
  }
}


void MPI_Send_chunk(int* arr, unsigned long length, int to, int rank)
{
#if DEBUG
  printf("%d MPI_S i %d %lu\n", rank, to, length);
#endif
  if (length < CHUNK_SIZE) 
  {
#if DEBUG
    printf("%d MPI_S no chunk %d %lu\n", rank, to, length);
#endif
    MPI_Send(arr, length, MPI_INT, to, 0, MPI_COMM_WORLD);
  }
  else
  {
    unsigned long num_bcasts = length / CHUNK_SIZE;
    unsigned long cur_off = 0;    
    for (int i = 0; i < num_bcasts; ++i)
    {            
#if DEBUG
      printf("%d MPI_S chunk %d %d %lu\n", rank, to, i, cur_off);
#endif
      MPI_Send(&arr[cur_off], CHUNK_SIZE, MPI_INT, to, 0, MPI_COMM_WORLD);
      cur_off += CHUNK_SIZE;
    }
    unsigned long final_size = length - cur_off;
#if DEBUG
    printf("%d MPI_S chunk final %d %lu %lu\n", rank, to, final_size, cur_off);
#endif    
    assert(final_size > 0);
    MPI_Send(&arr[cur_off], final_size, MPI_INT, to, 0, MPI_COMM_WORLD);
  }
#if DEBUG
  printf("%d done\n", rank);
#endif
}


void MPI_Send_chunk(unsigned long* arr, unsigned long length, int to, int rank)
{
#if DEBUG
  printf("%d MPI_S ul %d %lu\n", rank, to, length);
#endif
  if (length < CHUNK_SIZE) 
  {
#if DEBUG
    printf("%d MPI_S no chunk %d %lu\n", rank, to, length);
#endif
    MPI_Send(arr, length, MPI_UNSIGNED_LONG, to, 0, MPI_COMM_WORLD);
  }
  else
  {
    unsigned long num_bcasts = length / CHUNK_SIZE;
    unsigned long cur_off = 0;    
    for (int i = 0; i < num_bcasts; ++i)
    {            
#if DEBUG
      printf("%d MPI_S chunk %d %d %lu\n", rank, to, i, cur_off);
#endif
      MPI_Send(&arr[cur_off], CHUNK_SIZE, MPI_UNSIGNED_LONG, to, 0, MPI_COMM_WORLD);
      cur_off += CHUNK_SIZE;
    }
    unsigned long final_size = length - cur_off;
#if DEBUG
    printf("%d MPI_S chunk final %d %lu %lu\n", rank, to, final_size, cur_off);
#endif    
    assert(final_size > 0);
    MPI_Send(&arr[cur_off], final_size, MPI_UNSIGNED_LONG, to, 0, MPI_COMM_WORLD);
  }
#if DEBUG
  printf("%d done\n", rank);
#endif
}


void MPI_Recv_chunk(float* arr, unsigned long length, int from, int rank)
{
#if DEBUG
  printf("%d MPI_R f %d %lu\n", rank, from, length);
#endif
  if (length < CHUNK_SIZE) 
  {    
#if DEBUG
  printf("%d MPI_R no chunk %d %lu\n", rank, from, length);
#endif
    MPI_Recv(arr, length, MPI_FLOAT, from, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }
  else
  {
    unsigned long num_bcasts = length / CHUNK_SIZE;
    unsigned long cur_off = 0;
    for (int i = 0; i < num_bcasts; ++i)
    {    
#if DEBUG
      printf("%d MPI_R chunk %d %d %lu\n", rank, from, i, cur_off);
#endif       
      MPI_Recv(&arr[cur_off], CHUNK_SIZE, MPI_FLOAT, from, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      cur_off += CHUNK_SIZE;
    }
    unsigned long final_size = length - cur_off;   
#if DEBUG
      printf("%d MPI_R chunk final %d %lu %lu\n", rank, from, final_size, cur_off);
#endif  
    assert(final_size > 0);
    MPI_Recv(&arr[cur_off], final_size, MPI_FLOAT, from, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }
#if DEBUG
  printf("%d done\n",rank);
#endif
}


void MPI_Recv_chunk(short* arr, unsigned long length, int from, int rank)
{
#if DEBUG
  printf("%d MPI_R s %d %lu\n", rank, from, length);
#endif
  if (length < CHUNK_SIZE) 
  {    
#if DEBUG
  printf("%d MPI_R no chunk %d %lu\n", rank, from, length);
#endif
    MPI_Recv(arr, length, MPI_SHORT, from, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }
  else
  {
    unsigned long num_bcasts = length / CHUNK_SIZE;
    unsigned long cur_off = 0;
    for (int i = 0; i < num_bcasts; ++i)
    {    
#if DEBUG
      printf("%d MPI_R chunk %d %d %lu\n", rank, from, i, cur_off);
#endif       
      MPI_Recv(&arr[cur_off], CHUNK_SIZE, MPI_SHORT, from, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      cur_off += CHUNK_SIZE;
    }
    unsigned long final_size = length - cur_off;   
#if DEBUG
      printf("%d MPI_R chunk final %d %lu %lu\n", rank, from, final_size, cur_off);
#endif  
    assert(final_size > 0);
    MPI_Recv(&arr[cur_off], final_size, MPI_SHORT, from, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }
#if DEBUG
  printf("%d done\n",rank);
#endif
}


void MPI_Recv_chunk(int* arr, unsigned long length, int from, int rank)
{
#if DEBUG
  printf("%d MPI_R i %d %lu\n", rank, from, length);
#endif
  if (length < CHUNK_SIZE) 
  {    
#if DEBUG
  printf("%d MPI_R no chunk %d %lu\n", rank, from, length);
#endif
    MPI_Recv(arr, length, MPI_INT, from, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }
  else
  {
    unsigned long num_bcasts = length / CHUNK_SIZE;
    unsigned long cur_off = 0;
    for (int i = 0; i < num_bcasts; ++i)
    {    
#if DEBUG
      printf("%d MPI_R chunk %d %d %lu\n", rank, from, i, cur_off);
#endif       
      MPI_Recv(&arr[cur_off], CHUNK_SIZE, MPI_INT, from, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      cur_off += CHUNK_SIZE;
    }
    unsigned long final_size = length - cur_off;   
#if DEBUG
      printf("%d MPI_R chunk final %d %lu %lu\n", rank, from, final_size, cur_off);
#endif  
    assert(final_size > 0);
    MPI_Recv(&arr[cur_off], final_size, MPI_INT, from, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }
#if DEBUG
  printf("%d done\n",rank);
#endif
}

void MPI_Recv_chunk(unsigned long* arr, unsigned long length, int from, int rank)
{
#if DEBUG
  printf("%d MPI_R ul %d %lu\n", rank, from, length);
#endif
  if (length < CHUNK_SIZE) 
  {    
#if DEBUG
  printf("%d MPI_R no chunk %d %lu\n", rank, from, length);
#endif
    MPI_Recv(arr, length, MPI_UNSIGNED_LONG, from, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }
  else
  {
    unsigned long num_bcasts = length / CHUNK_SIZE;
    unsigned long cur_off = 0;
    for (int i = 0; i < num_bcasts; ++i)
    {    
#if DEBUG
      printf("%d MPI_R chunk %d %d %lu\n", rank, from, i, cur_off);
#endif       
      MPI_Recv(&arr[cur_off], CHUNK_SIZE, MPI_UNSIGNED_LONG, from, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      cur_off += CHUNK_SIZE;
    }
    unsigned long final_size = length - cur_off;   
#if DEBUG
      printf("%d MPI_R chunk final %d %lu %lu\n", rank, from, final_size, cur_off);
#endif  
    assert(final_size > 0);
    MPI_Recv(&arr[cur_off], final_size, MPI_UNSIGNED_LONG, from, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }
#if DEBUG
  printf("%d done\n",rank);
#endif
}
