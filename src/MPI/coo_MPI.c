#include "coo_MPI.h"
#include <string.h>

void distribute_coo(coo_matrix_t *global_matrix, coo_matrix_t *local_matrix)
{
    i32 size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    if(rank == 0)
    {
        i32 rows_per_proc   = global_matrix->dim_x / size;
        i32 start_row       = rank * rows_per_proc;
        i32 end_row         = (rank == size -1) ? global_matrix->dim_x : start_row + rows_per_proc;
        usz nnz_global = global_matrix->nnz;

        usz send_count[size];
        memset(send_count, 0, size * sizeof(usz));
        for(usz i = 0; i < nnz_global; i++)
        {
            i32 owner = (global_matrix->row_index[i] * size) / global_matrix->dim_x;
            send_count[owner]++;
        }
        
        usz displ[size];
        memset(displ, 0, size * sizeof(usz));
        for(usz i = 1; i < size; i++)
        {
            displ[i] = displ[i-1] + send_count[i-1];
        }
    
        usz offsets[size];
        memset(offsets, 0, size * sizeof(usz));
        usz row_sorted[nnz_global];
        usz col_sorted[nnz_global];
        f64 val_sorted[nnz_global];
        for(usz i = 0; i < nnz_global; i++)
        {
            i32 owner = (global_matrix->row_index[i] * size) / global_matrix->dim_x;
            i32 pos = displ[owner] + offsets[owner]++;
            row_sorted[pos] = global_matrix->row_index[i];
            col_sorted[pos] = global_matrix->col_index[i];
            val_sorted[pos] = global_matrix->data[i];
        }

        for(usz i = 1; i < size; i++) 
        {
            MPI_Send(&(global_matrix->dim_x), 1             , MPI_UNSIGNED_LONG, i, 0, MPI_COMM_WORLD);
            MPI_Send(&(global_matrix->dim_y), 1             , MPI_UNSIGNED_LONG, i, 1, MPI_COMM_WORLD);
            MPI_Send(&send_count[i]         , 1             , MPI_UNSIGNED_LONG, i, 2, MPI_COMM_WORLD);
            MPI_Send(row_sorted + displ[i]  , send_count[i] , MPI_UNSIGNED_LONG, i, 3, MPI_COMM_WORLD);
            MPI_Send(col_sorted + displ[i]  , send_count[i] , MPI_UNSIGNED_LONG, i, 4, MPI_COMM_WORLD);
            MPI_Send(val_sorted + displ[i]  , send_count[i] , MPI_DOUBLE       , i, 5, MPI_COMM_WORLD);
        }
        
        allocate_COO(global_matrix->dim_x,global_matrix->dim_y ,send_count[0], local_matrix);

        for(usz i = 0; i < send_count[0]; i++)
        {
            local_matrix->data[i]       = val_sorted[i];
            local_matrix->row_index[i]  = row_sorted[i];
            local_matrix->col_index[i]  = col_sorted[i];
        }
    }
    else
    {  
        usz nnz_local   = 0;
        usz local_dim_x = 0;
        usz local_dim_y = 0;
        MPI_Recv(&local_dim_x   , 1, MPI_UNSIGNED_LONG, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&local_dim_y   , 1, MPI_UNSIGNED_LONG, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&nnz_local     , 1, MPI_UNSIGNED_LONG, 0, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        
        allocate_COO(local_dim_x, local_dim_y, nnz_local, local_matrix);
         
        MPI_Recv(local_matrix->row_index, nnz_local, MPI_UNSIGNED_LONG  , 0, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(local_matrix->col_index, nnz_local, MPI_UNSIGNED_LONG  , 0, 4, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(local_matrix->data     , nnz_local, MPI_DOUBLE         , 0, 5, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
}
