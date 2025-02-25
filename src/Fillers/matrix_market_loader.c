#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/mman.h>
#include <unistd.h>
#include <fcntl.h>

#include "matrix_market_loader.h"

// Credit : https://github.com/Rhythmicc/FastMatrixLoader/
void mm_load_coo(const ascii *filepath, coo_matrix_t *matrix)
{
    ascii *data;
    usz len = 0;

    usz rows, cols, nnz;

    i32 fd = open(filepath, O_RDONLY);
    if(fd < 0)
        return; 
    len = lseek(fd, 0, SEEK_END);
    
    ascii *mbuf = (ascii *)mmap(NULL, len, PROT_READ, MAP_PRIVATE, fd, 0);
    data = (ascii *)malloc(len);
    memcpy(data, mbuf, len);
    close(fd);
    munmap(mbuf, len);

    ascii *line = strtok(data, "\n");
    while(line && (line[0] == '%'))
        line = strtok(NULL, "\n"); 

    sscanf(line, "%zu%zu%zu", &rows, &cols, &nnz);
    
    // Symmetric
    len = len * 2;

    i32 *ia = (i32 *)malloc(len * sizeof(i32));
    i32 *ja = (i32 *)malloc(len * sizeof(i32)); 
    f64 *val = (f64 *)malloc(len * sizeof(f64));
    
    line = strtok(NULL, "\n");

    for(i32 idx = 0; idx < nnz; idx++)
    {
        sscanf(line, "%d%d%lg", ia + idx, ja + idx, val + idx);
        --ia[idx], --ja[idx];

        if(ia[idx] != ja[idx])
        {
            ++idx, ++nnz;
            usz tmp = ia[idx -1];
            ia[idx] = ja[idx -1];
            ja[idx] = tmp;

            val[idx] = val[idx -1];
        }

        if((line = strtok(NULL, "\n")) == NULL && (idx +1 < nnz))
        {
            free(data);
            free(ia);
            free(ja);
        }
    }
    
    allocate_COO(rows, cols, nnz, matrix);
    matrix->row_index = (usz*)ia;
    matrix->col_index = (usz *)ja;
    val = (f64 *)realloc(val, sizeof(f64) * nnz);
    matrix->data = val;
}
