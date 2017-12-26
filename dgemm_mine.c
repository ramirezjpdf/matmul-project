#include <stdio.h>
const char *dgemm_desc = "My awesome dgemm.";

#define MIN(A, B) A <= B ? A : B

#ifndef BLOCK_SIZE_MINE
#define BLOCK_SIZE_MINE ((int)64)
#endif

void square_dgemm(const int M, const double *A, const double *B, double *C)
{
    int jj, kk;
    int limitjj, limitkk;
    int j, k, i;
    for (jj = 0; jj < M; jj += BLOCK_SIZE_MINE)
    {
        for (kk = 0; kk < M; kk += BLOCK_SIZE_MINE)
        {
            limitjj = MIN((jj + BLOCK_SIZE_MINE), M);
            limitkk = MIN((kk + BLOCK_SIZE_MINE), M);

            for (j = jj; j < limitjj; ++j)
            {
                for (k = kk; k < limitkk; ++k)
                {
                    for (i = 0; i < M; ++i)
                    {
                        C[j * M + i] += A[k * M + i] * B[j * M + k];
                    }
                }
            }
        }
    }
}
