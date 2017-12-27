#include <stdio.h>
const char *dgemm_desc = "My awesome dgemm.";

#define MIN(A, B) A <= B ? A : B

#ifndef BLOCK_SIZE_JJJ
#define BLOCK_SIZE_JJJ ((int)16)
#endif

#ifndef BLOCK_SIZE_JJ
#define BLOCK_SIZE_JJ ((int)4)
#endif

#ifndef BLOCK_SIZE_KK
#define BLOCK_SIZE_KK ((int)16)
#endif

void square_dgemm(const int M, const double *A, const double *B, double *C)
{
    int jjj, jj, kk;
    int limitjjj, limitjj, limitkk;
    int j, k, i;
    for (jjj = 0; jjj < M; jjj += BLOCK_SIZE_JJJ)
    {
        for (kk = 0; kk < M; kk += BLOCK_SIZE_KK)
        {
            limitjjj = MIN((jjj + BLOCK_SIZE_JJJ), M);
            limitkk = MIN((kk + BLOCK_SIZE_KK), M);

            for (jj = jjj; jj < limitjjj; jj += BLOCK_SIZE_JJ)
            {
                limitjj = MIN(jj + BLOCK_SIZE_JJ, M);
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
}
