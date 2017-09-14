#include "type.h";
#include <math.h>


void orth(double *Q_loc, double *Y_loc, Info info)
{
    int i, j, k;

    int myN_loc = info.N_loc[info.intervalid];
    int M0 = info.M0;

    for (i = 0; i < myN_loc * M0; i++)
        Q_loc[i] = 0;

    double sum, sum_loc;
    double *mysum, *mysum_loc;
    mysum = calloc(M0, sizeof(double));
    mysum_loc = calloc(M0, sizeof(double));


    for (i = 0; i < M0; i++)
    {
        for (j = 0; j < myN_loc; j++)
            Q_loc[i * myN_loc + j] = Y_loc[i * myN_loc + j] - Q_loc[i * myN_loc + j] ;

        sum_loc = 0;
        for (j = 0; j < myN_loc; j++)
            sum_loc = sum_loc + Q_loc[i * myN_loc + j] * Q_loc[i * myN_loc + j];

        MPI_Allreduce(&sum_loc, &sum, 1, MPI_DOUBLE, MPI_SUM, info.intervalcomm);

        sum = sqrt(sum);

        for (j = 0; j < myN_loc; j++)
            Q_loc[i * myN_loc + j] = Q_loc[i * myN_loc + j] / sum;

        for (j = i + 1; j < M0; j++)
        {
            mysum_loc[j] = 0;
            mysum[j] = 0;
            for (k = 0; k < myN_loc; k++)
                mysum_loc[j] = mysum_loc[j] + Q_loc[i * myN_loc + k] * Y_loc[j * myN_loc + k];
        }

        MPI_Allreduce(mysum_loc + i + 1, mysum + i + 1, M0 - i - 1, MPI_DOUBLE, MPI_SUM, info.intervalcomm);

        for (j = i + 1; j < M0; j++)
        {
            for (k = 0; k < myN_loc; k++)
                Q_loc[j * myN_loc + k] = Q_loc[j * myN_loc + k] + mysum[j] * Q_loc[i * myN_loc + k];
        }
    }

    free(mysum);
    free(mysum_loc);


}

