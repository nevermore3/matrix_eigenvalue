/*************************************************
根据特征值分布ceig分割区间
*************************************************/
#include "type.h"
void DV_GetInterval(double *ceig, Info info)
{
    double *interval = info.interval;
    int myid = info.myid;
    int np = info.ninterval;
    int *DVtag = info.DVtag;
    if (myid < np)
    {
        int ninterval = info.ninterval;
        int i, id;
        id = (int)((myid + 1) * info.N * info.propotion * info.scal[0] / ninterval);
        interval[1] = ceig[id - 1];
        if (np > 1)
        {

            int mpitag = 1991;
            if (myid == 0)
            {
                interval[0] = info.mineig;

                if (DVtag[myid + 1] == 0)
                    MPI_Send(interval + 1, 1, MPI_DOUBLE, myid + 1, mpitag, MPI_COMM_WORLD);

            }
            else if (myid == np - 1)
            {
                MPI_Status status;
                if (DVtag[myid - 1] == 0)
                    MPI_Recv(interval, 1, MPI_DOUBLE, myid - 1, mpitag, MPI_COMM_WORLD, &status);

                while (interval[1] < interval[0] && id < info.history[1] )
                {
                    id = id + 1;
                    interval[1] = ceig[id - 1];
                }
                if (interval[1] < interval[0])
                    printf("***************************  error  !!! in DV_GetInterval\n");
            }
            else
            {
                MPI_Status status;
                if (DVtag[myid - 1] == 0)
                    MPI_Recv(interval, 1, MPI_DOUBLE, myid - 1, mpitag, MPI_COMM_WORLD, &status);

                while (interval[1] < interval[0] && id < info.history[1] )
                {
                    id = id + 1;
                    interval[1] = ceig[id - 1];
                }
                if (interval[1] < interval[0])
                    printf("***************************  error  !!! in DV_GetInterval\n");

                if (DVtag[myid + 1] == 0)
                    MPI_Send(interval + 1, 1, MPI_DOUBLE, myid + 1, mpitag, MPI_COMM_WORLD);
            }

        }
        printf("myid %d interval: \t %lf  %d  scal: %lf\n", info.myid, interval[1], id - 1, info.scal[0]);

    }

    MPI_Bcast(interval, 2, MPI_DOUBLE, 0, info.intervalcomm);
}
