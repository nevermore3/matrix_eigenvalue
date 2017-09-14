/****************************************************
获取谱分割的参数，判断终止条件
输入:
      int *neig
      int neigtarget
      int ninterval
输出:
      int *tag
*************************************************/
#include "type.h"
void DV_GetTag(int *tag, int neig_loc, int times, Info info)
{
    int *neig = info.neig;
    int *sneig = info.sneig;
    int *DVtag = info.DVtag;
    int myid = info.myid;
    int np = info.ninterval;
    int M0 = info.M0;

    int i, j;

    int mpitag = 1992;
    if (info.intervalid == 0)//跟进程修正结果
    {
        MPI_Status status;
        for (i = 0; i < np; i++)
        {
            if (i < myid && DVtag[i] == 0)
            {
                MPI_Recv(neig + i, 1, MPI_INT, i, mpitag, MPI_COMM_WORLD, &status);
                MPI_Send(&neig_loc, 1, MPI_INT, i, mpitag, MPI_COMM_WORLD);
            }
            if (i > myid && DVtag[i] == 0)
            {
                MPI_Send(&neig_loc, 1, MPI_INT, i, mpitag, MPI_COMM_WORLD);
                MPI_Recv(neig + i, 1, MPI_INT, i, mpitag, MPI_COMM_WORLD, &status);
            }
        }

        neig[myid] = neig_loc;

        int maxneig;
        double dist;
        sneig[0] = neig[0];
        maxneig = neig[0];
        dist = info.N * info.propotion / info.ninterval;
        if ((sneig[0] >= dist && sneig[0] < dist * 1.05) || times > 30)
            DVtag[0] = 1;

        for (i = 1; i < np; i++)
        {
            if (DVtag[i] == 0)
                sneig[i] = sneig[i - 1] + neig[i];
            else
                neig[i] = sneig[i] - sneig[i - 1];
            if (maxneig < neig[i])
                maxneig = neig[i];

            dist = (i + 1.0) * info.N * info.propotion / info.ninterval;
            if (maxneig < M0 && sneig[i] > dist && sneig[i] < dist + 8)
                DVtag[i] = 1;
        }
        printf("myid  %d tag %d  %d  %d  %d %d  %d  %d   \n", myid, DVtag[0] , DVtag[1] , DVtag[2] , DVtag[3],  DVtag[4] , DVtag[5] , DVtag[6]);
        printf("myid  %d sneig %d  %d  %d  %d %d  %d  %d \n", myid, sneig[0] , sneig[1] , sneig[2] , sneig[3], sneig[4] , sneig[5] , sneig[6]  );

    }

    MPI_Bcast(DVtag, info.ninterval, MPI_INT, 0, info.intervalcomm);

    *tag = DVtag[myid % np];

}
