#include "type.h"

void DV_End(Info info)
{
    int i;
    int *neig = info.neig;
    int *DVtag = info.DVtag;
    double *interval = info.interval;
    int ninterval = info.ninterval;
    int myid = info.myid;
    MPI_Barrier(MPI_COMM_WORLD);
    if (ninterval > 1)
    {
        if (myid < ninterval)
        {
            int id;

            for (i = 1; i < ninterval; i++)
                DVtag[0] = DVtag[0] + DVtag[i];

            int mpitag = 1991;
            MPI_Status status;
            if (myid == 0)
            {
                for (i = 1; i < ninterval; i++)
                    MPI_Recv(DVtag + i, 1, MPI_INT, i, mpitag, MPI_COMM_WORLD, &status);
            }
            else
            {
                MPI_Send(DVtag, 1, MPI_INT, 0, mpitag, MPI_COMM_WORLD);
            }

            mpitag = 1992;
            if (myid == 0)
            {
                int maxDVtag = DVtag[0];
                id = 0;
                for (i = 1; i < ninterval; i++)
                {
                    if (DVtag[i] > maxDVtag)
                    {
                        maxDVtag = DVtag[i];
                        id = i;
                    }
                }

                if (maxDVtag != ninterval)
                    printf("error in sumDVtag\n");
                for (i = 1; i < ninterval; i++)
                    MPI_Send(&id, 1, MPI_INT, i, mpitag, MPI_COMM_WORLD);
            }
            else
            {
                MPI_Recv( &id, 1, MPI_INT, 0, mpitag, MPI_COMM_WORLD, &status);
            }

            mpitag = 1993;
            if (myid == id)
            {
                for (i = 0; i < ninterval; i++)
                {
                    if (i != myid)
                        MPI_Send(neig, ninterval, MPI_INT, i, mpitag, MPI_COMM_WORLD);
                }
            }
            else
            {
                MPI_Recv(neig, ninterval, MPI_INT, id, mpitag, MPI_COMM_WORLD, &status);

            }


            mpitag = 1994 + myid;
            if (myid == 0)
            {
                MPI_Send(interval + 1, 1, MPI_DOUBLE, myid + 1, mpitag, MPI_COMM_WORLD);
            }
            else if (myid < ninterval - 1)
            {
                MPI_Recv(interval, 1, MPI_DOUBLE, myid - 1, mpitag - 1, MPI_COMM_WORLD, &status);

                MPI_Send(interval + 1, 1, MPI_DOUBLE, myid + 1, mpitag, MPI_COMM_WORLD);
            }
            else
            {
                MPI_Recv(interval, 1, MPI_DOUBLE, myid - 1, mpitag - 1, MPI_COMM_WORLD, &status);
            }
        }

    }

    MPI_Bcast(neig, info.ninterval, MPI_INT, 0, info.intervalcomm);
    MPI_Bcast(interval, 2, MPI_DOUBLE, 0, info.intervalcomm);

    if (myid == 0)
    {
        printf("neig:\n");
        for (i = 0; i < info.ninterval; i++)
            printf("%d \n", neig[i]);
        printf("\n");
    }

    // int number = info.ninterval;
    // int myid = info.myid;
    // int np = info.np;
    // int block, mpitag;
    // if (myid == 0)
    // {
    //     if (np > 1)
    //     {
    //         mpitag = 1;
    //         MPI_Send(neig, number, MPI_INT, mpitag, mpitag, MPI_COMM_WORLD);
    //         // printf("send  myid  %d  dest %d\n", myid, mpitag);
    //     }
    //     block = 1;
    //     while (block * 2 < np)
    //     {
    //         mpitag = block * 2;
    //         MPI_Send(neig, number, MPI_INT, mpitag, mpitag, MPI_COMM_WORLD);
    //         // printf("send  myid  %d  dest %d\n", myid, mpitag);
    //         block = block * 2;

    //     }
    // }
    // else
    // {
    //     MPI_Status status;
    //     block = 1;
    //     while (myid >= block * 2)
    //         block = block * 2;
    //     mpitag = myid % block;
    //     MPI_Recv(neig, number, MPI_INT, mpitag, myid, MPI_COMM_WORLD, &status);
    //     // printf("recv  myid  %d  source %d\n", myid, mpitag);

    //     block = block * 2;

    //     while (block + myid < np)
    //     {
    //         // if (myid == 1)
    //         // printf("~~~~~~~~~~~~~~~~~block %d np %d\n", block, np);
    //         mpitag = block + myid;
    //         MPI_Send(neig, number, MPI_INT, mpitag, mpitag, MPI_COMM_WORLD);
    //         // printf("send  myid  %d  dest %d\n", myid, mpitag);

    //         block = block * 2;
    //     }
    // }
}