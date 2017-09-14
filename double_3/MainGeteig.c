/*******************************************
用于计算本区间内特征值
*********************************************/
#include "type.h"
void Geteig(SPMultiply SPMul, LinearSystem LSystem, Info info)
{
    double mytime1, mytime2, mytime3, mytime4;
    mytime4 = MPI_Wtime();

    int neig = info.neig[info.myid % info.ninterval];
    info.Emin = info.interval[0];
    info.Emax = info.interval[1];
    info.M0 = (int)(1.5 * neig);

    int N = info.N;
    int M0 = info.M0;
    int N_loc = info.N_loc[info.intervalid];
    if (info.intervalid == 0)
        printf("myid  %d Max:%lf  Min: %lf  neig %d\n", info.myid, info.Emax, info.Emin, neig);

    mytime1 = MPI_Wtime();
    LS_InitLS(&LSystem, info);
    mytime2 = MPI_Wtime();
    if (info.myid == 0)
        printf("\n LS_InitLS %lf \n\n", mytime2 - mytime1);

    ZMUMPS_STRUC_C id;

    mytime1 = MPI_Wtime();
    Mumps_Init(&id, info);
    mytime2 = MPI_Wtime();
    if (info.myid == 0)
        printf("\n Mumps_Init %lf \n\n", mytime2 - mytime1);

    mytime1 = MPI_Wtime();
    Mumps_Analize(&id, LSystem, info);
    mytime2 = MPI_Wtime();
    if (info.myid == 0)
        printf("\n Mumps_Analize %lf \n\n", mytime2 - mytime1);

    mytime1 = MPI_Wtime();
    Mumps_Factor(&id);
    mytime2 = MPI_Wtime();
    if (info.myid == 0)
        printf("\n Mumps_Factor %lf \n\n", mytime2 - mytime1);

    double *lambda = calloc(M0, sizeof(double));
    double *Q_loc = calloc(N_loc * M0, sizeof(double));
    double *Y_loc = calloc(N_loc * M0, sizeof(double));
    Feast_SetRandomMatrix(Q_loc, info);


    double maxres = 1;
    int times = 0;

    while(times < 3)// (maxres >  0.000001 && times < 10) // (maxres > 0.0001 && maxres != 0)
    {
        mytime3 = MPI_Wtime();

        times++;
        info.times = times;
        mytime1 = MPI_Wtime();
        LS_InitRhs(LSystem, Q_loc, info);
        mytime2 = MPI_Wtime();
        if (info.myid == 0)
            printf("\n LS_InitRhs %lf \n\n", mytime2 - mytime1);

        mytime1 = MPI_Wtime();
        Mumps_Solve(&id);
        mytime2 = MPI_Wtime();
        if (info.myid == 0)
            printf("\n Mumps_Solve %lf \n\n", mytime2 - mytime1);

        mytime1 = MPI_Wtime();
        Feast_GetProjectionOp(Y_loc, LSystem, info);
        mytime2 = MPI_Wtime();
        if (info.myid == 0)
            printf("\n Feast_GetProjectionOp %lf \n\n", mytime2 - mytime1);
        //4. 投影&&正交化
        mytime1 = MPI_Wtime();
        Feast_Project(Q_loc, lambda, Y_loc, SPMul, info);
        mytime2 = MPI_Wtime();
        if (info.myid == 0)
            printf("\n Feast_Project %lf \n\n", mytime2 - mytime1);


        mytime1 = MPI_Wtime();
        Feast_GetRes(&maxres, lambda, Q_loc, SPMul, info);
        mytime2 = MPI_Wtime();
        if (info.myid == 0)
            printf("\n Feast_GetRes %lf \n\n", mytime2 - mytime1);

        if (info.intervalid == 0)
            printf("***************times  %d max: %.25f   time: %.5f\n", times, maxres, mytime2 - mytime3);

    }
    printf("^^^^^^^^^^^^^^myid %d   times  %d max: %.25f   time: %.5f\n", info.myid, times, maxres, mytime2 - mytime4);


    // //print the eig
    double *eig;
    double Emin = info.Emin;
    double Emax = info.Emax;

    if (info.intervalid == 0)
    {
        eig = calloc(neig, sizeof(double));
        int i, j;
        j = 0;
        //统计特征值的数目 => j
        for (i = 0; i < M0; i++)
        {
            if (lambda[i] < Emax && lambda[i] > Emin)
            {
                j = j + 1;
            }
        }

        //判断估计的特征值数目是否正确
        if (j != neig)
        {
            printf("error in neig\n");
            free(eig);
            neig = j;
            eig = calloc(neig, sizeof(double));
        }

        j = 0;
        for (i = 0; i < M0; i++)
        {
            if (lambda[i] < Emax && lambda[i] > Emin)
            {
                eig[j] = lambda[i];
                j = j + 1;
            }
        }

        //打印特征值信息
        for (i = 0; i < neig; i++)
            printf("i %d\t  eig:%.25f\n", i, eig[i]);

    }


    Mumps_End(&id);

    free(lambda);
    free(Y_loc);
    free(Q_loc);
}
