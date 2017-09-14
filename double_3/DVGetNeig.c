/*******************************************
用于计算区间内特征值的个数
*********************************************/
#include "type.h"
void LS_InitLS(LinearSystem *LSystem, Info info);
void LS_InitRhs(LinearSystem LSystem, double *Q_loc, Info info);
int DV_GetNeig(SPMultiply SPMul, LinearSystem LSystem, Info info)
{
    double mytime1, mytime2, mytime3, mytime4;

    int neig;
    info.Emin = info.interval[0];
    info.Emax = info.interval[1];
    info.M0 = (int)(1.9 * info.N * info.propotion / info.ninterval);
    if (info.myid < info.ninterval)
        printf("Emax %lf  Emin %lf M0 %d\n", info.Emax, info.Emin, info.M0);


    int N = info.N;
    int M0 = info.M0;
    int N_loc = SPMul.N_loc;

    //1. 初始化线性求解器右端项
    mytime1 = MPI_Wtime();
    LS_InitLS(&LSystem, info);
    mytime2 = MPI_Wtime();
    if (info.intervalid == 0)
        printf("\n LS_InitLS %lf \n\n", mytime2 - mytime1);

    //2. 初始化MUMPS
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
    // //3.初始化 *Y_loc, *Q_loc
    double *Y_loc, *Q_loc;
    Q_loc = calloc(N_loc * M0, sizeof(double));
    Y_loc = calloc(N_loc * M0, sizeof(double));
    Feast_SetRandomMatrix(Q_loc, info);

   if (info.myid == 0)
        printf("\n Feast_SetRandomMatrix %lf \n\n", mytime2 - mytime1); 
    
int errneig, neigpre = 0;
    errneig = 1;
    int times = 0;

    mytime4 = MPI_Wtime();
    while (errneig > 0 || (neig == 0 && times < 3))
    {
        times++;
        mytime3 = MPI_Wtime();

        //4.计算投影算�
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
        neig = Feast_GetNumber(Y_loc, info);
        mytime2 = MPI_Wtime();
        if (info.myid == 0)
            printf("\n Feast_GetNumber %lf \n\n", mytime2 - mytime1);


        errneig = abs(neig - neigpre);
        neigpre = neig;

        if (info.intervalid == 0)
            printf("myid  %d %lf %lf  neig: %d   time  %lf\n", info.myid, info.Emin, info.Emax, neig, mytime2 - mytime3);

        //2.正交化
        mytime1 = MPI_Wtime();
        orth( Q_loc, Y_loc, info);
        mytime2 = MPI_Wtime();
        if (info.myid == 0)
            printf("\n orth %lf \n\n", mytime2 - mytime1);

    }
    if (info.myid == 0)
        printf("DV_GetNeig  time  %lf\n", mytime2 - mytime4 );

    //结束
    Mumps_End(&id);
    free(Y_loc);
    free(Q_loc);
    if (info.mumpsid == 0)
        free(LSystem.rhs);

    return neig;


}
