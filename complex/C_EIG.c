#include "complex.h"

void Feast_SetRandomMatrix(mumps_double_complex *Q_loc, int length, int myid);
void orth(mumps_double_complex *Q_loc, mumps_double_complex *Y_loc, int N_loc, int M0);
void Feast_GetProjectionOp(mumps_double_complex *Y_loc, mumps_double_complex *Q_loc, SparseMatrixLoc_mumps C_loc,  Info info);
void Feast_Project(mumps_double_complex *lambda, mumps_double_complex *Y_loc, SPMultiply SPMul, Info info);
void Feast_GetRes(double *maxres, mumps_double_complex *lambda, mumps_double_complex *Q_loc, SPMultiply SPMul, Info info);
void test_orth(mumps_double_complex *Q, int N,  int M);
void test_Add(SparseMatrixLoc_mumps C_loc);
void Eig(SPMultiply SPMul, SparseMatrixLoc_mumps C_loc, Info info)
{
    int i, j;
    int myid = info.myid;
    int M0 = info.M0;
    int N_loc = info.N_loc[myid];
    printf(" M0   %d   N_loc  %d\n", M0, N_loc );
    C_loc.nrhs = M0;
    /********************  Init the matrix **************************************/
    int length = M0 * N_loc;
    mumps_double_complex *Y_loc = calloc(length, sizeof(mumps_double_complex));
    Feast_SetRandomMatrix(Y_loc, length, myid);
    // printf("\nY_loc: \n");
    // for (i = 0; i < length; i++)
    //     printf("\t  %d    %lf  %lf\n", i, Y_loc[i].r, Y_loc[i].i);

    /********************  orth the matrix **************************************/
    mumps_double_complex *Q_loc = calloc(length, sizeof(mumps_double_complex));

    // printf("\nQ_loc: \n");
    // for (i = 0; i < length; i++)
    //     printf("\t  %d    %lf  %lf\n", i, Q_loc[i].r, Q_loc[i].i);

    /**********************  test  orth  ***************************************/
    //test_orth(Q_loc, N_loc, M0);


    /********************  comput  **************************************/
    double maxres;
    mumps_double_complex *lambda = calloc(M0, sizeof(mumps_double_complex));
    for (i = 1; i < 10; i++)
    {
        orth( Q_loc, Y_loc, N_loc, M0);

        /********************  Feast_GetProjectionOp **************************************/
        printf("/**************  Feast_GetProjectionOp ******************************/\n");
        Feast_GetProjectionOp(Y_loc, Q_loc, C_loc,  info);

        /********************  Feast_Project **************************************/
        printf(" /********************  Feast_Project **************************************/\n");
        Feast_Project(lambda, Y_loc, SPMul, info);

        /********************  Feast_GetRes **************************************/
        printf(" /********************  Feast_GetRes **************************************/\n");
        Feast_GetRes(&maxres, lambda, Q_loc, SPMul, info);
    }


}

/************************************
生成随机矩阵
*****************************************/
void Feast_SetRandomMatrix(mumps_double_complex *Y_loc, int length, int myid)
{
    int three = 3;
    int iseed[4] = {myid, myid + 3, myid + 9, 1 };
    zlarnv_(&three, iseed, &length, Y_loc);
}


/************************************
矩阵施密特正交
*****************************************/
void orth(mumps_double_complex *Q_loc, mumps_double_complex *Y_loc, int N_loc, int M0)
{
    int i, j, k;

    double mymodle, modle;
    double *mysum, *mysum_loc;

    mysum = calloc(2 * M0, sizeof(double));
    mysum_loc = calloc(2 * M0, sizeof(double));


    for (i = 0; i < M0; i++)
    {
        for (j = 0; j < N_loc; j++)
        {
            Q_loc[i * N_loc + j].r = Y_loc[i * N_loc + j].r + Q_loc[i * N_loc + j].r;
            Q_loc[i * N_loc + j].i = Y_loc[i * N_loc + j].i + Q_loc[i * N_loc + j].i;
        }

        //单位化
        mymodle = 0;
        for (j = 0; j < N_loc; j++)
        {
            mymodle = mymodle + Q_loc[i * N_loc + j].r * Q_loc[i * N_loc + j].r + Q_loc[i * N_loc + j].i * Q_loc[i * N_loc + j].i;
        }
        MPI_Allreduce(&mymodle, &modle, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        modle = sqrt(modle);

        for (j = 0; j < N_loc; j++)
        {
            Q_loc[i * N_loc + j].r = Q_loc[i * N_loc + j].r / modle;
            Q_loc[i * N_loc + j].i = Q_loc[i * N_loc + j].i / modle;
        }

        //点积
        for (j = i + 1; j < M0; j++)
        {
            mysum_loc[2 * j] = 0; mysum_loc[2 * j + 1] = 0;
            mysum[2 * j] = 0; mysum[2 * j + 1] = 0;

            for (k = 0; k < N_loc; k++)
            {
                mysum_loc[2 * j] = mysum_loc[2 * j] + Q_loc[i * N_loc + k].r * Y_loc[j * N_loc + k].r + Q_loc[i * N_loc + k].i * Y_loc[j * N_loc + k].i;
                mysum_loc[2 * j + 1] = mysum_loc[2 * j + 1] + Q_loc[i * N_loc + k].r * Y_loc[j * N_loc + k].i - Q_loc[i * N_loc + k].i * Y_loc[j * N_loc + k].r;

            }
        }

        MPI_Allreduce(mysum_loc + 2 * i + 2, mysum + 2 * i + 2, 2 * (M0 - i - 1), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        for (j = i + 1; j < M0; j++)
        {
            for (k = 0; k < N_loc; k++)
            {
                Q_loc[j * N_loc + k].r = Q_loc[j * N_loc + k].r - mysum[2 * j] * Q_loc[i * N_loc + k].r + mysum[2 * j + 1] * Q_loc[i * N_loc + k].i;
                Q_loc[j * N_loc + k].i = Q_loc[j * N_loc + k].i - mysum[2 * j] * Q_loc[i * N_loc + k].i - mysum[2 * j + 1] * Q_loc[i * N_loc + k].r;
            }
        }
    }

}

/************************************
计算投影矩阵
*****************************************/
void Mumps_Init(ZMUMPS_STRUC_C *id);
void Mumps_Analize(ZMUMPS_STRUC_C *id, SparseMatrixLoc_mumps C_loc, int myid);
void Mumps_FactorSolve(ZMUMPS_STRUC_C *id);
void Mumps_End(ZMUMPS_STRUC_C *id);
void C_Add(SparseMatrixLoc_mumps C_loc, mumps_double_complex Ze);
void Gauss(int nGauss, int i, double *xe, double *we);
void test_Mumps_single(SparseMatrixLoc_mumps C_loc, mumps_double_complex *rhs, mumps_double_complex *result);
void Feast_GetProjectionOp(mumps_double_complex *Y_loc, mumps_double_complex *Q_loc, SparseMatrixLoc_mumps C_loc,  Info info)
{
    int myid = info.myid;
    int N = info.N;
    int M0 = info.M0;
    int *N_loc = info.N_loc;
    mumps_double_complex Emax = info.Emax;
    mumps_double_complex Emin = info.Emin;


    int i, j, k, nrow, mpitag = 1;
    MPI_Status status;
    MPI_Datatype MPI_complex;
    MPI_Type_contiguous(2, MPI_DOUBLE, &MPI_complex);
    MPI_Type_commit(&MPI_complex);

    int myN_loc = N_loc[myid];
    /******************************收集分布存储的正交矩阵*******************************/
    mumps_double_complex *Q;
    if (myid == 0)
        Q = calloc(N * M0, sizeof(mumps_double_complex));

    if (myid == 0)
    {
        for (i = 0; i < M0; i++)
            for (j = 0; j < myN_loc; j++)
                Q[i * N + j] = Q_loc[i * myN_loc + j];

        int offset, nrow;
        offset = myN_loc;
        mumps_double_complex *recvbuf;
        recvbuf = calloc(myN_loc * M0, sizeof(mumps_double_complex));

        for (i = 1; i < info.np; i++)
        {
            nrow = N_loc[i];
            MPI_Recv(recvbuf, nrow * M0, MPI_complex,  i, mpitag, MPI_COMM_WORLD, &status);
            for (j = 0; j < M0; j++)
                for (k = 0; k < nrow; k++)
                    Q[offset + j * N + k] = recvbuf[j * nrow + k];
            offset = offset + nrow;
        }
        free(recvbuf);
    }
    else
        MPI_Send(Q_loc, myN_loc * M0, MPI_complex, 0, mpitag, MPI_COMM_WORLD);

    /******************************计算投影矩阵Y*******************************/
    //右端项初始化
    mumps_double_complex *rhs;
    if (info.myid == 0)
    {
        rhs = calloc(N * M0, sizeof(mumps_double_complex));
        C_loc.rhs = rhs;
    }
    mumps_double_complex *Y;
    if (info.myid == 0)
        Y = calloc(N * M0, sizeof(mumps_double_complex));
    //mumps 求解器初始化
    ZMUMPS_STRUC_C id;
    Mumps_Init(&id);
    Mumps_Analize(&id, C_loc, myid);

    mumps_double_complex Ze, jac;

    mumps_double_complex Emid;
    Emid.r = 0.5 * (Emax.r + Emin.r);
    Emid.i = 0.5 * (Emax.i + Emin.i);

    double r;
    r = sqrt((Emax.r - Emin.r) * (Emax.r - Emin.r) + (Emax.i - Emin.i) * (Emax.i - Emin.i));
    for (i = 0; i < info.nGau; i++)
    {
        //右端项
        if (info.myid == 0)
            for (j = 0; j < N * M0; j++)
                rhs[j] = Q[j];

        //复数矩阵加减（ZeI-A）
        double xe, we, theta;
        Gauss(info.nGau,  i, &xe,  &we);
        theta = 3.1415926535898 * (xe + 1) / 2;

        mumps_double_complex Ze;
        Ze.r = 0.5 * Emid.r + 0.5 * r * cos(theta);
        Ze.i = 0.5 * Emid.i + 0.5 * r * sin(theta);

        C_Add(C_loc, Ze);

        //Mumps 求解
        Mumps_FactorSolve(&id);
        //test_Mumps_single(C_loc, Q, rhs);

        //将结果进行累加
        jac.r = 0.25 * we * r  * cos(theta);
        jac.i = 0.25 * we * r  * sin(theta);
        if (info.myid == 0)
        {
            for (j = 0; j < N * M0; j++)
            {
                Y[j].r = Y[j].r + jac.r * rhs[j].r - jac.i * rhs[j].i;
                Y[j].i = Y[j].i + jac.r * rhs[j].i + jac.i * rhs[j].r;
            }
        }


        //右端项
        if (info.myid == 0)
            for (j = 0; j < N * M0; j++)
                rhs[j] = Q[j];

        //复数矩阵加减（ZeI-A）
        Ze.i = -Ze.i;
        C_Add(C_loc, Ze);

        //Mumps 求解
        Mumps_FactorSolve(&id);
        // test_Mumps_single(C_loc, Q, rhs);
        //将结果进行累加
        jac.i = -jac.i;
        if (info.myid == 0)
        {
            for (j = 0; j < N * M0; j++)
            {
                Y[j].r = Y[j].r + jac.r * rhs[j].r - jac.i * rhs[j].i;
                Y[j].i = Y[j].i + jac.r * rhs[j].i + jac.i * rhs[j].r;
            }
        }

    }

    if (info.myid == 0)
        free(rhs);
    /******************************投影矩阵Y分布存储得Y_loc*******************************/
    if (info.myid == 0)
    {
        for (i = 0; i < M0; i++)
            for (j = 0; j < myN_loc; j++)
                Y_loc[i * myN_loc + j] = Y[i * N + j];

        int offset, nrow;
        offset = myN_loc;
        mumps_double_complex *sendbuf;
        sendbuf = calloc(myN_loc * M0, sizeof(mumps_double_complex));
        for (i = 1; i < info.np; i++)
        {
            nrow = N_loc[i];
            for (j = 0; j < M0; j++)
                for (k = 0; k < nrow; k++)
                    sendbuf[j * nrow + k] = Y[offset + j * N + k];
            MPI_Send(sendbuf,   nrow * M0, MPI_complex, i, mpitag, MPI_COMM_WORLD);
            offset = offset + nrow;
        }
        free(sendbuf);
    }
    else
    {
        MPI_Recv(Y_loc, myN_loc * M0, MPI_complex,  0, mpitag, MPI_COMM_WORLD, &status);
    }
    if (info.myid == 0)
        free(Y);
}


/************************************
根据计算得到的投影矩阵进行投影,并且将投影矩阵正交化
*****************************************/
void mul(mumps_double_complex *result, mumps_double_complex *Y_loc, SPMultiply SPMul, Info info);
void geteig( mumps_double_complex *A,  mumps_double_complex *B ,  mumps_double_complex *lambda, int N);
void Feast_Project(mumps_double_complex *lambda, mumps_double_complex *Y_loc, SPMultiply SPMul, Info info)
{
    MPI_Datatype MPI_complex;
    MPI_Type_contiguous(2, MPI_DOUBLE, &MPI_complex);
    MPI_Type_commit(&MPI_complex);

    int myid = info.myid;
    int M0 = info.M0;
    int i, j;
    int N_loc = info.N_loc[myid];


    char char_N = 'N', char_T = 'T', char_V = 'V', char_U = 'U';
    mumps_double_complex one, zero;
    one.r = 1; one.i = 0;
    zero.r = 0; zero.i = 0;
    int int_one = 1;

    mumps_double_complex *YT_loc = calloc(N_loc * M0, sizeof(mumps_double_complex));
    for (i = 0; i < N_loc * M0; i++)
    {
        YT_loc[i].r = Y_loc[i].r;
        YT_loc[i].i = -Y_loc[i].i;
    }
    /****************************** 投影Aq=YT*A*Y 和Bq=YT*Y *******************************/
    int error, lwork = N_loc * M0;
    mumps_double_complex *Aq_loc, *Aq, *Bq_loc, *Bq, *work_loc;
    Aq_loc = calloc(M0 * M0, sizeof(mumps_double_complex));
    Bq_loc = calloc(M0 * M0, sizeof(mumps_double_complex));
    work_loc = calloc(lwork, sizeof(mumps_double_complex));
    if (myid == 0)
    {
        Aq = calloc(M0 * M0, sizeof(mumps_double_complex));
        Bq = calloc(M0 * M0, sizeof(mumps_double_complex));
    }

    //Aq=Y^t A Y
    mul(work_loc, Y_loc, SPMul, info);
    zgemm_(&char_T, &char_N, &M0, &M0, &N_loc, &one, YT_loc, &N_loc, work_loc, &N_loc, &zero, Aq_loc, &M0); //Aq = Y^T * work
    MPI_Reduce(Aq_loc, Aq, 2 * M0 * M0, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    //Bq=Y^T * Y
    zgemm_(&char_T, &char_N, &M0, &M0, &N_loc, &one, YT_loc, &N_loc, Y_loc, &N_loc, &zero, Bq_loc, &M0);
    MPI_Reduce(Bq_loc, Bq, 2 * M0 * M0, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    /****************************** 计算Rize对Aq*X=λ*Bq*X *******************************/
    if (info.myid == 0)
    {
        geteig( Aq,  Bq , lambda , M0);
    }


    free(work_loc);
    free(Aq_loc);
    free(Bq_loc);
    free(YT_loc);
    free(Aq);
    free(Bq);


}

double distant(mumps_double_complex a, mumps_double_complex b)
{
    double c = (a.r - b.r) * (a.r - b.r) + (a.i - b.i) * (a.i - b.i);
    c = sqrt(c);
    return c;

}
void Feast_GetRes(double *maxres, mumps_double_complex *lambda, mumps_double_complex *Q_loc, SPMultiply SPMul, Info info)
{
    MPI_Datatype MPI_complex;
    MPI_Type_contiguous(2, MPI_DOUBLE, &MPI_complex);
    MPI_Type_commit(&MPI_complex);

    int myid = info.myid;
    int M0 = info.M0;
    int N_loc = info.N_loc[myid];

    mumps_double_complex mid;
    mid.r = (info.Emax.r + info.Emin.r) / 2;
    mid.i = (info.Emax.i + info.Emin.i) / 2;

    double offset;
    offset = distant(info.Emax, info.Emin) / 2;

    int i, j, k;
    mumps_double_complex *work_loc;
    double  *res_loc, *res, *x_loc, *x;
    work_loc = calloc(N_loc * M0, sizeof(mumps_double_complex));

    res_loc = calloc(M0, sizeof(double));
    x_loc = calloc(M0, sizeof(double));

    res = calloc(M0, sizeof(double));
    x = calloc(M0, sizeof(double));

    mul(work_loc, Q_loc, SPMul, info);
    MPI_Bcast(lambda, M0, MPI_complex, 0, MPI_COMM_WORLD);

    for (i = 0; i < M0; i++)
    {
        for (j = 0; j < N_loc; j++)
        {
            work_loc[i * N_loc + j].r =  work_loc[i * N_loc + j].r - lambda[i].r * Q_loc[i * N_loc + j].r + lambda[i].i * Q_loc[i * N_loc + j].i;
            work_loc[i * N_loc + j].i =  work_loc[i * N_loc + j].i - lambda[i].i * Q_loc[i * N_loc + j].r - lambda[i].r * Q_loc[i * N_loc + j].i;
        }
    }

    for (i = 0; i < M0; i++)
    {
        res_loc[i] = 0;

        x_loc[i] = 0;

        for (j = 0; j < N_loc; j++)
        {
            res_loc[i] = res_loc[i] + work_loc[i * N_loc + j].r * work_loc[i * N_loc + j].r + work_loc[i * N_loc + j].i * work_loc[i * N_loc + j].i;
            x_loc[i] = x_loc[i] + Q_loc[i * N_loc + j].r * Q_loc[i * N_loc + j].r + Q_loc[i * N_loc + j].i * Q_loc[i * N_loc + j].i;
        }
    }

    MPI_Reduce(res_loc, res, M0, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(x_loc,   x,   M0, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);


    if (info.myid == 0)
    {
        for (i = 0; i < M0; i++)
            res[i] = res[i] / x[i];

        int neig = 0;
        int maxid;
        *maxres = 0;
        for (i = 0; i < M0; i++)
        {
            if (distant(lambda[i], mid) < offset)
            {
                neig++;
                if (res[i] > *maxres)
                {
                    maxid = i;
                    *maxres = res[i];
                }
            }
        }
        //printf("\t\t maxres: %lf  neig %d %d \n", *maxres, neig, info.neig[info.myid % info.ninterval]);

        /*if (neig > info.neig[info.myid % info.ninterval])
        {
            for (i = 0; i < neig - info.neig[info.myid % info.ninterval]; i++)
            {
                res[maxid] = 0;
                lambda[maxid] = Emax + 1;

                *maxres = 0;
                for (j = 0; j < M0; j++)
                {
                    if (lambda[j] < Emax && lambda[j] > Emin)
                    {
                        if (res[j] > *maxres)
                        {
                            maxid = j;
                            *maxres = res[j];
                        }
                    }

                }
            }
            printf("\t\tcorect  maxres: %.25f\n", *maxres);
        }*/

     printf("\t\tcorect  maxres: %.25f\n", *maxres);

    }
    free(work_loc);
    free(res);
    free(x);
    free(res_loc);
    free(x_loc);


    MPI_Bcast(maxres, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

}


/**********************************************
测试矩阵正交性
*******************************************/
void test_orth(mumps_double_complex *Q, int N,  int M)
{
    char char_T = 'T', char_N = 'N';
    mumps_double_complex one, zero;
    one.r = 1; one.i = 0;
    zero.r = 0; zero.i = 0;
    mumps_double_complex *result = calloc(M * M, sizeof(mumps_double_complex));

    mumps_double_complex *QT    = calloc(N * M, sizeof(mumps_double_complex));
    int i, j;
    for (i = 0; i < M; i++)
    {
        for (j = 0; j < N; j++)
        {
            QT[i * N + j].r =  Q[i * N + j].r;
            QT[i * N + j].i = -Q[i * N + j].i;
        }
    }
    zgemm_(&char_T, &char_N, &M, &M, &N, &one, QT, &N, Q, &N, &zero, result, &M); //Aq = Y^T * work

    MPI_Datatype MPI_complex;
    MPI_Type_contiguous(2, MPI_DOUBLE, &MPI_complex);
    MPI_Type_commit(&MPI_complex);
    int myid, np;
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &np);
    if (myid == 0)
    {
        MPI_Status status;
        mumps_double_complex *recvbuf = calloc(M * M, sizeof(mumps_double_complex));
        for (i = 1; i < np; i++)
        {
            MPI_Recv(recvbuf, M * M, MPI_complex, i, i, MPI_COMM_WORLD, &status);
            for (j = 0; j < M * M; j++)
            {
                result[j].r = result[j].r + recvbuf[j].r;
                result[j].i = result[j].i + recvbuf[j].i;
            }
        }
    }
    else
    {
        MPI_Send(result, M * M, MPI_complex, 0, myid, MPI_COMM_WORLD);
    }

    if (myid == 0)
    {
        for (i = 0; i < M; i++)
        {
            for (j = 0; j < M; j++)
            {
                printf("\t %lf %lf\n", result[i * M + j].r, result[i * M + j].i);
            }
            printf("\n\n");
        }
    }
    free(result);
    free(QT);
}



void test_Add(SparseMatrixLoc_mumps C_loc)
{
    int N_loc = (C_loc).N_loc;
    int mynnz_loc = (C_loc).nnz_loc;
    int *ic_loc = (C_loc).ic_loc;
    int *jc_loc = (C_loc).jc_loc;
    mumps_double_complex *c_loc =  (C_loc).c_loc;
    mumps_double_complex *c_diag = (C_loc).c_diag;

    int i;
    printf(" /***************************** C_loc  ***********************************/\n");
    printf("\t mynnz_loc: %d\n", mynnz_loc);

    printf("\t ic_loc: ");
    for (i = 0; i < mynnz_loc; i++)
        printf("%d ", ic_loc[i]);
    printf("\n");

    printf("\t jc_loc: ");
    for (i = 0; i < mynnz_loc; i++)
        printf("%d ", jc_loc[i]);
    printf("\n");

    printf("\t c_loc: ");
    for (i = 0; i < mynnz_loc; i++)
    {
        printf("%lf  %lf\n", c_loc[i].r, c_loc[i].i);
        printf("\t        ");
    }
    printf("\n");

    printf("\t c_diag: ");
    for (i = 0; i < N_loc; i++)
    {
        printf("%lf  %lf\n", c_diag[i].r, c_diag[i].i);
        printf("\t         ");
    }
    printf("\n");
}


void test_Mumps_single(SparseMatrixLoc_mumps C_loc, mumps_double_complex *rhs, mumps_double_complex *result)
{
    int N = (C_loc).N_loc;
    int M = (C_loc).nrhs;
    int nnz = (C_loc).nnz_loc;
    int *ic = (C_loc).ic_loc;
    int *jc = (C_loc).jc_loc;
    mumps_double_complex *c =  (C_loc).c_loc;

    mumps_double_complex *myrhs = calloc(N * M, sizeof(mumps_double_complex));

    int *myic = calloc(N + 1, sizeof(int));
    int i, j;
    // printf("!!!!!!!!!!!!!!!!!!!!!!!!C: \n");
    for (i = 0; i < nnz; i++)
    {
        j = ic[i];
        myic[j] = myic[j] + 1;
        //    printf("\t   %d   %d   %lf  %lf\n", ic[i], jc[i], c[i].r, c[i].i);
    }
    myic[0] = 1;
    for (i = 1; i < N + 1; i++)
    {
        myic[i] = myic[i] + myic[i - 1];
    }

    char char_N = 'N';
    mumps_double_complex one, zero;
    one.r = 1; one.i = 0;
    zero.r = 0; zero.i = 0;
    char matdescra[6] = {'G', 'U', 'U', 'F'}; //矩阵A的类型
    mkl_zcsrmm (&char_N, &N, &M, &N, &one, matdescra, c, jc, myic, myic + 1, result, &N, &zero, myrhs, &N);

    mumps_double_complex error;
    error.r = 0; error.i = 0;

    for (i = 0; i < N * M; i++)
    {
        error.r = error.r + fabs(myrhs[i].r - rhs[i].r);
        error.i = error.i + fabs(myrhs[i].i - rhs[i].i);
    }
    printf("error   %lf  %lf  \n", error.r , error.i);
    // printf("result!!!!!!!!\n");
    // for (i = 0; i < N * M; i++)
    // {
    //     printf("\t %lf  %lf \n",result[i].r, result[i].i);
    // }
    // printf("rhs!!!!!!!!\n");
    // for (i = 0; i < N * M; i++)
    // {
    //     printf("\t %lf  %lf  %lf  %lf \n", myrhs[i].r, myrhs[i].i, rhs[i].r, rhs[i].i);
    // }

    free(myrhs);
    free(myic);
}


/********************************************二级程序**********************************/
/************************************
MUMPS 计算
*****************************************/
void Mumps_Init(ZMUMPS_STRUC_C *id)
{
    (*id).par = 1;
    (*id).sym = 0;
    (*id).comm_fortran = MPI_Comm_c2f(MPI_COMM_WORLD);

    (*id).job = -1;
    zmumps_c(id);
    if ((*id).info[0] != 0)
        printf("InitMumps info: %d %d \n", (*id).info[0], (*id).info[1]);
}
void Mumps_Analize(ZMUMPS_STRUC_C *id, SparseMatrixLoc_mumps C_loc, int myid)
{
    (*id).icntl[2] = -1;
    (*id).icntl[3] = 1;

    (*id).icntl[4] = 0;
    (*id).icntl[17] = 3;

    if (myid == 0)
    {
        (*id).rhs = C_loc.rhs;
        (*id).nrhs = C_loc.nrhs;
        (*id).lrhs = C_loc.N;
        (*id).n = C_loc.N;

    }

    (*id).nz_loc = C_loc.nnz_loc;
    (*id).irn_loc = C_loc.ic_loc;
    (*id).jcn_loc = C_loc.jc_loc;
    (*id).a_loc = C_loc.c_loc;
    (*id).job = 1;
    zmumps_c(id);

    if ((*id).info[0] != 0)
        printf("AnalizeMumps info: %d %d \n", (*id).info[0], (*id).info[1]);
}
void Mumps_FactorSolve(ZMUMPS_STRUC_C *id)
{
    (*id).job = 5;
    zmumps_c(id);
    if ((*id).info[0] != 0)
        printf("Mumps_FactorSolve info: %d %d \n", (*id).info[0], (*id).info[1]);
}
void Mumps_End(ZMUMPS_STRUC_C *id)
{
    (*id).job = -2;
    zmumps_c(id);
    if ((*id).info[0] != 0)
        printf("SolveMumps info: %d %d \n", (*id).info[0], (*id).info[1]);
}

/************************************
MUMPS 计算
*****************************************/
void C_Add(SparseMatrixLoc_mumps C_loc, mumps_double_complex Ze)
{
    int nnz_loc = C_loc.nnz_loc;

    int *ic_loc = C_loc.ic_loc;
    int *jc_loc = C_loc.jc_loc;
    mumps_double_complex *c_loc = C_loc.c_loc;
    mumps_double_complex *c_diag = C_loc.c_diag;

    int i;
    for (i = 0; i < nnz_loc; i++)
    {
        if (ic_loc[i] == jc_loc[i])
        {
            c_loc[i].r = Ze.r - c_diag[i].r;
            c_loc[i].i = Ze.i - c_diag[i].i;
        }
    }
}

/************************************
矩阵乘
*****************************************/
void mul(mumps_double_complex *result, mumps_double_complex *Y_loc, SPMultiply SPMul, Info info)
{
    int M0 = info.M0;
    int np = info.np;
    int myid = info.myid;


    int N_loc = SPMul.N_loc;

    int nnz_in = SPMul.nnz_in;
    int *ia_in = SPMul.ia_in;
    int *ja_in = SPMul.ja_in;
    mumps_double_complex *a_in = SPMul.a_in;

    int nnz_out = SPMul.nnz_out;
    int *ia_out = SPMul.ia_out;
    int *ja_out = SPMul.ja_out;
    mumps_double_complex *a_out = SPMul.a_out;

    int nsend = SPMul.nsend;
    int nrecv = SPMul.nrecv;
    int *recv_rowid = SPMul.recv_rowid;
    int *recv_count = SPMul.recv_count;
    int *send_rowid = SPMul.send_rowid;
    int *send_count = SPMul.send_count;

    MPI_Status status;
    MPI_Datatype MPI_complex;
    MPI_Type_contiguous(2, MPI_DOUBLE, &MPI_complex);
    MPI_Type_commit(&MPI_complex);

    //
    int i, j;

    mumps_double_complex *B_out;
    if (nrecv == 0)
        B_out = calloc(1, sizeof(mumps_double_complex));
    else
        B_out = calloc(nrecv * M0, sizeof(mumps_double_complex));

    int rowid, colid, idsend, idrecv, sendoffset, recvoffset, number, mpitag = 0;

    mumps_double_complex *sendbuf;
    sendbuf = calloc(N_loc * M0, sizeof(mumps_double_complex));
    idsend = 0;
    idrecv = 0;
    recvoffset = 0;
    sendoffset = 0;

    while (recv_count[idsend] == 0 && idsend < np)
        idsend++;
    while (send_count[idrecv] == 0 && idrecv < np)
        idrecv++;
    while (idsend < np || idrecv < np)
    {
        if (idsend == idrecv)
        {
            if (idrecv < myid)
            {
                MPI_Recv(B_out + recvoffset, recv_count[idrecv] * M0, MPI_complex, idrecv, mpitag, MPI_COMM_WORLD, &status);
                recvoffset = recvoffset + recv_count[idrecv] * M0;
                idrecv++;
            }

            if (idsend > myid)
            {
                for (i = 0; i < send_count[idsend]; i++)
                {
                    rowid = send_rowid[sendoffset];
                    for (j = 0; j < M0; j++)
                        sendbuf[i * M0 + j] = Y_loc[rowid + j * N_loc];
                    sendoffset++;
                }
                MPI_Send(sendbuf, send_count[idsend] * M0, MPI_complex, idsend, mpitag, MPI_COMM_WORLD);
                idsend++;
            }
        }
        else
        {
            if (idrecv < idsend || idsend == np)
            {
                MPI_Recv(B_out + recvoffset, recv_count[idrecv] * M0, MPI_complex, idrecv, mpitag, MPI_COMM_WORLD, &status);
                recvoffset = recvoffset + recv_count[idrecv] * M0;
                idrecv++;
            }
            else if (idsend < idrecv || idrecv == np)
            {
                for (i = 0; i < send_count[idsend]; i++)
                {
                    rowid = send_rowid[sendoffset];
                    for (j = 0; j < M0; j++)
                        sendbuf[i * M0 + j] = Y_loc[rowid + j * N_loc];
                    sendoffset++;
                }
                MPI_Send(sendbuf, send_count[idsend] * M0, MPI_complex, idsend, mpitag, MPI_COMM_WORLD);
                idsend++;
            }
        }
        while (recv_count[idsend] == 0 && idsend < np)
            idsend++;
        while (send_count[idrecv] == 0 && idrecv < np)
            idrecv++;
    }

    free(sendbuf);

    //内点相乘
    char char_N = 'N', char_T = 'T', char_V = 'V', char_U = 'U';
    mumps_double_complex one, zero;
    one.r = 1; one.i = 0;
    zero.r = 0; zero.i = 0;
    char matdescra[6] = {'G', 'U', 'U', 'F'}; //矩阵A的类型
    mkl_zcsrmm (&char_N, &N_loc, &M0, &N_loc, &one, matdescra, a_in, ja_in, ia_in, ia_in + 1, Y_loc, &N_loc, &zero, result, &N_loc);//work = A * Y;

    //外点相乘
    for (i = 0; i < nnz_out; i++)
    {
        rowid = ia_out[i];
        colid = ja_out[i];
        for (j = 0; j < M0; j++)
        {
            result[rowid + j * N_loc].r =  result[rowid + j * N_loc].r + a_out[i].r * B_out[colid * M0 + j].r;
            result[rowid + j * N_loc].i =  result[rowid + j * N_loc].i + a_out[i].i * B_out[colid * M0 + j].i;
        }
    }

    free(B_out);
}


/*****************************************
子空间特征值求解
****************************************/
void geteig( mumps_double_complex *A,  mumps_double_complex *B ,  mumps_double_complex *lambda, int N)
{
    char char_N = 'N', char_V = 'V';
    int i, j;
    mumps_double_complex *alpha, *beta;
    double *rwork;
    alpha = calloc(N, sizeof(mumps_double_complex));
    beta  = calloc(N, sizeof(mumps_double_complex));
    rwork = calloc(8 * N, sizeof(double));
    int error;

    int lwork = -1;
    mumps_double_complex zlwork;
    zgegv_(&char_N, &char_N, &N, A, &N, B, &N, alpha, beta, &zlwork, &N, &zlwork, &N, &zlwork, &lwork, rwork, &error);
    if (error != 0)
        printf("error in zgegv_  1 error %d\n", error);
    lwork = zlwork.r;


    mumps_double_complex *work;
    work  = calloc(lwork, sizeof(mumps_double_complex));


    zgegv_(&char_N, &char_N, &N, A, &N, B, &N, alpha, beta, &zlwork, &N, &zlwork, &N, work, &lwork, rwork, &error);
    if (error != 0)
        printf("error in dsygv_ error %d\n", error);


    double fenzi;
    // printf("alpha beta : \n");
    // for (i = 0; i < N; i++)
    // {
    //     printf("\t   i  %d   %lf %lf   %lf %lf\n", i, alpha[i].r, alpha[i].i, beta[i].r, beta[i].i);
    // }

    printf("lambda : \n");
    for (i = 0; i < N; i++)
    {
        fenzi = beta[i].r * beta[i].r + beta[i].i * beta[i].i;
        lambda[i].r = (alpha[i].r * beta[i].r + alpha[i].i * beta[i].i) / fenzi;
        lambda[i].i = (alpha[i].i * beta[i].r - alpha[i].r * beta[i].i) / fenzi;
        printf("\t   i  %d   %lf   %lf \n", i, lambda[i].r, lambda[i].i);
    }

    free(alpha);
    free(beta);
    free(rwork);
    free(work);
}
/************************************
高斯点
*****************************************/
void Gauss(int nGauss, int i, double *xe, double *we)
{
    switch (nGauss)
    {

    case (3):
        switch (i)
        {
        case (0):
            *xe = 0;
            *we = 0.888888888888888888889;
            break;
        case (1):
            *xe = sqrt(0.6);
            *we = 0.555555555555555555555;
            break;
        case (2):
            *xe = -sqrt(0.6);
            *we = 0.5555555555555555555555;
            break;
        }
        break;


    case (4):
        switch (i)
        {
        case (0):
            *xe = -0.339981043584856264792;
            *we = 0.652145154862546142644;
            break;
        case (1):
            *xe = 0.339981043584856264792;
            *we = 0.652145154862546142644;
            break;
        case (2):
            *xe = -0.861136311594052575248;
            *we = 0.347854845137453857383;
            break;
        case (3):
            *xe = 0.861136311594052575248;
            *we = 0.347854845137453857383;
            break;
        }
        break;


    case (5):
        switch (i)
        {
        case (0):
            *xe = 0;
            *we = 0.568888888888888888889;
            break;
        case (1):
            *xe = 0.538469310105683091018;
            *we = 0.4786287;
            break;
        case (2):
            *xe = -0.538469310105683091018;
            *we = 0.47862878;
            break;
        case (3):
            *xe = 0.906179845938663992811;
            *we = 0.2369269;
            break;
        case (4):
            *xe = -0.906179845938663992811;
            *we = 0.2369269;
            break;
        }
        break;


    case (6):
        switch (i)
        {
        case (0):
            *xe = -0.661209386466264513688;
            *we = 0.360761573048138607569;
            break;
        case (1):
            *xe = 0.661209386466264513688;
            *we = 0.360761573048138607569;
            break;
        case (2):
            *xe = -0.238619186083196908630;
            *we = 0.467913934572691047389;
            break;
        case (3):
            *xe = 0.238619186083196908630;
            *we = 0.467913934572691047389;
            break;
        case (4):
            *xe = -0.932469514203152027832;
            *we = 0.171324492379170345043;
            break;
        case (5):
            *xe = 0.932469514203152027832;
            *we = 0.171324492379170345043;
            break;
        }
        break;


    case (8):
        switch (i)
        {
        case (0):
            *xe = -0.183434642495649804936;
            *we = 0.362683783378361982976;
            break;
        case (1):
            *xe = 0.183434642495649804936;
            *we = 0.362683783378361982976;
            break;
        case (2):
            *xe = -0.525532409916328985830;
            *we = 0.313706645877887287338;
            break;
        case (3):
            *xe = 0.525532409916328985830;
            *we = 0.313706645877887287338;
            break;
        case (4):
            *xe = -0.796666477413626739567;
            *we = 0.222381034453374470546;
            break;
        case (5):
            *xe = 0.796666477413626739567;
            *we = 0.222381034453374470546;
            break;
        case (6):
            *xe = -0.960289856497536231661;
            *we = 0.101228536290376259154;
            break;
        case (7):
            *xe = 0.960289856497536231661;
            *we = 0.101228536290376259154;
            break;
        }
        break;


    case (10):
        switch (i)
        {
        case (0):
            *xe = -0.148874338981631210881;
            *we = 0.295524224714752870187;
            break;
        case (1):
            *xe = 0.148874338981631210881;
            *we = 0.295524224714752870187;
            break;
        case (2):
            *xe = -0.433395394129247190794;
            *we = 0.269266719309996355105;
            break;
        case (3):
            *xe = 0.433395394129247190794;
            *we = 0.269266719309996355105;
            break;
        case (4):
            *xe = -0.6794095682990244062070;
            *we = 0.219086362515982044000;
            break;
        case (5):
            *xe = 0.6794095682990244062070;
            *we = 0.219086362515982044000;
            break;
        case (6):
            *xe = -0.865063366688984510759;
            *we = 0.149451349150580593150;
            break;
        case (7):
            *xe = 0.865063366688984510759;
            *we = 0.149451349150580593150;
            break;
        case (8):
            *xe = -0.973906528517171720066;
            *we = 0.0666713443086881375920;
            break;
        case (9):
            *xe = 0.973906528517171720066;
            *we = 0.0666713443086881375920;
            break;
        }
        break;


    case (12):
        switch (i)
        {
        case (0):
            *xe = -0.125233408511468915478;
            *we = 0.249147045813402785006;
            break;
        case (1):
            *xe = 0.125233408511468915478;
            *we = 0.249147045813402785006;
            break;
        case (2):
            *xe = -0.367831498998180193757;
            *we = 0.233492536538354808758;
            break;
        case (3):
            *xe = 0.367831498998180193757;
            *we = 0.233492536538354808758;
            break;
        case (4):
            *xe = -0.587317954286617447312;
            *we = 0.203167426723065921743;
            break;
        case (5):
            *xe = 0.587317954286617447312;
            *we = 0.203167426723065921743;
            break;
        case (6):
            *xe = -0.769902674194304687059;
            *we = 0.160078328543346226338;
            break;
        case (7):
            *xe = 0.769902674194304687059;
            *we = 0.160078328543346226338;
            break;
        case (8):
            *xe = -0.904117256370474856682;
            *we = 0.106939325995318430960;
            break;
        case (9):
            *xe = 0.904117256370474856682;
            *we = 0.106939325995318430960;
            break;
        case (10):
            *xe = -0.981560634246719250712;
            *we = 0.0471753363865118271952;
            break;
        case (11):
            *xe = 0.981560634246719250712;
            *we = 0.0471753363865118271952;
            break;
        }
        break;


    case (16):
        switch (i)
        {
        case (0):
            *xe = -0.0950125098376374401877;
            *we = 0.189450610455068496287;
            break;
        case (1):
            *xe = 0.0950125098376374401877;
            *we = 0.189450610455068496287;
            break;
        case (2):
            *xe = -0.281603550779258913231;
            *we = 0.182603415044923588872;
            break;
        case (3):
            *xe = 0.281603550779258913231;
            *we = 0.182603415044923588872;
            break;
        case (4):
            *xe = -0.458016777657227386350;
            *we = 0.169156519395002538183;
            break;
        case (5):
            *xe = 0.458016777657227386350;
            *we = 0.169156519395002538183;
            break;
        case (6):
            *xe = -0.617876244402643748452;
            *we = 0.149595988816576732080;
            break;
        case (7):
            *xe = 0.617876244402643748452;
            *we = 0.149595988816576732080;
            break;
        case (8):
            *xe = -0.755404408355003033891;
            *we = 0.124628971255533872056;
            break;
        case (9):
            *xe = 0.755404408355003033891;
            *we = 0.124628971255533872056;
            break;
        case (10):
            *xe = -0.865631202387831743866;
            *we = 0.0951585116824927848073;
            break;
        case (11):
            *xe = 0.865631202387831743866;
            *we = 0.0951585116824927848073;
            break;
        case (12):
            *xe = -0.944575023073232576090;
            *we = 0.0622535239386478928628;
            break;
        case (13):
            *xe = 0.944575023073232576090;
            *we = 0.0622535239386478928628;
            break;
        case (14):
            *xe = -0.989400934991649932601;
            *we = 0.0271524594117540948514;
            break;
        case (15):
            *xe = 0.989400934991649932601;
            *we = 0.0271524594117540948514;
            break;
        }
        break;


    case (20):
        switch (i)
        {
        case (0):
            *xe = -0.0765265211334973337513;
            *we = 0.152753387130725850699;
            break;
        case (1):
            *xe = 0.0765265211334973337513;
            *we = 0.152753387130725850699;
            break;
        case (2):
            *xe = -0.227785851141645078076;
            *we = 0.149172986472603746785;
            break;
        case (3):
            *xe = 0.227785851141645078076;
            *we = 0.149172986472603746785;
            break;
        case (4):
            *xe = -0.373706088715419560662;
            *we = 0.142096109318382051326;
            break;
        case (5):
            *xe = 0.373706088715419560662;
            *we = 0.142096109318382051326;
            break;
        case (6):
            *xe = -0.510867001950827097985;
            *we = 0.131688638449176626902;
            break;
        case (7):
            *xe = 0.510867001950827097985;
            *we = 0.131688638449176626902;
            break;
        case (8):
            *xe = -0.636053680726515025467;
            *we = 0.118194531961518417310;
            break;
        case (9):
            *xe = 0.636053680726515025467;
            *we = 0.118194531961518417310;
            break;
        case (10):
            *xe = -0.746331906460150792634;
            *we = 0.101930119817240435039;
            break;
        case (11):
            *xe = 0.746331906460150792634;
            *we = 0.101930119817240435039;
            break;
        case (12):
            *xe = -0.839116971822218823420;
            *we = 0.0832767415767047487264;
            break;
        case (13):
            *xe = 0.839116971822218823420;
            *we = 0.0832767415767047487264;
            break;
        case (14):
            *xe = -0.912234428251325905857;
            *we = 0.0626720483341090635663;
            break;
        case (15):
            *xe = 0.912234428251325905857;
            *we = 0.0626720483341090635663;
            break;
        case (16):
            *xe = -0.963971927277913791287;
            *we = 0.0406014298003869413320;
            break;
        case (17):
            *xe = 0.963971927277913791287;
            *we = 0.0406014298003869413320;
            break;
        case (18):
            *xe = -0.993128599185094924776;
            *we = 0.0176140071391521183115;
            break;
        case (19):
            *xe = 0.993128599185094924776;
            *we = 0.0176140071391521183115;
            break;
        }
        break;


    case (24):
        switch (i)
        {
        case (0):
            *xe = -0.0640568928626056260827;
            *we = 0.127938195346752156976;
            break;
        case (1):
            *xe = 0.0640568928626056260827;
            *we = 0.127938195346752156976;
            break;
        case (2):
            *xe = -0.191118867473616309153;
            *we = 0.125837456346828296117;
            break;
        case (3):
            *xe = 0.191118867473616309153;
            *we = 0.125837456346828296117;
            break;
        case (4):
            *xe = -0.315042679696163374398;
            *we = 0.121670472927803391202;
            break;
        case (5):
            *xe = 0.315042679696163374398;
            *we = 0.121670472927803391202;
            break;
        case (6):
            *xe = -0.433793507626045138478;
            *we = 0.115505668053725601353;
            break;
        case (7):
            *xe = 0.433793507626045138478;
            *we = 0.115505668053725601353;
            break;
        case (8):
            *xe = -0.545421471388839535649;
            *we = 0.107444270115965634785;
            break;
        case (9):
            *xe = 0.545421471388839535649;
            *we = 0.107444270115965634785;
            break;
        case (10):
            *xe = -0.648093651936975569268;
            *we = 0.0976186521041138882720;
            break;
        case (11):
            *xe = 0.648093651936975569268;
            *we = 0.0976186521041138882720;
            break;
        case (12):
            *xe = -0.740124191578554364260;
            *we = 0.0861901615319532759152;
            break;
        case (13):
            *xe = 0.740124191578554364260;
            *we = 0.0861901615319532759152;
            break;
        case (14):
            *xe = -0.820001985973902921981;
            *we = 0.0733464814110803057346;
            break;
        case (15):
            *xe = 0.820001985973902921981;
            *we = 0.0733464814110803057346;
            break;
        case (16):
            *xe = -0.886415527004401034190;
            *we = 0.0592985849154367807461;
            break;
        case (17):
            *xe = 0.886415527004401034190;
            *we = 0.0592985849154367807461;
            break;
        case (18):
            *xe = -0.938274552002732758539;
            *we = 0.0442774388174198061695;
            break;
        case (19):
            *xe = 0.938274552002732758539;
            *we = 0.0442774388174198061695;
            break;
        case (20):
            *xe = -0.974728555971309498199;
            *we = 0.0285313886289336631809;
            break;
        case (21):
            *xe = 0.974728555971309498199;
            *we = 0.0285313886289336631809;
            break;
        case (22):
            *xe = -0.995187219997021360195;
            *we = 0.0123412297999871995469;
            break;
        case (23):
            *xe = 0.995187219997021360195;
            *we = 0.0123412297999871995469;
            break;
        }
        break;


    case (32):
        switch (i)
        {
        case (0):
            *xe = -0.0640568928626056260827;
            *we = 0.0965400885147278005666;
            break;
        case (1):
            *xe = 0.0640568928626056260827;
            *we = 0.0965400885147278005666;
            break;
        case (2):
            *xe = -0.144471961582796493484;
            *we = 0.0956387200792748594185;
            break;
        case (3):
            *xe = 0.144471961582796493484;
            *we = 0.0956387200792748594185;
            break;
        case (4):
            *xe = -0.239287362252137074544;
            *we = 0.0938443990808045656367;
            break;
        case (5):
            *xe = 0.239287362252137074544;
            *we = 0.0938443990808045656367;
            break;
        case (6):
            *xe = -0.331868602282127649782;
            *we = 0.0911738786957638847129;
            break;
        case (7):
            *xe = 0.331868602282127649782;
            *we = 0.0911738786957638847129;
            break;
        case (8):
            *xe = -0.421351276130635345353;
            *we = 0.0876520930044038111450;
            break;
        case (9):
            *xe = 0.421351276130635345353;
            *we = 0.0876520930044038111450;
            break;
        case (10):
            *xe = -0.506899908932229390044;
            *we = 0.0833119242269467552223;
            break;
        case (11):
            *xe = 0.506899908932229390044;
            *we = 0.0833119242269467552223;
            break;
        case (12):
            *xe = -0.587715757240762329066;
            *we = 0.0781938957870703064685;
            break;
        case (13):
            *xe = 0.587715757240762329066;
            *we = 0.0781938957870703064685;
            break;
        case (14):
            *xe = -0.663044266930215200960;
            *we = 0.0723457941088485062287;
            break;
        case (15):
            *xe = 0.663044266930215200960;
            *we = 0.0723464814110803057346;
            break;
        case (16):
            *xe = -0.732182118740289680412;
            *we = 0.0658222227763618468406;
            break;
        case (17):
            *xe = 0.732182118740289680412;
            *we = 0.0658222227763618468406;
            break;
        case (18):
            *xe = -0.794483795967942406965;
            *we = 0.0586840934785355471448;
            break;
        case (19):
            *xe = 0.794483795967942406965;
            *we = 0.0586840934785355471448;
            break;
        case (20):
            *xe = -0.849367613732569970160;
            *we = 0.0509980592623761761959;
            break;
        case (21):
            *xe = 0.849367613732569970160;
            *we = 0.0509980592623761761959;
            break;
        case (22):
            *xe = -0.896321155766052123971;
            *we = 0.0428358980222266806557;
            break;
        case (23):
            *xe = 0.896321155766052123971;
            *we = 0.0428358980222266806557;
            break;
        case (24):
            *xe = -0.934906075937739689159;
            *we = 0.0342738629130214331033;
            break;
        case (25):
            *xe = 0.934906075937739689159;
            *we = 0.0342738629130214331033;
            break;
        case (26):
            *xe = -0.964762255587506430761;
            *we = 0.0253920653092620594561;
            break;
        case (27):
            *xe = 0.964762255587506430761;
            *we = 0.0253920653092620594561;
            break;
        case (28):
            *xe = -0.985611511545268335400;
            *we = 0.0162743947309056706058;
            break;
        case (29):
            *xe = 0.985611511545268335400;
            *we = 0.0162743947309056706058;
            break;
        case (30):
            *xe = -0.997263861849481563534;
            *we = 0.00701861000947009660028;
            break;
        case (31):
            *xe = 0.997263861849481563534;
            *we = 0.00701861000947009660028;
            break;
        }
        break;


    case (48):
        switch (i)
        {
        case (0):
            *xe = -0.0323801709628693620343;
            *we = 0.0647376968126839225006;
            break;
        case (1):
            *xe = 0.0323801709628693620343;
            *we = 0.0647376968126839225006;
            break;
        case (2):
            *xe = -0.0970046992094626989322;
            *we = 0.0644661644359500822082;
            break;
        case (3):
            *xe = 0.0970046992094626989322;
            *we = 0.0644661644359500822082;
            break;
        case (4):
            *xe = -0.161222356068891718055;
            *we = 0.0639242385846481866207;
            break;
        case (5):
            *xe = 0.161222356068891718055;
            *we = 0.0639242385846481866207;
            break;
        case (6):
            *xe = -0.224763790394689061224;
            *we = 0.0631141922862540256548;
            break;
        case (7):
            *xe = 0.224763790394689061224;
            *we = 0.0631141922862540256548;
            break;
        case (8):
            *xe = -0.287362487355455576728;
            *we = 0.0620394231598926639029;
            break;
        case (9):
            *xe = 0.287362487355455576728;
            *we = 0.0620394231598926639029;
            break;
        case (10):
            *xe = -0.348755886292160738148;
            *we = 0.0607044391658938800517;
            break;
        case (11):
            *xe = 0.348755886292160738148;
            *we = 0.0607044391658938800517;
            break;
        case (12):
            *xe = -0.408686481990716729925;
            *we = 0.0591148396983956357477;
            break;
        case (13):
            *xe = 0.408686481990716729925;
            *we = 0.0591148396983956357477;
            break;
        case (14):
            *xe = -0.466902904750958404535;
            *we = 0.0572772921004032157044;
            break;
        case (15):
            *xe = 0.466902904750958404535;
            *we = 0.0572772921004032157044;
            break;
        case (16):
            *xe = -0.523160974722233033658;
            *we = 0.0551995036999841628676;
            break;
        case (17):
            *xe = 0.523160974722233033658;
            *we = 0.0551995036999841628676;
            break;
        case (18):
            *xe = -0.577224726083972703838;
            *we = 0.0528901894851936670964;
            break;
        case (19):
            *xe = 0.577224726083972703838;
            *we = 0.0528901894851936670964;
            break;
        case (20):
            *xe = -0.628867396776513624013;
            *we = 0.0503590355538544749590;
            break;
        case (21):
            *xe = 0.628867396776513624013;
            *we = 0.0503590355538544749590;
            break;
        case (22):
            *xe = -0.677872379632663905208;
            *we = 0.0476166584924904748267;
            break;
        case (23):
            *xe = 0.677872379632663905208;
            *we = 0.0476166584924904748267;
            break;
        case (24):
            *xe = -0.724034130923814654658;
            *we = 0.0446745608566942804201;
            break;
        case (25):
            *xe = 0.724034130923814654658;
            *we = 0.0446745608566942804201;;
            break;
        case (26):
            *xe = -0.767159032515740339276;
            *we = 0.0415450829434647492133;
            break;
        case (27):
            *xe = 0.767159032515740339276;
            *we = 0.0415450829434647492133;
            break;
        case (28):
            *xe = -0.807066204029442627087;
            *we = 0.0382413510658307063158;
            break;
        case (29):
            *xe = 0.807066204029442627087;
            *we = 0.0382413510658307063158;
            break;
        case (30):
            *xe = -0.843588261624393530704;
            *we = 0.0347772225647704388909;
            break;
        case (31):
            *xe = 0.843588261624393530704;
            *we = 0.0347772225647704388909;
            break;
        case (32):
            *xe = -0.876572020274247885885;
            *we = 0.0311672278327980889025;
            break;
        case (33):
            *xe = 0.876572020274247885885;
            *we = 0.0311672278327980889025;
            break;
        case (34):
            *xe = -0.905879136715569672805;
            *we = 0.0274265097083569482001;
            break;
        case (35):
            *xe = 0.905879136715569672805;
            *we = 0.0274265097083569482001;
            break;
        case (36):
            *xe = -0.931386690706554333107;
            *we = 0.0235707608393243791410;
            break;
        case (37):
            *xe = 0.931386690706554333107;
            *we = 0.0235707608393243791410;
            break;
        case (38):
            *xe = -0.952987703160430860724;
            *we = 0.0196161604573555278142;
            break;
        case (39):
            *xe = 0.952987703160430860724;
            *we = 0.0196161604573555278142;
            break;
        case (40):
            *xe = -0.970591592546247250472;
            *we = 0.0155793157229438487279;
            break;
        case (41):
            *xe = 0.970591592546247250472;
            *we = 0.0155793157229438487279;
            break;
        case (42):
            *xe = -0.984124583722826857765;
            *we = 0.0114772345792345394895;
            break;
        case (43):
            *xe = 0.984124583722826857765;
            *we = 0.0114772345792345394895;
            break;
        case (44):
            *xe = -0.993530172266350757526;
            *we = 0.00732755390127626210220;
            break;
        case (45):
            *xe = 0.993530172266350757526;
            *we = 0.00732755390127626210220;
            break;
        case (46):
            *xe = -0.998771007252426118580;
            *we = 0.00315334605230583863260;
            break;
        case (47):
            *xe = 0.998771007252426118580;
            *we = 0.00315334605230583863260;
            break;
        }
        break;
    }

}



