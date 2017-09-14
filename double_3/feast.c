#include "type.h"
/************************************
生成随机矩阵
*****************************************/
void Feast_SetRandomMatrix(double *Q_loc, Info info)
{
    int myN_loc = info.N_loc[info.intervalid];
    int *N_loc = info.N_loc;
    int N = info.N;
    int M0 = info.M0;

    double *Q;

    int i, j, k;
    int length = N * M0;
    //printf("aaaaaaaaaaaaaaaaaaaaaaaaa  intervalid %d \n", info.intervalid);
    if (info.intervalid == 0)
    {
        Q = calloc(N * M0, sizeof(double));

        int three = 3;
        int iseed[4] = { 0, 3, 9, 1 };
        dlarnv_(&three, iseed, &length, Q);

        for (i = 0; i < M0; i++)
            for (j = 0; j < myN_loc; j++)
                Q_loc[i * myN_loc + j] = Q[i * N + j];

        int offset, nrow;
        offset = myN_loc;

        double *sendbuf;
        sendbuf = calloc(myN_loc * M0, sizeof(double));

        for (i = 1; i < info.intervalnp; i++)
        {
            // printf("%d \n", type_N);
            nrow = N_loc[i];
            for (j = 0; j < M0; j++)
                for (k = 0; k < nrow; k++)
                    sendbuf[j * nrow + k] = Q[offset + j * N + k];

            MPI_Send(sendbuf, nrow * M0, MPI_DOUBLE, i, i, info.intervalcomm);
            offset = offset + nrow;
        }
        free(sendbuf);
        free(Q);

    }
    else
    {
        MPI_Status status;
        MPI_Recv(Q_loc, myN_loc * M0, MPI_DOUBLE,  0, info.intervalid, info.intervalcomm, &status);
    }
    //printf("bbbbbbbbbbbbbbbbbbbbbbbbbbb  intervalid %d \n", info.intervalid);


}


/************************************
计算投影矩阵
*****************************************/
void Feast_GetProjectionOp(double *Y_loc, LinearSystem linearsystem, Info info)
{
    int myN_loc = info.N_loc[info.intervalid];
    int *N_loc = info.N_loc;
    int N = info.N;
    int M0 = info.M0;

    int i, j, k, nrow, mpitag;
    MPI_Status status;

    for (i = 0; i < myN_loc * M0; i++)
        Y_loc[i] = 0;

    double *WY;
    WY = calloc(myN_loc * M0, sizeof(double));

    mpitag = 0;
    if (info.mumpsid == 0)
    {
        mumps_double_complex jac, *rhs;
        jac = linearsystem.jac;
        rhs = linearsystem.rhs;

        double *sendbuf;
        sendbuf = calloc((myN_loc + 1) * M0, sizeof(double));
        int sendlength, offset = 0;

        for (i = 0; i < info.intervalnp; i++)
        {
            if (i < info.nGau)
            {
                if (i < info.intervalid)
                {
                    MPI_Recv(WY,  myN_loc * M0, MPI_DOUBLE, i, mpitag, info.intervalcomm, &status);
                    for (j = 0; j < myN_loc * M0; j++)
                        Y_loc[j] = Y_loc[j] + WY[j];

                    nrow = N_loc[i];
                    sendlength = nrow * M0;
                    for (j = 0; j < nrow; j++)
                        for (k = 0; k < M0; k++)
                            sendbuf[j + k * nrow] = jac.r * rhs[offset + j + k * N].r - jac.i * rhs[offset + j + k * N].i;
                    MPI_Send(sendbuf, sendlength, MPI_DOUBLE, i, mpitag, info.intervalcomm);
                    // if (i == 2)
                    // printf("offset  %d sendbuf %lf %lf %lf\n", offset, sendbuf[0], sendbuf[1], sendbuf[2]);
                    offset = offset + nrow;
                }
                else if (i > info.intervalid)
                {
                    nrow = N_loc[i];
                    sendlength = nrow * M0;
                    for (j = 0; j < nrow; j++)
                        for (k = 0; k < M0; k++)
                            sendbuf[j + k * nrow] = jac.r * rhs[offset + j + k * N].r - jac.i * rhs[offset + j + k * N].i;
                    MPI_Send(sendbuf, sendlength, MPI_DOUBLE, i, mpitag, info.intervalcomm);
                    offset = offset + nrow;

                    MPI_Recv(WY,  myN_loc * M0, MPI_DOUBLE, i, mpitag, info.intervalcomm, &status);
                    for (j = 0; j < myN_loc * M0; j++)
                        Y_loc[j] = Y_loc[j] + WY[j];
                }
                else if (i == info.intervalid)
                {
                    nrow = N_loc[i];
                    sendlength = nrow * M0;
                    for (j = 0; j < nrow; j++)
                        for (k = 0; k < M0; k++)
                            Y_loc[j + k * nrow] = Y_loc[j + k * nrow] + jac.r * rhs[offset + j + k * N].r - jac.i * rhs[offset + j + k * N].i;
                    offset = offset + nrow;
                }
            }
            else
            {
                nrow = N_loc[i];
                sendlength = nrow * M0;
                for (j = 0; j < nrow; j++)
                    for (k = 0; k < M0; k++)
                        sendbuf[j + k * nrow] = jac.r * rhs[offset + j + k * N].r - jac.i * rhs[offset + j + k * N].i;
                MPI_Send(sendbuf, sendlength, MPI_DOUBLE, i, mpitag, info.intervalcomm);
                // if (i == 2)
                //     printf("offset  %d sendbuf %lf %lf %lf\n", offset, sendbuf[0], sendbuf[1], sendbuf[2]);
                offset = offset + nrow;
            }
        }
        free(sendbuf);
        if (offset != N)
            printf("error in mul  offset\n");
    }
    else
    {
        for (i = 0; i < info.nGau; i++)
        {
            MPI_Recv(WY,  myN_loc * M0, MPI_DOUBLE, i, mpitag, info.intervalcomm, &status);

            for (j = 0; j < myN_loc * M0; j++)
                Y_loc[j] = Y_loc[j] + WY[j];
        }
    }
    free(WY);
}


/************************************
根据计算得到的投影矩阵进行投影,并且将投影矩阵正交化
*****************************************/
void mul(double *result, double *Y_loc, SPMultiply SPMul, Info info);
void Feast_Project(double *Q_loc, double *lambda, double *Y_loc, SPMultiply SPMul, Info info)
{
    double mytime1, mytime2, mytime3, mytime4;

    int M0 = info.M0;
    int i, j;
    int N_loc = info.N_loc[info.intervalid];


    char char_N = 'N', char_T = 'T', char_V = 'V', char_U = 'U';
    double one = 1, zero = 0;
    int int_one = 1;
    char matdescra[6] = {'G', 'U', 'U', 'F'}; //矩阵A的类型

    int error, lwork = N_loc * M0;
    double *work_loc, *Aq_loc, *Aq, *Bq_loc, *Bq;
    work_loc = calloc(N_loc * M0, sizeof(double));
    Aq_loc = calloc(M0 * M0, sizeof(double));
    Bq_loc = calloc(M0 * M0, sizeof(double));

    Aq = calloc(M0 * M0, sizeof(double));
    Bq = calloc(M0 * M0, sizeof(double));

    //Aq=Y^t A Y
    mytime1 = MPI_Wtime();
    mul(work_loc, Y_loc, SPMul, info);
    //if (info.intervalid == 0)
    // printf("mul_loc  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf \n", work_loc[0], work_loc[1], work_loc[2], work_loc[3], work_loc[4], work_loc[5], work_loc[6], work_loc[7], work_loc[8], work_loc[9]);

    mytime4 = MPI_Wtime();
    dgemm_(&char_T, &char_N, &M0, &M0, &N_loc, &one, Y_loc, &N_loc, work_loc, &N_loc, &zero, Aq_loc, &M0); //Aq = Y^T * work
    mytime3 = MPI_Wtime();
    MPI_Reduce(Aq_loc, Aq, M0 * M0, MPI_DOUBLE, MPI_SUM, 0, info.intervalcomm);
    mytime2 = MPI_Wtime();
    if (info.myid == 0)
        printf("\t\t Aq\ttime %lf  time  %lf %lf %lf \n", mytime2 - mytime1, mytime4 - mytime1, mytime3 - mytime4, mytime2 - mytime3);


    //mytime1 = MPI_Wtime();
    dgemm_(&char_T, &char_N, &M0, &M0, &N_loc, &one, Y_loc, &N_loc, Y_loc, &N_loc, &zero, Bq_loc, &M0);//Bq=Y^T * Y
    MPI_Reduce(Bq_loc, Bq, M0 * M0, MPI_DOUBLE, MPI_SUM, 0, info.intervalcomm);
    mytime2 = MPI_Wtime();
    if (info.myid == 0)
        printf("\t\t ~~~~~Bq\ttime %lf \n", mytime2 - mytime1);

    if (info.intervalid == 0)
    {
        mytime1 = MPI_Wtime();
        dsygv_(&int_one, &char_V, &char_U, &M0, Aq, &M0, Bq, &M0, lambda, work_loc, &lwork, &error);//Aq * Aq' = lambda * Bq * Aq'
        if (error != 0)
            printf("error in dsygv_ error %d\n", error);
        mytime2 = MPI_Wtime();
        if (info.myid == 0)
            printf("\t\t\t GV\ttime %lf \n", mytime2 - mytime1);
    }

    mytime1 = MPI_Wtime();
    MPI_Bcast(Aq, M0 * M0, MPI_DOUBLE, 0, info.intervalcomm);
    //  printf("Aq %lf %lf %lf %lf %lf %lf\n", Aq[0], Aq[1], Aq[2], Aq[3], Aq[4], Aq[5]);


    dgemm_(&char_N, &char_N, &N_loc, &M0, &M0, &one, Y_loc, &N_loc, Aq, &M0, &zero, Q_loc, &N_loc); //Q = Y * Aq
    //  printf("Q_loc %lf %lf %lf %lf %lf %lf\n", Q_loc[0], Q_loc[1], Q_loc[2], Q_loc[3], Q_loc[4], Q_loc[5]);
    mytime2 = MPI_Wtime();
    if (info.myid == 0)
        printf("\t\t QQQQQ\ttime %lf \n", mytime2 - mytime1);

    free(work_loc);
    free(Aq_loc);
    free(Bq_loc);

    free(Aq);
    free(Bq);


}


/************************************
计算求解的特征值的误差
*****************************************/
void Feast_GetRes(double *maxres, double *lambda, double *Q_loc, SPMultiply SPMul, Info info)
{

    int N_loc = info.N_loc[info.intervalid];
    int M0 = info.M0;
    double Emax = info.Emax;
    double Emin = info.Emin;

    int i, j, k;

    double *work_loc, *res_loc, *res, *x_loc, *x;
    work_loc = calloc(N_loc * M0, sizeof(double));
    res_loc = calloc(M0, sizeof(double));
    x_loc = calloc(M0, sizeof(double));
    res = calloc(M0, sizeof(double));
    x = calloc(M0, sizeof(double));

    mul(work_loc, Q_loc, SPMul, info);
    MPI_Bcast(lambda, M0, MPI_DOUBLE, 0, info.intervalcomm);

    for (i = 0; i < M0; i++)
        for (j = 0; j < N_loc; j++)
            work_loc[i * N_loc + j] =  work_loc[i * N_loc + j] - lambda[i] * Q_loc[i * N_loc + j];

    for (i = 0; i < M0; i++)
    {
        res_loc[i] = 0;
        x_loc[i] = 0;
        for (j = 0; j < N_loc; j++)
        {
            res_loc[i] = res_loc[i] + fabs(work_loc[i * N_loc + j]);
            x_loc[i] = x_loc[i] + fabs(Q_loc[i * N_loc + j]);
        }
    }

    MPI_Reduce(res_loc, res, M0, MPI_DOUBLE, MPI_SUM, 0, info.intervalcomm);
    MPI_Reduce(x_loc,   x,   M0, MPI_DOUBLE, MPI_SUM, 0, info.intervalcomm);


    if (info.intervalid == 0)
    {
        for (i = 0; i < M0; i++)
            res[i] = res[i] / x[i];

        int neig = 0;
        int maxid;
        *maxres = 0;
        for (i = 0; i < M0; i++)
        {
            if (lambda[i] < Emax && lambda[i] > Emin)
            {
                neig++;
                if (res[i] > *maxres)
                {
                    maxid = i;
                    *maxres = res[i];
                }
            }
        }
        printf("\t\t maxres: %lf  neig %d %d \n", *maxres, neig, info.neig[info.myid % info.ninterval]);

        if (neig > info.neig[info.myid % info.ninterval])
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
        }

    }
    free(work_loc);
    free(res);
    free(x);
    free(res_loc);
    free(x_loc);


    MPI_Bcast(maxres, 1, MPI_DOUBLE, 0, info.intervalcomm);
}



/************************************
获取估计的特征值数目
*****************************************/
int  Feast_GetNumber(double *Y_loc, Info info)
{

    int int_one = 1;
    double one = 1, zero = 0;
    char char_T = 'T', char_N = 'N', char_U = 'U', char_V = 'V';

    int i, j, error;
    int N = info.N;
    int N_loc = info.N_loc[info.intervalid];
    int M0 = info.M0;
    int length = M0 * M0;

    double *Bq, *Bq_loc, *work, *lambda_B;
    Bq = calloc(M0 * M0, sizeof(double));
    Bq_loc = calloc(M0 * M0, sizeof(double));
    work = calloc(M0 * M0, sizeof(double));
    lambda_B = calloc(M0, sizeof(double));

    // printf("myid %d Y_loc  %lf %lf %lf \n", info.myid, Y_loc[0], Y_loc[4], Y_loc[6]);

    dgemm_(&char_T, &char_N, &M0, &M0, &N_loc, &one, Y_loc, &N_loc, Y_loc, &N_loc, &zero, Bq_loc, &M0);//Bq=Y^T * Y

    MPI_Reduce(Bq_loc, Bq, M0 * M0, MPI_DOUBLE, MPI_SUM, 0, info.intervalcomm);

    int neig;
    if (info.intervalid == 0)
    {
        for (i = 0; i < M0 * M0; i++)
            Bq_loc[i] = Bq[i];

        dsyev_(&char_N, &char_U, &M0, Bq_loc, &M0, lambda_B, work, &length, &error);//Bq * Bq' = lambda * Bq'
        if (error != 0)
            printf("error in dsyev_ error %d Mo  %d\n", error, M0);

        neig = 0;
        for (i = 0; i < M0; i++)
        {
            if (lambda_B[i] > 0.25)
                neig = neig + 1;
        }
    }


    MPI_Bcast(&neig, 1, MPI_INT, 0, info.intervalcomm);
    free(Bq);
    free(Bq_loc);
    free(work);
    free(lambda_B);

    return neig;


}