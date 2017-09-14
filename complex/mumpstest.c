#include "type.h"


int main(int argc, char *argv[])
{
    int myid, np,i;
    MPI_Init(&argc, &argv);   
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &np);

    SparseMatrixLoc A_loc;
    C_Init_Aloc(&A_loc, argc, argv);

    //右端项赋值
    int N = A_loc.N;
    if(myid == 0)
    {
        mumps_double_complex *rhs = calloc(N, sizeof(mumps_double_complex));
        for(i = 0; i < N; i++)
        {
            rhs[i].r = 1.0;
            rhs[i].i = 1.0;
        }
        A_loc.rhs = rhs;

    }
    
    
    //-----------------------------//
    ZMUMPS_STRUC_C id;
    Mumps_Init(&id);
    Mumps_Analize(&id, A_loc, myid);
    Mumps_FactorSolve(&id);
    Mumps_End(&id);


    if(myid == 0)
    {
        for(i = 0; i < 100; i++)
        {
            printf("rhs[%d].r is %lf\trhs[%d].i is %lf\n",i, rhs[i].r,i,rhs[i].i);
        }
    }


    MPI_Finalize();
    return 0;


}

   
/*********************************************
读入矩阵信息
***********************************************/
void C_Init_Aloc(SparseMatrixLoc *A_loc, int argc, char *argv[])
{
    int i, j, mpitag;

    //计算通信算子
    int myid, np;
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &np);

    /*******************************************************************************
                    ia  ja  a N nnz
    ********************************************************************************/
    int *ia, *ja, N, nnz;
    mumps_double_complex *a;
    if (myid == 0)
    {
        char *name = argv[1];
        FILE *fp;
        fp = fopen(name, "r");
        fscanf(fp, "%d%d%d\n", &N, &N, &nnz);

        a = calloc(nnz, sizeof(mumps_double_complex));
        ia = calloc(N + 1, sizeof(int));
        ja = calloc(nnz, sizeof(int));

        for (i = 0; i < nnz; i++)
        {
            fscanf(fp, "%d%d%lf%lf\n", &j, ja + i, &(a[i].r), &(a[i].i));
            ia[j] = ia[j] + 1;
        }
        fclose(fp);

        ia[0] = 1;
        for (i = 0; i < N; i++)
            ia[i + 1] = ia[i + 1] + ia[i];
    }

    MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);


    /*******************************************************************************
                  myN_loc  mynnz_loc  myfstrow
                  N_loc    nnz_loc    fstrow
    ********************************************************************************/
    int myN_loc, mynnz_loc, myfstrow, *N_loc, *nnz_loc, *fstrow;
    N_loc = calloc(np, sizeof(int));
    nnz_loc = calloc(np, sizeof(int));
    fstrow = calloc(np + 1, sizeof(int));
    for (i = 0; i < np; i++)
    {
        if (i < N % np)
        {
            N_loc[i] = N / np + 1;
            fstrow[i] = N_loc[i] * i;
        }
        else
        {
            N_loc[i] = N / np;
            fstrow[i] = N_loc[i] * i + N % np;
        }
    }
    fstrow[np] = fstrow[np - 1] + N_loc[np - 1];

    if (fstrow[np] != N)
        printf("error in fstrow\n");

    if (myid == 0)
    {
        int begin, end;
        for (i = 0; i < np; i++)
        {
            begin = fstrow[i];
            end = begin + N_loc[i];
            nnz_loc[i] = ia[end] - ia[begin];
        }

    }
    MPI_Bcast(nnz_loc, np, MPI_INT, 0, MPI_COMM_WORLD);

    myN_loc = N_loc[myid];
    mynnz_loc = nnz_loc[myid];
    myfstrow = fstrow[myid];
    printf("myid %d myN_loc %d mynnz_loc %d myfstrow %d \n", myid, myN_loc, mynnz_loc, myfstrow);


    /*******************************************************************************
                     ia_loc  ja_loc  a_loc
    ********************************************************************************/
    int *ia_loc, *ja_loc;
    mumps_double_complex *a_loc;

    ia_loc = calloc(myN_loc + 1, sizeof(int));
    ja_loc = calloc(mynnz_loc,   sizeof(int));
    a_loc  = calloc(mynnz_loc,   sizeof(mumps_double_complex));

    MPI_Datatype MPI_complex;
    MPI_Type_contiguous(2, MPI_DOUBLE, &MPI_complex);
    MPI_Type_commit(&MPI_complex);

    if (myid == 0)
    {
        for (i = 0; i < myN_loc + 1; i++)
            ia_loc[i] = ia[i];

        for (i = 0; i < mynnz_loc; i++)
        {
            ja_loc[i] = ja[i];
            a_loc[i] = a[i];
        }
        int offset = mynnz_loc;

        for (i = 1; i < np; i++)
        {
            mpitag = 1;
            MPI_Send(ia + fstrow[i], N_loc[i] + 1, MPI_INT,     i, mpitag, MPI_COMM_WORLD);
            mpitag = 2;
            MPI_Send(ja + offset,    nnz_loc[i],   MPI_INT,     i, mpitag, MPI_COMM_WORLD);
            mpitag = 3;
            MPI_Send(a + offset,     nnz_loc[i],   MPI_complex, i, mpitag, MPI_COMM_WORLD);

            offset = offset + nnz_loc[i];

        }
    }
    else
    {
        ia_loc = calloc(N_loc[myid] + 1, sizeof(int));
        ja_loc = calloc(nnz_loc[myid],   sizeof(int));
        a_loc  = calloc(nnz_loc[myid],   sizeof(mumps_double_complex));

        MPI_Status status;
        mpitag = 1;
        MPI_Recv(ia_loc, myN_loc + 1, MPI_INT,     0, mpitag, MPI_COMM_WORLD, &status);
        mpitag = 2;
        MPI_Recv(ja_loc, mynnz_loc,   MPI_INT,     0, mpitag, MPI_COMM_WORLD, &status);
        mpitag = 3;
        MPI_Recv(a_loc,  mynnz_loc,   MPI_complex, 0, mpitag, MPI_COMM_WORLD, &status);

        for (i = 1; i <= N_loc[myid]; i++)
        {
            ia_loc[i] = ia_loc[i] - ia_loc[0] + 1;
        }
        ia_loc[0] = 1;

    }
    if (myid == 0)
    {
        free(a);
        free(ia);
        free(ja);
    }

//    (*info).myid = myid;
//    (*info).np = np;
//    (*info).N_loc = N_loc;
//    (*info).fstrow = fstrow;
//    (*info).N = N;

    (*A_loc).N_loc = myN_loc;
    (*A_loc).nnz_loc = mynnz_loc;
    (*A_loc).fstrow = myfstrow;

    (*A_loc).ia_loc = ia_loc;
    (*A_loc).ja_loc = ja_loc;
    (*A_loc).a_loc = a_loc;

}


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
void Mumps_Analize(ZMUMPS_STRUC_C *id, SparseMatrixLoc C_loc, int myid)
{
    (*id).icntl[2] = -1;
    (*id).icntl[3] = 1;

    (*id).icntl[4] = 0;
    (*id).icntl[17] = 3;

    if (myid == 0)
    {
        (*id).rhs = C_loc.rhs;    //右端项数组（一维或者二维)
        (*id).nrhs = 1;           //the number of right hand vectors
        (*id).lrhs = C_loc.N;     //vector的长度
        (*id).n = C_loc.N;        //vector的长度

    }

    (*id).nz_loc = C_loc.nnz_loc;   //number of nozero-value
    (*id).irn_loc = C_loc.ia_loc;   //row index
    (*id).jcn_loc = C_loc.ja_loc;   //col index
    (*id).a_loc = C_loc.a_loc;      //nozero-value
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
