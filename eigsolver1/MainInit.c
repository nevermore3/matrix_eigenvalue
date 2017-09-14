#include <stdio.h>
#include "type.h"

void InitInfo(Info *info);
void InitData(SparseMatrix *A, Info *info, char *argv[]);
void InitMul_System(SPMultiply *SPMul, LinearSystem *LSystem, SparseMatrix A, Info *info);
/*********************************************
               初始化
***********************************************/
void Init(SPMultiply *SPMul, LinearSystem *LSystem, Info *info, int argc, char *argv[])
{
    InitInfo(info);

    SparseMatrix A;

    InitData(&A, info, argv);

    InitMul_System(SPMul, LSystem, A, info);

    if ((*info).myid == 0)
    {
        free(A.ia);
        free(A.ja);
        free(A.a);
    }
}


/*********************************************
           初始化info的信息
***********************************************/
void InitInfo(Info *info)
{
    int nGau = 4;
    (*info).nGau = nGau;
    (*info).interval = calloc(2, sizeof(double));

    //计算通信算子
    int myid, np, color;
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &np);

/*
    MPI_Comm mumpscomm;
    int mumpsid, mumpsnp;
    //why----------------------------//
    color = myid % (*info).nGau ;
    MPI_Comm_split(MPI_COMM_WORLD, color, myid, &mumpscomm);
    MPI_Comm_rank(mumpscomm, &mumpsid);
    MPI_Comm_size(mumpscomm, &mumpsnp);

    MPI_Comm hostcomm;
    int hostid, hostnp;
    MPI_Comm_split(MPI_COMM_WORLD, mumpsid, myid, &hostcomm);
    MPI_Comm_rank(hostcomm, &hostid);
    MPI_Comm_size(hostcomm, &hostnp);
*/

    (*info).myid = myid;
    (*info).np = np;
/*
    (*info).mumpscomm = mumpscomm;
    (*info).mumpsid = mumpsid;
    (*info).mumpsnp = mumpsnp;

    (*info).hostcomm = hostcomm;
    (*info).hostid = hostid;
    (*info).hostnp = hostnp;
*/
}


/*********************************************
读入矩阵信息
***********************************************/
void InitData(SparseMatrix *A, Info *info, char *argv[])
{
    int i, j, mpitag = 1;
    int myid = (*info).myid;

    /***************************** *a  *ia  *ja  N  nnz*********************************/
    int *ia, *ja,  N, nnz;
    double *a;
    if (myid == 0)
    {
        char *name = argv[1];
        FILE *fp;
        fp = fopen(name, "r");
        fscanf(fp, "%d%d%d\n", &N, &N, &nnz);
        a = calloc(nnz, sizeof(double));
        ia = calloc(N + 1, sizeof(int));
        ja = calloc(nnz, sizeof(int));

        for (i = 0; i < nnz; i++)
        {
            fscanf(fp, "%d%d%lf\n", &j, ja + i, a + i);
            ia[j] = ia[j] + 1;
        }
        fclose(fp);

        ia[0] = 1;
        for (i = 0; i < N; i++)
            ia[i + 1] = ia[i + 1] + ia[i];
    }
    MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);

    /***************************** ceig *********************************/
    (*A).N = N;
    (*A).nnz = nnz;
    (*A).ia = ia;
    (*A).ja = ja;
    (*A).a = a ;

    (*info).N = N;
}


/*********************************************
              初始化数据
***********************************************/
void InitMatrixLoc(SparseMatrixLoc *A_loc, SparseMatrix A, Info *info);
void InitMul(SPMultiply *SPMul, SparseMatrixLoc A_loc, Info *info);
void InitSystem(LinearSystem *LSystem, SparseMatrixLoc A_loc, Info *info);
void InitMul_System(SPMultiply *SPMul, LinearSystem *LSystem, SparseMatrix A, Info *info)
{

    SparseMatrixLoc A_loc;

    InitMatrixLoc(&A_loc, A, info);

    InitMul(SPMul, A_loc, info);

    InitSystem(LSystem, A_loc, info);

    free(A_loc.ia_loc);
    free(A_loc.ja_loc);
    free(A_loc.a_loc);
    free(A_loc.nnz_loc);
    free(A_loc.fstrow);

}

/*********************************************
         数据分布存储
***********************************************/
void InitMatrixLoc(SparseMatrixLoc *A_loc, SparseMatrix A, Info *info)
{
    int np = (*info).np;
    int myid = (*info).myid;

    int *ia = A.ia;
    int *ja = A.ja;
    double *a = A.a;
    int N = A.N;
    int nnz = A.nnz;

    int i, j, k;

    /****************** *infoN_loc, *infonnz_loc, *infofstrow **************************/
    int *infoN_loc = calloc(np, sizeof(int));
    int *infonnz_loc = calloc(np, sizeof(int));
    int *infofstrow = calloc(np + 1, sizeof(int));

    for (i = 0; i < np; i++)
    {
        if (i < N % np)
        {
            infoN_loc[i] = N / np + 1;
            infofstrow[i] = infoN_loc[i] * i;
        }
        else
        {
            infoN_loc[i] = N / np;
            infofstrow[i] = infoN_loc[i] * i + N % np;
        }
    }
    infofstrow[np] = infofstrow[np - 1] + infoN_loc[np - 1];
    if (myid == 0)
    {
        int begin, end;
        for (i = 0; i < np; i++)
        {
            begin = infofstrow[i];
            end = begin + infoN_loc[i];
            infonnz_loc[i] = ia[end] - ia[begin];
        }

    }
    MPI_Bcast(infonnz_loc, np, MPI_INT, 0, MPI_COMM_WORLD);

    /******************  N_loc, nnz_loc, fstrow **************************/
    int N_loc, nnz_loc, fstrow;
    N_loc = infoN_loc[myid];
    fstrow = infofstrow[myid];
    nnz_loc = infonnz_loc[myid];


    /******************   *ia_loc, *ja_loc, *a_loc  **************************/
    int *ia_loc, *ja_loc;
    double *a_loc;
    int mpitag;
    if (myid == 0)
    {
        ia_loc = calloc(N_loc + 1, sizeof(int));
        ja_loc = calloc(nnz_loc,   sizeof(int));
        a_loc  = calloc(nnz_loc,   sizeof(double));

        for (i = 0; i < N_loc + 1; i++)
            ia_loc[i] = ia[i];
        for (i = 0; i < nnz_loc; i++)
        {
            ja_loc[i] = ja[i];
            a_loc[i] = a[i];
        }
        int offset = nnz_loc;

        for (i = 1; i < np; i++)
        {
            mpitag = 1;
            MPI_Send(ia + infofstrow[i], infoN_loc[i] + 1, MPI_INT,    i, mpitag, MPI_COMM_WORLD);
            mpitag = 2;
            MPI_Send(ja + offset,        infonnz_loc[i],   MPI_INT,    i, mpitag, MPI_COMM_WORLD);
            mpitag = 3;
            MPI_Send(a + offset,         infonnz_loc[i],   MPI_DOUBLE, i, mpitag, MPI_COMM_WORLD);

            offset = offset + infonnz_loc[i];
        }
    }
    else
    {
        ia_loc = calloc(N_loc + 1, sizeof(int));
        ja_loc = calloc(nnz_loc,   sizeof(int));
        a_loc  = calloc(nnz_loc,   sizeof(double));

        MPI_Status status;
        mpitag = 1;
        MPI_Recv(ia_loc, N_loc + 1, MPI_INT,    0, mpitag, MPI_COMM_WORLD, &status);
        mpitag = 2;
        MPI_Recv(ja_loc, nnz_loc,   MPI_INT,    0, mpitag, MPI_COMM_WORLD, &status);
        mpitag = 3;
        MPI_Recv(a_loc,  nnz_loc,   MPI_DOUBLE, 0, mpitag, MPI_COMM_WORLD, &status);

        for (i = 1; i <= N_loc; i++)
        {
            ia_loc[i] = ia_loc[i] - ia_loc[0] + 1;
        }
        ia_loc[0] = 1;

    }

    (*A_loc).ia_loc = ia_loc;
    (*A_loc).ja_loc = ja_loc;
    (*A_loc).a_loc = a_loc;
    (*A_loc).N_loc = infoN_loc;
    (*A_loc).nnz_loc = infonnz_loc;
    (*A_loc).fstrow = infofstrow;
    (*info).N_loc = infoN_loc;


}

/*********************************************
            初始化乘法
***********************************************/
void InitMult_sortja(int N, int *ja, int *ia, double *a);
void InitMul(SPMultiply *SPMul, SparseMatrixLoc A_loc, Info *info)
{
    int np = (*info).np;
    int myid = (*info).myid;

    int myN_loc = A_loc.N_loc[myid];
    int mynnz_loc = A_loc.nnz_loc[myid];
    int myfstrow = A_loc.fstrow[myid];

    int *fstrow = A_loc.fstrow;
    int *ia_loc = A_loc.ia_loc;
    int *ja_loc = A_loc.ja_loc;
    double *a_loc = A_loc.a_loc;

    int i, j, begin, end, id, idin, idout, from, to;

    //统计内点和外点的个数
    /******************   nnz_in , nnz_out  **************************/
    int nnz_in , nnz_out;
    begin = myfstrow;
    end = myfstrow + myN_loc + 1;
    nnz_in = 0;
    nnz_out = 0;
    for (i = 0; i < mynnz_loc; i++)
    {
        if (ja_loc[i] > begin && ja_loc[i] < end)
            nnz_in++;
        else
            nnz_out++;
    }

    //构造行压缩存储的内点以及三元组存储的外点
    /******************   *ia_in, *ja_in, *ia_out, *ja_out *a_in, *a_out **************************/
    int *ia_in, *ja_in, *ia_out, *ja_out;
    double *a_in, *a_out;
    ia_in = calloc(myN_loc + 1, sizeof(int));
    ja_in = calloc(nnz_in,    sizeof(int));
    a_in =  calloc(nnz_in,    sizeof(double));

    if (nnz_out > 0)
    {
        ia_out = calloc(nnz_out, sizeof(int));
        ja_out = calloc(nnz_out, sizeof(int));
        a_out =  calloc(nnz_out, sizeof(double));
    }

    idin = 0;
    idout = 0;
    ia_in[0] = 1;
    for (i = 0; i < myN_loc; i++)
    {
        from = ia_loc[i] - 1;
        to = ia_loc[i + 1] - 1;

        for (j = from; j < to; j++)
        {
            if (ja_loc[j] > begin && ja_loc[j] < end)
            {
                ia_in[i + 1] = ia_in[i + 1] + 1;
                ja_in[idin] = ja_loc[j] - begin;
                a_in[idin] = a_loc[j];
                idin++;
            }
            else
            {
                ia_out[idout] = i;
                ja_out[idout] = ja_loc[j] - 1;
                a_out[idout] = a_loc[j];
                idout++;
            }
        }

        ia_in[i + 1] = ia_in[i + 1] + ia_in[i];
    }
    if (idin != nnz_in || idout != nnz_out)
        printf("!!!!!!!!!!!!error in calucation for  idin %d nnz_in: %d idout %d nnz_out: %d \n", idin, nnz_in, idout, nnz_out);

    InitMult_sortja( nnz_out, ja_out, ia_out, a_out);//外点按列排序


    //统计计算所需要的外部稠密阵信息————行数，行号以及所在进程号
    /******************   nrecv, *rowid_out, *recv_count **************************/
    int nrecv, *rowid_out, *recv_count;
    if (nnz_out == 0)
    {
        nrecv = 0;
        rowid_out = calloc(1, sizeof(int));
        recv_count = calloc(1, sizeof(int));
    }
    else
    {
        nrecv = 1;
        for (i = 1; i < nnz_out; i++)
            if (ja_out[i] > ja_out[i - 1])
                nrecv++;

        rowid_out = calloc(nrecv, sizeof(int));
        rowid_out[0] = ja_out[0];
        ja_out[0] = 0;
        id = 1;
        for (i = 1; i < nnz_out; i++)
        {
            if (ja_out[i] > rowid_out[id - 1])
            {
                rowid_out[id] = ja_out[i];
                id++;
                ja_out[i] = ja_out[i - 1] + 1;
            }
            else
            {
                ja_out[i] = ja_out[i - 1];
            }
        }
        if (id != nrecv)
            printf("error in caculate rowid_out \n");

        recv_count = calloc(np, sizeof(int));
        for (i = 0; i < nrecv; i++)
        {
            id = 0;
            while (rowid_out[i] >= fstrow[id + 1])
            {
                id++;
            }
            recv_count[id] = recv_count[id] + 1;
        }
        if (id > np)
            printf("error in caculate recv_count\n");

    }

    //记录需要发送给各个处理的行数
    /******************  *send_count **************************/
    int  *send_count;
    send_count = calloc(np, sizeof(int));
    int mpitag, sourceid, destid;
    MPI_Status status;
    for (i = 1; i < np; i++)
    {
        sourceid = (myid - i + np) % np;
        destid = (myid + i) % np;
        if (myid < np / 2 )
        {
            MPI_Recv(send_count + sourceid, 1, MPI_INT, sourceid, myid, MPI_COMM_WORLD, &status);
            MPI_Send(recv_count + destid,   1, MPI_INT, destid,   destid, MPI_COMM_WORLD);
        }

        else
        {
            MPI_Send(recv_count + destid,   1, MPI_INT, destid,   destid, MPI_COMM_WORLD);
            MPI_Recv(send_count + sourceid, 1, MPI_INT, sourceid, myid, MPI_COMM_WORLD, &status);
        }
    }


    //统计需要发送的次数
    /******************  nsend **************************/
    int nsend;
    nsend = 0;
    for (i = 0; i < np; i++)
        nsend = nsend + send_count[i];


    //记录需要发送的各个进程号以及行号
    /******************  *send_rowid **************************/
    int *send_rowid;
    if (nsend > 0)
        send_rowid = calloc(nsend, sizeof(int));
    else
        send_rowid = calloc(1, sizeof(int));


    int idsend, idrecv, sendoffset, recvoffset, number;
    idsend = 0;
    idrecv = 0;
    sendoffset = 0;
    recvoffset = 0;

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
                MPI_Recv(send_rowid + recvoffset, send_count[idrecv], MPI_INT, idrecv, myid, MPI_COMM_WORLD, &status);
                recvoffset = recvoffset + send_count[idrecv];
                idrecv++;
            }
            if (idsend > myid)
            {
                MPI_Send(rowid_out + sendoffset, recv_count[idsend], MPI_INT,   idsend, idsend, MPI_COMM_WORLD);
                sendoffset = sendoffset + recv_count[idsend];
                idsend++;
            }
        }
        else
        {
            if (idrecv < idsend || idsend == np)
            {
                MPI_Recv(send_rowid + recvoffset, send_count[idrecv], MPI_INT, idrecv, myid, MPI_COMM_WORLD, &status);
                recvoffset = recvoffset + send_count[idrecv];
                idrecv++;
            }
            else if (idsend < idrecv || idrecv == np)
            {
                MPI_Send(rowid_out + sendoffset, recv_count[idsend], MPI_INT,   idsend, idsend, MPI_COMM_WORLD);
                sendoffset = sendoffset + recv_count[idsend];
                idsend++;
            }
        }
        while (recv_count[idsend] == 0 && idsend < np)
            idsend++;
        while (send_count[idrecv] == 0 && idrecv < np)
            idrecv++;
    }

    if (nsend != recvoffset || nrecv != sendoffset)
        printf("error in caculate send_rowid\n");


    for (i = 0; i < nsend; i++)
        send_rowid[i] = send_rowid[i] - myfstrow;

    (*SPMul).N_loc = myN_loc;
    (*SPMul).fstrow = myfstrow;

    (*SPMul).nnz_in = nnz_in;
    (*SPMul).ia_in = ia_in;
    (*SPMul).ja_in = ja_in;
    (*SPMul).a_in = a_in;

    (*SPMul).nnz_out = nnz_out;
    (*SPMul).ia_out = ia_out;
    (*SPMul).ja_out = ja_out;
    (*SPMul).a_out = a_out;

    (*SPMul).nrecv = nrecv;
    (*SPMul).rowid_out = rowid_out;
    (*SPMul).recv_count = recv_count;

    (*SPMul).nsend = nsend;
    (*SPMul).send_count = send_count;
    (*SPMul).send_rowid = send_rowid;
}

/*********************************************
           初始化线性求解器
***********************************************/
int AddPre(int N_loc, int fstrow, int *ia_loc, int *ja_loc);
void InitSystem(LinearSystem *LSystem, SparseMatrixLoc A_loc, Info *info)
{
    int myid = (*info).myid;

    int myN_loc = A_loc.N_loc[myid];
    int A_nnz_loc = A_loc.nnz_loc[myid];
    int myfstrow = A_loc.fstrow[myid];
    int *ia_loc = A_loc.ia_loc;
    int *ja_loc = A_loc.ja_loc;
    double *a_loc = A_loc.a_loc;

    /******************************** mynnz_loc **************************************/
    int mynnz_loc;
    mynnz_loc = AddPre(myN_loc, myfstrow, ia_loc, ja_loc);

    /******************************** myic_loc myjc_loc myc_loc myc_diag  *********************/
    int *myic_loc = calloc(mynnz_loc, sizeof(int));
    int *myjc_loc = calloc(mynnz_loc, sizeof(int));
    double *myc_diag = calloc(myN_loc, sizeof(double));
    mumps_double_complex *myc_loc = calloc(mynnz_loc, sizeof(mumps_double_complex));


    int i, j, k, from, to, rowid, colid, tag, offset = 0;
    for (i = 0; i < myN_loc; i++)
    {
        rowid = i + 1 + myfstrow;

        from = ia_loc[i] - 1;
        to = ia_loc[i + 1] - 1;

        tag = 0;


        for (j = from; j < to; j++)
        {
            colid = ja_loc[j];
            if (colid < rowid)
            {
                myic_loc[offset] = rowid;
                myjc_loc[offset] = colid;
                myc_loc[offset].r = -a_loc[j];
                offset = offset + 1;
            }

            else if (colid == rowid)
            {
                myic_loc[offset] = rowid;
                myjc_loc[offset] = rowid;
                offset = offset + 1;
                myc_diag[i] = a_loc[j];
                tag = 1;
            }

            else if (colid > rowid)
            {
                if (tag == 0)
                {
                    myic_loc[offset] = rowid;
                    myjc_loc[offset] = rowid;
                    offset = offset + 1;
                    tag = 1;

                    myic_loc[offset] = rowid;
                    myjc_loc[offset] = colid;
                    myc_loc[offset].r = -a_loc[j];
                    offset = offset + 1;
                }
                else
                {
                    myic_loc[offset] = rowid;
                    myjc_loc[offset] = colid;
                    myc_loc[offset].r = -a_loc[j];
                    offset = offset + 1;
                }
            }
        }

        if (tag == 0)
        {
            myic_loc[offset] = rowid;
            myjc_loc[offset] = rowid;
            offset = offset + 1;
        }
    }
    if (offset != mynnz_loc)
        printf("error  %d  %d  %d \n", offset, mynnz_loc, k);

    /******************************** nnz_loc  N_loc *C_nnz_loc *C_N_loc  *********************/
    int *C_nnz_loc = calloc((*info).hostnp, sizeof(int));
    int *C_N_loc = calloc((*info).hostnp, sizeof(int));
    MPI_Allgather(&mynnz_loc, 1, MPI_INT, C_nnz_loc, 1, MPI_INT, (*info).hostcomm);
    MPI_Allgather(&myN_loc,   1, MPI_INT, C_N_loc,   1, MPI_INT, (*info).hostcomm);

    int nnz_loc = 0;
    int N_loc = 0;
    for (i = 0; i < (*info).hostnp; i++)
    {
        nnz_loc = nnz_loc + C_nnz_loc[i];
        N_loc = N_loc + C_N_loc[i];
    }


    /******************************** ic_loc jc_loc c_loc c_diag  *********************/
    int *ic_loc = calloc(nnz_loc, sizeof(int));
    int *jc_loc = calloc(nnz_loc, sizeof(int));
    double *c_diag = calloc(N_loc, sizeof(double));
    mumps_double_complex *c_loc = calloc(nnz_loc, sizeof(mumps_double_complex));

    MPI_Datatype MPI_complex;
    MPI_Type_contiguous(2, MPI_DOUBLE, &MPI_complex);
    MPI_Type_commit(&MPI_complex);

    int *displs = calloc((*info).hostnp, sizeof(int));
    for (i = 1; i < (*info).hostnp; i++)
    {
        displs[i] = displs[i - 1] + C_nnz_loc[i - 1];
    }
    MPI_Allgatherv(myic_loc, mynnz_loc, MPI_INT, ic_loc, C_nnz_loc, displs, MPI_INT, (*info).hostcomm);
    MPI_Allgatherv(myjc_loc, mynnz_loc, MPI_INT, jc_loc, C_nnz_loc, displs, MPI_INT, (*info).hostcomm);
    MPI_Allgatherv(myc_loc, mynnz_loc, MPI_complex, c_loc, C_nnz_loc, displs, MPI_complex, (*info).hostcomm);

    displs[0] = 0;
    for (i = 1; i < (*info).hostnp; i++)
    {
        displs[i] = displs[i - 1] + C_N_loc[i - 1];
    }
    MPI_Allgatherv(myc_diag, myN_loc, MPI_DOUBLE, c_diag, C_N_loc, displs, MPI_DOUBLE, (*info).hostcomm);
    MPI_Type_free(&MPI_complex);




    /*******************************************************************************
                           out
    ********************************************************************************/

    free(myic_loc);
    free(myjc_loc);
    free(myc_loc);
    free(myc_diag);

    free(displs);
    free(C_nnz_loc);
    free(C_N_loc);

    (*LSystem).nnz_loc = nnz_loc;
    (*LSystem).N_loc = N_loc;
    (*LSystem).N = (*info).N;
    (*LSystem).ic_loc = ic_loc;
    (*LSystem).jc_loc = jc_loc;
    (*LSystem).c_loc = c_loc;
    (*LSystem).c_diag = c_diag;


}












/********************************************************
四级程序                            排序程序
*******************************************************/
void InitMult_minisortja(int N, int block, int *ja, int *Wja, int *ia, int *Wia, double *a, double *Wa)
{
    int begin, id, ida, idb, endida, endidb;

    id = 0;
    for (begin = 0; begin < N; begin = begin + 2 * block)
    {
        ida = begin;
        idb = ida + block;
        endida = ida + block;
        endidb = idb + block;

        if (idb > N) idb = N;
        if (endida > N) endida = N;
        if (endidb > N) endidb = N;

        while (ida < endida || idb < endidb)
        {
            if (ida < endida && idb < endidb)
            {
                if (ja[ida] < ja[idb])
                {
                    Wja[id] = ja[ida];
                    Wia[id] = ia[ida];
                    Wa[id] = a[ida];

                    ida++;
                }
                else
                {
                    Wja[id] = ja[idb];
                    Wia[id] = ia[idb];
                    Wa[id] =  a[idb];
                    idb++;
                }
            }
            else
            {
                if (idb == endidb)
                {
                    Wja[id] = ja[ida];
                    Wia[id] = ia[ida];
                    Wa[id] = a[ida];

                    ida++;
                }
                else
                {
                    Wja[id] = ja[idb];
                    Wia[id] = ia[idb];
                    Wa[id] = a[idb];

                    idb++;
                }
            }

            id++;
        }

    }

    if (id != N)
        printf("error in minisort  id: %d\n ", id);

    if (endidb != N)
        printf("error minisort endidb: %d\n", endidb);

}
void InitMult_sortja(int N, int *ja, int *ia, double *a)
{
    int i, block, tag;
    int *Wja, *Wia;
    double *Wa;

    Wia = calloc(N, sizeof(int));
    Wja = calloc(N, sizeof(int));
    Wa = calloc(N, sizeof(double));

    tag = 0;
    for (block = 1; block < N; block = block * 2)
    {
        if ( tag == 1)
        {
            InitMult_minisortja(N, block, Wja, ja, Wia, ia, Wa, a);
            tag = 0;
        }
        else
        {
            InitMult_minisortja(N, block, ja, Wja, ia, Wia, a, Wa);
            tag = 1;
        }
    }

    if (tag == 1)
    {
        for (i = 0; i < N; i++)
        {
            a[i] = Wa[i];
            ja[i] = Wja[i];
            ia[i] = Wia[i];
        }
    }

    free(Wia);
    free(Wja);
    free(Wa);
}

/********************************************************
加法程序
*******************************************************/
int AddPre(int N_loc, int fstrow, int *ia_loc, int *ja_loc)
{
    int mynnz_loc = 0;
    int i, j, from, to, rowid, colid;
    for (i = 0; i < N_loc; i++)
    {
        rowid = i + 1 + fstrow;

        from = ia_loc[i] - 1;
        to = ia_loc[i + 1] - 1;
        for (j = from; j < to; j++)
        {
            colid = ja_loc[j];
            if (colid != rowid)
            {
                mynnz_loc = mynnz_loc + 1;
            }
        }
    }

    mynnz_loc = N_loc + mynnz_loc;
    return mynnz_loc;
}
