#include <stdio.h>
#include "complex.h"


void C_Init_Aloc(SparseMatrixLoc *A_loc, Info *info, int argc, char *argv[]);

void C_Init_SPMul(SPMultiply *SPMul, SparseMatrixLoc A_loc, Info info);

void C_Init_Cloc(SparseMatrixLoc_mumps *C_loc,  SparseMatrixLoc A_loc, Info info);

void freeA_loc(SparseMatrixLoc A_loc);

void C_Init(SPMultiply *SPMul, SparseMatrixLoc_mumps *C_loc, Info *info, int argc, char *argv[])
{
    SparseMatrixLoc A_loc;

    C_Init_Aloc(&A_loc, info, argc, argv);

    C_Init_SPMul(SPMul, A_loc, *info);

    C_Init_Cloc(C_loc, A_loc,  *info);

    freeA_loc(A_loc);
}





/*********************************************
读入矩阵信息
***********************************************/
void C_Init_Aloc(SparseMatrixLoc *A_loc, Info *info, int argc, char *argv[])
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
    //MPI_Barrier(MPI_COMM_WORLD);
    if (myid == 0)
    {
        free(a);
        free(ia);
        free(ja);
    }

    (*info).myid = myid;
    (*info).np = np;
    (*info).N_loc = N_loc;
    (*info).fstrow = fstrow;
    (*info).N = N;

    (*A_loc).N_loc = myN_loc;
    (*A_loc).nnz_loc = mynnz_loc;
    (*A_loc).fstrow = myfstrow;

    (*A_loc).ia_loc = ia_loc;
    (*A_loc).ja_loc = ja_loc;
    (*A_loc).a_loc = a_loc;

}


/*********************************************
矩阵乘
***********************************************/
void InitMult_sortja(int N, int *ja, int *ia, mumps_double_complex *a);
void C_Init_SPMul(SPMultiply *SPMul, SparseMatrixLoc A_loc, Info info)
{
    int myfstrow = A_loc.fstrow;
    int myN_loc = A_loc.N_loc;
    int mynnz_loc = A_loc.nnz_loc;

    int *ia_loc = A_loc.ia_loc;
    int *ja_loc = A_loc.ja_loc;
    mumps_double_complex *a_loc = A_loc.a_loc;

    int np = info.np;
    int myid = info.myid;
    int *fstrow = info.fstrow;


    int i, j, k, begin, end, from, to, id, idin, idout, sourceid, destid;
    MPI_Status status;
    /*******************************************************************************
    统计内点和外点的个数                nnz_in , nnz_out;
    ********************************************************************************/
    int nnz_in , nnz_out;
    begin = myfstrow;
    end = fstrow[myid + 1] + 1;
    nnz_in = 0;
    nnz_out = 0;
    for (i = 0; i < mynnz_loc; i++)
    {
        if (ja_loc[i] > begin && ja_loc[i] < end)
            nnz_in++;
        else
            nnz_out++;
    }

    /*******************************************************************************
    构造行压缩存储的内点以及三元组存储的外点   *ia_in, *ja_in, *ia_out, *ja_out *a_in, *a_out
    ********************************************************************************/
    int *ia_in, *ja_in, *ia_out, *ja_out;
    mumps_double_complex *a_in, *a_out;
    if (nnz_in > 0)
    {
        ia_in = calloc(myN_loc + 1, sizeof(int));
        ja_in = calloc(nnz_in,    sizeof(int));
        a_in =  calloc(nnz_in,    sizeof(mumps_double_complex));
    }
    if (nnz_out > 0)
    {
        ia_out = calloc(nnz_out, sizeof(int));
        ja_out = calloc(nnz_out, sizeof(int));
        a_out =  calloc(nnz_out, sizeof(mumps_double_complex));
    }

    if (nnz_out > 0 || nnz_in > 0)
    {
        idin = 0;
        idout = 0;
        ia_in[0] = 1;
        for (i = 0; i < myN_loc; i++)
        {
            from = ia_loc[i] - 1;
            to = ia_loc[i + 1] - 1;

            for (j = from; j < to; j++)
            {
                id++;
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
    }

    /*******************************************************************************
    统计计算所需要的外部稠密阵信息————行数，行号以及所在进程号   nrecv, *recv_rowid, *recv_count
    ********************************************************************************/
    int nrecv, *recv_rowid, *recv_count;
    if (nnz_out == 0)
    {
        nrecv = 0;
        recv_count = calloc(np, sizeof(int));
    }
    else
    {
        //nrecv
        nrecv = 1;
        for (i = 1; i < nnz_out; i++)
            if (ja_out[i] > ja_out[i - 1])
                nrecv++;

        //recv_rowid
        recv_rowid = calloc(nrecv, sizeof(int));
        recv_rowid[0] = ja_out[0];
        ja_out[0] = 0;
        j = 1;
        for (i = 1; i < nnz_out; i++)
        {
            if (ja_out[i] > recv_rowid[j - 1])
            {
                recv_rowid[j] = ja_out[i];
                j++;
                ja_out[i] = ja_out[i - 1] + 1;
            }
            else
            {
                ja_out[i] = ja_out[i - 1];
            }
        }
        if (j != nrecv)
            printf("error in caculate recv_rowid \n");

        //recv_count
        recv_count = calloc(np, sizeof(int));
        for (i = 0; i < nrecv; i++)
        {
            j = 0;
            while (recv_rowid[i] >= fstrow[j + 1])
            {
                j++;
            }
            recv_count[j] = recv_count[j] + 1;
        }
        if (j > np)
            printf("error in caculate recv_count\n");

    }


    /*******************************************************************************
    记录需要发送给各个处理的行数   send_count
    ********************************************************************************/
    int  *send_count;
    send_count = calloc(np, sizeof(int));
    int mpitag = 1;
    for (i = 1; i < np; i++)
    {
        sourceid = (myid - i + np) % np;
        destid = (myid + i) % np;
        if (myid < np / 2 )
        {
            MPI_Recv(send_count + sourceid, 1, MPI_INT, sourceid, mpitag, MPI_COMM_WORLD, &status);
            MPI_Send(recv_count + destid,   1, MPI_INT, destid,   mpitag, MPI_COMM_WORLD);
        }

        else
        {
            MPI_Send(recv_count + destid,   1, MPI_INT, destid,   mpitag, MPI_COMM_WORLD);
            MPI_Recv(send_count + sourceid, 1, MPI_INT, sourceid, mpitag, MPI_COMM_WORLD, &status);
        }
    }




    /*******************************************************************************
    统计需要发送的次数   nsend
    ********************************************************************************/
    int nsend;
    nsend = 0;
    for (i = 0; i < np; i++)
    {
        nsend = nsend + send_count[i];
    }

    /*******************************************************************************
     记录需要发送的各个进程号以及行号   send_count
     ********************************************************************************/
    int *send_rowid;
    if (nsend > 0)
        send_rowid = calloc(nsend, sizeof(int));

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
                MPI_Recv(send_rowid + recvoffset, send_count[idrecv], MPI_INT, idrecv, mpitag, MPI_COMM_WORLD, &status);
                recvoffset = recvoffset + send_count[idrecv];
                idrecv++;
            }
            if (idsend > myid)
            {
                MPI_Send(recv_rowid + sendoffset, recv_count[idsend], MPI_INT,   idsend, mpitag, MPI_COMM_WORLD);
                sendoffset = sendoffset + recv_count[idsend];
                idsend++;
            }
        }
        else
        {
            if (idrecv < idsend || idsend == np)
            {
                MPI_Recv(send_rowid + recvoffset, send_count[idrecv], MPI_INT, idrecv, mpitag, MPI_COMM_WORLD, &status);
                recvoffset = recvoffset + send_count[idrecv];
                idrecv++;
            }
            else if (idsend < idrecv || idrecv == np)
            {
                MPI_Send(recv_rowid + sendoffset, recv_count[idsend], MPI_INT, idsend, mpitag, MPI_COMM_WORLD);
                sendoffset = sendoffset + recv_count[idsend];
                idsend++;
            }
        }
        //}
        while (recv_count[idsend] == 0 && idsend < np)
            idsend++;
        while (send_count[idrecv] == 0 && idrecv < np)
            idrecv++;
    }

    if (nsend != recvoffset || nrecv != sendoffset)
        printf("error in caculate send_rowid\n");


    for (i = 0; i < nsend; i++)
        send_rowid[i] = send_rowid[i] - myfstrow;

    /*******************************************************************************
                                     out
    ********************************************************************************/
    (*SPMul).N_loc = myN_loc;
    (*SPMul).nnz_in = nnz_in;
    (*SPMul).nnz_out = nnz_out;
    (*SPMul).nrecv = nrecv;
    (*SPMul).nsend = nsend;

    (*SPMul).ia_in = ia_in;
    (*SPMul).ja_in = ja_in;
    (*SPMul).a_in = a_in;

    (*SPMul).ia_out = ia_out;
    (*SPMul).ja_out = ja_out;
    (*SPMul).a_out = a_out;

    (*SPMul).recv_rowid = recv_rowid;
    (*SPMul).recv_count = recv_count;

    (*SPMul).send_count = send_count;
    (*SPMul).send_rowid = send_rowid;

}

/*********************************************
矩阵加法
***********************************************/
int AddPre(int N_loc, int fstrow, int *ia_loc, int *ja_loc);
void C_Init_Cloc(SparseMatrixLoc_mumps *C_loc,  SparseMatrixLoc A_loc, Info info)
{
    int N_loc = A_loc.N_loc;
    int nnz_loc = A_loc.nnz_loc;
    int fstrow = A_loc.fstrow;
    int *ia_loc = A_loc.ia_loc;
    int *ja_loc = A_loc.ja_loc;
    mumps_double_complex *a_loc = A_loc.a_loc;

    /*******************************************************************************
                                         mynnz_loc;
    ********************************************************************************/
    int mynnz_loc;
    mynnz_loc = AddPre(N_loc, fstrow, ia_loc, ja_loc);

    /*******************************************************************************
                            ic_loc jc_loc c_loc c_diag
    ********************************************************************************/
    int *ic_loc;
    int *jc_loc;
    mumps_double_complex *c_loc;
    mumps_double_complex *c_diag;

    ic_loc = calloc(mynnz_loc, sizeof(int));
    jc_loc = calloc(mynnz_loc, sizeof(int));
    c_loc = calloc(mynnz_loc, sizeof(mumps_double_complex));
    c_diag = calloc(N_loc, sizeof(mumps_double_complex));

    int i, j, from, to, rowid, colid, tag, offset = 0;
    for (i = 0; i < N_loc; i++)
    {
        rowid = i + 1 + fstrow;

        from = ia_loc[i] - 1;
        to = ia_loc[i + 1] - 1;

        tag = 0;

        for (j = from; j < to; j++)
        {
            colid = ja_loc[j];
            if (colid < rowid)
            {
                ic_loc[offset] = rowid;
                jc_loc[offset] = colid;
                c_loc[offset].r = -a_loc[j].r;
                c_loc[offset].i = -a_loc[j].i;
                offset = offset + 1;
            }

            else if (colid == rowid)
            {
                ic_loc[offset] = rowid;
                jc_loc[offset] = rowid;
                offset = offset + 1;
                c_diag[i] = a_loc[j];
                tag = 1;
            }

            else if (colid > rowid)
            {
                if (tag == 0)
                {
                    ic_loc[offset] = rowid;
                    jc_loc[offset] = rowid;
                    offset = offset + 1;
                    tag = 1;

                    ic_loc[offset] = rowid;
                    jc_loc[offset] = colid;
                    c_loc[offset].r = -a_loc[j].r;
                    c_loc[offset].i = -a_loc[j].i;
                    offset = offset + 1;
                }
                else
                {
                    ic_loc[offset] = rowid;
                    jc_loc[offset] = colid;
                    c_loc[offset].r = -a_loc[j].r;
                    c_loc[offset].i = -a_loc[j].i;
                    offset = offset + 1;
                }
            }
        }

        
    }


    /*******************************************************************************
                           out
    ********************************************************************************/
    (*C_loc).nnz_loc = mynnz_loc;
    (*C_loc).N_loc = N_loc;
    (*C_loc).N = info.N;
    (*C_loc).ic_loc = ic_loc;
    (*C_loc).jc_loc = jc_loc;
    (*C_loc).c_loc = c_loc;
    (*C_loc).c_diag = c_diag;

}

/*********************************************
释放空间
***********************************************/
void freeA_loc(SparseMatrixLoc A_loc)
{
    free((A_loc).ia_loc);
    free((A_loc).ja_loc);
    free((A_loc).a_loc);

}

/********************************************************
排序程序
*******************************************************/
void InitMult_minisortja(int N, int block, int *ja, int *Wja, int *ia, int *Wia, mumps_double_complex *a, mumps_double_complex *Wa)
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

/********************************************************
排序程序
*******************************************************/
void InitMult_sortja(int N, int *ja, int *ia, mumps_double_complex *a)
{
    int i, block, tag;
    int *Wja, *Wia;
    mumps_double_complex *Wa;

    Wia = calloc(N, sizeof(int));
    Wja = calloc(N, sizeof(int));
    Wa = calloc(N, sizeof(mumps_double_complex));

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