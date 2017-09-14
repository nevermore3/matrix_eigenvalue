#include "type.h"

void mul(double *result, double *Y_loc, SPMultiply SPMul, Info info)
{
    int M0 = info.M0;
    int np = info.intervalnp;
    int myid = info.intervalid;


    int N_loc = SPMul.N_loc;
    int nnz_in = SPMul.nnz_in;
    int nnz_out = SPMul.nnz_out;
    int nsend = SPMul.nsend;
    int nrecv = SPMul.nrecv;

    int *ia_in = SPMul.ia_in;
    int *ja_in = SPMul.ja_in;
    double *a_in = SPMul.a_in;

    int *ia_out = SPMul.ia_out;
    int *ja_out = SPMul.ja_out;
    double *a_out = SPMul.a_out;

    int *rowid_out = SPMul.rowid_out;
    int *recv_count = SPMul.recv_count;
    int *send_rowid = SPMul.send_rowid;
    int *send_count = SPMul.send_count;

    int i, j;

    double *B_out;
    B_out = calloc(nrecv * M0, sizeof(double));

    int rowid, colid, idsend, idrecv, sendoffset, recvoffset, number, mpitag = 0;
    MPI_Status status;
    double *sendbuf;
    sendbuf = calloc(N_loc * M0, sizeof(double));
    idsend = 0;
    idrecv = 0;
    recvoffset = 0;
    sendoffset = 0;
    while (send_count[idsend] == 0 && idsend < np)
        idsend++;
    while (recv_count[idrecv] == 0 && idrecv < np)
        idrecv++;
    while (idsend < np || idrecv < np)
    {
        if (idsend == idrecv)
        {
            if (idrecv < myid)
            {
                MPI_Recv(B_out + recvoffset, recv_count[idrecv] * M0, MPI_DOUBLE, idrecv, mpitag, info.intervalcomm, &status);
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
                MPI_Send(sendbuf, send_count[idsend] * M0, MPI_DOUBLE, idsend, mpitag, info.intervalcomm);
                idsend++;
            }
        }
        else
        {
            if (idrecv < idsend || idsend == np)
            {
                MPI_Recv(B_out + recvoffset, recv_count[idrecv] * M0, MPI_DOUBLE, idrecv, mpitag, info.intervalcomm, &status);
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
                MPI_Send(sendbuf, send_count[idsend] * M0, MPI_DOUBLE, idsend, mpitag, info.intervalcomm);
                idsend++;
            }
        }
        while (send_count[idsend] == 0 && idsend < np)
            idsend++;
        while (recv_count[idrecv] == 0 && idrecv < np)
            idrecv++;
    }

    free(sendbuf);



    char char_N = 'N', char_T = 'T', char_V = 'V', char_U = 'U';
    double one = 1, zero = 0;
    int int_one = 1;
    char matdescra[6] = {'G', 'U', 'U', 'F'}; //矩阵A的类型

    mkl_dcsrmm (&char_N, &N_loc, &M0, &N_loc, &one, matdescra, a_in, ja_in, ia_in, ia_in + 1, Y_loc, &N_loc, &zero, result, &N_loc);//work = A * Y;

    for (i = 0; i < nnz_out; i++)
    {
        rowid = ia_out[i];
        colid = ja_out[i];
        for (j = 0; j < M0; j++)
            result[rowid + j * N_loc] =  result[rowid + j * N_loc] + a_out[i] * B_out[colid * M0 + j];
    }

    free(B_out);

}