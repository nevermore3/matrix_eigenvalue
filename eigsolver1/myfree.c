# include "type.h"
void freeSparseMatrix(SparseMatrix A)
{
    int *ia = A.ia;
    int *ja = A.ja;
    double *a = A.a;

    free(ia);
    free(ja);
    free(a);
}

void freeSparseMatrixLoc(SparseMatrixLoc A_loc)
{
    int *ia_loc = A_loc.ia_loc;
    int *ja_loc = A_loc.ja_loc;
    double *a_loc = A_loc.a_loc;

    free(ia_loc);
    free(ja_loc);
    free(a_loc);
}

void freeSPMultiply(SPMultiply SPMul)
{
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

    free(ia_in);
    free(ia_out);
    free(ja_in);
    free(ja_out);
    free(a_in);
    free(a_out);
    free(rowid_out);
    free(recv_count);
    free(send_count);
    free(send_rowid);
}

void freeLINEARSYSYTEM(LinearSystem linearsystem)
{
    int *ic_loc = linearsystem.ic_loc;
    int *jc_loc = linearsystem.jc_loc;
    mumps_double_complex *c_loc = linearsystem.c_loc;
    mumps_double_complex *rhs = linearsystem.rhs;

    int *nnz_loc = linearsystem.nnz_loc;
    int *offset = linearsystem.offset;

    free(ic_loc);
    free(jc_loc);
    free(c_loc);
    free(rhs);
    free(nnz_loc);
    free(offset);
}