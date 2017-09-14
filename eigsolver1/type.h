#include "mpi.h"
#include "zmumps_c.h"


typedef struct SPARSEMATRIX
{
    int *ia;
    int *ja;
    double *a;
    int N;
    int nnz;
} SparseMatrix;

typedef struct SPARSEMATRIXLOC
{
    int *ia_loc;
    int *ja_loc;
    double *a_loc;

    int *N_loc;
    int *nnz_loc;
    int *fstrow;
} SparseMatrixLoc;

typedef struct SPARSEMATRIXMULTIPLY
{
    int M0;
    int fstrow;
    int N_loc;

    int nnz_in;
    int *ia_in;
    int *ja_in;
    double *a_in;

    int nnz_out;
    int *ia_out;
    int *ja_out;
    double *a_out;

    double *B_out;

    int nrecv;
    int *rowid_out;
    int *recv_count;

    int nsend;
    int *send_rowid;
    int *send_count;
} SPMultiply;

typedef struct LINEARSYSYTEM
{
    int nnz_loc;
    int N_loc;
    int N;

    int *ic_loc;
    int *jc_loc;
    mumps_double_complex *c_loc;
    double *c_diag;

    mumps_double_complex *rhs;

    mumps_double_complex Ze;
    mumps_double_complex jac;

} LinearSystem;


typedef struct INFORMATION
{
    double *interval;
    int neig;
    int *N_loc;

    int N;
    int nGau;
    int M0;
    double Emax;
    double Emin;

    int times;

    int myid;
    int np;

    MPI_Comm hostcomm;
    int hostid;
    int hostnp;
    
    MPI_Comm mumpscomm;
    int mumpsid;
    int mumpsnp;

} Info;


