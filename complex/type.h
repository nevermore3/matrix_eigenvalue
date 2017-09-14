#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "zmumps_c.h"

typedef struct SPARSEMATRIX  //稀疏矩阵结构
{
    int *ia;   
    int *ja;  
    double *a;
    int N;  //多少行
    int nnz;   //非零元素的个数
} SparseMatrix;


typedef struct SPARSEMATRIXLOC   //稀疏矩阵在各个处理器上的数据结构  行压缩存储格式
{
    int *ia_loc;      
    int *ja_loc;
    mumps_double_complex *rhs
    mumps_double_complex *a_loc;
	int N_loc;           //每个进程的处理的行数
	int nnz_loc;         //每个进程处理的非零元的个数
    int fstrow;          //相对于全局矩阵每个进程处理行数的偏置值

    int N;
} SparseMatrixLoc;




typedef struct SPARSEMATRIXLOC_MUMPS
{
    mumps_double_complex *c_loc;
    double *c_diag;
    mumps_double_complex *rhs;    //mumps计算线性系统的右端项  

    int N;
    int nnz_loc;
    int N_loc;

    int *ic_loc;      //三元组形式的行
    int *jc_loc;       //三元组形式的列
	int max_in;
	int max_out;
    int  nrhs;        
	int *in;
	int *out;
}SparseMatrixLoc_mumps;



