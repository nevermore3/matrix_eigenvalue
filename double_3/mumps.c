#include "type.h"

/******************************************************
初始化MUMPS
******************************************************/
void Mumps_Init(ZMUMPS_STRUC_C *id, Info info)
{
    //MPI_Barrier(info.mumpscomm);

    (*id).par = 1;
    (*id).sym = 0;
    (*id).comm_fortran = MPI_Comm_c2f(info.mumpscomm);

    (*id).job = -1;
    zmumps_c(id);
    if ((*id).info[0] != 0)
        printf("InitMumps info: %d %d \n", (*id).info[0], (*id).info[1]);

}


/******************************************************
MUMPS分析阶段
******************************************************/
void Mumps_Analize(ZMUMPS_STRUC_C *id, LinearSystem LSystem, Info info)
{
    (*id).icntl[2] = -1;
    (*id).icntl[3] = 1;

    (*id).icntl[4] = 0;
    (*id).icntl[17] = 3;

    if (info.mumpsid == 0)
    {
        (*id).rhs = LSystem.rhs;
        (*id).nrhs = info.M0;
        (*id).lrhs = LSystem.N;
        (*id).n = LSystem.N;
    }

    (*id).nz_loc = LSystem.nnz_loc;
    (*id).irn_loc = LSystem.ic_loc;
    (*id).jcn_loc = LSystem.jc_loc;
    (*id).a_loc = LSystem.c_loc;


    (*id).job = 1;
    zmumps_c(id);

    if ((*id).info[0] != 0)
        printf("AnalizeMumps info: %d %d \n", (*id).info[0], (*id).info[1]);
}


/******************************************************
MUMPS分解阶段
******************************************************/
void Mumps_Factor(ZMUMPS_STRUC_C *id)
{
    // MPI_Barrier(info.mumpscomm);

    (*id).job = 2;
    zmumps_c(id);
    if ((*id).info[0] != 0)
        printf("FactorMumps info: %d %d \n", (*id).info[0], (*id).info[1]);
}


/******************************************************
MUMPS求解阶段
******************************************************/
void Mumps_Solve(ZMUMPS_STRUC_C *id)
{
    // MPI_Barrier(info.mumpscomm);

    (*id).job = 3;
    zmumps_c(id);
    if ((*id).info[0] != 0)
        printf("SolveMumps info: %d %d \n", (*id).info[0], (*id).info[1]);
}

/******************************************************
MUMPS结束
******************************************************/
void Mumps_End(ZMUMPS_STRUC_C *id)
{
    (*id).job = -2;
    zmumps_c(id);
    if ((*id).info[0] != 0)
        printf("SolveMumps info: %d %d \n", (*id).info[0], (*id).info[1]);
}
