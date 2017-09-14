#include "complex.h"

void C_Init(SPMultiply *SPMul, SparseMatrixLoc_mumps *C_loc, Info *info, int argc, char *argv[]);
void Eig(SPMultiply SPMul, SparseMatrixLoc_mumps C_loc, Info info);
int main(int argc, char *argv[])
{
    double mytime1, mytime2;
    MPI_Init(&argc, &argv);

    /**************  Init the program   *******************/
    Info info;
    SPMultiply SPMul;
    SparseMatrixLoc_mumps C_loc;


    mytime1 = MPI_Wtime();
    C_Init(&SPMul, &C_loc, &info, argc, argv);
    mytime2 = MPI_Wtime();

    printf("\n\t Init time %lf \n\n", mytime2 - mytime1);

    /**************  test the program   *******************/
    MPI_Barrier(MPI_COMM_WORLD);

    /**************  set the parameters ********************/
    info.nGau = 8;
    info.M0 = 4;
    (info.Emax).r = 3;   (info.Emin).i = 3;
    (info.Emin).r = 0;  (info.Emin).i = 0;

    /**************  eig solve program   *******************/
    Eig(SPMul, C_loc, info);

    /**************  end the program   *******************/
    MPI_Finalize();
    return 0;
}
