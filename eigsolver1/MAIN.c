#include "type.h"

void Init(SPMultiply *SPMul, LinearSystem *LSystem, Info *info, int argc, char *argv[]);
void Geteig(SPMultiply SPMul, LinearSystem LSystem, Info info);

int main(int argc, char *argv[])
{
    double mytime1, mytime2;
    /**************  Init the program   *******************/
    Info info;
    SPMultiply SPMul;
    LinearSystem LSystem;
    MPI_Init(&argc, &argv);

    int myid;
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);

    mytime1 = MPI_Wtime();
    Init(&SPMul, &LSystem, &info, argc, argv);

    MPI_Barrier(MPI_COMM_WORLD);
    mytime2 = MPI_Wtime();
    if (info.myid == 0)
        printf("\n Init time %lf \n\n", mytime2 - mytime1);

    /**************  comput  eigenvalue **************/
    info.neig = 102;
    info.interval[0] = 0.23;
    info.interval[1] = 2.98;
    
    MPI_Barrier(MPI_COMM_WORLD);
    mytime1 = MPI_Wtime();

    Geteig(SPMul, LSystem, info);

    MPI_Barrier(MPI_COMM_WORLD);
    mytime2 = MPI_Wtime();
    if (info.myid == 0)
        printf("\n Geteig time %lf \n\n", mytime2 - mytime1);

    MPI_Finalize();
    return 0;
}
