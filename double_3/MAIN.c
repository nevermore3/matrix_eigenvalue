#include "type.h"

void Init(SPMultiply *SPMul, LinearSystem *LSystem, Info *info, int argc, char *argv[]);
void Division(SPMultiply SPMul, LinearSystem LSystem, Info info);
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

    if (myid == 0)
        printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!MPI_Init\n");

    mytime1 = MPI_Wtime();
    Init(&SPMul, &LSystem, &info, argc, argv);
    MPI_Barrier(MPI_COMM_WORLD);
    mytime2 = MPI_Wtime();

    if (info.myid == 0)
        printf("\n Init time %lf \n\n", mytime2 - mytime1);

    /**************  call division function  **************/
    mytime1 = MPI_Wtime();
    Division(SPMul, LSystem, info);
    mytime2 = MPI_Wtime();

    MPI_Barrier(MPI_COMM_WORLD);
    if (info.myid == 0)
        printf("\n Division time %lf \n\n", mytime2 - mytime1);


    /**************  call division function  **************/
    //int neig[] = {1042,48,49,51};
    //info.neig = &neig[info.myid % info.ninterval];
    //double interva[] = {0.01, 0.656079};//0.002429808,0.028502863,0.04105427,0.052094756,0.062774786};
    //info.interval[0] = interva[info.myid % info.ninterval];
    //info.interval[1] = interva[info.myid % info.ninterval+1];
    MPI_Barrier(MPI_COMM_WORLD);
    mytime1 = MPI_Wtime();

    Geteig(SPMul, LSystem, info);

    mytime2 = MPI_Wtime();
    if (info.intervalid == 0)
        printf("\n  my interval  %d  Geteig time %lf \n\n", info.myid, mytime2 - mytime1);

    MPI_Barrier(MPI_COMM_WORLD);
    mytime2 = MPI_Wtime();
    if (info.myid == 0)
        printf("\n Geteig time %lf \n\n", mytime2 - mytime1);

    MPI_Finalize();
    return 0;
}
