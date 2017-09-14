out :MAIN.o C_MainInit.o C_EIG.o
    mpicc -o out MAIN.o MainInit.o MainGeteig.o linearsystem.o mumps.o feast.o mul.o  
     ~/Mumps/MUMPS_5.0.2/lib/libzmumps.a ~/Mumps/MUMPS_5.0.2/lib/libmumps_common.a ~/Mumps/MUMPS_5.0.2/lib/libpord.a -lmkl_scalapack_lp64 -lmkl_blacs_intelmpi_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -openmp -liomp5 -limf -lsvml -lirc -lifcore -lpthread -lm
MAIN.o : MAIN.c
    mpicc -O3 -I  ~/Mumps/MUMPS_5.0.2/include/  -c MAIN.c
MainInit.o: MainInit.c
    mpicc -O3 -I  ~/Mumps/MUMPS_5.0.2/include/  -c MainInit.c
MainGeteig.o :MainGeteig.c
    mpicc -O3 -I  ~/Mumps/MUMPS_5.0.2/include/  -c MainGeteig.c
linearsystem.o:linearsystem.c
    mpicc -O3 -I  ~/Mumps/MUMPS_5.0.2/include/  -c linearsystem.c
mumps.o:mumps.c
    mpicc -O3 -I  ~/Mumps/MUMPS_5.0.2/include/  -c mumps.c
feast.o:feast.c
    mpicc -O3 -I  ~/Mumps/MUMPS_5.0.2/include/  -c feast.c
mul.o:mul.c
    mpicc -O3 -I  ~/Mumps/MUMPS_5.0.2/include/  -c mul.c

##### -lzmumps_parmetis -lmumps_common_parmetis  -lpord_parmetis  -lmetis -lparmetis -lscalapack -llapack -lblacs -lblacsCinit -lblacsF77init -lblas -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core  -liomp5 -limf -lsvml -lirc -lifcore -lpthread -lm


clean :
    rm out *.o
 