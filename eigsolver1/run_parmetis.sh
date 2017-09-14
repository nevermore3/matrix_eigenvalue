mpicc -O3 -I ~/lib/inc/ -c MAIN.c
mpicc -O3 -I ~/lib/inc/ -c MainInit.c
mpicc -O3 -I ~/lib/inc/ -c MainGeteig.c


mpicc -O3 -I ~/lib/inc/ -c linearsystem.c
mpicc -O3 -I ~/lib/inc/ -c mumps.c
mpicc -O3 -I ~/lib/inc/ -c feast.c
mpicc -O3 -I ~/lib/inc/ -c mul.c

mpicc -o out_feast  MAIN.o MainInit.o MainGeteig.o linearsystem.o mumps.o feast.o mul.o  -L ~/lib/lib/ -lzmumps_parmetis -lmumps_common_parmetis  -lpord_parmetis  -lmetis -lparmetis -lscalapack -llapack -lblacs -lblacsCinit -lblacsF77init -lblas -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core  -liomp5 -limf -lsvml -lirc -lifcore -lpthread -lm


