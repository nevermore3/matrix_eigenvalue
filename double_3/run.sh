mpicc -O3 -I ~/lib/inc/ -c MAIN.c
mpicc -O3 -I ~/lib/inc/ -c MainInit.c
mpicc -O3 -I ~/lib/inc/ -c MainDivision.c
mpicc -O3 -I ~/lib/inc/ -c MainGeteig.c


mpicc -O3 -I ~/lib/inc/ -c DVGetCeig.c
mpicc -O3 -I ~/lib/inc/ -c DVGetInterval.c
mpicc -O3 -I ~/lib/inc/ -c DVGetNeig.c
mpicc -O3 -I ~/lib/inc/ -c DVGetTag.c
mpicc -O3 -I ~/lib/inc/ -c DVModify.c
mpicc -O3 -I ~/lib/inc/ -c DVEnd.c


mpicc -O3 -I ~/lib/inc/ -c linearsystem.c
mpicc -O3 -I ~/lib/inc/ -c mumps.c
mpicc -O3 -I ~/lib/inc/ -c feast.c
mpicc -O3 -I ~/lib/inc/ -c mul.c
mpicc -O3 -I ~/lib/inc/ -c orth.c


mpicc -o bsub_1000 MAIN.o MainInit.o MainDivision.o MainGeteig.o DVGetCeig.o DVGetInterval.o DVGetNeig.o DVGetTag.o DVModify.o DVEnd.o linearsystem.o mumps.o feast.o mul.o orth.o -L ~/lib/lib/ -lzmumps -lmumps_common -lpord -lscalapack -llapack -lblacs -lblacsCinit -lblacsF77init -lblas -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core  -liomp5 -limf -lsvml -lirc -lifcore -lpthread -lm

