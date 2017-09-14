mpiicc -O3 -I ~/lib/inc/ -c MAIN.c
mpiicc -O3 -I ~/lib/inc/ -c MainInit.c
mpiicc -O3 -I ~/lib/inc/ -c MainDivision.c
mpiicc -O3 -I ~/lib/inc/ -c MainGeteig.c
mpiicc -O3 -I ~/lib/inc/ -c DVGetCeig.c
mpiicc -O3 -I ~/lib/inc/ -c DVGetInterval.c
mpiicc -O3 -I ~/lib/inc/ -c DVGetNeig.c
mpiicc -O3 -I ~/lib/inc/ -c DVGetTag.c
mpiicc -O3 -I ~/lib/inc/ -c DVModify.c
mpiicc -O3 -I ~/lib/inc/ -c DVEnd.c
mpiicc -O3 -I ~/lib/inc/ -c linearsystem.c
mpiicc -O3 -I ~/lib/inc/ -c mumps.c
mpiicc -O3 -I ~/lib/inc/ -c feast.c
mpiicc -O3 -I ~/lib/inc/ -c mul.c
mpiicc -O3 -I ~/lib/inc/ -c orth.c
mpiicc -o out MAIN.o MainInit.o MainDivision.o MainGeteig.o DVGetCeig.o DVGetInterval.o DVGetNeig.o DVGetTag.o DVModify.o DVEnd.o linearsystem.o mumps.o feast.o mul.o orth.o -L ~/lib/lib_intelmpi/ -lzmumps -lmumps_common -lpord -lscalapack -llapack -lblacs -lblacsCinit -lblacsF77init -lblas -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core  -liomp5 -limf -lsvml -lirc -lifcore -lpthread -lm