mpicc -O3 -I ~/lib/inc/  -c MAIN.c
mpicc -O3 -I ~/lib/inc/  -c C_MainInit.c
mpicc -O3 -I ~/lib/inc/  -c C_EIG.c
mpicc -o out MAIN.o C_MainInit.o C_EIG.o -L ~/lib/lib/ -lzmumps_parmetis -lmumps_common_parmetis  -lpord_parmetis  -lmetis -lparmetis -lscalapack -llapack -lblacs -lblacsCinit -lblacsF77init -lblas -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core  -liomp5 -limf -lsvml -lirc -lifcore -lpthread -lm
