# Bash script to compile the file in serial
# echo "Compile the Debt Mobility program."

# gfortran ./src/mod_param.f90 ./src/mod_global.f90 ./src/mod_routines.f90 \
#      ./src/solve_system.f90 ./src/debt_main.f90 \
#      ./obj/minpack.o -llapack -latlas -lblas \
#      -o debt_main.out -Ofast -march=native -Wall

# echo "Compiled"

# Bash script to compile the file in parallel

echo "Compile the Debt Mobility program in paralle."

mpif90 ./src/mod_param.f90 ./src/mod_global.f90 ./src/mod_routines.f90 \
     ./src/solve_system.f90 ./src/debt_main.f90 \
     ./obj/minpack.o -llapack -latlas -lblas \
     -o debt_main_par.out -Ofast -march=native -Wall

echo "Compiled"