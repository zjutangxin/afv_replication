# Bash script to compile the file
echo "Compile the Debt Mobility program."

gfortran ./src/mod_param.f90 ./src/mod_global.f90 ./src/mod_routines.f90 ./src/solve_system.f90 ./src/debt_main.f90 -o debt_main.out -Ofast -march=native -Wall

echo "Compiled"