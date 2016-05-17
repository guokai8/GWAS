----------------------
Calulate the kinship
This program is writed for calculate the kinship,eigvector and eigvalue

1.Requirements

* To use the program,you should have Linux or UNIX-like platformi with g++

2.How to use

You can install the program just one step
Installation procedure:
* Step 1:
make
* then you could copy the binary file "kin_cal" anywhere you like
* Use the command:
* kin_cal infile parents number
* for calculating the kinship etc
* the result file are "kinship_res.csv","eigval.csv" and "eigvec.csv"
If you want remove the program just use "make clean"

3.Example

kin_cal example/example.csv 20

4.How to parallel it (If you have larger data)

* Step 1

Install OpenBLAS 
* Step 2

Install armadillo
* Step 3

g++ -o ../kin_cal src/main.cpp -I/where install armadillo/ -L/where the openblas/ -larmadillo -lopenblas -llapack
* Step 4
Just have fun!

