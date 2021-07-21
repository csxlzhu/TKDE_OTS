This code can be compiled with g++. If the compilation is not successful, please use the later version of g++ and c++11 or higher standard libraries.

For example: g++ "-std=c++11" -O3 -o run OVDO.cpp

The first parameter is the input file name.
The second parameter is the output file name.
The third parameter is the size of selection K.

For example: ./run input.txt output.txt 5

input:
n \\size of tree
x_i y_i \\ x_i is parent of y_i
...n-1 edges...
m \\size of I
x_i w_i \\ weight of x_i is w_i
...m nodes... 

output:
K-size selection set
summary score