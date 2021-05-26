all:
        g++ -std=c++14 -O3 -Ofast -mcpu=niagara2 -mfmaf main.cpp -o tsp