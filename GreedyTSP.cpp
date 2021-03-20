//
// Created by samue on 07/03/2021.
//

#include "GreedyTSP.h"
#include <queue>

struct Enlace    {
    int c1, c2;
    double coste;
};

// this is an strucure which implements the
// operator overloading (fuente https://www.geeksforgeeks.org/stl-priority-queue-for-structure-or-class/)
struct CompareCoste {
    bool operator()(Enlace const& e1, Enlace const& e2)
    {
        // return "true" if "p1" is ordered
        // before "p2", for example:
        return e1.coste < e2.coste;
    }
};

// Constructor
GreedyTSP::GreedyTSP(const std::vector<std::vector<int>> &matrizCoste)
    : mejorCoste(-1.0), coste(matrizCoste)
    {}
