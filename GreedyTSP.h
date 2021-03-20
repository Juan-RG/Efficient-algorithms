//
// Created by samue on 07/03/2021.
//

#ifndef P1GIT_GREEDYTSP_H
#define P1GIT_GREEDYTSP_H


#include <vector>

class GreedyTSP {
private:
    // Atributos privados
    double mejorCoste;
    std::vector<int> mejorCamino;
    std::vector<std::vector<int>> coste;
    // Métodos privados


public:
    // Métodos públicos
    GreedyTSP(const std::vector<std::vector<int>> &matrizCoste);
    void calcular();
};


#endif //P1GIT_GREEDYTSP_H
