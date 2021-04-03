#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <algorithm>
#include <vector>
#include <chrono>
#include <numeric>

using namespace std;

const bool DEBUG_INFO = true;

const bool DEBUG_INFO_EXTENDED = DEBUG_INFO & false;




struct nodo{
    double valor;
    std::vector<int> preludio;

    bool operator< (const nodo &a) {
        return valor < a.valor;
    }

};


struct arista {
    int nodo;
    double coste;
};

double obtenMejorAV(const std::vector<std::vector<double>>& costes, int dim, vector<int>& mejor) {
    vector<int> mejorCamino; mejorCamino.push_back(0);   // Nodo inicial
    double coste = 0, lastNode = 0;  // Coste inicial, nodo inicial
    struct arista mejorAristaLocal{-1, INT_MAX};

    cout << "Llamada a obtenMejorAV, dim " << dim << endl;

    for(int n = 0; n < dim - 1; n++)  {
        for(int i = 0; i < dim; i++) {
            // Descarta nodos repetidos y el propio nodo
            if (i != lastNode && find(mejorCamino.begin(), mejorCamino.end(), i) == mejorCamino.end()) {
                // Comprueba si es la arista con menor coste
                if (costes[lastNode][i] > 0 && costes[lastNode][i] < mejorAristaLocal.coste) {
                    mejorAristaLocal.nodo = i;
                    mejorAristaLocal.coste = costes[lastNode][i];
                }
            }
        }
        mejorCamino.push_back(lastNode = mejorAristaLocal.nodo);
        coste += mejorAristaLocal.coste;
        if (DEBUG_INFO) cout << "Anyadido nodo " << lastNode << " con coste " << mejorAristaLocal.coste << endl;

        // Reset valores
        mejorAristaLocal.coste = INT_MAX;
        mejorAristaLocal.nodo = -1;
    }


    mejor = mejorCamino;
    return coste + costes[lastNode][0];
}

// Lee y devuelve la matriz de coste, siendo tamanio la dimensión de la misma
std::vector<std::vector<double>> getMatriz(const string& datosEntrada, int& tamanio)   {
    std::vector<std::vector<double>> m;
    ifstream datos(datosEntrada);
    string linea;
    int filas = 0;
    while (getline(datos,linea)){
        m.emplace_back();
        stringstream s_stream(linea);
        string substr;
        while (s_stream >> substr){
            double numero = stod(substr);
            m[filas].push_back(numero);
        }
        filas++;
    }
    tamanio = filas;
    return m;
}

// Dado un camino, devuelve el coste
double valorarCamino(const std::vector<std::vector<double>>& costes, std::vector<int> camino) {
    double coste = 0;
    for (int j = 0; j < camino.size(); j++) {
        if(j == camino.size() - 1){
            if(costes[camino[j]][camino[0]] == 0) return -1;    // Camino cortado, valor 0
            else {
                coste = coste + costes[camino[j]][camino[0]];
            }
        } else {
            if(costes[camino[j]][camino[j + 1]] == 0) return -1;
            else {
                coste = coste + costes[camino[j]][camino[j + 1]];
            }
        }
    }
    return coste;
}

// Algoritmo recursivo de fuerza bruta
int obtenMejorRecursivo(const std::vector<std::vector<double>>& costes, int dim, vector<int>preludio, vector<int>& mejor)    {
    if(preludio.size() == dim)    {
        mejor = preludio;
        return valorarCamino(costes, preludio);
    }
    else {
        int minimo = INT_MAX;
        int aux;
        std::vector<int> auxCamino;
        for(int j = 0; j < dim; j++)    {
            if(find(preludio.begin(), preludio.end(), j) == preludio.end()) {   // Si no está en el camino
                vector<int> preludioSiguiente = preludio;
                preludioSiguiente.push_back(j);
                aux = obtenMejorRecursivo(costes, dim, preludioSiguiente, auxCamino);
                if (aux != -1 && aux < minimo) {   // No ha encontrado caminos validos con preludio [preludio]U[j]
                    minimo = aux;
                    mejor = auxCamino;
                    if(DEBUG_INFO)  {
                        cout << "Profundidad " << preludio.size() << endl;
                        cout << "Encontrado un mejor camino, coste = " << minimo << " :";
                        for (auto& i : mejor) {
                            cout << i << " -> ";
                        }
                        cout << mejor[0] << endl;
                    }
                }
            }
        }
        if(minimo != INT_MAX) return minimo;    // Devuelve el mínimo coste encontrado
        else return -1;     // No ha encontrado ningun camino válido con preludio [preludio]
    }
}
/*
int obtenMejor(const std::vector<std::vector<double>>& costes, int dim, vector<int>& mejor)    {
    int mejorCoste = INT_MAX;
    vector<int> preludio(1);
    preludio.at(0) = 0;
    int aux = obtenMejorRecursivo(costes, dim, preludio, mejor);
    if (aux != -1 && aux < mejorCoste) {   // No ha encontrado caminos validos con preludio [preludio]U[j]
        mejorCoste = aux;
        if(DEBUG_INFO)  {
            cout << "Encontrado el mejor camino, coste = " << mejorCoste << " :";
            for (auto& i : mejor) {
                cout << i << " -> ";
            }
            cout << mejor[0] << endl;
        }
    }
    return mejorCoste;
}
*/
double obtenMejorPermutaciones(const std::vector<std::vector<double>>& costes, int dim, vector<int>& mejor)    {
    double mejorCoste = INT_MAX, costeAux;
    vector<int> aux;
    for(int i = 0; i < dim; i++)
        aux.push_back(i);   // Crea el primer camino

    do {
        if(DEBUG_INFO_EXTENDED)  {
            cout << "nuevo camino ";
            for (auto& i : aux) {
                cout << i << " -> ";
            }
            cout << aux[0] << endl;
        }
        costeAux = valorarCamino(costes, aux);
        if(costeAux < mejorCoste && costeAux > 0)   {
            mejorCoste = costeAux;
            mejor = aux;
            if(DEBUG_INFO)  {
                cout << "Encontrado un mejor camino, coste = " << mejorCoste << " :";
                for (auto& i : mejor) {
                    cout << i << " -> ";
                }
                cout << mejor[0] << endl;
            }
        }
    } while (next_permutation( ++(aux.begin()), aux.end()));

    return mejorCoste;
}



int main() {
    string fichero = R"(..\a400.tsp)"; // Paso como argumento ?
    int filas;
    using std::chrono::high_resolution_clock;
    using std::chrono::duration_cast;
    using std::chrono::duration;
    using std::chrono::milliseconds;


    //rellenarMatriz(ciudades, filas);
    auto m = getMatriz(fichero, filas);


    //asigno -1 en el indice del recorrido
    std::vector<int> mejorCamino;

    // AV
    auto tInit = chrono::high_resolution_clock::now();
    double costeMinimo = obtenMejorAV(m, filas, mejorCamino);
    auto tEnd = chrono::high_resolution_clock::now();
    //int nMin = caminoMinimo(ciudades, caminos, &camino);
    chrono::duration<double, std::milli> ms_double = tEnd - tInit;
    cout << "Execution time: " << ms_double.count() << "ms" << endl;

    cout << "Mejor camino encontrado, coste = " << costeMinimo << endl;
    for (auto& i : mejorCamino) {
        cout << i << " -> ";
    }
    cout << mejorCamino.at(0) << endl;


    cout << endl << endl;


    // FB
    tInit = chrono::high_resolution_clock::now();
    costeMinimo = obtenMejorPermutaciones(m, filas, mejorCamino);
    tEnd = chrono::high_resolution_clock::now();
    //int nMin = caminoMinimo(ciudades, caminos, &camino);
    ms_double = tEnd - tInit;
    cout << "Execution time: " << ms_double.count() << "ms" << endl;

    cout << "Mejor camino encontrado, coste = " << costeMinimo << endl;
    for (auto& i : mejorCamino) {
        cout << i << " -> ";
    }
    cout << mejorCamino.at(0) << endl;

    return 0;
}



