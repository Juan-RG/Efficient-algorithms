#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <algorithm>
#include <vector>

using namespace std;

const bool DEBUG_INFO = true;

const bool DEBUG_INFO_EXTENDED = DEBUG_INFO & false;

/* START DEPRECATED
void rellenarMatriz(int ciudades[100][100], int& tamanio){    //toDo: Falta hacer que se pueda pasar el fichero por arv
    ifstream datos("C:\\Users\\samue\\Desktop\\AlgoritmiaBasica\\Practica\\a4.tsp");
    string linea;
    int filas = 0;
    while (getline(datos,linea)){
        stringstream s_stream(linea);
        string substr;
        int columna = 0;
        while (s_stream >> substr){
            int numero = stoi(substr);
            ciudades[filas][columna] = numero;
            columna++;
        }
        filas++;
    }
    tamanio = filas;
}




//metodo que recorre los caminos y compara cual es el menor
int caminoMinimo(int ciudades[100][100], vector<vector<int>> caminos, int *indiceCaminoMin) {
    int min = INT_MAX;

    for (int i = 0; i < caminos.size(); i++) {
        int coste = 0;
        for (int j = 0; j < caminos.at(i).size(); j++) {
            if(j == caminos.at(i).size() - 1){
                coste = coste + ciudades[caminos.at(i).at(j)][caminos.at(i).at(0)];
            }else{
                coste = coste + ciudades[caminos.at(i).at(j)][caminos.at(i).at(j + 1)];
            }
        }
        //Comprobamos si el es camino con < coste y asignamos indice
        if(coste < min){
            min = coste;
            *indiceCaminoMin = i;
        }

    }
    return min;
}

//metodo que recorre los caminos y compara cual es el menor
int caminoMinimo(std::vector<std::vector<int>> ciudades, vector<vector<int>> caminos, int *indiceCaminoMin) {
    int min = INT_MAX;

    for (int i = 0; i < caminos.size(); i++) {
        int coste = 0;
        for (int j = 0; j < caminos.at(i).size(); j++) {
            if(j == caminos.at(i).size() - 1){
                coste = coste + ciudades[caminos.at(i).at(j)][caminos.at(i).at(0)];
            }else{
                coste = coste + ciudades[caminos.at(i).at(j)][caminos.at(i).at(j + 1)];
            }
        }
        //Comprobamos si el es camino con < coste y asignamos indice
        if(coste < min){
            min = coste;
            *indiceCaminoMin = i;
        }

    }
    return min;
}


vector<vector<int>> generaCaminosRecursivo(int dim, vector<int>preludio)    {
    vector<vector<int>> caminos;
    if(preludio.size() == dim)    {
        caminos.push_back(preludio);
    }
    else {
        for(int j = 0; j < dim; j++)    {
            if(find(preludio.begin(), preludio.end(), j) == preludio.end())  {
                vector<int> preludioSiguiente = preludio;
                preludioSiguiente.push_back(j);
                vector<vector<int>> aux = generaCaminosRecursivo(dim, preludioSiguiente);
                caminos.insert(caminos.end(), aux.begin(), aux.end());
            }
        }
    }
    return caminos;
}

vector<vector<int>> generaCaminos(int dim)    {
    vector<vector<int>> caminos;
    vector<int> preludio(1);
    for(int j = 0; j < dim; j++)    {
        preludio.at(0) = j;
        vector<vector<int>> aux = generaCaminosRecursivo(dim, preludio);
        caminos.insert(caminos.end(), aux.begin(), aux.end());
    }
    return caminos;
}
END DEPRECATED */

// Lee y devuelve la matriz de coste, siendo tamanio la dimensión de la misma
std::vector<std::vector<int>> getMatriz(const string& datosEntrada, int& tamanio)   {
    std::vector<std::vector<int>> m;
    ifstream datos(datosEntrada);
    string linea;
    int filas = 0;
    while (getline(datos,linea)){
        m.emplace_back();
        stringstream s_stream(linea);
        string substr;
        while (s_stream >> substr){
            int numero = stoi(substr);
            m[filas].push_back(numero);
        }
        filas++;
    }
    tamanio = filas;
    return m;
}

// Dado un camino, devuelve el coste
int valorarCamino(const std::vector<std::vector<int>>& costes, std::vector<int> camino) {
    int coste = 0;
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
int obtenMejorRecursivo(const std::vector<std::vector<int>>& costes, int dim, vector<int>preludio, vector<int>& mejor)    {
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

int obtenMejor(const std::vector<std::vector<int>>& costes, int dim, vector<int>& mejor)    {
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

int obtenMejorPermutaciones(const std::vector<std::vector<int>>& costes, int dim, vector<int>& mejor)    {
    int mejorCoste = INT_MAX, costeAux;
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
        if(costeAux < mejorCoste)   {
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
    string fichero = R"(C:\Users\samue\Desktop\AlgoritmiaBasica\Practica\a4.tsp)"; // Paso como argumento ?
    int filas;
    //rellenarMatriz(ciudades, filas);
    auto m = getMatriz(fichero, filas);

    //asigno -1 en el indice del recorrido
    std::vector<int> mejorCamino;
    //int nMin = caminoMinimo(ciudades, caminos, &camino);
    int costeMinimo = obtenMejorPermutaciones(m, filas, mejorCamino);

    cout << "Mejor camino encontrado, coste = " << costeMinimo << endl;
    for (auto& i : mejorCamino) {
        cout << i << " -> ";
    }
    cout << mejorCamino.at(0);
    cout<<"\n";
    return 0;
}

