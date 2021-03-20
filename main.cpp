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


struct nodo{
    double valor;
    std::vector<int> preludio;

    bool operator< (const nodo &a) {
        return valor < a.valor;
    }

};


/**
 * revisar algo estoy haciendo mal
 * @param costes
 * @param dim
 * @param mejor
 * @return
 */
int obtenMejorAV(const std::vector<std::vector<double>>& costes, int dim, vector<int>& mejor) {


    /*
    vector<nodo> frontera;
    for (int i = 1; i < dim; ++i) {
                nodo pareja;
                pareja.valor = costes.at(0).at(i);
                pareja.preludio.push_back(i);
                frontera.push_back(pareja);
                sort(frontera.begin(), frontera.end());
    }


    while (frontera.at(0).preludio.size() != dim){

        nodo n = frontera.at(0);
        frontera.erase(frontera.cbegin());

        for (int i = 1; i < dim; ++i) {
            if(find(n.preludio.begin(), n.preludio.end(), i) == n.preludio.end()) {
                nodo newNodo;
                newNodo.valor = n.valor + costes.at(n.preludio.back()).at(i);
                for (int j = 0; j < n.preludio.size(); ++j) {                   //todo::cambiar a copy eficiente
                    newNodo.preludio.push_back(n.preludio.at(j));
                }
                newNodo.preludio.push_back(i);


                if (newNodo.preludio.size() == dim - 1){
                    newNodo.valor = newNodo.valor + costes.at(newNodo.preludio.back()).at(0);
                    newNodo.preludio.push_back(0);
                    //sort(frontera.begin(),frontera.end());
                }

            }
       }

        if (frontera.at(0).preludio.size() == dim - 1){
            frontera.at(0).valor = frontera.at(0).valor + costes.at(frontera.at(0).preludio.back()).at(0);
            frontera.at(0).preludio.push_back(0);
        }


       // sort(frontera.begin(),frontera.end());

    }*/

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
    string fichero = R"(..\a4.tsp)"; // Paso como argumento ?
    int filas;
    using std::chrono::high_resolution_clock;
    using std::chrono::duration_cast;
    using std::chrono::duration;
    using std::chrono::milliseconds;


    //rellenarMatriz(ciudades, filas);
    auto m = getMatriz(fichero, filas);


    //asigno -1 en el indice del recorrido
    std::vector<int> mejorCamino;
    auto tInit = chrono::high_resolution_clock::now();
    int Cminimo = obtenMejorAV(m, filas, mejorCamino);
    auto tEnd = chrono::high_resolution_clock::now();
    //int nMin = caminoMinimo(ciudades, caminos, &camino);
    chrono::duration<double, std::milli> ms_double = tEnd - tInit;
    cout << "Execution time: " << ms_double.count() << "ms" << endl;

    tInit = chrono::high_resolution_clock::now();
    double costeMinimo = obtenMejorPermutaciones(m, filas, mejorCamino);
    tEnd = chrono::high_resolution_clock::now();

    cout << "Mejor camino encontrado, coste = " << costeMinimo << endl;
    for (auto& i : mejorCamino) {
        cout << i << " -> ";
    }
    cout << mejorCamino.at(0) << endl;
     ms_double = tEnd - tInit;
    cout << "Execution time: " << ms_double.count() << "ms" << endl;


    return 0;
}



