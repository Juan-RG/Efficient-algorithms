#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <algorithm>
#include <vector>
#include <chrono>
#include <numeric>
#include <cmath>
#include <set>
#include <map>

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

struct arista {
    int nodo;
    double coste;
};

void obtenerMejor(const vector<std::vector<double>> &vector, int dim, std::vector<int>& vect, std::vector<std::vector<double>> &vector1);

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

void obtenerMejor(const vector<std::vector<double>> &datos, int dim, vector<int>& mejor ,std::vector<std::vector<double>>& matrizDatos) {

    for (int i = 0; i < dim; ++i) {
        for (int j = 0; j < dim; ++j) {
            matrizDatos[i][j] = datos[i][j];
        }
    }

    for (int i = 0; i < dim; ++i) {
        cout << i <<"  ";
        for (int j = 0; j < pow(2, dim - 1); ++j) {
            cout << matrizDatos[i][j] << " ";
        }
        cout << "\n";
    }
}

/*
double obtenMejorDynamic(const std::vector<std::vector<double>>& costes, int dim, vector<int>& mejor) {
    std::vector<std::vector<double>> costesConjuntos;
    for (int i = 0; i < dim; ++i) {
        costesConjuntos.emplace_back();
        for (int j = 0; j < pow(2,(dim - 1)); ++j) {
            costesConjuntos.at(i).emplace_back(-1);
        }
    }
    obtenerMejor(costes, dim, mejor,costesConjuntos);

    double masCorto, distancia, nat
    if (){ }
    return 1;
}
*/


/**                                                                                             INICIO
 * programacion dinamica
 *
 */
void insertarConjunto(set<int> S, map<std::set<int>,int>& mapaConjuntos, int max) {
    static int n = 0;
    //si el conjunto no existe le asigno su vvalor de N
    if (mapaConjuntos.find(S) == mapaConjuntos.end()){
        std::pair<std::set<int>, int> pareja(S,n++);
        mapaConjuntos.insert(pareja);
        if (n == max){
            n = 0;
        }
    }
}
/*
 * devuelvo el valor de matriz asignado al conjunto
 */
double calcularConjunto(set<int> S, map<std::set<int>, int>& mapaConjuntos) {
    return mapaConjuntos.find(S)->second;
}

double dynamicG(int i, set<int>& S,vector<vector<double>>& gtab ,vector<vector<double>>& costes, map<set<int>, int>& conjuntosVisitados, vector<vector<double>>& vertices){

    //asigno un valor entero de la matriz al conjunto ----Mejorable----
    insertarConjunto(S, conjuntosVisitados, gtab.at(0).size());

    //si conjunto vacio volvemos al vertice de inicio
    if (S.empty()){
        return costes[i][0];
    }
    //si valor diferente de -1 devolvemos el valor ya calculado del conjunto
    if (gtab.at(i).at(calcularConjunto(S,conjuntosVisitados)) != -1){
        return gtab.at(i).at(calcularConjunto(S,conjuntosVisitados));
    }

    //buscamos en todos posibles conjuntos y nos quedamos con el minimo
    double masCorto = INT_MAX;
    for(auto j : S){
        //copia para generar varios
        set<int> newS(S);
        //eliminamos J del conjunto
        newS.erase(j);
        //llamamos a la funcion recusivamente
        double distancia = costes[i][j] + dynamicG(j,newS, gtab, costes, conjuntosVisitados, vertices);
        //asignamos el minimo y el vertice al que hay que ir
        if (distancia < masCorto){
            masCorto = distancia;
            vertices.at(i).at(calcularConjunto(S,conjuntosVisitados)) = j;
        }
    }
    //asignamos el valor al conjunto
    gtab.at(i).at(calcularConjunto(S,conjuntosVisitados)) = masCorto;

    return masCorto;
}

/**                                                                                             Fin
 * fin programacion dinamica
 *
 */



int main() {
    string fichero = R"(..\a7.tsp)"; // Paso como argumento ?
    int filas;
    using std::chrono::high_resolution_clock;
    using std::chrono::duration_cast;
    using std::chrono::duration;
    using std::chrono::milliseconds;


    //rellenarMatriz(ciudades, filas);
    auto m = getMatriz(fichero, filas);


    //asigno -1 en el indice del recorrido
    std::vector<int> mejorCamino;
                    //   double costeMinimo = obtenMejorDynamic(m, filas, mejorCamino);              //todo:aclarar declaracion
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

    /**
     * Implementacion de programacion dinamica
     */

    vector<vector<double>> cities = getMatriz(fichero, filas);

    //matriz de [vertices][2^vertices - 1] --> creo que se podria reducir accediendo a costes con conjuntos de tamaño S == 1
    vector<vector<double>> gtab(cities.size());
    for (int i = 0; i < cities.size(); ++i) {
        gtab.at(i) = vector<double>((pow(2,cities.size())) - 1, -1);
    }
    //matriz de [vertices][2^vertices - 1] para guardar el camino de vuelta
    vector<vector<double>> vertices(cities.size());
    for (int i = 0; i < cities.size(); ++i) {
        vertices.at(i) = vector<double>((pow(2,cities.size())) - 1, -1);
    }

    //Genero el conjunto de posibles vertices para ir visitando
    set<int> S;
    for (int i = 1; i < cities.size(); ++i) {
        S.insert(i);
    }

    //gmapa de (conjuntos, indice) == conjuntosVisitados[number][S]
    std::map<set<int>, int> conjuntosVisitados;
    //llamo al metodo G desde el vertice 0
    cout << dynamicG(0, S, gtab, cities, conjuntosVisitados, vertices) << "\n";


    cout << "0 - ";
    //Recorro la matriz de minimos vertices para reconstruir el camino
    int conj = S.size();
    int vertice = 0;
    for (int i = 0; i < conj; ++i) {
        vertice = vertices.at(vertice).at(calcularConjunto(S,conjuntosVisitados));
        cout << vertice << " - ";
        S.erase(vertice);
    }
    cout << "0\n";

    //fin programacion dinamica
    return 0;
}



