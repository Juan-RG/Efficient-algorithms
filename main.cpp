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
#include <cfloat>
#include <queue>

using namespace std;

const bool DEBUG_INFO = false;

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

void ejecucionAV(int filas, const vector<std::vector<double>> &m, vector<int> &mejorCamino);

void ejecucionFB(int filas, const vector<std::vector<double>> &m, vector<int> &mejorCamino);

void ejecucionPDinamica(const vector<std::vector<double>> &m);

void ejecucionRamificacionYpoda(const vector<std::vector<double>> &m);

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



/**                                                                                            Inicio
 * Inicio Ramificacion y poda
 *
 */
/**
 * Comparo las tuplas por el coste
 * @param lhs
 * @param rhs
 * @return
 */
bool operator<=( const std::tuple<vector<int>, vector<std::vector<double>>, double>& lhs,
                 const std::tuple<vector<int>, vector<std::vector<double>>, double>& rhs ){
    return (std::get<2>(lhs) <= std::get<2>(rhs));
}

bool operator>=( const std::tuple<vector<int>, vector<std::vector<double>>, double>& lhs,
                 const std::tuple<vector<int>, vector<std::vector<double>>, double>& rhs ){
    return (std::get<2>(lhs) >= std::get<2>(rhs));
}


template <class T> struct less_equalH : binary_function <T,T,bool> {
    bool operator() (const T& x, const T& y) const {return x<=y;}
};

template <class T> struct greater_equalH : binary_function <T,T,bool> {
    bool operator() (const T& x, const T& y) const {return x>=y;}
};
/**
 * Metodo de reducir filas o columnas
 * !!!!!!Ojo con el int que tiene comportamientos malos por la precision de los doubles tendremos que ver en el hendrix¡¡¡¡¡¡¡
 * @param matriz
 * @return
 */
double reducirMatriz(std::vector<std::vector<double>>& matriz){

    double MAX = DBL_MAX;
    double total = 0;
    //** se podria marcar en el primer recorrido las columnas que ta tienen 0 para ahorranos el coste de pasar por todas columnas
    for (int i = 0; i < matriz.size(); ++i) {
        double minValue = DBL_MAX;
        for (int j = 0; j < matriz.size(); ++j) {
            if (matriz.at(i).at(j) < minValue){
                minValue = matriz.at(i).at(j);
            }

        }

        if ((minValue != MAX && minValue > 0)){
            total += minValue;
            for (int j = 0; j < matriz.size(); ++j) {
                if (j != i){
                    matriz.at(i).at(j) -= minValue;
                }
            }
        }
    }

    for (int i = 0; i < matriz.size(); ++i) {
        double minValue = DBL_MAX;
        for (int j = 0; j < matriz.size(); ++j) {
            if (matriz.at(j).at(i) < minValue){
                minValue = matriz.at(j).at(i);
            }

        }

        if ((minValue != MAX && minValue > 0)){
            total += minValue;
            for (int j = 0; j < matriz.size(); ++j) {
                if (j != i){
                    matriz.at(j).at(i) -= minValue;
                }
            }
        }
    }

    return total;
}

tuple<vector<int>, vector<std::vector<double>>, double> ramificacionPoda(const std::vector<std::vector<double>>& costes, int dim, vector<int>& mejor)    {
    std::vector<std::vector<double>> matrizNodoRaiz = costes;

    //Diagonales a infinito
    for (int i = 0; i < matrizNodoRaiz.size(); ++i) {
        for (int j = 0; j < matrizNodoRaiz.size(); ++j) {
            if(i == j){
                matrizNodoRaiz.at(i).at(j) = DBL_MAX;
            }
        }

    }

    //obtengo la matriz reducida del nodo raiz
    double coste = reducirMatriz(matrizNodoRaiz);

    vector<int> camino;
    camino.push_back(0);

    //Tupla de valores con el vector de camino, la matriz asociada a ese estado el coste del nodo
    tuple<vector<int>, vector<std::vector<double>>, double> nodoRaiz(camino, matrizNodoRaiz, coste);
    bool objetivo = true;

    //Priority quede para que ordene la cola por orden creciente como declaro arriba en greater_equal
    std::priority_queue<tuple<vector<int>, vector<std::vector<double>>, double>, vector<tuple<vector<int>, vector<std::vector<double>>, double>>,
            greater_equalH<tuple<vector<int>, vector<std::vector<double>>, double>>> explorados;

    //Ya los inserta el 1 el de menos valor
    explorados.push(nodoRaiz);

    //mostrar_queue(explorados);

    //valor del primer camino conseguido para realizar la poda
    double valorPoda = 0;
    bool asignacionPoda = 0;

    //Tupla de retorno
    tuple<vector<int>, vector<std::vector<double>>, double> solucion;


    while (!explorados.empty()){
        //saco el primer nodo de la cola
        tuple<vector<int>, vector<std::vector<double>>, double> nodo = explorados.top();
        explorados.pop();

        double coste = std::get<2>(nodo);
        bool podarNodo = 0;
        //compruebo si la poda esta activa y si hay que podar el nodo
        if (asignacionPoda && (coste > valorPoda)){
            podarNodo = 1;
        }
        //si no hay que podar seguimos generando niveles en el nodo
        if (!podarNodo){
            //Recorro desde 1 ya que el 0 es el 1 nodo
            for (int i = 1; i < costes.size(); ++i) {
                bool nuevoContenido = 0;
                //Nuevo camino de nodo y matriz asociada
                double costeTotalNodo = coste;
                vector<int> caminoNodo(get<0>(nodo));
                vector<std::vector<double>> matrizAsociada(get<1>(nodo));
                //Compruebo si el nodo ha sido visitado o no
                if(find(caminoNodo.begin(), caminoNodo.end(), i) == caminoNodo.end()) {
                    //Si no ha sido visitado indico que hay nuevo contenido para evitar nodos repetidos
                    nuevoContenido = 1;
                    //sumo el valor de coste de ir al nuevo nodo y lo añado al vector
                    costeTotalNodo += matrizAsociada.at(caminoNodo.back()).at(i);
                    caminoNodo.push_back(i);
                    //Realizo la funcion de reduccion de fila y columna de los nodos correspondientes y elimino la posibilidad de volver al primer 0 cuando corresponde
                    for (int j = 0; j < costes.size(); ++j) {
                        //accedo a la antepultima posicion para saber que filas hay que descartar
                        matrizAsociada.at(caminoNodo.rbegin()[1]).at(j) = DBL_MAX;
                        matrizAsociada.at(j).at(caminoNodo.back()) = DBL_MAX;
                        if (caminoNodo.size() < costes.size()){
                            matrizAsociada.at(caminoNodo.back()).at(0) = DBL_MAX;
                        }
                    }
                }
                //reduzco la matriz y lo sumo al coste del nodo
                int costeReducionNodo = reducirMatriz(matrizAsociada);
                //le sumo el coste de reduccion
                costeTotalNodo += costeReducionNodo;
                //si hay contenido lo añadimos a la lista de explorados o nodos vivos
                if (nuevoContenido){
                    tuple<vector<int>, vector<std::vector<double>>, double> nodoCoste(caminoNodo, matrizAsociada, costeTotalNodo);
                    explorados.push(nodoCoste);
                    //si es una solucion completa y es la primera asigno el coste al nodo * Esto se puede sustituir por una ejecucion del algoritmo voraz
                    if (caminoNodo.size() == costes.size()){
                        //si es la 1 vez
                        if (!asignacionPoda){
                            solucion = nodoCoste;
                            valorPoda = costeTotalNodo;
                            asignacionPoda = 1;
                            //si no actualizo cuando sea un nodo mejor
                        } else if (costeTotalNodo <= valorPoda){
                            valorPoda = costeTotalNodo;
                            solucion = nodoCoste;
                        }
                    }
                }
            }
        }
    }

    return solucion;


}

/**                                                                                            Fin
 * Fin Ramificacion y poda
 *
 */


int main() {
    string fichero = R"(..\a15.tsp)"; // Paso como argumento ?
    int filas;
    using std::chrono::high_resolution_clock;
    using std::chrono::duration_cast;
    using std::chrono::duration;
    using std::chrono::milliseconds;
    /**
   * algoritmio voraz
   */
    //rellenarMatriz(ciudades, filas);
    auto m = getMatriz(fichero, filas);

    //asigno -1 en el indice del recorrido
    std::vector<int> mejorCamino;
                    //   double costeMinimo = obtenMejorDynamic(m, filas, mejorCamino);              //todo:aclarar declaracion
  // AV
    auto tInit = chrono::high_resolution_clock::now();
    ejecucionAV(filas, m, mejorCamino);
    auto tEnd = chrono::high_resolution_clock::now();
    auto ms_double = tEnd - tInit;
    cout << "Execution time: " << ms_double.count() << "ms" << endl;
    cout << endl << endl;

    /**
       * fin algoritmo voraz
    */
    /**
   * fuerza bruta
   */
    // FB
    /*
    tInit = chrono::high_resolution_clock::now();
    ejecucionFB(filas, m, mejorCamino);
    tEnd = chrono::high_resolution_clock::now();
    ms_double = tEnd - tInit;
    cout << "Execution time: " << ms_double.count() << "ms" << endl;
    cout << endl << endl;
*/
    /**
     * fin fuerza bruta
     */

    /**
     * Implementacion de programacion dinamica
     */

    tInit = chrono::high_resolution_clock::now();
    ejecucionPDinamica(m);
    tEnd = chrono::high_resolution_clock::now();
    ms_double = tEnd - tInit;
    cout << "Execution time: " << ms_double.count() << "ms" << endl;
    cout << endl << endl;

    /**
     * fin Implementacion de programacion dinamica
     */

    /**
 * Implementacion Ramificacion y poda
 */

    tInit = chrono::high_resolution_clock::now();
    ejecucionRamificacionYpoda(m);
    tEnd = chrono::high_resolution_clock::now();
    ms_double = tEnd - tInit;
    cout << "Execution time: " << ms_double.count() << "ms" << endl;
    cout << endl << endl;

    /**
     * fin Implementacion de programacion dinamica
     */
    //fin programacion dinamica
    return 0;
}
void ejecucionRamificacionYpoda(const vector<std::vector<double>> &m) {
    vector<int> mejorCamino;
    tuple<vector<int>, vector<std::vector<double>>, double> solucion = ramificacionPoda(m, m.size(), mejorCamino);

    double costesFinal;
    vector<int> caminoSolcuion(get<0>(solucion));
    for (int i = 0; i < caminoSolcuion.size(); ++i) {
        cout << caminoSolcuion.at(i) << " ";
    }
    cout << "\n";
    for (int i = 1; i < caminoSolcuion.size(); ++i) {
        costesFinal += m.at(caminoSolcuion.at(i-1)).at(caminoSolcuion.at(i));
    }
    costesFinal += m.at(caminoSolcuion.back()).at(0);

    cout << "coste final "<< costesFinal<< "\n";
}
void ejecucionPDinamica(const vector<std::vector<double>> &m) {
    vector<vector<double>> cities = m;
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
    map<set<int>, int> conjuntosVisitados;
    //llamo al metodo G desde el vertice 0
    cout << dynamicG(0, S, gtab, cities, conjuntosVisitados, vertices) << "\n";


    cout << "0 -> ";
    //Recorro la matriz de minimos vertices para reconstruir el camino
    int conj = S.size();
    int vertice = 0;
    for (int i = 0; i < conj; ++i) {
        vertice = vertices.at(vertice).at(calcularConjunto(S,conjuntosVisitados));
        cout << vertice << " -> ";
        S.erase(vertice);
    }
    cout << "0\n";
}

void ejecucionFB(int filas, const vector<std::vector<double>> &m, vector<int> &mejorCamino) {

    auto costeMinimo = obtenMejorPermutaciones(m, filas, mejorCamino);

    cout << "Mejor camino encontrado, coste = " << costeMinimo << endl;
    for (auto& i : mejorCamino) {
        cout << i << " -> ";
    }
    cout << mejorCamino.at(0) << endl;
}

void ejecucionAV(int filas, const vector<std::vector<double>> &m, vector<int> &mejorCamino) {

    double costeMinimo = obtenMejorAV(m, filas, mejorCamino);
    //int nMin = caminoMinimo(ciudades, caminos, &camino);

    cout << "Mejor camino encontrado, coste = " << costeMinimo << endl;
    for (auto& i : mejorCamino) {
        cout << i << " -> ";
    }
    cout << mejorCamino.at(0) << endl;
}



