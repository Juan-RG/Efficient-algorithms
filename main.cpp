#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <algorithm>
#include <vector>
#include <chrono>
#include <cmath>
#include <set>
#include <map>
#include <cfloat>
#include <queue>

using namespace std;

const bool DEBUG_INFO = false;

const bool DEBUG_INFO_EXTENDED = DEBUG_INFO & false;


void ejecucionAV(int filas, const vector<std::vector<double>> &m, vector<int> &mejorCamino);

void ejecucionFB(int filas, const vector<std::vector<double>> &m, vector<int> &mejorCamino);

void ejecucionPDinamica(const vector<std::vector<double>> &m);

void ejecucionRamificacionYpoda(const vector<std::vector<double>> &m);


/**                                                 ***********************************
 *                                                  *  Algoritmo de algoritmo voraz   *
 *                                                  ***********************************
 */

/**
 * Estructura que representa una arista del problema TSP
 */
struct arista {
    int nodo;
    double coste;
};

/** Calcula el camino empleando el algoritmo voraz (búsquedas locales al último nodo expandido)
 *
 * @param costes: matriz de coste
 * @param dim: dimension del problema
 * @param mejor: referencia a un vector que contiene el mejor camino
 * @return mejor coste encontrado
 */
double obtenMejorAV(const std::vector<std::vector<double>>& costes, int dim, vector<int>& mejor) {
    vector<int> mejorCamino; mejorCamino.push_back(0);      // Nodo inicial
    double coste = 0, lastNode = 0;                         // Coste inicial, nodo inicial
    struct arista mejorAristaLocal{-1, INT_MAX};      //

    //cout << "Llamada a obtenMejorAV, dim " << dim << endl;
    // Ejecución del algoritmo voraz
    for(int n = 0; n < dim - 1; n++)  {
        // Recorre todos los nodos del problema
        for(int i = 0; i < dim; i++) {
            // Descarta nodos repetidos (que ya estén en el camino) y el propio nodo
            if (i != lastNode && find(mejorCamino.begin(), mejorCamino.end(), i) == mejorCamino.end()) {
                // Guarda la arista con menor coste
                if (costes[lastNode][i] > 0 && costes[lastNode][i] < mejorAristaLocal.coste) {
                    mejorAristaLocal.nodo = i;
                    mejorAristaLocal.coste = costes[lastNode][i];
                }
            }
        }

        // Tras encontrar la mejor arista, la anyade al camino
        mejorCamino.push_back(lastNode = mejorAristaLocal.nodo);
        // Actualiza coste
        coste += mejorAristaLocal.coste;
        if (DEBUG_INFO) cout << "Anyadido nodo " << lastNode << " con coste " << mejorAristaLocal.coste << endl;

        // Reset valores locales
        mejorAristaLocal.coste = INT_MAX;
        mejorAristaLocal.nodo = -1;
    }

    // Copia el mejor camino a la variable pasado por referencia
    mejor = mejorCamino;
    // Devuelve el coste del camino encontrado
    return coste + costes[lastNode][0];
}

/**                                                 ***************************************
 *                                                  *  fin Algoritmo de algoritmo voraz   *
 *                                                  ***************************************
 */

/**
 * Lee y devuelve la matriz de coste. Aquellas aristas no existentes se representarán como DBL_MAX (coste maximo)
 * @param datosEntrada Nombre del fichero del leer la matriz de coste
 * @param tamanio Paramétro referenciado en el que se devuelve adicionalmente la dimensión de la matriz
 * @return Matriz de coste leía
 */
std::vector<std::vector<double>> getMatriz(const string& datosEntrada, int& tamanio)   {
    std::vector<std::vector<double>> m;         // Matriz
    ifstream datos(datosEntrada);
    string linea;
    int filas = 0;
    while (getline(datos,linea)){
        m.emplace_back();
        stringstream s_stream(linea);
        string substr;
        // Por cada arista representada en la matriz...
        while (s_stream >> substr){
            // ...la lee...
            double numero = stod(substr);
            // .. y asigna el valor adecuado
            if(numero != 0) m[filas].push_back(numero); // Caso regular
            else m[filas].push_back(DBL_MAX);           // Caso arista no existe
        }
        filas++;
    }
    tamanio = filas;
    return m;   // Devuelve la matriz
}


/**
 * Dado un camino, devuelve el coste
 * @param costes Matriz de costes
 * @param camino Camino a valorar
 * @return Coste del camino en base a la matriz de coste
 */
double valorarCamino(const std::vector<std::vector<double>>& costes, std::vector<int> camino) {
    double coste = 0;
    // Contabiliza el coste total del camino
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
/**                                                 ********************************
 *                                                  *  Algoritmo de fuerza bruta   *
 *                                                  ********************************
 */
/**
 * Algoritmo recursivo de fuerza bruta (No se debe emplear este algoritmo dado que es menos eficiente y puede haber
 * problemas con la pila de llamadas).
 * @param costes Matriz de costes de referencia
 * @param dim Dimensión de la matriz de coste
 * @param preludio
 * @param mejor
 * @return Mejor coste de cualquier camino con el preludio pasado como parametro
 */
int obtenMejorRecursivo(const std::vector<std::vector<double>>& costes, int dim, vector<int>preludio, vector<int>& mejor)    {
    // Caso base, preludio ya es el camino completo
    if(preludio.size() == dim)    {
        mejor = preludio;
        // Devuelve el coste del camino
        return valorarCamino(costes, preludio);
    }
    // Caso recursivo, preludio no incluye todos los nodos (ciudades)
    else {
        int minimo = INT_MAX;
        int aux;
        std::vector<int> auxCamino;
        // Llama recursivamente a obtenMejorRecursivo con [preludio,nodo], siendo nodo to_do aquel nodo no incluido en
        // preludio (genera todos los posibles caminos).
        for(int j = 0; j < dim; j++)    {
            // Comprueba que i no esta en preludio
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
        // Devuelve el mejor
        if(minimo != INT_MAX) return minimo;    // Devuelve el mínimo coste encontrado
        // Si no ha encontrado nada, devuelve -1
        else return -1;     // No ha encontrado ningun camino válido con preludio [preludio]
    }
}


/**
 * Algoritmo de fuerza bruta implementado mediante permutaciones
 * @param costes Matriz de coste
 * @param dim Dimensión de la matriz de costes
 * @param mejor Mejor camino encotrado
 * @return Devuelve el mejor coste encontrado
 */
double obtenMejorPermutaciones(const std::vector<std::vector<double>>& costes, int dim, vector<int>& mejor)    {
    double mejorCoste = INT_MAX, costeAux;
    vector<int> aux;
    for(int i = 0; i < dim; i++)
        aux.push_back(i);   // Crea el primer camino: [0,1,2,...,dim-1]

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

/**                                                 ************************************
 *                                                  *  Fin Algoritmo de fuerza bruta   *
 *                                                  ************************************
 */

/**                                                 *****************************
 *                                                  *   programacion dinamica   *
 *                                                  *****************************
 */

/** Inserta conjuntos con su valor de coste
 *
 * @param S
 * @param mapaConjuntos
 * @param max
 */
void insertarConjunto(set<int> S, map<std::set<int>, int>& mapaConjuntos, int max) {
    static int n = 0;

    if (mapaConjuntos.find(S) == mapaConjuntos.end()){   //si el conjunto no existe le asigno su vvalor de N
        std::pair<std::set<int>, int> pareja(S, n++);
        mapaConjuntos.insert(pareja);
        n = n % max;
    }
}

/** Metodo que devuelve el coste asociado a un conjunto
 *
 * @param S
 * @param mapaConjuntos
 * @return mapaConjuntos[S]
 */
double calcularConjunto(set<int> S, map<std::set<int>, int>& mapaConjuntos) {
    return mapaConjuntos.find(S)->second;
}
/** Metodo de programacion dinamica en la que guardamos el coste de ir a los diferentes conjuntos para que a la hora
 *  de volver a generar el conjunto no tenga que volver a calcular los costes;
 *
 *
 * @param i vertice actual
 * @param S conjunto de vertices visitados
 * @param gtab Matriz de adyacente con el valor de coste del conjunto S
 * @param costes Matriz de costes entre vertices
 * @param conjuntosVisitados Mapa con los conjntos visitados
 * @param vertices matriz de vertices reccoridos
 * @return Valor del coste minimo
 */
double dynamicG(int i, set<int>& S, vector<vector<double>>& gtab, vector<vector<double>>& costes, map<set<int>,
        int>& conjuntosVisitados, vector<vector<double>>& vertices){

    //asigno un valor entero de la matriz al conjunto
    insertarConjunto(S, conjuntosVisitados, gtab.at(0).size());


    if (S.empty()){      //si conjunto vacío volvemos al vértice de inicio
        return costes[i][0];
    }
    //si valor diferente de -1 devolvemos el valor ya calculado del conjunto
    if (gtab.at(i).at(calcularConjunto(S,conjuntosVisitados)) != -1){
        return gtab.at(i).at(calcularConjunto(S,conjuntosVisitados));
    }


    double masCorto = INT_MAX;       //buscamos en todos posibles conjuntos y nos quedamos con el minimo
    for(auto j : S){
        set<int> newS(S);  //copia para generar varios
        newS.erase(j);  //eliminamos J del conjunto

        //llamamos a la funcion recursivamente
        double distancia = costes[i][j] + dynamicG(j,newS, gtab, costes, conjuntosVisitados, vertices);

        if (distancia < masCorto){  //asignamos el minimo y el vertice al que hay que ir
            masCorto = distancia;
            vertices.at(i).at(calcularConjunto(S,conjuntosVisitados)) = j;
        }
    }

    gtab.at(i).at(calcularConjunto(S,conjuntosVisitados)) = masCorto;  //asignamos el valor al conjunto

    return masCorto;
}

/**                                                 *****************************
 *                                                  * fin programacion dinamica *
 *                                                  *****************************
 */



/**                                                 *****************************
 *                                                  *   Ramificacion y poda     *
 *                                                  *****************************
 */


class Nodo   {
public:
    // ID del nodo (qué nodo es)
    int id;

    // Coste y camino
    double costeCamino;
    vector<int> camino;

    // Coste estimado del nodo
    double C;
    // Reducción efectuada sobre la matriz asociada
    double L;
    // Matriz de coste asociada (reducida)
    vector<vector<double>> matrizCosteAsociada;

    Nodo(const int id, const double coste, const vector<int> &camino, const vector<vector<double>>& matrizCoste)
            : id(id), costeCamino(coste), camino(camino), matrizCosteAsociada(matrizCoste)
    {
        this->L = this->reduceMatriz();
        this->C = this->L;
    }

    Nodo(const int id, const int idPadre, const double coste, const double costeEstimadoPadre, const vector<int> &camino, const vector<vector<double>>& matrizCoste)
            : id(id), costeCamino(coste), C(costeEstimadoPadre), camino(camino), matrizCosteAsociada(matrizCoste)
    {
        this->preparaMatriz(idPadre);
        this->L = this->reduceMatriz();
        this->C += this->L;
    }

    // Devuelve el coste estimado del nodo
    double getCosteEstimado() const   {
        return costeCamino + C;
    }

    // Reduce la matriz de coste asociada y devuelve L
    double reduceMatriz() {
        double aux, total = 0.0;

#ifdef DEBUG_REDUCCION
        cout << "Reduccion matriz: " << flush;
        for(auto &i : this->camino) cout << i << " ->";
        cout << endl;
        printMatriz(this->matrizCosteAsociada);
#endif

        // Reducción filas
        for(int i = 0; i < matrizCosteAsociada.size(); i++) {
            aux = DBL_MAX;
            // Busca el minimo valor
            for(int j = 0; j < matrizCosteAsociada.size() /*Es cuadrada*/; j++) {
                if(matrizCosteAsociada[i][j] < aux) aux = matrizCosteAsociada[i][j];
            }
            // Si se ha encontrado algo, reduce toda la fila y guarda el valor
            if(aux > 0.0 && aux < DBL_MAX) {
                total += aux;
                for (int j = 0; j < matrizCosteAsociada.size() /*Es cuadrada*/; j++) {
                    if(matrizCosteAsociada[i][j] < DBL_MAX) matrizCosteAsociada[i][j] -= aux;
                }
#ifdef DEBUG_REDUCCION
                cout << "Reduccion fila " << i << endl;
                printMatriz(this->matrizCosteAsociada);
#endif
            }
        }

        // Reduccion columnas
        for(int i = 0; i < matrizCosteAsociada.size(); i++) {
            aux = DBL_MAX;
            // Busca el minimo valor
            for(int j = 0; j < matrizCosteAsociada.size() /*Es cuadrada*/; j++) {
                if(matrizCosteAsociada[j][i] < aux) aux = matrizCosteAsociada[j][i];
            }
            // Si se ha encontrado algo, reduce toda la fila y guarda el valor
            if(aux > 0.0 && aux < DBL_MAX) {
                total += aux;
                for (int j = 0; j < matrizCosteAsociada.size() /*Es cuadrada*/; j++) {
                    if(matrizCosteAsociada[j][i] < DBL_MAX) matrizCosteAsociada[j][i] -= aux;
                }
#ifdef DEBUG_REDUCCION
                cout << "Reduccion columna " << i << endl;
                printMatriz(this->matrizCosteAsociada);
#endif
            }
        }

        // Devuelve el total reducido
        return total;
    }

    // Pone a infinito la fila del padre, la columna del nodo y el arco del nodo al padre
    void preparaMatriz(int idPadre)    {
        // Invalida los arcos a nodos anteriores para evitar ciclos
        for(auto &nodoAnterior : this->camino)
            matrizCosteAsociada[id][nodoAnterior] = DBL_MAX;

        // Invalida la fila del padre
        for(int j = 0; j < matrizCosteAsociada.size(); j++)
            matrizCosteAsociada[idPadre][j] = DBL_MAX;

        // Invalida la columna del nodo
        for (int i = 0; i < matrizCosteAsociada.size(); i++)
            matrizCosteAsociada[i][id] = DBL_MAX;
    }

    // Genera los hijos del nodo
    std::vector<Nodo> expande() const {
        // Vector de hijos
        std::vector<Nodo> hijos;

        for(int i = 0; i < matrizCosteAsociada.size(); i++) {
            // Si existe realmente un arco, genera el nodo hijo con arco (id,i) y lo añade al vector de hijos
            if(matrizCosteAsociada[this->id][i] < DBL_MAX)  {
                vector<int> caminoHijo = this->camino;
                caminoHijo.push_back(i);
                Nodo hijo(i, this->id, this->costeCamino + matrizCosteAsociada[this->id][i], this->C, caminoHijo, this->matrizCosteAsociada);
                hijos.push_back(hijo);
            }
        }
        // Devuelve los hijos generados
        return hijos;
    }

};

// Operadores de comparación para la clase Nodo
bool operator> (const Nodo &n1, const Nodo &n2) {
    return n1.getCosteEstimado() > n2.getCosteEstimado();
}

std::vector<int> TSPRamifPoda(std::vector<std::vector<double>>& matrizCoste, double &mejorCoste)   {
    // Nodo inicial
    // - El constructor ya se encarga de reducir la matriz y calcular L
    Nodo init(0, 0.0, vector<int>(1,0), matrizCoste);

    // Mejor coste estimado encontrado
    double mejorCosteEstimado = DBL_MAX;
    vector<int> mejorCamino;

    // El algoritmo se basa en el empled de una cola de prioridad
    std::priority_queue<Nodo, vector<Nodo>, std::greater<>> frontera;

    // Guarda el nodo inicial en la cola y comienza el algoritmo
    frontera.push(init);
    bool sigue = true;
    while(! frontera.empty() && sigue)   {
        // Saca el nodo con el mejor coste estimado
        Nodo mejorNodo = frontera.top(); frontera.pop();

        /*cout << "Saca nodo de frontera: ";
        for(auto &i : mejorNodo.camino)  {
            cout << i << "->";
        } cout << endl;
        */
        // Si saca una nodo que es peor que una hoja encontrada el algoritmo termina
        if(mejorNodo.C >= mejorCosteEstimado)   {
            sigue = false;
        }
        else    {
            // Si es un nodo final (hoja), es decir, ha pasado por todas las ciudades...
            if(mejorNodo.camino.size() == matrizCoste.size())   {
                //cout << "Nodo hoja encontrado (" << mejorNodo.id << ")" << endl;
                // .. si su coste estimado es menor que el actual, lo guarda como nuevo mejor coste
                if(mejorCosteEstimado > mejorNodo.getCosteEstimado()) {
                    mejorCosteEstimado = mejorNodo.getCosteEstimado();
                    mejorCamino = mejorNodo.camino;
                }
            }
                // Si no es nodo final, lo expande y guarda los hijos en la frontera.
            else    {
                //cout << "Expande nodo " << flush;
                //for(auto &i : mejorNodo.camino) cout << i << "->";
                //cout << " con coste estimado = " << mejorNodo.getCosteEstimado() << " (costeCamino="<<
                //    mejorNodo.costeCamino <<" C=" << mejorNodo.C << ")" << endl;
                //printMatriz(mejorNodo.matrizCosteAsociada);
                for(auto &i : mejorNodo.expande())  {
                    //cout << "Anyade nodo a la frontera " << i.id << endl;
                    frontera.push(i);
                }
            }
        }

        mejorCoste = mejorCosteEstimado;
    }

    return mejorCamino;
}

/**                                                 *****************************
 *                                                  * fin ramificación y poda   *
 *                                                  *****************************
 */

constexpr unsigned int str2int(const char* str, int h = 0)
{
    return !str[h] ? 5381 : (str2int(str, h+1) * 33) ^ str[h];
}

int main(int argc, char **argv) {
    if(argc != 3)   {
        cerr << "Modo de empleo: tsp -opt <nombre_fichero>" << endl;
        cerr << "Opciones disponibles:" << endl;
        cerr << "-fb  (Fuerza bruta)" << endl;
        cerr << "-aV  (Algoritmo voraz)" << endl;
        cerr << "-pd  (Programación dinámica)" << endl;
        cerr << "-rp  (Ramificación y poda)" << endl;
        cerr << "-all (Todos)" << endl;
        return -1;
    }
    else {

        string fichero(argv[2]);
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
        auto tInit = chrono::high_resolution_clock::now();
        auto tEnd = chrono::high_resolution_clock::now();
        auto ms_double = tEnd - tInit;
        std::vector<int> mejorCamino;


        switch (str2int(argv[1])) {
            case str2int("-fb") :
                cout << "Fuerza bruta" << endl;
                tInit = chrono::high_resolution_clock::now();
                ejecucionFB(filas, m, mejorCamino);
                tEnd = chrono::high_resolution_clock::now();
                ms_double = tEnd - tInit;
                cout << "Execution time: " << ms_double.count() << "ns" << endl;
                cout << endl << endl;
                break;
            case str2int("-av") :
                //asigno -1 en el indice del recorrido
                //std::vector<int> mejorCamino;
                // double costeMinimo = obtenMejorDynamic(m, filas, mejorCamino);
                // AV
                cout << "Algoritmo voraz" << endl;
                tInit = chrono::high_resolution_clock::now();
                ejecucionAV(filas, m, mejorCamino);
                tEnd = chrono::high_resolution_clock::now();
                ms_double = tEnd - tInit;
                cout  << "Execution time: " << ms_double.count() << "ns" << endl;
                cout << endl << endl;

                /**
                   * fin algoritmo voraz
                */
                break;
            case str2int("-pd") :
                /**
               * Implementacion de programacion dinamica
               */
                cout << "Programacion dinamica" << endl;
                tInit = chrono::high_resolution_clock::now();
                ejecucionPDinamica(m);
                tEnd = chrono::high_resolution_clock::now();
                ms_double = tEnd - tInit;
                cout << "Execution time: " << ms_double.count() << "ns" << endl;
                cout << endl << endl;

                /**
                 * fin Implementacion de programacion dinamica
                 */
                break;
            case str2int("-rp") :
                /**
                * Implementacion Ramificacion y poda
                */

                tInit = chrono::high_resolution_clock::now();
                cout << "Ramificacion y poda" << endl;
                ejecucionRamificacionYpoda(m);
                tEnd = chrono::high_resolution_clock::now();
                ms_double = tEnd - tInit;
                cout << "Execution time: " << ms_double.count() << "ns" << endl;
                cout << endl << endl;

                /**
                 * fin Implementacion de programacion dinamica
                 */
                break;
            case str2int("-all") :
                //asigno -1 en el indice del recorrido
                //std::vector<int> mejorCamino;
                // double costeMinimo = obtenMejorDynamic(m, filas, mejorCamino);
                // AV
                cout << "Algoritmo voraz" << endl;
                tInit = chrono::high_resolution_clock::now();
                ejecucionAV(filas, m, mejorCamino);
                tEnd = chrono::high_resolution_clock::now();
                ms_double = tEnd - tInit;
                cout << "Execution time: " << ms_double.count() << "ns" << endl;
                cout << endl << endl;

                /**
                   * fin algoritmo voraz
                */
                /**
               * fuerza bruta
               */
                // FB
                cout << "Fuerza bruta" << endl;
                tInit = chrono::high_resolution_clock::now();
                ejecucionFB(filas, m, mejorCamino);
                tEnd = chrono::high_resolution_clock::now();
                ms_double = tEnd - tInit;
                cout << "Execution time: " << ms_double.count() << "ns" << endl;
                cout << endl << endl;

                /**
                 * fin fuerza bruta
                 */

                /**
                 * Implementacion de programacion dinamica
                 */
                cout << "Programacion dinamica" << endl;
                tInit = chrono::high_resolution_clock::now();
                ejecucionPDinamica(m);
                tEnd = chrono::high_resolution_clock::now();
                ms_double = tEnd - tInit;
                cout << "Execution time: " << ms_double.count() << "ns" << endl;
                cout << endl << endl;

                /**
                 * fin Implementacion de programacion dinamica
                 */

                /**
             * Implementacion Ramificacion y poda
             */
                cout << "Ramificacion y poda" << endl;
                tInit = chrono::high_resolution_clock::now();
                ejecucionRamificacionYpoda(m);
                tEnd = chrono::high_resolution_clock::now();
                ms_double = tEnd - tInit;
                cout << "Execution time: " << ms_double.count() << "ns" << endl;
                cout << endl << endl;

                /**
                 * fin Implementacion de programacion dinamica
                 */
                //fin programacion dinamica
                break;

        }

        return 0;
    }
}

void ejecucionRamificacionYpoda(const vector<std::vector<double>> &m) {
    double mejorCoste;
    vector<vector<double>> mCoste = m;
    vector<int> solucion;
    solucion = TSPRamifPoda(mCoste, mejorCoste);

    for (auto &i : solucion) {
        cout << i << " -> ";
    }
    cout << "0" << endl;

    cout << "coste final "<< mejorCoste << endl;
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
    cout <<"Coste final: " << dynamicG(0, S, gtab, cities, conjuntosVisitados, vertices) << "\n";
    
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



