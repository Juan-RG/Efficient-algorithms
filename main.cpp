#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <algorithm>
#include <vector>

using namespace std;
void rellenarMatriz(int ciudades[100][100], int *fila);

int caminoMinimo(int ciudades[100][100], vector<vector<int>> recorridos, int *indiceMinCamino);

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





int main() {

    // Mirar si se puede mejorar mas fijo
    int ciudades[100][100];                     //toDo: Mirar a ver si podemos hacerlo de alguna otra forma
    int filas;
    rellenarMatriz(ciudades, &filas);
    vector<vector<int>> caminos = generaCaminos(filas);

   //asigno -1 en el indice del recorrido
    int camino = -1;
    int nMin = caminoMinimo(ciudades, caminos, &camino);

    for (int i = 0; i < caminos.at(camino).size(); ++i) {
        cout << caminos.at(camino).at(i) << " - ";
    }
    cout << caminos.at(camino).at(0);
    cout<<"\n";
    cout << "coste: " << nMin << "\n";
    return 0;


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

void rellenarMatriz(int ciudades[100][100], int *tamanio){    //toDo: Falta hacer que se pueda pasar el fichero por arv
    ifstream datos("datos");
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
    *tamanio = filas;
}