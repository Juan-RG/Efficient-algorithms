#include <iostream>
#include <fstream>
#include <random>
#include <iomanip>

using namespace std;

int MIN = 0;
int MAX = 100;

int main() {
    random_device rd;
    default_random_engine realAleatorio(rd());
    uniform_real_distribution<double> generador(MIN, MAX);



    int tamanio;
    cout << "tamanio x tamanio\n";
    cout << "tamanio tsp: \n";
    cin >> tamanio;

    ofstream fichero;

    string nombre = "a" + to_string(tamanio) + ".tsp";
    fichero.open(nombre);

    double numero;
    for (int i = 0; i < tamanio; ++i) {
        for (int j = 0; j < tamanio; ++j) {
            if (i == j){
                numero = 0;
            }else{
                numero = generador(realAleatorio);
            }
            fichero << setprecision(3) << numero << " ";
        }
        fichero << "\n";
    }

    return 0;
}
