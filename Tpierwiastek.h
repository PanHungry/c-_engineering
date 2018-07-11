#include <iostream>
#include <fstream>

using namespace std;

class Tpierwiastek
{
    string symbol;
    double radius;
    double pauling_electronegativity;
    int VEC;
    double ulamek_molowy;

public:
    Tpierwiastek(string nazwa, double promien, double pauling, double vec, double ci)
    {
        symbol = nazwa;
        radius = promien;
        pauling_electronegativity = pauling;
        VEC = vec;
        ulamek_molowy = ci;
    }
};
