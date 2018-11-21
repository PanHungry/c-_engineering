#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <sstream>
#include <algorithm>
#include <string>
#include <cstdlib>

using namespace std;

// Vectory denych wejściowych pierwiastków. Ładowanie danych jest pierwszą funkcją przy
// włączaniu programu z pliku dane_metali.txt

vector <string> symbol;
vector <double> promien;
vector <double> pauling;
vector <double> temperatura_m;
vector <double> vec;
vector <double> ulamek_molowy;
vector <int> liczba_atomowa;

vector <string> wyniki;

// Tablica entalpi par. Wpisywanie do tablicy jest drugą funkcją przy
// włączeniu programu z pliku entalpy_mix.txt

double tablica_entalpi_par[91][91];

void wpisuje_dane_pierwiastkow_do_vectorow() //funkcja wczytuje wartosci pierwiastkow do vectorów
{
    double suma_stechiometria=0;

        string symbol_atomu;
        int liczba_atomowa_atomu;
        double promien_atomu;
        double pauling_atomu;
        double vec_atomu;
        double tm_atomu;

        fstream plik;
        plik.open ("dane_metali.txt", ios::in);
        if(plik.good() == true)
        {
        for (int i=0; i<=92; i++)
        {
        plik >> symbol_atomu >> liczba_atomowa_atomu >> promien_atomu >> pauling_atomu >>vec_atomu >> tm_atomu;

        liczba_atomowa.push_back(liczba_atomowa_atomu);
        promien.push_back(promien_atomu);
        pauling.push_back(pauling_atomu);
        vec.push_back(vec_atomu);
        symbol.push_back(symbol_atomu);
        temperatura_m.push_back(tm_atomu);

        }

        plik.close();

        }
    else cout << "Nie udalo sie otworzyc pliku z danymi pierwiastkow!" << endl;
}

double licze_delta_chi(unsigned int ile,double ulamki[], double paulingi[])
{
   double suma_psi=0;

    for (int i=0; i<ile; i++)
    {
        suma_psi = suma_psi + paulingi[i]; //do obliczenia sredniej elektroujemnosci
    }

    double srednie_psi = suma_psi/ile; // srednia elektroujemnosc

    //////////////////////////////////////////////////////////////////////////////////

    double suma_pod_pierwiastkiem=0;

    for (int i=0; i<ile; i++)
    {
        suma_pod_pierwiastkiem = suma_pod_pierwiastkiem + ulamki[i]*pow((paulingi[i]-srednie_psi),2) ; // funckja pow(liczba,potega) podnosi liczbe do podanej potegi
    }

    return sqrt(suma_pod_pierwiastkiem);

}

double licze_VEC(unsigned int ile,double ulamki[],double vec[])
{
    double suma_ci_VECi=0;

   for (int i=0; i<ile; i++)
   {
      double ci_VECi=0;
      ci_VECi = ulamki[i] * vec[i];

      suma_ci_VECi = suma_ci_VECi + ci_VECi;
   }

   return suma_ci_VECi;
}

double licze_mala_delta(unsigned int ile, double ulamki[], double promien[])
{
    double suma_promieni=0;

    for (int i=0; i<ile; i++)
    {
        suma_promieni = suma_promieni + promien[i]; //do obliczenia sredniego r
    }

    double sredni_promien = suma_promieni/ile; // sredni promien ri

    ///////////////////////////////////////////////////

    double suma_pod_pierwiastkiem=0;

    for (int i=0; i<ile; i++)
    {
        suma_pod_pierwiastkiem = suma_pod_pierwiastkiem + ulamki[i]*pow(1-(promien[i]/sredni_promien),2) ; // funckja pow(liczba,potega) podnosi liczbe do podanej potegi
    }

    return 100*sqrt(suma_pod_pierwiastkiem);

}

double licze_delta_s(unsigned int ile, double ulamki[])
{
    const double R = 8.3144598;
    double suma=0;

    for (int i=0; i<ile; i++)
    {
        suma = suma + ulamki[i] * log(ulamki[i]);
    }

    return -1*R*suma;
}

double licze_temp_m(unsigned int ile, double ulamki[], double tm[]) // << do zrobienia, trzeba wpisac dane temperatur topnienia
{
   double suma_ci_tm=0;

   for (int i=0; i<ile; i++)
   {
      double ci_tmi=0;
      ci_tmi= tm[i] * ulamki[i];

      suma_ci_tm = suma_ci_tm + ci_tmi;
   }

   return suma_ci_tm;
}

double licze_entalpia_mix(unsigned int ile, double ulamki[], int liczby_atomowe_w_linijce[])
{

double entalpia_mix_pary=0;
double entalpia_mix_stopu=0;

    for (int i=0; i<ile-1; i++)
    {
    for (int j=i+1; j<=ile-1; j++)
        {
            entalpia_mix_pary = ulamki[i] * ulamki[j] * tablica_entalpi_par[liczby_atomowe_w_linijce[i]][liczby_atomowe_w_linijce[j]];
            entalpia_mix_stopu = entalpia_mix_stopu + entalpia_mix_pary;
        }
    }
    return 4*entalpia_mix_stopu;
}

void wczytuje_tabele_entalpi_par() // wczytuje dane entalpi par do tablicy tablica_entalpi_par
{
    fstream entalpie;
    entalpie.open("enthalpy_table.txt", ios::in);
    if(entalpie.good() == true)
    {
        for (int i=0; i<91; i++)
        {
            for (int j=0; j<91; j++)
            {
                entalpie >> tablica_entalpi_par[i][j];
            }
        }
    }

}

void wczytuje_liste_do_przeliczenia() // wczytuje liste do przeliczenia i wysyła linijkę do funkcji liczących
{

    fstream lista;
        lista.open ("lista_do_przeliczenia.txt", ios::in);
        if(lista.good() == true)
        {
            vector <string> linijka_do_przeliczenia_podzielona_spacjami;
            vector <string> symbole_w_linijce;
            vector <string> stechiometria_w_linijce_string;
        {
            while (!lista.eof())
            {
            string surowa_linijka_cala;

            string surowa_linijka;
            string references;
            string struktura;
            string processing_condition;

            getline(lista,surowa_linijka_cala);

            cout << "Pobrana linijka: " <<surowa_linijka_cala << endl;

            // Poniżej dzielę linijkę na: wzór, numer pracy, budowe...

            stringstream string_calej_linijki_w_strumien(surowa_linijka_cala);
            string_calej_linijki_w_strumien >> surowa_linijka >> processing_condition >> struktura >> references;

            wyniki.push_back(surowa_linijka);


            cout << "Processing Condition: " << processing_condition <<endl;
            cout << "Struktura: " <<struktura<<endl;
            cout << "References: " << references <<endl;


            //Poniżej przerobienie pobranej linijki bez spacji na taką ze spacjami i przyjującą 1 jako podstawową stechiometrię.

            int i=0;
            for (int i=0; i<surowa_linijka.size(); i++)
            {
                if(surowa_linijka[i+1]=='\0' && surowa_linijka[i]>=65)
                {
                    surowa_linijka += (" 1");
                }
                else if(surowa_linijka[i]>=97 && surowa_linijka[i+1]<58)
                {
                    surowa_linijka.insert(i+1, " ");
                }
                else if(surowa_linijka[i]<=57 && surowa_linijka[i]>=48 && surowa_linijka[i+1]>65)
                {
                    surowa_linijka.insert(i+1, " ");
                }
                else if(surowa_linijka[i]>=65 && surowa_linijka[i]<=90 && surowa_linijka[i+1]>=65 && surowa_linijka[i+1]<=90)
                {
                    surowa_linijka.insert(i+1, " 1 ");
                }
                else if(surowa_linijka[i]>=97 && surowa_linijka[i]<=122 && surowa_linijka[i+1]>=65 && surowa_linijka[i+1]<=90)
                {
                    surowa_linijka.insert(i+1, " 1 ");
                }
                else if(surowa_linijka[i]>=65 && surowa_linijka[i]<=90 && surowa_linijka[i+1]>=48 && surowa_linijka[i+1]<=57)
                {
                    surowa_linijka.insert(i+1, " ");
                }

            }

            cout <<"Po obrobce: " << surowa_linijka <<endl<<endl;

            //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

            stringstream ss(surowa_linijka);
            string symbol_lub_stechio;

            /*
            Poniżej podzielenie linijki na vectory symboli i stechiometrii
            */

            int parzystosc=0;
            while (getline(ss,symbol_lub_stechio, ' '))
            {
                    linijka_do_przeliczenia_podzielona_spacjami.push_back(symbol_lub_stechio);
            }

             for (int i=0; i<linijka_do_przeliczenia_podzielona_spacjami.size(); i=i+2)
            {
                symbole_w_linijce.push_back(linijka_do_przeliczenia_podzielona_spacjami[i]);
            }

             for (int i=1; i<linijka_do_przeliczenia_podzielona_spacjami.size(); i=i+2)
            {
                stechiometria_w_linijce_string.push_back(linijka_do_przeliczenia_podzielona_spacjami[i]);
            }

            vector <double> stechiometria_w_linijce_double(stechiometria_w_linijce_string.size());
            transform(stechiometria_w_linijce_string.begin(), stechiometria_w_linijce_string.end(),
                      stechiometria_w_linijce_double.begin(), [](string const &val) {return stod(val);});

            ///// obliczanie ułamka molowego ///////

            double suma_stechiometri=0;
            double ulamki_molowe_w_linijce[stechiometria_w_linijce_double.size()] {};

            for (int i=0; i<stechiometria_w_linijce_double.size(); i++)
            {
                suma_stechiometri = suma_stechiometri + stechiometria_w_linijce_double[i];
            }

            for (int i=0; i<stechiometria_w_linijce_double.size(); i++)
            {
                ulamki_molowe_w_linijce[i]=(stechiometria_w_linijce_double[i]/suma_stechiometri);
            }

            for (int i=0; i<stechiometria_w_linijce_double.size(); i++)
            {
            //cout << "Ulamek molowy skladnika " << i+1 << " wynosi: " << ulamki_molowe_w_linijce[i]<<endl;
            }

            ///// wczytywanie danych metalu z tablicy danych wyjsciowych //////

            int tablica_liczb_atomowych_w_linijce[symbole_w_linijce.size()] {};
            double tablica_promieni_atomow_w_linijce[symbole_w_linijce.size()] {};
            double tablica_paulingow_w_linijce[symbole_w_linijce.size()] {};
            double tablica_vec_w_linijce[symbole_w_linijce.size()] {};
            double tablica_tm_w_linijce[symbole_w_linijce.size()] {};

            for (int i=0; i<symbole_w_linijce.size(); i++)
            {
                for (int j=0; j<92; j++)
                {
                    if (symbole_w_linijce[i]==symbol[j])
                    {
                        tablica_liczb_atomowych_w_linijce[i]=liczba_atomowa[j];
                        tablica_promieni_atomow_w_linijce[i]=promien[j];
                        tablica_paulingow_w_linijce[i]=pauling[j];
                        tablica_vec_w_linijce[i]=vec[j];
                        tablica_tm_w_linijce[i]=temperatura_m[j];
                    }
                }
            }


            ///// WYNIKI

            wyniki.push_back(to_string(licze_VEC(symbole_w_linijce.size(), ulamki_molowe_w_linijce, tablica_vec_w_linijce)));
            wyniki.push_back(to_string(licze_mala_delta(symbole_w_linijce.size(), ulamki_molowe_w_linijce, tablica_promieni_atomow_w_linijce)));
            wyniki.push_back(to_string(licze_delta_chi(symbole_w_linijce.size(), ulamki_molowe_w_linijce, tablica_paulingow_w_linijce)));
            wyniki.push_back(to_string(licze_delta_s(symbole_w_linijce.size(), ulamki_molowe_w_linijce)));
            wyniki.push_back(to_string(licze_temp_m(symbole_w_linijce.size(), ulamki_molowe_w_linijce, tablica_tm_w_linijce)));
            wyniki.push_back(to_string(licze_entalpia_mix(symbole_w_linijce.size(), ulamki_molowe_w_linijce, tablica_liczb_atomowych_w_linijce)));
            wyniki.push_back(processing_condition);
            wyniki.push_back(struktura);
            wyniki.push_back(references);

            /*
            cout <<endl << "VEC linijki: " << licze_VEC(symbole_w_linijce.size(), ulamki_molowe_w_linijce, tablica_vec_w_linijce)<<endl;
            cout << "Mala delta linijki: " << licze_mala_delta(symbole_w_linijce.size(), ulamki_molowe_w_linijce, tablica_promieni_atomow_w_linijce)<<endl;
            cout << "Delta CHI: " <<licze_delta_chi(symbole_w_linijce.size(), ulamki_molowe_w_linijce, tablica_paulingow_w_linijce)<<endl;
            cout << "Delta S: " <<licze_delta_s(symbole_w_linijce.size(), ulamki_molowe_w_linijce)<<endl;
            cout << "Temperatura topnienia: "<<licze_temp_m(symbole_w_linijce.size(), ulamki_molowe_w_linijce, tablica_tm_w_linijce)<<endl;
            cout << "Entalpie par: "<< licze_entalpia_mix(symbole_w_linijce.size(), ulamki_molowe_w_linijce, tablica_liczb_atomowych_w_linijce) <<endl;

            cout << ".................................................................................................." <<endl<<endl;

            */

            //Sprzątanie po linijce, bo zaraz wczytujemy nową...

            linijka_do_przeliczenia_podzielona_spacjami.clear();
            symbole_w_linijce.clear();
            stechiometria_w_linijce_string.clear();
            stechiometria_w_linijce_double.clear();
            }

        }
}

}

void zapis_do_pliku()
{
    fstream result( "wyniki.txt", ios::out );
    if( result.good() )
    {
        for( int i = 0; i <= wyniki.size(); i++ )
        {
            result << wyniki[i] << '\t' ;

            if(i%10==9)
            {
               result <<endl;
            }

            result.flush();
        }
        result.close();
    }
}

int main()
{
    wpisuje_dane_pierwiastkow_do_vectorow(); // funkcja wczytujaca z plikow promienie, paulinga, vec i licząca ulamek molowy
    wczytuje_tabele_entalpi_par(); // funckja wczytujaca z plikow entalpie par
    wczytuje_liste_do_przeliczenia(); // funkcja wczytujaca liste zwiazkow do obliczenia
    zapis_do_pliku();


    /*

    dodać omege

    pisanie:
    1. plan ramowy (rozdzialy i podrozdzialy)

    1. stopy o wysokiej entropii
    1.1. opis termodynamiczny t-deltas itd TERMODYNAMIKA HEA
    1.2 kryteria empiryczne kto zaproponowal co wynika jakie sa dotychczas spojrzenie na werości KRYTERIA EMPIRYCZNE
    1.3. dlatego stopy istotne, krotki przeflad HEA , 4 core effects WŁAŚCIWOSCI HEA
    2. Cele pracy (weryfukacji kryteriow empirycznych, teraz globalne podejscie i sprawdzenie)
    3. Metodologia
    - ramowy plan, schemat blokowy
    - funkcje programu, jak licza
    4. wyniki
    5. dyskusja (jak sie ma kryteria do tego co zrobilismy)
    bez aneksu 40-60 stron
    literatura kilkadziesiąt pozycji, bo praca przeglądowa
    6. aneks



    WYKRESY:
    mala delta(omega)
    mala delta(Hmix)
    mala delta(VEC)
    z podziałem na FCC i BCC,


    */

}

