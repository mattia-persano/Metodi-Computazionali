#include<iostream>
#include<fstream>
#include<cmath>
#include<math.h>
#include<vector>
#include<array>

using namespace std;

void num(int N, double h, double *y, double e, double *c);

int main()
{

    int N;
    double L;
    double Emin;
    double Emax;
    double dE;
    double s;
    double a;
    double b;
    double V0;
    double k;

    cout << "Inserire la larghezza della buca L: " << endl;
    cin >> L;
    cout << "Inserire la spaziatura della buca N: " << endl;
    cin >> N;
    cout << "Inserire gli estremi dell'energia Emin e Emax e la spaziatura dE: " << endl;
    cin >> Emin >> Emax >> dE;
    cout << "Inserire gli estremi del gradino di potenziale a e b e il valore V0: " << endl;
    cin >> a >> b >> V0;
    cout << "Inserire la costante del potenziale armonico k: " << endl;
    cin >> k;
    cout << "Inserire il parametro di soglia s: " << endl;
    cin >> s;

    double h = L / N;
    int Nener = (Emax - Emin) / dE;
    if (Nener * h < (Emax - Emin))
    {
        Nener++;
    }

    // Creo array dinamici per le funzioni d'onda e i potenziali

    double *y1 = new double[N];
    double *y2 = new double[N];
    double *y3 = new double[N];
    double *v1 = new double[N];
    double *v2 = new double[N];
    double *v3 = new double[N];

    // Iniziamo definendo su ogni punto della griglia il potenziale corrispondente

    // Barriera piatta

    for (int i = 0; i < N; i++)
    {

        v1[i] = 0.;
    }

    // Gradino di potenziale

    for (int i = 0; i < N; i++)
    {

        if (i * h >= a && i * h <= b)
        {

            v2[i] = V0 * (-2.);
        }

        else
        {
            v2[i] = 0.;
        }
    }

    // Potenziale armonico

    for (int i = 0; i < N; i++)
    {

        v3[i] = k * pow(((h * i) - (L / 2.)), 2.) * (-2.);
    }

    // Check point 1: fino a quì tutto bene :D

    // creo un array dinamico per le energie

    double *ener = new double[Nener];

    for (int i = 0; i < Nener; i++)
    {

        ener[i] = (Emin + i * dE) * (-2.);
    }

    // Definisco un array dinamico le cui componenti saranno l' N-esimo elemento di y

    double *y1N = new double[Nener];
    double *y2N = new double[Nener];
    double *y3N = new double[Nener];

    // Inizializzando i primi elementi dei vettori possiamo far partire il ciclo da 2

    y1[0] = 0.;
    y2[0] = 0.;
    y3[0] = 0.;
    y1[1] = 1.;
    y2[1] = 1.;
    y3[1] = 1.;

    // Ora per ogni energia uso il propagatore di Numerov per ottenere le autofunzioni

    for (int i = 0; i <= Nener; i++)
    {

        num(N, h, y1, ener[i], v1);
        num(N, h, y2, ener[i], v2);
        num(N, h, y3, ener[i], v3);

        y1N[i] = y1[N - 1];
        y2N[i] = y2[N - 1];
        y3N[i] = y3[N - 1];
    }

    // Check point 2: fino a quì funonzia :D

    // Tramite metodo delle bisezione, cerco i valori dell'energia

    ofstream autov1("autovalori1.txt");
    ofstream autov2("autovalori2.txt");
    ofstream autov3("autovalori3.txt");

    int N1 = 0;
    int N2 = 0;
    int N3 = 0;

    for (int i = 0; i < Nener; i++)
    {

        // Buca piatta

        if (y1N[i] * y1N[i + 1] <= 0)
        {

            N1 = N1 + 1;

            double a = y1N[i];
            double b = y1N[i + 1];
            double c = y1N[i] + y1N[i + 1];
            double E1 = ener[i];
            double E2 = ener[i + 1];

            double NE;

            do
            {

                NE = (E1 + E2) / 2.;

                num(N, h, y1, NE, v1);

                if (y1[N - 1] * a <= 0)
                {

                    b = y1[N - 1];
                    E2 = NE;
                }

                else if (y1[N - 1] * b <= 0)
                {

                    a = y1[N - 1];
                    E1 = NE;
                }

            } while (abs(E1 - E2) >= s);

            autov1 << N1 << " " << NE * -0.5 << " " << 0.5 * pow(M_PI / L, 2.) * pow(N1, 2.) << endl;
        }

        // Gradino di potenziale

        if (y2N[i] * y2N[i + 1] <= 0)
        {

            N2++;

            double a = y2N[i];
            double b = y2N[i + 1];
            double c = y2N[i] + y2N[i + 1];
            double E1 = ener[i];
            double E2 = ener[i + 1];

            double NE;

            do
            {

                NE = (E1 + E2) / 2.;

                num(N, h, y2, NE, v2);

                if (y2[N - 1] * a <= 0)
                {

                    b = y2[N - 1];
                    E2 = NE;
                }

                else if (y2[N - 1] * b <= 0)
                {

                    a = y2[N - 1];
                    E1 = NE;
                }

            } while (abs(E1 - E2) >= s);

            autov2 << N2 << " " << NE * -0.5 << endl;
        }

        // Potenziale armonico

        if (y3N[i] * y3N[i + 1] <= 0)
        {

            N3++;

            double a = y3N[i];
            double b = y3N[i + 1];
            double c = y3N[i] + y3N[i + 1];
            double E1 = ener[i];
            double E2 = ener[i + 1];

            double NE;

            do
            {

                NE = (E1 + E2) / 2.;

                num(N, h, y3, NE, v3);

                if (y3[N - 1] * a <= 0)
                {

                    b = y3[N - 1];
                    E2 = NE;
                }

                else if (y3[N - 1] * b <= 0)
                {

                    a = y3[N - 1];
                    E1 = NE;
                }

            } while (abs(E1 - E2) >= s);

            autov3 << N3 << " " << NE * -0.5 << endl;
        }
    }

    autov1.close();
    autov2.close();
    autov3.close();

    // Check point 3: fino a qui tutto funonzia :D

    // Trovati gli autovalori, vanno trovate le corrispondenti autofunzioni normalizzate

    // Calcoliamo i coefficienti di normalizzazione imponendo che l'integrale dell'autofunzione sia uguale ad 1

    double *Y1 = new double[N];
    double *Y2 = new double[N];
    double *Y3 = new double[N];

    double W1 = 0.; // Coefficiente di normalizzazione al quadrato
    double W2 = 0.;
    double W3 = 0.;

    // Per calcolare la normalizzazione utilizziamo la regola di Simpson

    W1 = W1 + (h / 3.) * pow(y1[0], 2.);
    W2 = W2 + (h / 3.) * pow(y2[0], 2.);
    W3 = W3 + (h / 3.) * pow(y3[0], 2.);

    W1 = W1 + (h / 3.) * pow(y1[N - 1], 2.);
    W2 = W2 + (h / 3.) * pow(y2[N - 1], 2.);
    W3 = W3 + (h / 3.) * pow(y3[N - 1], 2.);

    for (int i = 1; i < N - 1; i = i + 2)
    {

        W1 = W1 + (4. * h / 3.) * pow(y1[i], 2.);
        W2 = W2 + (4. * h / 3.) * pow(y2[i], 2.);
        W3 = W3 + (4. * h / 3.) * pow(y3[i], 2.);

        W1 = W1 + (2. * h / 3.) * pow(y1[i + 1], 2.);
        W2 = W2 + (2. * h / 3.) * pow(y2[i + 1], 2.);
        W3 = W3 + (2. * h / 3.) * pow(y3[i + 1], 2.);
    }

    W1 = sqrt(W1);
    W2 = sqrt(W2);
    W3 = sqrt(W3);

    // Dividiamo ciascun valore dell'autofunzione per i rispettivi coefficienti di normalizzazione

    for (int i = 0; i < N; i++)
    {

        Y1[i] = y1[i] / W1;
        Y2[i] = y2[i] / W2;
        Y3[i] = y3[i] / W3;
    }

    // Creiamo i file in cui salvare le autofunzioni normalizzate

    ofstream pot1("Barriera.txt");
    ofstream auto1("Autofunzioni1.txt");
    ofstream pot2("Gradino.txt");
    ofstream auto2("Autofunzioni2.txt");
    ofstream pot3("Armonico.txt");
    ofstream auto3("Autofunzioni3.txt");

    for (int i = 0; i < N; i++)
    {

        pot1 << i * h << " " << v1[i] << endl;
        auto1 << i * h << " " << Y1[i] << endl;
        pot2 << i * h << " " << v2[i] << endl;
        auto2 << i * h << " " << Y2[i] << endl;
        pot3 << i * h << " " << v3[i] << endl;
        auto3 << i * h << " " << Y3[i] << endl;
    }

    pot1.close();
    auto1.close();
    pot2.close();
    auto2.close();
    pot3.close();
    auto3.close();

    delete[] y1;
    delete[] y2;
    delete[] y3;
    delete[] v1;
    delete[] v2;
    delete[] v3;
    delete[] ener;
    delete[] y1N;
    delete[] y2N;
    delete[] y3N;
    delete[] Y1;
    delete[] Y2;
    delete[] Y3;

    return 0;
}

void num(int N, double h, double *y, double e, double *c)
{

    double H = h * h;

    for (int i = 2; i < N; i++)
    {

        y[i] = 1. / (12. - (e - c[i]) * H) * ((y[i - 1] * (10. * H * (e - c[i - 1]) + 24.)) + (y[i - 2] * ((e - c[i - 2]) * H - 12.)));
    }
}