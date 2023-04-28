#include<iostream>
#include<fstream>
#include<cmath>
#include<math.h>
#include<vector>
#include<array>

using namespace std;

double num(double h, double phi1/*yi*/, double phi2 /*yi-1*/, double e, double c1 /*ci+1*/, double c2 /*ci*/, double c3/*ci-1*/);

int main ()
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

    cout << "Inserire la largezza della buca L: " << endl;
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

    double h=L/N;
    int Nener=(Emax-Emin)/dE;
    if(Nener*h<(Emax-Emin)){Nener++;}

    vector<double> y1;
    vector<double> y2;
    vector<double> y3;
    vector<double> v1;
    vector<double> v2;
    vector<double> v3;

    // Iniziamo definendo su ogni punto della griglia il potenziale corriespondente

    // Barriera piatta 

    for (int i = 0; i < N; i++)
    {
        v1[i]=0.;
        cout << v1[i] << endl;
    }
    
    // Gradino di potenziale

    for (int i = 0; i < N; i++)
    {
        
      if (i*h>=a && i*h<=b)
      {
            
            v2[i]=V0;

      }
      
      else{  v2[i]=0.; }
        
    }
    
    // Potenziale armonico

    for (int i = 0; i < N; i++)
    {
    
        v3[i]=k*pow(((h*i)-(L/2.)), 2);
    
    }
    
    // Creiamo un vettore in cui mettere tutti i valori di energia che scansioneremo

    vector<double> energie;

    for (int i = 0; i < Nener; i++)
    {
    
        energie[i] = Emin + (i * dE);

        cout << energie[i] << endl;
    
    }

    // Definisco un vettore le cui componenti saranno l' N-esimo elemento di y

    vector<double> y1N;
    vector<double> y2N;
    vector<double> y3N;

    // Inizializzando i primi elementi dei vettori possiamo far partire il ciclo da 2

    y1[0]=0.;
    y2[0]=0.;
    y3[0]=0.;
    y1[1]=1.;
    y2[1]=1.;
    y3[1]=1.;

    /*for (int i = 0; i < Nener; i++)
    {

            for (int j = 2; j < N; j++)
            {
                
                double y1[j]=num(h, y1[j-1], y1[j-2], energie[i], v1[j], v1[j-1], v1[j-2]);
                double y2[j]=num(h, y2[j-1], y2[j-2], energie[i], v2[j], v2[j-1], v2[j-2]);
                double y3[j]=num(h, y3[j-1], y3[j-2], energie[i], v3[j], v3[j-1], v3[j-2]);

            }

        y1N[i]=y1[N];
        y2N[i]=y2[N];
        y3N[i]=y3[N];

        cout << y1N[i] << " " << y2N[i] << " " << y3N[i] << endl;

    }*/
    
    

return 0;

}

double num(double h, double phi1/*yi*/, double phi2 /*yi-1*/, double e, double c1 /*ci+1*/, double c2 /*ci*/, double c3/*ci-1*/){

        double p=(1./(12.-(e-c1)*pow(h, 2.))) * (phi1*(10.*pow(h, 2.)*(e-c2)+24.)+phi2*((e-c3)*pow(h, 2.)-12.));

        return p;    
    
}  
