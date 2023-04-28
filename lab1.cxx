#include<iostream>
#include<fstream>
#include<math.h>
#include<string>
#include<cmath>
using namespace std;
int main () {

	

	int n;
	double h;
	double x0;
	
	double a,b;
	double S; 
	ofstream prova("graph.txt"); //creato un file di testo
	
	cout << "Quale programma vuoi usare:\n 1) Griglia di punti equispaziati con offset iniziale;\n 2) Griglia di punti in intervallo fissato." << endl;
	cin >> S;
	
	if(S==1){
		
		cout << "Inserire il numero di punti desiderato " << endl;
		cin >> n;
		cout << "Inserire il punto da cui partire " << endl;
		cin >> x0; 
		cout << "Inserire la spaziatura " << endl;
		cin >> h;
		double m=0; 
		double s=0;
		
		for(int i=0;i<n;i++){
		
			double m=x0+h*i;
			double s=sin(m);
			cout << m << "  " << s << endl;
			prova << m << "  "  << s << endl; //scrivo i valori calcolati nel file di testo
		
		}
	} 
	
	else if (S==2) { 
		
		cout << "Inserire il numero di punti " << endl;
		cin >> n;
		cout << "Inserire gli estremi dell'intervallo " << endl;
		cin >> a >> b;
		double m=0;
		double s=0;
	
		for (int i=0; i<=n; i++){
		
			double m=((b-a)/n)*i+a;
			double s= sin(m);
			cout << m << "  " << s  << endl;
			prova << m << "  " << s << endl;
			}
	
	
	prova.close();   //chiudo il file di testo
			
	}
	else {
		cout << "numero non valido" << endl;
	}

return 0;
}
