#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <time.h>
#include <string.h>
#include <complex>
#include <string>
using namespace std;


double **Makematriz(int fila,int columna)
{
double **matriz;
matriz = new double* [fila];

for(int j=0; j<fila; j++)
    {
        matriz[j] = new double [columna];
    }


return (matriz);
}

void Releaseme (double **v,int m)
{
   for(int t=0;t<m;t++)
      delete [] v[t];

   delete [] v;
}

int main(int argc, char *argv[])
{
  
 
 
    
double entrance;

ifstream fo("Tiempo_disociacion.txt");

double *numeros;
numeros=new double [5000000];
int contador=0;
double sumador=0;

while(fo>>entrance)
{

              ++contador;
              numeros[contador]=entrance;
              sumador=sumador+1;
       }
cout<<contador<<endl;

double promedio;

for(int i=1;i<contador+1;i++)
{

fo>>numeros[i];
//cout<<numeros[i]<<endl;
promedio=promedio+numeros[i];


}

promedio=promedio/sumador;
cout<<" Tiempo de vida promedio es="<<promedio<<endl;

    return 0;
}
