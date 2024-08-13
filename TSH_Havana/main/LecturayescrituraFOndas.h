#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <time.h>
#include <string.h>
#include <cstring>
#include <complex>
using namespace std;
double mA,mC,mB,mBC,mABC; ///estan en unidades de masa atomicas
double hbar;  /// Asi esta en las hojas de las constantes de Chuchi, esta en cm-1 * ps

complex<double> **coef;
double *ev; 
double posfondap,nivelesesenciales;
double pasogrilla,pasoinicialgrilla,rfinal,boundbroken;
int numerofondas,nptosfonda,numfondainicial,iteracionesciclo;
char *sacariteraciones;
int cicloiteraciones;
double **fondatotal;



////////////////////////////////////////////////////////////////
complex<double> **Makematriz2(int fila,int columna)
{
complex<double> **matriz;
matriz = new complex<double>* [fila];

for(int j=0; j<fila; j++)
    {
        matriz[j] = new complex<double> [columna];
    }


return (matriz);
}

void Releaseme2 (complex<double> **v,int m)
{
   for(int t=0;t<m;t++)
      delete [] v[t];

   delete [] v;
}

double **Makematriz3(int fila,int columna)
{
double **matriz;
matriz = new double* [fila];

for(int j=0; j<fila; j++)
    {
        matriz[j] = new double [columna];
    }


return (matriz);
}

void Releaseme3 (double **v,int m)
{
   for(int t=0;t<m;t++)
      delete [] v[t];

   delete [] v;
}


//////////////////////////////////////////////////////////////////

void lecturaprimera()
{
sacariteraciones=new char [5];
ifstream parini("parametros_iniciales");
parini>>mA>>mB>>mC>>hbar>>numerofondas>>numfondainicial>>nptosfonda>>pasogrilla>>pasoinicialgrilla>>nivelesesenciales>>boundbroken>>sacariteraciones>>rfinal;

ev=new double [numerofondas+1];

for(int i=0;i<numerofondas;i++)
parini>>ev[i];


cout<<"*****************************************************************************************************************************************|"<<endl;
cout<<"Las energias vibracionales son:"<<endl;

for(int i=0;i<numerofondas;i++)
cout<<"ev["<<i<<"]="<<ev[i]<<" corresponde a fonda"<<numfondainicial+i<<endl;


mBC=mB*mC/(mB+mC);
mABC=mA*(mB+mC)/(mA+mB+mC);

cout<<"\n mA="<<mA<<" | mB="<<mB<<"| mC="<<mC<<" | mBC="<<mBC<<" | mABC="<<mABC<<" | hbar="<<hbar<<" | punto inicial de la grilla="<<pasoinicialgrilla<<" | paso de la grilla="<<pasogrilla<<" | nivel de partida (en caso de utilizar factor de escalado)="<<nivelesesenciales<<endl;
 cicloiteraciones=atof(sacariteraciones);
cout<<"| numero de fondas="<<numerofondas<<"| numero de puntos de fonda="<<nptosfonda<<" | Primera fonda="<<numfondainicial<<" | tiempo de vida maximo="<<boundbroken<<"ps | iteraciones="<<cicloiteraciones<<" | R disociacion="<<rfinal<<"AA"<<endl;
cout<<"*****************************************************************************************************************************************|"<<endl;

fondatotal=Makematriz3(nptosfonda+1,numerofondas+1);


char pp[]="fonda";

FILE *output;
output=fopen("Valores","w+");


for(int i=numfondainicial;i<numerofondas+numfondainicial;i++)
{
 fprintf(output," %i \n",i);
}

fclose(output);

ifstream inpfon("Valores");


char *este;
este= new char [4];
char cambio[8];
strcpy(cambio,pp);
for(int i=0;i<numerofondas;i++)
{

inpfon>>este;
strcat(pp, este);
//cout<<pp<<endl;

fstream inpfonda(pp);

             for(int k=0;k<nptosfonda;k++)
                  {
                   for(int r=1;r<3;r++)
                  {

                   if(r==1)
                     {
                    inpfonda>>posfondap;
                     }else
                     {
                        inpfonda>>fondatotal[k][i];
                        //cout<<fondatotal[k][i]<<endl;
                     }

                    }
                       }



 
strcpy(pp,cambio);

}

}

//******************************************************************************
void normalizafondas()
{

double *probando,*normalizar;
probando=new double [nptosfonda];
normalizar=new double [numerofondas];

for(int ll=numerofondas-1;ll>=0;ll--) /// esta OK
   normalizar[ll]=1;
   

  double h2=pasogrilla,factor=h2/22.5,intf12,factor2;
  for(int pp=numerofondas-1;pp>=0;pp--)
  {
intf12=0;
for(int ll=0;ll<nptosfonda;ll+=4) /// esta OK
  {
intf12=intf12+7*fondatotal[ll][pp]*fondatotal[ll][pp]+32*fondatotal[ll+1][pp]*fondatotal[ll+1][pp]+12*fondatotal[ll+2][pp]*fondatotal[ll+2][pp]+32*fondatotal[ll+3][pp]*fondatotal[ll+3][pp]+7*fondatotal[ll+4][pp]*fondatotal[ll+4][pp];
}
 
 intf12=intf12*factor;
 
cout<<"La funcion de onda "<<pp<<" no esta normalizada, tienes que dividir por la raiz de este factor:= "<<intf12<<endl;
  
normalizar[pp]=sqrt(intf12);
 
 }
 
 for(int j=0;j<numerofondas;j++)
 {
   for(int i=0;i<nptosfonda;i++)
   {
    fondatotal[i][j]=fondatotal[i][j]/normalizar[j];
   }
 }
 
   for(int pp=numerofondas-1;pp>=0;pp--)
   {
   intf12=0;                            
             
    for(int ll=0;ll<nptosfonda;ll+=4) /// esta OK
  {
  intf12=intf12+7*fondatotal[ll][pp]*fondatotal[ll][pp]+32*fondatotal[ll+1][pp]*fondatotal[ll+1][pp]+12*fondatotal[ll+2][pp]*fondatotal[ll+2][pp]+32*fondatotal[ll+3][pp]*fondatotal[ll+3][pp]+7*fondatotal[ll+4][pp]*fondatotal[ll+4][pp];
}
 cout<<" F.onda "<<pp<<" normalizada a= "<<intf12*factor<<endl;
 } 
 
delete [] probando,normalizar;



}


//******************************************************************************

/*
char *mi_strcpy(char *cambio,int opcion)
{


   switch(opcion)
  {
       case 0:
    cambio="fonda10";
    break;
    
       case 1:
    cambio="fonda11";
    break;
    
       case 2:
    cambio="fonda12";
    break;
    
       case 3:
    cambio="fonda13";
    break;
    
       case 4:
    cambio="fonda14";
    break;
    
    
       case 5:
    cambio="fonda15";
    break;
    
    case 6:
     cambio="fonda16";

      break;

      case 7:
      cambio="fonda17";

      break;

      case 8:
        cambio="fonda18";

      break;
       case 9:
     cambio="fonda19";

      break;
       case 10:
      cambio="fonda20";

      break;
       case 11:
    cambio="fonda21";

      break;
       case 12:
       cambio="fonda22";

      break;
       case 13:
      cambio="fonda23";

      break;
       case 14:
     cambio="fonda24";

      break;
       case 15:
    cambio="fonda25";

      break;
       case 16:
      cambio="fonda26";

      break;
       case 17:
     cambio="fonda27";

      break;
       case 18:
      cambio="fonda28";

      break;
       case 19:
    cambio="fonda29";

      break;
        case 20:
    cambio="fonda30";

      break;
        case 21:
    cambio="fonda31";

      break;
        case 22:
    cambio="fonda32";

      break;
      

  }



   return cambio;
}

 inline double choice(double *fonda,int m)
{
   int ptos=288;

   switch(m)
  {
    case 0:

                for(int i=0;i<ptos;i++)
                  {
                  
                  fonda[i]=fo0[i];

                    }

      break;
  
    case 1:

                for(int i=0;i<ptos;i++)
                  {
                  
                  fonda[i]=fo1[i];

                    }

      break;

      case 2:
      for(int i=0;i<ptos;i++)
                  {

                  
                  fonda[i]=fo2[i];

                       }

      break;

      case 3:
        for(int i=0;i<ptos;i++)
                  {

                  
                  fonda[i]=fo3[i];

                       }

      break;

       case 4:
   for(int i=0;i<ptos;i++)
                  {

                  
                  fonda[i]=fo4[i];

                       }

      break;

       case 5:
     for(int i=0;i<ptos;i++)
                  {

                  
                  fonda[i]=fo5[i];

                       }

      break;

       case 6:
   for(int i=0;i<ptos;i++)
                  {

                  
                  fonda[i]=fo6[i];

                       }

      break;

       case 7:
  for(int i=0;i<ptos;i++)
                  {

                  
                  fonda[i]=fo7[i];

                       }

      break;

       case 8:
     for(int i=0;i<ptos;i++)
                  {

                  
                  fonda[i]=fo8[i];

                       }

      break;

       case 9:
    for(int i=0;i<ptos;i++)
                  {

                  
                  fonda[i]=fo9[i];

                       }

      break;

       case 10:
    for(int i=0;i<ptos;i++)
                  {

                  
                  fonda[i]=fo10[i];

                       }

      break;

       case 11:
     for(int i=0;i<ptos;i++)
                  {

                  
                  fonda[i]=fo11[i];

                       }

      break;

       case 12:
    for(int i=0;i<ptos;i++)
                  {

                  
                  fonda[i]=fo12[i];

                       }

      break;

       case 13:
 for(int i=0;i<ptos;i++)
                  {

                  
                  fonda[i]=fo13[i];

                       }

      break;

       case 14:
    for(int i=0;i<ptos;i++)
                  {

                  
                  fonda[i]=fo14[i];

                       }

      break;
        case 15:

                for(int i=0;i<ptos;i++)
                  {
                  
                  fonda[i]=fo15[i];

                    }

      break;
        case 16:

                for(int i=0;i<ptos;i++)
                  {
                  
                  fonda[i]=fo16[i];

                    }

      break;
        case 17:

                for(int i=0;i<ptos;i++)
                  {
                  
                  fonda[i]=fo17[i];

                    }

      break;
        case 18:

                for(int i=0;i<ptos;i++)
                  {
                  
                  fonda[i]=fo18[i];

                    }

      break;
        case 19:

                for(int i=0;i<ptos;i++)
                  {
                  
                  fonda[i]=fo19[i];

                    }

      break;
        case 20:

                for(int i=0;i<ptos;i++)
                  {
                  
                  fonda[i]=fo20[i];

                    }

      break;
        case 21:

                for(int i=0;i<ptos;i++)
                  {
                  
                  fonda[i]=fo21[i];

                    }

      break;
        case 22:

                for(int i=0;i<ptos;i++)
                  {
                  
                  fonda[i]=fo22[i];

                    }

      break;

  }


}

double  selecfonda2(int opcion,double *normalizacion)
{

int ptos=300;
 char *cambio;
cambio=new char [8];
double posicion;
fstream intp1(mi_strcpy(cambio,opcion));

switch(opcion)
  {
 case 15:
                  for(int i=0;i<ptos;i++)
                  {
                 for(int r=1;r<3;r++)
                  {
                   if(r==1)
                     {
                     intp1>>posicion;
                     }else
                     {
                        intp1>>fo15[i];
                        fo15[i]=fo15[i]/normalizacion[15];
                     }
                    }
                       }

      break;
      case 16:
                  for(int i=0;i<ptos;i++)
                  {
                  for(int r=1;r<3;r++)
                  {
                   if(r==1)
                     {
                     intp1>>posicion;
                     }else
                     {
                        intp1>>fo16[i];
                        fo16[i]=fo16[i]/normalizacion[16];
                     }
                    }
                       }

      break;
      case 17:
                  for(int i=0;i<ptos;i++)
                  {
                   for(int r=1;r<3;r++)
                  {
                   if(r==1)
                     {
                     intp1>>posicion;
                     }else
                     {
                        intp1>>fo17[i];
                        fo17[i]=fo17[i]/normalizacion[17];
                     }
                    }
                       }

      break;
      case 18:
                  for(int i=0;i<ptos;i++)
                  {

                   for(int r=1;r<3;r++)
                  {

                   if(r==1)
                     {
                     intp1>>posicion;
                     }else
                     {
                        intp1>>fo18[i];
                        fo18[i]=fo18[i]/normalizacion[18];
                     }

                    }


                       }

      break;
      case 19:
                  for(int i=0;i<ptos;i++)
                  {

                   for(int r=1;r<3;r++)
                  {

                   if(r==1)
                     {
                     intp1>>posicion;
                     }else
                     {
                        intp1>>fo19[i];
                        fo19[i]=fo19[i]/normalizacion[19];
                     }

                    }


                       }

      break;
      case 20:
                  for(int i=0;i<ptos;i++)
                  {

                   for(int r=1;r<3;r++)
                  {

                   if(r==1)
                     {
                     intp1>>posicion;
                     }else
                     {
                        intp1>>fo20[i];
                        fo20[i]=fo20[i]/normalizacion[20];
                     }

                    }


                       }

      break;
      case 21:
                  for(int i=0;i<ptos;i++)
                  {

                   for(int r=1;r<3;r++)
                  {

                   if(r==1)
                     {
                     intp1>>posicion;
                     }else
                     {
                        intp1>>fo21[i];
                        fo21[i]=fo21[i]/normalizacion[21];
                     }

                    }


                       }

      break;
      case 22:
                  for(int i=0;i<ptos;i++)
                  {

                   for(int r=1;r<3;r++)
                  {

                   if(r==1)
                     {
                     intp1>>posicion;
                     }else
                     {
                        intp1>>fo22[i];
                        fo22[i]=fo22[i]/normalizacion[22];
                     }

                    }


                       }

      break;

}
}

////////////////////////////////////////////////////////////////////////////////////////
double selecfonda(int opcion,double *normalizacion)
{
 int ptos=300;
 char *cambio;
cambio=new char [8];
double posicion;
fstream intp1(mi_strcpy(cambio,opcion));

 switch(opcion)
  {
    case 0:
                for(int i=0;i<ptos;i++)
                  {
                  for(int r=1;r<3;r++)
                  {
                     if(r==1)
                     {
                     intp1>>posicion;
                     }else
                     {
                        intp1>>fo0[i];
                        fo0[i]=fo0[i]/normalizacion[0];
                     }
                    }
                       }

      break;
    
    case 1:
                for(int i=0;i<ptos;i++)
                  {
                 for(int r=1;r<3;r++)
                  {
                     if(r==1)
                     {
                     intp1>>posicion;
                     }else
                     {
                        intp1>>fo1[i];
                        fo1[i]=fo1[i]/normalizacion[1];
                     }
                    }
                       }

      break;

      case 2:
     for(int i=0;i<ptos;i++)
                  {
                  for(int r=1;r<3;r++)
                  {
                if(r==1)
                     {
                     intp1>>posicion;
                     }else
                     {
                        intp1>>fo2[i];
                        fo2[i]=fo2[i]/normalizacion[2];
                     }
                   }
                       }

      break;

      case 3:
        for(int i=0;i<ptos;i++)
                  {
                for(int r=1;r<3;r++)
                  {
                   if(r==1)
                     {
                     intp1>>posicion;
                     }else
                     {
                        intp1>>fo3[i];
                        fo3[i]=fo3[i]/normalizacion[3];
                     }
                   }
                       }

      break;
       case 4:
     for(int i=0;i<ptos;i++)
                  {
                for(int r=1;r<3;r++)
                  {
                if(r==1)
                     {
                     intp1>>posicion;
                     }else
                     {
                        intp1>>fo4[i];
                        fo4[i]=fo4[i]/normalizacion[4];
                     }
                    }
                    }

      break;
       case 5:
                    for(int i=0;i<ptos;i++)
                  {

                   for(int r=1;r<3;r++)
                  {
                  if(r==1)
                    {
                     intp1>>posicion;
                     }else
                     {
                        intp1>>fo5[i];
                        fo5[i]=fo5[i]/normalizacion[5];
                     }
                    }
                       }

      break;
       case 6:
               for(int i=0;i<ptos;i++)
                  {
                   for(int r=1;r<3;r++)
                  {
                  if(r==1)
                     {
                     intp1>>posicion;
                     }else
                     {
                        intp1>>fo6[i];
                        fo6[i]=fo6[i]/normalizacion[6];
                     }
                   }
                       }

      break;
       case 7:
                  for(int i=0;i<ptos;i++)
                  {
                for(int r=1;r<3;r++)
                  {
               if(r==1)
                     {
                     intp1>>posicion;
                     }else
                     {
                        intp1>>fo7[i];
                        fo7[i]=fo7[i]/normalizacion[7];
                     }
                  }
                       }

      break;
       case 8:
                  for(int i=0;i<ptos;i++)
                  {
                for(int r=1;r<3;r++)
                  {
               if(r==1)
                     {
                     intp1>>posicion;
                     }else
                     {
                        intp1>>fo8[i];
                        fo8[i]=fo8[i]/normalizacion[8];
                     }
                   }
                      }

      break;
       case 9:
     for(int i=0;i<ptos;i++)
                  {
                for(int r=1;r<3;r++)
                  {
                 if(r==1)
                     {
                     intp1>>posicion;
                     }else
                     {
                        intp1>>fo9[i];
                        fo9[i]=fo9[i]/normalizacion[9];
                     }
                    }
                       }

      break;
       case 10:
                   for(int i=0;i<ptos;i++)
                  {
                for(int r=1;r<3;r++)
                  {
                   if(r==1)
                     {
                     intp1>>posicion;
                     }else
                     {
                        intp1>>fo10[i];
                        fo10[i]=fo10[i]/normalizacion[10];
                     }
                    }
                       }

      break;
       case 11:
                    for(int i=0;i<ptos;i++)
                  {
                   for(int r=1;r<3;r++)
                  {
                  if(r==1)
                     {
                     intp1>>posicion;
                     }else
                     {
                        intp1>>fo11[i];
                        fo11[i]=fo11[i]/normalizacion[11];
                     }
                    }
                       }

      break;
       case 12:
                   for(int i=0;i<ptos;i++)
                  {
                   for(int r=1;r<3;r++)
                  {
                   if(r==1)
                     {
                     intp1>>posicion;
                     }else
                     {
                        intp1>>fo12[i];
                        fo12[i]=fo12[i]/normalizacion[12];
                     }
                    }
                       }

      break;
       case 13:
                for(int i=0;i<ptos;i++)
                  {
                  for(int r=1;r<3;r++)
                  {

                   if(r==1)
                     {
                     intp1>>posicion;
                     }else
                     {
                        intp1>>fo13[i];
                        fo13[i]=fo13[i]/normalizacion[13];
                     }
                    }
                       }

      break;
       case 14:
                  for(int i=0;i<ptos;i++)
                  {
                  for(int r=1;r<3;r++)
                  {
                   if(r==1)
                     {
                     intp1>>posicion;
                     }else
                     {
                        intp1>>fo14[i];
                        fo14[i]=fo14[i]/normalizacion[14];
                     }
                    }
                       }

      break;
     
default :
selecfonda2(opcion,normalizacion);

  }
 delete [] cambio;
}
*/

