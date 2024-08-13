#include "LecturayescrituraFOndas.h"
#include "Potenciales.h"
#include "FSTU.h"
#include "Liouville.h"
#define Pi 3.141592654
#define hbarenergy 5.8064843
using namespace std;


//////////////////////////////////////////Para ordenar arrays//////////////////////
void intercambiar(double *x,double *y)
{
 double temp;
temp=*x;
*x=*y;
*y=temp;
return;
}
///////////////////////////////////////////algoritmo de la burbuja
void ordenararrays(double *numeros,int n)
{
for (int i=1;i<=n;i++)
{
 for(int j=i;j<=n;j++)
{
  if (*(numeros+j)<*(numeros+i))
intercambiar(numeros+i,numeros+j); 
}
}
return;
}
/////////////////////////////////////////////////////////////////////////


double queenergia(double w0[],int &nvibracional)
{
int nvfinal;
///////////////////////////////////////////Para analizar energeticamente en que nivel vibracional se encuentra///////////////////////
if(abs(MorseBr2(w0[1])+pow(w0[4],2)/(2*mBC))>=abs(ev[0]+ev[1])/2 && abs(MorseBr2(w0[1])+pow(w0[4],2)/(2*mBC))<=abs(ev[0]))
  {
  nvfinal=1;
 cout<<" nivel="<<nvfinal+15<<endl;
  }else
  {
    if(abs(MorseBr2(w0[1])+pow(w0[4],2)/(2*mBC))<=abs(ev[14]+ev[15])/2 && abs(MorseBr2(w0[1])+pow(w0[4],2)/(2*mBC))>=abs(ev[15]))
    {
    nvfinal=15;
cout<<" nivel="<<nvfinal+15<<endl;
    cout<<" Energia fuera de los limites estudiados"<<endl;
       return 0;
    }else
    {
      int serie=14;
  do
    {  
    nvfinal=0; /// Esto para saber en que nivel vibracional corresponde   
   if(abs(MorseBr2(w0[1])+pow(w0[4],2)/(2*mBC))<=abs(ev[serie]+ev[serie-1])/2 && abs(MorseBr2(w0[1])+pow(w0[4],2)/(2*mBC))>=abs(ev[serie]+ev[serie+1])/2)///ev[j] es el limi superior y ev[j+1] es el lim inf
   {
     nvfinal=serie;
cout<<" nivel="<<nvfinal+15<<endl;
   }
    --serie; 
    if(serie<0)
    {
     cout<<" Energia fuera de los limites estudiados"<<endl;
       return 0;
    }
  
    }while(nvfinal==0);
    }
}

nvibracional=nvfinal;
}

double funcion(double a,double b,int m, double *alpha,long double N,int  estadoinicial,int factordeescala)
{
int ncurvas=nivelesesenciales;  ///son estadoinicial curvas mas una curva 0 que es nv=15
double energiainicial,**v,**B,vel,newp,*potencial,*valorespaso;
int nvibracional,nvfinal;
potencial=new double[nptosfonda+1];
valorespaso=new double[nptosfonda+1];
///////////////////////////////Rescalado
int estadofinal,violacion,salto=0,conserva=0, internamente2, internamente;
float tpobi=0,tpobf;
///////////////////////////////////////////////////////////////////////Dinamica
 double delta=0,h,t,**w,**w0,*T,division=h/24,tiempodiso=3*h,RP,DRP,DR,DTH,TH,u_dwdR,TMP,intf12,intf12der,paso, h2=pasogrilla,mABC_1=1/mABC,mBC_1=1/mBC,factor=h2/22.5,energia,constante,h_2=2*h2,h_3=3*h2,h_4=4*h2; ////h2 es el paso de la;//El tiempo empieza en 3h porque los primeros pasos pertenecen a RK4, esto es debido al metodo de integracion
  int contador=0;
  w=Makematriz(m+1);
 w0=Makematriz(m+1);
T=new double [4];
    h=(b-a)/N; ///Step 1
    T[0]=a;
    for(int j=1; j<=m; j++)
    {
        w[0][j]=alpha[j];  ///Step 2, w[][], el primer termino significa  punto 1,2,3, el segundo: la ecuacion a integrar,
    }
    nvibracional=estadoinicial;
 int inicial=1;
///////////////////////////////////////////////////////////////////////////////
    FILE *salidatotal;
    salidatotal=fopen("Parametro_salida","a+");
    FILE *poblaciones;
    poblaciones=fopen("Poblaciones","a+");

  ///////////////////////////////Reservar memoria///////////////////////////////////// Si no hay espacio (me retorna 0) and its over
v=Makematriz(ncurvas+1);
B=Makematriz(ncurvas+1);
  //////////////////////////////////////////////////////////////////////////////////////////////////
    paso=pasoinicialgrilla;
         intf12der=0;
    for(int i=0;i<nptosfonda;i+=4) /// esta OK
  {
  intf12der=intf12der+7*fondatotal[i][nvibracional+factordeescala]*1/(paso*paso)*fondatotal[i][nvibracional+factordeescala]+32*fondatotal[i+1][nvibracional+factordeescala]*1/((paso+h2)*(paso+h2))*fondatotal[i+1][nvibracional+factordeescala]+12*fondatotal[i+2][nvibracional+factordeescala]*1/((paso+h_2)*(paso+h_2))*fondatotal[i+2][nvibracional+factordeescala]+32*fondatotal[i+3][nvibracional+factordeescala]*1/((paso+h_3)*(paso+h_3))*fondatotal[i+3][nvibracional+factordeescala]+7*fondatotal[i+4][nvibracional+factordeescala]*1/((paso+h_4)*(paso+h_4))*fondatotal[i+4][nvibracional+factordeescala];
paso+=h_4;
}
constante=(intf12der*factor)*mBC_1; /// <v|1/r^2|v>/mBc     
///////////////////////////////////////////////////////////////////////////////////
  RK4(w,T,m,h,nvibracional,constante,nvibracional+factordeescala,1);
   
     for (int i=4; i<=N; i++)       
      {
      if(salto==1)
      {
      cout<<"hay salto"<<" energiainicial="<<energiainicial<<" paso="<<i<<" "<<w0[1][1]<<" "<<w0[1][2]<<" "<<w0[1][3]<<" "<<w0[1][4]<<endl;
      cout<<" RK4 debido al salto"<<endl;
      T[0]=a; 
      for(int ss=m;ss>0;ss--)
       w[0][ss]=w0[1][ss];
       
        RK4(w,T,m,h,nvibracional,constante,nvibracional+factordeescala,1);
  tiempodiso=tiempodiso+3*h; ///los 3 paso del runge kutta
  conserva=1; ///Para que no calcule la conservacion de la energia despues de hacer el RK4
  i=4; ///Se reinicia el RK4
 if(factordeescala!=0)
{
cout<<" new factorrescalado"<<endl;
factordeescala=(nvibracional+factordeescala)-int(nivelesesenciales/2+1);///nuevo factor de escala//
nvibracional=nivelesesenciales/2+1;
}
  salto=0; ///Para que no vuelva a hacer el RK4
      }
      
       for(int j=1;j<=m;j++)
      {
         w0[1][j] = w[3][j]+(55.0*funciones(j,m,T[3],w[3],nvibracional,constante,nvibracional+factordeescala)-59.0*funciones(j,m,T[2],w[2],nvibracional,constante,nvibracional+factordeescala)+37.0*funciones(j,m,T[1],w[1],nvibracional,constante,nvibracional+factordeescala)-9.0*funciones(j,m,T[0],w[0],nvibracional,constante,nvibracional+factordeescala))*division; // predict W(I) //
      }
      for(int j=1;j<=m;j++) ///Las ecuaciones de integracion no depende explicitamente del tiempo
      {
       w0[1][j] = w[3][j]+(9.0*funciones(j,m,delta,w0[1],nvibracional,constante,nvibracional+factordeescala)+19.0*funciones(j,m,T[3],w[3],nvibracional,constante,nvibracional+factordeescala)-5.0*funciones(j,m,T[2],w[2],nvibracional,constante,nvibracional+factordeescala)+funciones(j,m,T[1],w[1],nvibracional,constante,nvibracional+factordeescala))*division; // correct W(I) //
      }
        if(tiempodiso>=boundbroken) ///el tiempo esta en ps
        {
            cout<<"molecula no disociada"<<endl;
            return 0;
        }
        vel=w0[1][3]*mABC_1;
     delta =a+ i*h;  ///Paso del tiempo del metodo Predictor Corrector
    tiempodiso=tiempodiso+h; //Tiempo real del sistema porque cuando salta se reinicia delta(definido arriba)
  internamente=nvibracional+factordeescala;
               ///////////////////////////////////Potenciales diabaticos/////////////////////////
paso=pasoinicialgrilla;  ///Las funciones de onda empiezan en 2 AA
for(int ciclo=0;ciclo<nptosfonda;ciclo++)
 {
   potencial[ciclo]=Morse(paso,w0[1][1],w0[1][2]); /// Me calcula una vez el potencial de Morse en esas coordenadas
   paso+=h2;
 }
for(int ii=ncurvas;ii>=0;ii--) /// Potenciales promedios de van der Waals
{   
  for(int jj=ncurvas;jj>=0;jj--)
  {         
intf12=0;
intf12der=0;
paso=pasoinicialgrilla;
  for(int ll=0;ll<nptosfonda;ll+=4) /// esta OK
  {
intf12 = intf12+7*fondatotal[ll][ii+factordeescala]*potencial[ll]*fondatotal[ll][jj+factordeescala]+32*fondatotal[ll+1][ii+factordeescala]*potencial[ll+1]*fondatotal[ll+1][jj+factordeescala]+12*fondatotal[ll+2][ii+factordeescala]*potencial[ll+2]*fondatotal[ll+2][jj+factordeescala]+32*fondatotal[ll+3][ii+factordeescala]*potencial[ll+3]*fondatotal[ll+3][jj+factordeescala]+7*fondatotal[ll+4][ii+factordeescala]*potencial[ll+4]*fondatotal[ll+4][jj+factordeescala];

intf12der=intf12der+7*fondatotal[ll][ii+factordeescala]*1/(paso*paso)*fondatotal[ll][jj+factordeescala]+32*fondatotal[ll+1][ii+factordeescala]*1/((paso+h2)*(paso+h2))*fondatotal[ll+1][jj+factordeescala]+12*fondatotal[ll+2][ii+factordeescala]*1/((paso+h_2)*(paso+h_2))*fondatotal[ll+2][jj+factordeescala]+32*fondatotal[ll+3][ii+factordeescala]*1/((paso+h_3)*(paso+h_3))*fondatotal[ll+3][jj+factordeescala]+7*fondatotal[ll+4][ii+factordeescala]*1/((paso+h_4)*(paso+h_4))*fondatotal[ll+4][jj+factordeescala];

paso+=h_4;
}
 B[ii][jj]=intf12der*factor*mBC_1*w0[1][4]*w0[1][4]*0.5;  
     B[jj][ii]=B[ii][jj];
 v[ii][jj]=intf12*factor;  
     v[jj][ii]=v[ii][jj];
   }
}
integracion(vel,w0[1][1],w0[1][2],nvibracional,estadofinal,ncurvas,v,B,violacion,newp,inicial,tiempodiso,w0[1],factordeescala);
if(violacion==1)
{
cout<<" Violacion de la norma"<<endl;
return 0;
}
inicial=0;
if(nvibracional!=estadofinal )
{
if(newp!=0)
{
cout<<"old pR="<<w0[1][3]<<" theta="<<w0[1][2]<<" constante old="<<constante<<endl;
w0[1][3]=newp; 
nvibracional=estadofinal;   /// Nuevo nivel vibracional 
       internamente=nvibracional+factordeescala;
  paso=pasoinicialgrilla;
      intf12der=0;
     for(int ll=0;ll<nptosfonda;ll+=4) /// esta OK
       {
       
       intf12der=intf12der+7*fondatotal[ll][nvibracional+factordeescala]*1/(paso*paso)*fondatotal[ll][nvibracional+factordeescala]+32*fondatotal[ll+1][nvibracional+factordeescala]*1/((paso+h2)*(paso+h2))*fondatotal[ll+1][nvibracional+factordeescala]+12*fondatotal[ll+2][nvibracional+factordeescala]*1/((paso+h_2)*(paso+h_2))*fondatotal[ll+2][nvibracional+factordeescala]+32*fondatotal[ll+3][nvibracional+factordeescala]*1/((paso+h_3)*(paso+h_3))*fondatotal[ll+3][nvibracional+factordeescala]+7*fondatotal[ll+4][nvibracional+factordeescala]*1/((paso+h_4)*(paso+h_4))*fondatotal[ll+4][nvibracional+factordeescala];
       
       paso+=h_4;
          }
      constante=(intf12der*factor)*mBC_1; /// <v|1/r^2|v>/mBc
      cout<<"new pR="<<w0[1][3]<<" theta="<<w0[1][2]<<" constante new= "<<constante<<endl;
salto=1;
cout<<"hay salto"<<" energiainicial="<<energiainicial<<" paso="<<i<<endl;
}
}
nvibracional=estadofinal;
if(i!=4 && conserva==0)
{
 if(abs(energiainicial-(v[nvibracional][nvibracional]+ev[internamente]+w0[1][3]*w0[1][3]*mABC_1*0.5+w0[1][4]*w0[1][4]*0.5*(constante+mABC_1/(w0[1][1]*w0[1][1]))))>1e-8)
 {
 cout<<"violacion de la energia"<<" paso="<<i<<endl;
  return 0;
 } 
else
{
internamente=nvibracional+factordeescala;
energiainicial=v[nvibracional][nvibracional]+ev[internamente]+w0[1][3]*w0[1][3]*mABC_1*0.5+w0[1][4]*w0[1][4]*0.5*(constante+mABC_1/(w0[1][1]*w0[1][1]));
if(salto==1)
cout<<energiainicial<<endl;
}
}else
{
energiainicial=v[nvibracional][nvibracional]+ev[internamente]+w0[1][3]*w0[1][3]*mABC_1*0.5+w0[1][4]*w0[1][4]*0.5*(constante+mABC_1/(w0[1][1]*w0[1][1]));
conserva=0;
}
 tpobf=tiempodiso;/// paso de tiempo para la poblacion
 if(abs(tpobf-tpobi)==0.5)
    {
    fprintf(poblaciones, " %1.4f %2.5f %i\n", tiempodiso,w0[1][1],internamente); ///Lo imrpime cuando termina el programa
       tpobi=tiempodiso; ///tpobi aquiere el paso temporal en que se imprimen los valores
    }
        if(w0[1][1]>=rfinal) ///esto es solamente para la ecuacion de R, esta en angstrom
        {
               internamente=nvibracional+factordeescala;
               fprintf(poblaciones, " %1.4f %2.5f %i\n", tiempodiso,w0[1][1],internamente); ///Lo imrpime cuando termina el programa
               fprintf(salidatotal," %2.13f %2.13f %2.13f %2.13f %2.13f  %2.13f %2.13f  %2.13f  %2.13f  %2.13f \n ",w0[1][1],w0[1][2],w0[1][3],w0[1][4],w0[1][3]*w0[1][3]*mABC_1,energiainicial,internamente,(sqrt(1+4*w0[1][4]*w0[1][4]/(hbarenergy*hbarenergy))-1)/2,w0[1][4]*w0[1][4]*constante/(2),tiempodiso); /// R,theta,pr,ptheta,KArBr2,Etotal,canal de salida, j,ErotBr2 ptheta^2*<v|1/r^2|v>/(2*m),tiempo_de_vida
               return (w0[1][3]*w0[1][3]*mABC_1*0.5+w0[1][4]*w0[1][4]*0.5*mABC_1/(w0[1][1]*w0[1][1]));
        }
        for (int k=1; k<=3; k++)    /* STEP 9 */
          {
          T[k-1] = T[k];
          for(int j=1;j<=m;j++)
          {
           w[k-1][j] = w[k][j];
          }
         }

      for(int j=1;j<=m;j++)
      {
        w[3][j] = w0[1][j];
      }
         T[3] = delta;    // STEP 10 //delta es el tiempo de integracion del metodo Predictor-Corrector

}
}

int main(int argc, char *argv[])
{
    srand(time(NULL));
    
  
int nv0=atoi(argv[2]);
cout<<nv0<<endl;


    double a,b;
    long double N;
    int m;
   int entrada;
  
cout<<argv[1]<<endl;

     
    double **matriz;
    matriz=new double* [5000];

    //  fstream in("Input.in");    ///Para permitir el flujo de entrada y salida
    //in>>entrada;

entrada=atoi(argv[3]);


    for(int j=0; j<5000; j++)
    {
        matriz[j] = new double [5000];


    }

  /// argv[1] da por ejemplo v=0/esque1-16.int, argv[2] da el nivel vibracional de partida nv0=16.....29, argv[3] da el numero de corrida 
         
    ifstream ci(argv[1]);
   

    FILE *matriz3;
    matriz3=fopen("Condiciones_iniciales.out","w+");

    for(int i=0; i<5000; i++)
    {

        for(int j=1; j<=7; j++)
        {

            ci>> matriz[i][j];

        }
    }
    fprintf(matriz3,"  r          R           theta    pr             pR         pt\n ");
    for(int i=0; i<5000; i++)
    {

        for(int j=1; j<=7; j++)
        {
            fprintf(matriz3," %f ",matriz[i][j]);
        }
        fprintf(matriz3," \n");
    }

    fclose(matriz3);

 if(entrada>=5000)
 {
   /*FILE *op;
           op=fopen("Input.in","w+");
           fprintf(op," %i\n",0);
          fclose(op);*/
          entrada=0;
 }




     double  *alpha; ///condiciones iniciales
      alpha=new double[5];


        alpha[1]=matriz[entrada][2];  // R
        alpha[2]=matriz[entrada][3];  // theta
        alpha[3]=matriz[entrada][5];  // pR
        alpha[4]=matriz[entrada][6];  // ptheta
              
      cout<<alpha[1]<<"  "<<alpha[2]<<" "<<alpha[3]<<" "<<alpha[4]<<endl;

Releaseme(matriz,5000);
///Lo que se hace es un sistema de 7 niveles (1,2,3,4,5,6,7) donde 4 representa el nivel en que se encuentra y los demas son 3 niveles para cada lado.
//10->0,11->1,12->2,13->3,14->4,15->5,16->6,17->7,18->8,19->9,20-10,21-11,22-12....,32->22 Asi esta guardado las Fondas.
/// el sistema empieza en unos de esos numeros de arriba, como quiero llevarlo a la notacion de 7 niveles (4 esta centrado) lo que hago es lo sgte:
////////////factordeescala=abs(estadoinical-4) ->estadoinical-10=nv0



/////////////////////////////////////////////////////////////////////Para comprobar que las F.ondas esten normalizadas
lecturaprimera();
normalizafondas();
  N= cicloiteraciones;
    m=4;///son 4 ecuaciones

    a=0;
    b=boundbroken;  ///Mi tiempo esta entre 0 y 600 ps
//cout<<b<<" "<<N<<endl;

//if(10.8>rfinal)
//cout<<"se cumple"<<rfinal<<endl;

int estadointermedio;

if(nivelesesenciales==0)
estadointermedio=0;
else
estadointermedio=nivelesesenciales/2+1;


//cout<<"aqui,,cmdkvfnb="<<estadointermedio<<" nivelesesenciales="<<nivelesesenciales<<endl;
/*double ppaso=pasoinicialgrilla;
for(int i=0;i<nptosfonda;i++)
{
cout<<ppaso<<" "<<fondatotal[i][0]<<" "<<fondatotal[i][10]<<" "<<fondatotal[i][22]<<endl;
ppaso=ppaso+pasogrilla;
}*/

int factordeescala,estadoinicial ;
  
  if(estadointermedio==0)
  {
  factordeescala=0;
  estadoinicial=abs(nv0-numfondainicial); 
  nivelesesenciales=numerofondas-1; ///Se incluye el 0 por esa razon elimino 1
  }else
  {
  factordeescala=abs(nv0-numfondainicial-estadointermedio);
  estadoinicial=estadointermedio;
  }
   cout<<" estado inicial :"<<nv0<<" internamente dentro del programa :"<<estadoinicial<<" factor de escala="<<factordeescala<<endl;

     cout<<"Numero de corrida "<<entrada<<endl;

if(nivelesesenciales!=0)
coef=Makematriz2(nivelesesenciales+1,nivelesesenciales+1);
else
coef=Makematriz2(numerofondas+1,numerofondas+1);

//////////////////////////////////////////////////////////////////////////////
 double **v,**B;
 double *potencial,*valorespaso;
potencial=new double[nptosfonda+1];
int ncurvas=nivelesesenciales;
 v=Makematriz(ncurvas+1);
 B=Makematriz(ncurvas+1);
 double RP,DRP,DR,DTH,TH,u_dwdR,TMP,intf12,intf12der,paso, h2=pasogrilla,mABC_1=1/mABC,mBC_1=1/mBC,factor=h2/22.5,energia; ////h2 es el paso de la
 double h_2=2*h2,h_3=3*h2,h_4=4*h2;
 paso=pasoinicialgrilla;  ///Las funciones de onda empiezan en 2 AA
for(int ciclo=0;ciclo<nptosfonda;ciclo++)
 {
   potencial[ciclo]=Morse(paso,alpha[1],alpha[2]); /// Me calcula una vez el potencial de Morse en esas coordenadas
   paso+=h2;
 }
for(int ii=ncurvas;ii>=0;ii--) /// Potenciales promedios de van der Waals
{   
  for(int jj=ncurvas;jj>=0;jj--)
  {         
intf12=0;
intf12der=0;
paso=pasoinicialgrilla;
 
  for(int ll=0;ll<nptosfonda;ll+=4) /// esta OK
  {
intf12 = intf12+7*fondatotal[ll][ii+factordeescala]*potencial[ll]*fondatotal[ll][jj+factordeescala]+32*fondatotal[ll+1][ii+factordeescala]*potencial[ll+1]*fondatotal[ll+1][jj+factordeescala]+12*fondatotal[ll+2][ii+factordeescala]*potencial[ll+2]*fondatotal[ll+2][jj+factordeescala]+32*fondatotal[ll+3][ii+factordeescala]*potencial[ll+3]*fondatotal[ll+3][jj+factordeescala]+7*fondatotal[ll+4][ii+factordeescala]*potencial[ll+4]*fondatotal[ll+4][jj+factordeescala];

intf12der=intf12der+7*fondatotal[ll][ii+factordeescala]*1/(paso*paso)*fondatotal[ll][jj+factordeescala]+32*fondatotal[ll+1][ii+factordeescala]*1/((paso+h2)*(paso+h2))*fondatotal[ll+1][jj+factordeescala]+12*fondatotal[ll+2][ii+factordeescala]*1/((paso+h_2)*(paso+h_2))*fondatotal[ll+2][jj+factordeescala]+32*fondatotal[ll+3][ii+factordeescala]*1/((paso+h_3)*(paso+h_3))*fondatotal[ll+3][jj+factordeescala]+7*fondatotal[ll+4][ii+factordeescala]*1/((paso+h_4)*(paso+h_4))*fondatotal[ll+4][jj+factordeescala];

paso+=h_4;
}
 B[ii][jj]=intf12der*factor*mBC_1; ///<v|1/r^2|v>  
     B[jj][ii]=B[ii][jj];
 v[ii][jj]=intf12*factor;  
     v[jj][ii]=v[ii][jj];
   }
}

int nvibracional=estadoinicial;
int internamente=nvibracional+factordeescala;
double energiainicial;
energiainicial=v[nvibracional][nvibracional]+ev[internamente]+alpha[3]*alpha[3]*mABC_1*0.5+alpha[4]*alpha[4]*0.5*(B[nvibracional][nvibracional]+mABC_1/(alpha[1]*alpha[1]));

cout<<" energia total de partida="<<energiainicial<<endl;
/////////////////////////////////////////////////////////////////////////////
 FILE *erest;
    erest=fopen("Erest","a+");  
 double ereste;
ereste=funcion(a,b,m,alpha,N,estadoinicial,factordeescala);
fprintf(erest, "%2.10f \n",ereste);


    FILE *opp;
    opp=fopen("final","w+");
           fprintf(opp," %i\n",entrada+1);  
           fclose(opp);
           

     delete [] alpha;



return 0;
}

