double funciones(int opcion,int m,double t,double x[],int estado,double constante,int estadoescalado)
{
    double r;
  
double mBC_1=1/mBC,mABC_1=1/mABC;

    if(opcion==1)
    {
    r=x[3]/mABC;
    
 return r;
    }
    
    if(opcion==2)
    {
     r=x[4]*(constante+mABC_1/(x[1]*x[1]));
     
 return r;
    }
    
if(opcion==3)
{
   double DR,TMP,paso,intf12;
    double h=pasogrilla,h2=2*h,h3=3*h,h4=4*h;//paso de la funcion de onda
DR=5e-6;
     TMP=x[1]+DR;
      DR=TMP-x[1];
 intf12=0;
paso=pasoinicialgrilla;///Las funciones de onda empiezan en 2 AA
double DR_1=1/DR;

  for(int i=0;i<nptosfonda;i+=4) /// esta OK
  {

intf12 = intf12+7*fondatotal[i][estadoescalado]*(0.5*(Morse(paso,TMP,x[2])-Morse(paso,x[1]-DR,x[2]))*DR_1)*fondatotal[i][estadoescalado]+32*fondatotal[i+1][estadoescalado]*(0.5*(Morse(paso+h,TMP,x[2])-Morse(paso+h,x[1]-DR,x[2]))*DR_1)*fondatotal[i+1][estadoescalado]+12*fondatotal[i+2][estadoescalado]*(0.5*(Morse(paso+h2,TMP,x[2])-Morse(paso+h2,x[1]-DR,x[2]))*DR_1)*fondatotal[i+2][estadoescalado]+32*fondatotal[i+3][estadoescalado]*(0.5*(Morse(paso+h3,TMP,x[2])-Morse(paso+h3,x[1]-DR,x[2]))*DR_1)*fondatotal[i+3][estadoescalado]+7*fondatotal[i+4][estadoescalado]*(0.5*(Morse(paso+h4,TMP,x[2])-Morse(paso+h4,x[1]-DR,x[2]))*DR_1)*fondatotal[i+4][estadoescalado];
paso+=h4;
}
//intf12=intf12*h/22.5;
r=x[4]*x[4]*mABC_1/(x[1]*x[1]*x[1])-(intf12*h/22.5);
 return r;
}

if(opcion==4)  ///Derivada con respecto a theta
{
  double DTH,TMP,paso,intf12;
    double h=pasogrilla,h2=2*h,h3=3*h,h4=4*h;//paso de la funcion de onda
DTH=5e-6;
  TMP=x[2]+DTH;
      DTH=TMP-x[2];
  intf12=0;
paso=pasoinicialgrilla;///Las funciones de onda empiezan en 2 AA
double DTH_1=1/DTH;
  for(int i=0;i<nptosfonda;i+=4) /// esta OK
  {

intf12 = intf12+7*fondatotal[i][estadoescalado]*(0.5*(Morse(paso,x[1],TMP)-Morse(paso,x[1],x[2]-DTH))*DTH_1)*fondatotal[i][estadoescalado]+32*fondatotal[i+1][estadoescalado]*(0.5*(Morse(paso+h,x[1],TMP)-Morse(paso+h,x[1],x[2]-DTH))*DTH_1)*fondatotal[i+1][estadoescalado]+12*fondatotal[i+2][estadoescalado]*(0.5*(Morse(paso+h2,x[1],TMP)-Morse(paso+h2,x[1],x[2]-DTH))*DTH_1)*fondatotal[i+2][estadoescalado]+32*fondatotal[i+3][estadoescalado]*(0.5*(Morse(paso+h3,x[1],TMP)-Morse(paso+h3,x[1],x[2]-DTH))*DTH_1)*fondatotal[i+3][estadoescalado]+7*fondatotal[i+4][estadoescalado]*(0.5*(Morse(paso+h4,x[1],TMP)-Morse(paso+h4,x[1],x[2]-DTH))*DTH_1)*fondatotal[i+4][estadoescalado];
paso+=h4;
}
r=-1*intf12*h/22.5;
 return r;
}

}

double Aleatorio (double maximo, double minimo)
{
	double aleatorio;
	aleatorio=(double)(maximo-minimo)*(double)rand()/RAND_MAX+minimo;
   return aleatorio;
}

///////////////////////////////////// Para crear y liberar apuntadores///////////////
 void Releaseme (double **v,int m)
{
   for(int t=0;t<m;t++)
      delete [] v[t];

   delete [] v;
}


 double **Makematriz(int m)
{
double **matriz;
matriz = new double* [m];

for(int j=0; j<m; j++)
    {
        matriz[j] = new double [m];
    }


return (matriz);
}

 
/////////////////////////////////////////////////////////////////////////////


int saltodow(double eta,double *probsalto,double &probdownfinal,int estadoinicial ) ///De nv hacia nv-1,-2,....
{
  if(probsalto[estadoinicial-1]>eta && (estadoinicial-1)>-1 ) /// Solo existen los estados de 1-m curvas diabaticas
  {
         
//  cout<<" salto de "<<estadoinicial<<" - "<<estadoinicial-1<<" R="<<posicion<<"  Probabilidad de salto "<<probsalto[estadoinicial-1]*100<<"%"<<endl;
  probdownfinal=probsalto[estadoinicial-1];
       return (estadoinicial-1);
  }else
  {
    for(int j=estadoinicial-2;j>=0;j--) // Si no salta para el nivel inmediato, para saber hacia que nivel salta, se incluye el 0 porque 0->nv=15
    {
      double suma=0,sumaanterior=0;

        if(estadoinicial!=j)
    {
     for(int r=estadoinicial-1;r>=j;r--)
    {
     suma+=probsalto[r];                         ///j=2, estadoinicial=4; suma=g43+g42
    }
      for(int r=estadoinicial-1;r>=j+1;r--)
    {
     sumaanterior+=probsalto[r];                   ///j=2, estadoinicial=4; suma=g43
    }
      if(sumaanterior<eta && eta<suma)
      {
          //  cout<<posicion<<" "<<v[j-1]<<" "<<v[j]<<endl;
     // cout<<" salto de "<<estadoinicial<<" - "<<j<<" R="<<posicion<<" Probabilidad de salto"<<probsalto[j]*100<<"%"<<endl;
 probdownfinal=suma;///la probabilidad donde vaya a saltar
        return j;
      }     
    }
  }
  return estadoinicial;
 }

 }

 int saltoup(double eta,double *probsalto,double &probupfinal,int estadoinicial,int m ) ///De nv hacia nv+1,+2,....
{
  if(probsalto[estadoinicial+1]>eta && (estadoinicial+1)<=m)
  {
 //cout<<" salto de "<<estadoinicial<<" - "<<estadoinicial+1<<" R="<<posicion<<" Probabilidad de salto "<<probsalto[estadoinicial+1]*100<<"%"<<endl;
  probupfinal=probsalto[estadoinicial+1];
       return (estadoinicial+1);

  }
  else
  {
    for(int j=estadoinicial+2;j<m;j++) // Si no salta para el nivel inmediato, para saber hacia que nivel salta
    {
      double sumapos=0,sumanterior=0;

        if(estadoinicial!=j)
    {
     for(int r=estadoinicial+1;r<=j;r++)  /// es la suma g12+g13+g14...... por ejemplo
    {
     sumapos+=probsalto[r];                   ///g12+g13<eta<g12+g13+g14
    }
         for(int r=j-1;r>estadoinicial;r--)  /// es la suma g12+g13 por ejemplo
    {
     sumanterior+=probsalto[r];
    }

      if(sumanterior<eta && eta<sumapos)
      {
     // cout<<" salto de "<<estadoinicial<<" - "<<j<<" R="<<posicion<<"  Probabilidad de salto "<<probsalto[j]*100<<"%"<<" suma="<<sumapos<<endl;
  probupfinal=sumapos;///la probabilidad donde vaya a saltar
        return j;
      }
    }

  }
  return estadoinicial;
 }


}

double Salto(double *probsalto,int m,double posicion,int  estadoinicial)
{
double eta,probdown=0,probup=0;
  eta=(double)(1-0)*(double)rand()/RAND_MAX+0;
int estadown=estadoinicial,estaup=estadoinicial; ///Estados iniciales, si existe la probabilidad de salto cambiara, sino se mantiene inicialmente

   if(probsalto[estadoinicial+1]>eta && (estadoinicial+1)<=m)
  {
 cout<<" salto de "<<estadoinicial<<" - "<<estadoinicial+1<<" R="<<posicion<<" Probabilidad de salto "<<probsalto[estadoinicial+1]*100<<"%"<<endl;
 cout<<"probup="<<probsalto[estadoinicial+1]<<endl;
       probup=probsalto[estadoinicial+1];
       estaup=(estadoinicial+1);
  }
  else
  {
    for(int j=estadoinicial+2;j<m;j++) // Si no salta para el nivel inmediato, para saber hacia que nivel salta
    {
      double sumapos=0,sumanterior=0;

        if(estadoinicial!=j)
    {
     for(int r=estadoinicial+1;r<=j;r++)  /// es la suma g12+g13+g14...... por ejemplo
    {
     sumapos+=probsalto[r];                   ///g12+g13<eta<g12+g13+g14
    }
         for(int r=j-1;r>estadoinicial;r--)  /// es la suma g12+g13 por ejemplo
    {
     sumanterior+=probsalto[r];
    }

      if(sumanterior<eta && eta<sumapos)
      {
      cout<<" salto de "<<estadoinicial<<" - "<<j<<" R="<<posicion<<"  Probabilidad de salto "<<probsalto[j]*100<<"%"<<" sumapos="<<sumapos<<" sumaanterior="<<sumanterior<<" eta="<<eta<<endl;
  cout<<"probup "<<probsalto[j]<<endl;
       probup=probsalto[j];
       estaup=j;
      }
    }

  }   
 }

 if(probsalto[estadoinicial-1]>eta && (estadoinicial-1)>-1 ) /// Solo existen los estados de 1-m curvas diabaticas
  {
         
  cout<<" salto de "<<estadoinicial<<" - "<<estadoinicial-1<<" R="<<posicion<<"  Probabilidad de salto "<<probsalto[estadoinicial-1]*100<<"%"<<endl;
  cout<<"prodown="<<probsalto[estadoinicial-1]<<endl;
        estadown=(estadoinicial-1);
        probdown=probsalto[estadoinicial-1];
  }else
  {
    for(int j=estadoinicial-2;j>=0;j--) // Si no salta para el nivel inmediato, para saber hacia que nivel salta, se incluye el 0 porque 0->nv=15
    {
      double suma=0,sumaanterior=0;

        if(estadoinicial!=j)
    {
     for(int r=estadoinicial-1;r>=j;r--)
    {
     suma+=probsalto[r];                         ///j=2, estadoinicial=4; suma=g43+g42
    }
      for(int r=estadoinicial-1;r>=j+1;r--)
    {
     sumaanterior+=probsalto[r];                   ///j=2, estadoinicial=4; suma=g43
    }
      if(sumaanterior<eta && eta<suma)
      {
          //  cout<<posicion<<" "<<v[j-1]<<" "<<v[j]<<endl;
      cout<<" salto de "<<estadoinicial<<" - "<<j<<" R="<<posicion<<" Probabilidad de salto"<<probsalto[j]*100<<"%"<<" sumapos="<<suma<<" sumaanterior="<<sumaanterior<<" eta="<<eta<<endl;
        cout<<"probdown "<<probsalto[j]<<endl;
        probdown=probsalto[j];
        estadown=j;

      }     
    }
  }
 }

//cout<<" estadodown: "<<estadown<<"  estadoup: "<<estaup<<" probup: "<<probup<<" probdown: "<<probdown<<endl;
if(probdown>probup)
{
cout<<"Se elige el estado:"<<estadown<<endl;
return estadown;
}

if(probup>probdown)
{
cout<<"Se elige el estado:"<<estaup<<endl;
return estaup;
}

if(probup==0 && probdown==0)
return estadoinicial;

if(probup==probdown)
{
cout<<"Ambos tienen la misma probabilidad, se elige random a cual de los estados saltara"<<endl;
double random=(double)(1-0)*(double)rand()/RAND_MAX+0;
if(random>0.5)
{
cout<<"Se elige el estado: "<<estaup<<endl;
return estaup;
}else
{
cout<<"Se elige el estado: "<<estadown<<endl;
return estadown;
}

}


}

double Rescalamiento2(int estadoinicial,int estfinal,int &estadofinal,double *probasalto,double **v,double **B,double vel,double posicion,double angulo,double &newp,complex<double> **coef2,int factordeescala)
{
int interno,interno2;
  cout<<"densidad inicial "<<real(coef[estadoinicial][estadoinicial])*real(coef[estadoinicial][estadoinicial])+imag(coef[estadoinicial][estadoinicial])*imag(coef[estadoinicial][estadoinicial])<<" densidad a la que salta "<<real(coef[estfinal][estfinal])*real(coef[estfinal][estfinal])+imag(coef[estfinal][estfinal])*imag(coef[estfinal][estfinal])<<endl;
      double RP,DRP,DR,DTH,TH,TMP,paso,intf12,intf12der,aa,bb,rescalado;
    double h2=pasogrilla,h_2=2*h2,h3=3*h2,h4=4*h2;//paso de la funcion de onda
     DR=5e-6;
     TMP=posicion+DR;
      DR=TMP-posicion;
 intf12=0;
  intf12der=0;
paso=pasoinicialgrilla;///Las funciones de onda empiezan en 2 AA
double DR_1=1/DR;
for(int i=0;i<nptosfonda;i+=4) /// esta OK
  {
intf12 = intf12+7*fondatotal[i][estadoinicial+factordeescala]*(0.5*(Morse(paso,TMP,angulo)-Morse(paso,posicion-DR,angulo))*DR_1)*fondatotal[i][estfinal+factordeescala]+32*fondatotal[i+1][estadoinicial+factordeescala]*(0.5*(Morse(paso+h2,TMP,angulo)-Morse(paso+h2,posicion-DR,angulo))*DR_1)*fondatotal[i+1][estfinal+factordeescala]+12*fondatotal[i+2][estadoinicial+factordeescala]*(0.5*(Morse(paso+h_2,TMP,angulo)-Morse(paso+h_2,posicion-DR,angulo))*DR_1)*fondatotal[i+2][estfinal+factordeescala]+32*fondatotal[i+3][estadoinicial+factordeescala]*(0.5*(Morse(paso+h3,TMP,angulo)-Morse(paso+h3,posicion-DR,angulo))*DR_1)*fondatotal[i+3][estfinal+factordeescala]+7*fondatotal[i+4][estadoinicial+factordeescala]*(0.5*(Morse(paso+h4,TMP,angulo)-Morse(paso+h4,posicion-DR,angulo))*DR_1)*fondatotal[i+4][estfinal+factordeescala];

 intf12der = intf12der+7*fondatotal[i][estfinal+factordeescala]*(0.5*(Morse(paso,TMP,angulo)-Morse(paso,posicion-DR,angulo))*DR_1)*fondatotal[i][estfinal+factordeescala]+32*fondatotal[i+1][estfinal+factordeescala]*(0.5*(Morse(paso+h2,TMP,angulo)-Morse(paso+h2,posicion-DR,angulo))*DR_1)*fondatotal[i+1][estfinal+factordeescala]+12*fondatotal[i+2][estfinal+factordeescala]*(0.5*(Morse(paso+h_2,TMP,angulo)-Morse(paso+h_2,posicion-DR,angulo))*DR_1)*fondatotal[i+2][estfinal+factordeescala]+32*fondatotal[i+3][estfinal+factordeescala]*(0.5*(Morse(paso+h3,TMP,angulo)-Morse(paso+h3,posicion-DR,angulo))*DR_1)*fondatotal[i+3][estfinal+factordeescala]+7*fondatotal[i+4][estfinal+factordeescala]*(0.5*(Morse(paso+h4,TMP,angulo)-Morse(paso+h4,posicion-DR,angulo))*DR_1)*fondatotal[i+4][estfinal+factordeescala];

paso+=h4;
}
intf12=intf12*h2/22.5;  ///<v'|derMorse|v> 
intf12der=-intf12der*h2/22.5;  ///F=-<v|derMorse|v> aqui calculo la fuerza
 ////////////////////////////////////////////////////Calculo del nuevo momento//////////////////////////////////////////
bb=vel*intf12;                    /// bbij
aa=intf12*intf12/(2*mABC);       /// aaij
interno=estadoinicial+factordeescala;
interno2=estfinal+factordeescala;
if((bb*bb+4*aa*((v[estadoinicial][estadoinicial]+ev[interno]+B[estadoinicial][estadoinicial])-(v[estfinal][estfinal]+ev[interno2]+B[estfinal][estfinal])))<0) ///Energeticamente imposible, el momento cambia de signo en la dinamica
{

if(intf12der*(vel*mABC)>=0)
{
newp=vel*mABC;
 cout<<"Salto frustado de "<<estadoinicial<<"-"<<estfinal<<" preescripcion +"<<endl;///Truhlar
}else
{
newp=-vel*mABC; ////  hopping_OM2 pag 6/30
 cout<<"Salto frustado de "<<estadoinicial<<"-"<<estfinal<<" preescripcion -"<<endl; ///Truhlar
}
   estfinal=estadoinicial; 
}
else
{
if(bb<0) /// Esto se hace para minimizar el salto,para que la variacion sea lo mas pequeña posible
{
   rescalado=(bb+sqrt(bb*bb+4*aa*((v[estadoinicial][estadoinicial]+ev[interno]+B[estadoinicial][estadoinicial])-(v[estfinal][estfinal]+ev[interno2]+B[estfinal][estfinal]))))/(2*aa);
   newp=vel*mABC-rescalado*intf12;
}else
{
  rescalado=(bb-sqrt(bb*bb+4*aa*((v[estadoinicial][estadoinicial]+ev[interno]+B[estadoinicial][estadoinicial])-(v[estfinal][estfinal]+ev[interno2]+B[estfinal][estfinal]))))/(2*aa);
  newp=vel*mABC-rescalado*intf12;
}
}
   estadofinal=estfinal;      
  if(estadofinal!= estadoinicial)
  cout<<"salto desde "<< estadoinicial<<" - "<< estadofinal<<endl;
  
   
}

void RK4(double **w,double *T,int m,double h,int estadoinicial,double constante,int estadoescala,int decido)
{
double **k1,**k2,**k3,**k4,**combinacion1,**combinacion2,**combinacion3;
 k1=Makematriz(7);
 k2=Makematriz(7);  ///porque son solamente 4 ecuaciones a integrar
 k3=Makematriz(7);
 k4=Makematriz(7);
combinacion1=Makematriz(7);
combinacion2=Makematriz(7);
combinacion3=Makematriz(7);

 for(int i=1; i<=3; i++)
    {

       for(int j=1; j<=m; j++)
        {
            k1[i-1][j]=h*funciones(j,m,T[i-1],w[i-1],estadoinicial,constante,estadoescala);

        }

        for(int j=1; j<=m; j++)
        {
            for(int p=1; p<=m; p++)
            {
                combinacion1[i-1][p]=w[i-1][p]+(k1[i-1][p])/2;
            }

            k2[i-1][j]=h*funciones(j,m,T[i-1]+h/2,combinacion1[i-1],estadoinicial,constante,estadoescala);

        }

        for(int j=1; j<=m; j++) ///Step 7
        {
            for(int p=1; p<=m; p++)
            {
                combinacion2[i-1][p]=w[i-1][p]+k2[i-1][p]/2;
            }

            k3[i-1][j]=h*funciones(j,m,T[i-1]+h/2,combinacion2[i-1],estadoinicial,constante,estadoescala);

        }
        for(int j=1; j<=m; j++) ///Step 8
        {
            for(int p=1; p<=m; p++)
            {
                combinacion3[i-1][p]=w[i-1][p]+k3[i-1][p];
            }

            k4[i-1][j]=h*funciones(j,m,T[i-1]+h,combinacion3[i-1],estadoinicial,constante,estadoescala);

        }

        for(int j=1; j<=m; j++) ///Step 9
        {
            w[i][j]=w[i-1][j]+(k1[i-1][j]+2*k2[i-1][j]+2*k3[i-1][j]+k4[i-1][j])/6;

        }


        T[i]=T[i-1]+h;
    }

if(decido==1)
{
    cout<<w[1][1]<<" "<<w[1][2]<<"  "<<w[1][3]<<"  "<<w[1][4]<<endl;
    cout<<w[2][1]<<" "<<w[2][2]<<"  "<<w[2][3]<<"  "<<w[2][4]<<endl;
    cout<<w[3][1]<<" "<<w[3][2]<<"  "<<w[3][3]<<"  "<<w[3][4]<<endl;

int nv0=estadoinicial; ///LO pongo para probar

double *potencial,*valorespaso;
potencial=new double[nptosfonda+1];
int estadointermedio;
nivelesesenciales=0;

if(nivelesesenciales==0)
estadointermedio=0;
else
estadointermedio=nivelesesenciales/2+1;

int factordeescala;

 if(estadointermedio==0)
  {
  factordeescala=0;
  nivelesesenciales=numerofondas-1; ///Se incluye el 0 por esa razon elimino 1
  }else
  {
  factordeescala=abs(nv0-numfondainicial-estadointermedio);
  }
  
 
   cout<<" estado inicial :"<<nv0<<" internamente dentro del programa :"<<estadoinicial<<" factor de escala="<<factordeescala<<endl;
int ncurvas=nivelesesenciales;
 double RP,DRP,DR,DTH,TH,u_dwdR,TMP,intf12,intf12der,paso, h2=pasogrilla,mABC_1=1/mABC,mBC_1=1/mBC,factor=h2/22.5,energia; ////h2 es el paso de la
 double **v,**B;
 v=Makematriz(ncurvas+1);
B=Makematriz(ncurvas+1);
 double h_2=2*h2,h_3=3*h2,h_4=4*h2;
 int inicial=1;

 for(int cc=1;cc<=3;cc++)
 {
  paso=pasoinicialgrilla;  ///Las funciones de onda empiezan en 2 AA
for(int ciclo=0;ciclo<nptosfonda;ciclo++)
 {
   potencial[ciclo]=Morse(paso,w[cc][1],w[cc][2]); /// Me calcula una vez el potencial de Morse en esas coordenadas
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
energiainicial=v[nvibracional][nvibracional]+ev[internamente]+w[cc][3]*w[cc][3]*mABC_1*0.5+w[cc][4]*w[cc][4]*0.5*(B[nvibracional][nvibracional]+mABC_1/(w[cc][1]*w[cc][1]));
cout<<" RKR4 energiainicial"<<cc<<"="<<energiainicial<<" nv"<<nvibracional<<" inter="<<internamente<<" factor de escala="<<factordeescala<<" "<<w[cc][3]<<" "<<w[cc][4]<<endl;
}



}
}

///Metodo FSTU forward (puede ser tambien backward, es decir se retrocede en el tiempo y se disminuye el paso)
double coeficientes(double vel,double posicion,double angulo,int estadoinicial,int &estadofinal,int m,double **v,double **B,double &newp,double tiempo,double energiagap,double tiempofrustrado,double factordeescala)  ///tiempofrustrado-> es el tiempo cuando sucede el salto frustrado y se activa el metodo FSTU
{
complex<double> **K,**combinacion1,**combinacion2,**combinacion3,**coef2;
K=Makematriz2(5,m+1);
coef2=Makematriz2(m+1,m+1);
combinacion1=Makematriz2(m+1,m+1);combinacion2=Makematriz2(m+1,m+1);combinacion3=Makematriz2(m+1,m+1);
complex<double> sumacompleja,cvder,expo,variable;
double *probasalto,aa,bb,rescalado,enelestado,hbar_1=1/hbar,t,h,a=tiempo;///tiempo inicial, tiempo final b, h=(b-a)/N
double sumapro,suma;
probasalto=new double [m+1];
int estfinal,N=1;      
h=boundbroken/(cicloiteraciones*10);      //// Para darle un pequeño intervalo de propagacion de un orden menor
    ///Step 1
complex<double> numeroi=complex<double> (0,-1);/// -i
      t=a;
      
  for(int j=m;j>=1;j--) ///actualizo la matriz
      {
        for(int l=m;l>=1;l--)  
        {         
           coef2[j][l]=coef[l][l];     
                  
            }
              }
              
              
 int interno1,interno2;             
 for(int i=1; i<=N; i++) ///Step 4
    {
     for(int diabata=m;diabata>=1;diabata--) /// Lo que hace es cambiar de fila, por el juego de C.iniciales
       {      
            interno1=diabata+factordeescala; 
            sumacompleja=complex<double> (0,0);    
          for(int ii=m;ii>=1;ii--) ////La suma se realiza por toda la fila
          {
          interno2=ii+factordeescala;
          real(expo)=cos((ev[interno1]-ev[interno2])*hbar_1*t);
          imag(expo)=sin((ev[interno1]-ev[interno2])*hbar_1*t);
          
           real(variable)=real(coef2[diabata][ii])*real(expo)-imag(expo)*imag(coef2[diabata][ii]);
       imag(variable)=imag(coef2[diabata][ii])*real(expo)+real(coef2[diabata][ii])*imag(expo);
            sumacompleja+=variable*(v[diabata][ii]+B[diabata][ii]);
           }   
       real(cvder)=real(sumacompleja)*real(numeroi)-imag(sumacompleja)*imag(numeroi);
       imag(cvder)=imag(sumacompleja)*real(numeroi)+real(sumacompleja)*imag(numeroi);
       K[1][diabata]=cvder*h*hbar_1;   
   }
     for(int j=m;j>=1;j--) ///actualizo la matriz
      {
        for(int l=m;l>=1;l--)  
        {         
           combinacion1[j][l]=coef2[l][l]+K[1][l]*0.5;     //lo unico que varia es el elemento diagonal las otras cosas se mantienen iguales                   
            }
              }
           
        for(int diabata=m;diabata>=1;diabata--) /// Lo que hace es cambiar de fila, por el juego de C.iniciales
       {
       interno1=diabata+factordeescala; 
       sumacompleja=complex<double> (0,0);
         for(int ii=m;ii>=1;ii--) ////La suma se realiza por toda la fila
          {
          interno2=ii+factordeescala;
          real(expo)=cos((ev[interno1]-ev[interno2])*hbar_1*(t+h*0.5));
          imag(expo)=sin((ev[interno1]-ev[interno2])*hbar_1*(t+h*0.5));
          
           real(variable)=real(combinacion1[diabata][ii])*real(expo)-imag(expo)*imag(combinacion1[diabata][ii]);
       imag(variable)=imag(combinacion1[diabata][ii])*real(expo)+real(combinacion1[diabata][ii])*imag(expo);
            sumacompleja+=variable*(v[diabata][ii]+B[diabata][ii]);
           }   
       real(cvder)=real(sumacompleja)*real(numeroi)-imag(sumacompleja)*imag(numeroi);
       imag(cvder)=imag(sumacompleja)*real(numeroi)+real(sumacompleja)*imag(numeroi);
       K[2][diabata]=cvder*h*hbar_1;     
       }  
           
for(int j=m;j>=1;j--)///actualizo la matriz
      {
        for(int l=m;l>=1;l--)  
        {         
           combinacion2[j][l]=coef2[l][l]+K[2][l]*0.5;;     
          }      
              }
           
       for(int diabata=m;diabata>=1;diabata--) /// Lo que hace es cambiar de fila, por el juego de C.iniciales
       {   
       interno1=diabata+factordeescala;
       sumacompleja=complex<double> (0,0);
        for(int ii=m;ii>=1;ii--) ////La suma se realiza por toda la fila
          {
           interno2=ii+factordeescala;
          real(expo)=cos((ev[interno1]-ev[interno2])*hbar_1*(t+h*0.5));
          imag(expo)=sin((ev[interno1]-ev[interno2])*hbar_1*(t+h*0.5));
          
           real(variable)=real(combinacion2[diabata][ii])*real(expo)-imag(expo)*imag(combinacion2[diabata][ii]);
       imag(variable)=imag(combinacion2[diabata][ii])*real(expo)+real(combinacion2[diabata][ii])*imag(expo);
            sumacompleja+=variable*(v[diabata][ii]+B[diabata][ii]);
           }   
       real(cvder)=real(sumacompleja)*real(numeroi)-imag(sumacompleja)*imag(numeroi);
       imag(cvder)=imag(sumacompleja)*real(numeroi)+real(sumacompleja)*imag(numeroi);
       K[3][diabata]=cvder*h*hbar_1;          
       }
 
        for(int j=m;j>=1;j--) ///actualizo la matriz
      {
        for(int l=m;l>=1;l--)  
        {         
           combinacion3[j][l]=coef2[l][l]+K[3][l];     
                  
            }
              }
       
        for(int diabata=m;diabata>=1;diabata--) /// Lo que hace es cambiar de fila, por el juego de C.iniciales
       {  
       interno1=diabata+factordeescala;
       sumacompleja=complex<double> (0,0);
         for(int ii=m;ii>=1;ii--) ////La suma se realiza por toda la fila
          {
          interno2=ii+factordeescala;
          real(expo)=cos((ev[interno1]-ev[interno2])*hbar_1*t);
          imag(expo)=sin((ev[interno1]-ev[interno2])*hbar_1*t);
          
           real(variable)=real(combinacion3[diabata][ii])*real(expo)-imag(expo)*imag(combinacion3[diabata][ii]);
       imag(variable)=imag(combinacion3[diabata][ii])*real(expo)+real(combinacion3[diabata][ii])*imag(expo);
            sumacompleja+=variable*(v[diabata][ii]+B[diabata][ii]);
           }   
       real(cvder)=real(sumacompleja)*real(numeroi)-imag(sumacompleja)*imag(numeroi);
       imag(cvder)=imag(sumacompleja)*real(numeroi)+real(sumacompleja)*imag(numeroi);
       K[4][diabata]=cvder*h*hbar_1;  
    }   
      
     for(int diabata=m;diabata>=1;diabata--) /// Lo que hace es cambiar de fila, por el juego de C.iniciales 
         coef2[diabata][diabata]=coef2[diabata][diabata]+(K[1][diabata]+2.0*K[2][diabata]+2.0*K[3][diabata]+K[4][diabata])/6.0;

///Pongo 2.0 para ue no me de problemas con el compilador   
   
          sumapro=0;
         for(int pp=m;pp>=1;pp--) /// La suma es la densidad que es la traza igual a 1
       sumapro+=real(coef2[pp][pp])*real(coef2[pp][pp])+imag(coef2[pp][pp])*imag(coef2[pp][pp]); ///(a+bi)(a-bi)=a^2+b^2 ->cv_x_cv*
                              
         
    if(abs(sumapro-1)>1e-8) ///Para la conservacion de la norma
  {
  cout<<"sumapr="<<sumapro<<endl;  
    return 0;
  }
 /////////////////////////////////////////////////Probabilidad de salto hacia otro estado
enelestado=2*h*hbar_1/(real(coef2[estadoinicial][estadoinicial])*real(coef2[estadoinicial][estadoinicial])+imag(coef2[estadoinicial][estadoinicial])*imag(coef2[estadoinicial][estadoinicial]));  ////2(deltat)/(hbar*cv*cv)

double termA,termB;
interno1=estadoinicial+factordeescala;
for(int allestados=m;allestados>=1;allestados--)
{
interno2=allestados+factordeescala;
if(allestados!=estadoinicial)
{
///Imag(cv*_x_cv'_x_exp(i(Ev-Ev')t/hbar))=cos(t(Ev-Ev')/hbar)*[real(cv')*imag(cv*)+real(cv*)*imag(cv')]+sen(t(Ev-Ev')/hbar)*[real(cv')*real(cv*)-imag(cv')*imag(cv*)]
 
 termA=cos((ev[interno1]-ev[interno2])*hbar_1*t)*(real(coef[allestados][allestados])*imag(conj(coef[estadoinicial][estadoinicial]))+imag(coef[allestados][allestados])*real(coef[estadoinicial][estadoinicial]));

 termB=sin((ev[interno1]-ev[interno2])*hbar_1*t)*(real(coef[allestados][allestados])*real(coef[estadoinicial][estadoinicial])-imag(coef[allestados][allestados])*imag(conj(coef[estadoinicial][estadoinicial])));
 
//cout<<termA<<" "<<termB<<"  "<<(termA+termB)*(v[estadoinicial][allestados])*enelestado<<endl;
 probasalto[allestados]=(termA+termB)*(v[estadoinicial][allestados])*enelestado;
 
 if(probasalto[allestados]<0  || probasalto[allestados]<1e-15)
 probasalto[allestados]=0;
 
}else
{
probasalto[allestados]=0;
}
}  
///////////////////////////////////////////////////////////////////
  estfinal=Salto(probasalto,m,posicion,estadoinicial);
       if(estadoinicial!=estfinal)   //// Porque puede haber salto frustrado 
    {   
    cout<<"Antes del rescalamiento FSTU"<<estadoinicial<<" estadofinal="<<estadofinal<<endl;
    Rescalamiento2(estadoinicial,estfinal,estadofinal,probasalto,v,B,vel,posicion,angulo,newp,coef2,factordeescala);
    if(estadoinicial!=estadofinal && abs(tiempo-tiempofrustrado)<hbar/(2*abs(energiagap)))
    {
    cout<<"FSTU funciono"<<endl;
    cout<<"diferencia tiempo="<<abs(tiempo-tiempofrustrado)<<" energiagap="<<energiagap<<endl;
    estadofinal=estadofinal;
    estadoinicial=estadoinicial;
    }else
    {
    cout<<"FSTU no funciono"<<endl;
    estadofinal=estadoinicial;
    }
 cout<<" FSTU->estadoinicial="<<estadoinicial<<" estadofinal="<<estadofinal<<endl;
        return 0;
}else
{
//cout<<coef[estadoinicial][estadoinicial]<<" "<<coef[estfinal][estfinal]<<endl;
//cout<<t<<" "<<real(coef2[0][0])<<" "<<real(coef2[1][1])<<" "<<real(coef2[2][2])<<" "<<real(coef2[3][3])<<" "<<real(coef2[4][4])<<" "<<real(coef2[5][5])<<" "<<real(coef2[6][6])<<" "<<real(coef2[7][7])<<" "<<real(coef2[8][8])<<" "<<real(coef2[9][9])<<" "<<real(coef2[10][10])<<" "<<real(coef2[11][11])<<" "<<real(coef2[12][12])<<" "<<real(coef2[13][13])<<" "<<real(coef2[14][14])<<" "<<real(coef2[15][15])<<endl;
// cout<<"No hay probabilidad de  salto en el FSTU"<<endl;
}
/////////////////////////////////////////////////////////////////////////       
  t=a+i*(h);
 }
 

    newp=0;          
estadofinal=estadoinicial;

 for(int t=0;t<5;t++)
      delete [] K[t];
 
delete [] K;
delete [] probasalto;
  Releaseme2 (combinacion1,m+1);
    Releaseme2 (combinacion2,m+1);
      Releaseme2 (combinacion3,m+1);
}


double funcion2(double a,double b,int m, double *alpha,long double N,int  estadoinicial,int &estfinal,double energiagap,double &newp,int factordeescala)
{
int ncurvas=nivelesesenciales;  ///son estadoinicial curvas mas una curva 0 que es nv=14
double energiainicial;
//////////////////////////////////////////////// Potenciales diabaticos
double **v,**B;
int nvibracional,nvfinal;
double*potencial,*valorespaso;
potencial=new double[nptosfonda+1];
valorespaso=new double[nptosfonda+1];
///////////////////////////////Rescalado
int estadofinal,violacion,salto=0,conserva=0;
double vel;
///////////////////////////////////////////////////////////////////////Dinamica
 double delta=0,h,t,**w,**w0,*T;
  int contador=0;
  w=Makematriz(m+1);
 w0=Makematriz(m+1);
T=new double [4];
    h=(b-a)/N; ///Step 1
    T[0]=a;
double division=h/24,tiempodiso=3*h;//El tiempo empieza en 3h porque los primeros pasos pertenecen a RK4, esto es debido al metodo de integracion
    for(int j=1; j<=m; j++)
    {
        w[0][j]=alpha[j];  ///Step 2, w[][], el primer termino significa  punto 1,2,3, el segundo: la ecuacion a integrar,
    }
      nvibracional=estadoinicial;
    estadofinal=nvibracional;  //son iguales hasta que cambie
    
 double RP,DRP,DR,DTH,TH,u_dwdR,TMP,intf12,intf12der,paso, h2=pasogrilla,mABC_1=1/mABC,mBC_1=1/mBC,factor=h2/22.5,energia,constante; ////h2 es el paso de la
 double h_2=2*h2,h_3=3*h2,h_4=4*h2;
 int inicial=1;

 ///////////////////////////////Reservar memoria///////////////////////////////////// Si no hay espacio (me retorna 0) and its over
if((v=Makematriz(ncurvas+1))==NULL)
return 0;
if((B=Makematriz(ncurvas+1))==NULL)
return 0;
/////////////////////////////////////////////////////////////////////////////////
 paso=pasoinicialgrilla;
         intf12der=0;
    for(int i=0;i<nptosfonda;i+=4) /// esta OK
  {
  intf12der=intf12der+7*fondatotal[i][nvibracional+factordeescala]*1/(paso*paso)*fondatotal[i][nvibracional+factordeescala]+32*fondatotal[i+1][nvibracional+factordeescala]*1/((paso+h2)*(paso+h2))*fondatotal[i+1][nvibracional+factordeescala]+12*fondatotal[i+2][nvibracional+factordeescala]*1/((paso+h_2)*(paso+h_2))*fondatotal[i+2][nvibracional+factordeescala]+32*fondatotal[i+3][nvibracional+factordeescala]*1/((paso+h_3)*(paso+h_3))*fondatotal[i+3][nvibracional+factordeescala]+7*fondatotal[i+4][nvibracional+factordeescala]*1/((paso+h_4)*(paso+h_4))*fondatotal[i+4][nvibracional+factordeescala];
paso+=h_4;
}
constante=(intf12der*factor)*mBC_1; /// <v|1/r^2|v>/mBc  
///////////////////////////////////////////////////////////////////////////////////
  RK4(w,T,m,h,nvibracional,constante,nvibracional+factordeescala,0);
   
     for (int i=4; i<=N; i++)       /* STEP 7 */  /* T0, W0 will be used in place of t, w resp. */
      {
         
       for(int j=1;j<=m;j++)
      {
         w0[1][j] = w[3][j]+(55.0*funciones(j,m,T[3],w[3],nvibracional,constante,nvibracional+factordeescala)-59.0*funciones(j,m,T[2],w[2],nvibracional,constante,nvibracional+factordeescala)+37.0*funciones(j,m,T[1],w[1],nvibracional,constante,nvibracional+factordeescala)-9.0*funciones(j,m,T[0],w[0],nvibracional,constante,nvibracional+factordeescala))*division; // predict W(I) //
      }
      for(int j=1;j<=m;j++) ///Las ecuaciones de integracion no depende explicitamente del tiempo
      {
       w0[1][j] = w[3][j]+(9.0*funciones(j,m,delta,w0[1],nvibracional,constante,nvibracional+factordeescala)+19.0*funciones(j,m,T[3],w[3],nvibracional,constante,nvibracional+factordeescala)-5.0*funciones(j,m,T[2],w[2],nvibracional,constante,nvibracional+factordeescala)+funciones(j,m,T[1],w[1],nvibracional,constante,nvibracional+factordeescala))*division; // correct W(I) //
      }
       
        vel=w0[1][3]*mABC_1;
     delta =a+ i*h;  ///Paso del tiempo del metodo Predictor Corrector
    tiempodiso=tiempodiso+h; //Tiempo real del sistema porque cuando salta se reinicia delta(definido arriba)
         ///////////////////////////////////Potenciales diabaticos/////////////////////////
paso=pasoinicialgrilla;  ///Las funciones de onda empiezan en 2 AA
for(int ciclo=0;ciclo<nptosfonda;ciclo++)
 {
   potencial[ciclo]=Morse(paso,w0[1][1],w0[1][2]); /// Me calcula una vez el potencial de Morse en esas coordenadas
   paso+=h2;
 }
for(int ii=ncurvas;ii>=1;ii--) /// Potenciales promedios de van der Waals
{   
  for(int jj=ncurvas;jj>=1;jj--)
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
cout<<a+tiempodiso<<" "<<w0[1][1]<<" "<<w0[1][2]<<" "<<w0[1][3]<<" "<<w0[1][4]<<endl;
coeficientes(vel,w0[1][1],w0[1][2],nvibracional,estadofinal,ncurvas,v,B,newp,tiempodiso+a,energiagap,a,factordeescala);///El FSTU, tiempo del salto frustrado-t del paso reducido, en este caso seria-> abs((tiempodiso+a)-a)
if(nvibracional!=estadofinal)
{
cout<<"hubo salto"<<endl;
estfinal=estadofinal;
return 0;
}
//cout<<tiempodiso<<" "<<w0[1][1]<<" "<<w0[1][3]<<" "<<w0[1][2]<<"  "<<w0[1][4]<<" "<<v[2][2]<<" "<<B[2][2]<<endl;
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
estfinal=nvibracional;

 Releaseme(v,ncurvas+1);
 Releaseme(B,ncurvas+1);
 delete [] potencial,valorespaso;
}




