

double Rescalamiento(int estadoinicial,int estfinal,int &estadofinal,double *probasalto,double **v,double **B,double vel,double posicion,double angulo,double &newp,double &energiagap,int factordeescala)
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
energiagap=v[estfinal][estfinal]-v[estadoinicial][estadoinicial]-vel*vel*mABC*0.5; ///Estudios_teoricos_complejos_moleculares.pdf pag 42/310
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

void probasaltofun(int estadoinicial,int m,double t,double h,double *probasalto,double **v,double factordeescala)
{
double hbar_1=1/hbar;
double enelestado=2*h*hbar_1/(real(coef[estadoinicial][estadoinicial])*real(coef[estadoinicial][estadoinicial])+imag(coef[estadoinicial][estadoinicial])*imag(coef[estadoinicial][estadoinicial]));  ////2(deltat)/(hbar*cv*cv)

double termA,termB;
int interno1,interno2;
interno1=estadoinicial+factordeescala;

int limiteinf;

if(factordeescala==0)
limiteinf=0;
else
limiteinf=1;

for(int allestados=m;allestados>=limiteinf;allestados--)
{
if(allestados!=estadoinicial)
{
interno2=allestados+factordeescala;
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

//cout<<"estado="<<allestados<<" "<<probasalto[allestados]<<endl;
}  

}

int integracion(double vel,double posicion,double angulo,int estadoinicial,int &estadofinal,int m,double **v,double **B,int &violacion,double &newp,int inicial,double tiempo,double *w0,int factordeescala)
{
complex<double> **K,**combinacion1,**combinacion2,**combinacion3;
K=Makematriz2(5,m+1);
int estadointermedio;
int limiteinf;

if(factordeescala==0)
{
estadointermedio=estadoinicial;
limiteinf=0;
}
else
{
estadointermedio=nivelesesenciales/2+1;
limiteinf=1;
}

combinacion1=Makematriz2(m+1,m+1);combinacion2=Makematriz2(m+1,m+1);combinacion3=Makematriz2(m+1,m+1);
complex<double> sumacompleja,cvder,expo,variable;
double *probasalto,aa,bb,rescalado,enelestado,hbar_1=1/hbar,t,h,a=tiempo;///tiempo inicial, tiempo final b, h=(b-a)/N
double sumapro,suma;
probasalto=new double [m+1];
int estfinal,N=1;      
h=boundbroken/cicloiteraciones;
//cout<<"paso Liouville="<<h<<" boundbroken="<<boundbroken<<" cicloiteraciones="<<cicloiteraciones<<" division="<<boundbroken/cicloiteraciones<<endl;     //// Para darle un pequeño intervalo de propagacion
    ///Step 1
complex<double> numeroi=complex<double> (0,-1);/// -i
   //cout<<v[1][1]<<" "<<v[2][2]<<" "<<v[3][3]<<endl;
    t=a;     
if(inicial==1) ///Esto es solo para la primera iteracion
 {      ///las filas muestran las condiciones iniciales para cada v (1-> ncurvas) por ejemplo: 0 0 .... 1 para v1, para v2 0 0 ......1 , etc
   for(int j=m;j>=limiteinf;j--)
      {
        for(int l=m;l>=limiteinf;l--)  ///Tengo que dejarlo asi
        {
         if(l==estadointermedio)  ///// Me imprime un 1 en la columna que esta, 4 siempre va ser el estado en que se encuentra la curva
          {
           coef[j][l]=complex<double>(1,0);
          }else
          {
          coef[j][l]=complex<double>(0,0);
          coef[l][j]=coef[j][l];
          }          
            }
              }
 }else
 {
 
  for(int j=m;j>=limiteinf;j--)  ///Para que me actualize la matriz de cv
      {
        for(int l=m;l>=limiteinf;l--)  
        {         
           coef[j][l]=coef[l][l];     
                 
            }
              }
 }

 int interno1,interno2;
 for(int i=1; i<=N; i++) ///Step 4
    {
     for(int diabata=m;diabata>=limiteinf;diabata--) /// Lo que hace es cambiar de fila, por el juego de C.iniciales
       {            
            sumacompleja=complex<double> (0,0);   
            interno1=diabata+factordeescala; 
          for(int ii=m;ii>=limiteinf;ii--) ////La suma se realiza por toda la fila
          {
                interno2=ii+factordeescala;
          real(expo)=cos((ev[interno1]-ev[interno2])*hbar_1*t);
          imag(expo)=sin((ev[interno1]-ev[interno2])*hbar_1*t);
          
           real(variable)=real(coef[diabata][ii])*real(expo)-imag(expo)*imag(coef[diabata][ii]);
       imag(variable)=imag(coef[diabata][ii])*real(expo)+real(coef[diabata][ii])*imag(expo);
            sumacompleja+=variable*(v[diabata][ii]+B[diabata][ii]);
           }   
       real(cvder)=real(sumacompleja)*real(numeroi)-imag(sumacompleja)*imag(numeroi);
       imag(cvder)=imag(sumacompleja)*real(numeroi)+real(sumacompleja)*imag(numeroi);
       K[1][diabata]=cvder*h*hbar_1;   
   }
    
     for(int j=m;j>=limiteinf;j--) ///actualizo la matriz
      {
        for(int l=m;l>=limiteinf;l--)  
        {         
           combinacion1[j][l]=coef[l][l]+K[1][l]*0.5;     //lo unico que varia es el elemento diagonal las otras cosas se mantienen iguales                   
            }
              }
           
        for(int diabata=m;diabata>=limiteinf;diabata--) /// Lo que hace es cambiar de fila, por el juego de C.iniciales
       {
       sumacompleja=complex<double> (0,0);
       interno1=diabata+factordeescala;
         for(int ii=m;ii>=limiteinf;ii--) ////La suma se realiza por toda la fila
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
           
for(int j=m;j>=limiteinf;j--)///actualizo la matriz
      {
        for(int l=m;l>=limiteinf;l--)  
        {         
           combinacion2[j][l]=coef[l][l]+K[2][l]*0.5;;     
          }      
              }
           
       for(int diabata=m;diabata>=limiteinf;diabata--) /// Lo que hace es cambiar de fila, por el juego de C.iniciales
       {   
        interno1=diabata+factordeescala;
       sumacompleja=complex<double> (0,0);
        for(int ii=m;ii>=limiteinf;ii--) ////La suma se realiza por toda la fila
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
 
        for(int j=m;j>=limiteinf;j--) ///actualizo la matriz
      {
        for(int l=m;l>=limiteinf;l--)  
        {         
           combinacion3[j][l]=coef[l][l]+K[3][l];     
                  
            }
              }
       
        for(int diabata=m;diabata>=limiteinf;diabata--) /// Lo que hace es cambiar de fila, por el juego de C.iniciales
       {  
       interno1=diabata+factordeescala;
       sumacompleja=complex<double> (0,0);
         for(int ii=m;ii>=limiteinf;ii--) ////La suma se realiza por toda la fila
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
      
     for(int diabata=m;diabata>=limiteinf;diabata--) /// Lo que hace es cambiar de fila, por el juego de C.iniciales 
         coef[diabata][diabata]=coef[diabata][diabata]+(K[1][diabata]+2.0*K[2][diabata]+2.0*K[3][diabata]+K[4][diabata])/6.0;

///Pongo 2.0 para ue no me de problemas con el compilador   
     
          sumapro=0;
         for(int pp=m;pp>=limiteinf;pp--) /// La suma es la densidad que es la traza igual a 1
       sumapro+=real(coef[pp][pp])*real(coef[pp][pp])+imag(coef[pp][pp])*imag(coef[pp][pp]); ///(a+bi)(a-bi)=a^2+b^2 ->cv_x_cv*
                              
     
    if(abs(sumapro-1)>1e-8) ///Para la conservacion de la norma
  {
  cout<<sumapro<<endl;
    violacion=1;
    return 0;
  }
 /////////////////////////////////////////////////Probabilidad de salto hacia otro estado
//cout<<tiempo<<" "<<coef[1][1]<<" "<<coef[2][2]<<" "<<coef[4][4]<<" "<<coef[7][7]<<endl;
probasaltofun(estadoinicial,m,t,h,probasalto,v,factordeescala);
///////////////////////////////////////////////////////////////////
  estfinal=Salto(probasalto,m,posicion,estadoinicial); 
       if(estadoinicial!=estfinal)   //// Porque puede haber salto frustrado 
    {   
    double energiagap;
   // cout<<coef[estadoinicial][estadoinicial]<<" "<<coef[estfinal][estfinal]<<endl;
    Rescalamiento(estadoinicial,estfinal,estadofinal,probasalto,v,B,vel,posicion,angulo,newp,energiagap,factordeescala);
  if(estadoinicial==estadofinal) ///Si hay salto frustrado
  {
  cout<<"se activa FSTU"<<endl;
      //  cout<<"t="<<tiempo<<" "<<w0[1]<<" "<<w0[3]<<" "<<w0[2]<<"  "<<w0[4]<<"  "<<v[2][2]<<" "<<B[2][2]<<endl;
        funcion2(tiempo,tiempo+1e-4,4, w0,10,estadoinicial,estadofinal,energiagap,newp,factordeescala);
        cout<<"se desactiva FSTU"<<endl;
       cout<<"estadoinicial="<<estadoinicial<<" estadofinal="<<estadofinal<<endl;
    }else
    {  
   
  if(factordeescala==0)
      {  
        for(int j=m;j>=limiteinf;j--) ///Reset a los coeficientes
      {
        for(int l=m;l>=limiteinf;l--)  ///Tengo que dejarlo asi
        {
         if(l==estfinal )  ///// Me imprime un 1 en la columna que esta
          {
           coef[j][l]=complex<double>(1,0);
          }else
          {
          coef[j][l]=complex<double>(0,0);
          coef[l][j]=coef[j][l];
          }          
            }
              }  
      }else
      {
        for(int hk=m;hk>=limiteinf;hk--)
        {
          if(hk==estadointermedio)
         coef[hk][hk]=1;
          else
        coef[hk][hk]=0;
             }    
      } 
    }
 
        return 0;
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

