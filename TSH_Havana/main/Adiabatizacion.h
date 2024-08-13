//**********************************Metodo de Jacobi*********************************
void jacobi (int orden, double **a, double *d)
 {
   double **v;
v = new double* [orden+1];

   for(int j=0; j<orden+1; j++)
    {
        v[j] = new double [orden+1];
    }
 
   // Define las variables de tipo entero

   long i, j, ip, iq, nrot;

// Define las variables de tipo doble

   double *b, *z;

   b = new double [orden];

   z = new double [orden];

   b = new double [orden]; z = new double [orden];

   double c, g, h, s, sm, t, tau, theta, tresh;

   // Inicializa a la matriz identidad

   for (ip = 0; ip < orden; ip++) {

      for (iq = 0; iq < orden; iq++) {

         v[ip][iq] = 0;

      }

      v[ip][ip] = 1;

   }

// Inicializa b y d a la diagonal de a

   for (ip = 0; ip < orden; ip++) {

      b[ip] = a[ip][ip];

      d[ip] = b[ip];

      z[ip] = 0;

   }

   nrot = 0;
    for (i = 0; i < 50; i++) {

      sm = 0;

      for (ip = 0; ip < orden - 1; ip++) {

         for (iq = ip + 1; iq < orden; iq++) {

            sm +=fabs(a[ip][iq]);

         }

      }

      if (sm == 0) break;

      if (i < 4)

         tresh = 0.2*sm/(orden*orden);

      else

         tresh = 0.0;

      for (ip =0; ip < orden -1; ip++) {

         for (iq = ip + 1; iq < orden; iq++) {

            g = 100.0*fabs(a[ip][iq]);

            if(i>4 && (double)(fabs(d[ip])+g) == (double)fabs(d[ip])

               && (double)(fabs(d[iq])+g) == (double)fabs(d[iq]))

               a[ip][iq] = 0.0;

            else if (fabs(a[ip][iq]) > tresh) {

               h = d[iq] - d[ip];

               if ((double)(fabs(h)+g) == (double)fabs(h))

                  t = (a[ip][iq])/h;   // t = 1/(2theta)

                   else {

                  theta = 0.5*h/(a[ip][iq]);

                  t = 1.0/(fabs(theta)+sqrt(1.0+theta*theta));

                  if(theta < 0.0) t = -t;

               }

                  c = 1.0/sqrt(1+t*t);

                  s = t*c;

                  tau = s/(1.0+c);

                  h = t*a[ip][iq];

                  z[ip] -=h;

                  z[iq] +=h;

                  d[ip] -=h;

                  d[iq] +=h;

                  a[ip][iq] = 0.0;
 // Varía desde 0 hasta  ip - 1

               for (j =0; j < ip; j++) {

                  g = a[j][ip];

                  h = a[j][iq];

                  a[j][ip] = g - s*(h+g*tau);

                  a[j][iq] = h + s*(g-h*tau);

               }

// Varía desde ip+1 hasta  iq - 1

               for (j =ip+1; j < iq; j++) {

                  g = a[ip][j];

                  h = a[j][iq];

                  a[ip][j] = g - s*(h+g*tau);

                  a[j][iq] = h + s*(g-h*tau);

               }

               for (j =iq+1; j < orden; j++) {

                  g = a[ip][j];

                  h = a[iq][j];

                  a[ip][j] = g - s*(h+g*tau);

                  a[iq][j] = h + s*(g-h*tau);

               }


               for (j =0; j < orden; j++) {

                  g = v[j][ip];

                  h = v[j][iq];

                  v[j][ip] = g - s*(h+g*tau);

                  v[j][iq] = h + s*(g-h*tau);
   }

               ++(nrot);

            }

         }

      }

         for (ip = 0; ip < orden; ip++) {

            b[ip] = b[ip]+z[ip];

            d[ip] = b[ip];

            z[ip] = 0.0;

         }

   }

   delete [] b, z;
 for(int t=0;t<orden+1;t++)
      delete [] v[t];

   delete [] v;
    

  
}



void change(double *d,int i,int curvas)
{
    double cambio1,cambio2;
  int j;


  do{
  j=i+1;


  do
  {
      cambio1=d[i];
  cambio2=d[j];

   if(d[i]<=d[j] )
    {

       d[j]=cambio1;
      d[i]=cambio2;

        }else
        {
          cambio1=d[i];
          cambio2=d[j];
        }

   ++j;

 }while(j<curvas);
++i;
}while(i<(curvas-1));



}

