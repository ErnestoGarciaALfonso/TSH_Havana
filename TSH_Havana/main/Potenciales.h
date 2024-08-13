//////////////////////////////////////////////Potenciales tipo Morse
double Morse(double r,double R,double cita) ///  Ar-Br2 (B)
{
    double U_vwd,R_pri_1,R_pri_2;
     double D=114;
    double alpha=1.8;
    double R_medio=4.1;


         R_pri_1=sqrt(pow(R,2)+pow(r,2)/4+r*R*cos(cita));
         R_pri_2=sqrt(pow(R,2)+pow(r,2)/4-r*R*cos(cita));

    U_vwd=D*(exp(-1*2*alpha*(R_pri_1-R_medio))-2*exp(-1*alpha*(R_pri_1-R_medio)))+D*(exp(-1*2*alpha*(R_pri_2-R_medio))-2*exp(-1*alpha*(R_pri_2-R_medio)));
    return(U_vwd);
}


double MorseBr2(double r)/// Br-Br (B)
{
    double Dr=3788;
    double alphar=2.045;
    double r_medio=2.667;
    double   udwrBr2;
    udwrBr2=Dr*(exp(-1*2*alphar*(r-r_medio))-2*exp(-1*alphar*(r-r_medio)));
    return (udwrBr2);
}




