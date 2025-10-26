#ifndef FONCTION_H
#define FONCTION_H


double* v(double x, double y);
double u0(double x, double y);







double alpha( double x,double y);
double beta( double x,double y);
double gama( double x,double y);
    


int modulo_pos(int k,int d);

double u_exacte(double x, double y ,double t);
double error2(double* sol);



#endif // FONCTION_H