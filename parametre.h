// parametre.h

#ifndef PARAMETRE_H
#define PARAMETRE_H

extern int space_scheme, time_scheme;
extern double dx, dy, xmin, xmax, ymin, ymax, Tf, CFL;
extern double dt_imp;
extern int Nx, Ny, cas;

void initialiser_parametres();



#endif // PARAMETRE_H