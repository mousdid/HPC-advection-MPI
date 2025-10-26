#include "parametre.h"
#include <stdio.h>
#include </usr/lib/x86_64-linux-gnu/openmpi/include/mpi.h>
 
int space_scheme=1, time_scheme=2;
double dx, dy, xmin, xmax, ymin, ymax, Tf, CFL;
double dt_imp=0.1;
int Nx=50, Ny=50, cas=2;

void initialiser_parametres() {
    switch (cas) {
        case 1:
            xmin=-1.0;
            xmax=1.0;
            ymin=-0.5;
            ymax=0.5;
            Tf=2.0;///normalement c 2
            dx=(xmax-xmin)/Nx;
            dy=(ymax-ymin)/Ny;
            break;
        case 2:
            xmin=-1.0;
            xmax=1.0;
            ymin=-1.0;
            ymax=1.0;
            Tf=2;
            dx=(xmax-xmin)/Nx;
            dy=(ymax-ymin)/Ny;
            break;
        case 3:
            xmin=-1;
            xmax=1;
            ymin=-1;
            ymax=1;
            Tf=1;
            dx=(xmax-xmin)/Nx;
            dy=(ymax-ymin)/Ny;
            break;
        default:
            printf("choisir un cas!");
            break;
    }
}