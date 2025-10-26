#include "stdio.h"
#include "stdlib.h"
#include "fonction.h"
#include "algebre.h"
#include "BiCGstab.h"
#include "parametre.h"
#include "math.h"
#include <mpi.h>
extern int space_scheme, time_scheme;
extern double dx, dy, xmin, xmax, ymin, ymax, Tf, CFL;
extern int Nx, Ny, cas;

void charge_c(int me, int n, int np, int *iBeg, int *iEnd);



int main(int argc, char* argv[]) {
    initialiser_parametres();
    int n=Nx*Ny;
    // Find the CFL
    double max_vx = fabs(v(xmin, ymin)[0]); // Initialization with the bottom left corner
    double current_vx;
    for (int I = 0; I < Nx * Ny; I++) {
        current_vx = fabs(v(maillage(I, 0), maillage(I, 1))[0]);
        if (current_vx > max_vx) {
            max_vx = current_vx;
        }
    }

    double max_vy = fabs(v(xmin, ymin)[1]); // Initialization with the bottom left corner
    double current_vy;
    for (int I = 0; I < Nx * Ny; I++) {
        current_vy = fabs(v(maillage(I, 0), maillage(I, 1))[1]);
        if (current_vy > max_vy) {
            max_vy = current_vy;
        }
    }

    CFL = (max_vx / dx) + (max_vy / dy);
    double dt;
    dt = 0.9 / CFL;

    // Initialize MPI
    int nproc, rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    
    printf("i'm proc %d and i get dt=%f\n",rank, dt);

    int iBeg, iEnd;
    charge_c(rank, n, nproc, &iBeg, &iEnd);
    printf("ibeg=%d",iBeg);
    int taille_me = iEnd - iBeg + 1;
    int taille_me2 = taille_me + 2 * Nx;
    double *u_bas = (double *)malloc(Nx * sizeof(double));
    double *u_haut = (double *)malloc(Nx * sizeof(double));
    double *vect_u_me = (double *)malloc(taille_me * sizeof(double));
    double *vect_u = (double *)malloc(taille_me2 * sizeof(double));
    double *vect_un = (double *)malloc(taille_me2 * sizeof(double));
    double *U_global;


    if (rank == 0){
        U_global = (double*)malloc(Nx * Ny * sizeof(double));
    }

    for (int i = iBeg; i < iEnd; i++) {
        vect_u[i - iBeg + Nx] = u0(maillage(i, 0), maillage(i, 1));
        printf("%f\n",vect_u[i - iBeg + Nx]);
    }
    for (int i = 0; i < Nx; i++) {
        u_bas[i] = u0(maillage(i + iBeg - Nx, 0), maillage(i + iBeg - Nx, 1));
        u_haut[i] = u0(maillage(i + iBeg + taille_me, 0), maillage(i + iBeg + taille_me, 1));
    }
    for (int i = 0; i < Nx; i++) {
        vect_u[i] = u_bas[i];
        vect_u[taille_me + i] = u_haut[i];
    }

    if (rank == 0) {
            MPI_Send(&vect_u_me[0], Nx, MPI_DOUBLE, nproc -1, 0, MPI_COMM_WORLD);
            MPI_Send(&vect_u_me[taille_me-Nx], Nx, MPI_DOUBLE, rank +1 , 0, MPI_COMM_WORLD);
            MPI_Recv(&u_haut[0], Nx, MPI_DOUBLE, rank +1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&u_bas[0], Nx, MPI_DOUBLE, nproc -1 , 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
    else if (rank == nproc -1) {
            MPI_Send(&vect_u_me[0], Nx, MPI_DOUBLE, rank -1, 0, MPI_COMM_WORLD);
            MPI_Send(&vect_u_me[taille_me-Nx], Nx, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
            MPI_Recv(&u_haut[0], Nx, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&u_bas[0], Nx, MPI_DOUBLE, rank -1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
    else {
            MPI_Send(&vect_u_me[0], Nx, MPI_DOUBLE, rank -1, 0, MPI_COMM_WORLD);
            MPI_Send(&vect_u_me[taille_me-Nx], Nx, MPI_DOUBLE, rank +1, 0, MPI_COMM_WORLD);
            MPI_Recv(&u_haut[0], Nx, MPI_DOUBLE, rank +1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&u_bas[0], Nx, MPI_DOUBLE, rank -1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    for (int I = 0; I< Nx ; ++I) { 
            vect_u[I]=u_bas[I];
    }
    
    for (int I=iBeg; I<= iEnd; ++I) {
            vect_u[I-iBeg+Nx] = vect_u_me[I-iBeg];
    }

    for (int I = taille_me +Nx ; I< taille_me + 2*Nx   ; ++I) {
            vect_u[I]=u_haut[I - taille_me - Nx];
    }

    for (int I = iBeg; I <= iEnd; ++I){
        vect_un[I-iBeg] = produit_MV(vect_u, iBeg, iEnd)[I-iBeg];
        printf("proc %d , vect[%d] = %f\n",rank,I-iBeg,vect_un[I-iBeg+Nx]);
    }


    free(u_bas);
    free(u_haut);
    free(vect_un);
    free(vect_u);
    free(vect_u_me); 
    if (rank == 0){
        free(U_global);
    }

    MPI_Finalize();
    return 0;
}

void charge_c(int me,int n,int np ,int *iBeg, int *iEnd)
{
  int r;
    r=n%np;
    if (me<r) {
        *iBeg=me*(n/np+1);
        *iEnd=(me+1)*(n/np+1)-1;
    }
    else{
        *iBeg=r*(n/np+1)+(me-r)*(n/np);
        *iEnd=*iBeg+(n/np)-1;
    }
}