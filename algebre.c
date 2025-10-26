#include "algebre.h"
#include "math.h"
#include "parametre.h"
#include <stdlib.h> 
#include "fonction.h"
#include <stdio.h>
#include </usr/lib/x86_64-linux-gnu/openmpi/include/mpi.h>
#include "charge.h"
extern int space_scheme, time_scheme;
extern double dx, dy, xmin, xmax, ymin, ymax, Tf, CFL;
extern int Nx, Ny, cas;


double* produit_MV( double* vecteur) {

    


    int nproc,rank;
    MPI_Comm_size(MPI_COMM_WORLD,&nproc);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);

    
    int iBeg,iEnd;
    charge_c(rank,Ny*Nx,nproc ,&iBeg, &iEnd);
    int taille_me=iEnd-iBeg+1;

    double* resultat = (double*)malloc(taille_me* sizeof(double)); 
     double* vecteur_apres = (double*)malloc(Nx* sizeof(double)); 
     double *vecteur_avant = (double*)malloc(Nx* sizeof(double)); 
    // Calcul du produit matrice-vecteur
    if (vecteur== NULL) {
        
        return NULL;
    }
    else{

            switch(space_scheme) {


                
                case 1://centré

 { 


                MPI_Send(&vecteur[0],(Nx),MPI_DOUBLE,modulo_pos(rank-1,nproc),0,MPI_COMM_WORLD);
                MPI_Send(&vecteur[taille_me-Nx],(Nx),MPI_DOUBLE,modulo_pos(rank+1,nproc),0,MPI_COMM_WORLD);
                MPI_Recv(&vecteur_apres[0],(Nx),MPI_DOUBLE,modulo_pos(rank+1,nproc),0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
                MPI_Recv(&vecteur_avant[0],(Nx),MPI_DOUBLE,modulo_pos(rank-1,nproc),0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
             
                 int k=0 ;
                if (iBeg==0){k=0;}
                else if (modulo_pos(iBeg, Nx) == 0) {k=(iBeg/Nx-1); } 
                else{k=iBeg/Nx;} ;

                                                                  

 
                    for (int I=iBeg;I<iEnd+1;I++){
                                                   
                       
                        if ((I%Nx==0 ) && (I>0) ){
                            k=k+1;
                        }
                    
                        int I_plus,I_moins,I_plus_x,I_moins_x;

                        double val_I_plus,val_I_moins,val_I_plus_x,val_I_moins_x;



                        I_plus=modulo_pos(I+1,Nx)+k*Nx;
                        I_moins=modulo_pos(I-1,Nx)+k*Nx;
                        I_plus_x=modulo_pos(I+Nx,Nx*Ny);
                        I_moins_x=modulo_pos(I+(Ny-1)*Nx,Nx*Ny);
                        

                       

                        
                        ///recherhce d'indice dans les tableaux spécifiques

                    get_true_value(I_plus, iBeg, iEnd, vecteur, vecteur_avant, vecteur_apres, &val_I_plus);
                    get_true_value(I_moins, iBeg, iEnd, vecteur, vecteur_avant, vecteur_apres, &val_I_moins);
                    get_true_value(I_plus_x, iBeg, iEnd, vecteur, vecteur_avant, vecteur_apres, &val_I_plus_x);
                    get_true_value(I_moins_x, iBeg, iEnd, vecteur, vecteur_avant, vecteur_apres, &val_I_moins_x);
         
 
                      resultat[I-iBeg]=0.5*alpha(maillage(I,0),maillage(I,1))*val_I_moins-0.5*alpha(maillage(I,0),maillage(I,1))*val_I_plus-0.5*beta(maillage(I,0),maillage(I,1))*val_I_plus_x+0.5*beta(maillage(I,0),maillage(I,1))*val_I_moins_x;



                    }

                    break;
                   }

     


   
    
                case 2://Upwind

                   { 



MPI_Send(&vecteur[0],(Nx),MPI_DOUBLE,modulo_pos(rank-1,nproc),0,MPI_COMM_WORLD);
MPI_Send(&vecteur[taille_me-Nx],(Nx),MPI_DOUBLE,modulo_pos(rank+1,nproc),0,MPI_COMM_WORLD);
MPI_Recv(&vecteur_apres[0],(Nx),MPI_DOUBLE,modulo_pos(rank+1,nproc),0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
MPI_Recv(&vecteur_avant[0],(Nx),MPI_DOUBLE,modulo_pos(rank-1,nproc),0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
             
                 int k=0 ;
                if (iBeg==0){k=0;}
                else if (modulo_pos(iBeg, Nx) == 0) {k=(iBeg/Nx-1); } 
                else{k=iBeg/Nx;} ;

                                                                  

 
                    for (int I=iBeg;I<iEnd+1;I++){
                                                   
                       
                        if ((I%Nx==0 ) && (I>0) ){
                            k=k+1;
                        }
                    
                        int I_plus,I_moins,I_plus_x,I_moins_x;

                        double val_I_plus,val_I_moins,val_I_plus_x,val_I_moins_x,val_I;



                        I_plus=modulo_pos(I+1,Nx)+k*Nx;
                        I_moins=modulo_pos(I-1,Nx)+k*Nx;
                        I_plus_x=modulo_pos(I+Nx,Nx*Ny);
                        I_moins_x=modulo_pos(I+(Ny-1)*Nx,Nx*Ny);
                        

                       

                        
                        ///recherhce d'indice dans les tableaux spécifiques

                    get_true_value(I_plus, iBeg, iEnd, vecteur, vecteur_avant, vecteur_apres, &val_I_plus);
                    get_true_value(I_moins, iBeg, iEnd, vecteur, vecteur_avant, vecteur_apres, &val_I_moins);
                    get_true_value(I_plus_x, iBeg, iEnd, vecteur, vecteur_avant, vecteur_apres, &val_I_plus_x);
                    get_true_value(I_moins_x, iBeg, iEnd, vecteur, vecteur_avant, vecteur_apres, &val_I_moins_x);
                    get_true_value(I, iBeg, iEnd, vecteur, vecteur_avant, vecteur_apres, &val_I);
                   



                        if( (alpha(maillage(I,0),maillage(I,1))>=0.)&&((beta(maillage(I,0),maillage(I,1)))>=0.))// vx>0 Vy>0
                        {
                            


                        
                        resultat[I-iBeg]=alpha(maillage(I,0),maillage(I,1))*val_I_moins+gama(maillage(I,0),maillage(I,1))*val_I+beta(maillage(I,0),maillage(I,1))*val_I_moins_x;

                        }




                        else if( (alpha(maillage(I,0),maillage(I,1))<0.)&&((beta(maillage(I,0),maillage(I,1)))>=0.))// vx<0 Vy>0
                        {
                            double s=alpha(maillage(I,0),maillage(I,1))-beta(maillage(I,0),maillage(I,1));
                        resultat[I-iBeg]=-alpha(maillage(I,0),maillage(I,1))*val_I_plus+s*val_I+beta(maillage(I,0),maillage(I,1))*val_I_moins_x;
                        }


                        else if( (alpha(maillage(I,0),maillage(I,1))>=0.)&&((beta(maillage(I,0),maillage(I,1)))<0.))// vx>0 Vy<0
                        {
                            
                       double s=-alpha(maillage(I,0),maillage(I,1))+beta(maillage(I,0),maillage(I,1));
                            resultat[I-iBeg]=alpha(maillage(I,0),maillage(I,1))*val_I_moins+s*val_I-beta(maillage(I,0),maillage(I,1))*val_I_plus_x;
                        }
                        else if  ( (alpha(maillage(I,0),maillage(I,1))<0.)&&((beta(maillage(I,0),maillage(I,1)))<0.))//vx<0 et vy<0
                        {

                          resultat[I-iBeg]=-alpha(maillage(I,0),maillage(I,1))*val_I_plus-gama(maillage(I,0),maillage(I,1))*val_I-beta(maillage(I,0),maillage(I,1))*val_I_plus_x;
                        }





                    }

                    break;
                   }

     }

    free( vecteur_apres); 
     free(vecteur_avant); 
   
    }
    return resultat; // Retourne le pointeur vers le vecteur résultat





}
int indexe(int i,int j){ // focntion qui donne l'indice de l'element stocké en fonction de i et j
    return (j-1)*Nx+i-1;
} 
int couple(int I,int axis){ //focntion inverse de index qui donne i ou j selon axe choisi
    int result;
    if (axis==0) 
    {
        result= I%Nx+1;
    }
    else if (axis==1)
    {
        result= I/Nx+1;
    }
    return result;
}

double maillage(int I,int axis){ // fonction qui donne xi ou xj selon axe choisi
    double result;
    if (axis==0) 
    {
         result=xmin+couple(I,0)*dx;
    }
    else if (axis==1)
    {
        result=ymin+couple(I,1)*dy;
    }
    return result;
    
}

double produitScalaire(double *A, double *B, int taille) {
    double produit = 0;
    for (int i = 0; i < taille; i++) {
        produit += A[i] * B[i];
    }
    return produit;
}
double norm(double *A, int taille) {
    double produit = 0.0;
    double resultat ;
    for (int i = 0; i < taille; i++) {
        produit += A[i] * A[i];
    }
    resultat=sqrt(produit);
    return resultat;
}

void copierTableau(double *source, double *destination, int taille) {
    for (int i = 0; i < taille; i++) {
        destination[i] = source[i];
    }
}

void sommeTableaux(double *A, double *B, double *C, int taille) {
    for (int i = 0; i < taille; i++) {
        C[i] = A[i] + B[i];
    }
}

void differenceTableaux(double *A, double *B, double *C, int taille) {
    for (int i = 0; i < taille; i++) {
        C[i] = A[i] - B[i];
    }
}

void scalaireMultiplieTableau(double scalaire, double *vecteur, double *resultat, int taille) {
    for (int i = 0; i < taille; i++) {
        resultat[i] = scalaire * vecteur[i];
    }
}
/////!!!!!problem d'indice
void get_true_index(int I, int iBeg, int iEnd, int *i_loc) {

     int nproc,rank;
    MPI_Comm_size(MPI_COMM_WORLD,&nproc);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);

    if (rank==0)
    {
 if ((I <= iEnd) && (I >= iBeg)) {
        *i_loc = I - iBeg;
    } else if (I > iEnd) {
        *i_loc = I - iEnd - 1;
    } else if((I<Nx*Ny)&&(I>Nx*Ny-Nx-1) ){
        *i_loc = I - Nx*Ny +Nx;
    }
    }
    else if(rank==nproc-1)
    {
 if ((I <= iEnd) && (I >= iBeg)) {
        *i_loc = I - iBeg;
    } else if ((I>= 0)&&(I<Nx ) ){
        *i_loc = I;
    } else if (I < iBeg) {
        *i_loc = I - iBeg + Nx + 1;
    }
    }
    else{

    if ((I <= iEnd) && (I >= iBeg)) {
        *i_loc = I - iBeg;
    } else if (I > iEnd) {
        *i_loc = I - iEnd - 1;
    } else if (I < iBeg) {
        *i_loc = I - iBeg + Nx + 1;
    }
    }

    
}

void get_true_value(int I, int iBeg, int iEnd, double* vecteur, double* vecteur_avant, double* vecteur_apres, double* value) {

    int i_loc;
    get_true_index(I, iBeg, iEnd, &i_loc);

    int nproc,rank;
    MPI_Comm_size(MPI_COMM_WORLD,&nproc);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);

    if (rank==0)
    {
 if ((I <= iEnd) && (I >= iBeg)) {
        *value = vecteur[i_loc];
    } else if (I > iEnd) {
        *value = vecteur_apres[i_loc];
    }else if((I<Nx*Ny)&&(I>Nx*Ny-Nx-1) ){
        *value = vecteur_avant[i_loc];
    }
    }
    else if(rank==nproc-1)
    {
 if ((I <= iEnd) && (I >= iBeg)) {
        *value = vecteur[i_loc];
    } else if ((I>= 0)&&(I<Nx ) ){
       *value = vecteur_apres[i_loc];
    } else if (I < iBeg) {
        *value = vecteur_avant[i_loc];
    }
    }
    else{

    if ((I <= iEnd) && (I >= iBeg)) {
        *value = vecteur[i_loc];
    } else if (I > iEnd) {
       *value = vecteur_apres[i_loc];
    } else if (I < iBeg) {
        *value = vecteur_avant[i_loc];
    }
    }



    
    
}