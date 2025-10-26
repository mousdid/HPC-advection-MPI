#include "stdio.h"
#include "stdlib.h"
#include "fonction.h"
#include "algebre.h"
#include "BiCGstab.h"
#include "parametre.h"
#include "math.h"
#include </usr/lib/x86_64-linux-gnu/openmpi/include/mpi.h>
#include "charge.h"
extern int space_scheme, time_scheme;
extern double dx, dy, xmin, xmax, ymin, ymax, Tf, CFL;
extern int Nx, Ny, cas;

int main(int argc, char* argv[]) {

    
    MPI_Init(&argc,&argv);
int nproc,rank;
    MPI_Comm_size(MPI_COMM_WORLD,&nproc);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    
   

   

    int iBeg,iEnd;
   
  
charge_c(rank,Nx*Ny,nproc ,&iBeg, &iEnd);
int taille_me=iEnd-iBeg+1;
    initialiser_parametres();
  
    // testing

    //  double dt_test=0.01;
    //  double* x0 = (double*)malloc(Nx*Ny* sizeof(double));//ce vecteur permet l'affichage d'une colonne de la matrice
    // for (int i=0 ;i<Nx*Ny ;i++)
    // {   
    //     if (i==3){
    //         x0[i]=1;
    //     }
    //     else{
    //     x0[i]=0.;
    //     }
        
    // }
    
    // for (int i=0 ;i<Nx*Ny ;i++)
    // {   
    //     if (space_scheme==1){//A=I-dt*B
    //         //printf("produit implicit: %f\n",x0[i]-dt_test*produit_MV(x0)[i]);
            
    //     }
    //     else if(space_scheme==2){//A=I+dt*B
    //         //printf("produit explicit: %f\n",x0[i]+dt_test*produit_MV(x0)[i]);
           
    //     }
    // }








//if (rank ==0)
//{






    //initialisation
 double* vect_u_me = (double*)malloc(taille_me* sizeof(double));
    for (int i=iBeg ;i<iEnd+1;i++)
    {
        vect_u_me[i-iBeg]=u0(maillage(i,0),maillage(i,1));/////local
   
    }
     double* vect_un_me = (double*)malloc(taille_me * sizeof(double));///// solution_temp




    switch(time_scheme) {
    case 1: //Explicit
        {
            //Trouver la cfl


            //trouver V_x max
            double max_vx; // Déclaration de max_vx au début du bloc de code
            double current_vx;
            max_vx = fabs(v(xmin, ymin)[0]); // Initialisation avec le coin en bas à gauche
           for (int I=0;I<Nx*Ny;I++){
            
           
                    current_vx = fabs(v( maillage(I, 0), maillage( I, 1))[0]);
                    
                    if (current_vx > max_vx) {
                        max_vx = current_vx;
                    }
                }
            
            //trouver V_y max
           double max_vy; // Déclaration de max_vx au début du bloc de code
            double current_vy;
            max_vy = fabs(v(xmin, ymin)[1]); // Initialisation avec le coin en bas à gauche
           for (int I=0;I<Nx*Ny;I++){
            
           
                    current_vy = fabs(v( maillage(I, 0), maillage( I, 1))[1]);
                   
                    if (current_vy > max_vy) {
                        max_vy = current_vy;
                    }
                }
            

            CFL = (max_vx / dx) + (max_vy / dy);
          
           
            double dt = 0.9 / CFL;
            printf("dt=%f\n",dt);
            
        ///ecriture initiale
            char filenam[120]; // Définissez la taille en fonction du format du nom du fichier
            sprintf(filenam, "sol.%d_proc_%d.dat", 0,rank); // Formatage du nom du fichier
            FILE* file_explicit = fopen(filenam, "w"); // Ouvrir le fichier en écriture
            if (file_explicit == NULL) {
                printf("Erreur lors de la création du fichier %s\n", filenam);
            return 1;
            }
            //stocker le vecteur initial
            for (int i = iBeg; i < iEnd+1; i++) {
                fprintf(file_explicit, "%f %f %f\n", maillage(i, 0), maillage(i, 1),vect_u_me[i-iBeg]);

            }
            fprintf(file_explicit, "\n"); // Saut de ligne pour passer à l'itération de temps suivante
            fclose(file_explicit); // Fermer le fichier
            




            // Boucle temporelle pour l'explicite
            double T_explicit = 0;
            int Nt=0;
            while (T_explicit < Tf) {
                // Mise à jour du vecteur de résolution pour la prochaine étape de temps
                copierTableau(produit_MV(vect_u_me), vect_un_me, taille_me);
                scalaireMultiplieTableau(dt, vect_un_me, vect_un_me, taille_me);
                sommeTableaux(vect_u_me, vect_un_me, vect_u_me, taille_me);
                T_explicit += dt; // Incrémentation du temps
                Nt+=1;
               
                    
                


///ecriture en temps final

    
                


                char filename[120]; // Définissez la taille en fonction du format du nom du fichier
                sprintf(filename, "sol.%d_proc_%d.dat", Nt,rank); // Formatage du nom du fichier
                FILE* file_explicit_final = fopen(filename, "w"); // Ouvrir le fichier en écriture
                if (file_explicit_final == NULL) {
                    printf("Erreur lors de la création du fichier %s\n", filename);
                return 1;
                }
                for (int i = iBeg; i < iEnd+1; i++) {//le proc 0 qui se charge de la solution
                    fprintf(file_explicit_final, "%f %f %f\n", maillage(i, 0), maillage(i, 1),vect_u_me[i-iBeg]);
                }
                fprintf(file_explicit_final, "\n"); // Saut de ligne pour passer à l'itération de temps suivante
                fclose(file_explicit_final); // Fermer le fichier


            }

            

            

            
            
        break;
        }



    case 2: //Implicit
    {



     ///ecriture initiale
            char filenam[120]; // Définissez la taille en fonction du format du nom du fichier
            sprintf(filenam, "sol.%d_proc_%d.dat", 0,rank); // Formatage du nom du fichier
            FILE* file_explicit = fopen(filenam, "w"); // Ouvrir le fichier en écriture
            if (file_explicit == NULL) {
                printf("Erreur lors de la création du fichier %s\n", filenam);
            return 1;
            }
            //stocker le vecteur initial
            for (int i = iBeg; i < iEnd+1; i++) {
                fprintf(file_explicit, "%f %f %f\n", maillage(i, 0), maillage(i, 1),vect_u_me[i-iBeg]);

            }
            fprintf(file_explicit, "\n"); // Saut de ligne pour passer à l'itération de temps suivante
            fclose(file_explicit); // Fermer le fichier
            







 












            // Boucle temporelle pour l'explicite
            double T_explicit = 0;
            int Nt=0;
            while (T_explicit < Tf) {
                
                // Résolution du système d'équations pour la prochaine étape de temps
                bicgstab(vect_u_me, vect_un_me, taille_me, 1E-6, 5000);
                copierTableau(vect_un_me,vect_u_me,taille_me);
                T_explicit += dt_imp; // Incrémentation du temps
                Nt+=1;
               
                    
                



///ecriture en temps final

    
                


                char filename[120]; // Définissez la taille en fonction du format du nom du fichier
                sprintf(filename, "sol.%d_proc_%d.dat", Nt,rank); // Formatage du nom du fichier
                FILE* file_explicit_final = fopen(filename, "w"); // Ouvrir le fichier en écriture
                if (file_explicit_final == NULL) {
                    printf("Erreur lors de la création du fichier %s\n", filename);
                return 1;
                }
                for (int i = iBeg; i < iEnd+1; i++) {//le proc 0 qui se charge de la solution
                    fprintf(file_explicit_final, "%f %f %f\n", maillage(i, 0), maillage(i, 1),vect_u_me[i-iBeg]);
                }
                fprintf(file_explicit_final, "\n"); // Saut de ligne pour passer à l'itération de temps suivante
                fclose(file_explicit_final); // Fermer le fichier
            }  
            break;  
            }
       
    }
   
    free(vect_u_me);
    free(vect_un_me);




 
MPI_Finalize();








   
    return 0;
    }




































