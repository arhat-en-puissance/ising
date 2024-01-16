/*
Ce fichier réutilise le contenu du fichier evolution_mcmc.c dans la fonction metro()
D'autres fonctions sont communes.
Il permet de stocker la valeur de l'aimantation par site moyenne en fonction du champ extérieur et cela pour différentes températures.
Ce programme a tourné huit heures pour parcourir toutes les températures.
Si vous souhaitez le tester je vous conseille de supprimer la boucle sur la température
et de tester seulement pour une valeur comme T = 1 par exemple.
*/




#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>

// Variables


FILE *fmat_h; //état du réseau après metropolis pour chaque pas de champ extérieur
FILE *faim; //fichier aimantation par site en fonction du champ extérieur qui augmente
FILE *faim2; //fichier aimantation par site en fonction du champ extérieur qui diminue

const int mat_size = 30; // taille de l'echantillon de matériau simulé
const long int N = 1e8;//nombre d'itérations de la boucle
const double J = 1; // Constante couple -> énergie d'interaction

// Prototypes de fonctions

void fprint_mat(int *mat);
int gen_rd_mat(int *mat);
double energy(int *mat, double h, double J);
void metro(int *mat, double h, double J, double T);
double aimantation_par_site(int *mat);



// Main

int main(){


    int N_T = 26;
    double T_min = 0.25, T_max = 3.5;
    double pas_T = (T_max - T_min)/N_T;

    int N_h = 64;
    double h_min = -2, h_max = 2;
    double pas_h = (h_max-h_min)/N_h;

    srand(time(NULL)); // On initialise la fonction rand()
    fmat_h =fopen("mat_h.dat", "w");
    faim = fopen("aimantation_h.dat", "w");
    faim2 = fopen("aimantation_h2.dat", "w");


    int *mat;


    for(int i_T = 0; i_T<=N_T;i_T++){

        int nd = 100+i_T;
        double T = T_min + i_T*pas_T;

        printf("T = %f\n",T);
        mat = calloc(mat_size*mat_size, sizeof(int));

        int n_up = gen_rd_mat(mat);

        fprintf(faim, "%f \t %d", T, nd);
        fprintf(faim, "\n");

        //Augmentation du champ extérieur

        for(int i_h=0;i_h<=N_h;i_h++){
            double h = h_min+i_h*pas_h;

            double ener;

            metro(mat, h,J,T);
            fprintf(faim, "%f \t %f", h, aimantation_par_site(mat));
            fprintf(faim, "\n");


        }


        //Diminution du champ extérieur

        for(int i_h=0;i_h<=N_h;i_h++){
            double h = h_max-i_h*pas_h;

            printf("h = %f\n", h);



            double ener;

            metro(mat, h,J,T);
            fprintf(faim2, "%f \t %f", h, aimantation_par_site(mat));
            fprintf(faim2, "\n");

        }
    }

    free(mat);

    fclose(faim);
    fclose(faim2);

    return 0;

}








// Fonctions


void metro(int *mat, double h, double J, double T){


    // Boucle principale, Monte Carlo

    long int accept_1 = 0, accept_2 = 0, refus = 0; //Compteurs de matrices acceptées et refusées
    for(long int n=0; n<N;n++){


        int i = rand()%(mat_size);
        int j = rand()%(mat_size);

        double delta_e = 2* (*(mat+i+j*mat_size)) * (
                                    J*(  *(mat+(i+1)%mat_size+j*mat_size)  // la somme sur les plus proches voisin de [i][j]
                                    + *(mat+(i-1+mat_size)%mat_size+j*mat_size)
                                    + *(mat+i+((j+1)%mat_size)*mat_size)
                                    + *(mat+i+((j-1+mat_size)%mat_size)*mat_size))

                                    +h);


        if(delta_e < 0){
                    (*(mat+i+j*mat_size))*=-1;
                    accept_1++;
                }
        else{//delta>=0
            double proba = exp(-delta_e/T);
            double tirage =  (float)rand() / (float)RAND_MAX;
            if(tirage<proba){
                (*(mat+i+j*mat_size))*=-1;
                accept_2++;
                }
            else{
                refus++;
            }
        }

    }

}







int gen_rd_mat(int *mat){// génère un matrice de 1 et -1 de taille mat_size
    int n = 0;
    for(int i = 0; i<=mat_size-1;i++){//On s'arrête avant la dernière colonne et ligne qui servent à forcer la périodicité
        for(int j = 0; j<=mat_size-1;j++){
            int s = (rand()%2)*2-1;
            *(mat+i+mat_size*j)= (rand()%2)*2-1;
            if(s==1){
                n++;
            }

        }

    }
    return n;
}


double energy(int *mat, double h, double J){//Calcule l'énergie du matériau
    double en = 0;
    for(int i=0; i<=mat_size-1;i++){
        for(int j=0; j<=mat_size-1;j++){
            en += -h* ( *(mat+i+j*mat_size));
            en += -(J/2)* ( *(mat+i+j*mat_size)) * ( *(mat+(i+1)%mat_size+mat_size*j)
                      + *(mat+(i-1+mat_size)%mat_size+mat_size*j)
                      + *(mat+i+((j+1)%mat_size)*mat_size)
                      + *(mat+i+((j-1+mat_size)%mat_size)*mat_size));
        }
    }
    return en;
}



double aimantation_par_site(int *mat){

    double aimantation = 0;

    for(int i = 0; i<mat_size;i++){
        for(int j = 0; j<mat_size;j++){
            aimantation += *(mat+i+mat_size*j);
        }
    }
    aimantation /= mat_size*mat_size;
    return aimantation;
}

