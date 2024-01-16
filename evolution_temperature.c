/*
Ce fichier réutilise le contenu du fichier evolution_mcmc.c dans la fonction metro()
D'autres fonctions sont communes.
Il permet de stocker les grandeurs intéressantes en fonction de la température.
*/





#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>

// Variables


FILE *fmat_t; //état du réseau après metropolis pour chaque pas de température
FILE *faim; //fichier aimatation par site en fonction de la température
FILE *fenps;//fichier energie par site en fonction de la température
FILE *fsus; //fichier susceptibilité magnétique en fonction de la température

const int mat_size = 30; // taille de l'echantillon de matériau simulé
const long int N = 1e8;//nombre d'itérations de la boucle
const double h = 0; // champ extérieur -> énergie sur site
const double J = 0.5; // Constante couple -> énergie d'interaction


// Prototypes de fonctions


void fprint_mat(int *mat);
int gen_rd_mat(int *mat);
double energy(int *mat, double h, double J);
void metro(int *mat, double h, double J, double T, double * aimantations);
double aimantation_par_site(int *mat);




// Main

int main(){




    int N_T = 20;
    double T_min = 0.1, T_max = 5.1;
    double pas_T = (T_max-T_min)/N_T;


    fmat_t =fopen("mat_T.dat", "w");
    faim = fopen("aimantation_T.dat", "w");
    fenps = fopen("energie_T.dat", "w");
    fsus = fopen("susceptibilite_T.dat","w");

    for(int i_T=0;i_T<=N_T;i_T++){
        double T = T_min+i_T*pas_T;

        printf("T = %f\n", T);

        int *mat;

        mat = calloc(mat_size*mat_size, sizeof(int));

        int n_up = gen_rd_mat(mat); //Attribue une valeur de -1 et 1 aux éléments de mat en renvoie le nombre de 1


        srand(time(NULL)); // On initialise la fonction rand()

        double* aimantations;
        aimantations = calloc(990, sizeof(double));
        double susceptibilite_mag = 0, aimantation_moyenne = 0;

        for(int k = 0; k<=990;k++){
            aimantation_moyenne += *(aimantations+k);
        }

        aimantation_moyenne /= 990;

        for(int k = 0; k<=990;k++){
            susceptibilite_mag += (*(aimantations+k)-aimantation_moyenne) * (*(aimantations+k)-aimantation_moyenne);
        }

        susceptibilite_mag /= 990;

        metro(mat, h,J,T, aimantations);
        fprintf(faim, "%f \t %f", T, aimantation_moyenne);
        fprintf(faim, "\n");

        fprintf(fenps, "%f \t %f", T, energy(mat, h, J)/(mat_size*mat_size)); //énergie totale divisée par le nombre de sites
        fprintf(fenps, "\n");

        fprintf(fsus, "%f \t %f", T, susceptibilite_mag/T);
        fprintf(fsus, "\n");

        for(int i = 0; i<=mat_size-1;i++){//écrit la matrice la matrice dans un fichier à chaque centième du processus d'équlibrage par métropolis
            for(int j = 0; j<=mat_size-1;j++){//les matrices sont séparées par python
                fprintf(fmat_t, "%d\t", *(mat+i+mat_size*j));
            }
            fprintf(fmat_t, "\n");
        }

        free(mat);
    }





    return 0;

}








// Fonctions



void metro(int *mat, double h, double J, double T, double * aimantations){


    int test = 0;

    // Boucle principale, Monte Carlo

    long int accept_1 = 0, accept_2 = 0, refus = 0; //Compteurs de matrices acceptées et refusées
    for(long int n=0; n<N;n++){


        int i = rand()%(mat_size);
        int j = rand()%(mat_size);

        double delta_e = 2* (*(mat+i+j*mat_size)) * (// la somme sur les plus proches voisin de [i][j]
                                    J*(  *(mat+(i+1)%mat_size+j*mat_size) // en tenant compté de la périodicité
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
        if(n%(N/10)==0){
            printf("%ld/%d\n", n*10/N +1, 10);

        }

        int k = n*1000/N;

        if(n%(N/1000)==0 && k>=10){
            *(aimantations+k-10) += aimantation_par_site(mat);
            test ++;
            printf("k = %d \n", k);

        }

    }




    fprint_mat(mat);

    printf("test = %d \n", test);
}




void fprint_mat(int *mat){// Enregistre une matrice dans un fichier
    FILE *mat_mag;
    mat_mag = fopen("mat_mag.dat","w");
    for(int i = 0; i<=mat_size-1;i++){
        for(int j = 0; j<=mat_size-1;j++){
            fprintf(mat_mag, "%d\t", *(mat+i+mat_size*j));
        }
        fprintf(mat_mag, "\n");
    }
    fclose(mat_mag);
}



int gen_rd_mat(int *mat){// génère un matrice de 1 et -1 de taille mat_size
    int n = 0;
    for(int i = 0; i<=mat_size-1;i++){
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
            en += -(J/2)* ( *(mat+i+j*mat_size)) * ( *(mat+(i+1)%mat_size+mat_size*j)// C'est ici qu'intervient la périodicité aux bords
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
