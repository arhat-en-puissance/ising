/*
Ce fichier contient l'algorithme établi suivant la méthode MCMC
Il permet les matrices au cours du temps et donc de visualiser l'évolution des domaines de Weiss.
*/






#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>

// Variables

FILE *fmatev; //fichier de l'évolution de la matrice au cours des itérations de l'algorithme de métropolis


const int mat_size = 100; // taille de l'echantillon de matériau simulé
const long int N = 1e8;//nombre d'itérations de la boucle
const double h = 0; // champ extérieur -> énergie sur site
const double J = 1; // Constante couple -> énergie d'interaction
const double T = 0.5; // Température


// Prototypes de fonctions


void fprint_mat(int matr_size, __int8 mat[matr_size][matr_size]);
int gen_rd_mat(int mat_size, __int8 mat[mat_size][mat_size]);
double energy(int m_size, __int8 mat[m_size][m_size], double h, double J);



// Main

int main(){


    fmatev = fopen("matev2.dat", "w");

    __int8 mat[mat_size][mat_size];

    int n_up = gen_rd_mat(mat_size, mat); //Attribue une valeur de -1 et 1 aux éléments de mat en renvoie le nombre de 1


    srand(time(NULL)); // On initialise la fonction rand() avec une seed qui dépend de l'heure

    double ener;



    // Boucle principale, Monte Carlo

    long int accept_1 = 0, accept_2 = 0, refus = 0; //Compteurs de matrices acceptées et refusées


    for(long int n=0; n<N;n++){

        int i = rand()%(mat_size);
        int j = rand()%(mat_size);

        double delta_e = 2*mat[i][j] * ( // la somme sur les plus proches voisin de [i][j]
                                    J*(mat[(i+1)%mat_size][j]// Ici il n'y a plus le /2 car on ne compte chaque couple bien qu'une fois
                                    + mat[(i-1+mat_size)%mat_size][j] //on tient compte de la périodicité
                                    + mat[i][(j+1)%mat_size]
                                    + mat[i][(j-1+mat_size)%mat_size])

                                    +h);


        if(delta_e < 0){
                    mat[i][j]*=-1;
                    accept_1++;
                    ener+=delta_e;
                }
        else{//delta>=0
            double proba = exp(-delta_e/T);
            double tirage =  (float)rand() / (float)RAND_MAX;
            if(tirage<proba){
                mat[i][j]*=-1;
                accept_2++;
                ener+=delta_e;
                }
            else{
                refus++;
            }
        }



        if(n%(N/100)==0){
            printf("ntot = %ld, avancee = %ld\n",N*n, n*100/N);

            for(int i = 0; i<=mat_size-1;i++){//écrit la matrice la matrice dans un fichier à chaque centième du processus d'équlibrage par métropolis
                for(int j = 0; j<=mat_size-1;j++){//les matrices sont séparées par python
                    fprintf(fmatev, "%d\t", mat[i][j]);
                }
                fprintf(fmatev, "\n");
            }
        }


    }



    fclose(fmatev);

    printf("accept_1 = %ld, accept_2 = %ld, refus = %ld, total = %ld \n", accept_1, accept_2, refus, accept_1+accept_2+refus);

    fprint_mat(mat_size, mat);

    return 0;

}








// Fonctions


void fprint_mat(int matr_size, __int8 mat[matr_size][matr_size]){// Enregistre une matrice dans un fichier
    FILE *mat_mag;
    mat_mag = fopen("mat_mag.dat","w");
    for(int i = 0; i<=matr_size-1;i++){
        for(int j = 0; j<=matr_size-1;j++){
            fprintf(mat_mag, "%d\t", mat[i][j]);
        }
        fprintf(mat_mag, "\n");
    }
    fclose(mat_mag);
}



int gen_rd_mat(int mat_size, __int8 mat[mat_size][mat_size]){// génère un matrice de 1 et -1 de taille mat_size
    int n = 0;
    for(int i = 0; i<=mat_size-1;i++){
        for(int j = 0; j<=mat_size-1;j++){
            __int8 s = (rand()%2)*2-1;
            mat[i][j]= s;
            if(s==1){
                n++;
            }

        }

    }
    return n;
}


double energy(int m_size, __int8 mat[m_size][m_size], double h, double J){//Calcule l'énergie du matériau
    double en = 0;
    for(int i=0; i<=m_size-1;i++){
        for(int j=0; j<=m_size-1;j++){
            en += -h*mat[i][j];
            en += -(J/2)*mat[i][j] * (mat[(i+1)%m_size][j]//On divise par 2 pour ne compter chaque couple qu'une fois
                      + mat[(i-1+m_size)%m_size][j]//On tient compte de la périodidicté
                      + mat[i][(j+1)%m_size]
                      + mat[i][(j-1+m_size)%m_size]);
        }
    }
    return en;
}
