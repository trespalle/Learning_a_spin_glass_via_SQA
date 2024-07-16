#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define NormRANu (2.3283063671E-10F)
#define normNorm (4.656612873E-10F) ///1/2^31, para que r1279 nunca valga 1
#define NBITM1 31
#define dim 2
#define z 4


unsigned int larueda[256], number = 0;
unsigned char ind_ran, i1, i2, i3;
double AMIGOTE, PERESEUDO;
double bond_penalty;

int secuencia_rand[2048], index1[2048], index2[2048], iaux;


void ini_ran(int seed);
double Rand(void);

void crea_lista_vecinos(int** vecinos, int L, int N);
void muestra_lista(int **vecinos, int N);
void ini_r1279();
double r1279(void);

int sign(int n);
void read_D(char *namefile, int i_seed, char **D, int C, int N);
void comprobacion(char **D, int C, int N);
void create_couplings_list(int**vecinos, int N, int zmax, int**JJ, int M);
void tabulate_lambda (int **lambda, char **D, int C, int N, int **JJ, int M_times_R_guess, int R);
void tabulate_FUNK (double **FUNK, int zmax, double temp);
double compute_PLH( double **FUNK, int zmax, char **D, int **lambda, int C, int N, int R);
void MC_STEP(int **JJ, int Mmax, int *columna_i_new, int *columna_j_new, char **D, int **lambda, int C, int N, double **FUNK, int zmax, double teff, int R, double acoplo_2, double acoplo_3, double acoplo_4, double PLH_strength, int **all_coups, int *M_times_R_guess_dir, int M);
void print_couplings_mean(int **JJ, int M);
void create_name(char *namefile, int i_seed);
void gamma_stats(double *gamma_vector, int n_seeds);
void lee_input(char *namefile, int *N_dir, int *C_dir, int *n_temps_dir, int **temps_dir, int *n_seeds_dir, double *gamma_ini_dir, double*gamma_fin_dir, int *tau_dir, int *first_seed_dir, int *R_dir, int *n_wl_dir, double *teff_dir, double *PLH_strength_dir);
void change_g_name(char *name_g_file, int i_seed);
void read_geometry(char *namegfile, int *M_dir, int *zmax_dir, int **vecinos, int ***JJ_0_dir, int N);
void update_name(char *namefile, int temp);
void save_results(int *temps, double *av_gamma, double *av_gamma2, double *av_sensit, double *av_sensit2, double *av_precis, double *av_precis2, double *av_F1, double *av_F1_2, int n_temps, int n_seeds, int first_seed);
void tabulate_acoplos(double *acoplo_2_list, double *acoplo_3_list, double *acoplo_4_list, int tau, double gamma_ini, double delta_gamma, double temp, int R);
void make_the_guess(int **JJ, int **all_coups, int Mmax, int R);
double compute_success(int **JJ_0, int **JJ, int M);
void append_new_coupling(int **all_coups, int **JJ, int *M_times_R_guess_dir, int thechosen, int coupling, int tau);
void delete_coupling(int **all_coups, int **JJ, int *M_guess_dir, int thechosen, int N, int Mmax);
void WORLDLINE_FLIP(int **JJ, int Mmax, int **columna_i_new_wl, int **columna_j_new_wl, char **D, int **lambda, int C, int N, double **FUNK, int zmax, double teff, int R, double acoplo_2, double acoplo_3, double acoplo_4, double PLH_strength, int **all_coups, int *M_times_R_guess_dir, int M);
void random_initial_couplings(int **all_coups, int **JJ, int Mmax, int *Mguess_dir, double alpha, int R);
double GAMMA(int **real_all_coups, int **all_coups, int Mmax, int **JJ_0, int M, int N);
double F1_score(int **real_all_coups, int Mmax, int **all_coups, double *sens_dir, double *precis_dir);
void save_trotter_copies( int **JJ, int M_times_R_guess, char *filename);






int main()
{
    ini_ran(123456789);
    ini_r1279();
    double teff, gamma, delta_gamma, gamma_ini, gamma_fin, PLH_strength;
    int mu, i, j, k, e, N, C, M, R, zmax, tau, i_seed, n_seeds, n_wl, index;
    int n_temps, i_temps, first_seed, Mmax, M_times_R_guess;
    int *integer_temps;
    double aux_gamma, aux_success, aux_precis;
    double temperatura;


    ///Leemos los parámetros
    lee_input("input_UG_SQA.txt", &N, &C, &n_temps, &integer_temps, &n_seeds, &gamma_ini, &gamma_fin, &tau, &first_seed, &R, &n_wl, &teff, &PLH_strength);
    Mmax = N*(N-1)/2;

    printf("N = %d\nC = %d\nn_temps =%d\n", N, C, n_temps);
    for(i=0; i<n_temps; i++) printf("%d ", integer_temps[i]);
    printf("\nn_seeds = %d\ngamma_ini = %lf\ngamma_fin = %lf\ntau = %d\nfirst_seed = %d\nR = %d\nn_wl = %d\nteff = %lf\nPLH_strength = %lf\n", n_seeds, gamma_ini, gamma_fin, tau, first_seed, R, n_wl, teff, PLH_strength);
    printf("bond pen = %lf\n", bond_penalty);



    ///Generamos las principales matrices y vectores:
    int **vecinos = (int**) malloc(N*sizeof(int*));
    int **JJ_0; /// couplings reales (los que queremos inferir)
    int **JJ; /// variables libres de nuestro simulated annealing

    JJ = (int**) malloc(Mmax*R*sizeof(int*));
    for(i=0; i<Mmax*R; i++) JJ[i] = (int*) malloc(3*sizeof(int));
    int **all_coups = (int**) malloc(Mmax*R*sizeof(int*));
    for(i=0; i<Mmax*R; i++) all_coups[i] = (int*) malloc(3*sizeof(int));
    int **real_all_coups = (int**) malloc(Mmax*sizeof(int*));
    for(i=0; i<Mmax; i++) real_all_coups[i] = (int*) malloc(3*sizeof(int));


    k=0;
    for(i=0; i<N; i++) for(j=(i+1); j<N; j++)
    {
        for(e=0; e<R; e++)
        {
            all_coups[k+Mmax*e][1] = i;
            all_coups[k+Mmax*e][2] = j;
            all_coups[k+Mmax*e][0] = -1; ///por defecto esta desactivado
        }
        k++;
    }

    char **D = (char**) malloc(C*sizeof(char*)); ///aquí se guardan las configuraciones sampleadas
    for(i=0; i<C; i++) D[i] = (char*) malloc(N*sizeof(char));
    int **lambda = (int**) malloc(C*sizeof(int*)); /// el argumento de la tangente hiperbolica
    for(i=0; i<C; i++) lambda[i] = (int*) malloc(N*R*sizeof(int));
    double **FUNK = (double**) malloc(2*sizeof(double*)); ///log(...)
    for(i=0; i<2; i++) FUNK[i] = (double*) malloc((2*N+1)*sizeof(double));


    int *columna_i_new = (int*) malloc(C*sizeof(int));
    int *columna_j_new = (int*) malloc(C*sizeof(int));
    int **columna_i_new_wl = (int**) malloc(R*sizeof(int*));
    for(i=0; i<R; i++) columna_i_new_wl[i] = (int*) malloc(C*sizeof(int));
    int **columna_j_new_wl = (int**) malloc(R*sizeof(int*));
    for(i=0; i<R; i++) columna_j_new_wl[i] = (int*) malloc(C*sizeof(int));

    double *av_gamma = (double*) malloc(n_temps*sizeof(double));
    double *av_gamma2 = (double*) malloc(n_temps*sizeof(double));
    double *av_sensit = (double*) malloc(n_temps*sizeof(double));
    double *av_sensit2 = (double*) malloc(n_temps*sizeof(double));
    double *av_precis = (double*) malloc(n_temps*sizeof(double));
    double *av_precis2 = (double*) malloc(n_temps*sizeof(double));
    double *av_F1 = (double*) malloc(n_temps*sizeof(double));
    double *av_F1_2 = (double*) malloc(n_temps*sizeof(double));
    double *acoplo_2_list = (double*) malloc(tau*sizeof(double));
    double *acoplo_3_list = (double*) malloc(tau*sizeof(double));
    double *acoplo_4_list = (double*) malloc(tau*sizeof(double));

    delta_gamma = (gamma_ini-gamma_fin)/tau;
    tabulate_acoplos(acoplo_2_list, acoplo_3_list, acoplo_4_list, tau, gamma_ini, delta_gamma, teff, R);
    //for(i=0; i<tau; i++) printf("%lf\n", acoplo_2_list[i]);


    ///inicialicamos a cero el vector averaged_gamma(T):
    for(i=0; i<n_temps; i++) av_gamma[i] = av_gamma2[i] = 0.0;
    for(i=0; i<n_temps; i++) av_sensit[i] = av_sensit2[i] = 0.0;
    for(i=0; i<n_temps; i++) av_precis[i] = av_precis2[i] = 0.0;
    for(i=0; i<n_temps; i++) av_F1[i] = av_F1_2[i] = 0.0;


    double amigo;


    ///nombres de archivos;
    char namefile[] = "t000/D_000.dat";
    char namegfile[] = "GEOMETRY/G_000.dat";
    //printf("%d\n", n_seeds);




    for(i_seed = 0; i_seed<n_seeds; i_seed++)
    {
        create_name(namefile, (first_seed + i_seed));

        //printf("%s\n", namefile);

        ///leemos la geometria:
        change_g_name(namegfile, first_seed + i_seed);
        read_geometry(namegfile, &M, &zmax, vecinos, &JJ_0, N);
        //muestra_lista(vecinos, N);





        ///bucle en temperaturas:
        for(i_temps = 0; i_temps<n_temps; i_temps++)
        {
            ///leemos las configuraciones sampleadas (matriz D)
            update_name(namefile, integer_temps[i_temps]);
            printf("%s\n", namefile);
            read_D(namefile, i_seed, D, C, N);
            //comprobacion(D, C, N);
            //for(i=0; i<N; i++) printf("%d\n", D[1][i]);

            for(i=0; i<Mmax*R; i++) all_coups[i][0] = -1;///por defecto esta desactivado
            ///initial random guess:
            random_initial_couplings(all_coups, JJ, Mmax, &M_times_R_guess, 0.956, R);



            ///tabulamos lambda y FUNK
            tabulate_lambda(lambda, D, C, N, JJ, M_times_R_guess, R);
            //for(i=0; i<N; i++) printf(" %d\n", lambda[15][i]);
            temperatura = integer_temps[i_temps]*0.01;
            tabulate_FUNK(FUNK, N, temperatura);
            //for(i=0; i<(2*zmax+1); i++) printf("%lf %lf\n", FUNK[0][i], FUNK[1][i]);


            ///SIMULACION DE MONTECARLO:

            gamma = gamma_ini ;
            PERESEUDO = 0.0;
            PERESEUDO = compute_PLH(FUNK, N, D, lambda, C, N, R);



            for(i=0; i<tau; i++)
            {

                //double pseudo;

                if(i%n_wl==0) WORLDLINE_FLIP(JJ, Mmax, columna_i_new_wl, columna_j_new_wl, D, lambda, C, N, FUNK, N, teff, R,acoplo_2_list[i], acoplo_3_list[i], acoplo_4_list[i], PLH_strength, all_coups, &M_times_R_guess, M);
                else MC_STEP(JJ, Mmax, columna_i_new, columna_j_new, D, lambda, C, N, FUNK, N, teff, R, acoplo_2_list[i], acoplo_3_list[i], acoplo_4_list[i], PLH_strength, all_coups, &M_times_R_guess, M);
                gamma -= delta_gamma;
                //printf("gamma = %lf\n", gamma);
                //teff-=delta_teff;
                //printf("%lf %lf %lf\n", acoplo_2_list[i], acoplo_3_list[i], acoplo_4_list[i]);
                //print_couplings_mean(JJ, M);
                //pseudo = compute_PLH(FUNK, N, D, lambda, C, N, R);
                //printf("pseudo = %lf\n", pseudo);
                //printf("PSEUDO = %lf\n", PERESEUDO);

            }




            make_the_guess(JJ, all_coups, Mmax, R);
            //print_couplings_mean(JJ, M);
            //save_trotter_copies(JJ, M_times_R_guess, "trotter_copies.txt");

            aux_gamma = GAMMA(real_all_coups, all_coups, Mmax, JJ_0, M, N);
            av_gamma[i_temps] += aux_gamma;
            av_gamma2[i_temps] += aux_gamma*aux_gamma;
            aux_gamma = F1_score(real_all_coups, Mmax, all_coups, &aux_success, &aux_precis);
            av_sensit[i_temps] += aux_success;
            av_sensit2[i_temps] += aux_success*aux_success;
            av_precis[i_temps] += aux_precis;
            av_precis2[i_temps] += aux_precis*aux_precis;
            av_F1[i_temps] += aux_gamma;
            av_F1_2[i_temps] += aux_gamma*aux_gamma;


        }

        ///liberamos los vectores que se crean de iteracion a iteracion:
        for(i=0; i<N; i++) free(vecinos[i]);
        for(i=0; i<M; i++) free(JJ_0[i]);
        free(JJ_0);


    }


    save_results(integer_temps, av_gamma, av_gamma2, av_sensit, av_sensit2, av_precis, av_precis2, av_F1, av_F1_2, n_temps, n_seeds , first_seed);
    //for(i=0; i<n_temps; i++) av_gamma[i]/=n_seeds;
    //for(i=0; i<n_temps; i++) printf("%lf\n", av_gamma[i]);


    for(i=0; i<C; i++)
    {
        free(D[i]);
        free(lambda[i]);
    }
    free(D);
    free(lambda);
    for(i=0; i<2; i++) free(FUNK[i]);
    free(FUNK);
    free(columna_i_new);
    free(columna_j_new);
    for(i=0;i<R; i++)
    {
        free(columna_i_new_wl[i]);
        free(columna_j_new_wl[i]);
    }
    free(columna_i_new_wl);
    free(columna_j_new_wl);

    free(integer_temps);
    free(av_gamma);
    free(av_gamma2);
    free(vecinos);
    free(acoplo_2_list);
    free(acoplo_3_list);
    free(acoplo_4_list);

    for(i=0; i<Mmax*R; i++)
    {
        free(JJ[i]);
        free(all_coups[i]);
    }
    for(i=0; i<Mmax; i++) free(real_all_coups[i]);
    free(JJ);
    free(all_coups);
    free(real_all_coups);
    free(av_sensit);
    free(av_sensit2);
    free(av_precis);
    free(av_precis2);
    free(av_F1);
    free(av_F1_2);



    return 0;
}


void ini_ran(int seed)
{
    int INI, FACTOR, SUM, i;
    INI = seed;
    FACTOR = 67397;
    SUM = 7364893;

    for(i=0; i<256; i++)
    {
        INI = (INI*FACTOR + SUM);
        larueda[i] = INI;
    }
    ind_ran = i1 = i2 = i3 = 0;
}

double Rand(void)
{
    double r;
    i1 = ind_ran-24;
    i2 = ind_ran-55;
    i3 = ind_ran-61;
    larueda[ind_ran] = larueda[i1] + larueda[i2];
    number = (larueda[ind_ran]^larueda[i3]);
    ind_ran++;
    r = number*NormRANu;
    return r;
}

void crea_lista_vecinos(int** vecinos, int L, int N)
{
    int x, y, index, aux;
    for(x=0; x<L; x++) for(y=0; y<L; y++)
    {
        index = x + L*y;
        vecinos[index][0] = index + L; ///vecino de arriba
        vecinos[index][1] = index + 1; ///vecino de la derecha
        vecinos[index][2] = index - L; ///vecino de abajo
        vecinos[index][3] = index -1; ///vecino de la izquierda
    }

    for(index=0; index<N; index++)
    {
        //printf("%d\n",vecinos[index][2]);
        //printf("%lf\n", floor(((double)vecinos[index][2])/N)*N);
        vecinos[index][0] -= floor(((double)vecinos[index][0])/N)*N;
        vecinos[index][2] = vecinos[index][2] - floor(((double)vecinos[index][2])/N)*N;
        //printf("%d\n",vecinos[index][2]);


        aux = vecinos[index][1]-floor((double)index/L)*L;
        aux-= floor((double)aux/L)*L;
        vecinos[index][1] = floor((double)index/L)*L + aux;

        aux = vecinos[index][3]-floor((double)index/L)*L;
        aux-= floor((double)aux/L)*L;
        vecinos[index][3] = floor((double)index/L)*L + aux;

    }
}

void muestra_lista(int **vecinos, int N)
{
    int i,j;
    for(i=0; i<N; i++)
    {
        printf("%d   ", i);
        for(j=0; j<z; j++) printf("%d ", vecinos[i][j]);
        printf("\n");
    }
}

int sign(int n)
{
    if(n>0) return 1;
    if(n<0) return -1;
    if(n==0) return 0;
}



void ini_r1279()
{
    int i,j, one_bit;

    iaux = 0;
    for(i=0; i<2048; i++)
    {
        index1[i] = (i-1279)&2047;
        index2[i] = (i-418)&2047;

        secuencia_rand[i] = 0;
        for(j=0; j<=NBITM1; j++)
        {
            one_bit = 0;
            if(Rand()>0.5) one_bit = 1;

            one_bit = one_bit << j; ///we shift j places to the left
            secuencia_rand[i] = secuencia_rand[i]|one_bit; ///hacemos un or bit a bit
        }
        secuencia_rand[i] = 2*secuencia_rand[i] + 1;
    }
}

double r1279(void)
{
    int numerin;
    iaux = (iaux+1)&2047;
    secuencia_rand[iaux] = secuencia_rand[index1[iaux]]*secuencia_rand[index2[iaux]];
    numerin = (secuencia_rand[iaux]>>1)&2147483647;
    return ((double)numerin)*normNorm;
}

void read_D(char *namefile, int i_seed, char **D, int C, int N)
{
    int mu, i;
    FILE *f;

    f = fopen(namefile, "rb");
    if(f==NULL) exit(i_seed);

    for(mu=0; mu<C; mu++)
    {
        fread(D[mu], sizeof(char), N, f);
    }
    fclose(f);
}

void comprobacion(char **D, int C, int N)
{
    int mu, i;
    double suma;

    for(mu=0; mu<C; mu++)
    {
        suma = 0.0;
        for(i=0; i<N; i++) suma += D[mu][i]*1.0;
        printf("%d %lf\n", mu, suma/N);
    }
}

void create_couplings_list(int**vecinos, int N, int zmax, int**JJ, int M)
{
    int i, j, k;

    k=0;
    for(i=0; i<N; i++)
    {
        j=0;
        while((j<zmax)&&(vecinos[i][j]!=i))
        {
            if(abs(vecinos[i][j])>i) ///we only consider ordered pairs
            {
                JJ[k][1] = i;
                JJ[k][2] = abs(vecinos[i][j]);
                if(vecinos[i][j]>0) JJ[k][0] = 1;
                else JJ[k][0] = -1;
                //printf("El coupling %d vale %d\n", k, JJ[k][0]);
                k++;

            }
            j++;
        }
    }
    //printf("k = %d and M = %d\n", k, M);

}

void tabulate_lambda (int **lambda, char **D, int C, int N, int **JJ, int M_times_R_guess, int R)
{
    int mu, i, k, tau, twoN_3, a, b;
    twoN_3 = 2*N-3;

    for(mu=0; mu<C; mu++) for(i=0; i<N*R; i++) lambda[mu][i] = 0;

    for(i=0; i<M_times_R_guess; i++)
    {
        a = JJ[i][1];
        b = JJ[i][2];
        tau = abs(JJ[i][0])-1;
        if (tau<0) printf("tau = %d\n", tau);
        for(mu=0; mu<C; mu++)
        {
            lambda[mu][tau*N+JJ[i][1]] += sign(JJ[i][0])*D[mu][JJ[i][2]];
            lambda[mu][tau*N+JJ[i][2]] += sign(JJ[i][0])*D[mu][JJ[i][1]];
        }
    }

}

void tabulate_FUNK (double **FUNK, int zmax, double temp)
{
    int i, j , aux_spin;
    double eps;
    for(i=0; i<2; i++)
    {
        aux_spin = -1 +2*i;
        for(j=-zmax; j<=zmax; j++)
        {
            eps = 0.5*(1.0 + aux_spin*tanh(((double)j)/temp));
            if(eps<1.0e-100) eps = 1.0e-100;
            FUNK[i][j+zmax] = log(eps);
        }
    }
}


double compute_PLH( double **FUNK, int zmax, char **D, int **lambda, int C, int N, int R)
{
    double suma, PLH = 0.0;
    int mu, i, aux_ind, tau;

    for(tau=0; tau<R; tau++) for(i=0; i<N; i++)
    {
        suma = 0.0;
        for(mu=0; mu<C; mu++)
        {
            aux_ind = (1+D[mu][i])/2;
            suma += FUNK[aux_ind][lambda[mu][tau*N+i]+zmax];
        }
        PLH += suma/C;
    }
    return PLH;
}

void MC_STEP(int **JJ, int Mmax, int *columna_i_new, int *columna_j_new, char **D, int **lambda, int C, int N, double **FUNK, int zmax, double teff, int R, double acoplo_2, double acoplo_3, double acoplo_4, double PLH_strength, int **all_coups, int *M_times_R_guess_dir, int M)
{
    char flag = 0;
    int e, k, i, j, mu, coupling, i_1, j_1, tau, thechosen, new_coupling, delta_coup, delta_coup_2;
    int prev_tau, next_tau, prev_coup, next_coup;
    double DELTA_PLH, omega, DELTA_K2, DELTA_K3, DELTA_K4;
    double suma, pen_term;

    for(e=0; e<Mmax*R; e++)
    {
        do
        {
            k = Mmax*Rand();
        } while(k==Mmax);
        do
        {
            tau = R*Rand();
        } while(tau==R);
        thechosen = tau*Mmax + k;
        prev_tau = tau-1;
        if(prev_tau<0) prev_tau = R-1;
        next_tau = tau+1;
        if(next_tau==R) next_tau=0;
        if(all_coups[prev_tau*Mmax + k][0]<0) prev_coup = 0;
        else prev_coup = sign(JJ[all_coups[prev_tau*Mmax + k][0]][0]);
        if(all_coups[next_tau*Mmax + k][0]<0) next_coup = 0;
        else next_coup = sign(JJ[all_coups[next_tau*Mmax + k][0]][0]);


        i = all_coups[thechosen][1];
        j = all_coups[thechosen][2];
        if(all_coups[thechosen][0]<0) coupling = 0;
        else coupling = sign(JJ[all_coups[thechosen][0]][0]);
        //printf("Coupling = %d\", coupling);

        ///hasta aquí lo que hemos hecho es seleccionar un coupling al azar y ver cuanto valen sus copias de Trotter contiguas


        ///PROPONEMOS UN CAMBIO:
        if(Rand()>0.5) new_coupling = coupling +1;
        else new_coupling = coupling -1;
        if(new_coupling == -2) new_coupling = 1;
        if(new_coupling == 2) new_coupling = -1;
        delta_coup = new_coupling - coupling;
        delta_coup_2 = new_coupling*new_coupling -coupling*coupling;

        for(mu=0; mu<C; mu++)
        {
            columna_i_new[mu] = lambda[mu][tau*N + i] +delta_coup*D[mu][j];
            columna_j_new[mu] = lambda[mu][tau*N + j] +delta_coup*D[mu][i];
        } ///(solo cambian dos columnas de la matriz lambda al cambiar el coupling J_ij

        /*printf("Seleccionamos la trotter copy %d que vale %d\n", thechosen, coupling);
        printf("que es entre los espines %d y %d\n", i, j);
        printf("tau = %d, k = %d, prev_tau = %d; next_tau = %d\n", tau, k, prev_tau, next_tau);
        printf("prev_coup = %d, next_coup = %d, new_coup = %d\n", prev_coup, next_coup, new_coupling);
        printf("lambda_i = %d, lambda_j = %d\n", lambda[15][tau*N +i], lambda[15][tau*N+j]);
        printf("columna_new_i = %d, columna_new_j = %d\n", columna_i_new[15], columna_j_new[15]);*/


        ///CÁLCULO DE DELTA_PLH:
        DELTA_PLH = 0.0;
        for(mu=0; mu<C; mu++)
        {

            i_1 = (D[mu][i]+1)/2;
            j_1 = (D[mu][j]+1)/2;

            DELTA_PLH += (FUNK[i_1][columna_i_new[mu]+zmax] + FUNK[j_1][columna_j_new[mu]+zmax]);
            DELTA_PLH -= (FUNK[i_1][lambda[mu][tau*N + i]+zmax] + FUNK[j_1][lambda[mu][tau*N + j]+zmax]);
        }
        AMIGOTE = DELTA_PLH/C;
        DELTA_PLH/=(C*R);
        DELTA_PLH*=PLH_strength;
        /*if(coupling*new_coupling!=0)*/// printf("DELTA_PLH vale %lf\n", DELTA_PLH);
        DELTA_K2 = acoplo_2*delta_coup*(prev_coup + next_coup);
        DELTA_K3 = acoplo_3*delta_coup_2;
        DELTA_K4 = acoplo_4*delta_coup_2*(prev_coup*prev_coup + next_coup*next_coup);
        //if(coupling*new_coupling!=0) printf("acoplo2 = %lf, DELTA_K2 = %lf\n", acoplo_2, DELTA_K2);
        DELTA_PLH -= (DELTA_K2+DELTA_K3+DELTA_K4);
        //if(coupling*new_coupling!=0)*/ printf("DELTA_PLH final vale %lf\n", DELTA_PLH);

        ///el termino de penalty afecta a DELTA_PLH:
        if(new_coupling*coupling==0) ///si se ha creado o destruido algun enlace
        {
            if((*M_times_R_guess_dir)<M*R) pen_term = (abs(new_coupling)-abs(coupling))*bond_penalty; ///aquí hacemos que los procesos de creacion aumenten PLH
            if((*M_times_R_guess_dir)>M*R) pen_term = -(abs(new_coupling)-abs(coupling))*bond_penalty; ///aquí hacemos que los procesos de destruccion disminuyan PLH
            if((*M_times_R_guess_dir)==M*R) pen_term = -bond_penalty; ///si  ya estas en el target, tanto crear como destruir enlaces se penaliza
            DELTA_PLH += pen_term/R;
        }

        ///VEMOS SI ACEPTAMOS O NO (ALGORITMO DE METROPOLIS)
        flag = 0;
        if(DELTA_PLH>0.0) flag = 1;
        else
        {
            omega = Rand();
            //if(coupling*new_coupling!=0)printf("El factor Boltzmann vale %lf y omega vale %lf\n",exp(DELTA_PLH/teff), omega);
            if(omega<exp(DELTA_PLH/teff)) flag = 1;

        }

        //if(coupling*new_coupling==0) flag = 0;


        ///LO QUE SUCEDE SI ACEPTAMOS:
        if(flag == 1)
        {
            PERESEUDO += AMIGOTE;
            //printf("Aceptamos el cambio\n");
            ///hay 3 casos: creacion, destruccion y flip

            ///CREACION:
            if(coupling==0) append_new_coupling(all_coups, JJ, M_times_R_guess_dir, thechosen, new_coupling, tau);

            ///DESTRUCCION
            if(new_coupling==0) delete_coupling(all_coups, JJ, M_times_R_guess_dir, thechosen, N, Mmax);

            ///FLIP
            if(new_coupling*coupling!=0) JJ[all_coups[thechosen][0]][0] = -JJ[all_coups[thechosen][0]][0];


            /*suma = lambda[15][tau*N+i];
            tabulate_lambda(lambda, D, C, N, JJ, *M_times_R_guess_dir, R);
            if(lambda[15][tau*N+i]!=columna_i_new[15])
            {
                printf("Vamos a ver que ha pasado aqui. Teniamos un coupling, la copia %d del bond %d, que valia %d. Era entre %d y %d\n", tau, k, coupling, i, j);
                printf("Hemos cambiado el coupling a %d. Lambda(15; tau = %d, i = %d) valia %d. Computamos cuanto valdra ahora y nos sale que %d\n", new_coupling, tau, i, suma, columna_i_new[15]);
                printf("Esto es porque el espin (15,%d) vale %d\n", j, D[15][j]);
                printf("No obstante, al aplicar los cambios y volver a tabular lambda, obtenemos que Lambda(15; tau = %d, i = %d) vale %d\n", tau, i, lambda[15][tau*N+i]);
                printf("%d vs %d\n", lambda[15][tau*N+i], columna_i_new[15]);

            }*/


            for(mu=0; mu < C; mu++)
            {
                lambda[mu][tau*N+i] = columna_i_new[mu];
                lambda[mu][tau*N+j] = columna_j_new[mu];
            }



            //printf("Delta PLH real = %lf\n", suma - compute_PLH(FUNK, zmax, D, lambda, C, N));
            //printf("Supuesto delta PLH = %lf\n", DELTA_PLH);
            //tabulate_lambda(lambda, D, C, N, JJ, *M_times_R_guess_dir, R);
            //printf("lambda_i = %d, lambda_j = %d\n", lambda[15][tau*N +i], lambda[15][tau*N+j]);

        }
        //else if(coupling*new_coupling!=0) printf("No aceptamos el cambio\n");
        //printf("\n");

    }
    //printf("number of nonzero Trotter copies = %d\n", *M_times_R_guess_dir);

}

void print_couplings_mean(int **JJ, int M)
{
    double suma = 0.0;
    int k;
    for(k=0; k<M; k++) suma += (double)JJ[k][0];
    printf("El promedio de los couplings es %lf\n", suma/M);

}


void create_name(char *namefile, int i_seed)
{
    int i, dummy;
    dummy = i_seed;
    for(i=0; i<3; i++)
    {
        namefile[9-i] = dummy%10;
        dummy = (dummy-namefile[9-i])/10;
        namefile[9-i]+=48;
    }
    //printf("%s\n", namefile);
}

void gamma_stats(double *gamma_vector, int n_seeds)
{
    double suma, suma2;
    suma = suma2 = 0.0;
    int i;

    for(i=0; i<n_seeds; i++)
    {
        suma+=gamma_vector[i];
        suma2+=gamma_vector[i]*gamma_vector[i];
    }
    suma/=n_seeds;
    suma2/=n_seeds;
    suma2 = suma2 - suma*suma;
    suma2*=(((double)n_seeds)/((double)(n_seeds-1)));
    printf("gamma = %lf pm %lf", suma, sqrt(suma/n_seeds));
}

void lee_input(char *namefile, int *N_dir, int *C_dir, int *n_temps_dir, int **temps_dir, int *n_seeds_dir, double *gamma_ini_dir, double*gamma_fin_dir, int *tau_dir, int *first_seed_dir, int *R_dir, int *n_wl_dir, double *teff_dir, double *PLH_strength_dir)
{
    int i;
    FILE *f;

    f = fopen(namefile, "rt");
    if(f==NULL) exit(1);

    fscanf(f, "%d\n%d\n%d\n", N_dir, C_dir, n_temps_dir);
    //printf("%d\n", *n_temps_dir);

    (*temps_dir) = (int*)malloc((*n_temps_dir)*sizeof(int));
    for(i=0; i<(*n_temps_dir); i++)
    {
        fscanf(f, "%d ", &(*temps_dir)[i]);
    }

    fscanf(f, "\n%d\n%lf\n%lf\n%d\n%d\n%d\n%d\n%lf\n%lf\n%lf", n_seeds_dir, gamma_ini_dir, gamma_fin_dir, tau_dir, first_seed_dir, R_dir, n_wl_dir, teff_dir, PLH_strength_dir, &bond_penalty);
    fclose(f);
}

void change_g_name(char *name_g_file, int i_seed)
{
    int dummy, i;
    dummy = i_seed;
    for(i=0; i<3; i++)
    {
        name_g_file[13-i] = dummy%10;
        dummy =  (dummy-name_g_file[13-i])/10;
        name_g_file[13-i] += 48; ///convertimos a codigo ASCII
    }
}

void read_geometry(char *namegfile, int *M_dir, int *zmax_dir, int **vecinos, int ***JJ_0_dir, int N)
{
    int i,j;
    FILE *f;

    f = fopen(namegfile, "rb");
    if(f==NULL) exit(3);

    fread(M_dir, sizeof(int), 1, f);
    fread(zmax_dir, sizeof(int), 1, f);

    for(i=0; i<N; i++) vecinos[i] = (int*) malloc((*zmax_dir)*sizeof(int));
    for(i=0; i<N; i++) fread(vecinos[i], sizeof(int), (*zmax_dir), f);
    fclose(f);

    (*JJ_0_dir) = (int**) malloc((*M_dir)*sizeof(int*)); /// couplings reales (los que queremos inferir)
    for(i=0; i<(*M_dir); i++) (*JJ_0_dir)[i] = (int*) malloc(3*sizeof(int));

    create_couplings_list(vecinos, N, *zmax_dir, *JJ_0_dir, *M_dir); ///aquí le estamos pasando la matriz de vecinos "signada" (los vecinos aparecen con un - delante si la relacion es antiferromagnetica)
    ///des-signamos la matriz de vecinos:
    for(i=0; i<N; i++) for(j=0; j<(*zmax_dir); j++) vecinos[i][j] = abs(vecinos[i][j]);

}

void update_name(char *namefile, int temp)
{
    int i, dummy;

    dummy = temp;
    for(i=0; i<3; i++)
    {
        namefile[3-i] = dummy%10;
        dummy =  (dummy - namefile[3-i])/10;
        namefile[3-i] += 48; ///pasamos a codigo ASCII
    }

}

void save_results(int *temps, double *av_gamma, double *av_gamma2, double *av_sensit, double *av_sensit2, double *av_precis, double *av_precis2, double *av_F1, double *av_F1_2, int n_temps, int n_seeds, int first_seed)
{
    int i, final_seed, dummy;
    char namefile[] = "av_gamma_vs_T_000_to_000_bp_0000.txt";

    final_seed = first_seed + n_seeds -1;

    dummy = first_seed;
    for(i=0; i<3; i++)
    {
        namefile[16-i] = dummy%10;
        dummy = (dummy - namefile[16-i])/10;
        namefile[16-i] += 48;
    }

    dummy = final_seed;
    for(i=0; i<3; i++)
    {
        namefile[23-i] = dummy%10;
        dummy = (dummy - namefile[23-i])/10;
        namefile[23-i] += 48;
    }

    dummy = (bond_penalty+0.0003)*1000.0;
    for(i=0; i<4; i++)
    {
        namefile[31-i] = dummy%10;
        dummy = (dummy - namefile[31-i])/10;
        namefile[31-i] += 48;
    }

    FILE *f = fopen(namefile, "wt");
    if(f==NULL) exit(first_seed);
    for(i=0; i<n_temps; i++) fprintf(f,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %d\n", temps[i]*0.01, av_gamma[i], av_gamma2[i], av_sensit[i], av_sensit2[i], av_precis[i], av_precis2[i], av_F1[i], av_F1_2[i], n_seeds);
    fclose(f);
    printf("%s", namefile);
}

void tabulate_acoplos(double *acoplo_2_list, double *acoplo_3_list, double *acoplo_4_list, int tau, double gamma_ini, double delta_gamma, double temp, int R)
{
    int i;
    double gamma;
    gamma = gamma_ini;
    for(i=0; i<tau; i++)
    {
        acoplo_2_list[i] = temp*log(tanh(0.5*gamma/temp/R));
        acoplo_3_list[i] = -2.0*temp*log(tanh(gamma/temp/R)/sqrt(2.0));
        acoplo_4_list[i] = temp*log(tanh(gamma/temp/R));
        gamma -= delta_gamma;
    }
}

void make_the_guess(int **JJ, int **all_coups, int Mmax, int R)
{
    int k, j, tau, n_zero, n_pos, n_neg, coupling;

    for(k=0; k<Mmax; k++)
    {
        n_zero = n_pos = n_neg = 0;
        //printf("Estudiamos las guesses del coupling %d\n", k);
        for(tau=0; tau<R; tau++)
        {
            if(all_coups[k+tau*Mmax][0]<0) n_zero++;
            else
            {
                coupling = JJ[all_coups[k+tau*Mmax][0]][0];
                if(coupling>0) n_pos++;
                else n_neg++;
            }

        }

        if(n_pos>n_zero)
        {
            if(n_pos>n_neg) all_coups[k][0] = 1;
            if(n_pos==n_neg) all_coups[k][0] = 0;
            if(n_pos<n_neg) all_coups[k][0] = -1;
        }
        if(n_pos<=n_zero)
        {
            if(n_neg>n_zero) all_coups[k][0] = -1;
            else all_coups[k][0] = 0;
        }

    }
}

double compute_success(int **JJ_0, int **JJ, int M)
{
    int k;
    int n_aciertos = 0;
    double success;
    for(k=0; k<M; k++) if(JJ_0[k][0]==JJ[k][0]) n_aciertos++;
    success = (n_aciertos*1.0)/M;
    return success;
}

void append_new_coupling(int **all_coups, int **JJ, int *M_times_R_guess_dir, int thechosen, int coupling, int tau)
{
    int our_M_x_R = *M_times_R_guess_dir;
    JJ[our_M_x_R][0] = coupling*(tau+1);
    JJ[our_M_x_R][1] = all_coups[thechosen][1];
    JJ[our_M_x_R][2] = all_coups[thechosen][2];
    all_coups[thechosen][0] = our_M_x_R;
    our_M_x_R++;
    *M_times_R_guess_dir = our_M_x_R;
}

void delete_coupling(int **all_coups, int **JJ, int *M_guess_dir, int thechosen, int N, int Mmax)
{
    int i, j, a, b, k, position, tau, our_M, twoN_3;
    our_M = *M_guess_dir;
    position = all_coups[thechosen][0];
    twoN_3 = 2*N-3;



    for(i=position; i<our_M; i++) for(j=0; j<3; j++) JJ[i][j] = JJ[i+1][j];
    all_coups[thechosen][0] = -1; ///coupling desactivado
    our_M--;
    ///actualizamos las posiciones de todos los demas:
    for(i=position; i<our_M; i++)
    {
        a = JJ[i][1];
        b = JJ[i][2];
        k = (-a*a +twoN_3*a -2 +2*b)/2;
        tau = abs(JJ[i][0])-1;
        all_coups[tau*Mmax + k][0]--;
    }

    *M_guess_dir = our_M;
}

void WORLDLINE_FLIP(int **JJ, int Mmax, int **columna_i_new_wl, int **columna_j_new_wl, char **D, int **lambda, int C, int N, double **FUNK, int zmax, double teff, int R, double acoplo_2, double acoplo_3, double acoplo_4, double PLH_strength, int **all_coups, int *M_times_R_guess_dir, int M)
{
    char flag = 0;
    int e, k, i, j, mu, coupling, i_1, j_1, tau, thechosen, new_coupling, delta_coup, delta_coup_2;
    int prev_tau, next_tau, prev_coup, next_coup, elevator;
    double DELTA_PLH, omega, DELTA_K2, DELTA_K3, DELTA_K4;
    double suma, pen_term;

    for(e=0; e<Mmax; e++)
    {
        do
        {
            k = Mmax*r1279();
        } while(k==Mmax);


        DELTA_PLH = DELTA_K2 = DELTA_K3 = DELTA_K4 = 0.0;
        i = all_coups[k][1];
        j = all_coups[k][2];

        ///PROPONEMOS UN CAMBIO:
        if(Rand()>0.5) elevator = 1;
        else elevator = -1;

        for(tau=0; tau<R; tau++)
        {
            thechosen = tau*Mmax + k;
            prev_tau = tau-1;
            if(prev_tau<0) prev_tau = R-1;
            next_tau = tau+1;
            if(next_tau==R) next_tau=0;
            if(all_coups[prev_tau*Mmax + k][0]<0) prev_coup = 0;
            else prev_coup = sign(JJ[all_coups[prev_tau*Mmax + k][0]][0]);
            if(all_coups[next_tau*Mmax + k][0]<0) next_coup = 0;
            else next_coup = sign(JJ[all_coups[next_tau*Mmax + k][0]][0]);

            if(all_coups[thechosen][0]<0) coupling = 0;
            else coupling = sign(JJ[all_coups[thechosen][0]][0]);

            ///hasta aquí lo que hemos hecho es seleccionar una replica del coupling seleccionado al azar
            //printf("Seleccionamos el coupling numero %d que vale %d\n", thechosen, JJ[thechosen][0]);
            //printf("que es entre los espines %d y %d\n", i, j);

            new_coupling = coupling + elevator; ///el mismo cambio para todas las Trotter copies
            if(new_coupling == -2) new_coupling = 1;
            if(new_coupling == 2) new_coupling = -1;
            delta_coup = new_coupling - coupling;
            delta_coup_2 = new_coupling*new_coupling -coupling*coupling;

            for(mu=0; mu<C; mu++)
            {
                columna_i_new_wl[tau][mu] = lambda[mu][tau*N + i] +delta_coup*D[mu][j];
                columna_j_new_wl[tau][mu] = lambda[mu][tau*N + j] +delta_coup*D[mu][i];
            } ///(solo cambian dos columnas de la matriz lambda al cambiar el coupling J_ij

            ///CÁLCULO DE DELTA_PLH:

            for(mu=0; mu<C; mu++)
            {

                i_1 = (D[mu][i]+1)/2;
                j_1 = (D[mu][j]+1)/2;

                DELTA_PLH += (FUNK[i_1][columna_i_new_wl[tau][mu]+zmax] + FUNK[j_1][columna_j_new_wl[tau][mu]+zmax]);
                DELTA_PLH -= (FUNK[i_1][lambda[mu][tau*N + i]+zmax] + FUNK[j_1][lambda[mu][tau*N + j]+zmax]);
            }

            DELTA_K2 += acoplo_2*delta_coup*(prev_coup + next_coup);
            DELTA_K3 += acoplo_3*delta_coup_2;
            DELTA_K4 += acoplo_4*delta_coup_2*(prev_coup*prev_coup + next_coup*next_coup);

        }
        AMIGOTE = DELTA_PLH/C;
        DELTA_PLH/=(C*R);
        DELTA_PLH*=PLH_strength;
        DELTA_PLH -= (DELTA_K2+DELTA_K3+DELTA_K4);

        ///el termino de penalty afecta a DELTA_PLH:
        if(new_coupling*coupling==0) ///si se ha creado o destruido algun enlace
        {
            if((*M_times_R_guess_dir)<M*R) pen_term = (abs(new_coupling)-abs(coupling))*bond_penalty; ///aquí hacemos que los procesos de creacion aumenten PLH
            if((*M_times_R_guess_dir)>M*R) pen_term = -(abs(new_coupling)-abs(coupling))*bond_penalty; ///aquí hacemos que los procesos de destruccion disminuyan PLH
            if((*M_times_R_guess_dir)==M*R) pen_term = -bond_penalty; ///si  ya estas en el target, tanto crear como destruir enlaces se penaliza
            DELTA_PLH += pen_term/R;
        }


        ///VEMOS SI ACEPTAMOS O NO (ALGORITMO DE METROPOLIS)
        flag = 0;
        if(DELTA_PLH>0.0) flag = 1;
        else
        {
            omega = r1279();
            //printf("El factor Boltzmann vale %lf y omega vale %lf\n",exp(DELTA_PLH/teff), omega);
            if(omega<exp(DELTA_PLH/teff)) flag = 1;

        }


        ///LO QUE SUCEDE SI ACEPTAMOS:
        if(flag == 1)
        {
            //suma = compute_PLH(FUNK, zmax,  D, lambda, C, N);
            //printf("Aceptamos el cambio %lf\n");
            PERESEUDO+=AMIGOTE;
            for(tau=0; tau<R; tau++)
            {
                thechosen = tau*Mmax + k;
                if(all_coups[thechosen][0]<0) coupling = 0;
                else coupling = sign(JJ[all_coups[thechosen][0]][0]);
                new_coupling = coupling + elevator; ///el mismo cambio para todas las Trotter copies
                if(new_coupling == -2) new_coupling = 1;
                if(new_coupling == 2) new_coupling = -1;

                ///hay 3 casos: creacion, destruccion y flip

                ///CREACION:
                if(coupling==0) append_new_coupling(all_coups, JJ, M_times_R_guess_dir, thechosen, new_coupling, tau);

                ///DESTRUCCION
                if(new_coupling==0) delete_coupling(all_coups, JJ, M_times_R_guess_dir, thechosen, N, Mmax);

                ///FLIP
                if(new_coupling*coupling!=0) JJ[all_coups[thechosen][0]][0] = -JJ[all_coups[thechosen][0]][0];

                JJ[thechosen][0] = -JJ[thechosen][0]; ///cambiamos el coupling
                for(mu=0; mu<C; mu++)
                {
                    lambda[mu][tau*N+i] = columna_i_new_wl[tau][mu];
                    lambda[mu][tau*N+j] = columna_j_new_wl[tau][mu];
                }

            }

            //printf("Delta PLH real = %lf\n", suma - compute_PLH(FUNK, zmax, D, lambda, C, N));
            //printf("Supuesto delta PLH = %lf\n", DELTA_PLH);
        }
        //else printf("No aceptamos el cambio\n");


    }
}

void random_initial_couplings(int **all_coups, int **JJ, int Mmax, int *Mguess_dir, double alpha, int R)
{
    ///alpha es la probabilidad de que un coupling valga 0. Se recomienda fijar alta, porque normalmente M/Mmax es pequeño

    int i, j, our_M, tau, k;

    for(i=0; i<Mmax*R; i++) for(j=0; j<3; j++) JJ[i][j] = 0;
    our_M = 0;

    for(tau=0; tau<R; tau++) for(k=0; k<Mmax; k++)
    {
        i = tau*Mmax + k;
        if(Rand()<alpha) all_coups[i][0] = -1; ///está desactivado
        else
        {
            all_coups[i][0] = our_M; ///posicion del coupling en la lista de couplings activos
            JJ[our_M][1] = all_coups[i][1];
            JJ[our_M][2] = all_coups[i][2];
            if(Rand()>0.5) JJ[our_M][0] = tau+1;
            else JJ[our_M][0] = -(tau+1);
            our_M++;
        }
    }
    *Mguess_dir = our_M;
}

double GAMMA(int **real_all_coups, int **all_coups, int Mmax, int **JJ_0, int M, int N)
{
    int i, j, k, index, twoN_3;
    double suma, diff;

    twoN_3 = 2*N-3;

    ///Primero, construimos "real_all_coups":
    k=0;
    for(i=0; i<N; i++) for(j=i+1; j<N; j++)
    {
        real_all_coups[k][1] = i;
        real_all_coups[k][2] = j;
        real_all_coups[k][0] = 0; ///desactivado por defecto
        k++;
    }


    for(k=0; k<M; k++)
    {
        i = JJ_0[k][1];
        j = JJ_0[k][2];
        index = (-i*i +twoN_3*i -2 +2*j)/2;
        real_all_coups[index][0] = JJ_0[k][0];
    }


    ///Ahora, calculamos GAMMA:

    suma = 0.00;
    for(k=0; k<Mmax; k++)
    {
        diff = (real_all_coups[k][0] - all_coups[k][0])*1.0;
        suma+=diff*diff;
    }

    suma/=(M);

    return suma;
}

double F1_score(int **real_all_coups, int Mmax, int **all_coups, double *sens_dir, double *precis_dir)
{
    int real_coup, guess_coup;
    int k;
    int n_correct_nonzero_predictions, n_nonzero_bonds, n_nonzero_predictions;


    n_nonzero_bonds = n_nonzero_predictions = n_correct_nonzero_predictions = 0;
    for(k=0; k<Mmax; k++)
    {
        real_coup = real_all_coups[k][0];
        guess_coup = all_coups[k][0];

        if(real_coup!=0) n_nonzero_bonds++;
        if(guess_coup!=0) n_nonzero_predictions++;
        if((real_coup==guess_coup)&&(real_coup!=0)) n_correct_nonzero_predictions++;
        /*if((real_coup!=0)&&(guess_coup==0))
        {
            printf("Problema con el link %d\n", k);
            printf("real_coup = %d, guess = %d\n", real_coup, guess_coup);
        }
        if((real_coup==0)&&(guess_coup!=0))
        {
            printf("Problema con el link %d\n", k);
            printf("real_coup = %d, guess = %d\n", real_coup, guess_coup);
        }*/


    }
    printf("n correct nonzero predictions = %d, n_nonzero_bonds = %d, n_nonzero_predictions = %d\n", n_correct_nonzero_predictions, n_nonzero_bonds, n_nonzero_predictions);

    *sens_dir = (n_correct_nonzero_predictions*1.00)/n_nonzero_bonds;
    *precis_dir = (n_correct_nonzero_predictions*1.00)/n_nonzero_predictions;
    //printf("n_nonzero_bonds = %d, n_nonzero_predictions = %d, n_aciertos = %d, s = %lf\n", n_nonzero_bonds, n_nonzero_predictions, n_correct_nonzero_predictions, *sens_dir);
    return 2.0*(*sens_dir)*(*precis_dir)/((*sens_dir)+(*precis_dir));
}



void save_trotter_copies( int **JJ, int M_times_R_guess, char *filename)
{
    int i,j;
    FILE *f = fopen(filename, "wt");
    if(f==NULL) exit(14);
    for(i=0; i<M_times_R_guess; i++)
    {
        for(j=0; j<3; j++) fprintf(f, "%d ", JJ[i][j]);
        fprintf(f, "\n");
    }
    fclose(f);


}



















