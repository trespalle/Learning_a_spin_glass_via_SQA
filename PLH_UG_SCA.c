#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define NormRANu (2.3283063671E-10F)
#define normNorm (4.656612873E-10F) ///1/2^31, para que r1279 nunca valga 1
#define NBITM1 31
#define dim 2
#define z 4

#define bond_penalty 0.010F
#define M_aver 200


unsigned int larueda[256], number = 0;
unsigned char ind_ran, i1, i2, i3;

unsigned int secuencia_rand[2048], index1[2048], index2[2048], iaux;


void ini_ran(int seed);
double Rand(void);

void crea_lista_vecinos(int** vecinos, int L, int N);
void muestra_lista(int **vecinos, int N);
void ini_r1279();
double r1279(void);

void read_D(char *namefile, int i_seed, char **D, int C, int N);
void comprobacion(char **D, int C, int N);
void create_couplings_list(int**vecinos, int N, int zmax, int**JJ, int M);
void tabulate_lambda (int **lambda, char **D, int C, int N, int **JJ, int M);
void tabulate_FUNK (double **FUNK, int zmax, double temp);
double compute_PLH( double **FUNK, int zmax, char **D, int **lambda, int C, int N);
void MC_STEP(int **JJ, int M, int *columna_i_new, int *columna_j_new, char **D, int **lambda, int C, int N, double **FUNK, int zmax, double teff, int **all_coups, int *Mguess_dir, int M_true);
void print_couplings_mean(int **JJ, int M);
void create_name(char *namefile, int i_seed);
void gamma_stats(double *gamma_vector, int n_seeds);
void lee_input(char *namefile, int *N_dir, int *C_dir, int *n_temps_dir, int **temps_dir, int *n_seeds_dir, double *teff_dir, int *tau_dir, int *first_seed_dir);
void change_g_name(char *name_g_file, int i_seed);
void read_geometry(char *namegfile, int *M_dir, int *zmax_dir, int **vecinos, int ***JJ_0_dir, int N);
void update_name(char *namefile, int temp);
void save_results(int *temps, double *av_gamma, double *av_gamma2, double *av_sensit, double *av_sensit2, double *av_precis, double *av_precis2, double *av_F1, double *av_F1_2, int n_temps, int n_seeds, int first_seed);
void random_initial_couplings(int **all_coups, int **JJ, int Mmax, int *Mguess_dir, double alpha);
void append_new_coupling(int **all_coups, int **JJ, int *M_guess_dir, int k, int coupling);
void delete_coupling(int **all_coups, int **JJ, int *M_guess_dir, int k, int N);
double GAMMA(int **real_all_coups, int Mmax, int M, int **all_coups, int Mguess, int **JJ_0, int **JJ, int N);
double F1_score(int **real_all_coups, int Mmax, int **all_coups, int **JJ_0, int **JJ, double *sens_dir, double *precis_dir);



int main()
{
    ini_ran(423456789);
    ini_r1279();
    double teff, delta_teff, teff_ini;
    int mu, i, j, k, N, C, M, Mguess, Mmax, zmax, tau, i_seed, n_seeds;
    int n_temps, i_temps, first_seed;
    int *integer_temps;
    double aux_gamma, aux_sensit, aux_precis, aux_F1, suma;
    double temperatura;



    ///Leemos los parámetros
    lee_input("input_SCA.txt", &N, &C, &n_temps, &integer_temps, &n_seeds, &teff_ini, &tau, &first_seed);

    /*printf("N = %d\nC = %d\nn_temps =%d\n", N, C, n_temps);
    for(i=0; i<n_temps; i++) printf("%d ", integer_temps[i]);
    printf("\nn_seeds = %d\nteff_ini = %lf\ntau = %d\nfirst_seed = %d\n", n_seeds, teff_ini, tau, first_seed);*/

    Mmax = N*(N-1)/2;



    ///Generamos las principales matrices y vectores:
    int **vecinos = (int**) malloc(N*sizeof(int*));
    int **JJ_0; /// couplings reales (los que queremos inferir)
    int **JJ; /// variables libres de nuestro simulated annealing
    JJ = (int**) malloc(Mmax*sizeof(int*));
    for(i=0; i<Mmax; i++) JJ[i] = (int*) malloc(3*sizeof(int));
    int **all_coups = (int**) malloc(Mmax*sizeof(int*));
    for(i=0; i<Mmax; i++) all_coups[i] = (int*) malloc(3*sizeof(int));
    int **real_all_coups = (int**) malloc(Mmax*sizeof(int*));
    for(i=0; i<Mmax; i++) real_all_coups[i] = (int*) malloc(3*sizeof(int));


    k=0;
    for(i=0; i<N; i++) for(j=(i+1); j<N; j++)
    {
        all_coups[k][1] = i;
        all_coups[k][2] = j;
        all_coups[k][0] = -1; ///por defecto esta desactivado
        k++;
    }



    char **D = (char**) malloc(C*sizeof(char*)); ///aquí se guardan las configuraciones sampleadas
    for(i=0; i<C; i++) D[i] = (char*) malloc(N*sizeof(char));
    int **lambda = (int**) malloc(C*sizeof(int*)); /// el argumento de la tangente hiperbolica
    for(i=0; i<C; i++) lambda[i] = (int*) malloc(N*sizeof(int));
    double **FUNK = (double**) malloc(2*sizeof(double*)); ///log(...)

    int **auxlambda = (int**) malloc(C*sizeof(int*)); /// el argumento de la tangente hiperbolica
    for(i=0; i<C; i++) auxlambda[i] = (int*) malloc(N*sizeof(int));

    int *columna_i_new = (int*) malloc(C*sizeof(int));
    int *columna_j_new = (int*) malloc(C*sizeof(int));
    double *av_gamma = (double*) malloc(n_temps*sizeof(double));
    double *av_gamma2 = (double*) malloc(n_temps*sizeof(double));
    double *av_sensit = (double*) malloc(n_temps*sizeof(double));
    double *av_sensit2 = (double*) malloc(n_temps*sizeof(double));
    double *av_precis = (double*) malloc(n_temps*sizeof(double));
    double *av_precis2 = (double*) malloc(n_temps*sizeof(double));
    double *av_F1 = (double*) malloc(n_temps*sizeof(double));
    double *av_F1_2 = (double*) malloc(n_temps*sizeof(double));

    ///inicialicamos a cero el vector averaged_gamma(T):
    for(i=0; i<n_temps; i++) av_gamma[i] = av_gamma2[i] = av_sensit[i] = av_sensit2[i] = 0.0;
    for(i=0; i<n_temps; i++) av_precis[i] = av_precis2[i] = av_F1[i] = av_F1_2[i] = 0.0;


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
        zmax = N;


        //muestra_lista(vecinos, N);

        ///creamos FUNK
        for(i=0; i<2; i++) FUNK[i] = (double*) malloc((2*zmax+1)*sizeof(double));


        ///bucle en temperaturas:
        for(i_temps = 0; i_temps<n_temps; i_temps++)
        {
            ///leemos las configuraciones sampleadas (matriz D)
            update_name(namefile, integer_temps[i_temps]);
            //printf("%s\n", namefile);
            read_D(namefile, i_seed, D, C, N);
            //comprobacion(D, C, N);
            //for(i=0; i<N; i++) printf("%d\n", D[1][i]);


            ///initial random guess:
            random_initial_couplings(all_coups, JJ, Mmax, &Mguess, 0.975); ///se inventa qué couplings están activados y cuanto valen

            ///tabulamos lambda y FUNK
            tabulate_lambda(lambda, D, C, N, JJ, Mguess);
            //for(i=0; i<N; i++) printf(" %d\n", lambda[15][i]);
            temperatura = integer_temps[i_temps]*0.01;
            tabulate_FUNK(FUNK, zmax, temperatura);
            //for(i=0; i<(2*zmax+1); i++) printf("%lf %lf\n", FUNK[0][i], FUNK[1][i]);

            ///SIMULACION DE MONTECARLO:

            teff = teff_ini ;
            delta_teff = teff_ini/tau;
            for(i=0; i<tau; i++)
            {
                MC_STEP(JJ, Mmax, columna_i_new, columna_j_new, D, lambda, C, N, FUNK, zmax, teff, all_coups, &Mguess, M);
                teff -= delta_teff;
                //printf("%d\n", Mguess);
                //printf("%lf\n", teff);
                //print_couplings_mean(JJ, M);
            }
            //print_couplings_mean(JJ, M);
            /*suma = 0.0;
            tabulate_lambda(auxlambda, D,C, N, JJ, Mguess);
            for(mu=0; mu<C; mu++) for(i=0; i<N; i++) suma+=fabs(auxlambda[mu][i]-lambda[mu][i]);
            printf("%lf\n",suma);*/



            aux_gamma = GAMMA(real_all_coups, Mmax, M, all_coups, Mguess, JJ_0, JJ, N);
            aux_F1 = F1_score(real_all_coups, Mmax, all_coups, JJ_0, JJ, &aux_sensit, &aux_precis);

            av_gamma[i_temps] += aux_gamma;
            av_gamma2[i_temps] += aux_gamma*aux_gamma;
            av_sensit[i_temps] += aux_sensit;
            av_sensit2[i_temps] += aux_sensit*aux_sensit;
            av_precis[i_temps] += aux_precis;
            av_precis2[i_temps] += aux_precis*aux_precis;
            av_F1[i_temps] += aux_F1;
            av_F1_2[i_temps] += aux_F1*aux_F1;



        }

        ///liberamos los vectores que se crean de iteracion a iteracion:
        for(i=0; i<N; i++) free(vecinos[i]);
        for(i=0; i<M; i++)
        {
            free(JJ_0[i]);
        }
        for(i=0; i<2; i++) free(FUNK[i]);
        free(JJ_0);

    }


    save_results(integer_temps, av_gamma, av_gamma2, av_sensit, av_sensit2, av_precis, av_precis2, av_F1, av_F1_2, n_temps, n_seeds , first_seed);
    //for(i=0; i<n_temps; i++) av_gamma[i]/=n_seeds;
    //for(i=0; i<n_temps; i++) printf("%lf\n", av_gamma[i]);


    for(i=0; i<C; i++)
    {
        free(D[i]);
        free(lambda[i]);
        free(auxlambda[i]);
    }
    free(D);
    free(lambda);
    free(auxlambda);
    free(FUNK);
    free(columna_i_new);
    free(columna_j_new);

    for(i=0; i<Mmax; i++)
    {
        free(JJ[i]);
        free(all_coups[i]);
        free(real_all_coups[i]);
    }
    free(JJ);
    free(all_coups);
    free(real_all_coups);

    free(integer_temps);
    free(av_gamma);
    free(av_gamma2);
    free(vecinos);
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
    numerin = (secuencia_rand[iaux]>>1); ///&2147483647
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

void tabulate_lambda (int **lambda, char **D, int C, int N, int **JJ, int M)
{
    int mu, i, k;

    for(mu=0; mu<C; mu++) for(i=0; i<N; i++) lambda[mu][i] = 0;

    for(k=0; k<M; k++)
    {
        for(mu=0; mu<C; mu++)
        {
            lambda[mu][JJ[k][1]] += JJ[k][0]*D[mu][JJ[k][2]];
            lambda[mu][JJ[k][2]] += JJ[k][0]*D[mu][JJ[k][1]];
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

double compute_PLH( double **FUNK, int zmax, char **D, int **lambda, int C, int N)
{
    double suma, PLH = 0.0;
    int mu, i, aux_ind;

    for(i=0; i<N; i++)
    {
        suma = 0.0;
        for(mu=0; mu<C; mu++)
        {
            aux_ind = (1+D[mu][i])/2;
            suma += FUNK[aux_ind][lambda[mu][i]+zmax];
        }
        PLH += suma/C;
    }
    return PLH;
}

void MC_STEP(int **JJ, int M, int *columna_i_new, int *columna_j_new, char **D, int **lambda, int C, int N, double **FUNK, int zmax, double teff, int **all_coups, int *Mguess_dir, int M_true)
{
    char flag = 0;
    int e, k, i, j, mu, coupling, new_coupling, delta_coup, i_1, j_1;
    double DELTA_PLH, omega;
    double suma, pen_term;
    //printf("PEPE");

    for(e=0; e<M; e++)
    {
        do
        {
            k = M*Rand();
        } while(k==M);

        i = all_coups[k][1];
        j = all_coups[k][2];
        if(all_coups[k][0]<0) coupling = 0;
        else coupling = JJ[all_coups[k][0]][0];
        ///hasta aquí lo que hemos hecho es seleccionar un coupling al azar
        //printf("Seleccionamos el coupling numero %d que vale %d\n", k, coupling);
        //printf("que es entre los espines %d y %d\n", i, j);

        ///PROPONEMOS UN CAMBIO:
        if(Rand()>0.5) new_coupling = coupling +1;
        else new_coupling = coupling -1;
        if(new_coupling == -2) new_coupling = 1;
        if(new_coupling == 2) new_coupling = -1;
        delta_coup = new_coupling - coupling;
        //printf("Delta_coup vale %d\n", delta_coup);

        for(mu=0; mu<C; mu++)
        {
            columna_i_new[mu] = lambda[mu][i] +delta_coup*D[mu][j];
            columna_j_new[mu] = lambda[mu][j] +delta_coup*D[mu][i];
        } ///(solo cambian dos columnas de la matriz lambda al cambiar el coupling J_ij

        ///CÁLCULO DE DELTA_PLH:
        DELTA_PLH = 0.0;
        for(mu=0; mu<C; mu++)
        {

            i_1 = (D[mu][i]+1)/2;
            j_1 = (D[mu][j]+1)/2;

            DELTA_PLH += (FUNK[i_1][columna_i_new[mu]+zmax] + FUNK[j_1][columna_j_new[mu]+zmax]);
            DELTA_PLH -= (FUNK[i_1][lambda[mu][i]+zmax] + FUNK[j_1][lambda[mu][j]+zmax]);
            //if(e==4808) printf("Los indices funk nuevos son %d y %d, los viejos son %d y %d\n", columna_i_new[mu],columna_j_new[mu], lambda[mu][i], lambda[mu][j]);
            //if(e==4808) printf("(Recuerden que zmax = %d)\n", zmax);
            //if(e==4808) printf("Seguimiento: ahora DELTAPLH vale %lf\n", DELTA_PLH);
        }
        DELTA_PLH/=C;
        //printf("Iteracion %d: DELTA_PLH vale %lf,\n", e, DELTA_PLH);

        ///el termino de penalty afecta a DELTA_PLH:
        if(new_coupling*coupling==0) ///si se ha creado o destruido algun enlace
        {
            if((*Mguess_dir)<M_true) pen_term = (abs(new_coupling)-abs(coupling))*bond_penalty; ///aquí hacemos que los procesos de creacion aumenten PLH
            if((*Mguess_dir)>M_true) pen_term = -(abs(new_coupling)-abs(coupling))*bond_penalty; ///aquí hacemos que los procesos de destruccion disminuyan PLH
            if((*Mguess_dir)==M_true) pen_term = -bond_penalty; ///si  ya estas en el target, tanto crear como destruir enlaces se penaliza
            DELTA_PLH += pen_term;
        }
        //printf("Estamos en Mguess = %d y el man queria pasar de %d a %d, asi que le dije que el cambio en PLH seria %lf\n", *Mguess_dir, coupling, new_coupling, pen_term);

        ///VEMOS SI ACEPTAMOS O NO (ALGORITMO DE METROPOLIS)
        flag = 0;
        if(DELTA_PLH>0.0) flag = 1;
        else
        {
            omega = Rand();
            //printf("El factor Boltzmann vale %lf y omega vale %lf\n",exp(DELTA_PLH/teff), omega);
            if(omega<exp(DELTA_PLH/teff)) flag = 1;

        }

        ///LO QUE SUCEDE SI ACEPTAMOS:
        if(flag == 1)
        {
            //suma = compute_PLH(FUNK, zmax,  D, lambda, C, N);
            //printf("Aceptamos el cambio, teff = %lf\n", teff);
            ///hay 3 casos: creacion, destruccion y flip

            ///CREACION:
            if(coupling==0) append_new_coupling(all_coups, JJ, Mguess_dir, k, new_coupling);

            ///DESTRUCCION
            if(new_coupling==0) delete_coupling(all_coups, JJ, Mguess_dir, k, N);

            ///FLIP
            if(new_coupling*coupling!=0) JJ[all_coups[k][0]][0] = -JJ[all_coups[k][0]][0];


            for(mu=0; mu < C; mu++)
            {
                lambda[mu][i] = columna_i_new[mu];
                lambda[mu][j] = columna_j_new[mu];
            }
            //printf("Delta PLH real = %lf\n", compute_PLH(FUNK, zmax, D, lambda, C, N) -suma);
            //printf("Supuesto delta PLH = %lf\n", DELTA_PLH);
        }
        //else printf("No aceptamos el cambio\n");


    }
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

void lee_input(char *namefile, int *N_dir, int *C_dir, int *n_temps_dir, int **temps_dir, int *n_seeds_dir, double *teff_dir, int *tau_dir, int *first_seed_dir)
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

    fscanf(f, "\n%d\n%lf\n%d\n%d", n_seeds_dir, teff_dir, tau_dir, first_seed_dir);
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
    char namefile[] = "av_gamma_vs_T_000_to_000_bonpen001.txt";

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


    FILE *f = fopen(namefile, "wt");
    if(f==NULL) exit(first_seed);
    for(i=0; i<n_temps; i++) fprintf(f,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %d\n", temps[i]*0.01, av_gamma[i], av_gamma2[i], av_sensit[i], av_sensit2[i], av_precis[i], av_precis2[i], av_F1[i], av_F1_2[i], n_seeds);
    fclose(f);
}

void random_initial_couplings(int **all_coups, int **JJ, int Mmax, int *Mguess_dir, double alpha)
{
    ///alpha es la probabilidad de que un coupling valga 0. Se recomienda fijar alta, porque normalmente M/Mmax es pequeño

    int i, j, our_M;

    for(i=0; i<Mmax; i++) for(j=0; j<3; j++) JJ[i][j] = 0;
    our_M = 0;

    for(i=0; i<Mmax; i++)
    {
        if(Rand()<alpha) all_coups[i][0] = -1; ///está desactivado
        else
        {
            all_coups[i][0] = our_M; ///posicion del coupling en la lista de couplings activos
            JJ[our_M][1] = all_coups[i][1];
            JJ[our_M][2] = all_coups[i][2];
            if(Rand()>0.5) JJ[our_M][0] = 1;
            else JJ[our_M][0] = -1;
            our_M++;
        }
    }
    *Mguess_dir = our_M;
}

void append_new_coupling(int **all_coups, int **JJ, int *M_guess_dir, int k, int coupling)
{
    int our_M = *M_guess_dir;
    JJ[our_M][0] = coupling;
    JJ[our_M][1] = all_coups[k][1];
    JJ[our_M][2] = all_coups[k][2];
    all_coups[k][0] = our_M;
    our_M++;
    *M_guess_dir = our_M;
}

void delete_coupling(int **all_coups, int **JJ, int *M_guess_dir, int k, int N)
{
    int i, j, a, b, position,index, our_M, twoN_3;
    our_M = *M_guess_dir;
    position = all_coups[k][0];
    twoN_3 = 2*N-3;



    for(i=position; i<our_M; i++) for(j=0; j<3; j++) JJ[i][j] = JJ[i+1][j];
    all_coups[k][0] = -1; ///coupling desactivado
    our_M--;
    ///actualizamos las posiciones de todos los demas:
    for(i=position; i<our_M; i++)
    {
        a = JJ[i][1];
        b = JJ[i][2];
        index = (-a*a +twoN_3*a -2 +2*b)/2;
        all_coups[index][0]--;
    }

    *M_guess_dir = our_M;
}


double GAMMA(int **real_all_coups, int Mmax, int M, int **all_coups, int Mguess, int **JJ_0, int **JJ, int N)
{
    int i, j, k, index, twoN_3;
    double suma, diff, dummy;
    int real_coup, guess_coup;

    twoN_3 = 2*N-3;

    ///Primero, construimos "real_all_coups":
    k=0;
    for(i=0; i<N; i++) for(j=i+1; j<N; j++)
    {
        real_all_coups[k][1] = i;
        real_all_coups[k][2] = j;
        real_all_coups[k][0] = -1; ///desactivado por defecto
        k++;
    }


    for(k=0; k<M; k++)
    {
        i = JJ_0[k][1];
        j = JJ_0[k][2];
        index = (-i*i +twoN_3*i -2 +2*j)/2;
        //printf("%d\n", Mguess);
        real_all_coups[index][0] = k;
    }


    ///Ahora, calculamos GAMMA:

    suma = 0.00;
    for(k=0; k<M; k++)
    {
        i = JJ_0[k][1];
        j = JJ_0[k][2];
        index = (-i*i +twoN_3*i -2 +2*j)/2;
        if(all_coups[index][0]<0) diff = 1.0;
        else diff = (JJ_0[k][0] - JJ[all_coups[index][0]][0])*1.0;

        suma += diff*diff;
    }
    ///Acabamos de sumar la contribucion de los couplings que son no nulos en la lista original.
    ///Ahora vamos a sumar la constribucion de los couplings que no son nulos en la "lista estimada"
    ///y cuya contribucion no hemos contado ya:
    for(k=0; k<Mguess; k++)
    {
        i = JJ[k][1];
        j = JJ[k][2];
        index = (-i*i +twoN_3*i -2 +2*j)/2;
        if(real_all_coups[index][0]<0) suma+=1.0;
    }

    suma/=Mmax;
    //printf("Gamma calculado a mi manera es %lf\n", suma);

    /*///comprobacion:
    dummy = 0.0;
    for(k=0; k<Mmax; k++)
    {
        if(real_all_coups[k][0]<0) real_coup = 0;
        else real_coup = JJ_0[real_all_coups[k][0]][0];
        if(all_coups[k][0]<0) guess_coup = 0;
        else guess_coup = JJ[all_coups[k][0]][0];
        diff = (guess_coup-real_coup)*1.0;
        dummy += diff*diff;
    }

    dummy/=Mmax;
    printf("Gamma calculado de forma directa es %lf\n", dummy);*/
    return suma;
}

double F1_score(int **real_all_coups, int Mmax, int **all_coups, int **JJ_0, int **JJ, double *sens_dir, double *precis_dir)
{
    int real_coup, guess_coup;
    int k;
    int n_correct_nonzero_predictions, n_nonzero_bonds, n_nonzero_predictions;


    n_nonzero_bonds = n_nonzero_predictions = n_correct_nonzero_predictions = 0;
    for(k=0; k<Mmax; k++)
    {
        if(real_all_coups[k][0]<0) real_coup = 0;
        else real_coup = JJ_0[real_all_coups[k][0]][0];
        if(all_coups[k][0]<0) guess_coup = 0;
        else guess_coup = JJ[all_coups[k][0]][0];

        if(real_coup!=0) n_nonzero_bonds++;
        if(guess_coup!=0) n_nonzero_predictions++;
        if((real_coup==guess_coup)&&(real_coup!=0)) n_correct_nonzero_predictions++;

    }
    *sens_dir = (n_correct_nonzero_predictions*1.00)/n_nonzero_bonds;
    *precis_dir = (n_correct_nonzero_predictions*1.00)/n_nonzero_predictions;
    return 2.0*(*sens_dir)*(*precis_dir)/((*sens_dir)+(*precis_dir));
}
























