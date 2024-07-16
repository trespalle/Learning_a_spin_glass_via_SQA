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

int secuencia_rand[2048], index1[2048], index2[2048], iaux;


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
void MC_STEP(int **JJ, int M, int *columna_i_new, int *columna_j_new, char **D, int **lambda, int C, int N, double **FUNK, int zmax, double teff);
void print_couplings_mean(int **JJ, int M);
double GAMMA(int **JJ_0, int **JJ, int M);
void create_name(char *namefile, int i_seed);
void gamma_stats(double *gamma_vector, int n_seeds);
void lee_input(char *namefile, int *N_dir, int *C_dir, int *n_temps_dir, int **temps_dir, int *n_seeds_dir, double *teff_dir, int *tau_dir, int *first_seed_dir);
void change_g_name(char *name_g_file, int i_seed);
void read_geometry(char *namegfile, int *M_dir, int *zmax_dir, int **vecinos, int ***JJ_0_dir, int N);
void update_name(char *namefile, int temp);
void save_results(int *temps, double *av_gamma, double *av_gamma2, int n_temps, int n_seeds, int first_seed);


int main()
{
    ini_ran(123456789);
    ini_r1279();
    double teff, delta_teff, teff_ini;
    int mu, i, j, k, N, C, M, zmax, tau, i_seed, n_seeds;
    int n_temps, i_temps, first_seed;
    int *integer_temps;
    double aux_gamma;
    double temperatura;
    double amigo;


    ///We read the values of the parameters
    lee_input("input_SCA.txt", &N, &C, &n_temps, &integer_temps, &n_seeds, &teff_ini, &tau, &first_seed);

    printf("N = %d\nC = %d\nn_temps =%d\n", N, C, n_temps);
    for(i=0; i<n_temps; i++) printf("%d ", integer_temps[i]);
    printf("\nn_seeds = %d\nteff_ini = %lf\ntau = %d\nfirst_seed = %d\n", n_seeds, teff_ini, tau, first_seed);



    ///We create the main vectors and matrices:
    int **vecinos = (int**) malloc(N*sizeof(int*));
    int **JJ_0; /// real couplings (the ones we have to estimate)
    int **JJ; /// free variables of our simulated annealing

    char **D = (char**) malloc(C*sizeof(char*)); ///here we store the sampled equilibrium configurations
    for(i=0; i<C; i++) D[i] = (char*) malloc(N*sizeof(char));
    int **lambda = (int**) malloc(C*sizeof(int*)); /// argument of the hyperbolic tangent
    for(i=0; i<C; i++) lambda[i] = (int*) malloc(N*sizeof(int));
    double **FUNK = (double**) malloc(2*sizeof(double*)); ///log(...)

    int *columna_i_new = (int*) malloc(C*sizeof(int));
    int *columna_j_new = (int*) malloc(C*sizeof(int));
    double *av_gamma = (double*) malloc(n_temps*sizeof(double));
    double *av_gamma2 = (double*) malloc(n_temps*sizeof(double));

    ///we set to zero the averaged_gamma(T) vector:
    for(i=0; i<n_temps; i++) av_gamma[i] = av_gamma2[i] = 0.0;


    ///names of files;
    char namefile[] = "t000/D_000.dat";
    char namegfile[] = "GEOMETRY/G_000.dat";
    printf("%d\n", n_seeds);

    for(i_seed = 0; i_seed<n_seeds; i_seed++)
    {
        create_name(namefile, (first_seed + i_seed));

        //printf("%s\n", namefile);

        ///we read the geometry:
        change_g_name(namegfile, first_seed + i_seed);
        read_geometry(namegfile, &M, &zmax, vecinos, &JJ_0, N);
        JJ = (int**) malloc(M*sizeof(int*));
        for(i=0; i<M; i++) JJ[i] = (int*) malloc(3*sizeof(int));
        for(i=0; i<M; i++) for(j=1; j<3; j++) JJ[i][j] = JJ_0[i][j];

        //muestra_lista(vecinos, N);

        ///we create FUNK
        for(i=0; i<2; i++) FUNK[i] = (double*) malloc((2*zmax+1)*sizeof(double));


        ///temperature loop:
        for(i_temps = 0; i_temps<n_temps; i_temps++)
        {
            ///we read sampled configurations (D matrix)
            update_name(namefile, integer_temps[i_temps]);
            printf("%s\n", namefile);
            read_D(namefile, i_seed, D, C, N);
            //comprobacion(D, C, N);
            //for(i=0; i<N; i++) printf("%d\n", D[1][i]);


            ///initial random guess:
            for(i=0; i<M; i++)
            {
                if(Rand()>0.5) JJ[i][0] = 1;
                else JJ[i][0] = -1; ///random initial guess of the couplings
            }


            ///we tabulate lambda and FUNK
            tabulate_lambda(lambda, D, C, N, JJ, M);
            //for(i=0; i<N; i++) printf(" %d\n", lambda[15][i]);
            temperatura = integer_temps[i_temps]*0.01;
            tabulate_FUNK(FUNK, zmax, temperatura);
            //for(i=0; i<(2*zmax+1); i++) printf("%lf %lf\n", FUNK[0][i], FUNK[1][i]);

            ///MONTE CARLO SIMULATION:

            teff = teff_ini ;
            delta_teff = teff_ini/tau;
            for(i=0; i<tau; i++)
            {
                MC_STEP(JJ, M, columna_i_new, columna_j_new, D, lambda, C, N, FUNK, zmax, teff);
                teff -= delta_teff;
                amigo = compute_PLH(FUNK, zmax, D, lambda, C, N);
                printf("%lf\n", amigo);
                //print_couplings_mean(JJ, M);
            }
            //print_couplings_mean(JJ, M);

            aux_gamma = GAMMA(JJ_0, JJ, M);
            av_gamma[i_temps] += aux_gamma;
            av_gamma2[i_temps] += aux_gamma*aux_gamma;


        }

        ///we free those vectors that are created again at each iteration:
        for(i=0; i<N; i++) free(vecinos[i]);
        for(i=0; i<M; i++)
        {
            free(JJ[i]);
            free(JJ_0[i]);
        }
        for(i=0; i<2; i++) free(FUNK[i]);
        free(JJ_0);
        free(JJ);

    }


    save_results(integer_temps, av_gamma, av_gamma2, n_temps, n_seeds , first_seed);
    //for(i=0; i<n_temps; i++) av_gamma[i]/=n_seeds;
    //for(i=0; i<n_temps; i++) printf("%lf\n", av_gamma[i]);


    for(i=0; i<C; i++)
    {
        free(D[i]);
        free(lambda[i]);
    }
    free(D);
    free(lambda);
    free(FUNK);
    free(columna_i_new);
    free(columna_j_new);

    free(integer_temps);
    free(av_gamma);
    free(av_gamma2);
    free(vecinos);



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
    for(i=0; i<2; i++)
    {
        aux_spin = -1 +2*i;
        for(j=-zmax; j<=zmax; j++) FUNK[i][j+zmax] = log(0.5*(1.0 + aux_spin*tanh(((double)j)/temp)));
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

void MC_STEP(int **JJ, int M, int *columna_i_new, int *columna_j_new, char **D, int **lambda, int C, int N, double **FUNK, int zmax, double teff)
{
    char flag = 0;
    int e, k, i, j, mu, coupling, i_1, j_1;
    double DELTA_PLH, omega;
    double suma;

    for(e=0; e<M; e++)
    {
        do
        {
            k = M*r1279();
        } while(k==M);

        i = JJ[k][1];
        j = JJ[k][2];
        coupling =  JJ[k][0];
        ///up to now, we have chosen a coupling at random
        //printf("Seleccionamos el coupling numero %d que vale %d\n", k, JJ[k][0]);
        //printf("que es entre los espines %d y %d\n", i, j);
        for(mu=0; mu<C; mu++)
        {
            columna_i_new[mu] = lambda[mu][i] -2*coupling*D[mu][j];
            columna_j_new[mu] = lambda[mu][j] -2*coupling*D[mu][i];
        } ///(only two columns of the lambda matrix change when you change coupling J_ij)

        ///WE COMPUTE DELTA_PLH:
        DELTA_PLH = 0.0;
        for(mu=0; mu<C; mu++)
        {

            i_1 = (D[mu][i]+1)/2;
            j_1 = (D[mu][j]+1)/2;

            DELTA_PLH += (FUNK[i_1][columna_i_new[mu]+zmax] + FUNK[j_1][columna_j_new[mu]+zmax]);
            DELTA_PLH -= (FUNK[i_1][lambda[mu][i]+zmax] + FUNK[j_1][lambda[mu][j]+zmax]);
        }
        DELTA_PLH/=C;
        //printf("DELTA_PLH vale %lf\n", DELTA_PLH);

        ///LET'S SEE IF WE ACCEPT OR NOT (METROPOLIS ALGORITHM)
        flag = 0;
        if(DELTA_PLH>0.0) flag = 1;
        else
        {
            omega = r1279();
            //printf("El factor Boltzmann vale %lf y omega vale %lf\n",exp(DELTA_PLH/teff), omega);
            if(omega<exp(DELTA_PLH/teff)) flag = 1;

        }

        ///WHAT HAPPENS IF WE ACCEPT:
        if(flag == 1)
        {
            //suma = compute_PLH(FUNK, zmax,  D, lambda, C, N);
            //printf("Aceptamos el cambio, teff = %lf\n", teff);
            JJ[k][0] = -JJ[k][0]; ///we change the coupling
            for(mu=0; mu < C; mu++)
            {
                lambda[mu][i] = columna_i_new[mu];
                lambda[mu][j] = columna_j_new[mu];
            }
            //printf("Delta PLH real = %lf\n", suma - compute_PLH(FUNK, zmax, D, lambda, C, N));
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

double GAMMA(int **JJ_0, int **JJ, int M)
{
    double suma = 0.0;
    int k;

    for(k=0; k<M; k++) suma += (JJ_0[k][0]-JJ[k][0])*(JJ_0[k][0]-JJ[k][0]);
    return suma/M;
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
        name_g_file[13-i] += 48; ///convert to ASCII code
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

    (*JJ_0_dir) = (int**) malloc((*M_dir)*sizeof(int*));
    for(i=0; i<(*M_dir); i++) (*JJ_0_dir)[i] = (int*) malloc(3*sizeof(int));

    create_couplings_list(vecinos, N, *zmax_dir, *JJ_0_dir, *M_dir); ///the neighbours list is "signed" (nbrs appear with a - sign if the coupling is antiferromagnetic)
    ///we "unsign" the neigbours list:
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
        namefile[3-i] += 48; ///we convert to ASCII code
    }

}

void save_results(int *temps, double *av_gamma, double *av_gamma2, int n_temps, int n_seeds, int first_seed)
{
    int i, final_seed, dummy;
    char namefile[] = "av_gamma_vs_T_000_to_000.txt";

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
    for(i=0; i<n_temps; i++) fprintf(f,"%lf %lf %lf %d\n", temps[i]*0.01, av_gamma[i], av_gamma2[i], n_seeds);
    fclose(f);
}





















