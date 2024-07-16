#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>

void pilla_nombres_archivos(char *namefile, int *n_names_dir, int *long_max_dir, char ***names_inputs)
{
    FILE *f;
    int thegrefg, aitana, i;
    char dummy;
    int longmax = 0;

    f = fopen(namefile, "rt");
    if(f==NULL) exit(1);

    aitana = 0;
    while(feof(f)==0)
    {
        thegrefg = 0;
        while(((dummy = fgetc(f))!='\n')&&(dummy != EOF)) thegrefg++;
        if(thegrefg>longmax) longmax = thegrefg;
        aitana++;
    }
    fclose(f);

    *long_max_dir = longmax;
    *n_names_dir = aitana;

    printf("Numero de filas = %d\n", aitana);

    (*names_inputs) = (char**) malloc(aitana*sizeof(char*));
    for(i=0; i<aitana; i++) (*names_inputs)[i] = (char*) malloc((longmax+1)*sizeof(char));
    for(thegrefg = 0; thegrefg<aitana; thegrefg++) for(i=0; i<=longmax; i++) (*names_inputs)[thegrefg][i] = '\0';

    f = fopen(namefile, "rt");
    if(f==NULL) exit(2);

    aitana = 0;
    while(feof(f)==0)
    {
        i = 0;
        while(((dummy = fgetc(f))!='\n')&&(dummy != EOF))
        {
            (*names_inputs)[aitana][i] = dummy;
            printf("%c", dummy);
            i++;
        }
        aitana++;
        printf("\n");
    }
    fclose(f);


}



int main()
{

    char namefile[] = "names_inputs.txt";
    int longmax, n_names;
    char **names_inputs;
    int i, j, n_temps, counter;
    FILE *f;
    double gamma, gamma_2, sensit, sensit_2, precis, precis_2, F1, F1_2;


    n_temps = 15;

    double *av_gamma = (double*) malloc(n_temps*sizeof(double));
    double *av_gamma_2 = (double*) malloc(n_temps*sizeof(double));
    double *av_sensit = (double*) malloc(n_temps*sizeof(double));
    double *av_sensit_2 = (double*) malloc(n_temps*sizeof(double));
    double *av_precis = (double*) malloc(n_temps*sizeof(double));
    double *av_precis_2 = (double*) malloc(n_temps*sizeof(double));
    double *av_F1 = (double*) malloc(n_temps*sizeof(double));
    double *av_F1_2 = (double*) malloc(n_temps*sizeof(double));
    double *temps = (double*) malloc(n_temps*sizeof(double));
    int *n_mans = (int*) malloc(n_temps*sizeof(int));

    for(i=0; i<n_temps; i++) n_mans[i] = 0;
    for(i=0; i<n_temps; i++) av_gamma[i] = av_gamma_2[i] = temps[i] = 0.0;
    for(i=0; i<n_temps; i++) av_sensit[i] = av_sensit_2[i] = av_precis[i] = av_precis_2[i] = av_F1[i] = av_F1_2[i] = 0.0;



    pilla_nombres_archivos(namefile, &n_names, &longmax, &names_inputs);

    for(i=0; i<n_names; i++)
    {
        f = fopen(names_inputs[i], "rt");
        for(j=0; j<n_temps; j++)
        {
            fscanf(f, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %d\n", &temps[j], &gamma, &gamma_2, &sensit, &sensit_2, &precis, &precis_2, &F1, &F1_2, &counter);
            av_gamma[j] += gamma;
            av_gamma_2[j] += gamma_2;
            av_sensit[j] += sensit;
            av_sensit_2[j] += sensit_2;
            av_precis[j] += precis;
            av_precis_2[j] += precis_2;
            av_F1[j] += F1;
            av_F1_2[j] += F1_2;
            n_mans[j] += counter;
            printf("%d\n", counter);
        }
        fclose(f);
    }

    ///FINAL RESULTS

    printf("%lf %lf %lf %lf %lf %lf %lf %lf %lf %d\n", temps[0], av_gamma[0], av_gamma_2[0], av_sensit[0], av_sensit_2[0], av_precis[0], av_precis_2[0], av_F1[0], av_F1_2[0], n_mans[0]);

    for(i=0; i<n_temps; i++)
    {
        av_gamma[i]/=n_mans[i];
        av_gamma_2[i]/=n_mans[i];
        av_gamma_2[i] -= av_gamma[i]*av_gamma[i];
        av_gamma_2[i] = sqrt(av_gamma_2[i]/(n_mans[i]-1));

        av_sensit[i]/=n_mans[i];
        av_sensit_2[i]/=n_mans[i];
        av_sensit_2[i] -= av_sensit[i]*av_sensit[i];
        av_sensit_2[i] = sqrt(av_sensit_2[i]/(n_mans[i]-1));

        av_precis[i]/=n_mans[i];
        av_precis_2[i]/=n_mans[i];
        av_precis_2[i] -= av_precis[i]*av_precis[i];
        av_precis_2[i] = sqrt(av_precis_2[i]/(n_mans[i]-1));

        av_F1[i]/=n_mans[i];
        av_F1_2[i]/=n_mans[i];
        av_F1_2[i] -= av_F1[i]*av_F1[i];
        av_F1_2[i] = sqrt(av_F1_2[i]/(n_mans[i]-1));
    }

     printf("%lf %lf %lf %lf %lf %lf %lf %lf %lf %d\n", temps[0], av_gamma[0], av_gamma_2[0], av_sensit[0], av_sensit_2[0], av_precis[0], av_precis_2[0], av_F1[0], av_F1_2[0], n_mans[0]);


    f = fopen("final_results.txt", "wt");
    if(f==NULL) exit(3);
    for(i=0; i<n_temps; i++) fprintf(f, "%lf %lf %lf %lf %lf %lf %lf %lf %lf\n", temps[i], av_gamma[i], av_gamma_2[i], av_sensit[i], av_sensit_2[i], av_precis[i], av_precis_2[i], av_F1[i], av_F1_2[i]);
    fclose(f);






    for(i=0; i<n_names; i++) free(names_inputs[i]);
    free(names_inputs);
    free(av_gamma);
    free(av_gamma_2);
    free(temps);
    free(n_mans);
    free(av_sensit);
    free(av_sensit_2);
    free(av_precis);
    free(av_precis_2);
    free(av_F1);
    free(av_F1_2);


    return 0;
}






