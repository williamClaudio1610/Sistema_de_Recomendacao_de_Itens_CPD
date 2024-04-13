#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#define ALEATORIO ((double)random() / (double)RAND_MAX)

double** criar_matriz(int linha, int coluna, int val_inicio) {
    double** M = (double**) malloc(linha * sizeof(double*));
    for (int i = 0; i < linha; i++) {
        M[i] = (double*) malloc(coluna * sizeof(double));
        for (int j = 0; j < coluna; j++) {
            M[i][j] = val_inicio;
        }
    }
    return M;
}

void destruir_matriz(double** M, int linha) {
    for (int i = 0; i < linha; i++) {
        free(M[i]);
    }
    free(M);
}

void copiar_matriz(double** origem, double** destino, int linha, int coluna) {
    for (int i = 0; i < linha; i++) {
        for (int j = 0; j < coluna; j++) {
            destino[i][j] = origem[i][j];
        }
    }
}

void preenche_aleatorio_LR(double** L, double** R, int nU, int nI, int nF)
{
    srand(0);
    for (int i = 0; i < nU; i++)
        for (int j = 0; j < nF; j++)
            L[i][j] = ALEATORIO / (double)nF;
    for (int i = 0; i < nF; i++)
        for (int j = 0; j < nI; j++)
            R[i][j] = ALEATORIO / (double)nF;
}

void multiplicacao(double** L, double** R, double** B, int nU, int nI, int nF) {
    for (int i = 0; i < nU; i++) {
        for (int j = 0; j < nI; j++) {
            double soma = 0.0;
            for (int k = 0; k < nF; k++) {
                soma += L[i][k] * R[k][j];
            }
            B[i][j] = soma;
        }
    }
}

void gerarL(double** A, double**L, double**R, double** B, double** L_prox, int nU, int nI, int nF, double alfa) {
    double  entraL=0;
#pragma omp parallel for collapse(2)
    for (int i = 0; i < nF; i++) {
        for (int j = 0; j < nU; j++) {
            double somatorioL = 0.0;
            for (int k = 0; k < nI; k++) {
                if (A[j][k] > 0)
                    somatorioL += 2 * (A[j][k] - B[j][k]) * (-R[i][k]);
            }
            L_prox[j][i] = L[j][i] - alfa * somatorioL;
        }
    }
}

void gerarR(double** A, double**L, double**R, double** B, double** R_prox, int nU, int nI, int nF, double alfa) {
    double entraR=0;
#pragma omp parallel for collapse(2)
    for (int i = 0; i < nF; i++) {
        for (int j = 0; j < nI; j++) {
            double somatorioR= 0.0;
            for (int k = 0; k < nU; k++) {
                if (A[k][j] > 0)
                    somatorioR += 2 * (A[k][j] - B[k][j]) * (-L[k][i]);
            }
            R_prox[i][j] = R[i][j] - alfa * somatorioR;
        }
    }
}

int main(int argc, char** argv)
{
    if (argc != 2) {
        printf("Passe como argumento o nome do ficheiro dos dados de entrada\n");
        return -1;
    }

    double start, end, time;
    start = omp_get_wtime();

    FILE *arquivo, *arquivo1;

    arquivo = fopen(argv[1], "r");
    arquivo1 = fopen("output.txt", "w");
    int numIteracoes, numCarac, linha, coluna, numDifZero;
    double alfa;

    fscanf(arquivo, "%d", &numIteracoes);
    fscanf(arquivo, "%lf", &alfa);
    fscanf(arquivo, "%d", &numCarac);
    fscanf(arquivo, "%d%d%d", &linha, &coluna, &numDifZero);

    double** A = criar_matriz(linha, coluna, 0);
    double** B = criar_matriz(linha, coluna, 0);
    double** L = criar_matriz(linha, numCarac, 0);
    double** R = criar_matriz(numCarac, coluna, 0);

    double** L_prox = criar_matriz(linha, numCarac, 0);
    double** R_prox = criar_matriz(numCarac, coluna, 0);


    #pragma omp for
    for (int i = 1; i <= numDifZero; i++) {
        int y, x;
        double val;
        fscanf(arquivo, "%d%d%lf", &y, &x, &val);
        A[y][x] = val;
    }

    preenche_aleatorio_LR(L, R, linha, coluna, numCarac);

    multiplicacao(L, R, B, linha, coluna, numCarac);

    int itera = 0;

    while(itera < numIteracoes){
        gerarL(A, L, R, B, L_prox, linha, coluna, numCarac, alfa);
        gerarR(A, L, R, B, R_prox, linha, coluna, numCarac, alfa);
        copiar_matriz(L_prox, L, linha, numCarac);
        copiar_matriz(R_prox, R, numCarac, coluna);
        multiplicacao(L, R, B, linha, coluna, numCarac);
        itera++;
    }

    end = omp_get_wtime();
/*
    for (int i = 0; i < linha; i++) {
        for (int j = 0; j < coluna; j++) {
            printf("%lf ", B[i][j]);
        }
        printf("\n");
    }
*/

    for (int i = 0; i < linha; i++) {
        double maior = 0;
        int indice = 0;
        for (int j = 0; j < coluna; ++j) {
            if(B[i][j] > maior && A[i][j] == 0){
                maior = B[i][j];
                indice =  j;
            }
        }

        fprintf(arquivo1 ,"%d\n", indice);
    }

    //printf("\n%d", numIteracoes);
    time = end - start;
    printf("\nTempo total: %f\n", time);

    // Obter o número de threads utilizadas
    int num_threads = omp_get_max_threads();
    //printf("Numero de threads utilizadas: %d\n", num_threads);


    destruir_matriz(A, linha);
    destruir_matriz(B, linha);
    destruir_matriz(L, linha);
    destruir_matriz(R, numCarac);

    //destruir_matriz(B_prox, linha);
    destruir_matriz(L_prox, linha);
    destruir_matriz(R_prox, numCarac);



    fclose(arquivo);
    fclose(arquivo1);
}
