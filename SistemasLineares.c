#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>

#include "utils.h"
#include "SistemasLineares.h"


/*!
  \brief Essa função calcula a norma L2 do resíduo de um sistema linear 

  \param SL Ponteiro para o sistema linear
  \param x Solução do sistema linear
  \param res Valor do resíduo

  \return Norma L2 do resíduo.
*/
real_t normaL2Residuo(SistLinear_t *SL, real_t *x, real_t *res)
{
    real_t sum = 0.0f;
    for (int i = 0; i < SL->n; i++)
        sum += pow(res[i], 2.0f); 
    return sqrt(sum);
}


/*!
  \brief Método da Eliminação de Gauss

  \param SL Ponteiro para o sistema linear
  \param x ponteiro para o vetor solução
  \param tTotal time gasto pelo método

  \return código de erro. 0 em caso de sucesso.
*/
int eliminacaoGauss (SistLinear_t *SL, real_t *x, double *tTotal)
{
    SistLinear_t *clone = copiar_SL(SL);
    int i, k, j, result;
    double m, time = timestamp();
    
    // Triangula usando o pivoteamento parcial
    for (k = 0; k < clone->n - 1; k++){
        unsigned int max_index = k;
        float max = fabs(clone->A[k][k]);
        for (i = k + 1; i < clone->n; i++)
            if (fabs(clone->A[i][k]) > max)
                max_index = i;

        if (max_index != k){ // Se terminou em um índice diferente de onde começou, troca
            real_t aux;
            for (int c = 0; c < SL->n; c++){
                aux = clone->A[k][c];
                clone->A[k][c] = clone->A[max_index][c];
                clone->A[max_index][c] = aux;
            }

            aux = clone->b[k];
            clone->b[k] = clone->b[max_index];
            clone->b[max_index] = aux;
        }

        for (i = k + 1; i < clone->n; i++){
            m = clone->A[i][k] / clone->A[k][k];
            if (invalid(m)){
                fprintf(stderr, "Gauss-Jordan floating point error.\n");
                liberaSistLinear(clone);
                return -1;
            }
            clone->A[i][k] = 0.0f;
            for (j = k + 1; j < clone->n; j++){
                clone->A[i][j] -= (m * clone->A[k][j]);
                if (invalid(clone->A[i][j])){
                    fprintf(stderr, "Gauss-Jordan floating point error.\n");
                    liberaSistLinear(clone);
                    return -1;
                }
            }
                
            clone->b[i] -= m * clone->b[k];
            if (invalid(clone->b[i])){
                fprintf(stderr, "Gauss-Jordan floating point error.\n");
                liberaSistLinear(clone);
                return -1;
            }
        }
    }

    result = retrosubs(clone);
    if (result){
        liberaSistLinear(clone);
        return -1;
    }

    *tTotal = timestamp() - time;
    memcpy(x, clone->b, sizeof(real_t) * clone->n);

    liberaSistLinear(clone);

    return 0;
}


/*!
  \brief Método de Jacobi

  \param SL Ponteiro para o sistema linear
  \param x ponteiro para o vetor solução. Ao iniciar função contém
            valor inicial
  \param tTotal time gasto pelo método

  \return código de erro. Um nr positivo indica sucesso e o nr
          de iterações realizadas. Um nr. negativo indica um erro:
          -1 (não converge) -2 (sem solução)
*/
int gaussJacobi (SistLinear_t *SL, real_t *x, double *tTotal)
{
    if (!jacobi_converge(SL)){
        fprintf(stderr, "Gauss-Jacobi doesn't converge.\n");
        return -1;
    }
        
    real_t* curr_iter = calloc(SL->n, sizeof(real_t)); // Valores usados na atual iteração
    must_alloc(curr_iter, __func__);

    real_t* next_iter = malloc(SL->n * sizeof(real_t)); // Valores a serem usados na próxima iteração
    must_alloc(next_iter, __func__);
    for (int i = 0; i < SL->n; i++) next_iter[i] = FLT_MAX; // Inicia o vetor com valores muito diferentes da primeira iteração

    real_t* aux_iter = calloc(SL->n, sizeof(real_t)); // Vetor auxiliar para conservar os valores a serem comparados entre iterações
    must_alloc(aux_iter, __func__);

    void **ptrs = malloc(sizeof(void*) * 3);
    ptrs[0] = curr_iter;
    ptrs[1] = next_iter; 
    ptrs[2] = aux_iter;

    real_t prev_diff = FLT_MAX; // A maior diferença entre as iterações

    int iter, i, k;
    double sum, time = timestamp();
    // Enquanto forem muito diferentes e iter não ultrapassou o limite de iterações, itera
    for (iter = 0; iter < MAXIT && too_different(curr_iter, next_iter, SL->n, SL->erro); iter++){
        memcpy(curr_iter, aux_iter, sizeof(real_t) * SL->n);
        for (i = 0; i < SL->n; i++){
            sum = 0.0f;
            for (k = 0; k < SL->n; k++){
                if (k != i) // Impede que some o pivô
                    sum += SL->A[i][k] * curr_iter[k];
                if (invalid(sum)){
                    fprintf(stderr, "Gauss-Jacobi floating point error.\n");
                    free_these(ptrs, 3);
                    return -3;
                }
            }

            if (SL->A[i][i] == 0.0f && (SL->b[i] - sum) != 0.0f){
                fprintf(stderr, "Gauss-Jacobi no solution.\n");
                free_these(ptrs, 3);
                return -2;
            }
            else
                next_iter[i] = (SL->b[i] - sum) / SL->A[i][i];

            if (invalid(next_iter[i])){
                fprintf(stderr, "Gauss-Jacobi floating point error.\n");
                free_these(ptrs, 3);
                return -3;
            }
        }
        memcpy(aux_iter, next_iter, sizeof(real_t) * SL->n);
    }

    *tTotal = timestamp() - time;
    memcpy(x, next_iter, sizeof(real_t) * SL->n);
    
    free_these(ptrs, 3);

    return iter;
}


/*!
  \brief Método de Gauss-Seidel

  \param SL Ponteiro para o sistema linear
  \param x ponteiro para o vetor solução. Ao iniciar função contém
            valor inicial
  \param tTotal time gasto pelo método

  \return código de erro. Um nr positivo indica sucesso e o nr
          de iterações realizadas. Um nr. negativo indica um erro:
          -1 (não converge) -2 (sem solução)
  */
int gaussSeidel (SistLinear_t *SL, real_t *x, double *tTotal)
{
    if (!seidel_converge(SL)){
        fprintf(stderr, "Gauss-Seidel doesn't converge.\n");
        return -1;
    }

    real_t* prev_iter = malloc(SL->n * sizeof(real_t)); // Valores da iteração anterior
    must_alloc(prev_iter, __func__);
    for (int i = 0; i < SL->n; i++) prev_iter[i] = FLT_MAX; // Inicia o vetor com valores muito diferentes da primeira iteração

    real_t* curr_iter = calloc(SL->n, sizeof(real_t)); // Valores da iteração atual
    must_alloc(curr_iter, __func__);

    void **ptrs = malloc(sizeof(void*) * 2);
    ptrs[0] = prev_iter;
    ptrs[1] = curr_iter; 

    real_t prev_diff = FLT_MAX; // A maior diferença entre as iterações

    int iter, i, k;
    double sum, time = timestamp();
    // Enquanto forem muito diferentes e iter não ultrapassou o limite de iterações, itera
    for (iter = 0; iter < MAXIT && too_different(prev_iter, curr_iter, SL->n, SL->erro); iter++){
        memcpy(prev_iter, curr_iter, sizeof(real_t) * SL->n);
        for (i = 0; i < SL->n; i++){
            sum = 0.0f;
            for (k = 0; k < SL->n; k++){
                if (k != i)
                    sum += SL->A[i][k] * curr_iter[k];
                if (invalid(sum)){
                    fprintf(stderr, "Gauss-Seidel floating point error.\n");
                    free_these(ptrs, 2);
                    return -3;
                }
            }

            if (SL->A[i][i] == 0.0f && (SL->b[i] != 0.0f)){
                fprintf(stderr, "No solution.\n");
                free_these(ptrs, 2);
                return -2;
            }
            else
                curr_iter[i] = (SL->b[i] - sum) / SL->A[i][i];

            if (invalid(curr_iter[i])){
                fprintf(stderr, "Gauss-Seidel floating point error.\n");
                free_these(ptrs, 2);
                return -3;
            }
        }
    }

    *tTotal = timestamp() - time;
    memcpy(x, curr_iter, sizeof(real_t) * SL->n);

    free_these(ptrs, 2);

    return iter;
}


/*!
  \brief Método de Refinamento

  \param SL Ponteiro para o sistema linear
  \param x ponteiro para o vetor solução. Ao iniciar função contém
            valor inicial para início do refinamento
  \param tTotal time gasto pelo método

  \return código de erro. Um nr positivo indica sucesso e o nr
          de iterações realizadas. Um nr. negativo indica um erro:
          -1 (não converge) -2 (sem solução)
  */
int refinamento(SistLinear_t *SL, real_t *x, double *tTotal)
{
    real_t *prev_iter = malloc(sizeof(real_t) * SL->n);
    must_alloc(prev_iter, __func__);

    real_t norma = normaL2Residuo(SL, x, residue(SL, x));

    int iter, result;
    double time = timestamp();
    for (iter = 0; iter < MAXIT && norma > MAXNORMA; iter++){
        result = refine(SL, x);
        if (max_distance(prev_iter, x, SL->n) < SL->erro)
            return iter;
        memcpy(prev_iter, x, sizeof(real_t) * SL->n);
    
        if (result < 0)
            return result;

        norma = normaL2Residuo(SL, x, residue(SL, x));
    }

    *tTotal = timestamp() - time;

    return iter;
}


/*!
  \brief Alocaçao de memória 

  \param n tamanho do SL

  \return ponteiro para SL. NULL se houve erro de alocação
  */
SistLinear_t *alocaSistLinear(unsigned int n)
{
    SistLinear_t *new = (SistLinear_t*) malloc(sizeof(SistLinear_t));
    must_alloc(new, __func__);

    new->n = n;

    // Alocação de matriz contígua
    new->A = (real_t**) malloc(n * sizeof(real_t*));
    must_alloc(new->A, __func__);

    new->A[0] = (real_t*) calloc(n * n, sizeof(real_t));
    must_alloc(new->A[0], __func__);

    for (int i = 0; i < n; i++)
        new->A[i] = new->A[0] + i * n;

    new->b = (real_t*) calloc(n, sizeof(real_t));
    must_alloc(new->b, __func__);

    new->erro = (real_t) 0.0f;

    return new;
}


/*!
  \brief Liberaçao de memória 

  \param sistema linear SL
  */
void liberaSistLinear (SistLinear_t *SL)
{
    free(SL->A[0]);
    free(SL->A);
    free(SL->b);
    free(SL);
}


/*!
  \brief Leitura de SL a partir de Entrada padrão (stdin).

  \return sistema linear SL. NULL se houve erro (leitura ou alocação)
  */
SistLinear_t *lerSistLinear ()
{
    unsigned int n;
    real_t erro;
    fscanf(stdin, "%u\n%e\n", &n, &erro);

    SistLinear_t* SL = alocaSistLinear(n);
    SL->erro = erro;

    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            fscanf(stdin, "%f", &(SL->A[i][j]));

    for (int i = 0; i < n; i++)
        fscanf(stdin, "%f", &(SL->b[i]));

    return SL;
}


// Exibe SL na saída padrão
void prnSistLinear (SistLinear_t *SL)
{
    for (int i = 0; i < SL->n; i++)
        prnVetor(SL->A[i], SL->n);
}


// Exibe um vetor na saída padrão
void prnVetor (real_t *v, unsigned int n)
{
    for (int i = 0; i < n; i++)
        printf("%f ", v[i]);
    printf("\n");
}
