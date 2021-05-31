#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>

#include "utils.h"
#include "SistemasLineares.h"


// Cria uma cópia de um SL
SistLinear_t *copiar_SL(SistLinear_t* SL){
    SistLinear_t *new = alocaSistLinear(SL->n);

    memcpy(new->A[0], SL->A[0], sizeof(real_t) * SL->n * SL->n);
    memcpy(new->b, SL->b, sizeof(real_t) * SL->n);
    new->erro = SL->erro;

    return new;
}


// Retrosubstitui as variáveis do sistema para solucioná-lo
int retrosubs(SistLinear_t* SL){
    for (int i = SL->n - 1; i >= 0; i--){
        if (invalid(SL->A[i][i])){
            fprintf(stderr, "No solution.\n");
            return -1;
        }
        SL->b[i] /= SL->A[i][i];
        SL->A[i][i] = 1.0f;
        for (int j = i - 1; j >= 0; j--){
            SL->b[j] -= SL->A[j][i] * SL->b[i];
            SL->A[j][i] = 0.0f;
        }
    }
    return 0;
}


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
  \param tTotal tempo gasto pelo método

  \return código de erro. 0 em caso de sucesso.
*/
int eliminacaoGauss (SistLinear_t *SL, real_t *x, double *tTotal)
{
    SistLinear_t *clone = copiar_SL(SL);

    int i, k, j;
    double m, tempo = timestamp();
    // Triangula usando o pivoteamento parcial
    for (k = 0; k < clone->n - 1; k++){
        unsigned int max_index = k;
        float max = clone->A[k][k];
        for (i = k + 1; i < clone->n; i++)
            if (clone->A[i][k] > max)
                max_index = i;

        if (max_index != k){ // Se terminou em um índice diferente de onde começou, troca
            real_t* aux1 = clone->A[k];
            clone->A[k] = clone->A[max_index];
            clone->A[max_index] = clone->A[k];

            real_t aux2 = clone->b[k];
            clone->b[k] = clone->b[max_index];
            clone->b[max_index] = clone->b[k];
        }

        for (i = k + 1; i < clone->n; i++){
            if (invalid(clone->A[k][k])){
                fprintf(stderr, "No solution.\n"); //revisar
                return -1;
            }
            m = clone->A[i][k] / clone->A[k][k];
            clone->A[i][k] = 0.0f;
            for (j = k + 1; j < clone->n; j++)
                clone->A[i][j] -= (m * clone->A[k][j]);
            clone->b[i] -= m * clone->b[k];
        }
    }

    if (retrosubs(clone)) // Se falhar
        return -1;

    *tTotal = timestamp() - tempo;
    memcpy(x, clone->b, sizeof(real_t) * clone->n);

    liberaSistLinear(clone);

    return 0;
}


// Compara todos os elementos de dois vetores e ve se a diferença entre eles são muito diferentes (baseado no erro)
int too_different(real_t* prev, real_t* next, unsigned int n, real_t error){
    for (int i = 0; i < n; i++)
        if (fabs(next[i] - prev[i]) > error)
            return 1;
    return 0;
}


// Verifica se as iterações estão convergindo
int converging(real_t* prev, real_t* next, real_t prev_diff, unsigned int n){
    real_t dif, max_dif = fabs(next[0] - prev[0]);
    for (int i = 1; i < n; i++){
        dif = fabs(next[i] - prev[i]);
        if (dif > max_dif)
            max_dif = dif;
    }

    if (max_dif <= prev_diff)
        return 1;
        
    fprintf(stderr, "Not converging.\n");
    return 0;
}


/*!
  \brief Método de Jacobi

  \param SL Ponteiro para o sistema linear
  \param x ponteiro para o vetor solução. Ao iniciar função contém
            valor inicial
  \param tTotal tempo gasto pelo método

  \return código de erro. Um nr positivo indica sucesso e o nr
          de iterações realizadas. Um nr. negativo indica um erro:
          -1 (não converge) -2 (sem solução)
*/
int gaussJacobi (SistLinear_t *SL, real_t *x, double *tTotal)
{
    real_t* curr_iter = calloc(SL->n, sizeof(real_t)); // Valores usados na atual iteração
    must_alloc(curr_iter, __func__);

    real_t* next_iter = malloc(SL->n * sizeof(real_t)); // Valores a serem usados na próxima iteração
    must_alloc(next_iter, __func__);
    for (int i = 0; i < SL->n; i++) next_iter[i] = FLT_MAX; // Inicia o vetor com valores muito diferentes da primeira iteração

    real_t* aux_iter = calloc(SL->n, sizeof(real_t)); // Vetor auxiliar para conservar os valores a serem comparados entre iterações
    must_alloc(aux_iter, __func__);

    real_t prev_diff = FLT_MAX; // A maior diferença entre as iterações

    int iter, i, k;
    double sum, tempo = timestamp();
    // Enquanto forem muito diferentes e iter não ultrapassou o limite de iterações, itera
    for (iter = 0; iter < MAXIT && too_different(curr_iter, next_iter, SL->n, SL->erro); iter++){
        memcpy(curr_iter, aux_iter, sizeof(real_t) * SL->n);

        for (i = 0; i < SL->n; i++){
            sum = 0.0f;
            for (k = 0; k < SL->n; k++)
                sum += SL->A[i][k] * curr_iter[k];
            sum -= SL->A[i][i] * curr_iter[i];

            if (invalid(SL->A[i][i])){
                fprintf(stderr, "No solution.\n");
                return -2;
            }
            next_iter[i] = (SL->b[i] - sum) / SL->A[i][i];
        }

        if (!converging(curr_iter, next_iter, prev_diff, SL->n))
            return -1;

        memcpy(aux_iter, next_iter, sizeof(real_t) * SL->n);
    }

    *tTotal = timestamp() - tempo;
    memcpy(x, next_iter, sizeof(real_t) * SL->n);
    
    free(curr_iter);
    free(next_iter);
    free(aux_iter);

    return iter;
}


/*!
  \brief Método de Gauss-Seidel

  \param SL Ponteiro para o sistema linear
  \param x ponteiro para o vetor solução. Ao iniciar função contém
            valor inicial
  \param tTotal tempo gasto pelo método

  \return código de erro. Um nr positivo indica sucesso e o nr
          de iterações realizadas. Um nr. negativo indica um erro:
          -1 (não converge) -2 (sem solução)
  */
int gaussSeidel (SistLinear_t *SL, real_t *x, double *tTotal)
{
    real_t* prev_iter = malloc(SL->n * sizeof(real_t)); // Valores da iteração anterior
    must_alloc(prev_iter, __func__);
    for (int i = 0; i < SL->n; i++) prev_iter[i] = FLT_MAX; // Inicia o vetor com valores muito diferentes da primeira iteração

    real_t* curr_iter = calloc(SL->n, sizeof(real_t)); // Valores da iteração atual
    must_alloc(curr_iter, __func__);

    real_t prev_diff = FLT_MAX; // A maior diferença entre as iterações

    int iter, i, k;
    double sum, tempo = timestamp();
    // Enquanto forem muito diferentes e iter não ultrapassou o limite de iterações, itera
    for (iter = 0; iter < MAXIT && too_different(prev_iter, curr_iter, SL->n, SL->erro); iter++){
        memcpy(prev_iter, curr_iter, sizeof(real_t) * SL->n);
        for (i = 0; i < SL->n; i++){
            sum = 0.0f;
            for (k = 0; k < SL->n; k++)
                sum += SL->A[i][k] * curr_iter[k];
            sum -= SL->A[i][i] * curr_iter[i];
            if (invalid(SL->A[i][i])){
                fprintf(stderr, "No solution.\n");
                return -2;
            }
            curr_iter[i] = (SL->b[i] - sum) / SL->A[i][i];
        }

        if (!converging(prev_iter, curr_iter, prev_diff, SL->n))
            return -1;
    }

    *tTotal = timestamp() - tempo;

    // Escreva a solução em x
    memcpy(x, curr_iter, sizeof(real_t) * SL->n);

    free(prev_iter);
    free(curr_iter);

    return iter;
}


real_t distancia_max(real_t *a, real_t *b, unsigned int n){
    real_t diff, max = fabs(a[0] - b[0]);
    for (int i = 1; i < n; i++){
        diff = fabs(a[i] - b[i]);
        if (diff > max)
            max = diff;
    }
    return max;
}


int refinar(SistLinear_t *SL, real_t *x){
    SistLinear_t* clone = copiar_SL(SL);

    double time;
    
    real_t *res = residuo(SL, x);
    prnVetor(res, SL->n);

    real_t *w = malloc(SL->n * sizeof(real_t));
    must_alloc(w, __func__);
    
    memcpy(clone->b, res, sizeof(real_t) * SL->n);

    int result = eliminacaoGauss(clone, w, &time);
    
    if (result >= 0) // Caso tenha dado tudo certo
        for (int i = 0; i < SL->n; i++)
            x[i] = x[i] + w[i];

    free(res);
    free(w);

    return result;
}


/*!
  \brief Método de Refinamento

  \param SL Ponteiro para o sistema linear
  \param x ponteiro para o vetor solução. Ao iniciar função contém
            valor inicial para início do refinamento
  \param tTotal tempo gasto pelo método

  \return código de erro. Um nr positivo indica sucesso e o nr
          de iterações realizadas. Um nr. negativo indica um erro:
          -1 (não converge) -2 (sem solução)
  */
int refinamento (SistLinear_t *SL, real_t *x, double *tTotal)
{
    int iter;

    for (iter = 0; iter < 10; iter++){
        refinar(SL, x);
    }

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

    new->erro = (real_t) 0.0;

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
    for (int i = 0; i < SL->n; i++){
        for (int j = 0; j < SL->n; j++)
            printf("%f ", SL->A[i][j]);
        printf("\n");
    }
}


// Exibe um vetor na saída padrão
void prnVetor (real_t *v, unsigned int n)
{
    for (int i = 0; i < n; i++)
        printf("%f ", v[i]);
    printf("\n");
}
