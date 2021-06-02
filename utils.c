#include "utils.h"
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <string.h>
#include <stdlib.h>

/*  Retorna tempo em milisegundos

    Forma de uso:
 
    double tempo;
    tempo = timestamp();
    <trecho de programa do qual se deseja medir tempo>
    tempo = timestamp() - tempo;
*/

double timestamp(void)
{
  struct timeval tp;
  gettimeofday(&tp, NULL);
  return((double)(tp.tv_sec*1000.0 + tp.tv_usec/1000.0));
}


// Certifica que o ponteiro foi alocado
void must_alloc(void *ptr, const char *desc){
    if (!ptr){
        fprintf(stderr, "Memory allocation failure: %s\n", desc);
        exit(-1);
    }
}


// Cria uma cópia de um SL
SistLinear_t *copiar_SL(SistLinear_t* SL){
    SistLinear_t *new = alocaSistLinear(SL->n);

    memcpy(new->A[0], SL->A[0], sizeof(real_t) * SL->n * SL->n);
    memcpy(new->b, SL->b, sizeof(real_t) * SL->n);
    new->erro = SL->erro;

    return new;
}


// Checa se um número é inoperável
int invalid(real_t num){
    if (num == NAN)         return 1;
    if (num == INFINITY)    return 1;
    if (num == -INFINITY)   return 1;
    return 0;
}


// Retrosubstitui as variáveis do sistema para solucioná-lo
int retrosubs(SistLinear_t* SL){
    for (int i = SL->n - 1; i >= 0; i--){
        SL->b[i] /= SL->A[i][i];
        if (invalid(SL->b[i])){
            fprintf(stderr, "Retrosubs floating point failure.\n");
            return -1;
        }
        SL->A[i][i] = 1.0f;
        for (int j = i - 1; j >= 0; j--){
            SL->b[j] -= SL->A[j][i] * SL->b[i];
            if (invalid(SL->b[j])){
                fprintf(stderr, "Retrosubs floating point failure.\n");
                return -1;
            }
            SL->A[j][i] = 0.0f;
        }
    }
    return 0;
}


// Calcula o resíduo de um sistema linear e sua solução
real_t *residue(SistLinear_t *SL, real_t *x){
    real_t *res = malloc(SL->n * sizeof(real_t));
    must_alloc(res, __func__);

    int i, k;
    double ax;

    for (i = 0; i < SL->n; i++){
        ax = 0.0f;
        for (k = 0; k < SL->n; k++)
            ax += SL->A[i][k] * x[k];
        res[i] = SL->b[i] - ax;
    }
    
    return res;
}


void free_these(void **ptrs, unsigned int n){
    for (int i = 0; i < n; i++)
        free(ptrs[i]);
    free(ptrs);
}


int jacobi_converge(SistLinear_t *SL){
    double sum;
    for (int i = 0; i < SL->n; i++){
        sum = 0.0f;
        for (int j = 0; j < SL->n; j++){
            if (i != j)
                sum += SL->A[i][j];
        }
        if (fabs(sum) >= fabs(SL->A[i][i]))
            return 0;
    }
    return 1;
}


int seidel_converge(SistLinear_t *SL){
    double sum;
    real_t *betas = malloc(sizeof(real_t) * SL->n);
    int i, j;
    for (i = 0; i < SL->n; i++){
        sum = 0.0f;
        for (j = 0; j < i; j++)
            sum += betas[j] * SL->A[i][j];

        for (j = i + 1; j < SL->n; j++)
            sum += SL->A[i][j];
        sum = fabs(sum);
        if (sum > SL->A[i][i]){
            free(betas);
            return 0;
        }
        betas[i] = sum / fabs(SL->A[i][i]);
    }
    free(betas);
    return 1;
}


// Compara todos os elementos de dois vetores e ve se a diferença entre eles são muito diferentes (baseado no erro)
int too_different(real_t* prev, real_t* curr, unsigned int n, real_t error)
{
    for (int i = 0; i < n; i++)
        if (fabs(curr[i] - prev[i]) > error)
            return 1;
    return 0;
}


real_t max_distance(real_t *a, real_t *b, unsigned int n)
{
    real_t diff, max = fabs(a[0] - b[0]);
    for (int i = 1; i < n; i++){
        diff = fabs(a[i] - b[i]);
        if (diff > max)
            max = diff;
    }
    return max;
}


int refine(SistLinear_t *SL, real_t *x)
{
    SistLinear_t* clone = copiar_SL(SL);

    real_t *res = residue(SL, x);

    real_t *w = malloc(SL->n * sizeof(real_t));
    must_alloc(w, __func__);
    
    memcpy(clone->b, res, sizeof(real_t) * SL->n);

    double time;
    int result = eliminacaoGauss(clone, w, &time);
    
    if (result >= 0) // Caso tenha dado tudo certo
        for (int i = 0; i < SL->n; i++)
            x[i] = x[i] + w[i];

    free(res);
    free(w);

    return result;
}