#include "utils.h"
#include <stdio.h>
#include <math.h>
#include <float.h>

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


// Checa se um número é inoperável
int invalid(real_t num){
    if (num == 0.0f)        return 1;
    if (num == NAN)         return 1;
    if (num == INFINITY)    return 1;
    if (num == -INFINITY)   return 1;
    return 0;
}


real_t *residuo(SistLinear_t *SL, real_t *x){
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
