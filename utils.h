#ifndef __UTILS_H__
#define __UTILS_H__

#include <stdlib.h>
#include <sys/time.h>
#include "SistemasLineares.h"

/*  Retorna tempo em milisegundos

    Forma de uso:
 
    double tempo;
    tempo = timestamp();
    <trecho de programa do qual se deseja medir tempo>
    tempo = timestamp() - tempo;
*/

#define MAXNORMA 5.0f

double timestamp(void);

// Certifica que a memória foi alocada
void must_alloc(void *ptr, const char *desc);

// Libera toda a memória apontada por ptrs, inclusive ptrs 
void free_these(void **ptrs, unsigned int n);

// Cria uma cópia de um SistLinear_t
SistLinear_t *copiar_SL(SistLinear_t* SL);

// Retrosubstitui as variáveis do sistema para solucioná-lo
int retrosubs(SistLinear_t* SL);

// Verifica se um número é invalido
int invalid(real_t num);

// Calcula o resíduo de um sistema linear e sua solução
real_t *residue(SistLinear_t *SL, real_t *x);

// Verifica se o sistema converge utilizando iteações de Gauss-Jacobi
int jacobi_converge(SistLinear_t *SL);

// Verifica se o sistema converge utilizando iteações de Gauss-Seidel
int seidel_converge(SistLinear_t *SL);

// Retorna a distância máxima entre os elementos de um vetor
real_t max_distance(real_t *a, real_t *b, unsigned int n);

// Refina uma resultado
int refine(SistLinear_t *SL, real_t *x);

// Verifica se as duas soluções são muito diferentes
int too_different(real_t* prev, real_t* curr, unsigned int n, real_t error);

#endif // __UTILS_H__