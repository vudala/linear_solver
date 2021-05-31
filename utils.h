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

double timestamp(void);

real_t *residuo(SistLinear_t *SL, real_t *x);

int invalid(real_t num);

void must_alloc(void *ptr, const char *desc);

#endif // __UTILS_H__

