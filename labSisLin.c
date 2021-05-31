#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "utils.h"
#include "SistemasLineares.h"


#define MAXNORMA 5.0f

int main (){
    int result, counter = 1;
    double time;
    SistLinear_t *SL = NULL;
    real_t *x = NULL;
    real_t *res = NULL;
    double norma;
    
    while (!feof(stdin)){
        SL = lerSistLinear();
        x = malloc(sizeof(real_t) * SL->n);
        res = malloc(sizeof(real_t) * SL->n);
        must_alloc(res, __func__);

        printf("***** Sistema %i --> n = %i, erro: %f\n", counter, SL->n, SL->erro);

        result = eliminacaoGauss(SL, x, &time);
        if (result == 0){
            printf("===> Eliminação de Gauss: %1.10f ms\n", time);
            printf("--> X: ");
            prnVetor(x, SL->n);
            res = residuo(SL, x);
            norma = normaL2Residuo(SL, x, res);
            printf("--> Norma L2 do residuo: %f\n\n", norma);
            if (norma > MAXNORMA)
                printf("\n\n\n REFINAR \n\n\n");
        }
        else 
            printf("%i\n", result);

        result = gaussJacobi(SL, x, &time);
        if (result >= 0){
            printf("===> Jacobi: %1.10f ms --> %i iterações\n", time, result);
            printf("--> X: ");
            prnVetor(x, SL->n);
            res = residuo(SL, x);
            norma = normaL2Residuo(SL, x, res);
            printf("--> Norma L2 do residuo: %f\n\n", norma);

            if (norma > MAXNORMA){
                result = refinamento(SL, x, &time);
                res = residuo(SL, x);
                norma = normaL2Residuo(SL, x, res);

                if (result >= 0){
                    printf("===> Refinamento: %1.10f ms --> %i iterações\n", time, result);
                    printf("--> X: ");
                    prnVetor(x, SL->n);
                    printf("--> Norma L2 do residuo: %f\n\n", norma);
                }
            }

            
        }
        else 
            printf("%i\n", result);

        result = gaussSeidel(SL, x, &time);
        if (result >= 0){
            printf("===> Gauss-Seidel: %1.10f ms --> %i iterações\n", time, result);
            printf("--> X: ");
            prnVetor(x, SL->n);
            res = residuo(SL, x);
            norma = normaL2Residuo(SL, x, res);
            printf("--> Norma L2 do residuo: %f\n\n", norma);

            if (norma > MAXNORMA){
                result = refinamento(SL, x, &time);
                printf("%i\n", result);
                res = residuo(SL, x);
                norma = normaL2Residuo(SL, x, res);

                if (result >= 0){
                    printf("===> Refinamento: %1.10f ms --> %i iterações\n", time, result);
                    printf("--> X: ");
                    prnVetor(x, SL->n);
                    printf("--> Norma L2 do residuo: %f\n\n", norma);
                }
            }
        }
        else 
            printf("%i\n", result);

        liberaSistLinear(SL);
        free(x);
        free(res);

        counter++;
        getchar(); // Consome o \n

    }
    
    return 0;
}