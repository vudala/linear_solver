#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "utils.h"
#include "SistemasLineares.h"


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

        printf("***** Sistema %i --> n = %i, erro: %f\n", counter, SL->n, SL->erro);
        fprintf(stderr, "***** Sistema %i --> n = %i, erro: %f\n", counter, SL->n, SL->erro);

        result = eliminacaoGauss(SL, x, &time);
        if (result == 0){
            printf("===> Eliminação de Gauss: %1.10f ms\n--> X: ", time);
            prnVetor(x, SL->n);
            res = residue(SL, x);
            norma = normaL2Residuo(SL, x, res);
            free(res);

            printf("--> Norma L2 do residuo: %f\n\n", norma);
            
            if (norma > MAXNORMA){
                result = refinamento(SL, x, &time);
                res = residue(SL, x);
                norma = normaL2Residuo(SL, x, res);
                free(res);

                if (result >= 0){
                    printf("===> Refinamento: %1.10f ms --> %i iterações\n--> X: ", time, result);
                    prnVetor(x, SL->n);
                    printf("--> Norma L2 do residuo: %f\n\n", norma);
                }
            }
        }

        result = gaussJacobi(SL, x, &time);
        if (result >= 0){
            printf("===> Jacobi: %1.10f ms --> %i iterações\n--> X: ", time, result);
            prnVetor(x, SL->n);

            res = residue(SL, x);
            norma = normaL2Residuo(SL, x, res);
            free(res);

            printf("--> Norma L2 do residuo: %f\n\n", norma);

            if (norma > MAXNORMA){
                result = refinamento(SL, x, &time);
                res = residue(SL, x);
                norma = normaL2Residuo(SL, x, res);
                free(res);

                if (result >= 0){
                    printf("===> Refinamento: %1.10f ms --> %i iterações\n--> X: ", time, result);
                    prnVetor(x, SL->n);
                    printf("--> Norma L2 do residuo: %f\n\n", norma);
                }
            }
        }

        result = gaussSeidel(SL, x, &time);
        if (result >= 0){
            printf("===> Gauss-Seidel: %1.10f ms --> %i iterações\n--> X: ", time, result);
            prnVetor(x, SL->n);

            res = residue(SL, x);
            norma = normaL2Residuo(SL, x, res);
            free(res);

            printf("--> Norma L2 do residuo: %f\n\n", norma);

            if (norma > MAXNORMA){
                result = refinamento(SL, x, &time);
                res = residue(SL, x);
                norma = normaL2Residuo(SL, x, res);
                free(res);

                if (result >= 0){
                    printf("===> Refinamento: %1.10f ms --> %i iterações\n--> X: ", time, result);
                    prnVetor(x, SL->n);
                    printf("--> Norma L2 do residuo: %f\n\n", norma);
                }
            }
        }

        liberaSistLinear(SL);
        free(x);

        counter++;
        getchar(); // Consome o \n
    }
    
    return 0;
}