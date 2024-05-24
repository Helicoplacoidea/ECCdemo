#pragma once
#include "LargeNumber.h"

#define ECC_len 17

typedef struct {
    LN x[ECC_len];
    LN y[ECC_len];
} ECC_point;

typedef struct {
    LN x[ECC_len];
    LN y[ECC_len];
    LN z[ECC_len];
} ECC_Jaco_point;

void sub_mod(LN* res, LN* a, LN* b, uint32_t digits);                                       // res = (a - b) % p
void sub_mod_n(LN* res, LN* a, LN* b, uint32_t digits);                                     // res = (a - b) % n
void Point_neg(ECC_point P, ECC_point* P_neg, uint32_t digits);                             // P_neg = -P
void Point_Add(ECC_point P, ECC_point Q, ECC_point* R, uint32_t digits);                    // R = P + Q
void Point_Self_Add(ECC_point P, ECC_point* R, uint32_t digits);                            // R = P + P
int w_NAF(LN* k, uint32_t w, int ki[]);                                                     // 将k化为NAF表示
void Point_Mul(ECC_point P, ECC_point* R, LN* k, uint32_t digits);                          // R = kP
void Point_Mix_Add(ECC_Jaco_point Q, ECC_point P, ECC_Jaco_point* temp, uint32_t digits);   // temp(Jaco) = Q(Jaco) + P
void Point_Jaco_Self_Add(ECC_Jaco_point Q, ECC_Jaco_point* temp, uint32_t digits);          // temp(Jaco) = Q(temp) + Q(temp)
void Point_Mix_Mul(ECC_point P, ECC_point* R, LN* k, uint32_t digits);                      // R = kP(Mix)
void Point_NAF_Mul(ECC_point P, ECC_point* R, LN* k, uint32_t digits);                      // R = kP

void Precomputation(int e, int d, ECC_point P);
void Fixed_base_comb_method(ECC_point P, ECC_point* R, LN* k, uint32_t digits);             // R = kP