#include "EC_caculation.h"

extern LN p[ECC_len], p192[ECC_len], p64[ECC_len], one[ECC_len];
extern LN a[ECC_len], b[ECC_len], n[ECC_len];

extern ECC_point Pf[16];
extern ECC_point Pfe[16];

int aw[16][4] =
{ {0,0,0,0},{0,0,0,1},{0,0,1,0},{0,0,1,1},
{0,1,0,0},{0,1,0,1},{0,1,1,0},{0,1,1,1},
{1,0,0,0},{1,0,0,1},{1,0,1,0},{1,0,1,1},
{1,1,0,0},{1,1,0,1},{1,1,1,0},{1,1,1,1} };

void sub_mod(LN* res, LN* a, LN* b, uint32_t digits)
{
    LN temp[4 * ECC_len];
    if (ln_cmp(a, b, digits) >= 0) {
        ln_sub(temp, a, b, digits);
        ln_mod(res, temp, digits, p, digits);
    }
    else {
        ln_sub(temp, b, a, digits);
        ln_mod(res, temp, digits, p, digits);
        ln_sub(res, p, res, digits);
    }
}

void sub_mod_n(LN* res, LN* a, LN* b, uint32_t digits)
{
    LN temp[4 * ECC_len];
    if (ln_cmp(a, b, digits) >= 0) {
        ln_sub(temp, a, b, digits);
        ln_mod(res, temp, digits, n, digits);
    }
    else {
        ln_sub(temp, b, a, digits);
        ln_mod(res, temp, digits, n, digits);
        ln_sub(res, n, res, digits);
    }
}

void Point_neg(ECC_point P, ECC_point* P_neg, uint32_t digits)
{
    LN zero[ECC_len];
    ln_assign_zero(zero, digits);
    sub_mod(P_neg->y, zero, P.y, digits);
    ln_assign(P_neg->x, P.x, digits);
}

void Point_Add(ECC_point P, ECC_point Q, ECC_point* R, uint32_t digits)   // P != Q
{
    LN lamda[ECC_len], y_sub[ECC_len], x_sub[ECC_len], x_inv[ECC_len];
    sub_mod(y_sub, Q.y, P.y, digits);
    sub_mod(x_sub, Q.x, P.x, digits);
    mod_inv(x_inv, x_sub, p, digits);
    ln_mod_mul(lamda, y_sub, x_inv, p, digits);
    //print_ln(lamda, digits);
    LN x_sum[ECC_len], lamda2[2 * ECC_len], x_sub2[2 * ECC_len];
    ln_add(x_sum, Q.x, P.x, digits);
    ln_mul(lamda2, lamda, lamda, digits);
    sub_mod(R->x, lamda2, x_sum, digits);
    sub_mod(x_sub2, P.x, R->x, digits);
    ln_mul(x_sub2, lamda, x_sub2, digits);
    sub_mod(R->y, x_sub2, P.y, digits);
}

void Point_Self_Add(ECC_point P, ECC_point* R, uint32_t digits)
{
    LN lamda[ECC_len], x_temp[2 * ECC_len], y_temp[ECC_len], y_inv[ECC_len], three[ECC_len];
    ln_mul(x_temp, P.x, P.x, digits);
    assign_one_digit(three, 3, digits);
    ln_mul(x_temp, x_temp, three, digits);
    ln_add(x_temp, x_temp, a, digits);
    ln_add(y_temp, P.y, P.y, digits);
    mod_inv(y_inv, y_temp, p, digits);
    ln_mod_mul(lamda, x_temp, y_inv, p, digits);
    //print_ln(lamda, digits);
    LN x_sum[ECC_len], lamda2[2 * ECC_len], x_sub2[2 * ECC_len];
    ln_add(x_sum, P.x, P.x, digits);
    ln_mul(lamda2, lamda, lamda, digits);
    sub_mod(R->x, lamda2, x_sum, digits);
    sub_mod(x_sub2, P.x, R->x, digits);
    ln_mul(x_sub2, lamda, x_sub2, digits);
    sub_mod(R->y, x_sub2, P.y, digits);
}

int w_NAF(LN* k, uint32_t w, int ki[])
{
    LN one[ECC_len], K[ECC_len];
    assign_one_digit(one, 1, ECC_len);
    ln_assign(K, k, ECC_len);
    int w_exp = pow(2, (double)w);
    int half = w_exp / 2;
    int i = 0;
    LN temp[ECC_len];
    ln_assign_zero(temp, ECC_len);
    while (ln_cmp(K, one, ECC_len) >= 0) {
        if (K[0] % 2 == 1) {
            ki[i] = K[0] & (w_exp - 1);
            if (ki[i] >= half) {
                ki[i] = ki[i] - w_exp;
                assign_one_digit(temp, -ki[i], ECC_len);
                ln_add(K, temp, K, ECC_len);
            }
            else {
                assign_one_digit(temp, ki[i], ECC_len);
                ln_sub(K, K, temp, ECC_len);
            }
        }
        else {
            ki[i] = 0;
        }
        ln_r_shift(K, K, 1, ECC_len);
        i++;
    }
    //for (; i > 0; i--) {
    //    cout << ki[i-1] << " ";
    //}
    //cout << endl;
    return i - 1;
}

void Point_Mul(ECC_point P, ECC_point* R, LN* k, uint32_t digits)
{
    LN E[ECC_len];
    int m = ln_get_bits(k, digits);
    int num, cnt, esh, ei;
    ECC_point Q, temp;

    for (int i = m - 1; i >= 0; i--) {
        ln_assign(E, k, ECC_len);
        num = i / LN_BITS;
        cnt = i - num * LN_BITS;
        esh = E[num] >> (cnt);
        ei = esh & 1U;
        if (i == m - 1) {
            ln_assign(Q.x, P.x, digits);
            ln_assign(Q.y, P.y, digits);
        }
        else {
            Point_Self_Add(Q, &temp, digits);
            ln_assign(Q.x, temp.x, digits);
            ln_assign(Q.y, temp.y, digits);
        }
        if (ei == 1 && (i != (m - 1))) {

            Point_Add(Q, P, &temp, digits);
            ln_assign(Q.x, temp.x, digits);
            ln_assign(Q.y, temp.y, digits);
        }
    }
    ln_assign(R->x, Q.x, digits);
    ln_assign(R->y, Q.y, digits);
    //print_ln(R->x, digits);
    //print_ln(R->y, digits);
}

void Point_Mix_Add(ECC_Jaco_point Q, ECC_point P, ECC_Jaco_point* temp, uint32_t digits)
{

    LN s1[4 * ECC_len], s2[4 * ECC_len], s3[4 * ECC_len], s4[4 * ECC_len], s5[4 * ECC_len], s6[4 * ECC_len];
    LN s7[4 * ECC_len], s8[4 * ECC_len], s9[4 * ECC_len];

    LN tmp1[4 * ECC_len], tmp2[4 * ECC_len], tmp3[4 * ECC_len];

    ln_assign_zero(s1, 4 * digits);
    ln_assign_zero(s2, 4 * digits);
    ln_assign_zero(s3, 4 * digits);
    ln_assign_zero(s4, 4 * digits);
    ln_assign_zero(s5, 4 * digits);
    ln_assign_zero(s6, 4 * digits);
    ln_assign_zero(s7, 4 * digits);
    ln_assign_zero(s8, 4 * digits);
    ln_assign_zero(s9, 4 * digits);
    ln_assign_zero(tmp1, 4 * digits);
    ln_assign_zero(tmp2, 4 * digits);
    ln_assign_zero(tmp3, 4 * digits);

    ln_assign(s1, Q.x, digits);

    ln_mul(s2, Q.z, Q.z, digits);
    ln_mod_mul(s2, s2, P.x, p, digits);

    sub_mod(s3, s1, s2, digits);           //-C

    ln_assign(s4, Q.y, digits);

    ln_mul(tmp1, Q.z, Q.z, digits);
    ln_mul(tmp2, Q.z, P.y, digits);
    ln_mod_mul(s5, tmp1, tmp2, p, digits);

    sub_mod(s6, s4, s5, digits);            //-D

    ln_add(s7, s1, s2, digits);
    ln_mod(s7, s7, digits, p, digits);

    ln_add(s8, s4, s5, digits);
    ln_mod(s8, s8, digits, p, digits);

    ln_assign_zero(tmp1, 4 * digits);
    ln_assign_zero(tmp2, 4 * digits);
    ln_mod_mul(tmp1, s6, s6, p, digits);
    ln_mul(tmp2, s3, s3, digits);
    ln_mod_mul(tmp2, tmp2, s7, p, digits);
    sub_mod(temp->x, tmp1, tmp2, digits);

    //ln_assign_zero(tmp1, 4 * digits);
    //ln_assign_zero(tmp2, 4 * digits);
    //ln_mul(tmp1, s3, s3, digits);
    //ln_mod_mul(tmp1, tmp1, Q.x, p, digits);
    //sub_mod(tmp2, tmp1, temp->x, digits);
    //ln_mod_mul(tmp3, tmp2, s3, p, digits);
    //ln_sub(tmp3, p, tmp3, digits);
    //ln_assign_zero(tmp1, 4 * digits);
    //ln_assign_zero(tmp2, 4 * digits);
    //ln_mul(tmp1, s3, s3, digits);
    //ln_mod_mul(tmp1, tmp1, s3, p, digits);
    //ln_mod_mul(tmp2, tmp1, Q.y, p, digits);
    //sub_mod(temp->y, tmp2, tmp3, digits);

    ln_assign_zero(tmp1, 4 * digits);
    ln_assign_zero(tmp2, 4 * digits);
    ln_mul(tmp1, s3, s3, digits);
    ln_mod_mul(tmp1, tmp1, s7, p, digits);
    ln_add(tmp2, temp->x, temp->x, digits);
    sub_mod(s9, tmp1, tmp2, digits);

    ln_assign_zero(tmp1, 4 * digits);
    ln_assign_zero(tmp2, 4 * digits);
    ln_mod_mul(tmp1, s9, s6, p, digits);
    ln_mul(tmp2, s3, s3, digits);
    ln_mod_mul(tmp2, tmp2, s3, p, digits);
    ln_mod_mul(tmp2, tmp2, s8, p, digits);
    sub_mod(tmp3, tmp1, tmp2, 2 * digits);
    LN two[ECC_len], two_inv[ECC_len];
    assign_one_digit(two, 2, digits);
    mod_inv(two_inv, two, p, digits);
    ln_mod_mul(temp->y, tmp3, two_inv, p, digits);

    ln_mod_mul(temp->z, Q.z, s3, p, digits);
}

void Point_Jaco_Self_Add(ECC_Jaco_point Q, ECC_Jaco_point* temp, uint32_t digits)
{
    LN A[4 * ECC_len], B[4 * ECC_len], C[4 * ECC_len], D[4 * ECC_len];
    LN two[ECC_len], three[ECC_len], four[ECC_len], eight[ECC_len];
    LN s1[4 * ECC_len], s2[4 * ECC_len], s3[4 * ECC_len];
    ln_assign_zero(A, 4 * digits);
    ln_assign_zero(B, 4 * digits);
    ln_assign_zero(C, 4 * digits);
    ln_assign_zero(D, 4 * digits);
    ln_assign_zero(s1, 4 * digits);
    ln_assign_zero(s2, 4 * digits);
    ln_assign_zero(s3, 4 * digits);

    assign_one_digit(two, 2, digits);
    assign_one_digit(three, 3, digits);
    assign_one_digit(four, 4, digits);
    assign_one_digit(eight, 8, digits);
    ln_mul(s1, Q.x, four, digits);      //4X1
    ln_mul(s2, Q.y, Q.y, digits);       //Y1^2
    ln_mod_mul(A, s1, s2, p, digits);   //4X1Y1^2

    ln_assign_zero(s1, 4 * digits);
    ln_mod_mul(s1, s2, s2, p, digits);         //y1^4
    ln_mod_mul(B, s1, eight, p, digits);

    ln_assign_zero(s1, 4 * digits);
    ln_assign_zero(s2, 4 * digits);
    ln_mul(s1, Q.z, Q.z, digits);
    ln_mod_mul(s1, s1, s1, p, digits);
    ln_mod_mul(s2, a, s1, p, digits);
    ln_mul(s3, Q.x, Q.x, digits);
    ln_mod_mul(s3, three, s3, p, digits);
    ln_add(C, s3, s2, digits);
    ln_mod(C, C, digits, p, digits);

    ln_assign_zero(s1, 4 * digits);
    ln_assign_zero(s2, 4 * digits);
    ln_assign_zero(s3, 4 * digits);
    ln_mul(s1, C, C, digits);
    ln_mul(s2, two, A, digits);
    sub_mod(D, s1, s2, digits);

    ln_assign_zero(s1, 4 * digits);
    ln_assign_zero(s2, 4 * digits);
    ln_assign_zero(s3, 4 * digits);
    ln_assign(temp->x, D, digits);
    sub_mod(s1, A, D, digits);
    ln_mul(s2, s1, C, digits);
    sub_mod(temp->y, s2, B, digits);
    ln_assign_zero(s1, 4 * digits);
    ln_assign_zero(s2, 4 * digits);
    ln_assign_zero(s3, 4 * digits);
    ln_mul(s1, Q.y, Q.z, digits);
    ln_mod_mul(temp->z, s1, two, p, digits);
}

void Point_Mix_Mul(ECC_point P, ECC_point* R, LN* k, uint32_t digits)
{
    LN E[ECC_len];
    int m = ln_get_bits(k, digits);
    int num, cnt, esh, ei;
    ECC_Jaco_point Q, temp;

    //cout << endl;
    for (int i = m - 1; i >= 0; i--) {
        ln_assign(E, k, ECC_len);
        num = i / LN_BITS;
        cnt = i - num * LN_BITS;
        esh = E[num] >> (cnt);
        ei = esh & 1U;
        if (i == m - 1) {
            ln_assign(Q.x, P.x, digits);
            ln_assign(Q.y, P.y, digits);
            assign_one_digit(Q.z, 1, ECC_len);
        }
        else {
            Point_Jaco_Self_Add(Q, &temp, digits);
            ln_assign(Q.x, temp.x, digits);
            ln_assign(Q.y, temp.y, digits);
            ln_assign(Q.z, temp.z, digits);
        }
        if (ei == 1 && (i != (m - 1))) {
            Point_Mix_Add(Q, P, &temp, digits);
            ln_assign(Q.x, temp.x, digits);
            ln_assign(Q.y, temp.y, digits);
            ln_assign(Q.z, temp.z, digits);
        }
    }

    LN Z2[ECC_len], Z3[ECC_len], Z2_inv[ECC_len], Z3_inv[ECC_len];
    ln_mod_mul(Z2, Q.z, Q.z, p, digits);
    ln_mod_mul(Z3, Z2, Q.z, p, digits);
    mod_inv(Z2_inv, Z2, p, digits);
    mod_inv(Z3_inv, Z3, p, digits);
    ln_mod_mul(R->x, Q.x, Z2_inv, p, digits);
    ln_mod_mul(R->y, Q.y, Z3_inv, p, digits);
}

void Point_NAF_Mul(ECC_point P, ECC_point* R, LN* k, uint32_t digits)
{
    LN E[ECC_len];
    ECC_Jaco_point Q, temp;
    ECC_point Pi[8];

    LN three[ECC_len], five[ECC_len], seven[ECC_len];
    assign_one_digit(three, 3, digits);
    assign_one_digit(five, 5, digits);
    assign_one_digit(seven, 7, digits);
    Pi[1] = P;
    Point_Mul(P, &Pi[3], three, digits);
    Point_Mul(P, &Pi[5], five, digits);
    Point_Mul(P, &Pi[7], seven, digits);

    for (int i = 0; i < 8; i += 2) {
        ln_assign(Pi[i].x, Pi[i + 1].x, digits);
        ln_sub(Pi[i].y, p, Pi[i + 1].y, digits);
    }

    int ki[256] = { 0 };
    int len = w_NAF(k, 4, ki);
    int cnt = 0;

    for (int i = len; i >= 0; i--) {
        cnt = ki[i];
        if (i == len) {
            ln_assign(Q.x, Pi[cnt].x, digits);
            ln_assign(Q.y, Pi[cnt].y, digits);
            assign_one_digit(Q.z, 1, ECC_len);
        }
        else {
            Point_Jaco_Self_Add(Q, &temp, digits);
            ln_assign(Q.x, temp.x, digits);
            ln_assign(Q.y, temp.y, digits);
            ln_assign(Q.z, temp.z, digits);
        }
        if ((cnt != 0) && (i != len)) {
            if (cnt > 0) {
                Point_Mix_Add(Q, Pi[cnt], &temp, digits);
                ln_assign(Q.x, temp.x, digits);
                ln_assign(Q.y, temp.y, digits);
                ln_assign(Q.z, temp.z, digits);
            }
            else {
                Point_Mix_Add(Q, Pi[-cnt - 1], &temp, digits);
                ln_assign(Q.x, temp.x, digits);
                ln_assign(Q.y, temp.y, digits);
                ln_assign(Q.z, temp.z, digits);
            }
        }
    }
    if (ln_is_zero(Q.z, ECC_len)) {
        printf("结果为无穷远点，错误！");
        return;
    }
    LN Z2[ECC_len], Z3[ECC_len], Z2_inv[ECC_len], Z3_inv[ECC_len];
    ln_mod_mul(Z2, Q.z, Q.z, p, digits);
    ln_mod_mul(Z3, Z2, Q.z, p, digits);
    mod_inv(Z2_inv, Z2, p, digits);
    mod_inv(Z3_inv, Z3, p, digits);
    ln_mod_mul(R->x, Q.x, Z2_inv, p, digits);
    ln_mod_mul(R->y, Q.y, Z3_inv, p, digits);
}

void Precomputation(int e, int d, ECC_point P)
{
    //int cnt = 0;
    LN CNT[ECC_len];
    LN CNT_2[2 * ECC_len];
    ln_assign_zero(CNT, ECC_len);
    ln_assign_zero(CNT_2, 2 * ECC_len);
    LN exp_2[ECC_len];
    LN exp_2_e[ECC_len];
    ln_assign_zero(exp_2, ECC_len);
    ln_assign_zero(exp_2_e, ECC_len);
    ln_assign_2exp(exp_2_e, e, ECC_len);
    cout << endl;
    for (int i = 0; i < 16; i++) {
        for (int j = 3; j >= 0; j--) {
            if (j == 3) {
                assign_one_digit(exp_2, 1, ECC_len);
            }
            else {
                ln_assign_2exp(exp_2, (3 - j) * d, ECC_len);
            }
            if (aw[i][j] == 1) {
                ln_add(CNT, CNT, exp_2, ECC_len);
            }
            //cnt += pow(2, j * d) * aw[i][j];
        }
        //print_ln(CNT, ECC_len);
        ln_mul(CNT_2, CNT, exp_2_e, ECC_len);
        //print_ln(CNT_2, 2*ECC_len);
        if (ln_is_zero(CNT, ECC_len) == 1) {
            Pf[i].x[0] = -1;
            Pf[i].x[0] = -1;
            Pfe[i].x[0] = -1;
            Pfe[i].x[0] = -1;
        }
        else {
            Point_NAF_Mul(P, &Pf[i], CNT, ECC_len);
            Point_NAF_Mul(P, &Pfe[i], CNT_2, ECC_len);
            //Point_Mul(P, &Pf[i], CNT, ECC_len);
            //Point_Mul(P, &Pfe[i], CNT_2, ECC_len);
            //print_ln(Pf[i].x, ECC_len);
            //print_ln(Pf[i].y, ECC_len);
            //print_ln(Pfe[i].x, ECC_len);
            //print_ln(Pfe[i].y, ECC_len);
        }
        ln_assign_zero(CNT, ECC_len);
        ln_assign_zero(CNT_2, 2 * ECC_len);
        ln_assign_zero(exp_2, ECC_len);
    }
}

void Fixed_base_comb_method(ECC_point P, ECC_point* R, LN* k, uint32_t digits)
{
    int m = ln_get_bits(k, digits);
    int w = 4;
    //printf("m:%d\n", m);
    int d = 0, e = 0;
    if (m % w == 0) {
        d = m / w;
    }
    else {
        d = m / w + 1;
    }
    //printf("d:%d\n", d);
    if (d % 2 == 0) {
        e = d / 2;
    }
    else {
        e = d / 2 + 1;
    }
    ECC_Jaco_point Q, temp;
    ln_assign_zero(Q.x, digits);
    ln_assign_zero(Q.y, digits);
    ln_assign_zero(Q.z, digits);

    LN E[ECC_len];
    ln_assign(E, k, digits);
    int ki[256] = { 0 };
    int num, cnt, esh;
    for (int i = m - 1; i >= 0; i--) {
        num = i / LN_BITS;
        cnt = i - num * LN_BITS;
        esh = E[num] >> (cnt);
        ki[i] = esh & 1U;
        //cout << ki[i] << " ";
    }
    int cnt0 = 0;
    int cnt1 = 0;
    ECC_point P_sum;
    for (int i = e - 1; i >= 0; i--) {
        Point_Jaco_Self_Add(Q, &temp, digits);
        ln_assign(Q.x, temp.x, digits);
        ln_assign(Q.y, temp.y, digits);
        ln_assign(Q.z, temp.z, digits);

        cnt0 = ki[i] + 2 * ki[d + i] + 4 * ki[2 * d + i] + 8 * ki[3 * d + i];
        cnt1 = ki[i + e] + 2 * ki[d + i + e] + 4 * ki[2 * d + i + e] + 8 * ki[3 * d + i + e];

        if ((Pf[cnt0].x[0] == -1) && (Pfe[cnt1].x[0] != -1)) {
            ln_assign(P_sum.x, Pfe[cnt1].x, digits);
            ln_assign(P_sum.y, Pfe[cnt1].y, digits);
        }
        else if ((Pf[cnt0].x[0] != -1) && (Pfe[cnt1].x[0] == -1)) {
            ln_assign(P_sum.x, Pf[cnt0].x, digits);
            ln_assign(P_sum.y, Pf[cnt0].y, digits);
        }
        else if ((Pf[cnt0].x[0] == -1) && (Pfe[cnt1].x[0] == -1)) {
            continue;
        }
        else {
            Point_Add(Pf[cnt0], Pfe[cnt1], &P_sum, digits);
        }

        if (i == e - 1) {
            ln_assign(Q.x, P_sum.x, digits);
            ln_assign(Q.y, P_sum.y, digits);
            assign_one_digit(Q.z, 1, digits);
        }
        else {
            Point_Mix_Add(Q, P_sum, &temp, digits);
            ln_assign(Q.x, temp.x, digits);
            ln_assign(Q.y, temp.y, digits);
            ln_assign(Q.z, temp.z, digits);
        }

    }
    LN Z2[ECC_len], Z3[ECC_len], Z2_inv[ECC_len], Z3_inv[ECC_len];
    ln_mod_mul(Z2, Q.z, Q.z, p, digits);
    ln_mod_mul(Z3, Z2, Q.z, p, digits);
    mod_inv(Z2_inv, Z2, p, digits);
    mod_inv(Z3_inv, Z3, p, digits);
    ln_mod_mul(R->x, Q.x, Z2_inv, p, digits);
    ln_mod_mul(R->y, Q.y, Z3_inv, p, digits);
}