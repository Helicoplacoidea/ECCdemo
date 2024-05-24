// ECDSA签名及验签demo
#include <math.h>
#include <iostream>
#include "LargeNumber.h"
#include "PrimeGenaration.h"
#include "EC_caculation.h"

LN p[ECC_len], p192[ECC_len], p64[ECC_len], one[ECC_len];   //参数p
LN a[ECC_len], b[ECC_len], n[ECC_len];                      //参数a，b，n
ECC_point G;                                                //基点G
LN k[ECC_len];                                              //随机数k

ECC_point Pf[16];
ECC_point Pfe[16];


int main()
{
    double duration;
    LARGE_INTEGER fre = { 0 };
    LARGE_INTEGER startcnt = { 0 };
    LARGE_INTEGER endcnt = { 0 };

    QueryPerformanceFrequency(&fre);

    p[0] = 0x08f1dfc3;
    p[1] = 0x722edb8b;
    p[2] = 0x5c45517d;
    p[3] = 0x45728391;
    p[4] = 0xbf6ff7de;
    p[5] = 0xe8b92435;
    p[6] = 0x4c044f18;
    p[7] = 0x8542d69e;

    a[0] = 0x3937e498; a[1] = 0xec65228b; a[2] = 0x6831d7e0; a[3] = 0x2f3c848b; a[4] = 0x73bbfeff; a[5] = 0x2417842e; a[6] = 0xfa32c3fd; a[7] = 0x787968b4;

    b[0] = 0x27c5249a; b[1] = 0x6e12d1da; b[2] = 0xb16ba06e; b[3] = 0xf61d59a5; b[4] = 0x484bfe48; b[5] = 0x9cf84241; b[6] = 0xb23b0c84; b[7] = 0x63e4c6d3;

    n[0] = 0xc32e79b7; n[1] = 0x5ae74ee7; n[2] = 0x0485628d; n[3] = 0x29772063; n[4] = 0xbf6ff7dd; n[5] = 0xe8b92435; n[6] = 0x4c044f18; n[7] = 0x8542d69e;

    G.x[0] = 0x7fedd43d; G.x[1] = 0x4c4e6c14; G.x[2] = 0xadd50bdc; G.x[3] = 0x32220b3b;
    G.x[4] = 0xc3cc315e; G.x[5] = 0x746434eb; G.x[6] = 0x1b62eab6; G.x[7] = 0x421debd6;
    G.y[0] = 0xe46e09a2; G.y[1] = 0xa85841b9; G.y[2] = 0xbfa36ea1; G.y[3] = 0xe5d7fdfc;
    G.y[4] = 0x153b70c4; G.y[5] = 0xd47349d2; G.y[6] = 0xcbb42c07; G.y[7] = 0x0680512b;

    LN d[ECC_len];
    ln_assign_zero(d, ECC_len); //私钥
    d[0] = 0x20bb0da0; d[1] = 0x08ddbc29; d[2] = 0xb89463f2; d[3] = 0x34aa7f7c; d[4] = 0x3fbf3535; d[5] = 0x5e2efe28; d[6] = 0xa00637bd; d[7] = 0x1649ab77;
    //assign_one_digit(d, 7, ECC_len);
    printf("私钥:\n");
    print_ln(d, ECC_len);
    printf("\n");

    ECC_point Pb;   //公钥  
    //Point_Mul(G, &Pb, d, ECC_len);
    //Point_Mix_Mul(G, &Pb, d, ECC_len);
    Point_NAF_Mul(G, &Pb, d, ECC_len);
    printf("公钥:\n");
    print_ln(Pb.x, ECC_len);
    print_ln(Pb.y, ECC_len);

    uint8_t block[68];  //随机数k
    initialize_rand();
    generate_rand(block, 68);
    uint32_t len = strlen((const char*)block);
    str_to_ln(block, k, ECC_len, len);
    ln_mod(k, k, ECC_len, n, ECC_len);
    printf("生成随机数k:\n");
    print_ln(k, ECC_len);

    LN m[ECC_len];
    uint8_t M[ECC_len * 4] = { 0 };  //明文
    assign_one_digit(m, 20020808, ECC_len);
    printf("明文:\n");
    print_ln(m, ECC_len);

    ECC_point X1;
    LN C[2 * ECC_len];
    ln_assign_zero(C, 2 * ECC_len);

    Precomputation(32, 64, G);
    Fixed_base_comb_method(G, &X1, k, ECC_len);

    QueryPerformanceCounter(&startcnt);
    LN r[ECC_len], s[ECC_len], k_inv[ECC_len], H[ECC_len] = { 0 };  //签名
    ln_assign(r, X1.x, ECC_len);
    cout << "r:";
    print_ln(r, ECC_len);
    std::hash<char*> ptrHash;
    ln_to_str(M, m, ECC_len, 4 * ln_get_digits(m, ECC_len));
    cout << "hash:" << ptrHash((char*)M) << endl;
    H[0] = ptrHash((char*)M) & LN_MAX;
    H[1] = (ptrHash((char*)M) >> LN_BITS) & LN_MAX;
    LN temp[ECC_len];   //d*r
    ln_mod_mul(temp, d, r, n, ECC_len);
    sub_mod_n(s, H, temp, ECC_len);
    mod_inv(k_inv, k, n, ECC_len);
    ln_mod_mul(s, s, k_inv, n, ECC_len);

    QueryPerformanceCounter(&endcnt);
    duration = 1000000 * ((double)endcnt.QuadPart - (double)startcnt.QuadPart) / (double)fre.QuadPart;

    printf("签名耗时%f us\n", duration);


    QueryPerformanceCounter(&startcnt);
    ECC_point hashG, rQ, tempU, U;

    LN s_inv[ECC_len];
    mod_inv(s_inv, s, n, ECC_len);
    Point_NAF_Mul(G, &hashG, H, ECC_len);
    Point_NAF_Mul(Pb, &rQ, r, ECC_len);
    ln_sub(rQ.y, p, rQ.y, ECC_len);
    Point_Add(hashG, rQ, &tempU, ECC_len);
    Point_NAF_Mul(tempU, &U, s_inv, ECC_len);
    cout << "xU:";
    print_ln(U.x, ECC_len);
    QueryPerformanceCounter(&endcnt);
    duration = 1000000 * ((double)endcnt.QuadPart - (double)startcnt.QuadPart) / (double)fre.QuadPart;
    printf("验签耗时%f us\n", duration);
    system("pause");
    return 0;
}