// ElGamal加解密demo
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

    LN m[ECC_len];  //明文
    assign_one_digit(m, 965181398, ECC_len);
    printf("明文:\n");
    print_ln(m, ECC_len);

    ECC_point X1, X2;
    LN C[2*ECC_len];
    ln_assign_zero(C, 2 * ECC_len);

    clock_t start, finish;
    Precomputation(32, 64, G);
    Fixed_base_comb_method(G, &X1, k, ECC_len);
    start = clock();
    for (int i = 0; i < 100; i++) {
        Point_Mul(G, &X1, k, ECC_len);
        //Point_Mix_Mul(G, &X1, k, ECC_len);
        //Point_NAF_Mul(G, &X1, k, ECC_len);
        //Fixed_base_comb_method(G, &X1, k, ECC_len);
    }
    finish = clock();
    duration = finish - start;
    printf("点乘耗时%f ms\n", duration);


    Precomputation(32, 64, Pb);

    QueryPerformanceCounter(&startcnt);
    //Point_NAF_Mul(Pb, &X2, k, ECC_len);
    //Point_Mix_Mul(Pb, &X2, k, ECC_len);
    //Point_Mul(Pb, &X2, k, ECC_len);
    Fixed_base_comb_method(Pb, &X2, k, ECC_len);
    //print_ln(X2.x, ECC_len);
    //print_ln(X2.y, ECC_len);
    ln_mod_mul(C, m, X2.x, n, ECC_len);
    QueryPerformanceCounter(&endcnt);
    duration = 1000000 * ((double)endcnt.QuadPart - (double)startcnt.QuadPart) / (double)fre.QuadPart;

    printf("加密耗时%f us\n", duration);
    printf("加密得到密文:\n");
    print_ln(C, ECC_len);


    ECC_point X_2;

    Precomputation(32, 64, X1);
    //Point_Mix_Mul(X1, &X_2, d, ECC_len);
    //Point_NAF_Mul(X1, &X_2, d, ECC_len);
    QueryPerformanceCounter(&startcnt);
    Fixed_base_comb_method(X1, &X_2, d, ECC_len);
    //print_ln(X_2.x, ECC_len);
    //print_ln(X_2.y, ECC_len);
    LN x2_inv[ECC_len];
    LN res[2*ECC_len];
    ln_assign_zero(res, 2 * ECC_len);
    mod_inv(x2_inv, X_2.x, n, ECC_len);
    ln_mod_mul(res, x2_inv, C, n, ECC_len);
    QueryPerformanceCounter(&endcnt);
    duration = 1000000 * ((double)endcnt.QuadPart - (double)startcnt.QuadPart) / (double)fre.QuadPart;
    printf("解密耗时%f us\n", duration);
    printf("解密得到明文:\n");
    print_ln(res, ECC_len);
    system("pause");
    return 0;
}