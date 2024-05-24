// SM2加解密demo
#include "LargeNumber.h"
#include "PrimeGenaration.h"
#include "EC_caculation.h"
#include "SM3.h"
#pragma warning(disable:4996)
LN p[ECC_len];   //参数p
LN h[ECC_len];
LN a[ECC_len], b[ECC_len], n[ECC_len];                      //参数a，b，n
ECC_point G;                                                //基点G
LN k[ECC_len];                                              //随机数k

ECC_point Pf[16];
ECC_point Pfe[16];

void KDF(uint8_t* Z, int klen, uint8_t* K)
{
	uint32_t klenv;
	uint32_t ct = 1U;
	int v = 256;
	uint8_t input[1000] = { 0 };
	uint32_t Ha[100][ECC_len] = { 0 };
	strcpy((char*)input, (char*)Z);
	if (klen % v == 0) {
		klenv = klen / v;
	}
	else {
		klenv = klen / v + 1;
	}
	for (uint32_t ct = 1; ct <= klenv; ct++) {
		input[strlen((char*)Z)] = (ct&0xff000000)>>24;
		input[strlen((char*)Z) + 1] = (ct & 0x00ff0000) >> 16;
		input[strlen((char*)Z) + 2] = (ct & 0x0000ff00) >> 8;
		input[strlen((char*)Z) + 3] = (ct & 0x000000ff);
		SM3_encryption(input, strlen((char*)input) + 4, Ha[ct]);
		if (ct == klenv) {
			if (klen % v != 0) {
				int remain = klen % LN_BITS;
				Ha[ct][7 - klen / LN_BITS] = Ha[ct][7 - klen / LN_BITS] & (LN_MAX - (1 << (32 - remain)) + 1);
				ln_r_shift(Ha[ct], Ha[ct], 32 - remain, ECC_len);
				for (int i = 0; i < klen / LN_BITS + 1; i++) {
					Ha[ct][i] = Ha[ct][7 - klen / LN_BITS + i];
				}
				for (int i = klen / LN_BITS + 1; i < ECC_len; i++) {
					Ha[ct][i] = 0;
				}
			}
		}
		ln_to_str(K + 32 * (ct - 1), Ha[ct], ECC_len, 4 * ln_get_digits(Ha[ct], ECC_len));
	}
}

void SM2_encryption(uint8_t* M, ECC_point Pb, ECC_point *C1, LN* C2, LN* C3)
{
	printf("素数p:\n");
	print_ln(p, ECC_len);
	printf("系数a:\n");
	print_ln(a, ECC_len);
	printf("系数b:\n");
	print_ln(b, ECC_len);
	printf("坐标xG:\n");
	print_ln(G.x, ECC_len);
	printf("坐标yG:\n");
	print_ln(G.y, ECC_len);

	//ECC_point Pb;   //公钥  
	//Point_NAF_Mul(G, &Pb, d, ECC_len);
	printf("公钥:\n");
	print_ln(Pb.x, ECC_len);
	print_ln(Pb.y, ECC_len);

	uint8_t block[68];  //随机数k
	initialize_rand();
step:
	generate_rand(block, 68);
	uint32_t len = strlen((const char*)block);
	str_to_ln(block, k, ECC_len, len);
	ln_mod(k, k, ECC_len, n, ECC_len);
	printf("生成随机数k:\n");
	print_ln(k, ECC_len);
	//k[0] = 0x49dd7b4f, k[1] = 0x18e5388d, k[2] = 0x5546d490, k[3] = 0x8afa1742, k[4] = 0x3d957514, k[5] = 0x5b92fd6c, k[6] = 0x6ecfc2b9, k[7] = 0x4c62eefd;


	LN m[ECC_len];  //明文
	str_to_ln(M, m, ECC_len, strlen((char*)M));
	printf("明文:\n");
	print_ln(m, ECC_len);

	ECC_point X2, S;

	Precomputation(32, 64, G);
	Fixed_base_comb_method(G, C1, k, ECC_len);
	//Point_Mul(G, &C1, k, ECC_len);
	printf("C1坐标:\n");
	print_ln(C1->x, ECC_len);
	print_ln(C1->y, ECC_len);

	Point_NAF_Mul(Pb, &S, h, ECC_len);

	Precomputation(32, 64, Pb);
	//Point_NAF_Mul(Pb, &X2, k, ECC_len);
	//Point_Mix_Mul(Pb, &X2, k, ECC_len);
	//Point_Mul(Pb, &X2, k, ECC_len);
	Fixed_base_comb_method(Pb, &X2, k, ECC_len);

	LN t[ECC_len] = { 0 };
	uint8_t T[ECC_len * 4] = { 0 };
	uint8_t x_y[2 * ECC_len * 4] = { 0 };
	uint8_t x_M_y[10 * ECC_len * 4] = { 0 };
	//t[0] = 0x64491603, t[1] = 0xa379e902, t[2] = 0x71dfad8a, t[3] = 0xdae231b0, t[4] = 0x006e30;
	ln_to_str(x_y, X2.x, ECC_len, 4 * ln_get_digits(X2.x, ECC_len));
	ln_to_str(x_y + 32, X2.y, ECC_len, 4 * ln_get_digits(X2.y, ECC_len));
	KDF(x_y, 152, T);
	str_to_ln(T, t, ECC_len, 20);
	if (ln_is_zero(t, ECC_len)) {
		goto step;
	}

	for (int i = 0; i < ECC_len; i++)
		C2[i] = t[i] ^ m[i];

	printf("C2:\n");
	print_ln(C2, ECC_len);
	print_ln(X2.x, ECC_len);
	print_ln(m, ECC_len);
	print_ln(X2.y, ECC_len);

	ln_to_str(x_M_y, X2.x, ECC_len, 4 * ln_get_digits(X2.x, ECC_len));
	strcat((char*)(x_M_y + 32), (char*)M);
	ln_to_str(x_M_y + 32 + strlen((char*)M), X2.y, ECC_len, 4 * ln_get_digits(X2.y, ECC_len));

	SM3_encryption(x_M_y, strlen((char*)x_M_y), C3);
	printf("C3:\n");
	print_ln(C3, ECC_len);
}

void SM2_decryption(ECC_point C1, LN* C2, LN* C3, LN *d)
{
	LN right[ECC_len] = { 0 }, left[ECC_len] = { 0 };
	LN three[ECC_len];
	assign_one_digit(three, 3, ECC_len);
	ln_mod_exp_mary(right, C1.x, three, ECC_len, p, ECC_len);
	LN temp[ECC_len] = { 0 };
	ln_mod_mul(temp, C1.x, a, p, ECC_len);
	ln_add(temp, temp, b, ECC_len);
	ln_add(right, right, temp, ECC_len);
	ln_mod(right, right, ECC_len, p, ECC_len);
	ln_mod_mul(left, C1.y, C1.y, p, ECC_len);
	if (ln_cmp(right, left, ECC_len) != 0) {
		printf("C1不满足曲线！");
		return;
	}

	ECC_point X2;
	uint8_t x_y[2 * ECC_len * 4] = { 0 };
	uint8_t T[ECC_len * 4] = { 0 };
	LN t[ECC_len] = { 0 };
	Point_NAF_Mul(C1, &X2, d, ECC_len);
	ln_to_str(x_y, X2.x, ECC_len, 4 * ln_get_digits(X2.x, ECC_len));
	ln_to_str(x_y + 32, X2.y, ECC_len, 4 * ln_get_digits(X2.y, ECC_len));
	KDF(x_y, 152, T);
	str_to_ln(T, t, ECC_len, 20);
	if (ln_is_zero(t, ECC_len)) {
		printf("t为0！");
		return;
	}

	LN m_1[ECC_len];
	uint8_t M_1[100] = { 0 };
	for (int i = 0; i < ECC_len; i++)
		m_1[i] = t[i] ^ C2[i];
	ln_to_str(M_1, m_1, ECC_len, 4 * ln_get_digits(m_1, ECC_len));
	LN u[ECC_len] = { 0 };
	uint8_t x_M_1_y[10 * ECC_len * 4] = { 0 };
	ln_to_str(x_M_1_y, X2.x, ECC_len, 4 * ln_get_digits(X2.x, ECC_len));
	strcat((char*)(x_M_1_y + 32), (char*)M_1);
	ln_to_str(x_M_1_y + 32 + strlen((char*)M_1), X2.y, ECC_len, 4 * ln_get_digits(X2.y, ECC_len));

	SM3_encryption(x_M_1_y, strlen((char*)x_M_1_y), u);
	if (ln_cmp(u, C3, ECC_len) == 0) {
		printf("解密得到明文:%s", M_1);
	}
	else {
		printf("解密错误！");
	}
}

int main(void)
{
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

	ln_assign_zero(h, ECC_len);
	h[0] = 1;
	//椭圆曲线参数

	uint8_t M[] = "encryption standard";
	LN d[ECC_len];
	ln_assign_zero(d, ECC_len); //私钥
	d[0] = 0x20bb0da0; d[1] = 0x08ddbc29; d[2] = 0xb89463f2; d[3] = 0x34aa7f7c; d[4] = 0x3fbf3535; d[5] = 0x5e2efe28; d[6] = 0xa00637bd; d[7] = 0x1649ab77;
	ECC_point Pb;	//公钥
	ln_assign_zero(Pb.x, ECC_len);
	ln_assign_zero(Pb.y, ECC_len);
	Pb.x[0] = 0x628a7e0a; Pb.x[1] = 0x49a5cf70; Pb.x[2] = 0x581a0e48; Pb.x[3] = 0x0f7ba07e;
	Pb.x[4] = 0x67be491a; Pb.x[5] = 0xc1488afc; Pb.x[6] = 0xa8f3b508; Pb.x[7] = 0x435b39cc;
	Pb.y[0] = 0xba076a42; Pb.y[1] = 0x99ccf77b; Pb.y[2] = 0xdbadf453; Pb.y[3] = 0x01debb2c;
	Pb.y[4] = 0xc1cdf5fe; Pb.y[5] = 0x4c7895e2; Pb.y[6] = 0xf15feecb; Pb.y[7] = 0x75ddba78;

	ECC_point C1;
	ln_assign_zero(C1.x, ECC_len);
	ln_assign_zero(C1.y, ECC_len);
	LN C2[ECC_len] = { 0 }, C3[ECC_len] = { 0 };
	SM2_encryption(M, Pb, &C1, C2, C3);
	SM2_decryption(C1, C2, C3, d);
	system("pause");
	return 0;
}
