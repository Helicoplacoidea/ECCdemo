// SM2签名及验签demo
#include "LargeNumber.h"
#include "PrimeGenaration.h"
#include "EC_caculation.h"
#include "SM3.h"
#pragma warning(disable:4996)
LN p[ECC_len], p192[ECC_len], p64[ECC_len], one[ECC_len];   //参数p
LN a[ECC_len], b[ECC_len], n[ECC_len];                      //参数a，b，n
ECC_point G;                                                //基点G
LN k[ECC_len];                                              //随机数k

ECC_point Pf[16];
ECC_point Pfe[16];

void SM2_signature(uint8_t* M, LN* Z_A, LN* d, LN* r, LN* s)
{
	uint8_t temp[100] = { 0 };
	ln_to_str(temp, Z_A, ECC_len, 4 * ln_get_digits(Z_A, ECC_len));
	strcat((char*)temp, (char*)M);
	LN e[ECC_len] = { 0 };
	SM3_encryption(temp, strlen((char*)temp), e);
	printf("密码杂凑数值e:\n");
	print_ln(e, ECC_len);

	//k[0] = 0x1fb2f96f, k[1] = 0x260dbaae, k[2] = 0xdd72b727, k[3] = 0xc176d925, k[4] = 0x4817663f, k[5] = 0x94f94e93, k[6] = 0x385c175c, k[7] = 0x6cb28d99;
	uint8_t block[68];  //随机数k
	initialize_rand();
step:
	generate_rand(block, 68);
	uint32_t len = strlen((const char*)block);
	str_to_ln(block, k, ECC_len, len);
	ln_mod(k, k, ECC_len, n, ECC_len);
	printf("生成随机数k:\n");
	print_ln(k, ECC_len);

	ECC_point C1;
	Precomputation(32, 64, G);
	Fixed_base_comb_method(G, &C1, k, ECC_len);
	//Point_Mul(G, &C1, k, ECC_len);
	printf("C1坐标:\n");
	print_ln(C1.x, ECC_len);
	print_ln(C1.y, ECC_len);

	LN d_inv[ECC_len], one[ECC_len], temp1[ECC_len], temp2[ECC_len];
	ln_add(r, e, C1.x, ECC_len);
	ln_mod(r, r, ECC_len, n, ECC_len);
	LN r_k[ECC_len];
	ln_add(r_k, r, k, ECC_len);
	if (ln_is_zero(r, ECC_len) || ln_equal(r_k, n, ECC_len))
		goto step;

	printf("r:\n");
	print_ln(r, ECC_len);

	assign_one_digit(one, 1, ECC_len);
	ln_add(temp1, d, one, ECC_len);
	mod_inv(d_inv, temp1, n, ECC_len);

	ln_mul(temp2, r, d, ECC_len);
	sub_mod_n(temp2, k, temp2, ECC_len);
	ln_mod_mul(s, temp2, d_inv, n, ECC_len);
	if (ln_is_zero(s, ECC_len))
		goto step;
	printf("s:\n");
	print_ln(s, ECC_len);
}

void SM2_verification(uint8_t* M, LN* Z_A, ECC_point Pb, LN* r, LN* s)
{
	if (ln_cmp(r, n, ECC_len) >= 0) {
		printf("签名验证不通过！");
		return;
	}
	if (ln_cmp(s, n, ECC_len) >= 0) {
		printf("签名验证不通过！");
		return;
	}
	uint8_t temp[100] = { 0 };
	ln_to_str(temp, Z_A, ECC_len, 4 * ln_get_digits(Z_A, ECC_len));
	strcat((char*)temp, (char*)M);
	LN e[ECC_len] = { 0 };
	SM3_encryption(temp, strlen((char*)temp), e);
	printf("\n密码杂凑数值e:\n");
	print_ln(e, ECC_len);

	LN t[ECC_len] = { 0 };
	ln_add(t, r, s, ECC_len);
	ln_mod(t, t, ECC_len, n, ECC_len);
	if (ln_is_zero(t, ECC_len)) {
		printf("签名验证不通过！");
		return;
	}
	ECC_point C1, X1, X2;
	Point_NAF_Mul(G, &X1, s, ECC_len);
	Point_NAF_Mul(Pb, &X2, t, ECC_len);
	Point_Add(X1, X2, &C1, ECC_len);

	LN R[ECC_len] = { 0 };
	ln_add(R, e, C1.x, ECC_len);
	ln_mod(R, R, ECC_len, n, ECC_len);
	printf("R:\n");
	print_ln(R, ECC_len);

	if (ln_cmp(r, R, ECC_len) == 0) {
		printf("签名验证通过！");
		return;
	}
	else {
		printf("签名验证不通过！");
		return;
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
	//椭圆曲线参数

	uint8_t ID_A[] = "ALICE123@YAHOO.COM";
	uint8_t M[] = "message digest";
	uint8_t ENTL_A[2] = { 0 };
	ENTL_A[0] = 0x00, ENTL_A[1] = 0x90;

	LN d[ECC_len];
	ln_assign_zero(d, ECC_len); //私钥
	d[0] = 0x15897263; d[1] = 0x0c23661d; d[2] = 0x171b1b65; d[3] = 0x2a519a55; d[4] = 0x3dff7979; d[5] = 0x068c8d80; d[6] = 0xbd433c6c; d[7] = 0x128b2fa8;
	
	printf("阶n:\n");
	print_ln(n, ECC_len);
	printf("私钥d:\n");
	print_ln(d, ECC_len);
	
	ECC_point Pb;	//公钥
	ln_assign_zero(Pb.x, ECC_len);
	ln_assign_zero(Pb.y, ECC_len);
	Pb.x[0] = 0x4df2548a; Pb.x[1] = 0xe97c04ff; Pb.x[2] = 0xa5844495; Pb.x[3] = 0x02bb79e2;
	Pb.x[4] = 0x825be462; Pb.x[5] = 0x471bee11; Pb.x[6] = 0x8aa0f119; Pb.x[7] = 0x0ae4c779;
	Pb.y[0] = 0xb798e857; Pb.y[1] = 0xa9fe0c6b; Pb.y[2] = 0xa176d684; Pb.y[3] = 0x07353e53;
	Pb.y[4] = 0x17b7f16f; Pb.y[5] = 0x6352a73c; Pb.y[6] = 0x8f1cd4e1; Pb.y[7] = 0x7c0240f8;

	uint8_t sum[1000] = { 0 };
	sum[0] = ENTL_A[0];
	sum[1] = ENTL_A[1];
	strcat((char*)(sum + 2), (char*)ID_A);
	ln_to_str(sum + 2 + strlen((char*)ID_A), a, ECC_len, 4 * ln_get_digits(a, ECC_len));
	ln_to_str(sum + 2 + strlen((char*)ID_A) + 32, b, ECC_len, 4 * ln_get_digits(b, ECC_len));
	ln_to_str(sum + 2 + strlen((char*)ID_A) + 64, G.x, ECC_len, 4 * ln_get_digits(G.x, ECC_len));
	ln_to_str(sum + 2 + strlen((char*)ID_A) + 96, G.y, ECC_len, 4 * ln_get_digits(G.y, ECC_len));
	ln_to_str(sum + 2 + strlen((char*)ID_A) + 128, Pb.x, ECC_len, 4 * ln_get_digits(Pb.x, ECC_len));
	ln_to_str(sum + 2 + strlen((char*)ID_A) + 160, Pb.y, ECC_len, 4 * ln_get_digits(Pb.y, ECC_len));

	LN Z_A[ECC_len] = { 0 };
	SM3_encryption(sum, 2 + strlen((char*)ID_A) + 192, Z_A);

	LN r[ECC_len], s[ECC_len];
	SM2_signature(M, Z_A, d, r, s);
	SM2_verification(M, Z_A, Pb, r, s);
	system("pause");
	return 0;
}
