#include "LargeNumber.h"

#define shift_2_bits(x) (uint32_t)(((x) >> (LN_BITS - 2)) & 0x03)
#define shift_4_bits(x) (uint32_t)(((x) >> (LN_BITS - 4)) & 0x0f)

#define assign_one_digit(x, y, digits) \
	{                                  \
		ln_assign_zero(x, digits);     \
		x[0] = y;                      \
	}



LN ln_l_shift(LN* a, LN* b, uint32_t c, uint32_t digits)
{
	LN temp, carry;
	uint32_t cnt;
	if (c >= LN_BITS)
		return 0;
	carry = 0;
	cnt = LN_BITS - c;
	for (uint32_t i = 0; i < digits; i++)
	{
		temp = b[i];
		a[i] = (temp << c) | carry;
		carry = c ? (temp >> cnt) : 0;
	}
	return carry;
}

LN ln_r_shift(LN* a, LN* b, uint32_t c, uint32_t digits)
{
	LN temp, carry;
	uint32_t cnt;
	if (c >= LN_BITS)
		return 0;
	carry = 0;
	cnt = LN_BITS - c;
	int i = digits - 1;
	for (; i >= 0; i--)
	{
		temp = b[i];
		a[i] = (temp >> c) | carry;
		carry = c ? (temp << cnt) : 0;
	}
	return carry;
}

void str_to_ln(uint8_t* str, LN* ln, uint32_t digits, uint32_t size)
{
	LN temp;
	uint32_t i, u;
	int j;
	for (i = 0, j = size - 1; i < digits && j >= 0; i++)
	{
		temp = 0;
		for (u = 0; j >= 0 && u < LN_BITS; j--, u += 8)
		{
			temp |= ((LN)str[j]) << u;
		}
		ln[i] = temp;
	}
	for (; i < digits; i++)
	{
		ln[i] = 0;
	}
}

void ln_to_str(uint8_t* str, LN* ln, uint32_t digits, uint32_t size)
{
	LN temp;
	uint32_t i, u;
	int j;
	for (i = 0, j = size - 1; i < digits && j >= 0; i++)
	{
		temp = ln[i];
		for (u = 0; j >= 0 && u < LN_BITS; j--, u += 8)
		{
			str[j] = (uint8_t)(temp >> u);
		}
	}
	for (; j >= 0; j--)
	{
		str[j] = 0;
	}
	int cnt = 0;
	if ((str[0] == 0) && (str[1] != 0)) {
		cnt = 1;
	}
	else if ((str[0] == 0) && (str[1] == 0) && (str[2] != 0)) {
		cnt = 2;
	}
	else if ((str[0] == 0) && (str[1] == 0) && (str[2] == 0) && (str[3] != 0)) {
		cnt = 3;
	}
	for (j = 0; j < size - cnt; j++)
		str[j] = str[j + cnt];
	for (j = size - cnt; j < size; j++)
		str[j] = 0;
}

void ln_assign_zero(LN* a, uint32_t digits)
{
	for (uint32_t i = 0; i < digits; i++)
	{
		a[i] = 0;
	}
}

void ln_assign(LN* a, LN* b, uint32_t digits)
{
	for (uint32_t i = 0; i < digits; i++)
	{
		a[i] = b[i];
	}
}

void ln_assign_2exp(LN* a, uint32_t b, uint32_t digits)
{
	ln_assign_zero(a, digits);
	if (b >= (digits * LN_BITS))
	{
		return; // Out of range
	}

	a[b / LN_BITS] = (LN)1 << (b % LN_BITS);
}

static uint32_t ln_digit_bits(LN a)
{
	uint32_t i;
	for (i = 0; i < LN_BITS; i++)
	{
		if (a == 0)
			break;
		a >>= 1;
	}
	return i;
}

uint32_t ln_get_bits(LN* a, uint32_t digits)
{
	digits = ln_get_digits(a, digits);
	if (digits == 0)
		return 0;
	return (ln_digit_bits(a[digits - 1]) + LN_BITS * (digits - 1));
}

uint32_t ln_get_digits(LN* a, uint32_t digits)
{
	int i;
	uint32_t cnt = 0;
	for (i = digits - 1; i >= 0; i--)
	{
		if (a[i])
			break;
	}
	i++;
	return i;
}

int ln_cmp(LN* a, LN* b, uint32_t digits)
{
	for (int i = digits - 1; i >= 0; i--)
	{
		if (a[i] > b[i])
			return 1;
		if (a[i] < b[i])
			return -1;
	}
	return 0;
}

int ln_is_zero(LN* a, uint32_t digits)
{
	for (uint32_t i = 0; i < digits; i++)
	{
		if (a[i])
		{
			return 0;
		}
	}
	return 1;
}

int ln_is_one(LN* a, uint32_t digits)
{
	if (a[0] != 1)
		return 0;
	for (uint32_t i = 1; i < digits; i++)
	{
		if (a[i])
		{
			return 0;
		}
	}
	return 1;
}

LN ln_add(LN* a, LN* b, LN* c, uint32_t digits)
{
	LN carry = 0;
	LN temp;
	for (uint32_t i = 0; i < digits; i++)
	{
		temp = b[i] + carry;
		if (temp < carry)
		{
			// b+carry进位，只可能结果�?0
			temp = c[i];
		}
		else if ((temp += c[i]) < c[i])
		{
			// b+c+carry进位
			carry = 1;
		}
		else
		{
			carry = 0;
		}
		a[i] = temp;
	}
	return carry;
}

LN ln_sub(LN* a, LN* b, LN* c, uint32_t digits)
{
	LN borrow = 0;
	LN temp;
	for (uint32_t i = 0; i < digits; i++)
	{
		temp = b[i] - borrow;
		if (temp > (LN_MAX - borrow))
		{
			temp = LN_MAX - c[i];
		}
		else if ((temp -= c[i]) > (LN_MAX - c[i]))
		{
			borrow = 1;
		}
		else
		{
			borrow = 0;
		}
		a[i] = temp;
	}
	return borrow;
}

LN ln_sub_digit_mul(LN* a, LN* b, LN c, LN* d, uint32_t digits)
{
	dLN res;
	LN borrow, res_h, res_l;
	if (c == 0)
		return 0;

	borrow = 0;
	for (uint32_t i = 0; i < digits; i++)
	{
		res = d[i] * (dLN)c;
		res_h = (res >> LN_BITS) & LN_MAX;
		res_l = res & LN_MAX;
		a[i] = b[i] - borrow;
		if (a[i] > LN_MAX - borrow)
		{
			borrow = 1;
		}
		else
		{
			borrow = 0;
		}
		if ((a[i] -= res_l) > (LN_MAX - res_l))
		{
			borrow++;
		}
		borrow += res_h;
	}
	return borrow;
}

//void ln_mul(LN* a, LN* b, LN* c, uint32_t digits) // Blakley
//{												  // 有待优化
//	LN temp[2 * LN_MAX_DIGITS];
//	ln_assign_zero(temp, 2 * digits);
//	uint32_t b_digits, c_digits;
//	b_digits = ln_get_digits(b, digits);
//	c_digits = ln_get_digits(c, digits);
//
//	for (uint32_t i = 0; i < b_digits; i++)
//		temp[i + c_digits] += ln_add_digit_mul(&temp[i], &temp[i], b[i], c, c_digits);
//
//	ln_assign(a, temp, 2 * digits);
//}

void ln_mul(LN* c, LN* a, LN* b, uint32_t digits)
{
	LN t[2 * LN_MAX_DIGITS];
	ln_assign_zero(t, 2 * digits);

	LN C = 0, S = 0;
	dLN CS = 0;
	int s1 = ln_get_digits(a, digits);
	int s2 = ln_get_digits(b, digits);
	for (int i = 0; i < s2; i++) {
		C = 0;
		for (int j = 0; j < s1; j++) {
			CS = t[i + j] + a[j] * (dLN)b[i] + C;
			S = CS & LN_MAX;
			C = (CS >> LN_BITS) & LN_MAX;
			t[i + j] = S;
		}
		t[i + s1] = C;
	}
	ln_assign(c, t, 2 * digits);
}

void ln_div(LN* a, LN* b, LN* c, uint32_t cdigits, LN* d, uint32_t ddigits)
{
	dLN temp;
	LN ai, t, c_temp[2 * LN_MAX_DIGITS + 1], d_temp[LN_MAX_DIGITS];
	int i;
	uint32_t dddigits, shift;

	dddigits = ln_get_digits(d, ddigits);
	if (dddigits == 0)
		return;
	if (ln_is_zero(c, cdigits) == 1) {
		ln_assign_zero(a, cdigits);
		ln_assign_zero(b, dddigits);
		return;
	}

	shift = LN_BITS - ln_digit_bits(d[dddigits - 1]);
	ln_assign_zero(c_temp, dddigits);
	c_temp[cdigits] = ln_l_shift(c_temp, c, shift, cdigits);
	ln_l_shift(d_temp, d, shift, dddigits);
	t = d_temp[dddigits - 1];

	ln_assign_zero(a, cdigits);
	i = cdigits - dddigits;
	for (; i >= 0; i--)
	{
		if (t == LN_MAX)
		{
			ai = c_temp[i + dddigits];
		}
		else
		{
			temp = c_temp[i + dddigits - 1];
			temp += (dLN)c_temp[i + dddigits] << LN_BITS;
			ai = temp / (t + 1);
		}

		c_temp[i + dddigits] -= ln_sub_digit_mul(&c_temp[i], &c_temp[i], ai, d_temp, dddigits);
		while (c_temp[i + dddigits] || (ln_cmp(&c_temp[i], d_temp, dddigits) >= 0))
		{
			ai++;
			c_temp[i + dddigits] -= ln_sub(&c_temp[i], &c_temp[i], d_temp, dddigits);
		}
		a[i] = ai;
	}

	ln_assign_zero(b, ddigits);
	ln_r_shift(b, c_temp, shift, dddigits);
}

void ln_mod(LN* a, LN* b, uint32_t bdigits, LN* c, uint32_t cdigits)
{
	LN temp[2 * LN_MAX_DIGITS] = { 0 };
	ln_div(temp, a, b, bdigits, c, cdigits);
}

void ln_shift_mod(LN* a, LN* b, uint32_t bdigits, LN* c, uint32_t cdigits)
{
	int bbits, cbits, k;
	bbits = ln_get_bits(b, bdigits);
	cbits = ln_get_bits(c, cdigits);
	k = bbits - cbits;
	if (k < 0)
	{
		ln_assign(a, b, bdigits);
		return;
	}
	LN R1[LN_MAX_DIGITS], n[LN_MAX_DIGITS], k2[LN_MAX_DIGITS];

	if (k == 0)
	{
		if (ln_cmp(b, c, bdigits) == 1)
			ln_sub(a, b, c, bdigits);
		else if (ln_cmp(b, c, bdigits) == 0)
			ln_assign_zero(a, LN_MAX_DIGITS);
		else
			ln_assign(a, b, bdigits);
		return;
	}
	ln_assign(R1, b, bdigits);
	ln_assign_2exp(k2, k, LN_MAX_DIGITS);
	ln_mul(n, k2, c, LN_MAX_DIGITS);
	for (int i = 0; i <= k; i++)
	{
		if (ln_cmp(R1, n, LN_MAX_DIGITS) >= 0)
		{
			ln_sub(R1, R1, n, LN_MAX_DIGITS);
		}
		ln_r_shift(n, n, 1, LN_MAX_DIGITS);
	}
	ln_assign(a, R1, LN_MAX_DIGITS);
	return;
}

void ln_mod_mul(LN* a, LN* b, LN* c, LN* d, uint32_t digits)
{
	LN temp[2 * LN_MAX_DIGITS];
	ln_assign_zero(temp, 2 * LN_MAX_DIGITS);
	ln_mul(temp, b, c, digits);
	ln_mod(a, temp, 2 * digits, d, digits);
}

void ln_mod_exp(LN* a, LN* b, LN* c, uint32_t cdigits, LN* d, uint32_t ddigits)
{
	LN b_exp[3][LN_MAX_DIGITS], temp[LN_MAX_DIGITS], ci;
	uint32_t ci_bits, s;
	ln_assign(b_exp[0], b, ddigits);
	ln_mod_mul(b_exp[1], b, b_exp[0], d, ddigits);
	ln_mod_mul(b_exp[2], b, b_exp[1], d, ddigits);

	ln_assign_zero(temp, ddigits);
	temp[0] = 1;

	cdigits = ln_get_digits(c, cdigits);
	int i = cdigits - 1;
	for (; i >= 0; i--)
	{
		ci = c[i];
		ci_bits = LN_BITS;
		if (i == (int)(cdigits - 1))
		{
			while (!shift_2_bits(ci))
			{
				ci <<= 2;
				ci_bits -= 2;
			}
		}
		for (int j = 0; j < ci_bits; j += 2)
		{
			// print_ln(temp,LN_MAX_DIGITS);
			ln_mod_mul(temp, temp, temp, d, ddigits);
			// print_ln(temp,LN_MAX_DIGITS);
			ln_mod_mul(temp, temp, temp, d, ddigits);
			// print_ln(temp,LN_MAX_DIGITS);
			s = shift_2_bits(ci);
			if (s != 0)
			{
				ln_mod_mul(temp, temp, b_exp[s - 1], d, ddigits);
			}
			ci <<= 2;
		}
	}
	ln_assign(a, temp, ddigits);
}

void ln_mod_exp_binary(LN* b, LN* a, LN* e, uint32_t edigits, LN* m, uint32_t mdigits)
{
	LN temp[LN_MAX_DIGITS], p[LN_MAX_DIGITS], E[LN_MAX_DIGITS];
	ln_assign(temp, a, LN_MAX_DIGITS);

	int k = ln_get_bits(e, edigits);
	int ei, num, cnt, esh;
	for (int i = k - 2; i >= 0; i--)
	{
		num = i / LN_BITS;
		cnt = i - num * LN_BITS;
		ln_assign(E, e, LN_MAX_DIGITS);
		ln_mod_mul(p, temp, temp, m, mdigits);
		esh = E[num] >> (cnt);
		ei = esh & 1U;

		if (ei == 1)
		{
			ln_mod_mul(temp, p, a, m, mdigits);
		}
		else
		{
			ln_assign(temp, p, LN_MAX_DIGITS);
		}
		// ln_assign_zero(esh, LN_MAX_DIGITS);
	}
	ln_assign(b, temp, LN_MAX_DIGITS);
}

void ln_mod_exp_mary(LN* a, LN* b, LN* c, uint32_t cdigits, LN* d, uint32_t ddigits)
{
	LN b_exp[33][LN_MAX_DIGITS], temp[LN_MAX_DIGITS], ci;
	uint32_t ci_bits, s;
	ln_assign(b_exp[0], b, ddigits);
	for (int i = 1; i < 32; i++)
	{
		ln_mod_mul(b_exp[i], b, b_exp[i - 1], d, ddigits);
	}

	ln_assign_zero(temp, ddigits);
	temp[0] = 1;

	cdigits = ln_get_digits(c, cdigits);
	int i = cdigits - 1;
	for (; i >= 0; i--)
	{
		ci = c[i];
		ci_bits = LN_BITS;
		if (i == (int)(cdigits - 1))
		{
			while (!shift_4_bits(ci))
			{
				ci <<= 4;
				ci_bits -= 4;
			}
		}
		for (int j = 0; j < ci_bits; j += 4)
		{
			// print_ln(temp,LN_MAX_DIGITS);
			ln_mod_mul(temp, temp, temp, d, ddigits);
			// print_ln(temp,LN_MAX_DIGITS);
			ln_mod_mul(temp, temp, temp, d, ddigits);
			ln_mod_mul(temp, temp, temp, d, ddigits);
			ln_mod_mul(temp, temp, temp, d, ddigits);
			// print_ln(temp,LN_MAX_DIGITS);
			s = shift_4_bits(ci);
			if (s != 0)
			{
				ln_mod_mul(temp, temp, b_exp[s - 1], d, ddigits);
			}
			ci <<= 4;
		}
	}
	ln_assign(a, temp, ddigits);
}

void crt_ln_mod_exp(LN* c, LN* d, LN* p, LN* q, LN* res, uint32_t digits)
{
	LN d1[LN_MAX_DIGITS], d2[LN_MAX_DIGITS], C[LN_MAX_DIGITS], D[LN_MAX_DIGITS], p1[LN_MAX_DIGITS], q1[LN_MAX_DIGITS];
	LN one[LN_MAX_DIGITS];
	LN m1[LN_MAX_DIGITS], m2[LN_MAX_DIGITS], m22[LN_MAX_DIGITS], m[2 * LN_MAX_DIGITS];
	LN P[LN_MAX_DIGITS], Q[LN_MAX_DIGITS];
	assign_one_digit(one, 1, digits);
	ln_assign(C, c, digits);
	ln_assign(D, d, digits);
	ln_assign(P, p, digits);
	ln_assign(Q, q, digits);
	ln_assign(p1, p, digits);
	ln_assign(q1, q, digits);
	ln_sub(q1, q1, one, digits);
	ln_sub(p1, p1, one, digits);
	ln_mod(d1, D, digits, p1, digits);
	ln_mod(d2, D, digits, q1, digits);
	ln_mod_exp_mary(m1, C, d1, digits, P, digits);
	ln_mod_exp_mary(m2, C, d2, digits, Q, digits);
	LN temp[LN_MAX_DIGITS], p_1[LN_MAX_DIGITS];
	if (ln_cmp(m1, m2, digits) >= 0)
	{
		ln_sub(temp, m1, m2, digits);
	}
	else
	{
		ln_sub(temp, m2, m1, digits);
	}
	mod_inv(p_1, P, Q, digits);
	ln_mod_mul(m22, temp, p_1, Q, digits);
	ln_mul(m, m22, P, digits);
	ln_add(m, m, m1, digits);
	ln_assign(res, m, digits);
}

// void mod_inv(LN *a, LN *b, LN *c, uint32_t digits)	//欧拉定理直接求模逆，效率过低
// {
// 	LN temp[LN_MAX_DIGITS], one[LN_MAX_DIGITS], res[LN_MAX_DIGITS];
// 	assign_one_digit(temp, 1, LN_MAX_DIGITS);
// 	assign_one_digit(one, 1, LN_MAX_DIGITS);
// 	while (ln_cmp(temp, c, LN_MAX_DIGITS) == -1)
// 	{
// 		ln_mod_mul(res, b, temp, c, LN_MAX_DIGITS);
// 		if (ln_cmp(res, one, LN_MAX_DIGITS) == 0)
// 		{
// 			ln_assign(a, temp, LN_MAX_DIGITS);

// 			memset((uint8_t *)one, 0, sizeof(one));
// 			memset((uint8_t *)temp, 0, sizeof(temp));
// 			memset((uint8_t *)res, 0, sizeof(res));
// 			return;
// 		}
// 		ln_add(temp, temp, one, LN_MAX_DIGITS);
// 	}
// }

void mod_inv(LN* a, LN* b, LN* c, uint32_t digits)
{ // a=b-1modc
	// int m = -1, n;
	LN m[LN_MAX_DIGITS], n[LN_MAX_DIGITS], mod[LN_MAX_DIGITS], N[LN_MAX_DIGITS], B[LN_MAX_DIGITS], C[LN_MAX_DIGITS];

	ln_assign(N, c, digits);
	ln_assign(C, c, digits);
	ln_assign(B, b, digits);

	assign_one_digit(m, 1, LN_MAX_DIGITS);
	// int temp;
	LN temp1[LN_MAX_DIGITS], test[LN_MAX_DIGITS], temp2[2 * LN_MAX_DIGITS];
	// if (a < b) {
	// 	temp = a, a = b, b = temp;
	// }
	LN bi[500][LN_MAX_DIGITS], ai[500][LN_MAX_DIGITS];
	int cnt = 0;
	while (!ln_is_zero(m, digits))
	{

		// n = a / b;
		ln_div(n, mod, c, digits, b, digits);

		// temp = a, a = b;
		ln_assign(temp1, c, digits);
		ln_assign(c, b, digits);
		// m = temp - a * n;
		ln_mul(temp2, c, n, digits);
		ln_sub(m, temp1, temp2, digits);
		// b = m;
		ln_assign(b, m, digits);

		ln_assign(ai[cnt], n, digits);
		cnt++;
	}
	cnt--;
	assign_one_digit(bi[0], 1, digits);
	ln_assign(bi[1], ai[cnt - 1], digits);
	for (int i = 2; i <= cnt; i++)
	{
		ln_mul(bi[i], ai[cnt - i], bi[i - 1], digits);
		ln_add(bi[i], bi[i], bi[i - 2], digits);
	}
	// cout << "cnt:" << cnt << endl;
	if (cnt % 2 != 0)
	{
		ln_sub(a, N, bi[cnt], digits);
	}
	else
	{
		ln_assign(a, bi[cnt], digits);
	}
	ln_assign(b, B, digits);
	ln_assign(c, C, digits);
}

bool ln_gcd(LN* b, LN* c, uint32_t digits)
{
	// int m = -1, n;
	LN m[LN_MAX_DIGITS], n[LN_MAX_DIGITS], mod[LN_MAX_DIGITS], N[LN_MAX_DIGITS], B[LN_MAX_DIGITS], C[LN_MAX_DIGITS];

	ln_assign(N, c, digits);
	ln_assign(C, c, digits);
	ln_assign(B, b, digits);
	assign_one_digit(m, 1, LN_MAX_DIGITS);
	// int temp;
	LN temp1[LN_MAX_DIGITS], test[LN_MAX_DIGITS], temp2[2 * LN_MAX_DIGITS];
	// if (a < b) {
	// 	temp = a, a = b, b = temp;
	// }
	int cnt = 0;
	while (!ln_is_zero(m, digits))
	{
		// n = a / b;
		ln_div(n, mod, c, digits, b, digits);

		// temp = a, a = b;
		ln_assign(temp1, c, digits);
		ln_assign(c, b, digits);
		// m = temp - a * n;
		ln_mul(temp2, c, n, digits);
		ln_sub(m, temp1, temp2, digits);
		// b = m;
		ln_assign(b, m, digits);
		cnt++;
	}
	if (ln_is_zero(m, digits) && ln_is_one(c, digits))
	{
		ln_assign(b, B, digits);
		ln_assign(c, C, digits);
		return true;
	}
	else
	{
		ln_assign(b, B, digits);
		ln_assign(c, C, digits);
		return false;
	}
}

void print_ln(LN* a, uint32_t digits)
{
	uint8_t str[512];
	memset(str, 0, 512);
	ln_to_str(str, a, digits, 512);
	int i = 0;
	while (str[i] == 0)
		i++;
	if (i == 512)
	{
		printf("0\n");
		return;
	}
	for (; i < 512; i++)
		printf("%02X", str[i]);
	printf("\n");
}



void MonPro(LN* res, LN* a, LN* b, LN* r, LN* n, LN* n_1, uint32_t digits)
{
	LN t[2 * LN_MAX_DIGITS], m[LN_MAX_DIGITS], temp[LN_MAX_DIGITS], u[LN_MAX_DIGITS];
	LN mod[LN_MAX_DIGITS];
	ln_mul(t, a, b, digits);
	// mod_inv(n_1,n,r,digits);
	ln_mod_mul(m, t, n_1, r, digits); // shift
	ln_mul(temp, m, n, digits);
	ln_add(temp, temp, t, digits);
	ln_div(u, mod, temp, digits, r, digits);
	if (ln_cmp(u, n, digits) >= 0)
	{
		ln_sub(res, u, n, digits);
	}
	else
	{
		ln_assign(res, u, digits);
	}
	// ln_mod(res,u,digits,n,digits);
	return;
}

void MonMul(LN* res, LN* a, LN* b, LN* n, LN* n_1, LN* r, uint32_t digits)
{
	// // clock_t start,finish;
	// uint32_t bits = ln_get_bits(n, digits);
	// LN r[LN_MAX_DIGITS], a_1[LN_MAX_DIGITS], r_1[LN_MAX_DIGITS], n_1[LN_MAX_DIGITS], x[LN_MAX_DIGITS];
	// ln_assign_2exp(r, bits, digits);
	// // start = clock();
	// mod_inv(r_1, r, n, digits); // r_1
	// // finish = clock();
	// // double duration = finish - start;
	// // printf("模逆时间%f ms",duration);
	// LN temp[LN_MAX_DIGITS], mod[LN_MAX_DIGITS];
	// ln_mul(temp, r, r_1, digits);
	// ln_div(n_1, mod, temp, digits, n, digits);
	LN a_1[LN_MAX_DIGITS];
	ln_mod_mul(a_1, a, r, n, digits);
	MonPro(res, a_1, b, r, n, n_1, digits);
	return;
}

void ModExp(LN* res, LN* M, LN* e, LN* n, uint32_t digits)
{
	uint32_t bits = ln_get_bits(n, digits);
	LN r[LN_MAX_DIGITS], a_1[LN_MAX_DIGITS], r_1[LN_MAX_DIGITS], n_1[LN_MAX_DIGITS], x[LN_MAX_DIGITS];
	ln_assign_2exp(r, bits, digits);

	mod_inv(r_1, r, n, digits); // r_1

	LN temp[LN_MAX_DIGITS], mod[LN_MAX_DIGITS];
	ln_mul(temp, r, r_1, digits);
	ln_div(n_1, mod, temp, digits, n, digits);
	LN M_1[LN_MAX_DIGITS], x_1[LN_MAX_DIGITS];
	ln_mod_mul(M_1, M, r, n, digits);
	ln_mod(x_1, r, digits, n, digits);

	LN E[LN_MAX_DIGITS];
	ln_assign(temp, M, LN_MAX_DIGITS);

	int k = ln_get_bits(e, digits);
	int ei, num, cnt, esh;
	for (int i = k - 1; i >= 0; i--)
	{
		num = i / LN_BITS;
		cnt = i - num * LN_BITS;
		ln_assign(E, e, LN_MAX_DIGITS);
		// ln_mod_mul(p, temp, temp, m, mdigits);
		MonPro(x_1, x_1, x_1, r, n, n_1, digits);
		esh = E[num] >> (cnt);
		ei = esh & 1U;

		if (ei == 1)
		{
			// ln_mod_mul(temp, p, a, m, mdigits);
			MonPro(x_1, M_1, x_1, r, n, n_1, digits);
		}
	}
	LN one[LN_MAX_DIGITS];
	assign_one_digit(one, 1, LN_MAX_DIGITS);
	MonPro(res, x_1, one, r, n, n_1, digits);
}
