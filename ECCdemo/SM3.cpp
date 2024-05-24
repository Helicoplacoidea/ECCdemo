#include "SM3.h"

uint32_t w[68] = { 0 };
uint32_t w1[64] = { 0 };
uint32_t A = 0x7380166f, B = 0x4914b2b9, C = 0x172442d7, D = 0xda8a0600, E = 0xa96f30bc, F = 0x163138aa, G = 0xe38dee4d, H = 0xb0fb0e4e;
uint32_t V0 = 0x7380166f, V1 = 0x4914b2b9, V2 = 0x172442d7, V3 = 0xda8a0600, V4 = 0xa96f30bc, V5 = 0x163138aa, V6 = 0xe38dee4d, V7 = 0xb0fb0e4e;
uint32_t V00 = 0x7380166f, V11 = 0x4914b2b9, V22 = 0x172442d7, V33 = 0xda8a0600, V44 = 0xa96f30bc, V55 = 0x163138aa, V66 = 0xe38dee4d, V77 = 0xb0fb0e4e;
uint32_t SS1 = 0, SS2 = 0, TT1 = 0, TT2 = 0;
uint32_t T[4] = { 0x79cc4519,0x7a879d8a,0x7a879d8a,0x7a879d8a };


void str_to_int(uint8_t* input, uint32_t* plaintext)
{
	for (int i = 0; i < 64; i += 4) {
		plaintext[i / 4] = (((uint32_t)input[i] & 0x000000ff) << 24) |
			(((uint32_t)input[i + 1] & 0x000000ff) << 16) |
			(((uint32_t)input[i + 2] & 0x000000ff) << 8) |
			(((uint32_t)input[i + 3] & 0x000000ff));
	}
}

int add_padding(uint8_t* input, uint32_t* plaintext, uint64_t msglen, long long len)
{
	//int len = strlen((const char*)input);
	//int len = sizeof(input);
	if ((msglen * 8) % 512 <= 440 && (msglen * 8) % 512 > 0) {
		input[len] = 128;
		for (int i = len + 1; i < 60; i++) {
			input[i] = 0;
		}
		for (int i = 0; i < 56; i += 4) {
			plaintext[i / 4] = (((uint32_t)input[i] & 0x000000ff) << 24) |
				(((uint32_t)input[i + 1] & 0x000000ff) << 16) |
				(((uint32_t)input[i + 2] & 0x000000ff) << 8) |
				(((uint32_t)input[i + 3] & 0x000000ff));
		}
		plaintext[14] = ((msglen * 8) >> 32) & 0xffffffff;
		plaintext[15] = (msglen * 8) & 0xffffffff;
		return 0;
	}
	else if ((msglen * 8) % 512 >= 448 || (msglen * 8) % 512 == 0) {
		input[len] = 128;
		for (int i = len + 1; i < 124; i++) {
			input[i] = 0;
		}
		for (int i = 0; i < 120; i += 4) {
			plaintext[i / 4] = (((uint32_t)input[i] & 0x000000ff) << 24) |
				(((uint32_t)input[i + 1] & 0x000000ff) << 16) |
				(((uint32_t)input[i + 2] & 0x000000ff) << 8) |
				(((uint32_t)input[i + 3] & 0x000000ff));
		}
		plaintext[30] = ((msglen * 8) >> 32) & 0xffffffff;
		plaintext[31] = (msglen * 8) & 0xffffffff;
		return 1;
	}
}

void get_W(uint32_t* plaintext)
{
	int j = 0;
	for (; j < 16; j++) {
		w[j] = plaintext[j];
	}
	for (; j < 68; j++) {
		w[j] = P1(w[j - 16] ^ w[j - 9] ^ ROTL(w[j - 3], 15)) ^ ROTL(w[j - 13], 7) ^ w[j - 6];
	}
	for (j = 0; j < 64; j++) {
		w1[j] = w[j] ^ w[j + 4];
	}
}

void step_function()
{
	A = V0;
	B = V1;
	C = V2;
	D = V3;
	E = V4;
	F = V5;
	G = V6;
	H = V7;
	for (int j = 0; j < 64; j++) {
		SS1 = ROTL(ROTL(A, 12) + E + ROTL(T[j / 16], j), 7);
		SS2 = SS1 ^ ROTL(A, 12);
		if (j >= 0 && j <= 15) {
			TT1 = FF0 + D + SS2 + w1[j];
			TT2 = GG0 + H + SS1 + w[j];
		}
		else {
			TT1 = FF1 + D + SS2 + w1[j];
			TT2 = GG1 + H + SS1 + w[j];
		}
		D = C;
		C = ROTL(B, 9);
		B = A;
		A = TT1;
		H = G;
		G = ROTL(F, 19);
		F = E;
		E = P0(TT2);
		//printf("%08x %08x %08x %08x %08x %08x %08x %08x\n", A, B, C, D, E, F, G, H);
	}
	V0 = A ^ V0;
	V1 = B ^ V1;
	V2 = C ^ V2;
	V3 = D ^ V3;
	V4 = E ^ V4;
	V5 = F ^ V5;
	V6 = G ^ V6;
	V7 = H ^ V7;
}

void round(uint32_t* plaintext)
{
	get_W(plaintext);
	step_function();
	printf("%08x %08x %08x %08x %08x %08x %08x %08x\n", V0, V1, V2, V3, V4, V5, V6, V7);
}


//uint64_t msgsize(char* plainaddr)
//{
//	uint64_t size = 0;
//	FILE* fp = NULL;
//	if ((fp = fopen(plainaddr, "rb")) == NULL)
//	{
//		printf("open %s faied!", plainaddr);
//		exit(0);
//	}
//
//	fseek(fp, 0, SEEK_END);
//	size = ftell(fp);
//	fseek(fp, 0, SEEK_SET);
//	fclose(fp);
//	return size;
//}

//void SM3_encryption(char* inputfileaddr)
//{
//	uint8_t input[128] = { 0 };
//	uint32_t plaintext[64] = { 0 };
//	uint64_t msglen = msgsize(inputfileaddr);
//	printf("长度%lld\n", msglen);
//
//	FILE* fp;
//	if ((fp = fopen(inputfileaddr, "r")) == NULL)
//	{
//		perror("fail to read");
//		exit(1);
//	}
//	int len = 0;
//	uint64_t remain = msglen;
//	while (fgets((char*)input, MAX_LINE + 1, fp) != NULL) {
//		if (remain <= 64) {
//			int flag = add_padding(input, plaintext, msglen);
//			round(plaintext);
//			if (flag == 1) {
//				round(&plaintext[16]);
//			}
//		}
//		else {
//			str_to_int(input, plaintext);
//			round(plaintext);
//			remain = remain - 64;
//		}
//	}
//}

void SM3_encryption(uint8_t* input, uint64_t msglen, uint32_t* output)
{
	uint32_t plaintext[64] = { 0 };
	printf("长度%lld\n", msglen);

	int len = 0;
	long long remain = msglen;
	int i = 0;
	while (remain >= 0) {
		if (remain <= 64) {
			//printf("%d\n", input[64]);
			//printf("%d\n", input[65]);
			//printf("%d\n", input[66]);
			//printf("%d\n", input[67]);
			int flag = add_padding(input + 64 * i, plaintext, msglen, remain);

			round(plaintext);
			if (flag == 1) {
				round(&plaintext[16]);
			}
			remain = remain - 64;
			output[0] = V7, output[1] = V6, output[2] = V5, output[3] = V4,
			output[4] = V3, output[5] = V2, output[6] = V1, output[7] = V0;
			V0 = V00, V1 = V11, V2 = V22, V3 = V33, V4 = V44, V5 = V55, V6 = V66, V7 = V77;
		}
		else {
			str_to_int(input + 64 * i, plaintext);
			round(plaintext);
			remain = remain - 64;
		}
		i++;
	}
}