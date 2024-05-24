#pragma once
#include <stdint.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#define MAX_LINE 64
#define ROTL(x, n) (((x) << (n)) | ((x) >> (32 - n)))
#define P0(x) x ^ ROTL(x, 9) ^ ROTL(x, 17)
#define P1(x) x ^ ROTL(x, 15) ^ ROTL(x, 23)
#define FF0 (A ^ B ^ C)
#define FF1 ((A & B) | (A & C) | (B & C))
#define GG0 (E ^ F ^ G)
#define GG1 ((E & F) | (~E & G))

void str_to_int(uint8_t* input, uint32_t* plaintext);
int add_padding(uint8_t* input, uint32_t* plaintext, uint64_t msglen, long long len);
void get_W(uint32_t* plaintext);
void step_function();
void round(uint32_t* plaintext);
//uint64_t msgsize(char* plainaddr);
//void SM3_encryption(char* inputfileaddr);
void SM3_encryption(uint8_t* input, uint64_t msglen, uint32_t* output);