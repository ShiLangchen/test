#include <cassert> 

#include <iostream>
#include <fstream>  
#include <stdio.h>
#include <ctime>
#include <set>
#include <chrono> 
using namespace std;

//LOW_MC_ALGORITHM

int ppp = 0;//
int ans = 1;//

#define round 2
// #define block_size 60
unsigned char plain_glo[block_size / 3] = { 0 };	//global m
unsigned char cipher_glo[block_size / 3] = { 0 };	//global c

unsigned char plain[block_size / 3] = { 0};		//m
unsigned char cipher[block_size / 3] = { 0 };	//c

// state
unsigned sta_l1[block_size] = { 0 };	//l1
unsigned sta_k2[block_size] = { 0 };	//k2
unsigned sta_s2[block_size] = { 0 };	//s2
unsigned sta_l2[block_size] = { 0 };	//l2
unsigned sta_k3[block_size] = { 0 };	//k3

#pragma region LOW_MC_ALGORITHM

void encryption_withkey(unsigned char* plaineg, unsigned char* ciphereg, unsigned char* guess_key);
// void exhusted_search()
// {
// 	unsigned char all_key[block_size / 3] = { 0 };
// 	unsigned char pp[block_size / 3], cc[block_size / 3];
// 	for (size_t j = 0; j < ((unsigned long long)1 << block_size - 1); j++)
// 	{
// 		for (size_t i = 0; i < block_size / 3; i++)
// 		{
// 			all_key[i] = (j >> i * 3) & 7;
// 			pp[i] = cc[i];
// 			//cc[i] = 0;
// 		}
// 		encryption_withkey(pp, cc, all_key);
// 	}
// 	// printf("\n%d\n", cc[0]);
// }
void inv_martix(unsigned i_round);	//�����任 inverse
void chushihua_juzhen(unsigned i);	//init_matrix

unsigned char key[round + 1][block_size / 3];


unsigned roundmatrix[round][block_size][block_size];		//for encryption
unsigned inv_roundmatrix[round][block_size][block_size];	//for decryption


unsigned Sbox[8] =
{ 0x00, 0x01, 0x03, 0x06, 0x07, 0x04, 0x05, 0x02 };
unsigned invSbox[8] =
{ 0x00, 0x01, 0x07, 0x02, 0x05, 0x06, 0x03, 0x04 };

// Gaussian elimination to compute inverse matrix
void inv_martix(unsigned i_round)
{
	// Initialize the matrix as an identity matrix (I)
	for (size_t i = 0; i < block_size; i++)
	{
		for (size_t j = 0; j < block_size; j++)
		{
			inv_roundmatrix[i_round][i][j] = 0;
		}
	}
	for (size_t i = 0; i < block_size; i++)
	{
		inv_roundmatrix[i_round][i][i] = 1;
	}

	unsigned tempfinal[block_size][block_size];
	for (size_t i = 0; i < block_size; i++)
	{
		for (size_t j = 0; j < block_size; j++)
		{
			tempfinal[i][j] = roundmatrix[i_round][i][j];
		}
	}
	
	// Gaussian elimination
	// Row exchange
	int exchange;
	for (size_t i = 0; i < block_size; i++)
	{
		if (tempfinal[i][i] == 0)
		{
			for (size_t j = i + 1; j < block_size; j++)
			{
				if (tempfinal[j][i] == 1)
				{
					for (size_t k = 0; k < block_size; k++)
					{
						exchange = tempfinal[i][k];
						tempfinal[i][k] = tempfinal[j][k];
						tempfinal[j][k] = exchange;
					}
					for (size_t k = 0; k < block_size; k++)
					{
						exchange = inv_roundmatrix[i_round][i][k];
						inv_roundmatrix[i_round][i][k] = inv_roundmatrix[i_round][j][k];
						inv_roundmatrix[i_round][j][k] = exchange;
					}
					j = block_size;
				}
				else if (j == block_size - 1)
				{	// Irreversible ���󲻿���
					//printf("\nrank??%d\n", i);
					chushihua_juzhen(i_round);
					return;
				}
			}
			if (i == block_size - 1)
			{	// Irreversible
				//printf("\nrank??%d\n", i);
				chushihua_juzhen(i_round);
				return;
			}

		}
		// elimination ��Ԫ����
		for (size_t j = i + 1; j < block_size; j++)
		{
			if (tempfinal[j][i] == 1)
			{
				for (size_t k = 0; k < block_size; k++)
				{
					tempfinal[j][k] ^= tempfinal[i][k];
				}
				for (size_t k = 0; k < block_size; k++)
				{
					inv_roundmatrix[i_round][j][k] ^= inv_roundmatrix[i_round][i][k];
				}
			}
		}
	}
	// reverse elimination ������Ԫ
	for (int i = block_size - 1; i > -1; i--)
	{
		for (size_t j = 0; j < i; j++)
		{
			if (tempfinal[j][i] == 1)
			{
				for (size_t k = 0; k < block_size; k++)\
				{
					tempfinal[j][k] ^= tempfinal[i][k];
				}
				for (size_t k = 0; k < block_size; k++)
				{
					inv_roundmatrix[i_round][j][k] ^= inv_roundmatrix[i_round][i][k];
				}
			}
		}
	}
}

// initialize key for a 0-7 cycle
unsigned suijizhoingzi = 0;
void chushihua_miyao()
{
	unsigned seed = std::chrono::high_resolution_clock::now().time_since_epoch().count() + suijizhoingzi++;
    srand(seed);

	for (size_t i = 0; i < round + 1; i++)
	{
		for (size_t j = 0; j < block_size / 3; j++)
		{
			key[i][j] = j & 7;		//0~7
			// key[i][j] = rand() % 8;
		}
	}
}


void chushihua_juzhen(unsigned i)
{
	unsigned seed = std::chrono::high_resolution_clock::now().time_since_epoch().count() + suijizhoingzi++;
    srand(seed);
	// srand(time(0) + suijizhoingzi++);
	for (size_t j = 0; j < block_size; j++)
	{
		for (size_t k = 0; k < block_size; k++)
		{
			roundmatrix[i][j][k] = 0;
		}
	}
	for (size_t j = 0; j < block_size; j++)
	{
		for (size_t k = 0; k < block_size; k++)
		{
			roundmatrix[i][j][k] = rand() % 2;	// 0/1
		}
		//roundmatrix[i][j][j] = 1;
	}
	inv_martix(i);
}

void encryption(unsigned char* plain, unsigned char* cipher)
{
	for (size_t i = 0; i < round; i++)
	{
		// XOR: m + k
		for (size_t j = 0; j < block_size / 3; j++)
		{
			plain[j] ^= key[i][j];
		}
		unsigned temp_bit[block_size] = { 0 };

		// Sbox: S(m)
		for (size_t j = 0; j < block_size / 3; j++)
		{
			plain[j] = Sbox[plain[j]];
		}

		// matrix multiplication ����˷����� plain ����ȡ���أ����� roundmatrix[i] ���г˷�����
		for (size_t j = 0; j < block_size; j++)
		{
			for (size_t k = 0; k < block_size; k++)
			{
				temp_bit[j] += ((plain[k / 3] & (1 << (2 - k % 3))) >> (2 - k % 3)) & roundmatrix[i][j][k];
			}
		}

		// ÿ�������غϳ�һ���ֽ�
		for (size_t j = 0; j < block_size / 3; j++)
		{
			plain[j] = (temp_bit[j * 3 + 2] % 2) + ((temp_bit[j * 3 + 1] % 2) << 1) + ((temp_bit[j * 3] % 2) << 2);
		}
	}

	// final round: c = m + k
	for (size_t j = 0; j < block_size / 3; j++)
	{
		plain[j] ^= key[round][j];
		cipher[j] = plain[j];
	}
}

// ʹ�ø����Ĳ²���Կ guess_key
void encryption_withkey(unsigned char* plaineg, unsigned char* ciphereg, unsigned char* guess_key)
{

	for (size_t i = 0; i < round; i++)
	{
		for (size_t j = 0; j < block_size / 3; j++)
		{
			plaineg[j] ^= guess_key[j];
		}

		unsigned temp_bit[block_size] = { 0 };
		for (size_t j = 0; j < block_size / 3; j++)
		{
			plaineg[j] = Sbox[plaineg[j]];
		}

		for (size_t j = 0; j < block_size; j++)
		{
			for (size_t k = 0; k < block_size; k++)
			{
				temp_bit[j] += ((plaineg[k / 3] & (1 << (2 - k % 3))) >> (2 - k % 3)) & roundmatrix[i][j][k];
			}
		}

		for (size_t j = 0; j < block_size / 3; j++)
		{
			plaineg[j] = (temp_bit[j * 3 + 2] % 2) + ((temp_bit[j * 3 + 1] % 2) << 1) + ((temp_bit[j * 3] % 2) << 2);
		}
	}

	for (size_t j = 0; j < block_size / 3; j++)
	{
		plaineg[j] ^= guess_key[j];
		ciphereg[j] = plaineg[j];
	}
}

void decryption(unsigned char* cipher, unsigned char* plain)
{
	for (int i = round; i > 0; i--)
	{
		for (size_t j = 0; j < block_size / 3; j++)
		{
			cipher[j] ^= key[i][j];
		}

		unsigned temp_bit[block_size] = { 0 };
		for (size_t j = 0; j < block_size; j++)
		{
			for (size_t k = 0; k < block_size; k++)
			{
				temp_bit[j] += ((cipher[k / 3] & (1 << (2 - k % 3))) >> (2 - k % 3)) & inv_roundmatrix[i - 1][j][k];
			}
		}

		for (size_t j = 0; j < block_size / 3; j++)
		{
			cipher[j] = (temp_bit[j * 3 + 2] % 2) + ((temp_bit[j * 3 + 1] % 2) << 1) + ((temp_bit[j * 3] % 2) << 2);
		}

		for (size_t j = 0; j < block_size / 3; j++)
		{
			cipher[j] = invSbox[cipher[j]];
		}
	}
	for (size_t j = 0; j < block_size / 3; j++)
	{
		cipher[j] ^= key[0][j];
		plain[j] = cipher[j];
	}
}

void pre_work()
{
	for (size_t i = 0; i < block_size / 3; i++)
	{
		plain_glo[i] = plain[i];
	}

	chushihua_miyao();
	for (size_t i = 0; i < round; i++)
	{
		chushihua_juzhen(i);
	}

	encryption(plain, cipher);

	for (size_t i = 0; i < block_size / 3; i++)
	{
		cipher_glo[i] = cipher[i];
	}
	decryption(cipher, plain);
}

#pragma endregion
void creat_s_matrix(char Ma[block_size][block_size + block_size * block_size / 2 + 1])
{
	char temp_ma[block_size][block_size + block_size * block_size / 2 + 1] = { 0 };
	//z2 z1 z0 \\z2 z1 z0
	for (size_t j = 0; j < block_size / 3; j++)
	{
		//????s??????????
		char x2[block_size + block_size * block_size / 2 + 1] = { 0 }, x1[block_size + block_size * block_size / 2 + 1] = { 0 }, x0[block_size + block_size * block_size / 2 + 1] = { 0 };
		for (size_t k = 0; k < block_size + block_size * block_size / 2 + 1; k++)
		{
			x2[k] ^= Ma[j * 3][k];		//x2
			x1[k] ^= Ma[j * 3 + 1][k];	//x1
			x0[k] ^= Ma[j * 3 + 2][k];	//x0
		}
		
		//z2 = x2 + x0x1
		for (size_t k = 0; k < block_size + block_size * block_size / 2 + 1; k++)
		{
			temp_ma[j * 3][k] ^= Ma[j * 3][k];//x2
		}
		if (x0[block_size + block_size * block_size / 2] == 1)
		{
			for (size_t k = 0; k < block_size + block_size * block_size / 2; k++)
			{
				temp_ma[j * 3][k] ^= Ma[j * 3 + 1][k];//x1
			}
		}

		if (x1[block_size + block_size * block_size / 2] == 1)
		{
			for (size_t k = 0; k < block_size + block_size * block_size / 2; k++)
			{
				temp_ma[j * 3][k] ^= Ma[j * 3 + 2][k];//x0
			}
		}

		if (x1[block_size + block_size * block_size / 2] == 1 && x0[block_size + block_size * block_size / 2] == 1)
		{
			temp_ma[j * 3][block_size + block_size * block_size / 2] ^= 1;//x0
		}

		// x0*x1
		for (size_t k = 0; k < block_size; k++){
			for (size_t l = 0; l < block_size; l++){			
				if (x0[k] == 1 && x1[l] == 1)
				{
					int snum, bignum;
					if (k == l)
					{
						temp_ma[j * 3][k] ^= 1;
					}else{
						(k > l) ? (snum = l, bignum = k) : (snum = k, bignum = l);
						//????????????
						temp_ma[j * 3][block_size + (2 * block_size - 1 - snum) * snum / 2 + (bignum - snum - 1)] ^= 1;
					}
				}
			}
		}

		// z1 = x1 + x2 + x0x2
		if (x0[block_size + block_size * block_size / 2] == 1)
		{
			for (size_t k = 0; k < block_size + block_size * block_size / 2; k++)
			{
				temp_ma[j * 3 + 1][k] ^= Ma[j * 3][k];//x2
			}
		}
		if (x2[block_size + block_size * block_size / 2] == 1)
		{
			for (size_t k = 0; k < block_size + block_size * block_size / 2; k++)
			{
				temp_ma[j * 3 + 1][k] ^= Ma[j * 3 + 2][k];//x0
			}
		}
		if (x2[block_size + block_size * block_size / 2] == 1 && x0[block_size + block_size * block_size / 2] == 1)
		{
			temp_ma[j * 3 + 1][block_size + block_size * block_size / 2] ^= 1;//x0
		}
		for (size_t k = 0; k < block_size + block_size * block_size / 2 + 1; k++)
		{
			temp_ma[j * 3 + 1][k] ^= Ma[j * 3 + 1][k];//x1
			temp_ma[j * 3 + 1][k] ^= Ma[j * 3][k];//x2
		}

		//x0*x2
		for (size_t k = 0; k < block_size; k++)
		{
			for (size_t l = 0; l < block_size; l++)
			{
				if (x0[k] == 1 && x2[l] == 1)
				{
					int snum, bignum;
					if (k == l)
					{
						temp_ma[j * 3 + 1][k] ^= 1;
					}else{									
						(k > l) ? (snum = l, bignum = k) : (snum = k, bignum = l);
						//????????????
						temp_ma[j * 3 + 1][block_size + (2 * block_size - 1 - snum) * snum / 2 + (bignum - snum - 1)] ^= 1;
					}
				}
			}
		}

		//z0 = x0 + x1 + x2 + x1x2
		if (x1[block_size + block_size * block_size / 2] == 1)
		{
			for (size_t k = 0; k < block_size + block_size * block_size / 2; k++)
			{
				temp_ma[j * 3 + 2][k] ^= Ma[j * 3][k];//x2
			}
		}

		if (x2[block_size + block_size * block_size / 2] == 1)
		{
			for (size_t k = 0; k < block_size + block_size * block_size / 2; k++)
			{
				temp_ma[j * 3 + 2][k] ^= Ma[j * 3 + 1][k];//x0
			}
		}

		if (x2[block_size + block_size * block_size / 2] == 1 && x1[block_size + block_size * block_size / 2] == 1)
		{
			temp_ma[j * 3 + 2][block_size + block_size * block_size / 2] ^= 1;//x0
		}
		for (size_t k = 0; k < block_size + block_size * block_size / 2 + 1; k++)
		{
			temp_ma[j * 3 + 2][k] ^= Ma[j * 3 + 1][k];	//x1
			temp_ma[j * 3 + 2][k] ^= Ma[j * 3][k];		//x2
			temp_ma[j * 3 + 2][k] ^= Ma[j * 3 + 2][k];	//x0
		}

		// x1*x2
		for (size_t k = 0; k < block_size; k++)
		{
			for (size_t l = 0; l < block_size; l++)
			{
				if (x1[k] == 1 && x2[l] == 1)
				{
					int snum, bignum;
					if (k == l)
					{

						temp_ma[j * 3 + 2][k] ^= 1;
					}else{
						(k > l) ? (snum = l, bignum = k) : (snum = k, bignum = l);
						//????????????
						temp_ma[j * 3 + 2][block_size + (2 * block_size - 1 - snum) * snum / 2 + (bignum - snum - 1)] ^= 1;
					}
				}
			}
		}
	}

	// 
	for (size_t i = 0; i < block_size; i++)
	{
		for (size_t j = 0; j < block_size + block_size * block_size / 2 + 1; j++)
		{
			Ma[i][j] = temp_ma[i][j];
		}
	}
}



void creat_inv_s_matrix(char Ma[block_size][block_size + block_size * block_size / 2 + 1])
{
	char temp_ma[block_size][block_size + block_size * block_size / 2 + 1] = { 0 };
	//z2 z1 z0 \\z2 z1 z0
	for (size_t j = 0; j < block_size / 3; j++)
	{
		char z2[block_size + block_size * block_size / 2 + 1] = { 0 }, z1[block_size + block_size * block_size / 2 + 1] = { 0 }, z0[block_size + block_size * block_size / 2 + 1] = { 0 };
		for (size_t k = 0; k < block_size + block_size * block_size / 2 + 1; k++)
		{
			z2[k] ^= Ma[j * 3][k];		//x2
			z1[k] ^= Ma[j * 3 + 1][k];	//x1
			z0[k] ^= Ma[j * 3 + 2][k];	//x0
		}

		//x2 = z1 + z2 + z0z1
		for (size_t k = 0; k < block_size + block_size * block_size / 2 + 1; k++)
		{
			temp_ma[j * 3][k] ^= Ma[j * 3][k];//x2
			temp_ma[j * 3][k] ^= Ma[j * 3 + 1][k];//x2
		}
		//z2 = x2 + x0x1
		if (z0[block_size + block_size * block_size / 2] == 1)
		{
			for (size_t k = 0; k < block_size + block_size * block_size / 2; k++)
			{
				temp_ma[j * 3][k] ^= Ma[j * 3 + 1][k];//x1
			}
		}

		if (z1[block_size + block_size * block_size / 2] == 1)
		{
			for (size_t k = 0; k < block_size + block_size * block_size / 2; k++)
			{
				temp_ma[j * 3][k] ^= Ma[j * 3 + 2][k];//x0
			}
		}

		if (z1[block_size + block_size * block_size / 2] == 1 && z0[block_size + block_size * block_size / 2] == 1)
		{
			temp_ma[j * 3][block_size + block_size * block_size / 2] ^= 1;//x0
		}

		//z0*z1
		for (size_t k = 0; k < block_size; k++)
		{
			for (size_t l = 0; l < block_size; l++)
			{
				if (z0[k] == 1 && z1[l] == 1)
				{
					int snum, bignum;
					if (k == l)
					{
						temp_ma[j * 3][k] ^= 1;
					}
					else
					{
						(k > l) ? (snum = l, bignum = k) : (snum = k, bignum = l);
						temp_ma[j * 3][block_size + (2 * block_size - 1 - snum) * snum / 2 + (bignum - snum - 1)] ^= 1;
					}
				}
			}
		}

		//	x1 = z1 + z0z2
		if (z0[block_size + block_size * block_size / 2] == 1)
		{
			for (size_t k = 0; k < block_size + block_size * block_size / 2; k++)
			{
				temp_ma[j * 3 + 1][k] ^= Ma[j * 3][k];//x2
			}
		}

		if (z2[block_size + block_size * block_size / 2] == 1)
		{
			for (size_t k = 0; k < block_size + block_size * block_size / 2; k++)
			{
				temp_ma[j * 3 + 1][k] ^= Ma[j * 3 + 2][k];//x0
			}
		}

		if (z2[block_size + block_size * block_size / 2] == 1 && z0[block_size + block_size * block_size / 2] == 1)
		{
			temp_ma[j * 3 + 1][block_size + block_size * block_size / 2] ^= 1;//x0
		}
		for (size_t k = 0; k < block_size + block_size * block_size / 2 + 1; k++)
		{
			temp_ma[j * 3 + 1][k] ^= Ma[j * 3 + 1][k];//x1
		}

		// z0*z2
		for (size_t k = 0; k < block_size; k++)
		{
			for (size_t l = 0; l < block_size; l++)
			{
				if (z0[k] == 1 && z2[l] == 1)
				{
					int snum, bignum;
					if (k == l)
					{
						temp_ma[j * 3 + 1][k] ^= 1;
					}
					else
					{
						(k > l) ? (snum = l, bignum = k) : (snum = k, bignum = l);
						temp_ma[j * 3 + 1][block_size + (2 * block_size - 1 - snum) * snum / 2 + (bignum - snum - 1)] ^= 1;
					}
				}
			}
		}

		//	x0 = z0 + z1 + z2 + z1z2
		if (z1[block_size + block_size * block_size / 2] == 1)
		{
			for (size_t k = 0; k < block_size + block_size * block_size / 2; k++)
			{
				temp_ma[j * 3 + 2][k] ^= Ma[j * 3][k];//x2
			}
		}

		if (z2[block_size + block_size * block_size / 2] == 1)
		{
			for (size_t k = 0; k < block_size + block_size * block_size / 2; k++)
			{
				temp_ma[j * 3 + 2][k] ^= Ma[j * 3 + 1][k];//x0
			}
		}

		if (z2[block_size + block_size * block_size / 2] == 1 && z1[block_size + block_size * block_size / 2] == 1)
		{
			temp_ma[j * 3 + 2][block_size + block_size * block_size / 2] ^= 1;//x0
		}
		for (size_t k = 0; k < block_size + block_size * block_size / 2 + 1; k++)
		{
			temp_ma[j * 3 + 2][k] ^= Ma[j * 3 + 1][k];//x1
			temp_ma[j * 3 + 2][k] ^= Ma[j * 3][k];//x2
			temp_ma[j * 3 + 2][k] ^= Ma[j * 3 + 2][k];//x0
		}

		//z1*z2
		for (size_t k = 0; k < block_size; k++)
		{
			for (size_t l = 0; l < block_size; l++)
			{
				if (z1[k] == 1 && z2[l] == 1)
				{
					int snum, bignum;
					if (k == l)
					{
						temp_ma[j * 3 + 2][k] ^= 1;
					}
					else
					{
						(k > l) ? (snum = l, bignum = k) : (snum = k, bignum = l);
						temp_ma[j * 3 + 2][block_size + (2 * block_size - 1 - snum) * snum / 2 + (bignum - snum - 1)] ^= 1;
					}
				}
			}
		}
	}

	//
	for (size_t i = 0; i < block_size; i++)
	{
		for (size_t j = 0; j < block_size + block_size * block_size / 2 + 1; j++)
		{
			Ma[i][j] = temp_ma[i][j];
		}
	}
}

// inverse linear translation
void creat_invl_matrix(char Ma[block_size][block_size + block_size * block_size / 2 + 1], int round_k)
{
	char temp_ma[block_size][block_size + block_size * block_size / 2 + 1] = { 0 };

	for (size_t j = 0; j < block_size; j++)
	{
		for (size_t k = 0; k < block_size; k++)
		{
			if (inv_roundmatrix[round_k][j][k] == 1)
			{
				for (size_t i = 0; i < block_size + block_size * block_size / 2 + 1; i++)
				{
					temp_ma[j][i] ^= Ma[k][i];  //??j???��???��??????????k????????
				}
			}
		}
	}

	for (size_t j = 0; j < block_size; j++)
	{
		for (size_t k = 0; k < block_size + block_size * block_size / 2 + 1; k++)
		{
			Ma[j][k] = temp_ma[j][k];
		}
	}
}

// linear translation
void creat_l_matrix(char Ma[block_size][block_size + block_size * block_size / 2 + 1], int round_k)
{
	char temp_ma[block_size][block_size + block_size * block_size / 2 + 1] = { 0 };

	for (size_t j = 0; j < block_size; j++)
	{
		for (size_t k = 0; k < block_size; k++)
		{
			if (roundmatrix[round_k][j][k] == 1)
			{
				for (size_t i = 0; i < block_size + block_size * block_size / 2 + 1; i++)
				{
					temp_ma[j][i] ^= Ma[k][i];
				}
			}
		}
	}
	for (size_t j = 0; j < block_size; j++)
	{
		for (size_t k = 0; k < block_size + block_size * block_size / 2 + 1; k++)
		{
			Ma[j][k] = temp_ma[j][k];
		}
	}
}

char state_alpha[block_size][block_size + block_size * block_size / 2 + 1] = { 0 };
void juzhen_01()
{
	pre_work();

	// I Diagonal set
	for (size_t i = 0; i < block_size; i++)
	{
		state_alpha[i][i] ^= 1;
	}
	creat_s_matrix(state_alpha);
	creat_l_matrix(state_alpha, 0);
	
	// I
	for (size_t i = 0; i < block_size; i++)
	{
		state_alpha[i][i] ^= 1;
	}

	char state_beta[block_size][block_size + block_size * block_size / 2 + 1] = { 0 };

	for (size_t i = 0; i < block_size / 3; i++)
	{
		if (((cipher_glo[i] >> 2) & 1) == 1)
		{
			state_beta[i * 3][block_size + block_size * block_size / 2] = 1;
		}
		if (((cipher_glo[i] >> 1) & 1) == 1)
		{
			state_beta[i * 3 + 1][block_size + block_size * block_size / 2] = 1;
		}
		if (((cipher_glo[i]) & 1) == 1)
		{
			state_beta[i * 3 + 2][block_size + block_size * block_size / 2] = 1;
		}
	}

	// Diagonal set
	for (size_t i = 0; i < block_size; i++)
	{
		state_beta[i][i] ^= 1;
	}

	creat_invl_matrix(state_beta, 1);
	creat_inv_s_matrix(state_beta);

	for (size_t i = 0; i < block_size; i++)
	{
		for (size_t j = 0; j < block_size + block_size * block_size / 2 + 1; j++)
		{
			state_alpha[i][j] ^= state_beta[i][j];
		}
	}
}

void m2xcnf(){
    cout << "p cnf " << block_size + block_size * (block_size - 1) / 2 +1 << " " << block_size + 3 * block_size * (block_size - 1) / 2 +1 << endl;

    bool apd = 0;

	set<int> occur;

    for (int i = 0; i < block_size;i++){
        cout << "x";
        for (int j = 0; j < block_size + block_size * (block_size-1) / 2;j++){
            if(state_alpha[i][j]==1){
                cout << j + 1 << " ";
				if (j >= block_size && occur.find(j+1)==occur.end()){
					occur.insert(j + 1);
				}
			}
        }
        if(state_alpha[i][block_size+block_size*block_size/2]==0){
            cout << block_size + block_size * (block_size - 1) / 2 + 1<<" ";
            apd = 1;
        }
        cout << "0" << endl;
    }
    if (apd){
        cout<<block_size + block_size * (block_size - 1) / 2 + 1<<" 0"<<endl;
    }else{
        cout<<"-"<<block_size + block_size * (block_size - 1) / 2 + 1<<" 0"<<endl;
    }
    
    int k = block_size+1;
    for (int i = 1; i <= block_size - 1;i++){
        for (int j = i + 1;j<=block_size; j++)
        {
			if(occur.find(k)!=occur.end()){
				cout <<"-"<< k << " " << i << " 0" << endl;
            	cout <<"-"<< k << " " << j << " 0" << endl;	
            	cout << k << " -" << i <<" -" << j << " 0" << endl;
			}
            k++;
        }
    }
}

void m2cnf(){
	int num_v = block_size +(block_size + 1)*(block_size+ block_size * (block_size - 1) / 2);
	int num_c = 3 * block_size * (block_size - 1) / 2 +block_size*(2+4*(block_size+ block_size * (block_size - 1) / 2));
	cout << "p cnf " << num_v << " " << num_c << endl;

    bool apd = 0;
    
    int k = block_size+1;
    for (int i = 1; i <= block_size - 1;i++){
        for (int j = i + 1;j<=block_size; j++)
        {
			// if(occur.find(k)!=occur.end()){
				cout << i << " " <<"-"<< k << " 0"  << endl;
            	cout << j << " " <<"-"<< k << " 0" << endl;	
            	cout << k << " -" << i <<" -" << j << " 0" << endl;
			// }
            k++;
        }
    }

	for (int i = 0; i < block_size;i++){
		if(state_alpha[i][block_size+block_size*block_size/2]==1){
            cout<< k <<" 0"<<endl;
        }else{
			cout<<"-"<< k <<" 0"<<endl;
		}
        for (int j = 0; j < block_size + block_size * (block_size-1) / 2;j++){
            if(state_alpha[i][j]==1){
				//x=j+1,k=k,k'=k+1
				cout << j+1 << " -"<< k <<" "<< k+1 << " 0"  << endl;
				cout << "-"<< j+1 << " -"<< k <<" -"<< k+1 << " 0"  << endl;
				// cout << "-"<< j+1 << " "<< k << " 0"  << endl;
				// cout << "-"<< k+1 << " "<< k << " 0"  << endl;
                cout << j+1 << " -"<< k+1 <<" "<< k << " 0"  << endl;
				cout << "-"<< j+1 << " "<< k <<" "<< k+1 << " 0"  << endl;
				k++;
			}
        }   
		cout<<"-"<< k <<" 0"<<endl;
		k++;
    }
}

void writecnf(string filename){
	// ???????????
	filename.append(".cnf");
	ofstream outfile(filename, std::ofstream::out | std::ofstream::trunc);

	// ??????????  
    if (!outfile.is_open()) {  
        std::cerr << "Unable to open file for writing!" << std::endl;  
        return; 
    }  

	int num_v=block_size +(block_size + 1)*(block_size+ block_size * (block_size - 1) / 2);
	int num_c=3 * block_size * (block_size - 1) / 2 +block_size*(2+4*(block_size+ block_size * (block_size - 1) / 2));
	

    outfile << "p cnf " << num_v << " " << num_c << endl;

    bool apd = 0;
    
    int k = block_size+1;
    for (int i = 1; i <= block_size - 1;i++){
        for (int j = i + 1;j<=block_size; j++)
        {
			outfile <<"-"<< k << " " << i << " 0" << endl;
            outfile <<"-"<< k << " " << j << " 0" << endl;	
            outfile << k << " -" << i <<" -" << j << " 0" << endl;
			
            k++;
        }
    }

	for (int i = 0; i < block_size;i++){
		if(state_alpha[i][block_size+block_size*block_size/2]==1){
            outfile<< k <<" 0"<<endl;
        }else{
			outfile<<"-"<< k <<" 0"<<endl;
		}
        for (int j = 0; j < block_size + block_size * (block_size-1) / 2;j++){
            if(state_alpha[i][j]==1){
				//x=j+1,k=k,k'=k+1
				outfile << j+1 << " -"<< k <<" "<< k+1 << " 0"  << endl;
				outfile << "-"<< j+1 << " -"<< k <<" -"<< k+1 << " 0"  << endl;
                outfile << j+1 << " -"<< k+1 <<" "<< k << " 0"  << endl;
				outfile << "-"<< j+1 << " "<< k <<" "<< k+1 << " 0"  << endl;
				k++;
			}
        }   
		outfile<<"-"<< k <<" 0"<<endl;
		k++;
    }
	// ????  
    outfile.close();  
}

void writexnf(string filename){
	// ???????????
	filename.append(".xnf");
	ofstream outfile(filename, std::ofstream::out | std::ofstream::trunc);

	// ??????????  
    if (!outfile.is_open()) {  
        std::cerr << "Unable to open file for writing!" << std::endl;  
        return; 
    }  

    outfile << "p cnf " << block_size + block_size * (block_size - 1) / 2 +1 << " " << block_size + 3 * block_size * (block_size - 1) / 2 +1 << endl;

    bool apd = 0;

	set<int> occur;

    for (int i = 0; i < block_size;i++){
        outfile << "x";
        for (int j = 0; j < block_size + block_size * (block_size-1) / 2;j++){
            if(state_alpha[i][j]==1){
                outfile << j + 1 << " ";
				if (j >= block_size && occur.find(j+1)==occur.end()){
					occur.insert(j + 1);
				}
			}
        }
        if(state_alpha[i][block_size+block_size*block_size/2]==0){
            outfile << block_size + block_size * (block_size - 1) / 2 + 1<<" ";
            apd = 1;
        }
        outfile << "0" << endl;
    }
    if (apd){
        outfile<<block_size + block_size * (block_size - 1) / 2 + 1<<" 0"<<endl;
    }else{
        outfile<<"-"<<block_size + block_size * (block_size - 1) / 2 + 1<<" 0"<<endl;
    }
    
    int k = block_size+1;
    for (int i = 1; i <= block_size - 1;i++){
        for (int j = i + 1;j<=block_size; j++)
        {
			if(occur.find(k)!=occur.end()){
				outfile <<"-"<< k << " " << i << " 0" << endl;
            	outfile <<"-"<< k << " " << j << " 0" << endl;	
            	outfile << k << " -" << i <<" -" << j << " 0" << endl;
			}
            k++;
        }
    }
	// ????  
    outfile.close();  
}

int main(int argc, char* argv[])
{
	juzhen_01();

	string filename= argv[1];

	// auto start = chrono::high_resolution_clock::now();  
	// exhusted_search();
	// auto end = chrono::high_resolution_clock::now();  
  
    // // ??????????  
    // chrono::duration<double> duration = end - start;  
     
	// std::ofstream outputFile("compare2Exhusted.txt", std::ios::app); // ?????????  
    // if (outputFile.is_open()) {  
	// 	outputFile << "Exhusted: " << duration.count() << " s" << std::endl; 
    //     outputFile.close(); // ????  
    // } else {  
    //     // std::cout << "???????? 'output.txt'." << std::endl;  
    // }  

	for (size_t i = 0; i < block_size; i++)
	{
		for (size_t j = 0; j < block_size + block_size * block_size / 2 + 1; j++)
		{
			printf("%d ", state_alpha[i][j]);
		}
		printf(" \n");
	}

	writecnf(filename);
	// writexnf(filename);

	// m2cnf();
	// g++ equation2xcnf.cpp -o e2x.o
    // ./e2x.o >example27.cnf
    // minisat example27.cnf >result.log
}