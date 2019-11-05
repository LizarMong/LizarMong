#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>

#include "LizarMong_KEM.h"

// add(rdtsc)//
uint64_t start_cycle1, finish_cycle1, start_cycle2, finish_cycle2, start_cycle3, finish_cycle3, cycles1, cycles2, cycles3;
long long rdtsc(void)
{
  unsigned long long result;
  __asm__ volatile(".byte 15;.byte 49;shlq $32,%%rdx;orq %%rdx,%%rax"
    : "=a" (result) ::  "%rdx");
  return result;
}
// add(rdtsc)//


unsigned char pk[SEED_LEN+LWE_N];
unsigned char sk[LWE_N+LWE_N/8];

void Keygen_ring() {
	int res = 0;
	start_cycle3 = rdtsc();
	for (int l = 0; l < iter; ++l) {
		res = Keygen(pk, sk);
		if (res == 3) {
			printf("Fault detected\n");
			break;
		}
	}
	start_cycle3 = rdtsc() - start_cycle3;

	printf("    Keygen Cycles:  %lu \n",  start_cycle3/iter);
}

void EncDecTest_RING() {
	// Set a messages
	unsigned char c[LWE_N+LWE_N+size_of_delta] = {0,};
	unsigned char shared_k1[MESSAGE_LEN] = {0,};
	unsigned char shared_k2[MESSAGE_LEN] = {0,};
	int i, res_enc = 0, res_dec = 0;

	cycles1=0;
	cycles2=0;
	elapsed1 = 0;
	elapsed2 = 0;

	for (int l = 0; l < iter; ++l) {
		for (i = 0; i < testnum; i++) {
			start_cycle1 = rdtsc();
			res_enc = Enc(c, shared_k1, pk);
			finish_cycle1 = rdtsc();
			cycles1 += (finish_cycle1 - start_cycle1);

			if (res_enc == 3) {
				printf("Fault detected\n");
				goto BREAK;
			}

			start_cycle2 = rdtsc();
			res_dec = Dec(shared_k2, c, sk, pk);
			finish_cycle2 = rdtsc();
			cycles2 += (finish_cycle2 - start_cycle2);
			
			if (res_dec == 1) {
				printf("    Decryption Validity Error Type 1 : c3 components\n");
				goto BREAK;
			}

			if (res_dec == 2) {
				printf("    Decryption Validity Error Type 2 : c2 components\n");
				goto BREAK;
			}
		}
		
		// Correctness check
		for (i = 0; i < MESSAGE_LEN; ++i) {
			if (shared_k1[i] != shared_k2[i]) {
				printf("Correctness Error, %d\n", i);
				break;
			}
		}
		if (i < MESSAGE_LEN) break;
	}
BREAK:
	printf("    Enc Cycles: %lu \n",  cycles1/iter/testnum);
	printf("    Dec Cycles: %lu \n",  cycles2/iter/testnum);
}


void main() {
	printf("\n  //////////////////////////////////////////////////////////////////\n\n");
	printf("\t\t"PARAMNAME" Parameter\n\n");
	printf("    LWE dimension: %d, \t\tLWR dimension: %d\n", LWE_N, LWE_N);
	printf("    Plaintext dimension: %d, \t\tPlaintext Modulus: %d bits\t\n", (MESSAGE_LEN*8), 1);
	printf("    Public Key modulus: %d bits, \tCiphertext modulus: %d bits\t\n\n", LOG_Q, LOG_P);
	printf("  //////////////////////////////////////////////////////////////////\n\n");
	printf("\t\t\tPerformance Test\n\n");

	// Key Generation
	Keygen_ring();

	// Enc and Dec
	EncDecTest_RING();
}
