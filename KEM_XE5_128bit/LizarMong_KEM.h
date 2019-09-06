#include <stdint.h>
#include "params.h"

#define iter 1000		// iteration number for keygen & EncDec test
#define testnum 100		// repeatetion number of Enc Dec procedure in a single iteration

#define sft (sizeof(size_t) * 4 - 1)

clock_t start, finish, elapsed1, elapsed2;

int Keygen(unsigned char *pk, unsigned char *sk);

int Enc(unsigned char *c, unsigned char *shared_k, const unsigned char *pk);

int Dec(unsigned char *shared_key, unsigned char *c, const unsigned char *sk, const unsigned char *pk);


