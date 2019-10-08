#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>

#include "LizarMong_KEM.h"
#include "randombytes.h"
#include "fips202.h"
#include "xef.h"


int Keygen(unsigned char *pk, unsigned char *sk){
	unsigned char pk_a[LWE_N];
	unsigned char pk_b[LWE_N*2]={0,};
	unsigned char seed_a[SEED_LEN];
	uint16_t sk_s[HS];
	int i, j;

//////// Gen poly a //////// 
	randombytes(seed_a, SEED_LEN);	
	shake256(pk_a, LWE_N, seed_a, SEED_LEN);

//////// gen poly s ////////
	memset(sk, 0, LWE_N+LWE_N/8); // difference PKE Algorithm because random u.
	unsigned char seed_s[HS*4];
	unsigned int sk_random_idx;
	int hw=0, count = 0;

	randombytes(seed_s, HS*4);

	while (hw < HS) {
		sk_random_idx = seed_s[count++]; 
		sk_random_idx <<= 8;
		sk_random_idx ^= seed_s[count++];
		sk_random_idx &= (LWE_N - 1); // (seed1) || (seed2) 
		if (sk[sk_random_idx] == 0) {
			sk[sk_random_idx] = (seed_s[count++] & 0x02) - 1;
			hw++;
		}
		if (count >= HS*4 - 3) {
			randombytes(seed_s, HS*4);
			count = 0;
		}
	}
	
	if (hw != HS) { // fault detecting
		return 3;
	}

//////// gen s_idx //////// 
	int neg_start = 0, back_position = HS;	
	for (i = 0; i < LWE_N; ++i) {
		if (sk[i] == 0x01){sk_s[neg_start++] = i;}
		if (sk[i] == 0xff){sk_s[--back_position] = i;}
	}

/* ADD compare PKE algorithm */

//////// gen u //////// 
	unsigned char u[LWE_N/8];
	randombytes(u, LWE_N/8);

//////// concat sk = sk || u //////// 
	for (i=LWE_N; i<LWE_N+LWE_N/8; ++i){sk[i] = u[i-LWE_N];}

/*		END of the ADD		*/

//////// Initialize b as an error polynomial e //////// 
	unsigned char b0, b1, tmp2[LWE_N/4];
	randombytes(tmp2,LWE_N/4);	// tmp2[0]'s 0, 1bit = pk_b[0], tmp2[0]'s 2,3bit = pk_b[1]and so on
	for(j=0; j<LWE_N/4; ++j){ 	// Centered Binomial Distribution
		b0 = tmp2[j] & 0x01;
		tmp2[j] = tmp2[j] >> 1;
		b1 = tmp2[j] & 0x01;
		pk_b[j*4+0] = b0 -b1;
		tmp2[j] = tmp2[j] >> 1;
		b0 = tmp2[j] & 0x01;
		tmp2[j] = tmp2[j] >> 1;
		b1 = tmp2[j] & 0x01;
		pk_b[j*4+1] = b0 -b1;
		tmp2[j] = tmp2[j] >> 1;
		b0 = tmp2[j] & 0x01;
		tmp2[j] = tmp2[j] >> 1;
		b1 = tmp2[j] & 0x01;
		pk_b[j*4+2] = b0 -b1;
		tmp2[j] = tmp2[j] >> 1;
		b0 = tmp2[j] & 0x01;
		tmp2[j] = tmp2[j] >> 1;
		b1 = tmp2[j] & 0x01;
		pk_b[j*4+3] = b0 -b1;
		tmp2[j] = tmp2[j] >> 1;		
	}
	
	if (j != (LWE_N/4)) { // fault detecting
		return 3;
	}

//////// mult a*s and add e(pk_b) //////// 
	unsigned char startindex; // random byte for SCA countermeasure
	randombytes(&startindex, 1);
	for (i = 0; i < HS; ++i) {
		uint16_t deg = sk_s[(startindex + i) & 0x7F];
		uint16_t branch = (2 * ((((startindex + i) & 0x7F) - neg_start) >> sft & 0x1) - 1);
		for (int j = 0; j < LWE_N; ++j) {pk_b[deg + j] -= branch * pk_a[j];}
	}
	for (j = 0; j < LWE_N; ++j) {pk_b[j] -= pk_b[LWE_N + j];}

//////// Concat seed_genA || pk_b ////////
	memcpy(pk, seed_a, SEED_LEN);
	memcpy(pk+SEED_LEN, pk_b, LWE_N);

	return 0;
}



int Enc(unsigned char *c, unsigned char *shared_k, const unsigned char *pk){ // PKE's m = shared_k
	int i, j;
	unsigned char c1h_a[LWE_N*2]={0,};
	unsigned char c1h_b[LWE_N*2]={0,};	
	unsigned char *hash = NULL;
	unsigned char *hash_t = NULL; // difference PKE Algorithm because shared_k hash function.


//////// Generate a random polynomial delta //////// 
	unsigned char delta[size_of_delta] = { 0, };
	int sum = 0;
	randombytes(delta, size_of_delta);
	for (i = 0; i < size_of_delta; ++i) {
		sum += delta[i];
	}
	if (sum == 0) { // fault detecting
		return 3;
	}

//////// Set d = H'(delta) and Concat c //////// 
	hash = calloc(size_of_delta, sizeof(unsigned char));
	shake256(hash, size_of_delta, delta, size_of_delta);	// H'(delta) is LWE_N bit = LWE_N/8 bytes.
	memcpy(c+LWE_N+LWE_N, hash, size_of_delta);			// Concat C with d

//////// Set r = H(delta) //////// 
	unsigned char r[LWE_N]={0,};
	uint16_t r_idx[HR];
	unsigned int r_random_idx; 
	int hw=0, count = 0;

	hash = calloc(HR*4, sizeof(unsigned char));
	shake256(hash, HR*4, delta, size_of_delta);		// **Check1 같은 delta로 가능??**

	while (hw < HR) {
		r_random_idx = hash[count++]; 
		r_random_idx <<= 8;	
		r_random_idx ^= hash[count++];
		r_random_idx = r_random_idx & (LWE_N - 1); // (seed1) || (seed2) 
		if (r[r_random_idx] == 0) {
			r[r_random_idx] = (hash[count++] & 0x02) - 1;
			hw++;
		}
		if (count >= (HR*4 - 3)) { 
			printf("New seed hash!\n");
			shake256(hash, HR*4, hash, HR*4);
			count = 0;
		}
	}

//////// Generate r_idx //////// 
	int neg_start = 0, back_position = HR;

	for (i = 0; i < LWE_N; ++i) {
		if (r[i] == 0x01){r_idx[neg_start++] = i;}
		else if (r[i] == 0xff){r_idx[--back_position] = i;}
	}

//////// Encoding delta using Error Correcting Code //////// 
	unsigned char delta_hat[LWE_N / 8]={0,}; // delta_hat is result of Encoding(delta)

	memcpy(delta_hat, delta, size_of_delta);
	xe5_234_compute(delta_hat);				 // 256bit delta + 234bit redundancy = 490bit (512-490 = 22bit is already padding in code_line number 157)

//////// Parse seed_a||pk_b from pk and Make pk_a //////// 
	unsigned char pk_a[LWE_N];
	unsigned char pk_b[LWE_N];
	unsigned char seed_a[SEED_LEN];

	memcpy(seed_a, pk, SEED_LEN);
	shake256(pk_a, LWE_N, seed_a, SEED_LEN);

	memcpy(pk_b, pk+SEED_LEN, LWE_N);

//////// Initialize c1h_b as q/2 * delta_hat ////////  		
	for (i = 0; i < LWE_N / 8; ++i) {
		for (j = 0; j < 8; ++j) {
			c1h_b[8 * i + j] = (delta_hat[i] >> j) << _8_LOG_T;	//Shift each bit of delta_hat to MSB and save in c1h_b
		}
	}


//////// Compute a * r and b * r, and then add to c1h_a and c1h_b, respectively. //////// 
	unsigned char startindex; // random byte for SCA countermeasure
	randombytes(&startindex, 1);
	for(i = 0; i < HR; ++i){
		uint16_t branch = (2 * ((((startindex + i) & 0x7F) - neg_start) >> sft & 0x1) - 1);
		uint16_t deg = r_idx[(startindex + i) & 0x7F];
		for(j = 0; j < LWE_N; ++j){
			c1h_a[deg+j] += branch * pk_a[j];
			c1h_b[deg+j] += branch * pk_b[j];
		}
	}
	for(j = 0; j < LWE_N; ++j){
		c1h_a[j] -= c1h_a[LWE_N+j];
		c1h_b[j] -= c1h_b[LWE_N+j];
	}

//////// Send c1h_a and c1h_b from mod q to mod p and mod k //////// 

	for (i=0; i< LWE_N; ++i) {
		c[i] = ((c1h_a[i] + 0x02) & 0xfc);
		c[LWE_N + i] = ((c1h_b[i] + 0x08) & 0xf0);
	}

/* ADD compare PKE algorithm */

//////// G(c1, d, delta) //////// 
	hash_t=calloc((LWE_N+LWE_N+size_of_delta+size_of_delta), sizeof(unsigned char));
	memcpy(hash_t, c, (LWE_N+LWE_N+size_of_delta));	// c = c1||d
	memcpy(hash_t+(LWE_N+LWE_N+size_of_delta), delta, size_of_delta);
	shake256(shared_k, MESSAGE_LEN, hash_t, LWE_N+LWE_N+size_of_delta+size_of_delta);
	free(hash_t);

/*		END of the ADD		*/

	free(hash);
	return 0;
}




int Dec(unsigned char *shared_k, unsigned char *c, const unsigned char *sk, const unsigned char *pk){
	int res = 0;
	int i, j;

	unsigned char *hash = NULL;
	unsigned char *hash_t = NULL;

	unsigned char delta_hat[LWE_N / 8]={0,};
	unsigned char delta[size_of_delta]={0,};
	unsigned char c1h_a[LWE_N*2] = { 0, };
	unsigned char c1h_b[LWE_N*2] = { 0, };
	unsigned char decomp_delta[LWE_N*2]={0,};

//////// Initialize delta as c1h_b //////// 
	memcpy(decomp_delta, c+LWE_N, LWE_N);
	memcpy(c1h_a, c, LWE_N);

//////// Omit the task of changing mod_k to mod_p. because the data_type is unsigned_char. //////// 


//////// gen s_idx ////////
	uint16_t sk_s[HS];
	int neg_start = 0, back_position = HS;	
	for (i = 0; i < LWE_N; ++i) {
		if (sk[i] == 0x01){sk_s[neg_start++] = i;}
		if (sk[i] == 0xff){sk_s[--back_position] = i;}
	}

//////// Compute delta (delta(= c1b) + c1(=c1a) * s) //////// 
	unsigned char startindex; // random byte for SCA countermeasure
	randombytes(&startindex, 1);
	for(i = 0; i < HS; ++i){
		uint16_t branch = (2 * ((((startindex + i) & 0x7F) - neg_start) >> sft & 0x1) - 1);
		uint16_t deg = sk_s[(startindex + i) & 0x7F];
		for(int j = 0; j < LWE_N; ++j){
			decomp_delta[deg+j] += branch * c1h_a[j];
	    }
	}
	for(j = 0; j < LWE_N; ++j){
		decomp_delta[j] -= decomp_delta[LWE_N+j];
	}

//////// Compute delta = 2/p * delta //////// 
	for (i = 0; i < LWE_N; ++i) {
		decomp_delta[i] += 0x40;
		decomp_delta[i] >>= _8_LOG_T;
	}
//////// Set delta_hat //////// 
	for (i = 0; i < LWE_N/8; ++i) {
		for (j = 0; j < 8; ++j) {
			uint8_t a = (decomp_delta[8 * i + j]) << j;
			delta_hat[i] ^= a;  
		}
	}

//////// Decoding delta_hat using Error Correcting Code //////// 
	xe5_234_compute(delta_hat);
	xe5_234_fixerr(delta_hat);
	memcpy(delta, delta_hat, size_of_delta);

//////// Set d = H'(delta) //////// 
	//hash = calloc(size_of_delta, sizeof(unsigned char));
	unsigned char d[size_of_delta]={0,};
	shake256(d, size_of_delta, delta, size_of_delta); // H'(delta) is LWE_N bit = LWE_N/8 bytes.
	//memcpy(d, hash, LWE_N/8);

//////// check d==d //////// 
	for (i = 0; i < size_of_delta; ++i) {
		//printf("%04x, %04x %d \n",c[LWE_N+LWE_N+i],d[i],i);
		if (c[LWE_N+LWE_N+i] != d[i]){
			printf("delta error!!!!");		

			// G(c1, d, u)
			hash_t=calloc((LWE_N+LWE_N+size_of_delta+LWE_N/8), sizeof(unsigned char));
			memcpy(hash_t, c, LWE_N+LWE_N+size_of_delta);
			memcpy(hash_t+LWE_N+LWE_N+size_of_delta, sk+LWE_N, LWE_N/8);

			shake256(shared_k, MESSAGE_LEN, hash_t, LWE_N+LWE_N+size_of_delta+LWE_N/8);

			free(hash);
			free(hash_t);
			return res=1;
		}
	}

//////// Set r = H(delta) //////// 
	unsigned char r[LWE_N] = { 0, };
	unsigned int r_idx[HR];
	unsigned int r_random_idx; 
	int hw=0, count = 0;

	hash = calloc(HR*4, sizeof(unsigned char));
	shake256(hash, HR*4, delta, size_of_delta);	//check2 **같은 delta로 가능??**

	while (hw < HR) {
		r_random_idx = hash[count++]; 
		r_random_idx <<= 8;	
		r_random_idx ^= hash[count++];
		r_random_idx = r_random_idx & (LWE_N - 1); // (seed1) || (seed2) 
		if (r[r_random_idx] == 0) {
			r[r_random_idx] = (hash[count++] & 0x02) - 1;
			hw++;
		}
		if (count >= (HR*4 - 3)) { 
			printf("DEC_New seed hash!\n");
			shake256(hash, HR*4, delta, size_of_delta);/////////////
			count = 0;
		}
	}
	
//////// Generate r_idx //////// 
	neg_start = 0;
	back_position = HR;

	for (i = 0; i < LWE_N; ++i) {
		if (r[i] == 0x01){r_idx[neg_start++] = i;}
		else if (r[i] == 0xff){r_idx[--back_position] = i;}
	}

//////// Encoding delta using Error Correcting Code ////////
	memset(delta_hat, 0, LWE_N/8);
	memcpy(delta_hat, delta, size_of_delta);
	xe5_234_compute(delta_hat);


//////// Parse seed_a||pk_b from pk & Make pk_a //////// 
	unsigned char pk_a[LWE_N];
	unsigned char pk_b[LWE_N];
	unsigned char seed_a[SEED_LEN];

	memcpy(seed_a, pk, SEED_LEN);
	shake256(pk_a, LWE_N, seed_a, SEED_LEN);

	memcpy(pk_b, pk+SEED_LEN, LWE_N);

//////// Initialize c1h_b as q/2 * delta_hat ////////  	
	for (i = 0; i < LWE_N / 8; ++i) {
		for (j = 0; j < 8; ++j) {
			c1h_b[8 * i + j] = (delta_hat[i] >> j) << _8_LOG_T;// Shift each bit of delta_hat to MSB and save in c1h_b
		}
	}

//////// Compute a * r and b * r, and then add to c2h_a and c2h_b, respectively. //////// 
	memset(c1h_a, 0, LWE_N*2);
	randombytes(&startindex, 1);
	for(i = 0; i < HR; ++i){
		uint16_t branch = (2 * ((((startindex + i) & 0x7F) - neg_start) >> sft & 0x1) - 1);
		uint16_t deg = r_idx[(startindex + i) & 0x7F];
		for(j = 0; j < LWE_N; ++j){
			c1h_a[deg+j] += branch * pk_a[j];
			c1h_b[deg+j] += branch * pk_b[j];
		}
	}
	for(j = 0; j < LWE_N; ++j){
		c1h_a[j] -= c1h_a[LWE_N+j];
		c1h_b[j] -= c1h_b[LWE_N+j];
	}

//////// Send c2h_a and c2h_b from mod q to mod p and mod k //////// 
	for (i = 0; i < LWE_N; ++i) {
		c1h_a[i] = (c1h_a[i] + 0x02) & 0xfc;
		c1h_b[i] = (c1h_b[i] + 0x08) & 0xf0;
	}

//////// check c == c_hat ? //////// 

	for (i = 0; i < LWE_N; ++i) {
		if ((c[i] != ((c1h_a[i] + 0x02) & 0xfc)) || (c[LWE_N + i] != ((c1h_b[i] + 0x08) & 0xf0))){
			printf("c1h_a or c1h_b error!!!!");		

			// G(c1, d, u)
			hash_t=calloc((LWE_N+LWE_N+size_of_delta+LWE_N/8), sizeof(unsigned char));
			memcpy(hash_t, c, LWE_N+LWE_N+size_of_delta);
			memcpy(hash_t+LWE_N+LWE_N+size_of_delta, sk+LWE_N, LWE_N/8);

			shake256(shared_k, MESSAGE_LEN, hash_t, LWE_N+LWE_N+size_of_delta+LWE_N/8);

			free(hash_t);
			free(hash);
			return res = 2;
		}
	}

//////// G(c1, d, delta) //////// 
	hash_t=calloc((LWE_N+LWE_N+size_of_delta+size_of_delta), sizeof(unsigned char));
	memcpy(hash_t, c, LWE_N+LWE_N+size_of_delta);
	memcpy(hash_t+LWE_N+LWE_N+size_of_delta, delta, size_of_delta);
	shake256(shared_k, MESSAGE_LEN, hash_t, LWE_N+LWE_N+size_of_delta+size_of_delta);

	free(hash_t);
	free(hash);

	return res;
}
