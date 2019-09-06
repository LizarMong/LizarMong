/*
 * Copyright (c) 2018, PQShield
 * Markku-Juhani O. Saarinen
 */

//  Generic prototypes for error correction code

#ifndef _XEF_H_
#define _XEF_H_

#include <stdint.h>
#include <stddef.h>

//  Parametrized versions. f = 0..5, number of errors fixed

//  Computes the parity code, XORs it at the end of payload
//  len = payload (bytes). Returns (payload | xef) length in *bits*.
size_t xef_compute(void *block, size_t len, unsigned f);

//  Fixes errors based on parity code. Call xef_compute() first to get delta.
//  len = payload (bytes). Returns (payload | xef) length in *bits*.
size_t xef_fixerr(void *block, size_t len, unsigned f);


// specific code from optimized implementations

void xe5_234_compute(void *block);
void xe5_234_fixerr(void *block);

#endif /* _XEF_H_ */
