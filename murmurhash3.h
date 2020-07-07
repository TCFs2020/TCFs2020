// MurmurHash3 was written by Austin Appleby, and is placed in the public
// domain. The author hereby disclaims copyright to this source code.

#ifndef MURMURHASH3_H_INCLUDED
#define MURMURHASH3_H_INCLUDED

#include <stdint.h>

// MurmurHash3 produces 32-bit hash values on x86 platform.
uint32_t MurmurHash3_32(const void * key, int len, uint32_t seed);

// MurmurHash2 produces 64-bit hash values on x86 platform.
uint64_t MurmurHash2_64(const void * key, int len, uint64_t seed);

// MurmurHash3 produces 128-bit hash values on x86 platform.
void MurmurHash3_128(const void * key, int len, uint32_t seed, void * out);

#endif // MURMURHASH3_H_INCLUDED
