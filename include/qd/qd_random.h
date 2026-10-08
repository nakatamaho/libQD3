/*
 * include/qd/qd_random.h
 *
 * Random number source shared by every libQD3 random function: ddrand,
 * tdrand, qdrand, eddrand, ds/ts/qs_real::rand, the c_*_rand C functions
 * and the Fortran random_number interfaces.
 *
 * The generator is xoshiro256** seeded through splitmix64.  It does not use
 * std::rand, so its width and its sequences are the same on every platform.
 * The initial state equals qd_srand(0x9E3779B97F4A7C15).  Calls are
 * serialized by a mutex.
 */
#ifndef _QD_QD_RANDOM_H
#define _QD_QD_RANDOM_H

#include <stdint.h>
#include <qd/qd_config.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Reseed the shared generator. */
QD_API void qd_srand(uint64_t seed);

/* Next 64 uniformly distributed random bits. */
QD_API uint64_t qd_rand_u64(void);

#ifdef __cplusplus
}
#endif

#endif /* _QD_QD_RANDOM_H */
