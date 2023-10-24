#define ITERATIONS	1

#include "Helper.h"

void Try357(const uint64_t limit);

tpPrime SoA0(const tpPrime limit, bool sieve[], void*, void*);
tpPrime SoA1(const tpPrime limit, bool sieve[], void*, void*);
tpPrime SoA2(const tpPrime limit, uint8_t sieve[], void*, void*);
//tpPrime SoA3(const tpPrime limit, uint8_t sieve[], void*, void*);

tpPrime SoA_T(const tpPrime limit, bool sieve[], void*, void*);
tpPrime SoA_T1(const tpPrime limit, uint8_t sieve[], void*, void*);
tpPrime SoA_T2(const tpPrime limit, uint8_t sieve[], void*, void*);
tpPrime SoA_T3(const tpPrime limit, uint8_t sieve[], void*, void*);
tpPrime SoA_T4(const tpPrime limit, uint8_t sieve[], void*, void*);

tpPrime SoA_I(const tpPrime limit, uint8_t sieve[], void*, void*);
tpPrime SoA_I1(const tpPrime limit, uint8_t sieve[], void*, void*);
tpPrime SoA_I2(const tpPrime limit, uint8_t sieve[], void*, void*);
tpPrime SoA_I3(const tpPrime limit, uint8_t sieve[], void*, void*);

tpPrime SoA_LP(const tpPrime limit, uint8_t sieve[], void*, void*);
tpPrime SoA_P(const tpPrime limit, uint8_t sieve[], void*, void*);
tpPrime SoA_FP(const tpPrime limit, uint8_t sieve[], void*, void*);

tpPrime SoA_S(const tpPrime limit, void*, void*, void*);

void SoA_Interval(void);

//constexpr tpPrime LIMIT = 69'000'000'000;
//constexpr tpPrime LIMIT = 100'000'000'000'000;

constexpr tpPrime LIMIT = 1'000'000'000;

unsigned long long test(void);

int main()
{
    std::locale mylocale("");   // get global locale 
    std::cout.imbue(mylocale);  // imbue global locale for thousands delimiter

    //test(); return 0;
    
    Try357(LIMIT);

    Try_Sieve<bool, int, int, LIMIT + 1>
        (LIMIT, "SoA Naive", &SoA0, false);
    Try_Sieve<bool, int, int, LIMIT + 1>
        (LIMIT, "SoA Standard", &SoA1, false);
    Try_Sieve<uint8_t, int, int, LIMIT / 24 + 1>
        (LIMIT, "SoA 1bit 6k", &SoA2, false);
    Try_Sieve<bool, int, int, LIMIT / 2 + 1>
        (LIMIT, "SoA w. tuples", &SoA_T, false);
    Try_Sieve<uint8_t, int, int, LIMIT / 24 + 1>
        (LIMIT, "SoA w. tuples - 1bit6k", &SoA_T1, false);
    Try_Sieve<uint8_t, int, int, LIMIT / 24 + 1>
        (LIMIT, "SoA w. pattern - 1bit6k", &SoA_T2, false);    
    Try_Sieve<uint8_t, int, int, 4 * (LIMIT / 96 + 2)>
        (LIMIT, "SoA w. pattern - 4 buffers", &SoA_T3, false);
    Try_Sieve<uint8_t, int, int, 4 * (LIMIT / 96 + 2)>
        (LIMIT, "SoA w. 4 buffers - optimized", &SoA_T4, false);
    Try_Sieve<uint8_t, int, int, 4 * (LIMIT / 96 + 2)>
        (LIMIT, "SoA incremental", &SoA_I, false);
    Try_Sieve<uint8_t, int, int, 4 * (LIMIT / 96 + 2)>
        (LIMIT, "SoA incr. - tink", &SoA_I1, false);
    Try_Sieve<uint8_t, int, int, 4 * (LIMIT / 96 + 2)>
        (LIMIT, "SoA incr. - tink2", &SoA_I2, false);
    Try_Sieve<uint8_t, int, int, 4 * (LIMIT / 96 + 2)>
        (LIMIT, "SoA incr. - tink3", &SoA_I3, false);
    Try_Sieve<uint8_t, int, int, 4 * (LIMIT / 96 + 2)>
        (LIMIT, "SoA light parallel", &SoA_LP, false);
    Try_Sieve<uint8_t, int, int, 4 * (LIMIT / 96 + 2)>
        (LIMIT, "SoA parallel", &SoA_P, false);
    Try_Sieve<uint8_t, int, int, 4 * (LIMIT / 96 + 2)>
        (LIMIT, "SoA full parallel", &SoA_FP, false);
    Try_Sieve<int, int, int>
        (LIMIT, "SoA segmented", &SoA_S, false);

    //SoA_Interval();
}

