#include "Helper.h"

#include <immintrin.h>

unsigned char bit_width(unsigned long long x) 
{
    return x == 0 ? 1 : (unsigned char)(64 - _tzcnt_u64(x)); //intrinsic for 'count trailing zeroes'
}

// implementation for all unsigned integer types
unsigned int_sqrt(const unsigned long long n) 
{
    unsigned char shift = bit_width(n);
    shift += shift & 1; // round up to next multiple of 2
    unsigned result = 0;
    do 
    {
        shift -= 2;
        result <<= 1; // leftshift the result to make the next guess
        result |= 1;  // guess that the next bit is 1
        result ^= result * result > (n >> shift); // revert if guess too high
    } while (shift != 0);
    return result;
}

unsigned long long test(void)
{
	cTimer tmr; tmr.Start();

	unsigned long long sum = 0, sum1 = 0;
	for (unsigned long long i = 1; i < 1'00'000'000; i++)
	{
        sum += int_sqrt(i);
        //sum1 += sqrt(i);
    }
	tmr.Stop(true, "sqrt");
    std::cout << sum << " | " << sum1;
	return sum;
}