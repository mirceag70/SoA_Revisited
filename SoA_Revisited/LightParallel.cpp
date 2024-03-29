#include "Helper.h"

//uint8_t* sieve1; uint8_t* sieve5; uint8_t* sieve7; uint8_t* sieve11;
//bool GetBit (tpPrime n, uint8_t sieve[])
//{
//	tpPrime q = (n / 12) / 8, r = (n / 12) % 8;
//	return (sieve[q] & BIT_MASK[r]);
//};
//void ResetBit(tpPrime n, uint8_t sieve[])
//{
//	tpPrime q = (n / 12) / 8, r = (n / 12) % 8;
//	sieve[q] &= BIT_RESET_MASK[r];
//};
//void FlipBit(tpPrime n, uint8_t sieve[])
//{
//	tpPrime q = (n / 12) / 8, r = (n / 12) % 8;
//	sieve[q] ^= BIT_MASK[r];
//};

//void Pattern11(const tpPrime start, const tpPrime stop, uint8_t sieve[])
//{
//	const tpPrime xmax = (tpPrime)floor(sqrt(stop / 4));
//	for (tpPrime x = 1, jmp = 0; x <= xmax; x += (1 + jmp), jmp = 1 - jmp)
//	{
//		//get in position
//		tpPrime y, n0 = 4 * x * x;
//		if (n0 < start)
//		{
//			const tpPrime yy = (tpPrime)ceil(sqrt(start - n0));
//			y = 6 * (yy / 6) + 3;
//			if (y < yy) y += 6;
//		}
//		else y = 3;
//		//sieve
//		for (; ; y += 6)
//		{
//			tpPrime n = n0 + (y * y);
//			if (n > stop) break;
//			FlipBit(n, sieve);
//		}
//	}
//};
//auto Pattern12(const tpPrime start, const tpPrime stop, uint8_t sieve[])
//{
//	//for x = 0 we can not start with y = 1
//	uint8_t step = 2;
//	const tpPrime yy = (tpPrime)ceil(sqrt(start));
//	tpPrime y = 6 * (yy / 6) + 5;
//	if (y < yy)
//	{
//		y += 2;
//		if (y < yy) y += 4;
//		else step = 4;
//	}
//	for (; ; y += step, step = 6 - step)
//	{
//		tpPrime n = y * y;
//		if (n > stop) break;
//		FlipBit(n, sieve);
//	}
//	const tpPrime xmax = (tpPrime)floor(sqrt(stop / 4));
//	for (tpPrime x = 3; x <= xmax; x += 3)
//	{
//		//get in position
//		uint8_t step = 4;
//		tpPrime y, n0 = 4 * x * x;
//		if (n0 < start)
//		{
//			const tpPrime yy = (tpPrime)ceil(sqrt(start - n0));
//			y = 6 * (yy / 6) + 1;
//			if (y < yy)
//			{
//				y += 4;
//				if (y < yy) y += 2;
//				else step = 2;
//			}
//		}
//		else y = 1;
//		//sieve
//		for (; ; y += step, step = 6 - step)
//		{
//			tpPrime n = (4 * x * x) + (y * y);
//			if (n > stop) break;
//			FlipBit(n, sieve);
//		}
//	}
//};
//void Pattern5(const tpPrime start, const tpPrime stop, uint8_t sieve[])
//{
//	const tpPrime xmax = (tpPrime)floor(sqrt(stop / 4));
//	for (tpPrime x = 1, jmp = 0; x <= xmax; x += 1 + jmp, jmp = 1 - jmp)
//	{
//		//get in position
//		uint8_t step = 4;
//		tpPrime y, n0 = 4 * x * x;
//		if (n0 < start)
//		{
//			const tpPrime yy = (tpPrime)ceil(sqrt(start - n0));
//			y = 6 * (yy / 6) + 1;
//			if (y < yy)
//			{
//				y += 4;
//				if (y < yy) y += 2;
//				else step = 2;
//			}
//		}
//		else y = 1;
//		//sieve
//		for (; ; y += step, step = 6 - step)
//		{
//			tpPrime n = n0 + (y * y);
//			if (n > stop) break;
//			FlipBit(n, sieve);
//		}
//	}
//};
//void Pattern7(const tpPrime start, const tpPrime stop, uint8_t sieve[])
//{
//	const tpPrime xmax = (tpPrime)floor(sqrt(stop / 3));
//	for (tpPrime x = 1; x <= xmax; x += 2)
//	{
//		//get in position
//		uint8_t step = 2;
//		tpPrime y, n0 = 3 * x * x;
//		if (n0 < start)
//		{
//			const tpPrime yy = (tpPrime)ceil(sqrt(start - n0));
//			y = 6 * (yy / 6) + 2;
//			if (y < yy)
//			{
//				y += 2;
//				if (y < yy) y += 4;
//				else step = 4;
//			}
//		}
//		else y = 2;
//		//sieve			
//		for (; ; y += step, step = 6 - step)
//		{
//			tpPrime n = n0 + (y * y);
//			if (n > stop) break;
//			FlipBit(n, sieve);
//		}
//	}
//};
//void Pattern111(const tpPrime start, const tpPrime stop, uint8_t sieve[])
//{
//	const tpPrime xmax = (tpPrime)floor(sqrt(stop / 2));
//	tpPrime x = (tpPrime)ceil(sqrt(start / 3));
//	if (x & 1) x++;
//	for (; x <= xmax; x += 2)
//	{
//		//get in position
//		uint8_t step = 4;
//		tpPrime y, n0 = 3 * x * x;
//		if (n0 > stop)
//		{
//			const tpPrime yy = (tpPrime)ceil(sqrt(n0 - stop));
//			y = 6 * (yy / 6) + 1;
//			if (y < yy)
//			{
//				y += 4;
//				if (y < yy) y += 2;
//				else step = 2;
//			}
//		}
//		else y = 1;
//		//sieve			
//		for (; x > y; y += step, step = 6 - step)
//		{
//			tpPrime n = n0 - (y * y);
//			if (n < start) break;
//			FlipBit(n, sieve);
//		}
//	}
//};
//void Pattern112(const tpPrime start, const tpPrime stop, uint8_t sieve[])
//{
//	const tpPrime xmax = (tpPrime)floor(sqrt(stop / 2));
//	tpPrime x = (tpPrime)ceil(sqrt(start / 3));
//	if (not (x & 1)) x++;
//	for (; x <= xmax; x += 2)
//	{
//		//get in position
//		uint8_t step = 2;
//		tpPrime y, n0 = 3 * x * x;
//		if (n0 > stop)
//		{
//			const tpPrime yy = (tpPrime)ceil(sqrt(n0 - stop));
//			y = 6 * (yy / 6) + 2;
//			if (y < yy)
//			{
//				y += 2;
//				if (y < yy) y += 4;
//				else step = 4;
//			}
//		}
//		else y = 2;
//		//sieve			
//		for (; x > y; y += step, step = 6 - step)
//		{
//			tpPrime n = n0 - (y * y);
//			if (n < start) break;
//			FlipBit(n, sieve);
//		}
//	}
//};

template<typename FN1, typename FN2>
void Sieve2Patterns(FN1 fn1, FN2 fn2, const tpPrime limit, uint8_t sieve[])
{
	tpPrime iter_min = 0;
	if (limit > STEP_SIZE)
		for (; iter_min < (limit - STEP_SIZE); iter_min += STEP_SIZE)
		{
			fn1(iter_min, iter_min + STEP_SIZE, sieve);
			fn2(iter_min, iter_min + STEP_SIZE, sieve);
		}
	fn1(iter_min, limit, sieve); fn2(iter_min, limit, sieve);
};

template<typename FN>
void Sieve1Pattern(FN fn1, const tpPrime limit, uint8_t sieve[])
{
	tpPrime iter_min = 0;
	if (limit > STEP_SIZE)
		for (; iter_min < (limit - STEP_SIZE); iter_min += STEP_SIZE)
		{
			fn1(iter_min, iter_min + STEP_SIZE, sieve);
		}
	fn1(iter_min, limit, sieve);
};

void Sieve1(const tpPrime limit) 
{ 
	Sieve2Patterns(Pattern11, Pattern12, limit, sieve1); 
	//std::cout << "S1! "; 
}
void Sieve5(const tpPrime limit) 
{ 
	Sieve1Pattern(Pattern5, limit, sieve5); 
	//std::cout << "S5! "; 
}
void Sieve7(const tpPrime limit) 
{ 
	Sieve1Pattern(Pattern7, limit, sieve7); 
	//std::cout << "S7! ";
}
void Sieve11(const tpPrime limit)
{ 
	Sieve2Patterns(Pattern111, Pattern112, limit, sieve11); 
	//std::cout << "S11! ";
}

void MarkSquares(uint8_t work_sieve[], const tpPrime limit, const unsigned rem)
{
	auto MarkOneSquare = [&](uint8_t sieve[], const tpPrime n)
	{
		if (GetBit(n, sieve))
			for (tpPrime i = rem * n * n; i < limit; i += 12 * n * n)
				ResetBit(i, work_sieve);
	};

	tpPrime kmax = (tpPrime)floor(sqrt(limit));

	MarkOneSquare(sieve5, 5);
	MarkOneSquare(sieve7, 7);
	MarkOneSquare(sieve11, 11);
	for (tpPrime k = 12; k < kmax; k += 12)
	{
		MarkOneSquare(sieve1, k + 1);
		MarkOneSquare(sieve5, k + 5);
		MarkOneSquare(sieve7, k + 7);
		MarkOneSquare(sieve11, k + 11);
	}
};

//inc. tinkering 3
tpPrime SoA_I3(const tpPrime limit, uint8_t sv[], void*, void*)
{
	const tpPrime svlen = (limit / 96) + 2;
	uint8_t* sieve1 = sv;			  uint8_t* sieve5 = sieve1 + svlen;
	uint8_t* sieve7 = sieve5 + svlen; uint8_t* sieve11 = sieve7 + svlen;

	// 2 and 3 are not generated by the alghorithm
	tpPrime i, numprimes = 2; AddPrime(2); AddPrime(3);
	for (i = 0; i < 4 * svlen; i++) sv[i] = false;

auto Sieve2 = [&]<typename FN1, typename FN2>(FN1 fn1, FN2 fn2, uint8_t sieve[])
{
	tpPrime iter_min = 0;
	if (limit > STEP_SIZE)
		for (; iter_min < (limit - STEP_SIZE); iter_min += STEP_SIZE)
		{
			fn1(iter_min, iter_min + STEP_SIZE, sieve);
			fn2(iter_min, iter_min + STEP_SIZE, sieve);
		}
	fn1(iter_min, limit, sieve); fn2(iter_min, limit, sieve);
};
auto Sieve1 = [&]<typename FN>(FN fn1, uint8_t sieve[])
{
	tpPrime iter_min = 0;
	if (limit > STEP_SIZE)
		for (; iter_min < (limit - STEP_SIZE); iter_min += STEP_SIZE)
		{
			fn1(iter_min, iter_min + STEP_SIZE, sieve);
		}
	fn1(iter_min, limit, sieve);
};

cTimer tmr;
tmr.Start();

Sieve2(Pattern11, Pattern12, sieve1);
//tmr.LapTime(true, "ptrn1");

//5 mod 12
Sieve1(Pattern5, sieve5);
//tmr.LapTime(true, "ptrn5");

//7 mod 12
Sieve1(Pattern7, sieve7);
//tmr.LapTime(true, "ptrn7");

Sieve2(Pattern111, Pattern112, sieve11);

tmr.LapTime(true, "ptrn11");

// secondary loop - remove multiple of squares
auto MarkSquares = [&](uint8_t work_sieve[], const unsigned rem)
{
	auto MarkOneSquare = [&](uint8_t sieve[], const tpPrime n)
	{
		//tpPrime istp = 12 * n * n;
		if (GetBit(n, sieve))
			for (tpPrime i = rem * n * n; i < limit; i += 12 * n * n)
				ResetBit(i, work_sieve);
	};

	tpPrime kmax = (tpPrime)floor(sqrt(limit));

	//for (tpPrime k = 12; k < kmax; k += 12)
	//	MarkOneSquare(sieve1, k + 1);
	//for (tpPrime k = 0; k < kmax; k += 12)
	//	MarkOneSquare(sieve5, k + 5);
	//for (tpPrime k = 0; k < kmax; k += 12)
	//	MarkOneSquare(sieve7, k + 7);
	//for (tpPrime k = 0; k < kmax; k += 12)
	//	MarkOneSquare(sieve11, k + 11);

	MarkOneSquare(sieve5, 5);
	MarkOneSquare(sieve7, 7);
	MarkOneSquare(sieve11, 11);
	for (tpPrime k = 12; k < kmax; k += 12)
	{
		MarkOneSquare(sieve1, k + 1);
		MarkOneSquare(sieve5, k + 5);
		MarkOneSquare(sieve7, k + 7);
		MarkOneSquare(sieve11, k + 11);
	}
};
MarkSquares(sieve1, 1); MarkSquares(sieve5, 5);
MarkSquares(sieve7, 7); MarkSquares(sieve11, 11);
tmr.LapTime(true, "squares");

// last loop - get primes
for (tpPrime k = 0; k < svlen; k++)
{
	uint8_t sv1 = sieve1[k], sv5 = sieve5[k], sv7 = sieve7[k], sv11 = sieve11[k];
	auto CountPrime = [&](const tpPrime n) { numprimes++; AddPrime(n); };
	auto CountPrimes = [&](const unsigned b)
	{
		tpPrime n = 12 * (8 * k + b);
		if (sv1 & BIT_MASK[b]) CountPrime(n + 1);
		if (sv5 & BIT_MASK[b]) CountPrime(n + 5);
		if (sv7 & BIT_MASK[b]) CountPrime(n + 7);
		if (sv11 & BIT_MASK[b]) CountPrime(n + 11);
	};
	CountPrimes(0); CountPrimes(1); CountPrimes(2); CountPrimes(3);
	CountPrimes(4); CountPrimes(5); CountPrimes(6); CountPrimes(7);
}

tmr.LapTime(true, "counting");
tmr.Stop(/*true, "squares & counting"*/);
return numprimes;
}

void CountOneChunk(const tpPrime start, const tpPrime stop, tpPrime* numprimes)
{
	auto CountPrime = [&](const tpPrime n) { (*numprimes)++; AddPrimeThreaded(n); };
	for (tpPrime k = start; k < stop; k++)
	{
		uint8_t sv1 = sieve1[k], sv5 = sieve5[k], sv7 = sieve7[k], sv11 = sieve11[k];

		tpPrime n = 12 * (8 * k);
		if (sv1 & BIT_MASK[0]) CountPrime(n + 1);
		if (sv5 & BIT_MASK[0]) CountPrime(n + 5);
		if (sv7 & BIT_MASK[0]) CountPrime(n + 7);
		if (sv11 & BIT_MASK[0]) CountPrime(n + 11);
		n = 12 * (8 * k + 1);
		if (sv1 & BIT_MASK[1]) CountPrime(n + 1);
		if (sv5 & BIT_MASK[1]) CountPrime(n + 5);
		if (sv7 & BIT_MASK[1]) CountPrime(n + 7);
		if (sv11 & BIT_MASK[1]) CountPrime(n + 11);
		n = 12 * (8 * k + 2);
		if (sv1 & BIT_MASK[2]) CountPrime(n + 1);
		if (sv5 & BIT_MASK[2]) CountPrime(n + 5);
		if (sv7 & BIT_MASK[2]) CountPrime(n + 7);
		if (sv11 & BIT_MASK[2]) CountPrime(n + 11);
		n = 12 * (8 * k + 3);
		if (sv1 & BIT_MASK[3]) CountPrime(n + 1);
		if (sv5 & BIT_MASK[3]) CountPrime(n + 5);
		if (sv7 & BIT_MASK[3]) CountPrime(n + 7);
		if (sv11 & BIT_MASK[3]) CountPrime(n + 11);
		n = 12 * (8 * k + 4);
		if (sv1 & BIT_MASK[4]) CountPrime(n + 1);
		if (sv5 & BIT_MASK[4]) CountPrime(n + 5);
		if (sv7 & BIT_MASK[4]) CountPrime(n + 7);
		if (sv11 & BIT_MASK[4]) CountPrime(n + 11);
		n = 12 * (8 * k + 5);
		if (sv1 & BIT_MASK[5]) CountPrime(n + 1);
		if (sv5 & BIT_MASK[5]) CountPrime(n + 5);
		if (sv7 & BIT_MASK[5]) CountPrime(n + 7);
		if (sv11 & BIT_MASK[5]) CountPrime(n + 11);
		n = 12 * (8 * k + 6);
		if (sv1 & BIT_MASK[6]) CountPrime(n + 1);
		if (sv5 & BIT_MASK[6]) CountPrime(n + 5);
		if (sv7 & BIT_MASK[6]) CountPrime(n + 7);
		if (sv11 & BIT_MASK[6]) CountPrime(n + 11);
		n = 12 * (8 * k + 7);
		if (sv1 & BIT_MASK[7]) CountPrime(n + 1);
		if (sv5 & BIT_MASK[7]) CountPrime(n + 5);
		if (sv7 & BIT_MASK[7]) CountPrime(n + 7);
		if (sv11 & BIT_MASK[7]) CountPrime(n + 11);
	}
}

//void CountOneChunkX(const tpPrime start, const tpPrime stop, tpPrime* numprimes)
//{
//	auto CountPrime = [&](const tpPrime n) { (*numprimes)++; AddPrime(n); };
//	for (tpPrime k = start; k < stop; k++)
//	{
//		uint8_t sv1 = sieve1[k], sv5 = sieve5[k], sv7 = sieve7[k], sv11 = sieve11[k];
//		auto CountPrimes = [&](const unsigned b)
//		{
//			tpPrime n = 12 * (8 * k + b);
//			if (sv1 & BIT_MASK[b]) CountPrime(n + 1);
//			if (sv5 & BIT_MASK[b]) CountPrime(n + 5);
//			if (sv7 & BIT_MASK[b]) CountPrime(n + 7);
//			if (sv11 & BIT_MASK[b]) CountPrime(n + 11);
//			//if (sv1 & BIT_MASK[b]) { (*numprimes)++; AddPrime(n + 1); }; 
//			//if (sv5 & BIT_MASK[b]) { (*numprimes)++; AddPrime(n + 5); }; 
//			//if (sv7 & BIT_MASK[b]) { (*numprimes)++; AddPrime(n + 7); }; 
//			//if (sv11 & BIT_MASK[b]) { (*numprimes)++; AddPrime(n + 11); }; 
//		};
//		CountPrimes(0); CountPrimes(1);
//		CountPrimes(2); CountPrimes(3);
//		CountPrimes(4); CountPrimes(5);
//		CountPrimes(6); CountPrimes(7);
//	}
//}

//light parallel
tpPrime SoA_LP(const tpPrime limit, uint8_t sv[], void*, void*)
{
	const tpPrime svlen = (limit / 96) + 2;
	sieve1 = sv; sieve5 = sieve1 + svlen;
	sieve7 = sieve5 + svlen; sieve11 = sieve7 + svlen;

	// 2 and 3 are not generated by the alghorithm
	tpPrime i, numprimes = 2; AddPrime(2); AddPrime(3);
	for (i = 0; i < 4 * svlen; i++) sv[i] = false;

	cTimer tmr;
	tmr.Start();

	// sieve quadratics
	{
		std::thread t1 = std::thread(&Sieve1, limit);
		std::thread t5 = std::thread(&Sieve5, limit);
		std::thread t7 = std::thread(&Sieve7, limit);

		Sieve11(limit);

		t7.join(); t5.join(); t1.join();
	}
	tmr.LapTime(true, "patterns");

	// sieve multiple of squares
	{
		std::thread t1 = std::thread(&MarkSquares, sieve1, limit, 1);
		std::thread t5 = std::thread(&MarkSquares, sieve5, limit, 5);
		std::thread t7 = std::thread(&MarkSquares, sieve7, limit, 7);

		MarkSquares(sieve11, limit, 11);

		t7.join(); t5.join(); t1.join();
	}
	tmr.LapTime(true, "squares");

	// get primes

	//CountOneChunk(0, svlen / 4, &numprimes);
	//CountOneChunk(svlen / 4, svlen / 2, &numprimes);
	//CountOneChunk(svlen / 2, 3 * svlen / 4, &numprimes);
	//CountOneChunk(3 * svlen / 4, svlen, &numprimes);
	{
		tpPrime n1 = 0, n2 = 0, n3 = 0;
		std::thread t1 = std::thread(&CountOneChunk, 0, svlen / 4, &n1);
		std::thread t5 = std::thread(&CountOneChunk, svlen / 4, svlen / 2, &n2);
		std::thread t7 = std::thread(&CountOneChunk, svlen / 2, 3 * svlen / 4, &n3);
		CountOneChunk(3 * svlen / 4, svlen, &numprimes);
		t7.join(); t5.join(); t1.join();
		numprimes += n1 + n2 + n3;
	}

	//for (tpPrime k = 0; k < svlen; k++)
	//{
	//	uint8_t sv1 = sieve1[k], sv5 = sieve5[k], sv7 = sieve7[k], sv11 = sieve11[k];
	//	auto CountPrime = [&](const tpPrime n) { numprimes++; AddPrime(n); };
	//	auto CountPrimes = [&](const unsigned b)
	//	{
	//		tpPrime n = 12 * (8 * k + b);
	//		if (sv1 & BIT_MASK[b]) CountPrime(n + 1);
	//		if (sv5 & BIT_MASK[b]) CountPrime(n + 5);
	//		if (sv7 & BIT_MASK[b]) CountPrime(n + 7);
	//		if (sv11 & BIT_MASK[b]) CountPrime(n + 11);
	//	};
	//	CountPrimes(0); CountPrimes(1); CountPrimes(2); CountPrimes(3);
	//	CountPrimes(4); CountPrimes(5); CountPrimes(6); CountPrimes(7);
	//}

	tmr.LapTime(true, "counting");
	tmr.Stop(/*true, "squares & counting"*/);
	return numprimes;
}

uint8_t svrp[4 * (ROOT_LIMIT / 96 + 2)];
uint8_t root_primes_gap[NO_ROOT_PRIMES];
unsigned firstIdx[6] = {0, 125'793, 239'115/*113'322*/, 348'510/*109'395*/, 455'379/*106'869*/, 560'685/*105'306*/};

inline thread_local unsigned idx_last_prime_t = 0;
inline thread_local unsigned last_prime_t = 0;

inline void AddRootPrime(unsigned prime)
{
	root_primes_gap[idx_last_prime_t++] = (uint8_t)(prime - last_prime_t);
	last_prime_t = prime;
}
void ProcessOneChunk(const unsigned svlen, const unsigned nt, const unsigned i, unsigned* numprimes)
{
	const unsigned chunksz = svlen / nt;
	const unsigned start = i * chunksz;
	const unsigned stop = (i == (nt - 1)) ? svlen : start + chunksz;

	idx_last_prime_t = firstIdx[i];

	auto CountPrime = [&](const unsigned n) { (*numprimes)++; AddRootPrime(n); };
	for (unsigned k = start; k < stop; k++)
	{
		uint8_t sv1 = sieve1[k], sv5 = sieve5[k], sv7 = sieve7[k], sv11 = sieve11[k];

		unsigned n = 12 * (8 * k);
		if (sv1 & BIT_MASK[0]) CountPrime(n + 1);
		if (sv5 & BIT_MASK[0]) CountPrime(n + 5);
		if (sv7 & BIT_MASK[0]) CountPrime(n + 7);
		if (sv11 & BIT_MASK[0]) CountPrime(n + 11);
		n = 12 * (8 * k + 1);
		if (sv1 & BIT_MASK[1]) CountPrime(n + 1);
		if (sv5 & BIT_MASK[1]) CountPrime(n + 5);
		if (sv7 & BIT_MASK[1]) CountPrime(n + 7);
		if (sv11 & BIT_MASK[1]) CountPrime(n + 11);
		n = 12 * (8 * k + 2);
		if (sv1 & BIT_MASK[2]) CountPrime(n + 1);
		if (sv5 & BIT_MASK[2]) CountPrime(n + 5);
		if (sv7 & BIT_MASK[2]) CountPrime(n + 7);
		if (sv11 & BIT_MASK[2]) CountPrime(n + 11);
		n = 12 * (8 * k + 3);
		if (sv1 & BIT_MASK[3]) CountPrime(n + 1);
		if (sv5 & BIT_MASK[3]) CountPrime(n + 5);
		if (sv7 & BIT_MASK[3]) CountPrime(n + 7);
		if (sv11 & BIT_MASK[3]) CountPrime(n + 11);
		n = 12 * (8 * k + 4);
		if (sv1 & BIT_MASK[4]) CountPrime(n + 1);
		if (sv5 & BIT_MASK[4]) CountPrime(n + 5);
		if (sv7 & BIT_MASK[4]) CountPrime(n + 7);
		if (sv11 & BIT_MASK[4]) CountPrime(n + 11);
		n = 12 * (8 * k + 5);
		if (sv1 & BIT_MASK[5]) CountPrime(n + 1);
		if (sv5 & BIT_MASK[5]) CountPrime(n + 5);
		if (sv7 & BIT_MASK[5]) CountPrime(n + 7);
		if (sv11 & BIT_MASK[5]) CountPrime(n + 11);
		n = 12 * (8 * k + 6);
		if (sv1 & BIT_MASK[6]) CountPrime(n + 1);
		if (sv5 & BIT_MASK[6]) CountPrime(n + 5);
		if (sv7 & BIT_MASK[6]) CountPrime(n + 7);
		if (sv11 & BIT_MASK[6]) CountPrime(n + 11);
		n = 12 * (8 * k + 7);
		if (sv1 & BIT_MASK[7]) CountPrime(n + 1);
		if (sv5 & BIT_MASK[7]) CountPrime(n + 5);
		if (sv7 & BIT_MASK[7]) CountPrime(n + 7);
		if (sv11 & BIT_MASK[7]) CountPrime(n + 11);
	}
}

//root primes generator
void SoA_LP_gen_root_primes(void)
{
	const unsigned svlen = (ROOT_LIMIT / 96) + 2;
	sieve1 = svrp; sieve5 = sieve1 + svlen;
	sieve7 = sieve5 + svlen; sieve11 = sieve7 + svlen;

	// 2 and 3 are not generated by the alghorithm
	for (unsigned i = 0; i < 4 * svlen; i++) svrp[i] = false;

	// sieve quadratics
	{
		std::thread t1 = std::thread(&Sieve1, ROOT_LIMIT);
		std::thread t5 = std::thread(&Sieve5, ROOT_LIMIT);
		std::thread t7 = std::thread(&Sieve7, ROOT_LIMIT);
		Sieve11(ROOT_LIMIT);
		t7.join(); t5.join(); t1.join();
	}

	// sieve multiple of squares
	{
		std::thread t1 = std::thread(&MarkSquares, sieve1, ROOT_LIMIT, 1);
		std::thread t5 = std::thread(&MarkSquares, sieve5, ROOT_LIMIT, 5);
		std::thread t7 = std::thread(&MarkSquares, sieve7, ROOT_LIMIT, 7);
		MarkSquares(sieve11, ROOT_LIMIT, 11);
		t7.join(); t5.join(); t1.join();
	}

	// get primes
	unsigned numprimes = 0;
	{
		std::vector<std::thread> tc;
		const unsigned nt = 6; unsigned np[nt]{ 0 };
		for (auto i : Range<nt>())
			tc.push_back(std::thread(&ProcessOneChunk, svlen, nt, i, np + i));
		for (auto& t : tc) t.join();
		for (auto i : Range<nt>()) numprimes += np[i];
		assert(numprimes == NO_ROOT_PRIMES);
		//unsigned n1 = 0, n2 = 0, n3 = 0;
		//std::thread t1 = std::thread(&ProcessOneChunk, svlen, 4, 0, &n1);
		//std::thread t2 = std::thread(&ProcessOneChunk, svlen, 4, 1, &n2);
		//std::thread t3 = std::thread(&ProcessOneChunk, svlen, 4, 2, &n3);
		//							  ProcessOneChunk( svlen, 4, 3, &numprimes);
		//t3.join(); t2.join(); t1.join();
		//numprimes += n1 + n2 + n3;
	}

	//fix values
	root_primes_gap[firstIdx[1]] = 30;
	root_primes_gap[firstIdx[2]] = 2;
	root_primes_gap[firstIdx[3]] = 36;
	root_primes_gap[firstIdx[4]] = 62;
	root_primes_gap[firstIdx[5]] = 30;

	// check root primes
	//unsigned last_p = 0; AddPrime(2); AddPrime(3);
	//for (unsigned i = 0; i < NO_ROOT_PRIMES; i++)
	//{
	//	unsigned p = last_p + root_primes_gap[i];
	//	AddPrime(p);
	//	last_p = p;
	//}
}
