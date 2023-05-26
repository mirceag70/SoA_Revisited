#include "Helper.h"

//extern uint8_t* sieve1; extern uint8_t* sieve5; extern uint8_t* sieve7; extern uint8_t* sieve11;
//bool GetBit(tpPrime n, uint8_t sieve[])
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

void PMarkSquares(uint8_t work_sieve[], const tpPrime limit, const unsigned rem)
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

void SieveQuadraticsChunk1(const tpPrime limit, const unsigned nt, const unsigned i)
{
	//cTimer tmr;	tmr.Start();
	tpPrime chunksz = limit / nt;
	tpPrime start = i * chunksz;
	tpPrime stop = (i == (nt-1)) ? limit : start + chunksz;
	Sieve2Patterns(Pattern11, Pattern12, start, stop, sieve1);
	//tmr.Stop(true, "1");
}
void SieveQuadraticsChunk5(const tpPrime limit, const unsigned nt, const unsigned i)
{
	//cTimer tmr;	tmr.Start();
	tpPrime chunksz = limit / nt;
	tpPrime start = i * chunksz;
	tpPrime stop = (i == (nt - 1)) ? limit : start + chunksz;
	Sieve1Pattern(Pattern5, start, stop, sieve5);
	//tmr.Stop(true, "5");
}
void SieveQuadraticsChunk7(const tpPrime limit, const unsigned nt, const unsigned i)
{
	//cTimer tmr;	tmr.Start();
	tpPrime chunksz = limit / nt;
	tpPrime start = i * chunksz;
	tpPrime stop = (i == (nt - 1)) ? limit : start + chunksz;
	Sieve1Pattern(Pattern7, start, stop, sieve7);
	//tmr.Stop(true, "7");
}
void SieveQuadraticsChunk11(const tpPrime limit, const unsigned nt, const unsigned i)
{
	//cTimer tmr;	tmr.Start();
	tpPrime chunksz = limit / nt;
	tpPrime start = i * chunksz;
	tpPrime stop = (i == (nt - 1)) ? limit : start + chunksz;
	Sieve2Patterns(Pattern111, Pattern112, start, stop, sieve11);
	//tmr.Stop(true, "11");
}

void ParallelSieve(const tpPrime limit)
{
	std::vector<std::thread> t1, t5, t7, t11;

	const unsigned nt = 8; 
	const unsigned nt1 = nt + 0; 
	const unsigned nt5 = nt + 0;
	const unsigned nt7 = nt + 0; 
	const unsigned nt11 = nt + 0;

	for (auto i : Range<nt1>())
		t1.push_back(std::thread(&SieveQuadraticsChunk1, limit, nt1, i));
	for (auto i : Range<nt5>())
		t5.push_back(std::thread(&SieveQuadraticsChunk5, limit, nt5, i));
	for (auto i : Range<nt7>())
		t7.push_back(std::thread(&SieveQuadraticsChunk7, limit, nt7, i));
	for (auto i : Range<nt11>())
		t11.push_back(std::thread(&SieveQuadraticsChunk11, limit, nt11, i));

	for (auto& t : t11) t.join();
	for (auto& t : t7) t.join();
	for (auto& t : t5) t.join();
	for (auto& t : t1) t.join();
}

void PCountOneChunk(const tpPrime svlen, const unsigned nt, const unsigned i, tpPrime *np)
{
	const tpPrime chunksz = svlen / nt;
	const tpPrime start = i * chunksz;
	const tpPrime stop = (i == (nt - 1)) ? svlen : start + chunksz;
	*np = 0;
	auto CountPrime = [&](const tpPrime n) { (*np)++; AddPrimeThreaded(n); };
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

tpPrime ParallelCount(const tpPrime svlen)
{
	std::vector<std::thread> tc;
	const unsigned nt = 400; tpPrime np[nt]{0};
	for (auto i : Range<nt>())
		tc.push_back(std::thread(&PCountOneChunk, svlen, nt, i, np+i));
	for (auto& t : tc) t.join();
	tpPrime numprimes = 0;
	for (auto i : Range<nt>()) numprimes += np[i];
	return numprimes;
}

//CPU parallel
tpPrime SoA_P(const tpPrime limit, uint8_t sv[], void*, void*)
{
	const tpPrime svlen = (limit / 96) + 2;
	sieve1 = sv; sieve5 = sieve1 + svlen;
	sieve7 = sieve5 + svlen; sieve11 = sieve7 + svlen;

	// 2 and 3 are not generated by the alghorithm
	for (tpPrime i = 0; i < 4 * svlen; i++) sv[i] = false;

	cTimer tmr;
	tmr.Start();

	// sieve quadratics
	ParallelSieve(limit);
	tmr.LapTime(/*true, "patterns"*/);

	// sieve multiple of squares
	{
		std::thread t1 = std::thread(&PMarkSquares, sieve1, limit, 1);
		std::thread t5 = std::thread(&PMarkSquares, sieve5, limit, 5);
		std::thread t7 = std::thread(&PMarkSquares, sieve7, limit, 7);

		PMarkSquares(sieve11, limit, 11);

		t7.join(); t5.join(); t1.join();
	}
	tmr.LapTime(true, "squares");

	// get primes
	tpPrime numprimes = 2 + ParallelCount(svlen);

	//AddPrime(2); AddPrime(3);
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

	tmr.LapTime(/*true, "counting"*/);
	tmr.Stop(/*true, "squares & counting"*/);
	return numprimes;
}

void SoA_LP_gen_root_primes(void);

extern uint8_t root_primes_gap[];

void RPMarkSquares(uint8_t work_sieve[], const tpPrime limit, const unsigned rem)
{
	tpPrime last_p = 0;
	for (unsigned i = 0; i < NO_ROOT_PRIMES; i++)
	{
		tpPrime p = last_p + root_primes_gap[i];
		tpPrime n = rem * p * p;
		if (n > limit) break;
		for (; n <= limit; n += 12 * p * p)
			ResetBit(n, work_sieve);
		last_p = p;
	}
};

void RPMarkSquaresChunk(uint8_t work_sieve[], const unsigned rem,
						const tpPrime start, const tpPrime stop)
{
	tpPrime last_p = 0;
	for (unsigned i = 0; i < NO_ROOT_PRIMES; i++)
	{
		tpPrime p = last_p + root_primes_gap[i]; 
		tpPrime psqr = p * p;
		tpPrime n = rem * psqr;
		if (n < start)
			n = (12 * ((start - n) / (12 * psqr)) + rem) * psqr;
		if (n > stop) 
			break;
		while (n < start) n += 12 * psqr;
		for (; n <= stop; n += 12 * psqr)
			ResetBit(n, work_sieve);
		last_p = p;
	}
}

void RPMarkSquaresS(uint8_t sieve[], const unsigned rem,
					const tpPrime limit, const unsigned nt, const unsigned i)
{
	//cTimer tmr;
	//tmr.Start();
	tpPrime chunksz = limit / nt;
	tpPrime start = i * chunksz;
	tpPrime stop = (i == (nt - 1)) ? limit : start + chunksz;

	tpPrime iter_min = start;
	if ((stop - start) > STEP_SIZE)
		for (; iter_min < (stop - STEP_SIZE); iter_min += STEP_SIZE)
			RPMarkSquaresChunk(sieve, rem, iter_min, iter_min + STEP_SIZE);
	RPMarkSquaresChunk(sieve, rem, iter_min, stop);
	//tmr.Stop(true, (rem==1) ? "1": ((rem == 5) ? "5": ((rem == 7) ? "7" : "11")));
};

void ParallelMark(const tpPrime limit)
{
	std::vector<std::thread> t1, t5, t7, t11;

	const unsigned nt = 6;
	const unsigned nt1 = nt + 0;
	const unsigned nt5 = nt + 0;
	const unsigned nt7 = nt + 0;
	const unsigned nt11 = nt + 0;

	for (auto i : Range<nt1>())
		t1.push_back(std::thread(&RPMarkSquaresS, sieve1, 1, limit, nt1, i));
	for (auto i : Range<nt5>())
		t5.push_back(std::thread(&RPMarkSquaresS, sieve5, 5, limit, nt5, i));
	for (auto i : Range<nt7>())
		t7.push_back(std::thread(&RPMarkSquaresS, sieve7, 7, limit, nt7, i));
	for (auto i : Range<nt11>())
		t11.push_back(std::thread(&RPMarkSquaresS, sieve11, 11, limit, nt11, i));

	for (auto& t : t11) t.join();
	for (auto& t : t7) t.join();
	for (auto& t : t5) t.join();
	for (auto& t : t1) t.join();
}

//full parallel
tpPrime SoA_FP(const tpPrime limit, uint8_t sv[], void*, void*)
{
	cTimer tmr;
	tmr.Start();
	SoA_LP_gen_root_primes();

	const tpPrime svlen = (limit / 96) + 2;
	sieve1 = sv; sieve5 = sieve1 + svlen;
	sieve7 = sieve5 + svlen; sieve11 = sieve7 + svlen;
	for (tpPrime i = 0; i < 4 * svlen; i++) sv[i] = false;

	//tmr.LapTime(/*true, "init"*/);


	// sieve quadratics
	ParallelSieve(limit);
	//tmr.LapTime(true, "patterns");

	// sieve multiple of squares
	ParallelMark(limit);
	//{
	//	std::thread t11 = std::thread(&RPMarkSquaresS, sieve1, 1, limit, 2, 0);
	//	std::thread t12 = std::thread(&RPMarkSquaresS, sieve1, 1, limit, 2, 1);
	//	std::thread t51 = std::thread(&RPMarkSquaresS, sieve5, 5, limit, 2, 0);
	//	std::thread t52 = std::thread(&RPMarkSquaresS, sieve5, 5, limit, 2, 1);
	//	std::thread t7 = std::thread(&RPMarkSquaresS, sieve7, 7, limit, 1, 0);

	//	RPMarkSquaresS(sieve11, 11, limit, 1, 0);

	//	t7.join(); t52.join(); t51.join(); t12.join(); t11.join();
	//}
	//tmr.LapTime(true, "squares");

	// get primes
	tpPrime numprimes = 2 + ParallelCount(svlen);

	//AddPrime(2); AddPrime(3);
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

	//tmr.LapTime(true, "counting");
	tmr.Stop(/*true, "squares & counting"*/);
	return numprimes;
}
