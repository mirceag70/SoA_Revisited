#include "Helper.h"

constexpr tpPrime SEGMENT_SIZE = 10'000'000'000ull;
constexpr tpPrime SVLEN = 2 + SEGMENT_SIZE/96;
uint8_t sv1[SVLEN], sv5[SVLEN], sv7[SVLEN], sv11[SVLEN];

tpPrime segment_start = 0;

//sievers

#define COUNTERS
#ifdef COUNTERS 
tpPrime flips = 0, positionings = 0;
#endif

void Pattern11S(const tpPrime start, const tpPrime stop, uint8_t sieve[])
{
	const tpPrime nmax = stop - segment_start;
	const tpPrime xmax = (tpPrime)(sqrt(stop / 4));
	for (tpPrime x = 1, jmp = 0; x <= xmax; x += (1 + jmp), jmp = 1 - jmp)
	{
		//get in position
		tpPrime y, n0 = 4 * x * x;
		if (n0 < start)
		{
			const tpPrime yy = (tpPrime)ceil(sqrt(start - n0));
			y = 6 * (yy / 6) + 3;
			if (y < yy) y += 6;
		}
		else y = 3;
#ifdef COUNTERS 
		positionings++;
#endif
		//sieve
		for (tpPrime n = n0 + (y * y) - segment_start; n <= nmax; n += 12 * y + 36, y += 6)
		{
			assert(n >= start - segment_start);
			FlipBit(n, sieve); 
#ifdef COUNTERS 
			flips++; 
#endif
		}		
		//tpPrime n = n0 + (y * y);
		//assert(n >= start);
		//n -= segment_start;
		//for (y = 0; n < nmax; n += 12 * y + 36, y += 6)
		//{
		//	assert(n >= start - segment_start);
		//	FlipBit(n, sieve);
		//}
	}
}; 
void Pattern12S(const tpPrime start, const tpPrime stop, uint8_t sieve[])
{
	const tpPrime nmax = stop - segment_start;
	const tpPrime xmax = (tpPrime)(sqrt(stop / 4));
	{//for x = 0 we can not start with y = 1
		uint8_t step = 2;
		const tpPrime yy = (tpPrime)ceil(sqrt(start));
		tpPrime y = 6 * (yy / 6) + 5;
		if (y < yy)
		{
			y += 2;
			if (y < yy) y += 4;
			else step = 4;
		}
#ifdef COUNTERS 
		positionings++;
#endif
		for (tpPrime n = y * y - segment_start; n <= nmax; n += step * (2 * y + step), y += step, step = 6 - step)
		{
			assert(n >= start - segment_start);
			FlipBit(n, sieve);
#ifdef COUNTERS 
			flips++;
#endif
		}
		//	for (; ; y += step, step = 6 - step)
		//	{
		//		tpPrime n = y * y;
		//		if (n > stop) break;
		//		assert(n >= start);
		//		FlipBit(n - segment_start, sieve);
		//	}	
	}
	for (tpPrime x = 3; x <= xmax; x += 3)
	{
		//get in position
		uint8_t step = 4;
		tpPrime y, n0 = 4 * x * x;
		if (n0 < start)
		{
			const tpPrime yy = (tpPrime)ceil(sqrt(start - n0));
			y = 6 * (yy / 6) + 1;
			if (y < yy)
			{
				y += 4;
				if (y < yy) y += 2;
				else step = 2;
			}
		}
		else y = 1;
#ifdef COUNTERS 
		positionings++;
#endif
		//sieve
		for (tpPrime n = n0 + (y * y) - segment_start; n <= nmax; n += step * (2 * y + step), y += step, step = 6 - step)
		{
			assert(n >= start - segment_start);
			FlipBit(n, sieve);
#ifdef COUNTERS 
			flips++;
#endif
		}
		//		for (; ; y += step, step = 6 - step)
		//		{
		//			tpPrime n = (4 * x * x) + (y * y);
		//			if (n > stop) break;
		//			FlipBit(n - segment_start, sieve);
		//		}	
	}
};
void SieveQuadraticsChunk1S(const tpPrime sstart, const tpPrime sstop, const unsigned nt, const unsigned i)
{
	//cTimer tmr;	tmr.Start();
	tpPrime chunksz = (sstop - sstart) / nt;
	tpPrime start = sstart + i * chunksz;
	tpPrime stop = (i == (nt - 1)) ? sstop : start + chunksz;
	Sieve2Patterns(Pattern11S, Pattern12S, start, stop, sieve1);
	//tmr.Stop(true, "1");
}

void Pattern5S(const tpPrime start, const tpPrime stop, uint8_t sieve[])
{
	const tpPrime nmax = stop - segment_start;
	const tpPrime xmax = (tpPrime)(sqrt(stop / 4));
	for (tpPrime x = 1, jmp = 0; x <= xmax; x += 1 + jmp, jmp = 1 - jmp)
	{
		//get in position
		uint8_t step = 4;
		tpPrime y, n0 = 4 * x * x;
		if (n0 < start)
		{
			const tpPrime yy = (tpPrime)ceil(sqrt(start - n0));
			y = 6 * (yy / 6) + 1;
			if (y < yy)
			{
				y += 4;
				if (y < yy) y += 2;
				else step = 2;
			}
		}
		else y = 1;
#ifdef COUNTERS 
		positionings++;
#endif
		//sieve
		for (tpPrime n = n0 + (y * y) - segment_start; n <= nmax; n += step * (2 * y + step), y += step, step = 6 - step)
		{
			assert(n >= start - segment_start);
			FlipBit(n, sieve);
#ifdef COUNTERS 
			flips++;
#endif
		}
		//for (; ; y += step, step = 6 - step)
		//{
		//	tpPrime n = n0 + (y * y);
		//	if (n > stop) break;
		//	assert(n >= segment_start);
		//	FlipBit(n - segment_start, sieve);
		//}
	}
};
void SieveQuadraticsChunk5S(const tpPrime sstart, const tpPrime sstop, const unsigned nt, const unsigned i)
{
	//cTimer tmr;	tmr.Start();
	tpPrime chunksz = (sstop - sstart) / nt;
	tpPrime start = sstart + i * chunksz;
	tpPrime stop = (i == (nt - 1)) ? sstop : start + chunksz;
	Sieve1Pattern(Pattern5S, start, stop, sieve5);
	//tmr.Stop(true, "5");
}

void Pattern7S(const tpPrime start, const tpPrime stop, uint8_t sieve[])
{
	const tpPrime nmax = stop - segment_start;
	const tpPrime xmax = (tpPrime)(sqrt(stop / 3));
	for (tpPrime x = 1; x <= xmax; x += 2)
	{
		//get in position
		uint8_t step = 2;
		tpPrime y, n0 = 3 * x * x;
		if (n0 < start)
		{
			const tpPrime yy = (tpPrime)ceil(sqrt(start - n0));
			y = 6 * (yy / 6) + 2;
			if (y < yy)
			{
				y += 2;
				if (y < yy) y += 4;
				else step = 4;
			}
		}
		else y = 2;
#ifdef COUNTERS 
		positionings++;
#endif
		//sieve			
		for (tpPrime n = n0 + (y * y) - segment_start; n <= nmax; n += step * (2 * y + step), y += step, step = 6 - step)
		{
			assert(n >= start - segment_start);
			FlipBit(n, sieve);
#ifdef COUNTERS 
			flips++;
#endif
		}
		//for (; ; y += step, step = 6 - step)
		//{
		//	tpPrime n = n0 + (y * y);
		//	if (n > stop) break;
		//	assert(n >= segment_start);
		//	FlipBit(n - segment_start, sieve);
		//}
	}
};
void SieveQuadraticsChunk7S(const tpPrime sstart, const tpPrime sstop, const unsigned nt, const unsigned i)
{
	//cTimer tmr;	tmr.Start();
	tpPrime chunksz = (sstop - sstart) / nt;
	tpPrime start = sstart + i * chunksz;
	tpPrime stop = (i == (nt - 1)) ? sstop : start + chunksz;
	Sieve1Pattern(Pattern7S, start, stop, sieve7);
	//tmr.Stop(true, "7");
}

void Pattern111S(const tpPrime start, const tpPrime stop, uint8_t sieve[])
{
	const tpPrime nmin = start - segment_start;
	const tpPrime xmax = (tpPrime)(sqrt(stop / 2));
	tpPrime x = (tpPrime)(sqrt(start / 3));
	if (x & 1) x++;
	for (; x <= xmax; x += 2)
	{
		//get in position
		uint8_t step = 4;
		tpPrime y, n0 = 3 * x * x;
		if (n0 > stop)
		{
			const tpPrime yy = (tpPrime)ceil(sqrt(n0 - stop));
			y = 6 * (yy / 6) + 1;
			if (y < yy)
			{
				y += 4;
				if (y < yy) y += 2;
				else step = 2;
			}
		}
		else y = 1;
#ifdef COUNTERS 
		positionings++;
#endif
		//sieve			
		n0 -= y * y;
		if(n0 >= start)
		for (tpPrime n = n0 - segment_start; (x > y) and (n >= nmin); y += step, step = 6 - step)
		{
			FlipBit(n, sieve);
#ifdef COUNTERS 
			flips++;
#endif
			tpPrime dn = step * (2 * y + step);
			if (n >= dn)
				n -= dn;
			else
				break;
		}
		//for (tpPrime n = n0 - (y * y); (x > y) and (n >= start); y += step, step = 6 - step)
		//{
		//	FlipBit(n - segment_start, sieve);
		//	n -= step * (2 * y + step);
		//}
	}
};
void Pattern112S(const tpPrime start, const tpPrime stop, uint8_t sieve[])
{
	const tpPrime nmin = start - segment_start;
	const tpPrime xmax = (tpPrime)(sqrt(stop / 2));
	tpPrime x = (tpPrime)(sqrt(start / 3));
	if (not (x & 1)) x++;
	for (; x <= xmax; x += 2)
	{
		//get in position
		uint8_t step = 2;
		tpPrime y, n0 = 3 * x * x;
		if (n0 > stop)
		{
			const tpPrime yy = (tpPrime)ceil(sqrt(n0 - stop));
			y = 6 * (yy / 6) + 2;
			if (y < yy)
			{
				y += 2;
				if (y < yy) y += 4;
				else step = 4;
			}
		}
		else y = 2;
#ifdef COUNTERS 
		positionings++;
#endif
		//sieve			
		n0 -= y * y;
		if (n0 >= start)
		for (tpPrime n = n0 - segment_start; (x > y) and (n >= nmin); y += step, step = 6 - step)
		{
			FlipBit(n, sieve);
#ifdef COUNTERS 
			flips++;
#endif
			tpPrime dn = step * (2 * y + step);
			if (n >= dn)
				n -= dn;
			else
				break;
		}
		//for (; x > y; y += step, step = 6 - step)
		//{
		//	tpPrime n = n0 - (y * y);
		//	if (n < start) break;
		//	assert(n >= segment_start);
		//	FlipBit(n - segment_start, sieve);
		//}
	}
};
void SieveQuadraticsChunk11S(const tpPrime sstart, const tpPrime sstop, const unsigned nt, const unsigned i)
{
	//cTimer tmr;	tmr.Start();
	tpPrime chunksz = (sstop - sstart) / nt;
	tpPrime start = sstart + i * chunksz;
	tpPrime stop = (i == (nt - 1)) ? sstop : start + chunksz;
	Sieve2Patterns(Pattern111S, Pattern112S, start, stop, sieve11);
	//tmr.Stop(true, "11");
}

void SegmentSieve(const tpPrime start, const tpPrime stop)
{
	std::vector<std::thread> t1, t5, t7, t11;

	const unsigned nt = 8;
	const unsigned nt1 = nt + 0;
	const unsigned nt5 = nt + 0;
	const unsigned nt7 = nt + 0;
	const unsigned nt11 = nt + 0;

	for (auto i : Range<nt1>())
		t1.push_back(std::thread(&SieveQuadraticsChunk1S, start, stop, nt1, i));
	for (auto i : Range<nt5>())
		t5.push_back(std::thread(&SieveQuadraticsChunk5S, start, stop, nt5, i));
	for (auto i : Range<nt7>())
		t7.push_back(std::thread(&SieveQuadraticsChunk7S, start, stop, nt7, i));
	for (auto i : Range<nt11>())
		t11.push_back(std::thread(&SieveQuadraticsChunk11S, start, stop, nt11, i));

	for (auto& t : t11) t.join();
	for (auto& t : t7) t.join();
	for (auto& t : t5) t.join();
	for (auto& t : t1) t.join();
}

//markers
extern uint8_t root_primes_gap[];

void MarkSquaresChunk(uint8_t sieve[], const unsigned rem,
	const tpPrime start, const tpPrime stop)
{
	tpPrime last_p = 0;
	for (unsigned i = 0; i < NO_ROOT_PRIMES; i++)
	{
		tpPrime p = last_p + root_primes_gap[i];
		tpPrime psqr = p * p;
		tpPrime n = rem * psqr;
		if (n < start)
			n = (12 * (((start - n) / psqr) / 12) + rem) * psqr;
		if (n > stop) 
			break;
		while (n < start) n += 12 * psqr;
		for (; n <= stop; n += 12 * psqr)
			ResetBit(n - segment_start, sieve);
		last_p = p;
	}
}

void MarkSquaresS(uint8_t sieve[], const unsigned rem,
	const tpPrime sstart, const tpPrime sstop, const unsigned nt, const unsigned i)
{
	//cTimer tmr;
	//tmr.Start();
	tpPrime chunksz = (sstop - sstart) / nt;
	tpPrime start = sstart + i * chunksz;
	tpPrime stop = (i == (nt - 1)) ? sstop : start + chunksz;

	tpPrime iter_min = start;
	if ((stop - start) > STEP_SIZE)
		for (; iter_min < (stop - STEP_SIZE); iter_min += STEP_SIZE)
			MarkSquaresChunk(sieve, rem, iter_min, iter_min + STEP_SIZE);
	MarkSquaresChunk(sieve, rem, iter_min, stop);
	//tmr.Stop(true, (rem==1) ? "1": ((rem == 5) ? "5": ((rem == 7) ? "7" : "11")));
};

void SegmentMark(const tpPrime start, const tpPrime stop)
{
	std::vector<std::thread> t1, t5, t7, t11;

	const unsigned nt = 6;
	const unsigned nt1 = nt + 0;
	const unsigned nt5 = nt + 0;
	const unsigned nt7 = nt + 0;
	const unsigned nt11 = nt + 0;

	for (auto i : Range<nt1>())
		t1.push_back(std::thread(&MarkSquaresS, sieve1, 1, start, stop, nt1, i));
	for (auto i : Range<nt5>())
		t5.push_back(std::thread(&MarkSquaresS, sieve5, 5, start, stop, nt5, i));
	for (auto i : Range<nt7>())
		t7.push_back(std::thread(&MarkSquaresS, sieve7, 7, start, stop, nt7, i));
	for (auto i : Range<nt11>())
		t11.push_back(std::thread(&MarkSquaresS, sieve11, 11, start, stop, nt11, i));

	for (auto& t : t11) t.join();
	for (auto& t : t7) t.join();
	for (auto& t : t5) t.join();
	for (auto& t : t1) t.join();
}

//counters
void CountOneChunk(const unsigned nt, const unsigned i, tpPrime* np)
{
	const tpPrime chunksz = SVLEN / nt;
	const tpPrime start = i * chunksz;
	const tpPrime stop = (i == (nt - 1)) ? SVLEN : start + chunksz;
	*np = 0;
	auto CountPrime = [&](const tpPrime n) { (*np)++; AddPrimeThreaded(n); };
	for (tpPrime k = start; k < stop; k++)
	{
		uint8_t sv1 = sieve1[k], sv5 = sieve5[k], sv7 = sieve7[k], sv11 = sieve11[k];

		tpPrime n = segment_start + 12 * (8 * k);
		if (sv1 & BIT_MASK[0]) CountPrime(n + 1);
		if (sv5 & BIT_MASK[0]) CountPrime(n + 5);
		if (sv7 & BIT_MASK[0]) CountPrime(n + 7);
		if (sv11 & BIT_MASK[0]) CountPrime(n + 11);
		n = segment_start + 12 * (8 * k + 1);
		if (sv1 & BIT_MASK[1]) CountPrime(n + 1);
		if (sv5 & BIT_MASK[1]) CountPrime(n + 5);
		if (sv7 & BIT_MASK[1]) CountPrime(n + 7);
		if (sv11 & BIT_MASK[1]) CountPrime(n + 11);
		n = segment_start + 12 * (8 * k + 2);
		if (sv1 & BIT_MASK[2]) CountPrime(n + 1);
		if (sv5 & BIT_MASK[2]) CountPrime(n + 5);
		if (sv7 & BIT_MASK[2]) CountPrime(n + 7);
		if (sv11 & BIT_MASK[2]) CountPrime(n + 11);
		n = segment_start + 12 * (8 * k + 3);
		if (sv1 & BIT_MASK[3]) CountPrime(n + 1);
		if (sv5 & BIT_MASK[3]) CountPrime(n + 5);
		if (sv7 & BIT_MASK[3]) CountPrime(n + 7);
		if (sv11 & BIT_MASK[3]) CountPrime(n + 11);
		n = segment_start + 12 * (8 * k + 4);
		if (sv1 & BIT_MASK[4]) CountPrime(n + 1);
		if (sv5 & BIT_MASK[4]) CountPrime(n + 5);
		if (sv7 & BIT_MASK[4]) CountPrime(n + 7);
		if (sv11 & BIT_MASK[4]) CountPrime(n + 11);
		n = segment_start + 12 * (8 * k + 5);
		if (sv1 & BIT_MASK[5]) CountPrime(n + 1);
		if (sv5 & BIT_MASK[5]) CountPrime(n + 5);
		if (sv7 & BIT_MASK[5]) CountPrime(n + 7);
		if (sv11 & BIT_MASK[5]) CountPrime(n + 11);
		n = segment_start + 12 * (8 * k + 6);
		if (sv1 & BIT_MASK[6]) CountPrime(n + 1);
		if (sv5 & BIT_MASK[6]) CountPrime(n + 5);
		if (sv7 & BIT_MASK[6]) CountPrime(n + 7);
		if (sv11 & BIT_MASK[6]) CountPrime(n + 11);
		n = segment_start + 12 * (8 * k + 7);
		if (sv1 & BIT_MASK[7]) CountPrime(n + 1);
		if (sv5 & BIT_MASK[7]) CountPrime(n + 5);
		if (sv7 & BIT_MASK[7]) CountPrime(n + 7);
		if (sv11 & BIT_MASK[7]) CountPrime(n + 11);
	}
}

tpPrime SegmentCount(const tpPrime start, const tpPrime stop)
{
	std::vector<std::thread> tc;
	const unsigned nt = 400; tpPrime np[nt]{ 0 };
	for (auto i : Range<nt>())
		tc.push_back(std::thread(&CountOneChunk, nt, i, np + i));
	for (auto& t : tc) t.join();
	tpPrime numprimes = 0;
	for (auto i : Range<nt>()) numprimes += np[i];
	return numprimes;
}

//segmented
tpPrime SoA_SieveSegment(const tpPrime start, const tpPrime stop)
{
	//cTimer tmr;
	//tmr.Start();
	nln(); 	
#ifdef COUNTERS
	flips = positionings = 0;
#endif
	memset(sieve1, 0, SVLEN); memset(sieve5, 0, SVLEN);
	memset(sieve7, 0, SVLEN); memset(sieve11, 0, SVLEN);
	//tmr.LapTime(true, "reset");

	// sieve quadratics
	SegmentSieve(start, stop);
	//tmr.LapTime(true, "patterns");

	// sieve multiple of squares
	SegmentMark(start, stop);
	//tmr.LapTime(true, "squares");

	// get primes
	tpPrime numprimes = 0;
	numprimes = SegmentCount(start, stop);

#ifdef COUNTERS
	std::cout << flips << " | " << positionings;
#endif
	//for (tpPrime k = 0; k < SVLEN; k++)
	//{
	//	uint8_t sv1 = sieve1[k], sv5 = sieve5[k], sv7 = sieve7[k], sv11 = sieve11[k];
	//	auto CountPrime = [&](const tpPrime n) { numprimes++; AddPrime(n); };
	//	auto CountPrimes = [&](const unsigned b)
	//	{
	//		tpPrime n = segment_start + 12 * (8 * k + b);
	//		if (sv1 & BIT_MASK[b]) CountPrime(n + 1);
	//		if (sv5 & BIT_MASK[b]) CountPrime(n + 5);
	//		if (sv7 & BIT_MASK[b]) CountPrime(n + 7);
	//		if (sv11 & BIT_MASK[b]) CountPrime(n + 11);
	//	};
	//	CountPrimes(0); CountPrimes(1); CountPrimes(2); CountPrimes(3);
	//	CountPrimes(4); CountPrimes(5); CountPrimes(6); CountPrimes(7);
	//}

	//tmr.LapTime(true, "counting");
	//tmr.Stop(/*true, "squares & counting"*/);
	// std::cout << ".";
	return numprimes;
}

void SoA_LP_gen_root_primes(void);

tpPrime SoA_S(const tpPrime limit, void*, void*, void*)
{
	SoA_LP_gen_root_primes();

	sieve1 = sv1; sieve5 = sv5; sieve7 = sv7; sieve11 = sv11;
	AddPrime(2); AddPrime(3);

	tpPrime numprimes = 2;
	if(limit > SEGMENT_SIZE)
	for (segment_start = 0; segment_start < limit - SEGMENT_SIZE; segment_start += SEGMENT_SIZE)
	{
		numprimes += SoA_SieveSegment(segment_start, segment_start + SEGMENT_SIZE);
	}
	numprimes += SoA_SieveSegment(segment_start, limit);

	return numprimes;
}

constexpr tpPrime interval_base = 100'000'000'000'000ull; //1e14
//constexpr tpPrime interval_len = 1'000'000'000'000; //1e12
constexpr tpPrime interval_len = 500'000'000'00ull; //1e9
constexpr tpPrime interval_start = (interval_base - interval_len);
constexpr tpPrime interval_end = interval_base;

void SoA_Interval(void)
{
	cTimer tmr;
	tmr.Start();

	SoA_LP_gen_root_primes();
	tmr.LapTime(true, "root primes"); nln();

	sieve1 = sv1; sieve5 = sv5; sieve7 = sv7; sieve11 = sv11;
	AddPrime(2); AddPrime(3);

	tpPrime numprimes = 2;
	for (segment_start = interval_start; segment_start < interval_end - SEGMENT_SIZE; segment_start += SEGMENT_SIZE)
	{
		numprimes += SoA_SieveSegment(segment_start, segment_start + SEGMENT_SIZE);
	}
	numprimes += SoA_SieveSegment(segment_start, interval_end);

	nln(); tmr.LapTime(true, "interval"); tmr.Stop(); nln();
	std::cout << numprimes << " between " << interval_start << " and " << interval_end;
}