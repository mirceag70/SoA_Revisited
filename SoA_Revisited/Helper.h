#pragma once

#include <iostream>
#include <iomanip>
#include <queue>
#include <cassert>
#include <array>
#include <thread>
#include <numeric>
#include <fstream>
#include <chrono>
#include <limits>
#include <functional>
#include <typeinfo>
#include <string_view>

class cTimer
{
private:

	enum class InternalStates { Initial, Iddle, Counting, Pause, Stop };

	InternalStates timer_state = InternalStates::Initial;

	long time_accumulated = 0;
	clock_t mark_start = 0;
	mutable clock_t mark_lap = 0;

	constexpr long TimeDiff2Miliseconds(clock_t const tmDiff) const
	{
		return tmDiff * (1'000 / CLOCKS_PER_SEC);
	}

	void OutputTime(long const tm, std::string_view const msg = "") const
	{
		std::cout << " ( ";
		if (not msg.empty()) std::cout << "[" << msg << "] ";
		std::cout << std::fixed << std::setprecision(0) << tm;
		std::cout << " ms ) ";
	}

public:

	cTimer(void) : timer_state(InternalStates::Iddle) {};

	void Start(void)
	{
		assert(timer_state == InternalStates::Iddle);
		timer_state = InternalStates::Counting;

		// mark down current clock
		mark_lap = mark_start = clock();	// use portable clock()
	}

	void Stop(bool const print_time = false, const char* msg = NULL)
	{
		switch (timer_state)
		{
		case InternalStates::Counting:
			time_accumulated += TimeDiff2Miliseconds(clock() - mark_start);
			mark_start = 0;
			[[fallthrough]];
		case InternalStates::Pause:
			timer_state = InternalStates::Stop;
			break;
		default:
			assert(false);
		}

		if (print_time)
			OutputTime(time_accumulated, msg ? msg : "StopTime");
	}

	long LapTime(const bool print_time = false, const char* msg = NULL) const
	{
		assert(timer_state == InternalStates::Counting);

		// get current clock
		clock_t tm = clock();
		long laptime = TimeDiff2Miliseconds(tm - mark_lap);
		//prepare for next lap
		mark_lap = tm;

		if (print_time)
			OutputTime(laptime, msg ? msg : "LapTime");

		return laptime;
	}

	long GetTime(bool print_time = false) const
	{
		clock_t tm = 0;

		switch (timer_state)
		{
		case InternalStates::Counting:
			tm = TimeDiff2Miliseconds(tm - mark_start);
			[[fallthrough]];
		case InternalStates::Pause:
		case InternalStates::Stop:
			tm += time_accumulated;
			break;
		default:
			assert(false);
		}

		if (print_time)
			OutputTime(tm, "Time");

		return tm;
	}

	long Pause(bool print_time = false)
	{
		assert(timer_state == InternalStates::Counting);
		timer_state = InternalStates::Pause;

		time_accumulated += TimeDiff2Miliseconds(clock() - mark_start);
		mark_start = 0;

		if (print_time)
			OutputTime(time_accumulated, "PauseTime");

		return time_accumulated;
	}

	void Resume(void)
	{
		assert(timer_state == InternalStates::Pause);
		timer_state = InternalStates::Counting;

		mark_lap = mark_start = clock();
	}

	void Reset(void)
	{
		assert(timer_state == InternalStates::Stop);
		timer_state = InternalStates::Iddle;

		time_accumulated = mark_lap = mark_start = 0;
	}
};

constexpr std::size_t default_iterations = 5;

inline const void nln(bool delimiter = false)
{
	std::cout << std::endl;
	if (delimiter)
		std::cout << "\t--------------------------" << std::endl;
}

template<std::size_t range_length, typename INTEGRAL = int>
consteval auto Range(const int offset = 0)
{
	static_assert(std::is_integral<INTEGRAL>::value,
		"Integral type required for Range.");

	std::array<INTEGRAL, range_length> arr{};
	for (int i = 0; i < range_length; i++)
		arr[i] = i + offset;
	return arr;
}

template<class NUMBER>
constexpr double Average(std::vector<NUMBER> vec)
{
	static_assert(std::is_arithmetic<NUMBER>::value,
		"Numeric type required for Average.");

	return std::accumulate(vec.begin(), vec.end(), 0.0) / vec.size();
}

inline uint64_t PI_Nmax(uint64_t Nmax)
{
	double logNmax = log(Nmax);
	return (uint64_t)(Nmax / logNmax * (1 + 1.2762 / logNmax));
}

// scalar data type used for prime numbers
//go with 32b for small experiments
#define USE_64_BITS_PRIMES
//we need 64b for larger numbers		
#ifndef USE_64_BITS_PRIMES
typedef		uint32_t		tpPrime;
#else
typedef		uint64_t		tpPrime;
#endif // !USE_64_BITS

class cChecker
{
	std::ifstream file_with_primes;
	tpPrime max_value = 0;

public:
	tpPrime check_next_prime(tpPrime nGenerated)
	{
		if (nGenerated == 2)
			BackToZero();

		tpPrime nRead;
		file_with_primes >> nRead;
		if (nRead <= max_value)
			return nRead == nGenerated;
		else
			return true;
	}

	void BackToZero(void)
	{
		file_with_primes.seekg(0, file_with_primes.beg);
	}

	cChecker(std::string file_name = "C:/50MilPrimes.txt", tpPrime maxPrime = 982'451'653) :
		max_value(maxPrime)
	{
		file_with_primes.open(file_name, std::ifstream::in);
		assert(file_with_primes.is_open());
	}

	~cChecker() { file_with_primes.close(); }
};

inline uint8_t idx_last_primes = 0;
inline std::array<uint64_t, 256> last_primes;
inline void AddPrime(uint64_t prime)
{ 
#ifdef _DEBUG
	static cChecker ckr;
	assert(ckr.check_next_prime(prime));
#endif // _DEBUG

	last_primes[idx_last_primes++] = prime; 
}

inline thread_local uint8_t idx_last_primes_t = 0;
inline thread_local std::array<uint64_t, 256> last_primes_t;
inline void AddPrimeThreaded(uint64_t prime)
{
	last_primes_t[idx_last_primes_t++] = prime;
}

#ifndef ITERATIONS
#define ITERATIONS	5
#endif 

template<typename INTEGRAL1, typename INTEGRAL2, typename INTEGRAL3,
		std::size_t size1 = 0, std::size_t size2 = 0, std::size_t size3 = 0>
void Try_Sieve(const uint64_t LIMIT, std::string message,
	std::function<uint64_t(uint64_t, INTEGRAL1[], INTEGRAL2[], INTEGRAL3[])> Sieve, bool show_all = true)
{
	static_assert(std::is_integral<INTEGRAL1>::value, "Integral type required for INTEGRAL1.");
	static_assert(std::is_integral<INTEGRAL2>::value, "Integral type required for INTEGRAL2.");
	static_assert(std::is_integral<INTEGRAL3>::value, "Integral type required for INTEGRAL3.");

	INTEGRAL1* v1 = NULL; INTEGRAL2* v2 = NULL; INTEGRAL3* v3 = NULL;

	if (size1 > 0) v1 = new INTEGRAL1[size1];
	if (size2 > 0) v2 = new INTEGRAL2[size2];
	if (size3 > 0) v3 = new INTEGRAL3[size3];

	nln(true);
	cTimer tmr;
	std::vector<decltype(tmr.LapTime())> times;

	tmr.Start();
	std::cout << " - " << message << " - ";
	for (auto i : Range<ITERATIONS>())
	{

		auto numPrimes = Sieve(LIMIT, v1, v2, v3);

		if (i == 0 or show_all)
		{
			nln();
			std::cout << numPrimes << " primes up to " << LIMIT;
			times.push_back(tmr.LapTime(true));
		}
		else
		{
			times.push_back(tmr.LapTime(false));
		}
	}
	nln(true);
	tmr.Stop();
	std::cout << "Average compute time: " << Average(times);

	if (v1) delete[] v1;
	if (v2) delete[] v2;
	if (v3) delete[] v3;
}

constexpr uint8_t BIT_MASK[]	   = {  1,  2,  4,  8,  16,  32,  64,  128 };
constexpr uint8_t BIT_RESET_MASK[] = { ~1, ~2, ~4, ~8, ~16, ~32, ~64, ~128 };

inline uint64_t idx2no(uint64_t idx) { return 3 * idx + 5 - (idx & 1); };
inline uint64_t no2idx(uint64_t no)  { return no / 3 - 1; };

inline bool GetBit(tpPrime n, uint8_t sieve[])
{
	tpPrime q = (n / 12) / 8, r = (n / 12) % 8;
	return (sieve[q] & BIT_MASK[r]);
};
inline void ResetBit(tpPrime n, uint8_t sieve[])
{
	tpPrime q = (n / 12) / 8, r = (n / 12) % 8;
	sieve[q] &= BIT_RESET_MASK[r];
};
inline void FlipBit(tpPrime n, uint8_t sieve[])
{
	tpPrime q = (n / 12) / 8, r = (n / 12) % 8;
	sieve[q] ^= BIT_MASK[r];
};

inline uint8_t* sieve1; inline uint8_t* sieve5; inline uint8_t* sieve7; inline uint8_t* sieve11;

inline void Pattern11(const tpPrime start, const tpPrime stop, uint8_t sieve[])
{
	const tpPrime xmax = (tpPrime)floor(sqrt(stop / 4));
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
		//sieve
		for (; ; y += 6)
		{
			tpPrime n = n0 + (y * y);
			if (n > stop) break;
			FlipBit(n, sieve);
		}
	}
};
inline void Pattern12(const tpPrime start, const tpPrime stop, uint8_t sieve[])
{
	//for x = 0 we can not start with y = 1
	uint8_t step = 2;
	const tpPrime yy = (tpPrime)ceil(sqrt(start));
	tpPrime y = 6 * (yy / 6) + 5;
	if (y < yy)
	{
		y += 2;
		if (y < yy) y += 4;
		else step = 4;
	}
	for (; ; y += step, step = 6 - step)
	{
		tpPrime n = y * y;
		if (n > stop) break;
		FlipBit(n, sieve);
	}
	const tpPrime xmax = (tpPrime)floor(sqrt(stop / 4));
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
		//sieve
		for (; ; y += step, step = 6 - step)
		{
			tpPrime n = (4 * x * x) + (y * y);
			if (n > stop) break;
			FlipBit(n, sieve);
		}
	}
};
inline void Pattern5(const tpPrime start, const tpPrime stop, uint8_t sieve[])
{
	const tpPrime xmax = (tpPrime)floor(sqrt(stop / 4));
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
		//sieve
		for (; ; y += step, step = 6 - step)
		{
			tpPrime n = n0 + (y * y);
			if (n > stop) break;
			FlipBit(n, sieve);
		}
	}
};
inline void Pattern7(const tpPrime start, const tpPrime stop, uint8_t sieve[])
{
	const tpPrime xmax = (tpPrime)floor(sqrt(stop / 3));
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
		//sieve			
		for (; ; y += step, step = 6 - step)
		{
			tpPrime n = n0 + (y * y);
			if (n > stop) break;
			FlipBit(n, sieve);
		}
	}
};
inline void Pattern111(const tpPrime start, const tpPrime stop, uint8_t sieve[])
{
	const tpPrime xmax = (tpPrime)floor(sqrt(stop / 2));
	tpPrime x = (tpPrime)ceil(sqrt(start / 3));
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
		//sieve			
		for (; x > y; y += step, step = 6 - step)
		{
			tpPrime n = n0 - (y * y);
			if (n < start) break;
			FlipBit(n, sieve);
		}
	}
};
inline void Pattern112(const tpPrime start, const tpPrime stop, uint8_t sieve[])
{
	const tpPrime xmax = (tpPrime)floor(sqrt(stop / 2));
	tpPrime x = (tpPrime)ceil(sqrt(start / 3));
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
		//sieve			
		for (; x > y; y += step, step = 6 - step)
		{
			tpPrime n = n0 - (y * y);
			if (n < start) break;
			FlipBit(n, sieve);
		}
	}
};

constexpr unsigned ROOT_LIMIT = 10'000'000u;
constexpr unsigned NO_ROOT_PRIMES = 664'577u;

constexpr unsigned STEP_SIZE = 65'000'000;

template<typename FN1, typename FN2>
void Sieve2Patterns(FN1 fn1, FN2 fn2, const tpPrime start, const tpPrime stop, uint8_t sieve[])
{
	tpPrime iter_min = start;
	if ((stop - start) > STEP_SIZE)
		for (; iter_min < (stop - STEP_SIZE); iter_min += STEP_SIZE)
		{
			fn1(iter_min, iter_min + STEP_SIZE, sieve);
			fn2(iter_min, iter_min + STEP_SIZE, sieve);
		}
	fn1(iter_min, stop, sieve); fn2(iter_min, stop, sieve);
};

template<typename FN>
void Sieve1Pattern(FN fn1, const tpPrime start, const tpPrime stop, uint8_t sieve[])
{
	tpPrime iter_min = start;
	if ((stop - start) > STEP_SIZE)
		for (; iter_min < (stop - STEP_SIZE); iter_min += STEP_SIZE)
		{
			fn1(iter_min, iter_min + STEP_SIZE, sieve);
		}
	fn1(iter_min, stop, sieve);
};
