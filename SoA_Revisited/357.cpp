#include "Helper.h"

constexpr auto szM = 500'000;
bool vM[szM];

uint64_t Sieve357(const uint64_t Nmax, uint64_t IP[], uint64_t IQ[], int32_t JQ[], bool initialize = false)
{
    uint64_t numPrimes = 0;     //for this call
    for (int i = 0; i < szM; i++) vM[i] = true;

    static uint64_t nstart, nend;
    static uint32_t nq, stepq;
    static int32_t iq;      //idx of the last root prime generated so far
    if (initialize)
    {
        vM[1] = false;
        JQ[0] = 3;  IQ[0] = 9;
        iq = 0; stepq = 2; nq = 5;
        nstart = 0;
        nend = szM;
        numPrimes = 1; AddPrime(2); // account for 2
    }

    uint64_t nsqrt = (uint64_t)ceil(sqrt(nstart + szM));
    for (; nq <= nsqrt; nq += stepq, stepq = 6 - stepq)
    {   //generate rest of root primes for this iteration
        bool isprime = true;
        for (int i = 2; i <= iq; i++)
        {
            if ((nq / JQ[i]) * JQ[i] == nq)
            {
                isprime = false;
                break;
            }
        }
        if (isprime)
        {
            if (nq > 65500)
                iq = iq;
            iq++;
            JQ[iq] = nq;
            IQ[iq] = ((uint64_t)nq) * nq;
        }
    }

    //strike out all composites
    for (int i = 0; i <= iq; i++)
    {
        uint64_t n = IQ[i];
        uint32_t stp = 2 * JQ[i];
        while (n < nend)
        {
            vM[n - nstart] = false;
            n += stp;
        }
        IQ[i] = n;
    }
    //move primes in IP
    for (int i = 1; i < szM; i += 2)
        if (vM[i])
            if ((i + nstart) > Nmax)
                break;
            else
            {
                AddPrime(IP[numPrimes++] = nstart + i);
            }

    nstart = nend;
    nend += szM;

    return numPrimes;
}

void Try357(const uint64_t LIMIT)
{
    nln();
    uint64_t sz = PI_Nmax(szM);
    uint64_t* vulPrimes = new uint64_t[sz + 1];
    std::cout << sz << " : ";
    sz = PI_Nmax((uint64_t)ceil(sqrt(LIMIT)));
    std::cout << sz;
    uint64_t* IQ = new uint64_t[sz + 1];
    int32_t* JQ = new int32_t[sz + 1];

    nln(true);
    cTimer tmr;
    std::vector<decltype(tmr.LapTime())> times;

    tmr.Start();
    std::cout << " - Sieve 357 - ";
    for (auto i : Range<1>())
    {
        nln();

        uint64_t numPrimes = Sieve357(LIMIT, vulPrimes, IQ, JQ, true);
        for (uint64_t start = szM; start < LIMIT; start += szM)
            numPrimes += Sieve357(LIMIT, vulPrimes, IQ, JQ);

        std::cout << numPrimes << " primes up to " << LIMIT;
        times.push_back(tmr.LapTime(true));
    }
    nln(true);
    tmr.Stop();
    std::cout << "Average compute time: " << Average(times);

    nln(true);
    times.clear();

    delete[] vulPrimes;
    delete[] IQ;
    delete[] JQ;
}
