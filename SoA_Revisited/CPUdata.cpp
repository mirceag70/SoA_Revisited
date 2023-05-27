#define _WIN32_DCOM
#include <iostream>
#include <comdef.h>
#include <Wbemidl.h>

#pragma comment(lib, "wbemuuid.lib")

//This code uses the WMI API to connect to the root\WMI namespace and execute a query to retrieve the CurrentTemperature property from the MSAcpi_ThermalZoneTemperature class.
HRESULT GetCpuTemperature(LPLONG pTemperature)
{
    if (pTemperature == NULL)
        return E_INVALIDARG;

    *pTemperature = -1;
    HRESULT ci = CoInitialize(NULL);
    HRESULT hr = CoInitializeSecurity(NULL, -1, NULL, NULL, RPC_C_AUTHN_LEVEL_DEFAULT, RPC_C_IMP_LEVEL_IMPERSONATE, NULL, EOAC_NONE, NULL);
    if (SUCCEEDED(hr))
    {
        IWbemLocator* pLocator;
        hr = CoCreateInstance(CLSID_WbemAdministrativeLocator, NULL, CLSCTX_INPROC_SERVER, IID_IWbemLocator, (LPVOID*)&pLocator);
        if (SUCCEEDED(hr))
        {
            IWbemServices* pServices;
            BSTR ns = SysAllocString(L"root\\WMI");
            hr = pLocator->ConnectServer(ns, NULL, NULL, NULL, 0, NULL, NULL, &pServices);
            pLocator->Release();
            SysFreeString(ns);
            if (SUCCEEDED(hr))
            {
                BSTR query = SysAllocString(L"SELECT * FROM MSAcpi_ThermalZoneTemperature");
                BSTR wql = SysAllocString(L"WQL");
                IEnumWbemClassObject* pEnum;
                hr = pServices->ExecQuery(wql, query, WBEM_FLAG_RETURN_IMMEDIATELY | WBEM_FLAG_FORWARD_ONLY, NULL, &pEnum);
                SysFreeString(wql);
                SysFreeString(query);
                pServices->Release();
                if (SUCCEEDED(hr))
                {
                    IWbemClassObject* pObject;
                    ULONG returned;
                    hr = pEnum->Next(WBEM_INFINITE, 1, &pObject, &returned);
                    pEnum->Release();
                    if (SUCCEEDED(hr))
                    {
                        BSTR temp = SysAllocString(L"CurrentTemperature");
                        VARIANT v;
                        VariantInit(&v);
                        hr = pObject->Get(temp, 0, &v, NULL, NULL);
                        pObject->Release();
                        SysFreeString(temp);
                        if (SUCCEEDED(hr))
                        {
                            *pTemperature = V_I4(&v);
                        }
                        VariantClear(&v);
                    }
                }
            }
        }
    }
    if (ci == S_OK)
    {
        CoUninitialize();
    }
    return hr;
}
//The temperature value is returned in tenths of degrees Kelvin and needs to be converted to degrees Celsius1.
double GetCPUTemp(void)
{
    LONG temp;
    GetCpuTemperature(&temp);
    double tempCelsius = temp;
    tempCelsius = temp / 10 - 273.15;
    return tempCelsius;
}


#include <windows.h>
#include <vector>
#include <winternl.h>

#pragma comment(lib, "Ntdll.lib")

typedef struct _SYSTEM_PROCESSOR_PERFORMANCE_INFORMATION_R
{
    LARGE_INTEGER IdleTime;
    LARGE_INTEGER KernelTime;
    LARGE_INTEGER UserTime;
    LARGE_INTEGER DpcTime;
    LARGE_INTEGER InterruptTime;
    ULONG InterruptCount;
} SYSTEM_PROCESSOR_PERFORMANCE_INFORMATION_R;

static long long toInteger(LARGE_INTEGER const& integer)
{
#ifdef INT64_MAX // Does the compiler natively support 64-bit integers?
    return integer.QuadPart;
#else
    return (static_cast<long long>(integer.HighPart) << 32) | integer.LowPart;
#endif
}

//This code uses the NtQuerySystemInformation function to retrieve the SYSTEM_PROCESSOR_PERFORMANCE_INFORMATION structure for each core of the CPU.
//The structure contains information about the idle time, kernel time, and user time for each core.
//The code calculates the CPU usage per core by subtracting the idle time from the kernel time and user time and dividing by the total time1
class CPU
{
public:
    uint64_t prev_idle = 0;
    uint64_t prev_ker = 0;
    uint64_t prev_user = 0;
    uint64_t cur_idle = 0;
    uint64_t cur_ker = 0;
    uint64_t cur_user = 0;

    double get()
    {
        SYSTEM_PROCESSOR_PERFORMANCE_INFORMATION_R* a = new SYSTEM_PROCESSOR_PERFORMANCE_INFORMATION_R[4]; // 4 is the total of CPU (4 cores)
        NtQuerySystemInformation(SystemProcessorPerformanceInformation, a, sizeof(SYSTEM_PROCESSOR_PERFORMANCE_INFORMATION_R) * 4, NULL);

        prev_idle = cur_idle;
        prev_ker = cur_ker;
        prev_user = cur_user;

        cur_idle = 0;
        cur_ker = 0;
        cur_user = 0;

        // 4 is the total of CPU (4 cores)
        // Sum up the SYSTEM_PROCESSOR_PERFORMANCE_INFORMATION_R array so I can get the utilization from all of the CPU
        for (int i = 0; i < 4; ++i)
        {
            SYSTEM_PROCESSOR_PERFORMANCE_INFORMATION_R b = a[i];
            cur_idle += toInteger(b.IdleTime);
            cur_ker += toInteger(b.KernelTime);
            cur_user += toInteger(b.UserTime);
        }

        std::cout << "Cur idle " << cur_idle << '\n';
        std::cout << "Cur ker " << cur_ker << '\n';
        std::cout << "Cur user " << cur_user << '\n';

        uint64_t delta_idle = cur_idle - prev_idle;
        uint64_t delta_kernel = cur_ker - prev_ker;
        uint64_t delta_user = cur_user - prev_user;
        std::cout << "Delta idle " << delta_idle << '\n';
        std::cout << "Delta ker " << delta_kernel << '\n';
        std::cout << "Delta user " << delta_user << '\n';

        uint64_t total_sys = delta_kernel + delta_user;
        uint64_t kernel_total = delta_kernel - delta_idle;

        delete[] a;
        // return (total_sys - delta_idle) * 100.0 / total_sys;
        return (kernel_total + delta_user) * 100.0 / total_sys;
    }
};

void CPUtilizatoin()
{
    CPU a;
    std::cout << "starting" << '\n';
    while (1)
    {
        std::cout << a.get() << '\n';
        Sleep(1000);
    }
}

