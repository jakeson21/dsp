#pragma once
#include <complex>
#include <string>

template <typename T>
class IIR_Bandpass
{
public:
    IIR_Bandpass() = delete;
    IIR_Bandpass(double w0, double bw);
    virtual ~IIR_Bandpass(){}

    T operator()(T x);

    std::string to_string();

protected:
    T y[3]; // {n, n-1, n-2}
    T x[3]; // {n, n-1, n-2}

    double w0;
    double bw;
    double alpha;
    double beta;
    double gain;
};
