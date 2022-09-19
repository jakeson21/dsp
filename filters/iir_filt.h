#pragma once
#include <complex>
#include <string>
#include <vector>


template <typename T>
class IIR_Filter
{
public:
    IIR_Filter();
    IIR_Filter(const std::vector<double>& b, const std::vector<double>& a);
    virtual ~IIR_Filter(){}

    T operator()(T x);

    virtual std::string to_string();
    std::vector<double> num() const { return b; }
    std::vector<double> den() const { return a; }

protected:
    size_t N = 0; // filter order
    std::vector<T> y; // {z-N, z-N-1, ..., z-2, z-1, 1}
    std::vector<T> x; // {z-N, z-N-1, ..., z-2, z-1, 1}
    std::vector<double> b; // numerator
    std::vector<double> a; // denominator
};

void IIR_Bandpass(double w0, double bw, std::vector<double>& b, std::vector<double>& a);
void IIR_Bandstop(double w0, double bw, std::vector<double>& b, std::vector<double>& a);
void IIR_Lowpass(double w0, std::vector<double>& b, std::vector<double>& a);
void IIR_Highpass(double w0, std::vector<double>& b, std::vector<double>& a);
