#include <iir_filt.h>
#include <cmath>
#include <exception>
#include <complex>
#include <sstream>
#include <string.h>
#include <cassert>

template <typename T>
IIR_Filter<T>::IIR_Filter(const std::vector<double>& ib, const std::vector<double>& ia)
{
    assert(ib.size() == ia.size());
    this->b.assign(ib.begin(), ib.end());
    this->a.assign(ia.begin(), ia.end());

    N = this->b.size()-1;
    y.resize(N+1);
    x.resize(N+1);
}



template <typename T>
T IIR_Filter<T>::operator()(T input)
{
    static const size_t n = N;
    memcpy(&y[0], &y[1], y.size()*sizeof(T));
    memcpy(&x[0], &x[1], x.size()*sizeof(T));

    x[n] = input;
    y[n] = T(this->b[0])*x[n];
    for (size_t k=1; k<=N; k++)
    {
        y[n] += T(this->b[k])*x[n-k] - T(this->a[k])*y[n-k];
    }
    return y[n];
}

template <>
std::complex<float> IIR_Filter<std::complex<float>>::operator()(std::complex<float> input)
{
    static const size_t n = N;
    memcpy(&y[0], &y[1], n*sizeof(std::complex<float>));
    memcpy(&x[0], &x[1], n*sizeof(std::complex<float>));

    x[n] = input;
    y[n] = float(this->b[0])*x[n];
    for (size_t k=1; k<=N; k++)
    {
        y[n] += float(this->b[k])*x[n-k] - float(this->a[k])*y[n-k];
    }
    return y[n];
}

template <typename T>
std::string IIR_Filter<T>::to_string()
{
    std::stringstream ss;
    ss << "\nN = " << N << std::endl;
    size_t k=0;
    ss << "\nb = {";
    for (; k<N; k++)
    {
        ss << this->b[k] << ", ";
    } 
    ss << this->b[k] << "}";

    k=0;
    ss << "\na = {";
    for (; k<N; k++)
    {
        ss << this->a[k] << ", ";
    } 
    ss << this->a[k] << "}";

    return ss.str();
}

void IIR_Bandpass(double omega, double bandwidth, std::vector<double>& b, std::vector<double>& a)
{
    if (abs(omega)>1 || abs(bandwidth)>1)
    {
        throw std::invalid_argument("abs(omega)>1 or abs(bandwidth)>1");
    }

    omega *= M_PI;
    bandwidth *= M_PI;
    double alpha = (1.0 - std::sin(bandwidth)) / std::cos(bandwidth);
    double beta = std::cos(omega);
    double gain = (1.0 - alpha) / 2.0;

    b = {gain, 0.0, -gain};
    a = {1.0, -beta*(1.0 + alpha), alpha};
}

void IIR_Bandstop(double omega, double bandwidth, std::vector<double>& b, std::vector<double>& a)
{
    if (abs(omega)>1 || abs(bandwidth)>1)
    {
        throw std::invalid_argument("abs(omega)>1 or abs(bandwidth)>1");
    }

    omega *= M_PI;
    bandwidth *= M_PI;
    double alpha = (1.0 - std::sin(bandwidth)) / std::cos(bandwidth);
    double beta = std::cos(omega);
    double gain = (1.0 + alpha) / 2.0;

    b = {gain, -gain*2*beta, gain};
    a = {1.0, -beta*(1.0 + alpha), alpha};
}

void IIR_Lowpass(double omega, std::vector<double>& b, std::vector<double>& a)
{
    if (abs(omega)>1)
    {
        throw std::invalid_argument("abs(omega)>1");
    }

    omega *= M_PI;
    double alpha = (1.0 - std::sin(omega)) / std::cos(omega);
    double gain = (1.0 - alpha) / 2.0;

    b = {gain, gain, 0.0};
    a = {1.0, -alpha, 0.0};
}

void IIR_Highpass(double omega, std::vector<double>& b, std::vector<double>& a)
{
    if (abs(omega)>1)
    {
        throw std::invalid_argument("abs(omega)>1");
    }

    omega *= M_PI;
    double alpha = (1.0 - std::sin(omega)) / std::cos(omega);
    double gain = (1.0 + alpha) / 2.0;

    b = {gain, -gain, 0.0};
    a = {1.0, -alpha, 0.0};
}

template class IIR_Filter<float>;
template class IIR_Filter<double>;
template class IIR_Filter<std::complex<float>>;
template class IIR_Filter<std::complex<double>>;