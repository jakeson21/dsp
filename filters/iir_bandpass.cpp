#include <iir_bandpass.h>
#include <cmath>
#include <exception>
#include <complex>
#include <sstream>

template <typename T>
IIR_Bandpass<T>::IIR_Bandpass(double omega, double bandwidth)
{
    if (abs(omega)>1 || abs(bandwidth)>1)
    {
        throw std::invalid_argument("abs(omega)>1 or abs(bandwidth)>1");
    }

    omega *= M_PI;
    bandwidth *= M_PI;
    this->alpha = (1.0 - std::sin(bandwidth)) / std::cos(bandwidth);
    this->beta = std::cos(omega);
    this->gain = (1.0 - this->alpha) / 2.0;

    y[0] = T(0);
    y[1] = T(0);
    y[2] = T(0);
    x[0] = T(0);
    x[1] = T(0);
    x[2] = T(0);
}

template <typename T>
std::string IIR_Bandpass<T>::to_string()
{
    std::stringstream ss;
    ss << "alpha = " << alpha << ", ";
    ss << "beta = " << beta << ", ";
    ss << "gain = " << gain << ", ";
    return ss.str();
}

template <typename T>
T IIR_Bandpass<T>::operator()(T input)
{
    static const size_t n = 2;
    y[n-2] = y[n-1];
    y[n-1] = y[n];
    x[n-2] = x[n-1];
    x[n-1] = x[n];
    x[n] = input;

    // 1-a         1 - z^-2
    // ----*------------------------- = Hbp(z)
    //  2    1 - B(1+a)z^-1 + az^-2
    //
    // y[n] = G * (x[n] - x[n-2]) / ( B(1 + a)y[n-1] - ay[n-2] )
    y[n] = this->gain * (x[n] - x[n-2]) 
        + (this->beta * (T(1.0) + this->alpha) * y[n-1]) - this->alpha*y[n-2];
    return y[n];
}

template <>
std::complex<float> IIR_Bandpass<std::complex<float>>::operator()(std::complex<float> input)
{
    static const size_t n = 2;
    y[n-2] = y[n-1];
    y[n-1] = y[n];
    x[n-2] = x[n-1];
    x[n-1] = x[n];
    x[n] = input;

    // 1-a         1 - z^-2
    // ----*------------------------- = Hbp(z)
    //  2    1 - B(1+a)z^-1 + az^-2
    //
    // y[n] = G * (x[n] - x[n-2]) / ( B(1 + a)y[n-1] - ay[n-2] )
    y[n] = float(this->gain) * (x[n] - x[n-2]) 
        + (float(this->beta) * (1.0f + float(this->alpha)*y[n-1])) - float(this->alpha)*y[n-2];
    return y[n];
}

template class IIR_Bandpass<float>;
template class IIR_Bandpass<double>;
template class IIR_Bandpass<std::complex<float>>;
template class IIR_Bandpass<std::complex<double>>;
