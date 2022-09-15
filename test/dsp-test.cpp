#include <iir_bandpass.h>
#include <iostream>
#include <iomanip>      // std::setprecision

int main()
{
    
    IIR_Bandpass<double> bp(0.25, 0.025);
    std::cout << bp.to_string() << std::endl;
    auto y = bp(0.0);
    std::cout << "0 : " << std::setprecision(16) << y << std::endl;
    for (size_t n=1; n<=100; n++)
    {
        y = bp(1.0);
        std::cout << n << " : " << std::setprecision(16) << y << std::endl;
    }

    return EXIT_SUCCESS;
}

/*
from scipy import signal
import matplotlib.pyplot as plt
import numpy as np

fs = 1.0  # Sample frequency (Hz)
f0 = 0.25  # Frequency to be retained (Hz)
Q = 10.0  # Quality factor
# Design peak filter
b, a = signal.iirpeak(f0, Q, fs)
s = signal.dlti(b, a, dt=1)
t, y = signal.dsetp(s, n=101)
print(y)
plt.step(t, np.squeeze(y))
plt.grid()
plt.xlabel('n [samples]')
plt.ylabel('Amplitude')
plt.show()
*/