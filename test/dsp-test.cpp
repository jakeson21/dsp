#include <iir_filt.h>
#include <iostream>
#include <iomanip>      // std::setprecision
#include <vector>
#include <cmath>

#include <sigpack.h>

void do_filter(const std::vector<double>& num, const std::vector<double>& den)
{
    IIR_Filter<double> filt(num, den);
    size_t N = 101;
    std::vector<double> y;
    y.reserve(N);

    std::cout << filt.to_string() << std::endl;

    // Send step through
    y.push_back(filt(0.0));
    for (size_t n=1; n<N; n++)
    {
        auto yy = filt(1.0);
        y.push_back(yy);

        std::cout << yy << " ";
    }
    std::cout << std::endl;

    using namespace arma;
    using namespace sp;
    gplot gp;
    gp.window("Filter plots", 10, 10, 700, 700);
    arma::vec Mg(1024);
    arma::vec Ph(1024);
    arma::vec Y = arma::conv_to<arma::vec>::from(y);
    arma::vec b = arma::conv_to<arma::vec>::from(filt.num());
    arma::vec a = arma::conv_to<arma::vec>::from(filt.den());
    Mg = 20 * log10(freqz(b, a, 1024));
    Ph = phasez(b, a, 1024);

    // Plot
    gp.set_term("qt");
    gp.grid_on();
    gp.send2gp("set multiplot layout 3, 1");
    gp.plot_add(Y, "Step response");
    gp.plot_show();

    gp.send2gp("set xtics (\"0\" 1,\"0.5\" 512,\"1\" 1024)");
    gp.ylabel("Magnitude [dB]");
    gp.xlabel("Frequency [f/Fs]");
    gp.plot_add(Mg, "Filter");
    // gp.send2gp("set logscale xz");
    gp.plot_show();

    gp.send2gp("set xtics (\"0\" 1,\"0.5\" 512,\"1\" 1024)");
    gp.ylabel("Phase [rad]");
    gp.xlabel("Frequency [f/Fs]");
    gp.plot_add(Ph, "Filter");
    gp.plot_show();
    gp.send2gp("unset multiplot");
}



int main(int argc, char *argv[])
{
    double w0 = 0.25, bw = 0.025;
    if(argc==1)
    {
        printf("\nNo Command Line Argument Passed - Using defaults for w0=%f and Bw=%f.\n\n", w0, bw);
    }
    else if(argc==3)
    {
        w0 = atof(argv[1]);
        bw = atof(argv[2]);
        printf("\nCommand Line Argument Passed - Using w0=%f and Bw=%f.\n\n", w0, bw);
    }
    else
    {
        printf("\nInvalid number of Command Line Argument Passed.\n");
        return EXIT_FAILURE;
    }

    std::vector<double> b;
    std::vector<double> a;

    IIR_Bandpass(w0, bw, b, a);
    do_filter(b,a);
    IIR_Bandstop(w0, bw, b, a);
    do_filter(b,a);
    IIR_Lowpass(w0, b, a);
    do_filter(b,a);
    IIR_Highpass(w0, b, a);
    do_filter(b,a);
    
    return EXIT_SUCCESS;
}

/*
from scipy import signal
import matplotlib.pyplot as plt
import numpy as np

fs = 2.0  # Sample frequency (Hz)
f0 = 0.25  # Frequency to be retained (Hz)
Q = 10.0  # Quality factor
# Design peak filter
b, a = signal.iirpeak(f0, Q, fs)
s = signal.dlti(b, a, dt=1)
t, y = signal.dstep(s, n=101)
print(y)
plt.step(t, np.squeeze(y))
plt.grid()
plt.xlabel('n [samples]')
plt.ylabel('Amplitude')
plt.show()

w, h = signal.freqz(b, a, 1024)
fig, ax1 = plt.subplots()
ax1.set_title('Digital filter frequency response')
ax1.plot(w, 20 * np.log10(abs(h)), 'b')
ax1.set_ylabel('Amplitude [dB]', color='b')
ax1.set_xlabel('Frequency [rad/sample]')
ax2 = ax1.twinx()
angles = np.unwrap(np.angle(h))
ax2.plot(w, angles, 'g')
ax2.set_ylabel('Angle (radians)', color='g')
ax2.grid(True)
ax2.axis('tight')
plt.show()
*/