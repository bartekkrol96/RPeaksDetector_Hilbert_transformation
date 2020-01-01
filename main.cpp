#include "libalglib/stdafx.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "libalglib/fasttransforms.h"
#include <complex>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <wfdb/wfdb.h>

using namespace alglib;
using namespace std;

vector<double> hilbert(vector<double> input)
{
    double arr[input.size()];
    copy(input.begin(), input.end(), arr);
    real_1d_array x;
    x.setcontent(input.size(), arr);
    complex_1d_array f;
    fftr1d(x, f);
//    cout << "Fourier" << endl;
//    printf("%s\n", f.tostring(3).c_str());

    vector<std::complex<double>> fourier;
    for (int i=0; i<f.length(); i++)
    {
        std::complex<double> mycomplex(f[i].x, f[i].y);
        fourier.push_back(mycomplex);
    }

    int N = fourier.size();
    for(int i=(int)N/2+1;i<=N; i++)
    {
        fourier[i] = 0;
    }
    for(int i=1;i<(int)N/2; i++)
    {
        fourier[i] = std::complex<double> ((2*fourier[i].real()), (2*fourier[i].imag()));
    }

    alglib::complex array2[N];
    for(int i=0; i<N; i++)
    {
        alglib::complex a;
        a.x = fourier[i].real();
        a.y = fourier[i].imag();
        array2[i] = a;
    }
    complex_1d_array x2;
    x2.setcontent(N, array2);
    fftc1dinv(x2);
//    cout << "Fourier odwrotny" << endl;
//    printf("%s\n", x2.tostring(3).c_str());

    vector<double> hilbert;
    for (int i=0; i<x2.length(); i++)
    {
        std::complex<double> mycomplex(x2[i].x, x2[i].y);
        hilbert.push_back(mycomplex.imag());
//        cout << hilbert[i] << endl;
    }
    return hilbert;
}

vector<double> gradient(vector<double> input)
{
    // TODO napisac funckje liczaca pochodna
    vector<double> gradient;
    double dx = 0.001;
    double grad[input.size()];
    grad[0] = (input[1] - input[0]) / dx;
    for(int i=1; i<input.size()-1; i++)
    {
        grad[i] = (input[i+1] - input[i-1]) / (2*dx);  // for i in [1,N-2]
    }
    grad[input.size()] = (input[input.size()-1] - input[input.size()-2]) / dx;

    for(int j=0; j<input.size(); j++)
    {
        gradient.push_back(grad[j]);
    }
  return gradient;
}

double set_threshold(vector<double> segment, double prev_max)
{
    double max_of_segment  = *max_element(segment.begin(), segment.end());
    transform(segment.begin(), segment.end(), segment.begin(), [](double x){return x*x;});
    double average = accumulate(segment.begin(), segment.end(), 0.0)/segment.size();
    double rms = sqrt(average);
    double ratio = rms/max_of_segment;
    double threshold;
    cout << "Srednia : " <<average <<endl;
    cout << "RMS : " <<rms <<endl;
    cout << "Ratio : "<<ratio <<endl;
    cout << "Max : "<<max_of_segment<<endl;
    if (ratio > 0.18)
    {
        threshold = 0.39*max_of_segment;
    }
    else
    {
        threshold = 1.6*rms;
    }
//    if ((max_of_segment)>2*prev_max)
//    {
//        threshold = 0.39*prev_max;
//    }
    return threshold;
}

int main() {

    vector<double> segment{-0.145,-0.145,-0.145,
                           -0.145,-0.145,-0.145,-0.145,-0.145,-0.12,-0.135,-0.145,-0.15,-0.16,-0.155,-0.16,-0.175,
                           -0.18,-0.185,-0.17,-0.155,-0.175,-0.18,-0.19,-0.18,-0.155,-0.135,-0.155,-0.19,-0.205,
                           -0.235,-0.225,-0.245,-0.25,-0.26,-0.275,-0.275,-0.275,-0.265,-0.255,-0.265,-0.275,-0.29,
                           -0.29,-0.29,-0.29,-0.285,-0.295,-0.305,-0.285,-0.275,-0.275,-0.28,-0.285,-0.305,-0.29,
                           -0.3,-0.28,-0.29,-0.3,-0.315,-0.32,-0.335,-0.36,-0.385,-0.385,-0.405,-0.455,-0.485,-0.485,
                           -0.425,-0.33,-0.22,-0.07,0.12,0.375,0.62,0.78,0.84,0.765,0.52,0.17,-0.165,-0.365,-0.435,
                           -0.425,-0.37,-0.33,-0.325,-0.335,-0.345,-0.33,-0.325,-0.315,-0.31,-0.32,-0.335,-0.34,-0.325,
                           -0.345,-0.335,-0.33,-0.335,-0.33,-0.325,-0.33,-0.33,-0.345,-0.355,-0.335,-0.325,-0.305,
                           -0.32,-0.32,-0.33,-0.34,-0.335,-0.34,-0.345,-0.355,-0.355,-0.34,-0.33,-0.33,-0.33,-0.34,
                           -0.35,-0.325,-0.325,-0.33,-0.33,-0.335,-0.335,-0.34,-0.33,-0.34,-0.35,-0.355,-0.35,-0.345,
                           -0.33,-0.32,-0.335,-0.33,-0.345,-0.33,-0.335,-0.335,-0.345,-0.345,-0.355,-0.34,-0.34,-0.335,-0.33,
                           -0.35,-0.35,-0.345,-0.335,-0.335,-0.335,-0.35,-0.355,-0.355,-0.345,-0.345,-0.335,-0.35,
                           -0.36,-0.36,-0.36,-0.365,-0.36,-0.37,-0.385,-0.37,-0.36,-0.355,-0.36,-0.375,-0.375,-0.365,-0.365,
                           -0.36,-0.36,-0.365,-0.37,-0.355,-0.33,-0.325,-0.325,-0.335,-0.34,-0.315,-0.3,-0.3,-0.29,
                           -0.295,-0.29,-0.285,-0.275,-0.255,-0.25,-0.25,-0.265,-0.255,-0.245,-0.23,-0.245,-0.245,
                           -0.255,-0.255,-0.24,-0.25,-0.255,-0.245,-0.255,-0.25,-0.25,-0.265,-0.26,-0.26,-0.265,-0.27,-0.265,-0.26,-0.275,-0.28,
                           -0.29,-0.275,-0.27,-0.26,-0.28,-0.28,-0.285,-0.275,
                           -0.275,-0.265,-0.27,-0.285,-0.29,-0.28,-0.275,-0.285,-0.28,-0.3,-0.3,-0.305,-0.295,-0.3,-0.31};

    vector<double> test_gradient{1.0, 2.0, 4.0, 7.0, 11.0, 16.0};
    vector<double> moj_test;
    moj_test = gradient(test_gradient);
    cout << "Test gradeint : " <<endl;
    // TODO expexted: array([1. , 1.5, 2.5, 3.5, 4.5, 5. ])
    for (int i=0; i<moj_test.size(); i++)
    {
        cout<< (moj_test[i]) << endl;
    }

    vector<double> derivative;
    derivative = gradient(segment);
    vector<double> output;
    output = hilbert(derivative);

    vector<int> r_peaks;
    vector<double> above_threshold;
    double prev_max=0;
    double threshold = set_threshold(output, prev_max);

    cout << "Prog : " <<threshold << endl;
    copy_if(output.begin(), output.end(),
            back_inserter(above_threshold),[&threshold](double n){return n>threshold;});


    cout << "powyzej progu : " << above_threshold.size() <<endl;
    cout << "Powzyej progu : "<<endl;
    for (int i=0; i<above_threshold.size(); i++)
    {
        cout << above_threshold[i] << endl;
    }

    vector<int> above_threshold_indexes;
    for (int j=0; j<above_threshold.size(); j++)
    {
        cout << "taka wartosc szukamy : " << above_threshold[j] <<endl;
        for (int i=0; i<output.size(); ++i)
        {
            if (output[i]==above_threshold[j])
            {
                cout << "Index : " << i << endl;
                above_threshold_indexes.push_back(i);
            }
        }
    }


    double window = (int)(360*0.12);
    vector<int> r_peaks_in_segment;

    // todo nalozyc okna na te powyzej progu
    for(int i=0; i<above_threshold.size(); i++)
    {
        return 0;
    }

    // todo segmentacja i globalne: prev_max i r_peaks in signal





    return 0;
}