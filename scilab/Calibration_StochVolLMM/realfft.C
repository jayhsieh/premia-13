#include "realfft.h"

// "realfft.C" use the definition of direct fourier transform with positive in exponential and without any term outside the integral, then the definition of the inverse fourier tranform with negative in exponential and with 1/(2.PI)

void realfastfouriertransform(double *a, int tnn, int inversefft);//for fft put "fault" in inversefft, for inverse fft put "true" in inversefft



void realfastfouriertransform(double *a, int tnn, int inversefft)
{
    double twr;
    double twi;
    double twpr;
    double twM_PI;
    double twtemp;
    double ttheta;
    int i;
    int i1;
    int i2;
    int i3;
    int i4;
    double c1;
    double c2;
    double h1r;
    double h1i;
    double h2r;
    double h2i;
    double wrs;
    double wis;
               
    int nn;
    int ii;
    int jj;
    int n;
    int mmax;
    int m;
    int j;
    int istep;
    int isign;
    double wtemp;
    double wr;
    double wpr;
    double wM_PI;
    double wi;
    double theta;
    double tempr;
    double temM_PI;
 double ibidon;

    if( tnn==1 )
    {
        return;
    }
    if( !inversefft )
    {
        ttheta = 2*M_PI/tnn;
        c1 = 0.5;
        c2 = -0.5;
    }
    else
    {
        ttheta = 2*M_PI/tnn;
        c1 = 0.5;
        c2 = 0.5;
        ttheta = -ttheta;
  ibidon = sin(0.5*ttheta);
        twpr = -2.0*ibidon*ibidon;
        twM_PI = sin(ttheta);
        twr = 1.0+twpr;
        twi = twM_PI;
        for(i = 2; i <= tnn/4+1; i++)
        {
            i1 = i+i-2;
            i2 = i1+1;
            i3 = tnn+1-i2;
            i4 = i3+1;
            wrs = twr;
            wis = twi;


            h1r = c1*(a[i1]+a[i3]);
            h1i = c1*(a[i2]-a[i4]);
            h2r = -c2*(a[i2]+a[i4]);
            h2i = c2*(a[i1]-a[i3]);
            a[i1] = h1r+wrs*h2r-wis*h2i;
            a[i2] = h1i+wrs*h2i+wis*h2r;
            a[i3] = h1r-wrs*h2r+wis*h2i;
            a[i4] = -h1i+wrs*h2i+wis*h2r;
            twtemp = twr;
            twr = twr*twpr-twi*twM_PI+twr;
            twi = twi*twpr+twtemp*twM_PI+twi;
        }
        h1r = a[0];
        a[0] = c1*(h1r+a[1]);
        a[1] = c1*(h1r-a[1]);
    }
    if( inversefft )
    {
        isign = -1;
    }
    else
    {
        isign = 1;
    }
    n = tnn;
    nn = tnn/2;
    j = 1;
    for(ii = 1; ii <= nn; ii++)
    {
        i = 2*ii-1;
        if( j>i )
        {
            tempr = a[j-1];
            temM_PI = a[j];
            a[j-1] = a[i-1];
            a[j] = a[i];
            a[i-1] = tempr;
            a[i] = temM_PI;
        }
        m = n/2;
        while(m>=2&&j>m)
        {
            j = j-m;
            m = m/2;
        }
        j = j+m;
    }
    mmax = 2;
    while(n>mmax)
    {
        istep = 2*mmax;
        theta = 2*M_PI/(isign*mmax);
  ibidon = sin(0.5*theta);
        wpr = -2.0*ibidon*ibidon;
    //    wpr = -2.0*sqr(sin(0.5*theta));
        wM_PI = sin(theta);
        wr = 1.0;
        wi = 0.0;
        for(ii = 1; ii <= mmax/2; ii++)
        {
            m = 2*ii-1;
            for(jj = 0; jj <= (n-m)/istep; jj++)
            {
                i = m+jj*istep;
                j = i+mmax;
                tempr = wr*a[j-1]-wi*a[j];
                temM_PI = wr*a[j]+wi*a[j-1];
                a[j-1] = a[i-1]-tempr;
                a[j] = a[i]-temM_PI;
                a[i-1] = a[i-1]+tempr;
                a[i] = a[i]+temM_PI;
            }
            wtemp = wr;
            wr = wr*wpr-wi*wM_PI+wr;
            wi = wi*wpr+wtemp*wM_PI+wi;
        }
        mmax = istep;
    }
    if( inversefft )
    {
        for(i = 1; i <= 2*nn; i++)
        {
            a[i-1] = a[i-1]/nn;
        }
    }
    if( !inversefft )
    {
        ibidon = sin(0.5*ttheta);
        twpr = -2.0*ibidon*ibidon;
 // twpr = -2.0*sqr(sin(0.5*ttheta));
        twM_PI = sin(ttheta);
        twr = 1.0+twpr;
        twi = twM_PI;
        for(i = 2; i <= tnn/4+1; i++)
        {
            i1 = i+i-2;
            i2 = i1+1;
            i3 = tnn+1-i2;
            i4 = i3+1;
            wrs = twr;
            wis = twi;
            h1r = c1*(a[i1]+a[i3]);
            h1i = c1*(a[i2]-a[i4]);
            h2r = -c2*(a[i2]+a[i4]);
            h2i = c2*(a[i1]-a[i3]);
            a[i1] = h1r+wrs*h2r-wis*h2i;
            a[i2] = h1i+wrs*h2i+wis*h2r;
            a[i3] = h1r-wrs*h2r+wis*h2i;
            a[i4] = -h1i+wrs*h2i+wis*h2r;
            twtemp = twr;
            twr = twr*twpr-twi*twM_PI+twr;
            twi = twi*twpr+twtemp*twM_PI+twi;
        }
        h1r = a[0];
        a[0] = h1r+a[1];
        a[1] = h1r-a[1];
    }
}



