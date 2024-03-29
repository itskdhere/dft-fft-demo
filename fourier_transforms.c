#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <complex.h>
#include <math.h>
#include <time.h>

#define PI atan2f(1, 1) * 4
#define SAMPLE 4096

size_t n = SAMPLE;
float in[SAMPLE];
float complex out[SAMPLE];
// float out_sin[SAMPLE];
// float out_cos[SAMPLE];

void genSig(char sigType[4], float sigFreq)
{
    float inTmp[n];

    if (strcmp(sigType, "sin") == 0)
    {
        for (size_t i = 0; i < n; ++i)
        {
            float t = (float)i / (float)n;
            inTmp[i] = sin(2 * PI * t * sigFreq);
        }
    }
    else if (strcmp(sigType, "cos") == 0)
    {
        for (size_t i = 0; i < n; ++i)
        {
            float t = (float)i / (float)n;
            inTmp[i] = cos(2 * PI * t * sigFreq);
        }
    }
    else
    {
        printf("\nError: Invalid Signal Type!\n");
        exit(1);
    }

    for (size_t i = 0; i < n; ++i)
    {
        in[i] += inTmp[i];
    }
}

void dft(void)
{
    for (size_t f = 0; f < n; ++f)
    {
        out[f] = 0;
        // out_sin[f] = 0;
        // out_cos[f] = 0;
        for (size_t i = 0; i < n; ++i)
        {
            float t = (float)i / (float)n;
            out[f] += in[i] * cexp(2 * I * PI * t * f);
            // out_sin[f] += in[i] * sinf(2 * PI * t * f);
            // out_cos[f] += in[i] * cosf(2 * PI * t * f);
        }
    }
}

void fft(float in[], float complex out[], size_t n, size_t stride)
{
    if (n <= 1)
    {
        out[0] = in[0];
        return;
    }
    fft(in, out, n / 2, stride * 2);
    fft(in + stride, out + n / 2, n / 2, stride * 2);

    for (size_t k = 0; k < n / 2; ++k)
    {
        float t = (float)k / (float)n;
        float complex v = cexp(2 * I * PI * t) * out[k + n / 2];
        float complex e = out[k];
        out[k] = e + v;
        out[k + n / 2] = e - v;
    }
}

void display(float complex sig[])
{
    printf("\nSin Component(s): ");
    for (size_t f = 0; f < n; ++f)
    {
        float im = cimagf(sig[f]);
        float re = crealf(sig[f]);
        if (im >= 1)
            printf("%02zu Hz\t", f);
        // printf("%02zu\t%0.3f\t%0.3f\n", f, im, re);
    }

    printf("\nCos Component(s): ");
    for (size_t f = 0; f < n; ++f)
    {
        float im = cimagf(sig[f]);
        float re = crealf(sig[f]);
        if (re >= 1)
            printf("%02zu Hz\t", f);
        // printf("%02zu\t%0.3f\t%0.3f\n", f, im, re);
    }
}

int main(void)
{
    clock_t start, end;
    int sigNum;

    printf("\nNumber Of Signals To Be Generated: ");
    scanf("%d", &sigNum);

    for (int i = 1; i <= sigNum; i++)
    {
        char sigType[4];
        float sigFreq;
        fflush(stdin);
        printf("\nSignal %d Type ('sin' or 'cos'): ", i);
        fgets(sigType, 4, stdin);
        printf("Signal %d Frequency (in Hz): ", i);
        scanf("%f", &sigFreq);
        genSig(sigType, sigFreq);
    }

    printf("\n\nComputing Fourier Transforms...");
    start = clock();
    dft();
    end = clock();
    double cpu_time_used_dft = ((double)(end - start)) / CLOCKS_PER_SEC;

    printf("\n\nDiscrete Fourier Transform:");
    display(out);

    start = clock();
    fft(in, out, n, 1);
    end = clock();
    double cpu_time_used_fft = ((double)(end - start)) / CLOCKS_PER_SEC;

    printf("\n\nFast Fourier Transform:");
    display(out);

    printf("\n\n\nTime Taken For DFT [O(N^2)]: %f\n", cpu_time_used_dft);
    printf("Time Taken For FFT [O(n log(n))]: %f\n", cpu_time_used_fft);
    return 0;
}
