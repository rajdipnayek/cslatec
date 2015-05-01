#include <math.h>
#include "math_c99.h"

#ifdef __cplusplus
extern "C" {
#endif

double d_mod(const double *x, const double *y)
{
    return fmod(*x, *y);
}

void pow_ci(float_complex *r, const float_complex *z, const int *i)
{
    const double_complex s = cpow(
        cmplx(crealf(*z), cimagf(*z)),
        cmplx((double)*i, 0.)
    );
    *r = cmplxf(creal(s), cimag(s));
}

double pow_dd(const double *x, const double *y)
{
    return pow(*x, *y);
}

double pow_di(const double *x, const int *i)
{
    return pow(*x, (double)*i);
}

int pow_ii(const int *i, const int *j)
{
    return (int)pow((double)*i, (double)*j);
}

double pow_ri(const float *x, const int *i)
{
    return pow((double)*x, (double)*i);
}

void r_cnjg(float_complex *r, const float_complex *z)
{
    *r = conjf(*z);
}

double r_imag(const float_complex *z)
{
    return cimagf(*z);
}

double r_int(const float *x)
{
    return truncf(*x);
}

double r_lg10(const float *x)
{
    return log10f(*x);
}

double r_mod(const float *x, const float *y)
{
    return fmodf(*x, *y);
}

double r_sign(const float *x, const float *y)
{
    return copysignf(*x, *y);
}

void c_div(complex *r, const complex *z, const complex *w)
{
    *r = cdivf(*z, *w);
}

void c_cos(complex *r, const complex *z)
{
    *r = ccosf(*z);
}

void c_exp(complex *r, const complex *z)
{
    *r = cexpf(*z);
}

void c_log(complex *r, const complex *z)
{
    *r = clogf(*z);
}

void c_sin(complex *r, const complex *z)
{
    *r = csinf(*z);
}

void c_sqrt(complex *r, const complex *z)
{
    *r = csqrtf(*z);
}

#ifdef __cplusplus
}
#endif
