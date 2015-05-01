#ifndef NO_COMPLEX
#include <complex.h>

typedef float  _Complex float_complex;
typedef double _Complex double_complex;

static inline float_complex cdivf(float_complex z, float_complex w)
{
    return z / w;
}

static inline float_complex cmplxf(float x, float y)
{
    return x + y * _Complex_I;
}

static inline double_complex cmplx(double x, double y)
{
    return x + y * _Complex_I;
}

#else

typedef struct { float real, imag; } float_complex;

static float_complex cmplxf(float x, float y)
{
    const float_complex z = {x, y};
    return z;
}

static double_complex cmplx(double x, double y)
{
    const double_complex z = {x, y};
    return z;
}

static double creal(double_complex z)
{
    return z.real;
}

static float crealf(float_complex z)
{
    return z.real;
}

static double cimag(double_complex z)
{
    return z.imag;
}

static float cimagf(float_complex z)
{
    return z.imag;
}

static float_complex cdivf(float_complex z, float_complex w)
{
    /* not sure if this is 100% correct */
    const float w2 = w.real * w.real + w.imag * w.imag
    const float_complex r = {
        (z.real * w.real + z.imag * w.imag) / w2,
        (z.imag * w.real - z.real * w.imag) / w2
    };
    return r;
}

static float_complex conjf(float_complex z)
{
    const float_complex r = {z.real, -z.imag};
    return r;
}

#endif

#ifdef NO_COPYSIGN
static double copysign(double x, double y)
{
    return fabs(x) * (y >= 0 ? 1 : -1);
}
#endif

#ifdef NO_COPYSIGNF
static float copysignf(float x, float y)
{
    return (float)copysign((double)x, (double)y);
}
#endif

#ifdef NO_FMODF
static float fmodf(float x, float y)
{
    return (float)fmod((double)x, (double)y);
}
#endif

#ifdef NO_FLOORF
static float floorf(float x)
{
    return (float)floor((double)x);
}
#endif

#ifdef NO_LOG10F
static float log10f(float x)
{
    return (float)log10((double)x);
}
#endif

#ifdef NO_TRUNC
static double trunc(double x)
{
    return x >= 0 ? floorf(x) : -floor(-x);
}
#endif

#ifdef NO_TRUNCF
static float truncf(float x)
{
    return (float)trunc((double)x);
}
#endif
