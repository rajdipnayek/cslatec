#ifndef F2C_INCLUDE
#define F2C_INCLUDE
#include <math.h>
#include <string.h>
#ifdef __cplusplus
extern "C" {
#endif
#define TRUE_ (1)
#define FALSE_ (0)

typedef int integer;
typedef char *address;
typedef float real;
typedef double doublereal;
typedef struct { float r, i; } complex;
typedef int logical;

#ifdef f2c_i2
/* for -i2 */
typedef short ftnlen;
#else
typedef int ftnlen;
#endif

#ifdef __cplusplus
#define unknown_params ...
#else
#define unknown_params
#endif
typedef int    (*U_fp)(unknown_params), (*S_fp)(unknown_params);
typedef double (*D_fp)(unknown_params), (*E_fp)(unknown_params);
#undef unknown_params

/* I/O stuff */

/*external read, write*/
typedef struct
{	ftnlen cierr;
	ftnlen ciunit;
	ftnlen ciend;
	char *cifmt;
	ftnlen cirec;
} cilist;

/*internal read, write*/
typedef struct
{	ftnlen icierr;
	char *iciunit;
	ftnlen iciend;
	char *icifmt;
	ftnlen icirlen;
	ftnlen icirnum;
} icilist;

/*open*/
typedef struct
{	ftnlen oerr;
	ftnlen ounit;
	char *ofnm;
	ftnlen ofnmlen;
	char *osta;
	char *oacc;
	char *ofm;
	ftnlen orl;
	char *oblnk;
} olist;

/*close*/
typedef struct
{	ftnlen cerr;
	ftnlen cunit;
	char *csta;
} cllist;

/*rewind, backspace, endfile*/
typedef struct
{	ftnlen aerr;
	ftnlen aunit;
} alist;

#define abs(x) ((x) >= 0 ? (x) : -(x))
#define dabs(x) (double)abs(x)
#define min(a,b) ((a) <= (b) ? (a) : (b))
#define max(a,b) ((a) >= (b) ? (a) : (b))
#define dmin(a,b) (double)min(a,b)
#define dmax(a,b) (double)max(a,b)
#define bit_test(a,b)	((a) >> (b) & 1)
#define bit_clear(a,b)	((a) & ~((unsigned)1 << (b)))
#define bit_set(a,b)	((a) |  ((unsigned)1 << (b)))

double f__cabs(double, double);
char *F77_aloc(int Len, const char *whence);
void sig_die(const char*, int);
void _uninit_f2c(void *x, int type, long len);

/*
 * Public functions in libF77
 */

void c_cos(complex *r, complex *z);
void c_div(complex *c, complex *a, complex *b);
void c_exp(complex *r, complex *z);
void c_log(complex *r, complex *z);
void c_sin(complex *r, complex *z);
void c_sqrt(complex *r, complex *z);

int i_indx(char *a, char *b, ftnlen la, ftnlen lb);

int lbit_bits(int a, int b, int len);
int lbit_shift(int a, int b);
int lbit_cshift(int a, int b, int len);

void pow_ci(complex *p, complex *a, int *b);
double pow_dd(double *ap, double *bp);
double pow_di(double *ap, int *bp);
int pow_ii(int *ap, int *bp);
double pow_ri(float *ap, int *bp);

double r_abs(float *x);
double r_acos(float *x);
double r_asin(float *x);
double r_atan(float *x);
double r_atn2(float *x, float *y);
void r_cnjg(complex *r, complex *z);
double r_cos(float *x);
double r_cosh(float *x);
double r_dim(float *a, float *b);
double r_exp(float *x);
double r_imag(complex *z);
double r_int(float *x);
double r_lg10(float *x);
double r_log(float *x);
double r_mod(float *x, float *y);
double r_nint(float *x);
double r_sign(float *a, float *b);
double r_sin(float *x);
double r_sinh(float *x);
double r_sqrt(float *x);
double r_tan(float *x);
double r_tanh(float *x);

int s_cat(char *lp, char *rpp[], ftnlen rnp[], ftnlen *np, ftnlen ll);
int s_cmp(const char *a0, const char *b0, ftnlen la, ftnlen lb);
int s_paus(char *s, ftnlen n);
int s_stop(char *s, ftnlen n);

/* inline functions */

static inline double c_abs(const complex *z) { return hypot(z->r, z->i); }
static inline double d_abs(const double *x) { return fabs(*x); }
static inline double d_acos(const double *x) { return acos(*x); }
static inline double d_acosh(const double *x) { return acosh(*x); }
static inline double d_asin(const double *x) { return asin(*x); }
static inline double d_asinh(const double *x) { return asinh(*x); }
static inline double d_atan(const double *x) { return atan(*x); }
static inline double d_atanh(const double *x) { return atanh(*x); }
static inline double d_atn2(const double *x, double *y) { return atan2(*x, *y); }

static inline double d_cos(const double *x) { return cos(*x); }
static inline double d_cosh(const double *x) { return cosh(*x); }

static inline double d_dim(const double *a, double *b)
{
  double d = (*a - *b);
  return (d > 0)? d : 0;
}

static inline double d_exp(const double *x) { return exp(*x); }

static inline double d_int(const double *x) {
  double y = *x;
  return (y < 0)? floor(y) : -floor(-y);
}

static inline double d_lg10(const double *x) { return log10(*x); }

static inline double d_log(const double *x) { return log(*x); }
static inline double d_nint(const double *x) { return round(*x); }
static inline double d_prod(const float *x, const float *y) { return ((double)*x) * ((double)*x); }
static inline double d_sin(const double *x) { return sin(*x); }
static inline double d_tan(const double *x) { return tan(*x); }
static inline double d_sinh(const double *x) { return sinh(*x); }
static inline double d_sqrt(const double *x) { return sqrt(*x); }
static inline double d_tanh(const double *x) { return tanh(*x); }

static inline double d_sign(const double *a, const double *b)
{
  double x = fabs(*a);
  return (*b >= 0 ? x : -x);
}

static inline double derfc_(const double *x) { return erfc(*x); }
static inline double derf_(const double *x) { return erf(*x); }
static inline double erf_(const float *x) { return erf((double)(*x)); }
static inline double erfc_(const float *x) { return erfc((double)(*x)); }

static inline int i_abs(const int *x) { return abs(*x); }

static inline int i_dim(const int *a, const int *b)
{
  int d = (*a - *b);
  return (d > 0)? d : 0;
}

static inline int i_len(const char *s, ftnlen n) { return n; }
static inline int i_mod(const int *a, const int *b)
{
  return *a % *b;
}
static inline int i_nint(const float *x)
{
  return (int)round(*x);
}
static inline int i_dnnt(const double *x)
{
  return (int)round(*x);
}
static inline int i_sign(const int *a, const int *b)
{
  int x = abs(*a);
  return *b >= 0 ? x : -x;
}

static inline int s_copy(char *a, const char *b, ftnlen la, ftnlen lb)
{
  if (la <= lb) {
    memmove(a, b, la);
  } else {
    memmove(a, b, lb);
    memset(a, ' ', la - lb);
  }
  return 0;
}

static inline int i_sceiling(const float *r) {
  float x = *r;
  return ((int)(x) + ((x) > 0 && (x) != (int)(x)));
}

static inline int i_dceiling(const double *r) {
  double x = *r;
  return ((int)(x) + ((x) > 0 && (x) != (int)(x)));
}

/*************************************************************
 * LIBI77
 *
 * Public functions
 */

int c_dfe(cilist *a);
int c_due(cilist *a);
int c_sfe(cilist *a);
int c_sue(cilist *a);

int e_rdfe(void);
int e_rdue(void);
int e_rsfe(void);
int e_rsfi(void);
int e_rsle(void);
int e_rsli(void);
int e_rsue(void);
int e_wdfe(void);
int e_wdue(void);
int e_wsfi(void);
int e_wsfe(void);
int e_wsle(void);
int e_wsli(void);
int e_wsue(void);

void exit_(int *rc);

int f_back(alist *a);
int f_clos(cllist *a);
int f_end(alist *a);
void f_exit(void);
int f_open(olist *a);
int f_rew(alist *a);
int flush_(void);

int s_rdfe(cilist *a);
int s_rdue(cilist *a);
int s_rsfi(icilist *a);
int s_rsle(cilist *a);
int s_rsli(icilist *a);
int s_rsne(cilist *a);
int s_rsni(icilist *a);
int s_rsue(cilist *a);
int s_wdfe(cilist *a);
int s_wdue(cilist *a);
int s_wsfe(cilist *a);
int s_wsfi(icilist *a);
int s_wsle(cilist *a);
int s_wsli(icilist *a);
int s_wsne(cilist *a);
int s_wsni(icilist *a);
int s_wsue(cilist *a);

#ifdef __cplusplus
}
#endif
#endif
