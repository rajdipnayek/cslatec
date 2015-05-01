#ifndef F2C_INCLUDE
#define F2C_INCLUDE
#include <math.h>
#include <string.h>
#ifdef __cplusplus
extern "C" {
#endif
#define TRUE_ (1)
#define FALSE_ (0)

typedef int integer, logical, ftnlen;
typedef char *address;
typedef float real;
typedef double doublereal;
typedef struct { float r, i; } complex;

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
typedef struct {
    int cierr, ciunit, ciend;
    char *cifmt;
    int cirec;
} cilist;

/*internal read, write*/
typedef struct {
    int icierr;
    char *iciunit;
    int iciend;
    char *icifmt;
    int icirlen, icirnum;
} icilist;

/*open*/
typedef struct {
    int oerr, ounit;
    char *ofnm;
    int ofnmlen;
    char *osta, *oacc, *ofm;
    int orl;
    char *oblnk;
} olist;

/*close*/
typedef struct {
    int cerr, cunit;
    char *csta;
} cllist;

/*rewind, backspace, endfile*/
typedef struct {
    int aerr, aunit;
} alist;

#define abs(x) ((x) >= 0 ? (x) : -(x))
#define dabs(x) (double)abs(x)
#define min(a,b) ((a) <= (b) ? (a) : (b))
#define max(a,b) ((a) >= (b) ? (a) : (b))
#define dmin(a,b) (double)min(a,b)
#define dmax(a,b) (double)max(a,b)

void c_cos(complex *r, const complex *z);
void c_div(complex *c, const complex *a, const complex *b);
void c_exp(complex *r, const complex *z);
void c_log(complex *r, const complex *z);
void c_sin(complex *r, const complex *z);
void c_sqrt(complex *r, const complex *z);

int i_indx(char *a, char *b, int la, int lb);

void pow_ci(complex *p, const complex *a, const int *b);
double pow_dd(const double *ap, const double *bp);
double pow_di(const double *ap, const int *bp);
int pow_ii(const int *ap, const int *bp);
double pow_ri(const float *ap, const int *bp);

double r_abs(const float *x);
double r_acos(const float *x);
double r_asin(const float *x);
double r_atan(const float *x);
double r_atn2(const float *x, const float *y);
void r_cnjg(complex *r, const complex *z);
double r_cos(const float *x);
double r_cosh(const float *x);
double r_dim(float *a, float *b);
double r_exp(const float *x);
double r_imag(const complex *z);
double r_int(const float *x);
double r_lg10(const float *x);
double r_log(const float *x);
double r_mod(const float *x, const float *y);
double r_nint(const float *x);
double r_sign(const float *a, const float *b);
double r_sin(const float *x);
double r_sinh(const float *x);
double r_sqrt(const float *x);
double r_tan(const float *x);
double r_tanh(const float *x);

double d_mod(const double *x, const double *y);

int s_cat(char *lp, char **rpp, int *rnp, int *np, int ll);
int s_cmp(const char *a0, const char *b0, int la, int lb);
int s_stop(char *s, int n);

int e_rdue(void);
int e_rsfe(void);
int e_rsue(void);
int e_wdue(void);
int e_wsfe(void);
int e_wsfi(void);
int e_wsle(void);
int e_wsue(void);

int f_clos(cllist *a);
int f_end(alist *a);
int f_open(olist *a);
int f_rew(alist *a);

int s_rdue(cilist *a);
int s_rsfe(cilist *a);
int s_rsue(cilist *a);
int s_wdue(cilist *a);
int s_wsfe(cilist *a);
int s_wsfi(icilist *a);
int s_wsle(cilist *a);
int s_wsue(cilist *a);

int do_fio(int *number, char *ptr, int len);
int do_uio(int *number, char *ptr, int len);

/* inline functions */

static inline double c_abs(const complex *z) { return hypot(z->r, z->i); }
static inline double d_int(const double *x) {
    const double y = *x;
    return y < 0 ? floor(y) : -floor(-y);
}
static inline double d_lg10(const double *x) { return log10(*x); }
static inline double d_sign(const double *a, const double *b)
{
    const double x = fabs(*a);
    return *b >= 0 ? x : -x;
}

static inline double derfc_(const double *x) { return erfc(*x); }
static inline double derf_(const double *x) { return erf(*x); }
static inline double erf_(const float *x) { return erf((double)(*x)); }
static inline double erfc_(const float *x) { return erfc((double)(*x)); }
static inline int i_len(const char *s, int n) { return n; }
static inline int i_sign(const int *a, const int *b)
{
    const int x = abs(*a);
    return *b >= 0 ? x : -x;
}

static inline int s_copy(char *a, const char *b, int la, int lb)
{
    if (la <= lb) {
        memmove(a, b, la);
    } else {
        memmove(a, b, lb);
        memset(a, ' ', la - lb);
    }
    return 0;
}

static inline int i_sceiling(const float *r)
{
    return ((int)*r + (*r > 0 && *r != (int)*r));
}

static inline int i_dceiling(const double *r)
{
    return ((int)*r + (*r > 0 && *r != (int)*r));
}

#ifdef __cplusplus
}
#endif
#endif
