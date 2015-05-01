#ifndef F2C_TYPES_H
#define F2C_TYPES_H

typedef int integer;
typedef char *address;
typedef float real;
typedef double doublereal;
typedef struct { float r, i; } complex;
typedef int logical;

#ifdef f2c_i2
/* for -i2 */
typedef short flag;
#else
typedef int flag;
#endif
typedef flag ftnlen;

#ifdef __cplusplus
#define unknown_params ...
#else
#define unknown_params
#endif
typedef int (*U_fp)(unknown_params);
typedef integer (*I_fp)(unknown_params);
typedef real (*R_fp)(unknown_params);
typedef double (*D_fp)(unknown_params), (*E_fp)(unknown_params);
typedef void (*C_fp)(unknown_params);
typedef void (*Z_fp)(unknown_params);
typedef logical (*L_fp)(unknown_params);
typedef void (*H_fp)(unknown_params);
typedef int (*S_fp)(unknown_params);
#undef unknown_params

#endif
