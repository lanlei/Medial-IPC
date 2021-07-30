#pragma once
#ifndef DATA_CORE_H
#define DATA_CORE_H

#include <stdlib.h>
#include <math.h>

#define USE_DOUBLE_PRECISION

#ifdef USE_DOUBLE_PRECISION
typedef double qeal;
#define  GL_QEAL GL_DOUBLE
#define QEAL_ZERO 5e-14
#define MIN_VALUE 5e-14
#define QEAL_MIN DBL_MIN
#define QEAL_MAX DBL_MAX
#define IS_QEAL_ZERO(d) (abs(d) < MIN_VALUE)
#define Check_QEAL_ZERO(d) (abs(d) > MIN_VALUE ? d : 0)
#else  
typedef float qeal;
#define  GL_QEAL GL_QEAL
#define MIN_VALUE 3e-5
#define QEAL_MIN FLT_MIN
#define QEAL_MAX FLT_MAX
#define IS_QEAL_ZERO(d) (fabs(d) < MIN_VALUE)
#endif

#define IS_DOUBLE_ZERO(d) (abs(d) < MIN_VALUE)
#define IS_FLOAT_ZERO(d) (fabs(d) < 3e-5)

#define max2(x,y) ((x)>(y))?(x):(y);
#define min2(x,y) ((x)<(y))?(x):(y);
#define  max3(a,b,c)  (a>b?(a>c?a:c):(b>c?b:c))
#define  min3(a,b,c)  (a<b?(a<c?a:c):(b<c?b:c))

#define isInValidRange(x, l, h)(x > l ? (x <= h ? 1:0) : (0))

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define NORROW 0

#endif
