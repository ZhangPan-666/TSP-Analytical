#pragma once
#ifndef TSP_CORE_INCLUDED
#define TSP_CORE_INCLUDED
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <type_traits>
#include <vector>
/*
QD 2.3.24 double-double subset, consolidated with local includes.
Arithmetic and local IEEE/std::fma configuration are retained.

BSD 3-Clause License

Copyright (c) 2021, highprecision
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its
   contributors may be used to endorse or promote products derived from
   this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

This work was supported by the Director, Office of Science, Division
of Mathematical, Information, and Computational Sciences of the
U.S. Department of Energy under contract numbers DE-AC03-76SF00098 and
DE-AC02-05CH11231.

Copyright (c) 2003-2023, The Regents of the University of California,
through Lawrence Berkeley National Laboratory (subject to receipt of
any required approvals from U.S. Dept. of Energy) All rights reserved. 

By downloading or using this software you are agreeing to the modified
BSD license that is in file "LICENSE" in the main directory.
If you wish to use the software for commercial purposes
please contact the Technology Transfer Department at TTD@lbl.gov or
call 510-286-6457.


*/
#include <cmath>
#define QD_API
#define QD_INLINE 1
#define QD_HAVE_STD 1
#define QD_IEEE_ADD 1
#define QD_ISFINITE std::isfinite
#define QD_ISINF std::isinf
#define QD_ISNAN std::isnan
#define QD_FMA(a,b,c) std::fma(a,b,c)
#define QD_FMS(a,b,c) std::fma(a,b,-(c))

/*
 * include/inline.h
 *
 * This work was supported by the Director, Office of Science, Division
 * of Mathematical, Information, and Computational Sciences of the
 * U.S. Department of Energy under contract number DE-AC03-76SF00098.
 *
 * Copyright (c) 2000-2001
 *
 * This file contains the basic functions used both by double-double
 * and quad-double package.  These are declared as inline functions as
 * they are the smallest building blocks of the double-double and 
 * quad-double arithmetic.
 */
#ifndef _QD_INLINE_H
#define _QD_INLINE_H

#define _QD_SPLITTER 134217729.0               // = 2^27 + 1
#define _QD_SPLIT_THRESH 6.69692879491417e+299 // = 2^996

#ifdef QD_VACPP_BUILTINS_H
/* For VisualAge C++ __fmadd */
#include <builtins.h>
#endif

#include <cmath>
#include <limits>

namespace qd {

inline const double _d_nan = std::numeric_limits<double>::quiet_NaN();
inline const double _d_inf = std::numeric_limits<double>::infinity();

/*********** Basic Functions ************/
/* Computes fl(a+b) and err(a+b).  Assumes |a| >= |b|. */
inline double quick_two_sum(double a, double b, double &err) {
  double s = a + b;
  err = b - (s - a);
  return s;
}

/* Computes fl(a-b) and err(a-b).  Assumes |a| >= |b| */
inline double quick_two_diff(double a, double b, double &err) {
  double s = a - b;
  err = (a - s) - b;
  return s;
}

/* Computes fl(a+b) and err(a+b).  */
inline double two_sum(double a, double b, double &err) {
  double s = a + b;
  double bb = s - a;
  err = (a - (s - bb)) + (b - bb);
  return s;
}

/* Computes fl(a-b) and err(a-b).  */
inline double two_diff(double a, double b, double &err) {
  double s = a - b;
  double bb = s - a;
  err = (a - (s - bb)) - (b + bb);
  return s;
}

#ifndef QD_FMS
/* Computes high word and lo word of a */
inline void split(double a, double &hi, double &lo) {
  double temp;
  if (a > _QD_SPLIT_THRESH || a < -_QD_SPLIT_THRESH) {
    a *= 3.7252902984619140625e-09;  // 2^-28
    temp = _QD_SPLITTER * a;
    hi = temp - (temp - a);
    lo = a - hi;
    hi *= 268435456.0;          // 2^28
    lo *= 268435456.0;          // 2^28
  } else {
    temp = _QD_SPLITTER * a;
    hi = temp - (temp - a);
    lo = a - hi;
  }
}
#endif

/* Computes fl(a*b) and err(a*b). */
inline double two_prod(double a, double b, double &err) {
#ifdef QD_FMS
  double p = a * b;
  err = QD_FMS(a, b, p);
  return p;
#else
  double a_hi, a_lo, b_hi, b_lo;
  double p = a * b;
  split(a, a_hi, a_lo);
  split(b, b_hi, b_lo);
  err = ((a_hi * b_hi - p) + a_hi * b_lo + a_lo * b_hi) + a_lo * b_lo;
  return p;
#endif
}

/* Computes fl(a*a) and err(a*a).  Faster than the above method. */
inline double two_sqr(double a, double &err) {
#ifdef QD_FMS
  double p = a * a;
  err = QD_FMS(a, a, p);
  return p;
#else
  double hi, lo;
  double q = a * a;
  split(a, hi, lo);
  err = ((hi * hi - q) + 2.0 * hi * lo) + lo * lo;
  return q;
#endif
}

/* Computes the nearest integer to d. */
inline double nint(double d) {
  if (d == std::floor(d))
    return d;
  return std::floor(d + 0.5);
}

/* Computes the truncated integer. */
inline double aint(double d) {
  return (d >= 0.0) ? std::floor(d) : std::ceil(d);
}

/* These are provided to give consistent 
   interface for double with double-double and quad-double. */
inline void sincosh(double t, double &sinh_t, double &cosh_t) {
  sinh_t = std::sinh(t);
  cosh_t = std::cosh(t);
}

inline double sqr(double t) {
  return t * t;
}

inline double to_double(double a) { return a; }
inline int    to_int(double a) { return static_cast<int>(a); }

}

#endif /* _QD_INLINE_H */

/*
 * include/fpu.h
 *
 * This work was supported by the Director, Office of Science, Division
 * of Mathematical, Information, and Computational Sciences of the
 * U.S. Department of Energy under contract number DE-AC03-76SF00098.
 *
 * Copyright (c) 2001
 *
 * Contains functions to set and restore the round-to-double flag in the
 * control word of a x86 FPU.  The algorithms in the double-double and
 * quad-double package does not function with the extended mode found in
 * these FPU.
 */
#ifndef _QD_FPU_H
#define _QD_FPU_H


#ifdef __cplusplus
extern "C" {
#endif

/*
 * Set the round-to-double flag, and save the old control word in old_cw.
 * If old_cw is NULL, the old control word is not saved.
 */
QD_API void fpu_fix_start(unsigned int *old_cw);

/*
 * Restore the control word.
 */
QD_API void fpu_fix_end(unsigned int *old_cw);

#ifdef __cplusplus
}
#endif

#endif  /* _QD_FPU_H */

/*
 * include/bits.h
 *
 * This work was supported by the Director, Office of Science, Division
 * of Mathematical, Information, and Computational Sciences of the
 * U.S. Department of Energy under contract number DE-AC03-76SF00098.
 *
 * Copyright (c) 2000-2001
 *
 * This file defines various routines to get / set bits of a IEEE floating
 * point number.  This is used by the library for debugging purposes.
 */

#ifndef _QD_BITS_H
#define _QD_BITS_H

#include <iostream>

/* Returns the exponent of the double precision number.
   Returns INT_MIN is x is zero, and INT_MAX if x is INF or NaN. */
int get_double_expn(double x);

/* Prints 
     SIGN  EXPN  MANTISSA
   of the given double.  If x is NaN, INF, or Zero, this
   prints out the strings NaN, +/- INF, and 0.             */
void print_double_info(std::ostream &os, double x);


#endif  /* _QD_BITS_H */


/*
 * include/dd_real.h
 *
 * This work was supported by the Director, Office of Science, Division
 * of Mathematical, Information, and Computational Sciences of the
 * U.S. Department of Energy under contract number DE-AC03-76SF00098.
 *
 * Copyright (c) 2000-2007
 *
 * Double-double precision (>= 106-bit significand) floating point
 * arithmetic package based on David Bailey's Fortran-90 double-double
 * package, with some changes. See  
 *
 *   http://www.nersc.gov/~dhbailey/mpdist/mpdist.html
 *   
 * for the original Fortran-90 version.
 *
 * Overall structure is similar to that of Keith Brigg's C++ double-double
 * package.  See  
 *
 *   http://www-epidem.plansci.cam.ac.uk/~kbriggs/doubledouble.html
 *
 * for more details.  In particular, the fix for x86 computers is borrowed
 * from his code.
 *
 * Yozo Hida
 */

#ifndef _QD_DD_REAL_H
#define _QD_DD_REAL_H

#include <cmath>
#include <iostream>
#include <string>
#include <limits>

// Some compilers define isnan, isfinite, and isinf as macros, even for
// C++ codes, which cause havoc when overloading these functions.  We undef
// them here.
#ifdef isnan
#undef isnan
#endif

#ifdef isfinite
#undef isfinite
#endif

#ifdef isinf
#undef isinf
#endif

#ifdef max
#undef max
#endif

#ifdef min
#undef min
#endif

struct QD_API dd_real {
  double x[2];

  dd_real(double hi, double lo) { x[0] = hi; x[1] = lo; }
  dd_real() {x[0] = 0.0; x[1] = 0.0; }
  dd_real(double h) { x[0] = h; x[1] = 0.0; }
  dd_real(int h) {
    x[0] = (static_cast<double>(h));
    x[1] = 0.0;
  }

  dd_real (const char *s);
  explicit dd_real (const double *d) {
    x[0] = d[0]; x[1] = d[1];
  }

  static void error(const char *msg);

  double _hi() const { return x[0]; }
  double _lo() const { return x[1]; }

  static const dd_real _2pi;
  static const dd_real _pi;
  static const dd_real _3pi4;
  static const dd_real _pi2;
  static const dd_real _pi4;
  static const dd_real _e;
  static const dd_real _log2;
  static const dd_real _log10;
  static const dd_real _nan;
  static const dd_real _inf;

  static const double _eps;
  static const double _min_normalized;
  static const dd_real _max;
  static const dd_real _safe_max;
  static const int _ndigits;

  bool isnan() const { return QD_ISNAN(x[0]) || QD_ISNAN(x[1]); }
  bool isfinite() const { return QD_ISFINITE(x[0]); }
  bool isinf() const { return QD_ISINF(x[0]); }

  static dd_real add(double a, double b);
  static dd_real ieee_add(const dd_real &a, const dd_real &b);
  static dd_real sloppy_add(const dd_real &a, const dd_real &b);

  dd_real &operator+=(double a);
  dd_real &operator+=(const dd_real &a);

  static dd_real sub(double a, double b);

  dd_real &operator-=(double a);
  dd_real &operator-=(const dd_real &a);

  dd_real operator-() const;

  static dd_real mul(double a, double b);

  dd_real &operator*=(double a);
  dd_real &operator*=(const dd_real &a);

  static dd_real div(double a, double b);
  static dd_real sloppy_div(const dd_real &a, const dd_real &b);
  static dd_real accurate_div(const dd_real &a, const dd_real &b);
  
  dd_real &operator/=(double a);
  dd_real &operator/=(const dd_real &a);

  dd_real &operator=(double a);
  dd_real &operator=(const char *s);

  dd_real operator^(int n);
  static dd_real sqr(double d);

  static dd_real sqrt(double a);
  
  bool is_zero() const;
  bool is_one() const;
  bool is_positive() const;
  bool is_negative() const;

  static dd_real rand(void);

  void to_digits(char *s, int &expn, int precision = _ndigits) const;
  void write(char *s, int len, int precision = _ndigits, 
      bool showpos = false, bool uppercase = false) const;
  std::string to_string(int precision = _ndigits, int width = 0, 
      std::ios_base::fmtflags fmt = static_cast<std::ios_base::fmtflags>(0), 
      bool showpos = false, bool uppercase = false, char fill = ' ') const;
  int read(const char *s, dd_real &a);

  /* Debugging Methods */
  void dump(const std::string &name = "", std::ostream &os = std::cerr) const;
  void dump_bits(const std::string &name = "", 
                 std::ostream &os = std::cerr) const;

  static dd_real debug_rand();
};


namespace std {
  template <>
  class numeric_limits<dd_real> : public numeric_limits<double> {
  public:
    inline static double epsilon() { return dd_real::_eps; }
    inline static dd_real max() { return dd_real::_max; }
    inline static dd_real safe_max() { return dd_real::_safe_max; }
    inline static double min() { return dd_real::_min_normalized; }
    static const int digits = 104;
    static const int digits10 = 31;
  };
}

QD_API dd_real ddrand(void);
QD_API dd_real sqrt(const dd_real &a);

QD_API dd_real polyeval(const dd_real *c, int n, const dd_real &x);
QD_API dd_real polyroot(const dd_real *c, int n, 
    const dd_real &x0, int max_iter = 32, double thresh = 0.0);

QD_API inline bool isnan(const dd_real &a) { return a.isnan(); }
QD_API inline bool isfinite(const dd_real &a) { return a.isfinite(); }
QD_API inline bool isinf(const dd_real &a) { return a.isinf(); }

/* Computes  dd * d  where d is known to be a power of 2. */
QD_API dd_real mul_pwr2(const dd_real &dd, double d);

QD_API dd_real operator+(const dd_real &a, double b);
QD_API dd_real operator+(double a, const dd_real &b);
QD_API dd_real operator+(const dd_real &a, const dd_real &b);

QD_API dd_real operator-(const dd_real &a, double b);
QD_API dd_real operator-(double a, const dd_real &b);
QD_API dd_real operator-(const dd_real &a, const dd_real &b);

QD_API dd_real operator*(const dd_real &a, double b);
QD_API dd_real operator*(double a, const dd_real &b);
QD_API dd_real operator*(const dd_real &a, const dd_real &b);

QD_API dd_real operator/(const dd_real &a, double b);
QD_API dd_real operator/(double a, const dd_real &b);
QD_API dd_real operator/(const dd_real &a, const dd_real &b);

QD_API dd_real inv(const dd_real &a);

QD_API dd_real rem(const dd_real &a, const dd_real &b);
QD_API dd_real drem(const dd_real &a, const dd_real &b);
QD_API dd_real divrem(const dd_real &a, const dd_real &b, dd_real &r);

QD_API dd_real pow(const dd_real &a, int n);
QD_API dd_real pow(const dd_real &a, const dd_real &b);
QD_API dd_real npwr(const dd_real &a, int n);
QD_API dd_real sqr(const dd_real &a);

QD_API dd_real sqrt(const dd_real &a);
QD_API dd_real nroot(const dd_real &a, int n);

QD_API bool operator==(const dd_real &a, double b);
QD_API bool operator==(double a, const dd_real &b);
QD_API bool operator==(const dd_real &a, const dd_real &b);

QD_API bool operator<=(const dd_real &a, double b);
QD_API bool operator<=(double a, const dd_real &b);
QD_API bool operator<=(const dd_real &a, const dd_real &b);

QD_API bool operator>=(const dd_real &a, double b);
QD_API bool operator>=(double a, const dd_real &b);
QD_API bool operator>=(const dd_real &a, const dd_real &b);

QD_API bool operator<(const dd_real &a, double b);
QD_API bool operator<(double a, const dd_real &b);
QD_API bool operator<(const dd_real &a, const dd_real &b);

QD_API bool operator>(const dd_real &a, double b);
QD_API bool operator>(double a, const dd_real &b);
QD_API bool operator>(const dd_real &a, const dd_real &b);

QD_API bool operator!=(const dd_real &a, double b);
QD_API bool operator!=(double a, const dd_real &b);
QD_API bool operator!=(const dd_real &a, const dd_real &b);

QD_API dd_real nint(const dd_real &a);
QD_API dd_real floor(const dd_real &a);
QD_API dd_real ceil(const dd_real &a);
QD_API dd_real aint(const dd_real &a);

QD_API dd_real ddrand(void);

double to_double(const dd_real &a);
int    to_int(const dd_real &a);

QD_API dd_real exp(const dd_real &a);
QD_API dd_real ldexp(const dd_real &a, int exp);
QD_API dd_real log(const dd_real &a);
QD_API dd_real log10(const dd_real &a);

QD_API dd_real sin(const dd_real &a);
QD_API dd_real cos(const dd_real &a);
QD_API dd_real tan(const dd_real &a);
QD_API void sincos(const dd_real &a, dd_real &sin_a, dd_real &cos_a);

QD_API dd_real asin(const dd_real &a);
QD_API dd_real acos(const dd_real &a);
QD_API dd_real atan(const dd_real &a);
QD_API dd_real atan2(const dd_real &y, const dd_real &x);

QD_API dd_real sinh(const dd_real &a);
QD_API dd_real cosh(const dd_real &a);
QD_API dd_real tanh(const dd_real &a);
QD_API void sincosh(const dd_real &a, 
                      dd_real &sinh_a, dd_real &cosh_a);

QD_API dd_real asinh(const dd_real &a);
QD_API dd_real acosh(const dd_real &a);
QD_API dd_real atanh(const dd_real &a);

QD_API dd_real fabs(const dd_real &a);
QD_API dd_real abs(const dd_real &a);   /* same as fabs */

QD_API dd_real fmod(const dd_real &a, const dd_real &b);

QD_API std::ostream& operator<<(std::ostream &s, const dd_real &a);
QD_API std::istream& operator>>(std::istream &s, dd_real &a);
#ifdef QD_INLINE
#endif

#endif /* _QD_DD_REAL_H */


/*
 * include/dd_inline.h
 *
 * This work was supported by the Director, Office of Science, Division
 * of Mathematical, Information, and Computational Sciences of the
 * U.S. Department of Energy under contract number DE-AC03-76SF00098.
 *
 * Copyright (c) 2000-2001
 *
 * Contains small functions (suitable for inlining) in the double-double
 * arithmetic package.
 */
#ifndef _QD_DD_INLINE_H
#define _QD_DD_INLINE_H

#include <cmath>

#ifndef QD_INLINE
#define inline
#endif


/*********** Additions ************/
/* double-double = double + double */
inline dd_real dd_real::add(double a, double b) {
  double s, e;
  s = qd::two_sum(a, b, e);
  return dd_real(s, e);
}

/* double-double + double */
inline dd_real operator+(const dd_real &a, double b) {
  double s1, s2;
  s1 = qd::two_sum(a.x[0], b, s2);
  s2 += a.x[1];
  s1 = qd::quick_two_sum(s1, s2, s2);
  return dd_real(s1, s2);
}

/* double-double + double-double */
inline dd_real dd_real::ieee_add(const dd_real &a, const dd_real &b) {
  /* This one satisfies IEEE style error bound, 
     due to K. Briggs and W. Kahan.                   */
  double s1, s2, t1, t2;

  s1 = qd::two_sum(a.x[0], b.x[0], s2);
  t1 = qd::two_sum(a.x[1], b.x[1], t2);
  s2 += t1;
  s1 = qd::quick_two_sum(s1, s2, s2);
  s2 += t2;
  s1 = qd::quick_two_sum(s1, s2, s2);
  return dd_real(s1, s2);
}

inline dd_real dd_real::sloppy_add(const dd_real &a, const dd_real &b) {
  /* This is the less accurate version ... obeys Cray-style
     error bound. */
  double s, e;

  s = qd::two_sum(a.x[0], b.x[0], e);
  e += (a.x[1] + b.x[1]);
  s = qd::quick_two_sum(s, e, e);
  return dd_real(s, e);
}

inline dd_real operator+(const dd_real &a, const dd_real &b) {
#ifndef QD_IEEE_ADD
  return dd_real::sloppy_add(a, b);
#else
  return dd_real::ieee_add(a, b);
#endif
}

/* double + double-double */
inline dd_real operator+(double a, const dd_real &b) {
  return (b + a);
}


/*********** Self-Additions ************/
/* double-double += double */
inline dd_real &dd_real::operator+=(double a) {
  double s1, s2;
  s1 = qd::two_sum(x[0], a, s2);
  s2 += x[1];
  x[0] = qd::quick_two_sum(s1, s2, x[1]);
  return *this;
}

/* double-double += double-double */
inline dd_real &dd_real::operator+=(const dd_real &a) {
#ifndef QD_IEEE_ADD
  double s, e;
  s = qd::two_sum(x[0], a.x[0], e);
  e += x[1];
  e += a.x[1];
  x[0] = qd::quick_two_sum(s, e, x[1]);
  return *this;
#else
  double s1, s2, t1, t2;
  s1 = qd::two_sum(x[0], a.x[0], s2);
  t1 = qd::two_sum(x[1], a.x[1], t2);
  s2 += t1;
  s1 = qd::quick_two_sum(s1, s2, s2);
  s2 += t2;
  x[0] = qd::quick_two_sum(s1, s2, x[1]);
  return *this;
#endif
}

/*********** Subtractions ************/
/* double-double = double - double */
inline dd_real dd_real::sub(double a, double b) {
  double s, e;
  s = qd::two_diff(a, b, e);
  return dd_real(s, e);
}

/* double-double - double */
inline dd_real operator-(const dd_real &a, double b) {
  double s1, s2;
  s1 = qd::two_diff(a.x[0], b, s2);
  s2 += a.x[1];
  s1 = qd::quick_two_sum(s1, s2, s2);
  return dd_real(s1, s2);
}

/* double-double - double-double */
inline dd_real operator-(const dd_real &a, const dd_real &b) {
#ifndef QD_IEEE_ADD
  double s, e;
  s = qd::two_diff(a.x[0], b.x[0], e);
  e += a.x[1];
  e -= b.x[1];
  s = qd::quick_two_sum(s, e, e);
  return dd_real(s, e);
#else
  double s1, s2, t1, t2;
  s1 = qd::two_diff(a.x[0], b.x[0], s2);
  t1 = qd::two_diff(a.x[1], b.x[1], t2);
  s2 += t1;
  s1 = qd::quick_two_sum(s1, s2, s2);
  s2 += t2;
  s1 = qd::quick_two_sum(s1, s2, s2);
  return dd_real(s1, s2);
#endif
}

/* double - double-double */
inline dd_real operator-(double a, const dd_real &b) {
  double s1, s2;
  s1 = qd::two_diff(a, b.x[0], s2);
  s2 -= b.x[1];
  s1 = qd::quick_two_sum(s1, s2, s2);
  return dd_real(s1, s2);
}

/*********** Self-Subtractions ************/
/* double-double -= double */
inline dd_real &dd_real::operator-=(double a) {
  double s1, s2;
  s1 = qd::two_diff(x[0], a, s2);
  s2 += x[1];
  x[0] = qd::quick_two_sum(s1, s2, x[1]);
  return *this;
}

/* double-double -= double-double */
inline dd_real &dd_real::operator-=(const dd_real &a) {
#ifndef QD_IEEE_ADD
  double s, e;
  s = qd::two_diff(x[0], a.x[0], e);
  e += x[1];
  e -= a.x[1];
  x[0] = qd::quick_two_sum(s, e, x[1]);
  return *this;
#else
  double s1, s2, t1, t2;
  s1 = qd::two_diff(x[0], a.x[0], s2);
  t1 = qd::two_diff(x[1], a.x[1], t2);
  s2 += t1;
  s1 = qd::quick_two_sum(s1, s2, s2);
  s2 += t2;
  x[0] = qd::quick_two_sum(s1, s2, x[1]);
  return *this;
#endif
}

/*********** Unary Minus ***********/
inline dd_real dd_real::operator-() const {
  return dd_real(-x[0], -x[1]);
}

/*********** Multiplications ************/
/* double-double = double * double */
inline dd_real dd_real::mul(double a, double b) {
  double p, e;
  p = qd::two_prod(a, b, e);
  return dd_real(p, e);
}

/* double-double * (2.0 ^ exp) */
inline dd_real ldexp(const dd_real &a, int exp) {
  return dd_real(std::ldexp(a.x[0], exp), std::ldexp(a.x[1], exp));
}

/* double-double * double,  where double is a power of 2. */
inline dd_real mul_pwr2(const dd_real &a, double b) {
  return dd_real(a.x[0] * b, a.x[1] * b);
}

/* double-double * double */
inline dd_real operator*(const dd_real &a, double b) {
  double p1, p2;

  p1 = qd::two_prod(a.x[0], b, p2);
  p2 += (a.x[1] * b);
  p1 = qd::quick_two_sum(p1, p2, p2);
  return dd_real(p1, p2);
}

/* double-double * double-double */
inline dd_real operator*(const dd_real &a, const dd_real &b) {
  double p1, p2;

  p1 = qd::two_prod(a.x[0], b.x[0], p2);
  p2 += (a.x[0] * b.x[1] + a.x[1] * b.x[0]);
  p1 = qd::quick_two_sum(p1, p2, p2);
  return dd_real(p1, p2);
}

/* double * double-double */
inline dd_real operator*(double a, const dd_real &b) {
  return (b * a);
}

/*********** Self-Multiplications ************/
/* double-double *= double */
inline dd_real &dd_real::operator*=(double a) {
  double p1, p2;
  p1 = qd::two_prod(x[0], a, p2);
  p2 += x[1] * a;
  x[0] = qd::quick_two_sum(p1, p2, x[1]);
  return *this;
}

/* double-double *= double-double */
inline dd_real &dd_real::operator*=(const dd_real &a) {
  double p1, p2;
  p1 = qd::two_prod(x[0], a.x[0], p2);
  p2 += a.x[1] * x[0];
  p2 += a.x[0] * x[1];
  x[0] = qd::quick_two_sum(p1, p2, x[1]);
  return *this;
}

/*********** Divisions ************/
inline dd_real dd_real::div(double a, double b) {
  double q1, q2;
  double p1, p2;
  double s, e;

  q1 = a / b;

  /* Compute  a - q1 * b */
  p1 = qd::two_prod(q1, b, p2);
  s = qd::two_diff(a, p1, e);
  e -= p2;

  /* get next approximation */
  q2 = (s + e) / b;

  s = qd::quick_two_sum(q1, q2, e);

  return dd_real(s, e);
}

/* double-double / double */
inline dd_real operator/(const dd_real &a, double b) {

  double q1, q2;
  double p1, p2;
  double s, e;
  dd_real r;
  
  q1 = a.x[0] / b;   /* approximate quotient. */

  /* Compute  this - q1 * d */
  p1 = qd::two_prod(q1, b, p2);
  s = qd::two_diff(a.x[0], p1, e);
  e += a.x[1];
  e -= p2;
  
  /* get next approximation. */
  q2 = (s + e) / b;

  /* renormalize */
  r.x[0] = qd::quick_two_sum(q1, q2, r.x[1]);

  return r;
}

inline dd_real dd_real::sloppy_div(const dd_real &a, const dd_real &b) {
  double s1, s2;
  double q1, q2;
  dd_real r;

  q1 = a.x[0] / b.x[0];  /* approximate quotient */

  /* compute  this - q1 * dd */
  r = b * q1;
  s1 = qd::two_diff(a.x[0], r.x[0], s2);
  s2 -= r.x[1];
  s2 += a.x[1];

  /* get next approximation */
  q2 = (s1 + s2) / b.x[0];

  /* renormalize */
  r.x[0] = qd::quick_two_sum(q1, q2, r.x[1]);
  return r;
}

inline dd_real dd_real::accurate_div(const dd_real &a, const dd_real &b) {
  double q1, q2, q3;
  dd_real r;

  q1 = a.x[0] / b.x[0];  /* approximate quotient */

  r = a - q1 * b;
  
  q2 = r.x[0] / b.x[0];
  r -= (q2 * b);

  q3 = r.x[0] / b.x[0];

  q1 = qd::quick_two_sum(q1, q2, q2);
  r = dd_real(q1, q2) + q3;
  return r;
}

/* double-double / double-double */
inline dd_real operator/(const dd_real &a, const dd_real &b) {
#ifdef QD_SLOPPY_DIV
  return dd_real::sloppy_div(a, b);
#else
  return dd_real::accurate_div(a, b);
#endif
}

/* double / double-double */
inline dd_real operator/(double a, const dd_real &b) {
  return dd_real(a) / b;
}

inline dd_real inv(const dd_real &a) {
  return 1.0 / a;
}

/*********** Self-Divisions ************/
/* double-double /= double */
inline dd_real &dd_real::operator/=(double a) {
  *this = *this / a;
  return *this;
}

/* double-double /= double-double */
inline dd_real &dd_real::operator/=(const dd_real &a) {
  *this = *this / a;
  return *this;
}

/********** Remainder **********/
inline dd_real drem(const dd_real &a, const dd_real &b) {
  dd_real n = nint(a / b);
  return (a - n * b);
}

inline dd_real divrem(const dd_real &a, const dd_real &b, dd_real &r) {
  dd_real n = nint(a / b);
  r = a - n * b;
  return n;
}

/*********** Squaring **********/
inline dd_real sqr(const dd_real &a) {
  double p1, p2;
  double s1, s2;
  p1 = qd::two_sqr(a.x[0], p2);
  p2 += 2.0 * a.x[0] * a.x[1];
  p2 += a.x[1] * a.x[1];
  s1 = qd::quick_two_sum(p1, p2, s2);
  return dd_real(s1, s2);
}

inline dd_real dd_real::sqr(double a) {
  double p1, p2;
  p1 = qd::two_sqr(a, p2);
  return dd_real(p1, p2);
}


/********** Exponentiation **********/
inline dd_real dd_real::operator^(int n) {
  return npwr(*this, n);
}


/*********** Assignments ************/
/* double-double = double */
inline dd_real &dd_real::operator=(double a) {
  x[0] = a;
  x[1] = 0.0;
  return *this;
}

/*********** Equality Comparisons ************/
/* double-double == double */
inline bool operator==(const dd_real &a, double b) {
  return (a.x[0] == b && a.x[1] == 0.0);
}

/* double-double == double-double */
inline bool operator==(const dd_real &a, const dd_real &b) {
  return (a.x[0] == b.x[0] && a.x[1] == b.x[1]);
}

/* double == double-double */
inline bool operator==(double a, const dd_real &b) {
  return (a == b.x[0] && b.x[1] == 0.0);
}

/*********** Greater-Than Comparisons ************/
/* double-double > double */
inline bool operator>(const dd_real &a, double b) {
  return (a.x[0] > b || (a.x[0] == b && a.x[1] > 0.0));
}

/* double-double > double-double */
inline bool operator>(const dd_real &a, const dd_real &b) {
  return (a.x[0] > b.x[0] || (a.x[0] == b.x[0] && a.x[1] > b.x[1]));
}

/* double > double-double */
inline bool operator>(double a, const dd_real &b) {
  return (a > b.x[0] || (a == b.x[0] && b.x[1] < 0.0));
}

/*********** Less-Than Comparisons ************/
/* double-double < double */
inline bool operator<(const dd_real &a, double b) {
  return (a.x[0] < b || (a.x[0] == b && a.x[1] < 0.0));
}

/* double-double < double-double */
inline bool operator<(const dd_real &a, const dd_real &b) {
  return (a.x[0] < b.x[0] || (a.x[0] == b.x[0] && a.x[1] < b.x[1]));
}

/* double < double-double */
inline bool operator<(double a, const dd_real &b) {
  return (a < b.x[0] || (a == b.x[0] && b.x[1] > 0.0));
}

/*********** Greater-Than-Or-Equal-To Comparisons ************/
/* double-double >= double */
inline bool operator>=(const dd_real &a, double b) {
  return (a.x[0] > b || (a.x[0] == b && a.x[1] >= 0.0));
}

/* double-double >= double-double */
inline bool operator>=(const dd_real &a, const dd_real &b) {
  return (a.x[0] > b.x[0] || (a.x[0] == b.x[0] && a.x[1] >= b.x[1]));
}

/* double >= double-double */
inline bool operator>=(double a, const dd_real &b) {
  return (b <= a);
}

/*********** Less-Than-Or-Equal-To Comparisons ************/
/* double-double <= double */
inline bool operator<=(const dd_real &a, double b) {
  return (a.x[0] < b || (a.x[0] == b && a.x[1] <= 0.0));
}

/* double-double <= double-double */
inline bool operator<=(const dd_real &a, const dd_real &b) {
  return (a.x[0] < b.x[0] || (a.x[0] == b.x[0] && a.x[1] <= b.x[1]));
}

/* double <= double-double */
inline bool operator<=(double a, const dd_real &b) {
  return (b >= a);
}

/*********** Not-Equal-To Comparisons ************/
/* double-double != double */
inline bool operator!=(const dd_real &a, double b) {
  return (a.x[0] != b || a.x[1] != 0.0);
}

/* double-double != double-double */
inline bool operator!=(const dd_real &a, const dd_real &b) {
  return (a.x[0] != b.x[0] || a.x[1] != b.x[1]);
}

/* double != double-double */
inline bool operator!=(double a, const dd_real &b) {
  return (a != b.x[0] || b.x[1] != 0.0);
}

/*********** Micellaneous ************/
/*  this == 0 */
inline bool dd_real::is_zero() const {
  return (x[0] == 0.0);
}

/*  this == 1 */
inline bool dd_real::is_one() const {
  return (x[0] == 1.0 && x[1] == 0.0);
}

/*  this > 0 */
inline bool dd_real::is_positive() const {
  return (x[0] > 0.0);
}

/* this < 0 */
inline bool dd_real::is_negative() const {
  return (x[0] < 0.0);
}

/* Absolute value */
inline dd_real abs(const dd_real &a) {
  return (a.x[0] < 0.0) ? -a : a;
}

inline dd_real fabs(const dd_real &a) {
  return abs(a);
}

/* Round to Nearest integer */
inline dd_real nint(const dd_real &a) {
  double hi = qd::nint(a.x[0]);
  double lo;

  if (hi == a.x[0]) {
    /* High word is an integer already.  Round the low word.*/
    lo = qd::nint(a.x[1]);
    
    /* Renormalize. This is needed if x[0] = some integer, x[1] = 1/2.*/
    hi = qd::quick_two_sum(hi, lo, lo);
  } else {
    /* High word is not an integer. */
    lo = 0.0;
    if (std::abs(hi-a.x[0]) == 0.5 && a.x[1] < 0.0) {
      /* There is a tie in the high word, consult the low word 
         to break the tie. */
      hi -= 1.0;      /* NOTE: This does not cause INEXACT. */
    }
  }

  return dd_real(hi, lo);
}

inline dd_real floor(const dd_real &a) {
  double hi = std::floor(a.x[0]);
  double lo = 0.0;

  if (hi == a.x[0]) {
    /* High word is integer already.  Round the low word. */
    lo = std::floor(a.x[1]);
    hi = qd::quick_two_sum(hi, lo, lo);
  }

  return dd_real(hi, lo);
}

inline dd_real ceil(const dd_real &a) {
  double hi = std::ceil(a.x[0]);
  double lo = 0.0;

  if (hi == a.x[0]) {
    /* High word is integer already.  Round the low word. */
    lo = std::ceil(a.x[1]);
    hi = qd::quick_two_sum(hi, lo, lo);
  }

  return dd_real(hi, lo);
}

inline dd_real aint(const dd_real &a) {
  return (a.x[0] >= 0.0) ? floor(a) : ceil(a);
}

/* Cast to double. */
inline double to_double(const dd_real &a) {
  return a.x[0];
}

/* Cast to int. */
inline int to_int(const dd_real &a) {
  return static_cast<int>(a.x[0]);
}

/* Random number generator */
inline dd_real dd_real::rand() {
  return ddrand();
}

#endif /* _QD_DD_INLINE_H */

#include <string>

inline void append_expn(std::string &str, int expn);

/*
 * src/bits.cc
 *
 * This work was supported by the Director, Office of Science, Division
 * of Mathematical, Information, and Computational Sciences of the
 * U.S. Department of Energy under contract number DE-AC03-76SF00098.
 *
 * Copyright (c) 2000-2001
 *
 * Defines various routines to get / set bits of a IEEE floating point
 * number.  This used by the library for debugging purposes.
 */

#include <iostream>
#include <iomanip>
#include <cmath>
#include <climits>


#ifdef HAVE_IEEEFP_H
#include <ieeefp.h>
#endif

using std::setw;

inline int get_double_expn(double x) {
  if (x == 0.0)
    return INT_MIN;
  if (QD_ISINF(x) || QD_ISNAN(x))
    return INT_MAX;

  double y = std::abs(x);
  int i = 0;
  if (y < 1.0) {
    while (y < 1.0) {
      y *= 2.0;
      i++;
    }
    return -i;
  } else if (y >= 2.0) {
    while (y >= 2.0) {
      y *= 0.5;
      i++;
    }
    return i;
  }
  return 0;
}

inline void print_double_info(std::ostream &os, double x) {
  std::streamsize old_prec = os.precision(19);
  std::ios_base::fmtflags old_flags  = os.flags();
  os << std::scientific;

  os << setw(27) << x << ' ';
  if (QD_ISNAN(x) || QD_ISINF(x) || (x == 0.0)) {
    os << "                                                           ";
  } else {

    x = std::abs(x);
    int expn = get_double_expn(x);
    double d = std::ldexp(1.0, expn);
    os << setw(5) << expn << " ";
    for (int i = 0; i < 53; i++) {
      if (x >= d) {
        x -= d;
        os << '1';
      } else
        os << '0';
      d *= 0.5;
    }

    if (x != 0.0) {
      // should not happen
      os << " +trailing stuff";
    }
  }

  os.precision(old_prec);
  os.flags(old_flags);
}


/*
 * src/dd_const.cc
 *
 * This work was supported by the Director, Office of Science, Division
 * of Mathematical, Information, and Computational Sciences of the
 * U.S. Department of Energy under contract number DE-AC03-76SF00098.
 *
 * Copyright (c) 2000-2007
 */

inline const dd_real dd_real::_2pi = dd_real(6.283185307179586232e+00,
                                      2.449293598294706414e-16);
inline const dd_real dd_real::_pi = dd_real(3.141592653589793116e+00,
                                     1.224646799147353207e-16);
inline const dd_real dd_real::_pi2 = dd_real(1.570796326794896558e+00,
                                      6.123233995736766036e-17);
inline const dd_real dd_real::_pi4 = dd_real(7.853981633974482790e-01,
                                      3.061616997868383018e-17);
inline const dd_real dd_real::_3pi4 = dd_real(2.356194490192344837e+00,
                                       9.1848509936051484375e-17);
inline const dd_real dd_real::_e = dd_real(2.718281828459045091e+00,
                                    1.445646891729250158e-16);
inline const dd_real dd_real::_log2 = dd_real(6.931471805599452862e-01,
                                       2.319046813846299558e-17);
inline const dd_real dd_real::_log10 = dd_real(2.302585092994045901e+00,
                                        -2.170756223382249351e-16);
inline const dd_real dd_real::_nan = dd_real(qd::_d_nan, qd::_d_nan);
inline const dd_real dd_real::_inf = dd_real(qd::_d_inf, qd::_d_inf);

inline const double dd_real::_eps = 4.93038065763132e-32;  // 2^-104
inline const double dd_real::_min_normalized = 2.0041683600089728e-292;  // = 2^(-1022 + 53)
inline const dd_real dd_real::_max = 
    dd_real(1.79769313486231570815e+308, 9.97920154767359795037e+291);
inline const dd_real dd_real::_safe_max = 
    dd_real(1.7976931080746007281e+308, 9.97920154767359795037e+291);
inline const int dd_real::_ndigits = 31;



#include <cstdlib>

inline void append_expn(std::string &str, int expn) {
  int k;

  str += (expn < 0 ? '-' : '+');
  expn = std::abs(expn);

  if (expn >= 100) {
    k = (expn / 100);
    str += '0' + k;
    expn -= 100*k;
  }

  k = (expn / 10);
  str += '0' + k;
  expn -= 10*k;

  str += '0' + expn;
}


/*
 * src/fpu.cc
 *
 * This work was supported by the Director, Office of Science, Division
 * of Mathematical, Information, and Computational Sciences of the
 * U.S. Department of Energy under contract number DE-AC03-76SF00098.
 *
 * Copyright (c) 2000-2001
 *
 * Contains functions to set and restore the round-to-double flag in the
 * control word of a x86 FPU.
 */


#ifdef X86
#ifdef  _WIN32
#include <float.h>
#else

#ifdef HAVE_FPU_CONTROL_H
#include <fpu_control.h>
#endif

#ifndef _FPU_GETCW
#define _FPU_GETCW(x) asm volatile ("fnstcw %0":"=m" (x));
#endif

#ifndef _FPU_SETCW
#define _FPU_SETCW(x) asm volatile ("fldcw %0": :"m" (x));
#endif

#ifndef _FPU_EXTENDED
#define _FPU_EXTENDED 0x0300
#endif

#ifndef _FPU_DOUBLE
#define _FPU_DOUBLE 0x0200
#endif

#endif
#endif /* X86 */

extern "C" {

inline void fpu_fix_start(unsigned int *old_cw) {
#ifdef X86
#ifdef _WIN32
#ifdef __BORLANDC__
  /* Win 32 Borland C */
  unsigned short cw = _control87(0, 0);
  _control87(0x0200, 0x0300);
  if (old_cw) {
    *old_cw = cw;
  }
#else
  /* Win 32 MSVC */
  unsigned int cw = _control87(0, 0);
  _control87(0x00010000, 0x00030000);
  if (old_cw) {
    *old_cw = cw;
  }
#endif
#else
  /* Linux */
  volatile unsigned short cw, new_cw;
  _FPU_GETCW(cw);

  new_cw = (cw & ~_FPU_EXTENDED) | _FPU_DOUBLE;
  _FPU_SETCW(new_cw);
  
  if (old_cw) {
    *old_cw = cw;
  }
#endif
#endif
}

inline void fpu_fix_end(unsigned int *old_cw) {
#ifdef X86
#ifdef _WIN32

#ifdef __BORLANDC__
  /* Win 32 Borland C */
  if (old_cw) {
    unsigned short cw = (unsigned short) *old_cw;
    _control87(cw, 0xFFFF);
  }
#else
  /* Win 32 MSVC */
  if (old_cw) {
    _control87(*old_cw, 0xFFFFFFFF);
  }
#endif

#else
  /* Linux */
  if (old_cw) {
    int cw;
    cw = *old_cw;
    _FPU_SETCW(cw);
  }
#endif
#endif
}

#ifdef HAVE_FORTRAN

#define f_fpu_fix_start FC_FUNC_(f_fpu_fix_start, F_FPU_FIX_START)
#define f_fpu_fix_end   FC_FUNC_(f_fpu_fix_end,   F_FPU_FIX_END)

inline void f_fpu_fix_start(unsigned int *old_cw) {
  fpu_fix_start(old_cw);
}

inline void f_fpu_fix_end(unsigned int *old_cw) {
  fpu_fix_end(old_cw);
}

#endif

}
 

/*
 * src/dd_real.cc
 *
 * This work was supported by the Director, Office of Science, Division
 * of Mathematical, Information, and Computational Sciences of the
 * U.S. Department of Energy under contract number DE-AC03-76SF00098.
 *
 * Copyright (c) 2000-2007
 *
 * Contains implementation of non-inlined functions of double-double
 * package.  Inlined functions are found in dd_inline.h (in include directory).
 */
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cstring>
#include <iostream>
#include <iomanip>
#include <string>



#ifndef QD_INLINE
#endif

using std::cout;
using std::cerr;
using std::endl;
using std::ostream;
using std::istream;
using std::ios_base;
using std::string;
using std::setw;

/* This routine is called whenever a fatal error occurs. */
inline void dd_real::error(const char *msg) { 
  if (msg) { cerr << "ERROR " << msg << endl; }
}

/* Computes the square root of the double-double number dd.
   NOTE: dd must be a non-negative number.                   */
inline dd_real sqrt(const dd_real &a) {
  /* Strategy:  Use Karp's trick:  if x is an approximation
     to sqrt(a), then

        sqrt(a) = a*x + [a - (a*x)^2] * x / 2   (approx)

     The approximation is accurate to twice the accuracy of x.
     Also, the multiplication (a*x) and [-]*x can be done with
     only half the precision.
  */

  if (a.is_zero())
    return 0.0;

  if (a.is_negative()) {
    dd_real::error("(dd_real::sqrt): Negative argument.");
    return dd_real::_nan;
  }

  double x = 1.0 / std::sqrt(a.x[0]);
  double ax = a.x[0] * x;
  return dd_real::add(ax, (a - dd_real::sqr(ax)).x[0] * (x * 0.5));
}

/* Computes the square root of a double in double-double precision. 
   NOTE: d must not be negative.                                   */
inline dd_real dd_real::sqrt(double d) {
  return ::sqrt(dd_real(d));
}

/* Computes the n-th root of the double-double number a.
   NOTE: n must be a positive integer.  
   NOTE: If n is even, then a must not be negative.       */
inline dd_real nroot(const dd_real &a, int n) {
  /* Strategy:  Use Newton iteration for the function

          f(x) = x^(-n) - a

     to find its root a^{-1/n}.  The iteration is thus

          x' = x + x * (1 - a * x^n) / n

     which converges quadratically.  We can then find 
    a^{1/n} by taking the reciprocal.
  */

  if (n <= 0) {
    dd_real::error("(dd_real::nroot): N must be positive.");
    return dd_real::_nan;
  }

  if (n%2 == 0 && a.is_negative()) {
    dd_real::error("(dd_real::nroot): Negative argument.");
    return dd_real::_nan;
  }

  if (n == 1) {
    return a;
  } 
  if (n == 2) {
    return sqrt(a);
  }

  if (a.is_zero())
    return 0.0;

  /* Note  a^{-1/n} = exp(-log(a)/n) */
  dd_real r = abs(a);
  dd_real x = std::exp(-std::log(r.x[0]) / n);

  /* Perform Newton's iteration. */
  x += x * (1.0 - r * npwr(x, n)) / static_cast<double>(n);
  if (a.x[0] < 0.0)
    x = -x;
  return 1.0/x;
}

/* Computes the n-th power of a double-double number. 
   NOTE:  0^0 causes an error.                         */
inline dd_real npwr(const dd_real &a, int n) {
  
  if (n == 0) {
    if (a.is_zero()) {
      dd_real::error("(dd_real::npwr): Invalid argument.");
      return dd_real::_nan;
    }
    return 1.0;
  }

  dd_real r = a;
  dd_real s = 1.0;
  int N = std::abs(n);

  if (N > 1) {
    /* Use binary exponentiation */
    while (N > 0) {
      if (N % 2 == 1) {
        s *= r;
      }
      N /= 2;
      if (N > 0)
        r = sqr(r);
    }
  } else {
    s = r;
  }

  /* Compute the reciprocal if n is negative. */
  if (n < 0)
    return (1.0 / s);
  
  return s;
}

inline dd_real pow(const dd_real &a, int n) {
  return npwr(a, n);
}

inline dd_real pow(const dd_real &a, const dd_real &b) {
  return exp(b * log(a));
}

inline const int n_inv_fact = 15;
inline const double inv_fact[n_inv_fact][2] = {
  { 1.66666666666666657e-01,  9.25185853854297066e-18},
  { 4.16666666666666644e-02,  2.31296463463574266e-18},
  { 8.33333333333333322e-03,  1.15648231731787138e-19},
  { 1.38888888888888894e-03, -5.30054395437357706e-20},
  { 1.98412698412698413e-04,  1.72095582934207053e-22},
  { 2.48015873015873016e-05,  2.15119478667758816e-23},
  { 2.75573192239858925e-06, -1.85839327404647208e-22},
  { 2.75573192239858883e-07,  2.37677146222502973e-23},
  { 2.50521083854417202e-08, -1.44881407093591197e-24},
  { 2.08767569878681002e-09, -1.20734505911325997e-25},
  { 1.60590438368216133e-10,  1.25852945887520981e-26},
  { 1.14707455977297245e-11,  2.06555127528307454e-28},
  { 7.64716373181981641e-13,  7.03872877733453001e-30},
  { 4.77947733238738525e-14,  4.39920548583408126e-31},
  { 2.81145725434552060e-15,  1.65088427308614326e-31}
};

/* Exponential.  Computes exp(x) in double-double precision. */
inline dd_real exp(const dd_real &a) {
  /* Strategy:  We first reduce the size of x by noting that
     
          exp(kr + m * log(2)) = 2^m * exp(r)^k

     where m and k are integers.  By choosing m appropriately
     we can make |kr| <= log(2) / 2 = 0.347.  Then exp(r) is 
     evaluated using the familiar Taylor series.  Reducing the 
     argument substantially speeds up the convergence.       */  

  const double k = 512.0;
  const double inv_k = 1.0 / k;

  if (a.x[0] <= -709.0)
    return 0.0;

  if (a.x[0] >=  709.0)
    return dd_real::_inf;

  if (a.is_zero())
    return 1.0;

  if (a.is_one())
    return dd_real::_e;

  double m = std::floor(a.x[0] / dd_real::_log2.x[0] + 0.5);
  dd_real r = mul_pwr2(a - dd_real::_log2 * m, inv_k);
  dd_real s, t, p;

  p = sqr(r);
  s = r + mul_pwr2(p, 0.5);
  p *= r;
  t = p * dd_real(inv_fact[0][0], inv_fact[0][1]);
  int i = 0;
  do {
    s += t;
    p *= r;
    ++i;
    t = p * dd_real(inv_fact[i][0], inv_fact[i][1]);
  } while (std::abs(to_double(t)) > inv_k * dd_real::_eps && i < 5);

  s += t;

  s = mul_pwr2(s, 2.0) + sqr(s);
  s = mul_pwr2(s, 2.0) + sqr(s);
  s = mul_pwr2(s, 2.0) + sqr(s);
  s = mul_pwr2(s, 2.0) + sqr(s);
  s = mul_pwr2(s, 2.0) + sqr(s);
  s = mul_pwr2(s, 2.0) + sqr(s);
  s = mul_pwr2(s, 2.0) + sqr(s);
  s = mul_pwr2(s, 2.0) + sqr(s);
  s = mul_pwr2(s, 2.0) + sqr(s);
  s += 1.0;

  return ldexp(s, static_cast<int>(m));
}

/* Logarithm.  Computes log(x) in double-double precision.
   This is a natural logarithm (i.e., base e).            */
inline dd_real log(const dd_real &a) {
  /* Strategy.  The Taylor series for log converges much more
     slowly than that of exp, due to the lack of the factorial
     term in the denominator.  Hence this routine instead tries
     to determine the root of the function

         f(x) = exp(x) - a

     using Newton iteration.  The iteration is given by

         x' = x - f(x)/f'(x) 
            = x - (1 - a * exp(-x))
            = x + a * exp(-x) - 1.
           
     Only one iteration is needed, since Newton's iteration
     approximately doubles the number of digits per iteration. */

  if (a.is_one()) {
    return 0.0;
  }

  if (a.x[0] <= 0.0) {
    dd_real::error("(dd_real::log): Non-positive argument.");
    return dd_real::_nan;
  }

  dd_real x = std::log(a.x[0]);   /* Initial approximation */

  x = x + a * exp(-x) - 1.0;
  return x;
}

inline dd_real log10(const dd_real &a) {
  return log(a) / dd_real::_log10;
}

inline const dd_real _pi16 = dd_real(1.963495408493620697e-01,
                                     7.654042494670957545e-18);

/* Table of sin(k * pi/16) and cos(k * pi/16). */
inline const double sin_table [4][2] = {
  {1.950903220161282758e-01, -7.991079068461731263e-18},
  {3.826834323650897818e-01, -1.005077269646158761e-17},
  {5.555702330196021776e-01,  4.709410940561676821e-17},
  {7.071067811865475727e-01, -4.833646656726456726e-17}
};

inline const double cos_table [4][2] = {
  {9.807852804032304306e-01, 1.854693999782500573e-17},
  {9.238795325112867385e-01, 1.764504708433667706e-17},
  {8.314696123025452357e-01, 1.407385698472802389e-18},
  {7.071067811865475727e-01, -4.833646656726456726e-17}
};

/* Computes sin(a) using Taylor series.
   Assumes |a| <= pi/32.                           */
inline dd_real sin_taylor(const dd_real &a) {
  const double thresh = 0.5 * std::abs(to_double(a)) * dd_real::_eps;
  dd_real r, s, t, x;

  if (a.is_zero()) {
    return 0.0;
  }

  int i = 0;
  x = -sqr(a);
  s = a;
  r = a;
  do {
    r *= x;
    t = r * dd_real(inv_fact[i][0], inv_fact[i][1]);
    s += t;
    i += 2;
  } while (i < n_inv_fact && std::abs(to_double(t)) > thresh);

  return s;
}

inline dd_real cos_taylor(const dd_real &a) {
  const double thresh = 0.5 * dd_real::_eps;
  dd_real r, s, t, x;

  if (a.is_zero()) {
    return 1.0;
  }

  x = -sqr(a);
  r = x;
  s = 1.0 + mul_pwr2(r, 0.5);
  int i = 1;
  do {
    r *= x;
    t = r * dd_real(inv_fact[i][0], inv_fact[i][1]);
    s += t;
    i += 2;
  } while (i < n_inv_fact && std::abs(to_double(t)) > thresh);

  return s;
}

inline void sincos_taylor(const dd_real &a, 
                          dd_real &sin_a, dd_real &cos_a) {
  if (a.is_zero()) {
    sin_a = 0.0;
    cos_a = 1.0;
    return;
  }

  sin_a = sin_taylor(a);
  cos_a = sqrt(1.0 - sqr(sin_a));
}


inline dd_real sin(const dd_real &a) {  

  /* Strategy.  To compute sin(x), we choose integers a, b so that

       x = s + a * (pi/2) + b * (pi/16)

     and |s| <= pi/32.  Using the fact that 

       sin(pi/16) = 0.5 * sqrt(2 - sqrt(2 + sqrt(2)))

     we can compute sin(x) from sin(s), cos(s).  This greatly 
     increases the convergence of the sine Taylor series. */

  if (a.is_zero()) {
    return 0.0;
  }

  // approximately reduce modulo 2*pi
  dd_real z = nint(a / dd_real::_2pi);
  dd_real r = a - dd_real::_2pi * z;

  // approximately reduce modulo pi/2 and then modulo pi/16.
  dd_real t;
  double q = std::floor(r.x[0] / dd_real::_pi2.x[0] + 0.5);
  t = r - dd_real::_pi2 * q;
  int j = static_cast<int>(q);
  q = std::floor(t.x[0] / _pi16.x[0] + 0.5);
  t -= _pi16 * q;
  int k = static_cast<int>(q);
  int abs_k = std::abs(k);

  if (j < -2 || j > 2) {
    dd_real::error("(dd_real::sin): Cannot reduce modulo pi/2.");
    return dd_real::_nan;
  }

  if (abs_k > 4) {
    dd_real::error("(dd_real::sin): Cannot reduce modulo pi/16.");
    return dd_real::_nan;
  }

  if (k == 0) {
    switch (j) {
      case 0:
        return sin_taylor(t);
      case 1:
        return cos_taylor(t);
      case -1:
        return -cos_taylor(t);
      default:
        return -sin_taylor(t);
    }
  }

  dd_real u(cos_table[abs_k-1][0], cos_table[abs_k-1][1]);
  dd_real v(sin_table[abs_k-1][0], sin_table[abs_k-1][1]);
  dd_real sin_t, cos_t;
  sincos_taylor(t, sin_t, cos_t);
  if (j == 0) {
    if (k > 0) {
      r = u * sin_t + v * cos_t;
    } else {
      r = u * sin_t - v * cos_t;
    }
  } else if (j == 1) {
    if (k > 0) {
      r = u * cos_t - v * sin_t;
    } else {
      r = u * cos_t + v * sin_t;
    }
  } else if (j == -1) {
    if (k > 0) {
      r = v * sin_t - u * cos_t;
    } else if (k < 0) {
      r = -u * cos_t - v * sin_t;
    }
  } else {
    if (k > 0) {
      r = -u * sin_t - v * cos_t;
    } else {
      r = v * cos_t - u * sin_t;
    }
  }

  return r;
}

inline dd_real cos(const dd_real &a) {

  if (a.is_zero()) {
    return 1.0;
  }

  // approximately reduce modulo 2*pi
  dd_real z = nint(a / dd_real::_2pi);
  dd_real r = a - z * dd_real::_2pi;

  // approximately reduce modulo pi/2 and then modulo pi/16
  dd_real t;
  double q = std::floor(r.x[0] / dd_real::_pi2.x[0] + 0.5);
  t = r - dd_real::_pi2 * q;
  int j = static_cast<int>(q);
  q = std::floor(t.x[0] / _pi16.x[0] + 0.5);
  t -= _pi16 * q;
  int k = static_cast<int>(q);
  int abs_k = std::abs(k);

  if (j < -2 || j > 2) {
    dd_real::error("(dd_real::cos): Cannot reduce modulo pi/2.");
    return dd_real::_nan;
  }

  if (abs_k > 4) {
    dd_real::error("(dd_real::cos): Cannot reduce modulo pi/16.");
    return dd_real::_nan;
  }

  if (k == 0) {
    switch (j) {
      case 0:
        return cos_taylor(t);
      case 1:
        return -sin_taylor(t);
      case -1:
        return sin_taylor(t);
      default:
        return -cos_taylor(t);
    }
  }

  dd_real sin_t, cos_t;
  sincos_taylor(t, sin_t, cos_t);
  dd_real u(cos_table[abs_k-1][0], cos_table[abs_k-1][1]);
  dd_real v(sin_table[abs_k-1][0], sin_table[abs_k-1][1]);

  if (j == 0) {
    if (k > 0) {
      r = u * cos_t - v * sin_t;
    } else {
      r = u * cos_t + v * sin_t;
    }
  } else if (j == 1) {
    if (k > 0) {
      r = - u * sin_t - v * cos_t;
    } else {
      r = v * cos_t - u * sin_t;
    }
  } else if (j == -1) {
    if (k > 0) {
      r = u * sin_t + v * cos_t;
    } else {
      r = u * sin_t - v * cos_t;
    }
  } else {
    if (k > 0) {
      r = v * sin_t - u * cos_t;
    } else {
      r = - u * cos_t - v * sin_t;
    }
  }

  return r;
}

inline void sincos(const dd_real &a, dd_real &sin_a, dd_real &cos_a) {

  if (a.is_zero()) {
    sin_a = 0.0;
    cos_a = 1.0;
    return;
  }

  // approximately reduce modulo 2*pi
  dd_real z = nint(a / dd_real::_2pi);
  dd_real r = a - dd_real::_2pi * z;

  // approximately reduce module pi/2 and pi/16
  dd_real t;
  double q = std::floor(r.x[0] / dd_real::_pi2.x[0] + 0.5);
  t = r - dd_real::_pi2 * q;
  int j = static_cast<int>(q);
  int abs_j = std::abs(j);
  q = std::floor(t.x[0] / _pi16.x[0] + 0.5);
  t -= _pi16 * q;
  int k = static_cast<int>(q);
  int abs_k = std::abs(k);

  if (abs_j > 2) {
    dd_real::error("(dd_real::sincos): Cannot reduce modulo pi/2.");
    cos_a = sin_a = dd_real::_nan;
    return;
  }

  if (abs_k > 4) {
    dd_real::error("(dd_real::sincos): Cannot reduce modulo pi/16.");
    cos_a = sin_a = dd_real::_nan;
    return;
  }

  dd_real sin_t, cos_t;
  dd_real s, c;

  sincos_taylor(t, sin_t, cos_t);

  if (abs_k == 0) {
    s = sin_t;
    c = cos_t;
  } else {
    dd_real u(cos_table[abs_k-1][0], cos_table[abs_k-1][1]);
    dd_real v(sin_table[abs_k-1][0], sin_table[abs_k-1][1]);

    if (k > 0) {
      s = u * sin_t + v * cos_t;
      c = u * cos_t - v * sin_t;
    } else {
      s = u * sin_t - v * cos_t;
      c = u * cos_t + v * sin_t;
    }
  }

  if (abs_j == 0) {
    sin_a = s;
    cos_a = c;
  } else if (j == 1) {
    sin_a = c;
    cos_a = -s;
  } else if (j == -1) {
    sin_a = -c;
    cos_a = s;
  } else {
    sin_a = -s;
    cos_a = -c;
  }
  
}

inline dd_real atan(const dd_real &a) {
  return atan2(a, dd_real(1.0));
}

inline dd_real atan2(const dd_real &y, const dd_real &x) {
  /* Strategy: Instead of using Taylor series to compute 
     arctan, we instead use Newton's iteration to solve
     the equation

        sin(z) = y/r    or    cos(z) = x/r

     where r = sqrt(x^2 + y^2).
     The iteration is given by

        z' = z + (y - sin(z)) / cos(z)          (for equation 1)
        z' = z - (x - cos(z)) / sin(z)          (for equation 2)

     Here, x and y are normalized so that x^2 + y^2 = 1.
     If |x| > |y|, then first iteration is used since the 
     denominator is larger.  Otherwise, the second is used.
  */

  if (x.is_zero()) {
    
    if (y.is_zero()) {
      /* Both x and y is zero. */
      dd_real::error("(dd_real::atan2): Both arguments zero.");
      return dd_real::_nan;
    }

    return (y.is_positive()) ? dd_real::_pi2 : -dd_real::_pi2;
  } else if (y.is_zero()) {
    return (x.is_positive()) ? dd_real(0.0) : dd_real::_pi;
  }

  if (x == y) {
    return (y.is_positive()) ? dd_real::_pi4 : -dd_real::_3pi4;
  }

  if (x == -y) {
    return (y.is_positive()) ? dd_real::_3pi4 : -dd_real::_pi4;
  }

  dd_real r = sqrt(sqr(x) + sqr(y));
  dd_real xx = x / r;
  dd_real yy = y / r;

  /* Compute double precision approximation to atan. */
  dd_real z = std::atan2(to_double(y), to_double(x));
  dd_real sin_z, cos_z;

  if (std::abs(xx.x[0]) > std::abs(yy.x[0])) {
    /* Use Newton iteration 1.  z' = z + (y - sin(z)) / cos(z)  */
    sincos(z, sin_z, cos_z);
    z += (yy - sin_z) / cos_z;
  } else {
    /* Use Newton iteration 2.  z' = z - (x - cos(z)) / sin(z)  */
    sincos(z, sin_z, cos_z);
    z -= (xx - cos_z) / sin_z;
  }

  return z;
}

inline dd_real tan(const dd_real &a) {
  dd_real s, c;
  sincos(a, s, c);
  return s/c;
}

inline dd_real asin(const dd_real &a) {
  dd_real abs_a = abs(a);

  if (abs_a > 1.0) {
    dd_real::error("(dd_real::asin): Argument out of domain.");
    return dd_real::_nan;
  }

  if (abs_a.is_one()) {
    return (a.is_positive()) ? dd_real::_pi2 : -dd_real::_pi2;
  }

  return atan2(a, sqrt(1.0 - sqr(a)));
}

inline dd_real acos(const dd_real &a) {
  dd_real abs_a = abs(a);

  if (abs_a > 1.0) {
    dd_real::error("(dd_real::acos): Argument out of domain.");
    return dd_real::_nan;
  }

  if (abs_a.is_one()) {
    return (a.is_positive()) ? dd_real(0.0) : dd_real::_pi;
  }

  return atan2(sqrt(1.0 - sqr(a)), a);
}
 
inline dd_real sinh(const dd_real &a) {
  if (a.is_zero()) {
    return 0.0;
  }

  if (abs(a) > 0.05) {
    dd_real ea = exp(a);
    return mul_pwr2(ea - inv(ea), 0.5);
  }

  /* since a is small, using the above formula gives
     a lot of cancellation.  So use Taylor series.   */
  dd_real s = a;
  dd_real t = a;
  dd_real r = sqr(t);
  double m = 1.0;
  double thresh = std::abs((to_double(a)) * dd_real::_eps);

  do {
    m += 2.0;
    t *= r;
    t /= (m-1) * m;

    s += t;
  } while (abs(t) > thresh);

  return s;

}

inline dd_real cosh(const dd_real &a) {
  if (a.is_zero()) {
    return 1.0;
  }

  dd_real ea = exp(a);
  return mul_pwr2(ea + inv(ea), 0.5);
}

inline dd_real tanh(const dd_real &a) {
  if (a.is_zero()) {
    return 0.0;
  }

  if (std::abs(to_double(a)) > 0.05) {
    dd_real ea = exp(a);
    dd_real inv_ea = inv(ea);
    return (ea - inv_ea) / (ea + inv_ea);
  } else {
    dd_real s, c;
    s = sinh(a);
    c = sqrt(1.0 + sqr(s));
    return s / c;
  }
}

inline void sincosh(const dd_real &a, dd_real &s, dd_real &c) {
  if (std::abs(to_double(a)) <= 0.05) {
    s = sinh(a);
    c = sqrt(1.0 + sqr(s));
  } else {
    dd_real ea = exp(a);
    dd_real inv_ea = inv(ea);
    s = mul_pwr2(ea - inv_ea, 0.5);
    c = mul_pwr2(ea + inv_ea, 0.5);
  }
}

inline dd_real asinh(const dd_real &a) {
  return log(a + sqrt(sqr(a) + 1.0));
}

inline dd_real acosh(const dd_real &a) {
  if (a < 1.0) {
    dd_real::error("(dd_real::acosh): Argument out of domain.");
    return dd_real::_nan;
  }

  return log(a + sqrt(sqr(a) - 1.0));
}

inline dd_real atanh(const dd_real &a) {
  if (abs(a) >= 1.0) {
    dd_real::error("(dd_real::atanh): Argument out of domain.");
    return dd_real::_nan;
  }

  return mul_pwr2(log((1.0 + a) / (1.0 - a)), 0.5);
}

inline dd_real fmod(const dd_real &a, const dd_real &b) {
  dd_real n = aint(a / b);
  return (a - b * n);
}

inline dd_real ddrand() {
  static const double m_const = 4.6566128730773926e-10;  /* = 2^{-31} */
  double m = m_const;
  dd_real r = 0.0;
  double d;

  /* Strategy:  Generate 31 bits at a time, using lrand48 
     random number generator.  Shift the bits, and reapeat
     4 times. */

  for (int i = 0; i < 4; i++, m *= m_const) {
//    d = lrand48() * m;
    d = std::rand() * m;
    r += d;
  }

  return r;
}

/* polyeval(c, n, x)
   Evaluates the given n-th degree polynomial at x.
   The polynomial is given by the array of (n+1) coefficients. */
inline dd_real polyeval(const dd_real *c, int n, const dd_real &x) {
  /* Just use Horner's method of polynomial evaluation. */
  dd_real r = c[n];
  
  for (int i = n-1; i >= 0; i--) {
    r *= x;
    r += c[i];
  }

  return r;
}

/* polyroot(c, n, x0)
   Given an n-th degree polynomial, finds a root close to 
   the given guess x0.  Note that this uses simple Newton
   iteration scheme, and does not work for multiple roots.  */
inline dd_real polyroot(const dd_real *c, int n, 
    const dd_real &x0, int max_iter, double thresh) {
  dd_real x = x0;
  dd_real f;
  dd_real *d = new dd_real[n];
  bool conv = false;
  int i;
  double max_c = std::abs(to_double(c[0]));
  double v;

  if (thresh == 0.0) thresh = dd_real::_eps;

  /* Compute the coefficients of the derivatives. */
  for (i = 1; i <= n; i++) {
    v = std::abs(to_double(c[i]));
    if (v > max_c) max_c = v;
    d[i-1] = c[i] * static_cast<double>(i);
  }
  thresh *= max_c;

  /* Newton iteration. */
  for (i = 0; i < max_iter; i++) {
    f = polyeval(c, n, x);

    if (abs(f) < thresh) {
      conv = true;
      break;
    }
    x -= (f / polyeval(d, n-1, x));
  }
  delete [] d;

  if (!conv) {
    dd_real::error("(dd_real::polyroot): Failed to converge.");
    return dd_real::_nan;
  }

  return x;
}


/* Constructor.  Reads a double-double number from the string s
   and constructs a double-double number.                         */
inline dd_real::dd_real(const char *s) {
  if (dd_real::read(s, *this)) {
    dd_real::error("(dd_real::dd_real): INPUT ERROR.");
    *this = dd_real::_nan;
  }
}

inline dd_real &dd_real::operator=(const char *s) {
  if (dd_real::read(s, *this)) {
    dd_real::error("(dd_real::operator=): INPUT ERROR.");
    *this = dd_real::_nan;
  }
  return *this;
}

/* Outputs the double-double number dd. */
inline ostream &operator<<(ostream &os, const dd_real &dd) {
  bool showpos = (os.flags() & ios_base::showpos) != 0;
  bool uppercase =  (os.flags() & ios_base::uppercase) != 0;
  return os << dd.to_string(os.precision(), os.width(), os.flags(), 
      showpos, uppercase, os.fill());
}

/* Reads in the double-double number a. */
inline istream &operator>>(istream &s, dd_real &a) {
  char str[255];
  s >> str;
  a = dd_real(str);
  return s;
}

inline void dd_real::to_digits(char *s, int &expn, int precision) const {
  int D = precision + 1;  /* number of digits to compute */

  dd_real r = abs(*this);
  int e;  /* exponent */
  int i, d;

  if (x[0] == 0.0) {
    /* this == 0.0 */
    expn = 0;
    for (i = 0; i < precision; i++) s[i] = '0';
    return;
  }

  /* First determine the (approximate) exponent. */
  e = to_int(std::floor(std::log10(std::abs(x[0]))));

  if (e < -300) {
    r *= dd_real(10.0) ^ 300;
    r /= dd_real(10.0) ^ (e + 300);
  } else if (e > 300) {
    r = ldexp(r, -53);
    r /= dd_real(10.0) ^ e;
    r = ldexp(r, 53);
  } else {
    r /= dd_real(10.0) ^ e;
  }

  /* Fix exponent if we are off by one */
  if (r >= 10.0) {
    r /= 10.0;
    e++;
  } else if (r < 1.0) {
    r *= 10.0;
    e--;
  }

  if (r >= 10.0 || r < 1.0) {
    dd_real::error("(dd_real::to_digits): can't compute exponent.");
    return;
  }

  /* Extract the digits */
  for (i = 0; i < D; i++) {
    d = static_cast<int>(r.x[0]);
    r -= d;
    r *= 10.0;

    s[i] = static_cast<char>(d + '0');
  }

  /* Fix out of range digits. */
  for (i = D-1; i > 0; i--) {
    if (s[i] < '0') {
      s[i-1]--;
      s[i] += 10;
    } else if (s[i] > '9') {
      s[i-1]++;
      s[i] -= 10;
    }
  }

  if (s[0] <= '0') {
    dd_real::error("(dd_real::to_digits): non-positive leading digit.");
    return;
  }

  /* Round, handle carry */
  if (s[D-1] >= '5') {
    s[D-2]++;

    i = D-2;
    while (i > 0 && s[i] > '9') {
      s[i] -= 10;
      s[--i]++;
    }
  }

  /* If first digit is 10, shift everything. */
  if (s[0] > '9') { 
    e++; 
    for (i = precision; i >= 2; i--) s[i] = s[i-1]; 
    s[0] = '1';
    s[1] = '0';
  }

  s[precision] = 0;
  expn = e;
}

/* Writes the double-double number into the character array s of length len.
   The integer d specifies how many significant digits to write.
   The string s must be able to hold at least (d+8) characters.  
   showpos indicates whether to use the + sign, and uppercase indicates
   whether the E or e is to be used for the exponent. */
inline void dd_real::write(char *s, int len, int precision, 
    bool showpos, bool uppercase) const {
  string str = to_string(precision, 0, ios_base::scientific, showpos, uppercase);
  std::strncpy(s, str.c_str(), len-1);
  s[len-1] = 0;
}


inline void round_string(char *s, int precision, int *offset){
	/*
	 Input string must be all digits or errors will occur.
	 */

	int i;
	int D = precision ;

	/* Round, handle carry */
	  if (D>0 && s[D] >= '5') {
	    s[D-1]++;

	    i = D-1;
	    while (i > 0 && s[i] > '9') {
	      s[i] -= 10;
	      s[--i]++;
	    }
	  }

	  /* If first digit is 10, shift everything. */
	  if (s[0] > '9') {
	    // e++; // don't modify exponent here
	    for (i = precision; i >= 1; i--) s[i+1] = s[i];
	    s[0] = '1';
	    s[1] = '0';

	    (*offset)++ ; // now offset needs to be increased by one
	    precision++ ;
	  }

	  s[precision] = 0; // add terminator for array
}

inline string dd_real::to_string(int precision, int width, ios_base::fmtflags fmt, 
    bool showpos, bool uppercase, char fill) const {
  string s;
  bool fixed = (fmt & ios_base::fixed) != 0;
  bool sgn = true;
  int i, e = 0;

  if (isnan()) {
    s = uppercase ? "NAN" : "nan";
    sgn = false;
  } else {
    if (*this < 0.0)
      s += '-';
    else if (showpos)
      s += '+';
    else
      sgn = false;

    if (isinf()) {
      s += uppercase ? "INF" : "inf";
    } else if (*this == 0.0) {
      /* Zero case */
      s += '0';
      if (precision > 0) {
        s += '.';
        s.append(precision, '0');
      }
    } else {
      /* Non-zero case */
      int off = (fixed ? (1 + to_int(floor(log10(abs(*this))))) : 1);
      int d = precision + off;

      int d_with_extra = d;
      if(fixed)
    	  d_with_extra = std::max(60, d); // longer than the max accuracy for DD

      // highly special case - fixed mode, precision is zero, abs(*this) < 1.0
      // without this trap a number like 0.9 printed fixed with 0 precision prints as 0
      // should be rounded to 1.
      if(fixed && (precision == 0) && (abs(*this) < 1.0)){
    	  if(abs(*this) >= 0.5)
    		  s += '1';
    	  else
    		  s += '0';

    	  return s;
      }

      // handle near zero to working precision (but not exactly zero)
      if (fixed && d <= 0) {
        s += '0';
        if (precision > 0) {
          s += '.';
          s.append(precision, '0');
        }
      } else { // default

        char *t; //  = new char[d+1];
        int j;

        if(fixed){
        	t = new char[d_with_extra+1];
        	to_digits(t, e, d_with_extra);
        }
        else{
        	t = new char[d+1];
        	to_digits(t, e, d);
        }

        off = e + 1;

        if (fixed) {
          // fix the string if it's been computed incorrectly
          // round here in the decimal string if required
          round_string(t, d, &off);

          if (off > 0) {
            for (i = 0; i < off; i++) s += t[i];
            if (precision > 0) {
              s += '.';
              for (j = 0; j < precision; j++, i++) s += t[i];
            }
          } else {
            s += "0.";
            if (off < 0) s.append(-off, '0');
            for (i = 0; i < d; i++) s += t[i];
          }
        } else {
          s += t[0];
          if (precision > 0) s += '.';

          for (i = 1; i <= precision; i++)
            s += t[i];

        }
		delete [] t;
      }
    }

    // trap for improper offset with large values
    // without this trap, output of values of the for 10^j - 1 fail for j > 28
    // and are output with the point in the wrong place, leading to a dramatically off value
    if(fixed && (precision > 0)){
    	// make sure that the value isn't dramatically larger
    	double from_string = atof(s.c_str());

    	// if this ratio is large, then we've got problems
    	if( fabs( from_string / this->x[0] ) > 3.0 ){

    		int point_position;
    		char temp;

    		// loop on the string, find the point, move it up one
    		// don't act on the first character
    		for(i=1; i < s.length(); i++){
    			if(s[i] == '.'){
    				s[i] = s[i-1] ;
    				s[i-1] = '.' ;
    				break;
    			}
    		}

        	from_string = atof(s.c_str());
        	// if this ratio is large, then the string has not been fixed
        	if( fabs( from_string / this->x[0] ) > 3.0 ){
        		dd_real::error("Re-rounding unsuccessful in large number fixed point trap.") ;
        	}
    	}
    }


    if (!fixed && !isinf()) {
      /* Fill in exponent part */
      s += uppercase ? 'E' : 'e';
      append_expn(s, e);
    }
  }

  /* Fill in the blanks */
  int len = s.length();
  if (len < width) {
    int delta = width - len;
    if (fmt & ios_base::internal) {
      if (sgn)
        s.insert(static_cast<string::size_type>(1), delta, fill);
      else
        s.insert(static_cast<string::size_type>(0), delta, fill);
    } else if (fmt & ios_base::left) {
      s.append(delta, fill);
    } else {
      s.insert(static_cast<string::size_type>(0), delta, fill);
    }
  }

  return s;
}

/* Reads in a double-double number from the string s. */
inline int dd_real::read(const char *s, dd_real &a) {
  const char *p = s;
  char ch;
  int sign = 0;
  int point = -1;
  int nd = 0;
  int e = 0;
  bool done = false;
  dd_real r = 0.0;
  int nread;
  
  /* Skip any leading spaces */
  while (*p == ' ')
    p++;

  while (!done && (ch = *p) != '\0') {
    if (ch >= '0' && ch <= '9') {
      int d = ch - '0';
      r *= 10.0;
      r += static_cast<double>(d);
      nd++;
    } else {

      switch (ch) {

      case '.':
        if (point >= 0)
          return -1;        
        point = nd;
        break;

      case '-':
      case '+':
        if (sign != 0 || nd > 0)
          return -1;
        sign = (ch == '-') ? -1 : 1;
        break;

      case 'E':
      case 'e':
        nread = std::sscanf(p+1, "%d", &e);
        done = true;
        if (nread != 1)
          return -1;
        break;

      default:
        return -1;
      }
    }
    
    p++;
  }

  if (point >= 0) {
    e -= (nd - point);
  }

  if (e != 0) {
    r *= (dd_real(10.0) ^ e);
  }

  a = (sign == -1) ? -r : r;
  return 0;
}

/* Debugging routines */
inline void dd_real::dump(const string &name, std::ostream &os) const {
  std::ios_base::fmtflags old_flags = os.flags();
  std::streamsize old_prec = os.precision(19);
  os << std::scientific;

  if (name.length() > 0) os << name << " = ";
  os << "[ " << setw(27) << x[0] << ", " << setw(27) << x[1] << " ]" << endl;

  os.precision(old_prec);
  os.flags(old_flags);
}

inline void dd_real::dump_bits(const string &name, std::ostream &os) const {
  string::size_type len = name.length();
  if (len > 0) {
    os << name << " = ";
    len +=3;
  }
  os << "[ ";
  len += 2;
  print_double_info(os, x[0]);
  os << endl;
  for (string::size_type i = 0; i < len; i++) os << ' ';
  print_double_info(os, x[1]);
  os << " ]" << endl;
}

inline dd_real dd_real::debug_rand() { 

  if (std::rand() % 2 == 0)
    return ddrand();

  int expn = 0;
  dd_real a = 0.0;
  double d;
  for (int i = 0; i < 2; i++) {
    d = std::ldexp(static_cast<double>(std::rand()) / RAND_MAX, -expn);
    a += d;
    expn = expn + 54 + std::rand() % 200;
  }
  return a;
}
namespace tsp_native_math {
using Base = dd_real;
struct Real {
    Base value;
    Real() = default;
    Real(double v) : value(v) {}
    Real(int v) : value(double(v)) {}
    Real(const char *v) : value(v) {}
    explicit Real(const Base &v) : value(v) {}
    explicit operator double() const { return ::to_double(value); }
    Real &operator+=(const Real &b) {value+=b.value;return *this;}
    Real &operator-=(const Real &b) {value-=b.value;return *this;}
    Real &operator*=(const Real &b) {value*=b.value;return *this;}
    Real &operator/=(const Real &b) {value/=b.value;return *this;}
    bool isfinite() const {return value.isfinite();}
};
#define TA_NATIVE_BINARY(symbol) \
inline Real operator symbol(const Real &a,const Real &b) {return Real(a.value symbol b.value);} \
template<class T, typename std::enable_if<std::is_arithmetic<T>::value,int>::type =0> inline Real operator symbol(const Real &a,T b) {return Real(a.value symbol double(b));} \
template<class T, typename std::enable_if<std::is_arithmetic<T>::value,int>::type =0> inline Real operator symbol(T a,const Real &b) {return Real(double(a) symbol b.value);}
TA_NATIVE_BINARY(+)
TA_NATIVE_BINARY(-)
TA_NATIVE_BINARY(*)
TA_NATIVE_BINARY(/)
#undef TA_NATIVE_BINARY
#define TA_NATIVE_COMPARE(symbol) \
inline bool operator symbol(const Real &a,const Real &b) {return a.value symbol b.value;} \
template<class T, typename std::enable_if<std::is_arithmetic<T>::value,int>::type =0> inline bool operator symbol(const Real &a,T b) {return a.value symbol double(b);} \
template<class T, typename std::enable_if<std::is_arithmetic<T>::value,int>::type =0> inline bool operator symbol(T a,const Real &b) {return double(a) symbol b.value;}
TA_NATIVE_COMPARE(==)
TA_NATIVE_COMPARE(!=)
TA_NATIVE_COMPARE(<)
TA_NATIVE_COMPARE(>)
TA_NATIVE_COMPARE(<=)
TA_NATIVE_COMPARE(>=)
#undef TA_NATIVE_COMPARE
inline Real operator-(const Real &x) {return Real(-x.value);}
#define TA_NATIVE_UNARY(name) inline Real name(const Real &x) {return Real(::name(x.value));}
TA_NATIVE_UNARY(abs)
TA_NATIVE_UNARY(sqrt)
TA_NATIVE_UNARY(sin)
TA_NATIVE_UNARY(cos)
TA_NATIVE_UNARY(log)
TA_NATIVE_UNARY(atan)
TA_NATIVE_UNARY(atanh)
#undef TA_NATIVE_UNARY
inline Real atan2(const Real &y, const Real &x) {return Real(::atan2(y.value,x.value));}
inline Real pow(const Real &x, int n) {return Real(::pow(x.value,n));}
inline Real ldexp(const Real &x, int n) {return Real(::ldexp(x.value,n));}
inline bool isfinite(const Real &x) {return x.isfinite();}
inline Real log1p(const Real &x) {
    if (abs(x) > Real(".001")) return log(Real(1) + x);
    Real sum = -Real(1)/32;
    for (int k=31; k>=1; --k) sum = Real(k%2?1:-1)/k + x*sum;
    return x*sum;
}
inline Real asinh(const Real &x) {
    const bool negative = x < 0;
    const Real a = negative ? -x : x;
    const Real out = a > Real("1e100") ? log(a) + Real("0.6931471805599453094172321214581765680755001343602552541")
        : log1p(a + a*a/(1 + sqrt(1 + a*a)));
    return negative ? -out : out;
}
}

namespace tsp_math {
using Real = tsp_native_math::Real;
namespace working_math = tsp_native_math;
inline const Real PI("3.1415926535897932384626433832795028841971693993751");
struct Exp {
    int x, y, z;
};
inline constexpr std::array<Exp, 20> ex{{{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {2, 0, 0},
                                         {1, 1, 0}, {1, 0, 1}, {0, 2, 0}, {0, 1, 1}, {0, 0, 2},
                                         {3, 0, 0}, {2, 1, 0}, {2, 0, 1}, {1, 2, 0}, {1, 1, 1},
                                         {1, 0, 2}, {0, 3, 0}, {0, 2, 1}, {0, 1, 2}, {0, 0, 3}}};
inline constexpr std::uint32_t JET_FULL_MASK = (std::uint32_t(1) << 20) - 1;

inline thread_local std::uint32_t active_jet_mask = JET_FULL_MASK;
inline bool jet_active(int i) {
    return i >= 0 && i < 20 && ((active_jet_mask >> i) & std::uint32_t(1)) != 0;
}
struct JetScope {
    std::uint32_t previous;
    explicit JetScope(std::uint32_t mask) : previous(active_jet_mask) { active_jet_mask = mask; }
    ~JetScope() { active_jet_mask = previous; }
    JetScope(const JetScope &) = delete;
    JetScope &operator=(const JetScope &) = delete;
};

struct Prod {
    int i, j, k;
};
inline const auto products = [] {
    std::vector<Prod> v;
    for (int i = 0; i < 20; ++i)
        for (int j = 0; j < 20; ++j)
            for (int k = 0; k < 20; ++k)
                if (ex[i].x + ex[j].x == ex[k].x && ex[i].y + ex[j].y == ex[k].y &&
                    ex[i].z + ex[j].z == ex[k].z)
                    v.push_back({i, j, k});
    return v;
}();
inline const auto centered_square_products = [] {
    std::vector<Prod> v;
    for (auto t : products) if (t.i && t.j) v.push_back(t);
    return v;
}();
inline const auto centered_cube_products = [] {
    std::vector<Prod> v;
    for (auto t : products)
        if (ex[t.i].x + ex[t.i].y + ex[t.i].z == 2 &&
            ex[t.j].x + ex[t.j].y + ex[t.j].z == 1) v.push_back(t);
    return v;
}();

inline const auto inverse_products = [] {
    std::array<std::vector<Prod>,20> out;
    for(auto t:products) if(t.i) out[t.k].push_back(t);
    return out;
}();
struct J {
    std::array<Real, 20> a{};
    J() = default;
    J(double v) { a[0] = v; }
    J(const Real &v) { a[0] = v; }
    double value() const { return static_cast<double>(a[0]); }
    static J variable(const Real &v, int axis) {
        J r(v);
        if (jet_active(axis + 1)) r.a[axis + 1] = 1;
        return r;
    }
};

inline constexpr double CARLSON_TERMINAL_TOL = 1.0e-6;

inline double carlson_jet_tolerance(const J &x) {
    const Real large = Real(1000000) * (1 + ::tsp_math::working_math::abs(x.a[0]));
    for (int i = 1; i < 20; ++i) if (jet_active(i))
        if (::tsp_math::working_math::abs(x.a[i]) > large)
            return 1.0e-8;
    return CARLSON_TERMINAL_TOL;
}
inline J operator+(J a, const J &b) {
    for (int i = 0; i < 20; ++i) if (jet_active(i)) a.a[i] += b.a[i];
    return a;
}
inline J operator-(J a, const J &b) {
    for (int i = 0; i < 20; ++i) if (jet_active(i)) a.a[i] -= b.a[i];
    return a;
}
inline J operator-(J a) {
    for (int i = 0; i < 20; ++i) if (jet_active(i)) a.a[i] = -a.a[i];
    return a;
}
inline J operator*(const J &a, const J &b) {
    J c;
    bool constant_a = true, constant_b = true;
    for (int i = 1; i < 20; ++i) if (jet_active(i))
        if (a.a[i] != 0) { constant_a = false; break; }
    if (constant_a) {
        for (int i = 0; i < 20; ++i) if (jet_active(i)) c.a[i] = a.a[0] * b.a[i];
        return c;
    }
    for (int i = 1; i < 20; ++i) if (jet_active(i))
        if (b.a[i] != 0) { constant_b = false; break; }
    if (constant_b) {
        for (int i = 0; i < 20; ++i) if (jet_active(i)) c.a[i] = a.a[i] * b.a[0];
        return c;
    }
    for (auto t : products) if (jet_active(t.k) && jet_active(t.i) && jet_active(t.j))
        c.a[t.k] += a.a[t.i] * b.a[t.j];
    return c;
}
inline J operator*(J a, double b) {

    int exponent = 0;
    if (std::isfinite(b) && std::frexp(std::abs(b), &exponent) == 0.5) {
        for (int i = 0; i < 20; ++i) if (jet_active(i)) {
            a.a[i] = ::tsp_math::working_math::ldexp(a.a[i], exponent - 1);
            if (b < 0) a.a[i] = -a.a[i];
        }
        return a;
    }
    for (int i = 0; i < 20; ++i) if (jet_active(i)) a.a[i] *= b;
    return a;
}
inline J operator*(double b, J a) { return a * b; }
inline J operator/(J a, double b) {
    int exponent = 0;
    if (std::isfinite(b) && std::frexp(std::abs(b), &exponent) == 0.5) {
        for (int i = 0; i < 20; ++i) if (jet_active(i)) {
            a.a[i] = ::tsp_math::working_math::ldexp(a.a[i], 1 - exponent);
            if (b < 0) a.a[i] = -a.a[i];
        }
        return a;
    }
    for (int i = 0; i < 20; ++i) if (jet_active(i)) a.a[i] /= b;
    return a;
}
inline J operator+(J a, double b) {
    a.a[0] += b;
    return a;
}
inline J operator+(double b, J a) { return a + b; }
inline J operator-(J a, double b) {
    a.a[0] -= b;
    return a;
}
inline J operator-(double b, J a) { return -a + b; }
inline J operator*(J a, const Real &b) {
    for (int i = 0; i < 20; ++i) if (jet_active(i)) a.a[i] *= b;
    return a;
}
inline J operator*(const Real &b, J a) { return a * b; }
inline J compose(J x, Real f, Real d1, Real d2, Real d3) {
    x.a[0] = 0;

    J square, cube, out(f);
    for (auto t : centered_square_products)
        if (jet_active(t.k) && jet_active(t.i) && jet_active(t.j)) square.a[t.k] += x.a[t.i] * x.a[t.j];
    for (auto t : centered_cube_products)
        if (jet_active(t.k) && jet_active(t.i) && jet_active(t.j)) cube.a[t.k] += square.a[t.i] * x.a[t.j];
    for (int i = 1; i < 4; ++i) if (jet_active(i)) out.a[i] = x.a[i] * d1;
    for (int i = 4; i < 10; ++i) if (jet_active(i)) out.a[i] = x.a[i] * d1 + square.a[i] * d2;
    for (int i = 10; i < 20; ++i) if (jet_active(i)) out.a[i] = x.a[i] * d1 + square.a[i] * d2 + cube.a[i] * d3;
    return out;
}
inline J inverse_uncached(const J &x) {

    const Real v=x.a[0], magnitude=::tsp_math::working_math::abs(v);
    bool conditioned=magnitude>=Real("1e-8");
    if(conditioned) {
        const Real limit=10000*magnitude;
        for(int k=1;k<20;++k) if (jet_active(k))
            if(::tsp_math::working_math::abs(x.a[k])>limit) {conditioned=false;break;}
    }
    if(!conditioned) {
        const Real v2=v*v,v3=v2*v,v4=v3*v;
        return compose(x,1/v,-1/v2,1/v3,-1/v4);
    }
    J out(Real(1/x.a[0]));
    for(int k=1;k<20;++k) if (jet_active(k)) {
        Real sum=0;
        for(auto t:inverse_products[k]) sum+=x.a[t.i]*out.a[t.j];
        out.a[k]=-sum*out.a[0];
    }
    return out;
}
struct InverseJetEntry { J argument, value; };
struct InverseJetCache {
    std::array<InverseJetEntry, 32> entries{};
    std::size_t used = 0, replacement = 0;
};
inline thread_local InverseJetCache *active_inverse_cache = nullptr;
struct InverseJetScope {
    InverseJetCache cache;
    InverseJetCache *previous = active_inverse_cache;
    InverseJetScope() { active_inverse_cache = &cache; }
    ~InverseJetScope() { active_inverse_cache = previous; }
    InverseJetScope(const InverseJetScope &) = delete;
    InverseJetScope &operator=(const InverseJetScope &) = delete;
};
inline J inverse(const J &x) {
    if (!active_inverse_cache) return inverse_uncached(x);
    auto &cache = *active_inverse_cache;
    for (std::size_t k = 0; k < cache.used; ++k) {
        if (cache.entries[k].argument.a == x.a) {
            return cache.entries[k].value;
        }
    }

    const J out = inverse_uncached(x);
    const std::size_t k = cache.used < cache.entries.size()
        ? cache.used++ : cache.replacement++ % cache.entries.size();
    cache.entries[k].argument = x;
    cache.entries[k].value = out;
    return out;
}
inline J operator/(const J &a, const J &b) {
    bool constant = true;
    for (int i = 1; i < 20; ++i) if (jet_active(i))
        if (b.a[i] != 0) { constant = false; break; }
    if (constant) return a * Real(1 / b.a[0]);
    return a * inverse(b);
}
inline J &operator+=(J &a, const J &b) { return a = a + b; }
inline J &operator-=(J &a, const J &b) { return a = a - b; }
inline J &operator*=(J &a, const J &b) { return a = a * b; }
inline bool operator<(const J &a, const J &b) { return a.value() < b.value(); }
inline bool operator>(const J &a, const J &b) { return b < a; }
inline bool operator<=(const J &a, const J &b) { return !(a > b); }
inline bool operator>=(const J &a, const J &b) { return !(a < b); }
inline bool operator==(const J &a, const J &b) { return a.value() == b.value(); }
inline J sqrt(const J &x) {
    Real v = x.a[0], s = ::tsp_math::working_math::sqrt(v);
    if (v == 0) {
        bool zero = true;
        for (auto q : x.a)
            zero = zero && (q == 0);
        if (zero)
            return J(0);
    }
    if (!(v > 0))
        throw std::domain_error("nonpositive jet square root");
    J out(s);
    for(int k=1;k<20;++k)if(jet_active(k)) {
        Real sum=0;
        for(const auto&t:inverse_products[k])if(t.j)sum+=out.a[t.i]*out.a[t.j];
        out.a[k]=(x.a[k]-sum)/(2*s);
    }
    return out;
}
inline J log(const J &x) {
    Real v = x.a[0];
    if (!(v > 0))
        throw std::domain_error("nonpositive jet logarithm");
    return compose(x, ::tsp_math::working_math::log(v), 1 / v, -0.5 / (v * v), 1 / (3 * v * v * v));
}
inline J abs(const J &x) { return x.value() < 0 ? -x : x; }
inline J hypot(const J&a,const J&b) {
    if(std::abs(a.value())<std::abs(b.value()))return hypot(b,a);
    if(a.a[0]==0)return sqrt(a*a+b*b);
    const J ratio=b/a;
    return abs(a)*sqrt(1+ratio*ratio);
}
inline J hypot(const J&a,const J&b,const J&c) {
    if(std::abs(a.value())<std::abs(b.value()))return hypot(b,a,c);
    if(std::abs(a.value())<std::abs(c.value()))return hypot(c,b,a);
    if(a.a[0]==0)return sqrt(a*a+b*b+c*c);
    const J y=b/a,z=c/a;
    return abs(a)*sqrt(1+y*y+z*z);
}

struct TrigJetEntry {
    J argument, sine, cosine;
    Real s, c;
    bool sine_ready = false, cosine_ready = false;
};
struct TrigJetCache {
    std::array<TrigJetEntry, 32> entries{};
    std::size_t used = 0, replacement = 0;
};
inline thread_local TrigJetCache *active_trig_cache = nullptr;
struct TrigJetScope {
    TrigJetCache cache;
    TrigJetCache *previous = active_trig_cache;
    TrigJetScope() { active_trig_cache = &cache; }
    ~TrigJetScope() { active_trig_cache = previous; }
    TrigJetScope(const TrigJetScope &) = delete;
    TrigJetScope &operator=(const TrigJetScope &) = delete;
};
inline J trig_jet(const J &x, bool cosine) {
    TrigJetEntry *entry = nullptr;
    if (active_trig_cache) {
        auto &cache = *active_trig_cache;
        for (std::size_t k = 0; k < cache.used; ++k) {
            auto &candidate = cache.entries[k];
            if (candidate.argument.a == x.a) { entry = &candidate; break; }
        }
        if (!entry) {

            const Real s = ::tsp_math::working_math::sin(x.a[0]);
            const Real c = ::tsp_math::working_math::cos(x.a[0]);
            const std::size_t k = cache.used < cache.entries.size()
                ? cache.used++ : cache.replacement++ % cache.entries.size();
            entry = &cache.entries[k];
            entry->argument = x;
            entry->s = s;
            entry->c = c;
            entry->sine_ready = entry->cosine_ready = false;
        }
        bool &ready = cosine ? entry->cosine_ready : entry->sine_ready;
        J &out = cosine ? entry->cosine : entry->sine;
        if (!ready) {
            out = cosine ? compose(x, entry->c, -entry->s, -entry->c / 2, entry->s / 6)
                         : compose(x, entry->s, entry->c, -entry->s / 2, -entry->c / 6);
            ready = true;
        }
        return out;
    }
    const Real s = ::tsp_math::working_math::sin(x.a[0]), c = ::tsp_math::working_math::cos(x.a[0]);
    return cosine ? compose(x, c, -s, -c / 2, s / 6)
                  : compose(x, s, c, -s / 2, -c / 6);
}
inline J sin(const J &x) {
    return trig_jet(x, false);
}
inline J cos(const J &x) {
    return trig_jet(x, true);
}
inline J asinh(const J &x) {
    const Real v = x.a[0], d = 1 + v * v, s = ::tsp_math::working_math::sqrt(d);
    return compose(x, ::tsp_math::working_math::asinh(v), 1 / s, -v / (2 * d * s),
                   (2 * v * v - 1) / (6 * d * d * s));
}
inline J atanh(const J &x) {
    const Real v=x.a[0], d=1-v*v;
    return compose(x,::tsp_math::working_math::atanh(v),1/d,v/(d*d),(1+3*v*v)/(3*d*d*d));
}
inline J atan(const J &x) {
    Real v = x.a[0], d = 1 + v * v;
    return compose(x, ::tsp_math::working_math::atan(v), 1 / d, -v / (d * d), (3 * v * v - 1) / (3 * d * d * d));
}
inline J atan2(const J &y, const J &x) {
    J r;
    if (std::abs(x.value()) >= std::abs(y.value()))
        r = atan(y / x);
    else
        r = -atan(x / y);
    r.a[0] = ::tsp_math::working_math::atan2(y.a[0], x.a[0]);
    return r;
}
inline J pow(J x, int n) {
    if (n < 0)
        return inverse(pow(x, -n));
    J r(1);
    for (int k = 0; k < n; ++k)
        r = r * x;
    return r;
}
inline bool isfinite(const J &x) {
    for (const auto &v : x.a)
        if (!::tsp_math::working_math::isfinite(v))
            return false;
    return true;
}

inline int carlson_binary_scale(double largest) {
    if(std::isfinite(largest) && largest>0 && (largest>1e8 || largest<1e-8)) {
        int exponent=0;std::frexp(largest,&exponent);
        return exponent-exponent%2;
    }

    return 0;
}
inline J jet_binary_scale(J v,int exponent) {
    for(auto &c:v.a) c=::tsp_math::working_math::ldexp(c,exponent);
    return v;
}
inline J RF(J x, J y, J z) {
    if(const int e=carlson_binary_scale(std::max({std::abs(x.value()),std::abs(y.value()),std::abs(z.value())})))
        return jet_binary_scale(RF(jet_binary_scale(x,-e),jet_binary_scale(y,-e),jet_binary_scale(z,-e)),-e/2);
    const double tolerance = std::min({carlson_jet_tolerance(x),
        carlson_jet_tolerance(y), carlson_jet_tolerance(z)});
    for (int i = 0; i < 80; ++i) {
        J a = (x + y + z) / 3, dx = (a - x) / a, dy = (a - y) / a, dz = (a - z) / a;
        if (std::max({std::abs(dx.value()), std::abs(dy.value()), std::abs(dz.value())}) < tolerance) {
            J e2 = dx * dy - dz * dz, e3 = dx * dy * dz;
            return (1 + (e2 / 24 - Real(1) / 10 - e3 * (Real(3) / 44)) * e2 + e3 / 14) / sqrt(a);
        }
        J sx = sqrt(x), sy = sqrt(y), sz = sqrt(z), l = sx * (sy + sz) + sy * sz;
        x = (x + l) / 4;
        y = (y + l) / 4;
        z = (z + l) / 4;
    }
    throw std::runtime_error("RF duplication failed");
}
inline J RC(J x, J y) {
    if(const int e=carlson_binary_scale(std::max(std::abs(x.value()),std::abs(y.value()))))
        return jet_binary_scale(RC(jet_binary_scale(x,-e),jet_binary_scale(y,-e)),-e/2);
    const double tolerance = std::min(carlson_jet_tolerance(x), carlson_jet_tolerance(y));
    for (int i = 0; i < 80; ++i) {
        J a = (x + 2 * y) / 3, s = (y - a) / a;
        if (std::abs(s.value()) < tolerance)
            return (1 + s * s * (Real(3) / 10 + s * (Real(1) / 7 + s * (Real(3) / 8 + s * (Real(9) / 22))))) / sqrt(a);
        J l = 2 * sqrt(x) * sqrt(y) + y;
        x = (x + l) / 4;
        y = (y + l) / 4;
    }
    throw std::runtime_error("RC duplication failed");
}
inline J RD(J x, J y, J z) {
    if(const int e=carlson_binary_scale(std::max({std::abs(x.value()),std::abs(y.value()),std::abs(z.value())})))
        return jet_binary_scale(RD(jet_binary_scale(x,-e),jet_binary_scale(y,-e),jet_binary_scale(z,-e)),-3*e/2);
    const double tolerance = std::min({carlson_jet_tolerance(x),
        carlson_jet_tolerance(y), carlson_jet_tolerance(z)});
    J sum;
    double fac = 1;
    for (int i = 0; i < 80; ++i) {
        J a = (x + y + 3 * z) / 5, dx = (a - x) / a, dy = (a - y) / a, dz = (a - z) / a;
        if (std::max({std::abs(dx.value()), std::abs(dy.value()), std::abs(dz.value())}) < tolerance) {
            J ea = dx * dy, eb = dz * dz, ec = ea - eb, ed = ea - 6 * eb, ee = ed + 2 * ec;
            return 3 * sum + fac *
                                 (1 + ed * (-Real(3) / 14 + ed * (Real(9) / 88) - dz * ee * (Real(9) / 52)) +
                                  dz * (ee / 6 + dz * (-ec * (Real(9) / 22) + dz * ea * (Real(3) / 26)))) /
                                 (a * sqrt(a));
        }
        J sx = sqrt(x), sy = sqrt(y), sz = sqrt(z), l = sx * (sy + sz) + sy * sz;
        sum += fac / (sz * (z + l));
        fac *= 0.25;
        x = (x + l) / 4;
        y = (y + l) / 4;
        z = (z + l) / 4;
    }
    throw std::runtime_error("RD duplication failed");
}
inline J RJ(J x, J y, J z, J p) {
    if(const int e=carlson_binary_scale(std::max({std::abs(x.value()),std::abs(y.value()),std::abs(z.value()),std::abs(p.value())})))
        return jet_binary_scale(RJ(jet_binary_scale(x,-e),jet_binary_scale(y,-e),jet_binary_scale(z,-e),jet_binary_scale(p,-e)),-3*e/2);
    const double tolerance = std::min({carlson_jet_tolerance(x),
        carlson_jet_tolerance(y), carlson_jet_tolerance(z), carlson_jet_tolerance(p)});
    J sum;
    double fac = 1;
    for (int i = 0; i < 80; ++i) {
        J a = (x + y + z + 2 * p) / 5, dx = (a - x) / a, dy = (a - y) / a, dz = (a - z) / a, dp = (a - p) / a;
        if (std::max({std::abs(dx.value()), std::abs(dy.value()), std::abs(dz.value()),
                      std::abs(dp.value())}) < tolerance) {
            J ea = dx * (dy + dz) + dy * dz, eb = dx * dy * dz, ec = dp * dp, ed = ea - 3 * ec,
              ee = eb + 2 * dp * (ea - ec);
            return 3 * sum + fac *
                                 (1 + ed * (-Real(3) / 14 + ed * (Real(9) / 88) - ee * (Real(9) / 52)) +
                                  eb * (Real(1) / 6 + dp * (-Real(3) / 11 + dp * (Real(3) / 26))) +
                                  dp * ea * (Real(1) / 3 - dp * (Real(3) / 22)) - dp * ec / 3) /
                                 (a * sqrt(a));
        }
        J sx = sqrt(x), sy = sqrt(y), sz = sqrt(z), l = sx * (sy + sz) + sy * sz,
          al = p * (sx + sy + sz) + sx * sy * sz, be = p * (p + l) * (p + l);
        const J root_be=sqrt(p)*(p+l),rc_lambda=2*al*root_be+be;
        sum += fac * RC((al*al+rc_lambda)/4,(be+rc_lambda)/4);
        fac *= 0.25;
        x = (x + l) / 4;
        y = (y + l) / 4;
        z = (z + l) / 4;
        p = (p + l) / 4;
    }
    throw std::runtime_error("RJ duplication failed");
}
}

namespace tsp_native_math {
inline Real scalar_rf(Real x, Real y, Real z) {
    const double tolerance = 2e-7;
    for (int i = 0; i < 80; ++i) {
        Real a = (x + y + z) / 3, dx = (a - x) / a, dy = (a - y) / a, dz = (a - z) / a;
        if (std::max({std::abs(static_cast<double>(dx)), std::abs(static_cast<double>(dy)), std::abs(static_cast<double>(dz))}) < tolerance) {
            Real e2 = dx * dy - dz * dz, e3 = dx * dy * dz;
            return (1 + (e2 / 24 - Real(1) / 10 - e3 * (Real(3) / 44)) * e2 + e3 / 14) / sqrt(a);
        }
        Real sx = sqrt(x), sy = sqrt(y), sz = sqrt(z), l = sx * (sy + sz) + sy * sz;
        x = (x + l) / 4;
        y = (y + l) / 4;
        z = (z + l) / 4;
    }
    throw std::runtime_error("scalar_rf duplication failed");
}
inline Real scalar_rc(Real x, Real y) {
    const double tolerance = 2e-7;
    for (int i = 0; i < 80; ++i) {
        Real a = (x + 2 * y) / 3, s = (y - a) / a;
        if (std::abs(static_cast<double>(s)) < tolerance)
            return (1 + s * s * (Real(3) / 10 + s * (Real(1) / 7 + s * (Real(3) / 8 + s * (Real(9) / 22))))) / sqrt(a);
        Real l = 2 * sqrt(x) * sqrt(y) + y;
        x = (x + l) / 4;
        y = (y + l) / 4;
    }
    throw std::runtime_error("scalar_rc duplication failed");
}
inline Real scalar_rd(Real x, Real y, Real z) {
    const double tolerance = 2e-7;
    Real sum=0;
    double fac = 1;
    for (int i = 0; i < 80; ++i) {
        Real a = (x + y + 3 * z) / 5, dx = (a - x) / a, dy = (a - y) / a, dz = (a - z) / a;
        if (std::max({std::abs(static_cast<double>(dx)), std::abs(static_cast<double>(dy)), std::abs(static_cast<double>(dz))}) < tolerance) {
            Real ea = dx * dy, eb = dz * dz, ec = ea - eb, ed = ea - 6 * eb, ee = ed + 2 * ec;
            return 3 * sum + fac *
                                 (1 + ed * (-Real(3) / 14 + ed * (Real(9) / 88) - dz * ee * (Real(9) / 52)) +
                                  dz * (ee / 6 + dz * (-ec * (Real(9) / 22) + dz * ea * (Real(3) / 26)))) /
                                 (a * sqrt(a));
        }
        Real sx = sqrt(x), sy = sqrt(y), sz = sqrt(z), l = sx * (sy + sz) + sy * sz;
        sum += fac / (sz * (z + l));
        fac *= 0.25;
        x = (x + l) / 4;
        y = (y + l) / 4;
        z = (z + l) / 4;
    }
    throw std::runtime_error("scalar_rd duplication failed");
}
inline Real scalar_rj(Real x, Real y, Real z, Real p) {
    const double tolerance = 2e-7;
    Real sum=0;
    double fac = 1;
    for (int i = 0; i < 80; ++i) {
        Real a = (x + y + z + 2 * p) / 5, dx = (a - x) / a, dy = (a - y) / a, dz = (a - z) / a, dp = (a - p) / a;
        if (std::max({std::abs(static_cast<double>(dx)), std::abs(static_cast<double>(dy)), std::abs(static_cast<double>(dz)),
                      std::abs(static_cast<double>(dp))}) < tolerance) {
            Real ea = dx * (dy + dz) + dy * dz, eb = dx * dy * dz, ec = dp * dp, ed = ea - 3 * ec,
              ee = eb + 2 * dp * (ea - ec);
            return 3 * sum + fac *
                                 (1 + ed * (-Real(3) / 14 + ed * (Real(9) / 88) - ee * (Real(9) / 52)) +
                                  eb * (Real(1) / 6 + dp * (-Real(3) / 11 + dp * (Real(3) / 26))) +
                                  dp * ea * (Real(1) / 3 - dp * (Real(3) / 22)) - dp * ec / 3) /
                                 (a * sqrt(a));
        }
        Real sx = sqrt(x), sy = sqrt(y), sz = sqrt(z), l = sx * (sy + sz) + sy * sz,
          al = p * (sx + sy + sz) + sx * sy * sz, be = p * (p + l) * (p + l);
        sum += fac * scalar_rc(al * al, be);
        fac *= 0.25;
        x = (x + l) / 4;
        y = (y + l) / 4;
        z = (z + l) / 4;
        p = (p + l) / 4;
    }
    throw std::runtime_error("scalar_rj duplication failed");
}
}

namespace TSP { namespace detail {
using J=tsp_math::J;
using Real=tsp_math::Real;
struct EllipticValues {J F,E,P;};
template<int K> inline std::array<J,K+1> elliptic_moments(J theta) {
    std::array<J,K+1> I{};

    const Real s=sin(theta.a[0]),c=cos(theta.a[0]),s2=s*s,c2=c*c;
    I[0]=theta;
    Real odd=s,previous_power=1;
    auto lift=[&](int k,const Real &primitive) {
        const Real d1=previous_power*s2,d2=k*previous_power*s*c;
        const Real d3=Real(k)*previous_power*((2*k-1)*c2-s2)/3;
        I[k]=compose(theta,primitive,d1,d2,d3);
        previous_power=d1;
    };
    if (std::abs(theta.value())<=0.25) {

        constexpr int L=32;
        static const auto binom=[] {
            std::array<Real,L+1> out{};out[0]=1;
            for(int j=1;j<=L;++j)out[j]=out[j-1]*Real(2*j-1)/Real(2*j);
            return out;
        }();
        static const auto coefficients=[] {
            std::array<std::array<Real,L+1>,K+1> out{};
            for(int k=1;k<=K;++k)for(int j=0;j<=L;++j)
                out[k][j]=binom[j]/Real(2*k+2*j+1);
            return out;
        }();
        for (int k=1;k<=K;++k) {
            odd*=s2;
            Real h=coefficients[k][L];
            for (int j=L-1;j>=0;--j)
                h=coefficients[k][j]+s2*h;
            lift(k,odd*h);
        }
    } else {
        Real primitive=theta.a[0];
        for (int k=1;k<=K;++k) {
            primitive=primitive*(Real(2*k-1)/Real(2*k))-odd*c*(Real(1)/Real(2*k));
            lift(k,primitive);
            odd*=s2;
        }
    }
    return I;
}

template<int K> inline const std::array<std::array<Real,K+1>,2> &FEBinomialCoefficients() {
    static const auto values=[] {
        std::array<std::array<Real,K+1>,2> out{};
        out[0][0]=out[1][0]=1;
        for(int k=1;k<=K;++k) {
            out[0][k]=out[0][k-1]*Real(2*k-1)/Real(2*k);
            out[1][k]=out[1][k-1]*Real(2*k-3)/Real(2*k);
        }
        return out;
    }();
    return values;
}

template<int K,bool ComputePi=true> inline EllipticValues SmallParameterJets(J theta,J m,J n) {

    const auto I=elliptic_moments<K>(theta);
    J F=theta,E=theta,Pi=ComputePi?theta:J(0.0),mpower(1.0),h(1.0);
    const auto &coeff=FEBinomialCoefficients<K>();
    for (int k=1;k<=K;++k) {
        mpower*=m;
        const Real &cf=coeff[0][k];
        if constexpr(ComputePi) h=n*h+cf*mpower;
        F+=cf*mpower*I[k];
        E+=coeff[1][k]*mpower*I[k];
        if constexpr(ComputePi) Pi+=h*I[k];
    }
    return {F,E,Pi};
}

inline J EulerJet(const J&x){J out=x;for(int k=0;k<20;++k)out.a[k]*=tsp_math::ex[k].x+tsp_math::ex[k].y+tsp_math::ex[k].z;return out;}
inline void AssignDegree(J&out,const J&rhs,int d){for(int k=1;k<20;++k)if(tsp_math::ex[k].x+tsp_math::ex[k].y+tsp_math::ex[k].z==d)out.a[k]=rhs.a[k]/d;}
inline bool zero_jet(const J&x){for(const auto&a:x.a)if(a!=0)return false;return true;}
inline std::pair<J,J> EllipticFE(const J&u,const J&m) {
    struct Cache{bool valid=false;std::uint32_t mask=0;J u,m,F,E;};static thread_local Cache cache;
    if(cache.valid&&cache.mask==tsp_math::active_jet_mask&&cache.u.a==u.a&&cache.m.a==m.a)return {cache.F,cache.E};
    J F,E;
    if(std::abs(m.value())<.01) {
        const auto out=SmallParameterJets<16,false>(u,m,J(0));F=out.F;E=out.E;
    } else {
        const Real s=sin(u.a[0]),c=cos(u.a[0]),s2=s*s,B=1-m.a[0]*s2;
        const Real seed=s*tsp_native_math::scalar_rf(c*c,B,Real(1));
        F=J(seed);E=J(seed-m.a[0]*s*s2*tsp_native_math::scalar_rd(c*c,B,Real(1))/3);
        const J sj=sin(u),cj=cos(u),root=sqrt(1-m*sj*sj),em=EulerJet(m),et=EulerJet(u);
        const J mb=em/(2*m),ma=mb/(1-m),base=et/root-em*sj*cj/(2*(1-m)*root),base_e=et*root;
        const auto mask=tsp_math::active_jet_mask;
        for(int d=1;d<=3;++d) {
            if(!(mask>>(d==1?1:d==2?4:10)))break;
            tsp_math::JetScope scope(mask&(d==1?15u:d==2?1023u:tsp_math::JET_FULL_MASK));
            const J rf=base+ma*E-mb*F,re=base_e+mb*(E-F);AssignDegree(F,rf,d);AssignDegree(E,re,d);
        }
    }
    cache={true,tsp_math::active_jet_mask,u,m,F,E};return {F,E};
}
inline EllipticValues CarlsonElliptic(const J&u,const J&m,const J&n,bool complete) {
    const J s=complete?J(1):sin(u),c=complete?J(0):cos(u),s2=s*s;
    const J f=s*tsp_math::RF(c*c,1-m*s2,1);
    return {f,f-m*s*s2*tsp_math::RD(c*c,1-m*s2,1)/3,
        f+n*s*s2*tsp_math::RJ(c*c,1-m*s2,1,1-n*s2)/3};
}
inline EllipticValues EllipticSeeds(const J&u,const J&m,const J&n,bool complete) {
    if(std::abs((1-n).value())<.001) {
        if(m.value()>0)return CarlsonElliptic(u,m,n,complete);
        const auto fe=EllipticFE(u,m);
        const J s=complete?J(1):sin(u),c=complete?J(0):cos(u),s2=s*s;
        return {fe.first,fe.second,fe.first+n*s*s2*tsp_math::RJ(c*c,1-m*s2,1,1-n*s2)/3};
    }
    if(std::max(std::abs(m.value()),std::abs(n.value()))<.01){return SmallParameterJets<16>(u,m,n);}
    const auto fe=EllipticFE(u,m);const J F=fe.first,E=fe.second;
    if(zero_jet(n))return {F,E,F};
    const J s=sin(u),c=cos(u),s2=s*s,root=sqrt(1-m*s2),delta=n-m;
    if(std::abs(m.value())>.01&&tsp_native_math::abs(delta.a[0])<=Real(".001")*tsp_native_math::abs(m.a[0])) {
        constexpr int K=14;std::array<J,K+1>p{};const J inverse_den=1/(1-m*s2),ratio=s2*inverse_den;
        p[0]=(E-m*s*c/root)/(1-m);const J boundary=root*s*c,inv=1/(m*(m-1)),m2=m*m,twice=2*m;
        J previous2(0),previous=inverse_den;
        for(int k=1;k<=K;++k){const J h=previous*ratio,q=m2*h+twice*previous+(k>=2?previous2:J(0));
            J rhs=(k==1?E-F:J(0))+twice*p[k-1]-boundary*q;if(k>=2)rhs+=(2*k-3)*p[k-2];rhs+=2*(k-1)*(2*m-1)*p[k-1];
            p[k]=-rhs*inv/(2*k+1);previous2=previous;previous=h;}
        J out=p[K];for(int k=K-1;k>=0;--k)out=p[k]+delta*out;return {F,E,out};
    }
    if(std::abs(n.value())<.01&&std::abs(m.value())>.01) {
        constexpr int K=14;std::array<J,K+1>I{};I[0]=F;I[1]=(F-E)/m;
        J boundary=s*c*root;
        for(int k=1;k<K;++k){I[k+1]=(boundary-(2*k-1)*I[k-1]+2*k*(1+m)*I[k])/((2*k+1)*m);boundary*=s2;}
        J out=I[K];for(int k=K-1;k>=0;--k)out=I[k]+n*out;return {F,E,out};
    }
    if(std::abs(n.value())<1e-12||std::abs((m-n).value())<1e-14||std::abs((1-n).value())<1e-12)
        return CarlsonElliptic(u,m,n,complete);
    const Real ss=sin(u.a[0]),cc=cos(u.a[0]),ss2=ss*ss;
    J out(F.a[0]+n.a[0]*ss*ss2*tsp_native_math::scalar_rj(cc*cc,1-m.a[0]*ss2,Real(1),1-n.a[0]*ss2)/3);
    const J em=EulerJet(m),en=EulerJet(n),et=EulerJet(u),den=1-n*s2;
    const J bm=em/(2*(m-n)),bn=en/(2*(m-n)*(n-1)*n);
    const J base=et/(den*root)+bm*(-E/(m-1)+m*s*c/((m-1)*root))+
        bn*(n*E+(m-n)*F-n*n*root*s*c/den),coefficient=-bm+bn*(n*n-m);
    const auto mask=tsp_math::active_jet_mask;
    for(int d=1;d<=3;++d) {
        if(!(mask>>(d==1?1:d==2?4:10)))break;
        tsp_math::JetScope scope(mask&(d==1?15u:d==2?1023u:tsp_math::JET_FULL_MASK));
        AssignDegree(out,base+coefficient*out,d);
    }
    return {F,E,out};
}
}}
#include <fstream>
#include <string>

namespace TSP {
constexpr double PI = 3.141592653589793238462643383279502884;
using Vec = std::array<double,3>;
using Gravity = std::array<double,20>;
using Magnetic = std::array<double,10>;
struct Prism {
    std::array<Vec,3> vertices{};
    double r1=0, r2=0, density=1;
};
enum class Method { Automatic, Elliptic };
struct Options {
    Method method=Method::Elliptic;
    double relative_tolerance=1e-16, absolute_tolerance=1e-16;
};
struct Diagnostics {
    Method method=Method::Elliptic;
    std::size_t evaluations=0;
    double estimated_error=0;
};
namespace detail {
using tsp_math::J;
using tsp_math::Real;
using RV=std::array<Real,3>;
inline Vec cross(const Vec&a,const Vec&b) {return {a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]};}
inline double dot(const Vec&a,const Vec&b) {return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];}
inline Vec scale(Vec a,double s) {for(auto&v:a)v*=s;return a;}
inline Vec add(Vec a,const Vec&b) {for(int i=0;i<3;++i)a[i]+=b[i];return a;}
inline double norm(const Vec&a) {const double scale=std::max({std::abs(a[0]),std::abs(a[1]),std::abs(a[2])});
    if(scale==0)return 0;
    const double x=a[0]/scale,y=a[1]/scale,z=a[2]/scale;
    return scale*std::sqrt(x*x+y*y+z*z);}
inline Vec unit(Vec a) {double n=norm(a);if(!(n>0)||!std::isfinite(n))throw std::invalid_argument("invalid direction");return scale(a,1/n);}
inline Real rdot(const RV&a,const RV&b) {return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];}
inline RV rcross(const RV&a,const RV&b) {return {a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]};}
inline RV rscale(RV a,Real s) {for(auto&v:a)v*=s;return a;}
inline Real rnorm(const RV&a) {return tsp_native_math::sqrt(rdot(a,a));}
inline RV runit(RV a) {return rscale(a,1/rnorm(a));}
inline RV convert(const Vec&a) {return {Real(a[0]),Real(a[1]),Real(a[2])};}
inline J jdot(const RV&a,const std::array<J,3>&b) {return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];}
inline int index(int x,int y,int z) {
    for(int i=0;i<20;++i)if(tsp_math::ex[i].x==x&&tsp_math::ex[i].y==y&&tsp_math::ex[i].z==z)return i;
    throw std::logic_error("invalid derivative index");
}
inline int factorial(int n) {return n==3?6:n==2?2:1;}
struct Edge {
    RV a,b,N;Real gamma;
    Edge()=default;
    Edge(const RV&start,const RV&end):a(start) {
        N=runit(rcross(a,end));gamma=tsp_native_math::atan2(rnorm(rcross(a,end)),rdot(a,end));
        if(!(gamma>0&&gamma<tsp_math::PI-Real("1e-14")))throw std::invalid_argument("minor nonantipodal arcs required");
        b=rcross(N,a);
    }
};
struct Geometry {
    std::array<RV,3> v;
    std::array<Edge,3> edges;
    Real det;
    explicit Geometry(const Prism&p) {
        if(!(p.r1>0&&p.r2>p.r1)||!std::isfinite(p.r1)||!std::isfinite(p.r2)||!std::isfinite(p.density))
            throw std::invalid_argument("finite 0 < R1 < R2 and density required");
        for(int i=0;i<3;++i) {
            for(double x:p.vertices[i])if(!std::isfinite(x))throw std::invalid_argument("nonfinite vertex");
            if(norm(p.vertices[i])==0)throw std::invalid_argument("zero vertex");
            v[i]=runit(convert(p.vertices[i]));
        }
        det=rdot(v[0],rcross(v[1],v[2]));
        if(det<0) {std::swap_ranges(v[1].begin(),v[1].end(),v[2].begin());det=-det;}
        if(det<Real("1e-26"))throw std::invalid_argument("degenerate or excessively ill-conditioned triangle");
        const Real area=2*tsp_native_math::atan2(det,1+rdot(v[0],v[1])+rdot(v[1],v[2])+rdot(v[2],v[0]));
        if(!(area>0&&area<2*tsp_math::PI))throw std::invalid_argument("triangle must occupy an open hemisphere");
        for(int i=0;i<3;++i)edges[i]=Edge(v[i],v[(i+1)%3]);
    }
};
inline Real solid_angle(const Geometry&g) {
    return 2*tsp_native_math::atan2(g.det,1+rdot(g.v[0],g.v[1])+rdot(g.v[1],g.v[2])+rdot(g.v[2],g.v[0]));
}
template<int K> inline std::array<J,K+1> edge_moments(const J&xi,const J&eta,const Real&gamma) {
    std::array<J,K+1> out;out[0]=J(gamma);
    const Real s=tsp_native_math::sin(gamma),c=tsp_native_math::cos(gamma);
    out[1]=xi*s+eta*(1-c);
    const J end=xi*c+eta*s,derivative=-xi*s+eta*c,r2=xi*xi+eta*eta;
    J start_power(1),end_power(1);
    for(int j=2;j<=K;++j) {
        start_power*=xi;end_power*=end;
        out[j]=((j-1)*r2*out[j-2]-end_power*derivative+start_power*eta)/j;
    }
    return out;
}
inline J solid_harmonic_series(const Geometry&g,const Real&a,const std::array<J,3>&p,bool inner) {
    constexpr int K=32;
    const J r2=p[0]*p[0]+p[1]*p[1]+p[2]*p[2];
    const J radius=inner?J(1):tsp_math::sqrt(r2);
    std::array<J,3> q;for(int i=0;i<3;++i)q[i]=inner?p[i]/a:p[i]/r2;
    const J q2=q[0]*q[0]+q[1]*q[1]+q[2]*q[2];
    std::array<std::array<J,K+1>,3> moments;
    std::array<J,3> normal;
    for(int e=0;e<3;++e) {
        const auto&edge=g.edges[e];normal[e]=jdot(edge.N,q);
        moments[e]=edge_moments<K>(jdot(edge.a,q),jdot(edge.b,q),edge.gamma);
    }
    std::array<J,K/2+1> powers;powers[0]=J(1);
    for(int j=1;j<=K/2;++j)powers[j]=powers[j-1]*q2;
    const Real area=solid_angle(g);
    J total=area*(inner?(1-a*a)/2:(1-a*a*a)/3);
    const int last=inner&&p[0].a[0]==0&&p[1].a[0]==0&&p[2].a[0]==0?3:K;
    Real leading=1;
    for(int l=1;l<=last;++l) {
        leading*=Real(2*l-1)/l;
        Real coefficient=leading;
        J moment;
        for(int j=0;j<=(l-1)/2;++j) {
            const int degree=l-2*j;
            J boundary;for(int e=0;e<3;++e)boundary+=normal[e]*moments[e][degree-1];
            moment+=coefficient*degree*powers[j]*boundary;
            coefficient*=-Real((l-2*j)*(l-2*j-1))/Real(2*(j+1)*(2*l-2*j-1));
        }
        Real weight;
        if(!inner)weight=(1-tsp_native_math::pow(a,l+3))/(l+3);
        else if(l==1)weight=a*(1-a);
        else if(l==2)weight=-a*a*tsp_native_math::log(a);
        else weight=a*a*(1-tsp_native_math::pow(a,l-2))/(l-2);
        total+=weight*moment/(l*(l+1));
    }
    return total/radius;
}
inline void validate_point(const Vec&p) {for(double x:p)if(!std::isfinite(x))throw std::invalid_argument("nonfinite observer");}
inline bool inside_angle(const Geometry&g,const RV&n,bool strict=false) {
    for(const auto&e:g.edges)if(rdot(e.N,n)<=(strict?Real(0):Real("-1e-28")))return false;
    return true;
}
inline bool is_boundary(const Geometry&g,double r1,const Vec&p) {
    const RV q=convert(p);const Real R=rnorm(q);if(R==0)return false;
    const RV n=rscale(q,1/R);const Real eps("2e-14");
    if(R>=r1-eps&&R<=1+eps) {
        if((tsp_native_math::abs(R-r1)<eps||tsp_native_math::abs(R-1)<eps)&&inside_angle(g,n))return true;
        for(const auto&e:g.edges)if(tsp_native_math::abs(rdot(e.N,n))<eps) {
            Real t=tsp_native_math::atan2(rdot(n,e.b),rdot(n,e.a));
            if(t>=-eps&&t<=e.gamma+eps)return true;
        }
    }
    return false;
}

inline EllipticValues EvaluateElliptic(const J&u,const J&m,const J&nu) {
    auto elliptic=[&](const J&theta,bool complete) {
        return EllipticSeeds(complete?J(tsp_math::PI/2):theta,m,nu,complete);
    };
    const int periods=static_cast<int>(std::floor((u.value()+PI/2)/PI));
    J v=u-periods*tsp_math::PI;
    if(std::abs(std::abs(v.value())-PI/2)<1e-10) {
        const int sign=v.value()>0?1:-1;
        auto out=elliptic(v,true);J d=v-sign*tsp_math::PI/2;
        const J root=tsp_math::sqrt(1-m),den=1-nu;
        out.F=sign*out.F+d/root-m*tsp_math::pow(d,3)/(6*tsp_math::pow(root,3));
        out.E=sign*out.E+d*root+m*tsp_math::pow(d,3)/(6*root);
        out.P=sign*out.P+d/(den*root)-tsp_math::pow(d,3)*(nu/den+m/(2*(1-m)))/(3*den*root);
        if(periods) {auto c=elliptic(v,true);out.F+=2*periods*c.F;out.E+=2*periods*c.E;out.P+=2*periods*c.P;}
        return out;
    }
    auto out=elliptic(v,false);
    if(periods) {auto c=elliptic(v,true);out.F+=2*periods*c.F;out.E+=2*periods*c.E;out.P+=2*periods*c.P;}
    return out;
}
inline J RegularEdge(const Real&a,const J&phi,const J&h,const J&k,const J&R) {
    const J k2=k*k;
    const J one_minus_h=h.a[0]>=0?k2/(1+h):1-h;
    const J one_plus_h=h.a[0]>=0?1+h:k2/(1-h);
    const J sh=tsp_math::sin(phi/2),ch=tsp_math::cos(phi/2),ss=sh*sh,cc=ch*ch;
    const J c=h*(cc-ss),sn=2*sh*ch,distance=one_minus_h*cc+one_plus_h*ss;
    const J root_minus=h.a[0]>=0?k/tsp_math::sqrt(1+h):tsp_math::sqrt(1-h);
    const J root_plus=h.a[0]>=0?tsp_math::sqrt(1+h):k/tsp_math::sqrt(1-h);
    auto separation=[&](const Real&s) {
        const J scale=tsp_math::sqrt(2*s*R);
        return tsp_math::hypot(R-s,scale*root_minus*ch,scale*root_plus*sh);
    };
    const J D1=separation(a),D2=separation(Real(1));
    J logarithm;
    if(c.a[0]<0||R.a[0]<a)logarithm=tsp_math::log((1-R*c+D2)/(a-R*c+D1));
    else if(R.a[0]>=1)logarithm=tsp_math::log((D1+R*c-a)/(D2+R*c-1));
    else logarithm=tsp_math::log(1-R*c+D2)-tsp_math::log(R*R*distance*(1+c)/(D1+R*c-a));
    J out=k*tsp_math::pow(R,3)*h*sn*logarithm/2;
    auto endpoint=[&](const Real&s,const J&D) {
        const J z=s-R,t=s+R,d=tsp_math::hypot(z,tsp_math::sqrt(2*s*R)*root_minus),m=-4*s*R*h/(d*d);
        const J complement2=(z*z+2*s*R*one_plus_h)/(t*t);
        const J A=z*(z*z*(2*s+R)-3*tsp_math::pow(R,3)*k2)/12;
        const J B=t*(t*t*(2*s-R)+3*tsp_math::pow(R,3)*k2)/12;
        const J root_c=tsp_math::abs(ch),root_y=D/d,y=root_y*root_y;
        const J lambda=root_c*(root_y+1)+root_y;
        const J xx=(cc+lambda)/4,yy=(y+lambda)/4,zz=(1+lambda)/4;
        auto remainder=[&](const J&p,const J&root_p) {
            const J al=p*(root_c+root_y+1)+root_c*root_y,root_be=root_p*(p+lambda),be=root_be*root_be;
            const J rc_lambda=2*al*root_be+be;
            return sh*ss*(tsp_math::RC((al*al+rc_lambda)/4,(be+rc_lambda)/4)
                +tsp_math::RJ(xx,yy,zz,(p+lambda)/4)/12);
        };
        const J p1=cc+ss*z*z/(d*d),p2=cc+ss*complement2;
        const J S1=remainder(p1,tsp_math::hypot(ch,z*sh/d));
        const J S2=remainder(p2,tsp_math::hypot(ch,z*sh/t,tsp_math::sqrt(2*s*R)*root_plus*sh/t));
        EllipticValues fe;
        if(m.a[0]>0||m.a[0]<Real("-1e10")) {
            fe.F=sh*tsp_math::RF(xx,yy,zz);
            fe.E=fe.F-m*sh*ss*(1/(1+lambda)+tsp_math::RD(xx,yy,zz)/12);
        } else fe=EvaluateElliptic(phi/2,m,J(0));
        const J ab=tsp_math::abs(z);
        const J dp=tsp_math::atan2(-s*R*k*sn,ab*(D+ab)+2*s*R*one_minus_h*cc);
        const J dm=tsp_math::atan2(-s*R*k*sn,D*(t+D)+2*s*R*one_plus_h*cc);
        const J Cp=(z.a[0]>=0?1:-1)*(z*z*(2*s+R)-3*tsp_math::pow(R,3)*k2)/6;
        const J Cm=(t*t*(2*s-R)+3*tsp_math::pow(R,3)*k2)/6;
        return k*(R*R*d*fe.E/2-R*R*(3*R*R+s*s)*fe.F/(6*d)
                  -4*A*s*R*S1/tsp_math::pow(d,3)+4*B*s*R*(fe.F-complement2*S2)/(t*t*d))+Cp*dp+Cm*dm;
    };
    return out+endpoint(1,D2)-endpoint(a,D1);
}
inline J SmallProjectionEdge(const Edge&e,const Real&a,const J&R,const J&xi,const J&eta,const J&k,const J&jp,const J&jm) {
    constexpr int K=20;
    const auto moments=edge_moments<K>(xi,eta,e.gamma);
    auto primitive=[&](const Real&s) {
        std::array<J,K+1> d,z,logarithm,out;
        const J r2=R*R,b=2*s*R,A=r2+s*s,ratio=-b/A;
        d[0]=tsp_math::sqrt(A);z[0]=d[0]+s;logarithm[0]=tsp_math::log(z[0]);
        for(int j=1;j<=K;++j) {
            d[j]=d[j-1]*ratio*Real(3-2*j)/Real(2*j);
            z[j]=d[j]-(j==1?R:J(0));
            J sum;for(int i=1;i<j;++i)sum+=i*logarithm[i]*z[j-i];
            logarithm[j]=(z[j]-sum/j)/z[0];
        }
        out[0]=A*d[0]/3;
        for(int j=1;j<=K;++j)out[j]=(A*d[j]-b*d[j-1])/3
            +R*(s*d[j-1]-(j>=2?R*d[j-2]:J(0))+r2*(logarithm[j-1]-(j>=3?logarithm[j-3]:J(0))))/2;
        return out;
    };
    const auto upper=primitive(Real(1)),lower=primitive(a);
    std::array<J,K+1> kernel;J total;
    for(int j=0;j<=K;++j) {
        kernel[j]=upper[j]-lower[j]+(j>=2?kernel[j-2]:J(0));
        total+=(kernel[j]-(jp+(j%2?-jm:jm))/2)*moments[j];
    }
    return k*total;
}
inline std::array<Real,20> tensor_derivatives(const J&value) {
    if(!tsp_math::isfinite(value))throw std::runtime_error("nonfinite analytic tensor");
    std::array<Real,20>out{};
    for(int i=0;i<20;++i)if(tsp_math::jet_active(i)){const auto e=tsp_math::ex[i];out[i]=value.a[i]*(factorial(e.x)*factorial(e.y)*factorial(e.z));}return out;
}
struct EllipticKernel {
    Real r1;J R,inverseR,jp,jm;std::array<J,3>n;int sign_change;
    EllipticKernel(const Real&a,const RV&p):r1(a) {
        std::array<J,3>x;for(int i=0;i<3;++i)x[i]=J::variable(p[i],i);
        R=tsp_math::sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);inverseR=tsp_math::inverse(R);
        for(int i=0;i<3;++i)n[i]=x[i]*inverseR;
        auto pole=[&](const Real&s){const J z=s-R;const int sign=z.a[0]>=0?1:-1;return sign*(tsp_math::pow(z,3)/3+R*z*z/2);};
        jp=pole(Real(1))-pole(r1);jm=R*(1-r1*r1)/2+(1-r1*r1*r1)/3;
        sign_change=(R.a[0]<=1?1:-1)-(R.a[0]<=r1?1:-1);
    }
    void AddEdge(J&total,const Edge&e)const {
        const J k=jdot(e.N,n),xi=jdot(e.a,n),eta=jdot(e.b,n);
        if(std::hypot(xi.value(),eta.value())<.01){total+=SmallProjectionEdge(e,r1,R,xi,eta,k,jp,jm);return;}
        int phase=(xi*tsp_native_math::cos(e.gamma/2)+eta*tsp_native_math::sin(e.gamma/2)).a[0]>=0?1:-1;
        Real projection=tsp_native_math::atan2(eta.a[0],xi.a[0]);
        if(projection<0)projection+=2*tsp_math::PI;
        if(projection>e.gamma+Real("1e-24")&&projection<2*tsp_math::PI-Real("1e-24")&&std::abs(k.value())<1e-4&&
            std::min(std::abs(R.value()-1),std::abs(R.value()-static_cast<double>(r1)))<1e-5)phase=-1;
        const J h=phase*tsp_math::sqrt(xi*xi+eta*eta);
        J delta=tsp_math::atan2(phase*eta,phase*xi);
        delta+=2*tsp_math::PI*std::round((static_cast<double>(e.gamma)/2-delta.value())/(2*PI));
        const J lo=-delta,hi=J(e.gamma)-delta;
        total+=RegularEdge(r1,hi,h,k,R)-RegularEdge(r1,lo,h,k,R);
        if(sign_change) {
            const J ap=h.a[0]>=0?tsp_math::atan2(k*(1+h)*tsp_math::sin((hi-lo)/2),
                k*k*tsp_math::cos(hi/2)*tsp_math::cos(lo/2)+(1+h)*(1+h)*tsp_math::sin(hi/2)*tsp_math::sin(lo/2)):
                tsp_math::atan2(k*(1-h)*tsp_math::sin((hi-lo)/2),
                (1-h)*(1-h)*tsp_math::cos(hi/2)*tsp_math::cos(lo/2)+k*k*tsp_math::sin(hi/2)*tsp_math::sin(lo/2));
            total-=sign_change*tsp_math::pow(R,3)*k*k*ap/2;
        }
    }
};
inline std::array<Real,20> elliptic_tensor(const Geometry&g,const Real&r1,const RV&p) {
    const double radius=static_cast<double>(rnorm(p));
    if(radius<.03*static_cast<double>(r1)||radius>20){
        std::array<J,3>x;for(int i=0;i<3;++i)x[i]=J::variable(p[i],i);
        return tensor_derivatives(solid_harmonic_series(g,r1,x,radius<static_cast<double>(r1)));}
    const EllipticKernel kernel(r1,p);J total=solid_angle(g)*(kernel.jm-kernel.jp)/2;
    for(const auto&e:g.edges)kernel.AddEdge(total,e);return tensor_derivatives(total*kernel.inverseR);
}

inline bool harmonic_check(const std::array<Real,20>&a,bool interior,const Options&o) {
    auto check=[&](int i,int j,int k,const Real&target) {
        const Real budget=3*Real(o.absolute_tolerance)+Real(o.relative_tolerance)*
            (tsp_native_math::abs(a[i])+tsp_native_math::abs(a[j])+tsp_native_math::abs(a[k]));
        return tsp_native_math::abs(a[i]+a[j]+a[k]-target)<=budget;
    };
    if((tsp_math::active_jet_mask&1023)==1023&&!check(4,7,9,interior?-4*tsp_math::PI:Real(0)))return false;
    if(tsp_math::active_jet_mask==tsp_math::JET_FULL_MASK)for(auto ix:std::array<std::array<int,3>,3>{{{10,13,15},{11,16,18},{12,17,19}}})
        if(!check(ix[0],ix[1],ix[2],Real(0)))return false;
    return true;
}
template<class T> inline std::array<T,10> contract(const std::array<T,20>&g,const std::array<T,3>&m) {
    constexpr int indices[10][3]={{1,2,3},{4,5,6},{5,7,8},{6,8,9},{10,11,12},
        {11,13,14},{13,16,17},{12,14,15},{14,17,18},{15,18,19}};
    std::array<T,10> out{};for(int k=0;k<10;++k)if(tsp_math::jet_active(indices[k][0]))out[k]=m[0]*g[indices[k][0]]+m[1]*g[indices[k][1]]+m[2]*g[indices[k][2]];
    return out;
}
}

inline Vec Direction(double latitude,double longitude) {
    if(!std::isfinite(latitude)||!std::isfinite(longitude)||std::abs(latitude)>90)throw std::invalid_argument("invalid latitude/longitude");
    double f=latitude*PI/180,l=std::remainder(longitude,360.0)*PI/180;
    return {std::cos(f)*std::cos(l),std::cos(f)*std::sin(l),std::sin(f)};
}
inline Prism MakePrism(double Fai1,double Fai2,double Fai3,double Lamda1,double Lamda2,double Lamda3,double R1,double R2,double density=1) {
    Prism p{{Direction(Fai1,Lamda1),Direction(Fai2,Lamda2),Direction(Fai3,Lamda3)},R1,R2,density};detail::Geometry validate(p);return p;
}
namespace detail {
using WideGravity=std::array<Real,20>;
using WideAxes=std::array<RV,3>;
template<std::size_t N> inline std::array<double,N> rounded(const std::array<Real,N>&a) {
    std::array<double,N> out;
    for(std::size_t k=0;k<N;++k){out[k]=static_cast<double>(a[k]);if(!std::isfinite(out[k]))throw std::overflow_error("nonfinite field after dimensional scaling");}
    return out;
}
inline WideAxes local_axes(double latitude,double longitude) {
    Direction(latitude,longitude);
    const Real f=Real(latitude)*tsp_math::PI/180,l=Real(std::remainder(longitude,360.0))*tsp_math::PI/180;
    const Real sf=tsp_native_math::sin(f),cf=tsp_native_math::cos(f),sl=tsp_native_math::sin(l),cl=tsp_native_math::cos(l);
    return {{{-sf*cl,-sf*sl,cf},{-sl,cl,Real(0)},{cf*cl,cf*sl,sf}}};
}
inline WideAxes local_axes(const RV&p) {
    const Real r=rnorm(p),h=tsp_native_math::sqrt(p[0]*p[0]+p[1]*p[1]);
    if(h==0)return local_axes(p[2]<0?-90.:p[2]>0?90.:0.,0.);
    return {{{-p[2]*p[0]/(r*h),-p[2]*p[1]/(r*h),h/r},{-p[1]/h,p[0]/h,Real(0)},{p[0]/r,p[1]/r,p[2]/r}}};
}
inline WideGravity rotate(const WideGravity&a,const WideAxes&axes) {
    WideGravity out{};out[0]=a[0];
    for(int k=1;k<20;++k)if(tsp_math::jet_active(k)) {
        auto e=tsp_math::ex[k];int target[3],degree=0;
        for(int j=0;j<e.x;++j)target[degree++]=0;for(int j=0;j<e.y;++j)target[degree++]=1;for(int j=0;j<e.z;++j)target[degree++]=2;
        const int count=degree==1?3:degree==2?9:27;
        for(int code=0;code<count;++code){int v=code,ex[3]={};Real factor(1);
            for(int j=0;j<degree;++j){int axis=v%3;v/=3;++ex[axis];factor*=axes[target[j]][axis];}
            out[k]+=factor*a[index(ex[0],ex[1],ex[2])];}
    }
    return out;
}
inline WideGravity gravity_tensor(const Prism&prism,const Vec&point,const Options&options,Diagnostics*diagnostics) {
    validate_point(point);
    if(!(options.relative_tolerance>0&&options.absolute_tolerance>0)||!std::isfinite(options.relative_tolerance)||!std::isfinite(options.absolute_tolerance))
        throw std::invalid_argument("invalid analytic options");
    Geometry g(prism);const double length=prism.r2,r1=prism.r1/length;const Vec p=scale(point,1/length);
    if(is_boundary(g,r1,p))throw std::domain_error("full third-order tensors are undefined on a physical boundary (2e-14 relative guard)");
    auto a=elliptic_tensor(g,Real(prism.r1)/Real(length),rscale(convert(point),1/Real(length)));
    const double R=norm(p);const bool interior=R>r1&&R<1&&inside_angle(g,runit(convert(p)),true);
    if(!harmonic_check(a,interior,options))throw std::runtime_error("elliptic harmonic residual failed");
    for(int k=0;k<20;++k)if(tsp_math::jet_active(k)){auto e=tsp_math::ex[k];a[k]*=Real(prism.density)*tsp_native_math::pow(Real(length),2-e.x-e.y-e.z);}
    if(diagnostics)*diagnostics=Diagnostics{};return a;
}
}
inline Gravity GravityCartesian(const Prism&prism,const Vec&point,const Options&options={},Diagnostics*diagnostics=nullptr) {
    return detail::rounded(detail::gravity_tensor(prism,point,options,diagnostics));
}
inline Gravity GravityLocalCartesian(const Prism&prism,const Vec&point,const Options&options={},Diagnostics*diagnostics=nullptr) {
    const auto field=detail::gravity_tensor(prism,point,options,diagnostics);
    return detail::rounded(detail::rotate(field,detail::local_axes(detail::convert(point))));
}
inline Magnetic MagneticCartesian(Prism prism,const Vec&point,const Vec&magnetization,const Options&options={},Diagnostics*diagnostics=nullptr) {
    detail::validate_point(magnetization);prism.density=1;
    return detail::rounded(detail::contract(detail::gravity_tensor(prism,point,options,diagnostics),detail::convert(magnetization)));
}
inline std::array<Vec,3> LocalAxes(double latitude,double longitude) {
    const auto axes=detail::local_axes(latitude,longitude);
    return {{detail::rounded(axes[0]),detail::rounded(axes[1]),detail::rounded(axes[2])}};
}
inline Gravity RotateGravity(const Gravity&a,const std::array<Vec,3>&axes) {
    detail::WideGravity field;for(int k=0;k<20;++k)field[k]=detail::Real(a[k]);
    return detail::rounded(detail::rotate(field,{{detail::convert(axes[0]),detail::convert(axes[1]),detail::convert(axes[2])}}));
}
namespace work_detail {
inline constexpr int GravityOrder[20]={0,1,2,3,4,5,7,6,8,9,10,11,12,14,13,16,17,15,18,19};
using OutputPointers=std::array<double*,20>;
using MagneticOutputPointers=std::array<double*,10>;
inline void require(const double*p,int n) {if(n&&p==nullptr)throw std::invalid_argument("null array pointer");}
inline void counts(int m,int n) {if(m<0||n<0)throw std::invalid_argument("negative model/coordinate count");}
template<std::size_t N> inline void check_outputs(const std::array<double*,N>&p,int n) {for(auto v:p)require(v,n);}
inline Vec observation(double lat,double lon,double r) {if(!std::isfinite(r)||r<0)throw std::invalid_argument("invalid observation radius");return detail::scale(Direction(lat,lon),r);}
inline detail::RV source_magnetization_wide(const Prism&p,double mx,double my,double mz) {
    detail::RV center{};for(const auto&v:p.vertices){const auto u=detail::runit(detail::convert(v));for(int k=0;k<3;++k)center[k]+=u[k];}
    const auto axes=detail::local_axes(center);detail::RV out{};
    for(int k=0;k<3;++k)out[k]=mx*axes[0][k]+my*axes[1][k]+mz*axes[2][k];return out;
}
inline Vec source_magnetization(const Prism&p,double mx,double my,double mz) {return detail::rounded(source_magnetization_wide(p,mx,my,mz));}
template<std::size_t N> inline auto coordinate(int models,const std::array<const double*,N>&m,double lat,double lon,double r) {
    static_assert(N==9||N==11);constexpr int fields=N==9?20:10;
    const auto point=observation(lat,lon,r);const auto axes=detail::local_axes(lat,lon);std::array<detail::Real,fields>sum{};
    for(int j=0;j<models;++j) {
        Prism prism{{Direction(m[0][j],m[3][j]),Direction(m[1][j],m[4][j]),Direction(m[2][j],m[5][j])},m[6][j],m[7][j]};
        if constexpr(N==9)prism.density=m[8][j];
        const auto field=detail::rotate(detail::gravity_tensor(prism,point,{},nullptr),axes);
        if constexpr(N==9){for(int k=0;k<fields;++k)if(tsp_math::jet_active(k))sum[k]+=field[k];}
        else {
            const auto mg=source_magnetization_wide(prism,m[8][j],m[9][j],m[10][j]);detail::RV ml{};
            for(int k=0;k<3;++k)ml[k]=detail::rdot(axes[k],mg);
            const auto value=detail::contract(field,ml);for(int k=0;k<fields;++k)sum[k]+=value[k];
        }
    }
    return detail::rounded(sum);
}
struct Input {
    int models=0,points=0;std::array<std::vector<double>,11>model;std::array<std::vector<double>,3>coordinates;
    explicit Input(const char*path,bool magnetic) {
        std::ifstream in(path,std::ios::binary|std::ios::ate);if(!in)throw std::runtime_error(std::string("cannot open ")+path);
        const auto size=in.tellg();if(size<8)throw std::runtime_error("truncated parameter file");in.seekg(0);
        auto count=[&](int&v){std::int32_t n;in.read(reinterpret_cast<char*>(&n),4);if(!in||n<0)throw std::runtime_error("invalid count");v=n;};
        count(models);const int fields=magnetic?11:9;const std::uint64_t bytes=std::uint64_t(fields)*models*8;
        if(bytes>static_cast<std::uint64_t>(size)-8)throw std::runtime_error("model count exceeds file size");
        auto read=[&](std::vector<double>&v,int n){v.resize(n);if(n)in.read(reinterpret_cast<char*>(v.data()),std::streamsize(n)*8);
            if(!in)throw std::runtime_error("truncated parameter array");for(double x:v)if(!std::isfinite(x))throw std::invalid_argument("nonfinite parameter");};
        for(int i=0;i<fields;++i)read(model[i],models);count(points);
        if(8+bytes+std::uint64_t(points)*24!=static_cast<std::uint64_t>(size))throw std::runtime_error("parameter size mismatch");
        for(auto&v:coordinates)read(v,points);
    }
};
template<std::size_t N> inline void write_result(const char*path,const Input&input,const std::array<const double*,N>&values) {
    std::ofstream out(path,std::ios::binary|std::ios::trunc);if(!out)throw std::runtime_error("cannot create result");
    const std::int32_t n=input.points;out.write(reinterpret_cast<const char*>(&n),4);const auto bytes=std::streamsize(n)*8;
    for(const auto&v:input.coordinates)if(bytes)out.write(reinterpret_cast<const char*>(v.data()),bytes);
    std::vector<double>empty;
    if(std::find(values.begin(),values.end(),nullptr)!=values.end())empty.assign(n,std::numeric_limits<double>::quiet_NaN());
    for(auto v:values)if(bytes)out.write(reinterpret_cast<const char*>(v?v:empty.data()),bytes);
    out.close();if(!out)throw std::runtime_error("result write failed");
}
}

inline void TSP_GravityEstimate(int NumDensityModel,const double*Fai1,const double*Fai2,const double*Fai3,
    const double*Lamda1,const double*Lamda2,const double*Lamda3,const double*R1,const double*R2,const double*rho,
    int NumCoordinates,const double*Latitude,const double*Longitude,const double*Radius,
    double*V,double*Vx,double*Vy,double*Vz,double*Vxx,double*Vxy,double*Vyy,double*Vzx,double*Vzy,double*Vzz,
    double*Vxxx,double*Vxxy,double*Vxxz,double*Vxyz,double*Vyyx,double*Vyyy,double*Vyyz,double*Vzzx,double*Vzzy,double*Vzzz) {
    using namespace work_detail;counts(NumDensityModel,NumCoordinates);
    const std::array<const double*,9>m{Fai1,Fai2,Fai3,Lamda1,Lamda2,Lamda3,R1,R2,rho};for(auto v:m)require(v,NumDensityModel);
    for(auto v:{Latitude,Longitude,Radius})require(v,NumCoordinates);
    const OutputPointers out{V,Vx,Vy,Vz,Vxx,Vxy,Vzx,Vyy,Vzy,Vzz,Vxxx,Vxxy,Vxxz,Vyyx,Vxyz,Vzzx,Vyyy,Vyyz,Vzzy,Vzzz};check_outputs(out,NumCoordinates);
    for(int i=0;i<NumCoordinates;++i){auto a=coordinate(NumDensityModel,m,Latitude[i],Longitude[i],Radius[i]);for(int k=0;k<20;++k)out[k][i]=a[k];}
}
inline void TSP_MagneticEstimate(int NumMagnetizationModel,const double*Fai1,const double*Fai2,const double*Fai3,
    const double*Lamda1,const double*Lamda2,const double*Lamda3,const double*R1,const double*R2,const double*Mx,const double*My,const double*Mz,
    int NumCoordinates,const double*Latitude,const double*Longitude,const double*Radius,
    double*V,double*Vx,double*Vy,double*Vz,double*Vxx,double*Vxy,double*Vyy,double*Vzx,double*Vzy,double*Vzz) {
    using namespace work_detail;counts(NumMagnetizationModel,NumCoordinates);
    const std::array<const double*,11>m{Fai1,Fai2,Fai3,Lamda1,Lamda2,Lamda3,R1,R2,Mx,My,Mz};for(auto v:m)require(v,NumMagnetizationModel);
    for(auto v:{Latitude,Longitude,Radius})require(v,NumCoordinates);
    const MagneticOutputPointers out{V,Vx,Vy,Vz,Vxx,Vxy,Vyy,Vzx,Vzy,Vzz};check_outputs(out,NumCoordinates);
    for(int i=0;i<NumCoordinates;++i){auto a=coordinate(NumMagnetizationModel,m,Latitude[i],Longitude[i],Radius[i]);for(int k=0;k<10;++k)out[k][i]=a[k];}
}
inline void TSP_VzEstimate(int models,const double*f1,const double*f2,const double*f3,const double*l1,const double*l2,const double*l3,
    const double*r1,const double*r2,const double*rho,int n,const double*lat,const double*lon,const double*r,double*vz) {
    tsp_math::JetScope scope(15);
    work_detail::counts(models,n);const std::array<const double*,9>m{f1,f2,f3,l1,l2,l3,r1,r2,rho};for(auto p:m)work_detail::require(p,models);
    for(auto p:std::array<const double*,4>{lat,lon,r,vz})work_detail::require(p,n);
    for(int i=0;i<n;++i)vz[i]=work_detail::coordinate(models,m,lat[i],lon[i],r[i])[3];
}
inline void TSP_MagneticVectorEstimate(int models,const double*f1,const double*f2,const double*f3,const double*l1,const double*l2,const double*l3,
    const double*r1,const double*r2,const double*mx,const double*my,const double*mz,int n,const double*lat,const double*lon,const double*r,double*bx,double*by,double*bz) {
    tsp_math::JetScope scope(1023);
    work_detail::counts(models,n);const std::array<const double*,11>m{f1,f2,f3,l1,l2,l3,r1,r2,mx,my,mz};for(auto p:m)work_detail::require(p,models);
    for(auto p:std::array<const double*,6>{lat,lon,r,bx,by,bz})work_detail::require(p,n);
    for(int i=0;i<n;++i){auto v=work_detail::coordinate(models,m,lat[i],lon[i],r[i]);bx[i]=v[1];by[i]=v[2];bz[i]=v[3];}
}
}
#endif

#ifndef TSP_PROGRESS_INCLUDED
#define TSP_PROGRESS_INCLUDED
#include <atomic>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <mutex>
#ifdef _WIN32
#include <io.h>
#else
#include <unistd.h>
#endif
#ifdef TSP_MPI
#include <mpi.h>
#endif
#include <cstdarg>
#ifdef _OPENMP
#include <omp.h>
#endif
namespace tsp_detail {
inline void Stage(const char*label,const char*phase,const char*message) {
    std::printf("%s | %s %s\n",label,phase,message);std::fflush(stdout);
}
struct ParameterReport {
    const char*label;
    ParameterReport(const char*name,int kind,int mode,int ranks):label(name) {
        std::printf("\n%s +------------------------------------------------------------------\n%s | RUN PARAMETERS\n",label,label);
        Line("Field / outputs","%s",kind==0?"Gravity | all 20 components":kind==1?"Magnetics | all 10 components":kind==2?"Gravity | Vz only":"Magnetics | vector only");
        int threads=1;
#ifdef _OPENMP
        if(mode==1)threads=omp_get_max_threads();
#endif
        Line("Execution","%s | MPI ranks=%d | thread limit/rank=%d",mode==2?"MPI":mode==1?"OpenMP":"Serial",ranks,threads);
    }
    void Line(const char*key,const char*format,...) const {
        std::printf("%s | %-18s : ",label,key);va_list args;va_start(args,format);std::vprintf(format,args);va_end(args);std::putchar('\n');
    }
    void Range(const char*key,const double*data,std::size_t count,std::size_t stride=1) const {
        if(!count){Line(key,"n/a (empty)");return;}
        double lo=data[0],hi=lo;
        for(std::size_t k=1;k<count;++k){lo=std::min(lo,data[k*stride]);hi=std::max(hi,data[k*stride]);}
        Line(key,"%.12g .. %.12g",lo,hi);
    }
    void End() const {
        std::printf("%s +------------------------------------------------------------------\n",label);std::fflush(stdout);
    }
};

#ifdef TSP_MPI
inline void MPICheck(int error,const char*operation) {
    if(error==MPI_SUCCESS)return;
    char message[MPI_MAX_ERROR_STRING];int length=0;MPI_Error_string(error,message,&length);
    throw std::runtime_error(std::string(operation)+": "+std::string(message,length));
}
// A blocking barrier turns one missing rank into an unbounded wait. Guard the
// arrival with a timeout so a stall aborts the job instead of holding the nodes.
inline double SyncTimeout() {
    static const double seconds=[]() {const char*p=std::getenv("TSP_MPI_SYNC_TIMEOUT");return p?std::atof(p):1800.;}();
    return seconds;
}
inline void SyncBarrier(const char* where) {
    MPI_Request request;MPICheck(MPI_Ibarrier(MPI_COMM_WORLD,&request),where);
    const double limit=SyncTimeout(),started=MPI_Wtime();
    for(;;) {int complete=0;MPICheck(MPI_Test(&request,&complete,MPI_STATUS_IGNORE),where);
        if(complete)return;
        const double waited=MPI_Wtime()-started;
        if(limit>0&&waited>limit) {int rank=0;MPI_Comm_rank(MPI_COMM_WORLD,&rank);
            std::fprintf(stderr,"%s | rank %d: barrier not reached after %.0f s (TSP_MPI_SYNC_TIMEOUT=%.0f)\n",
                where,rank,waited,limit);std::fflush(stderr);MPI_Abort(MPI_COMM_WORLD,3);}}
}
inline std::pair<int,int> OwnedRange(int total,int rank,int ranks) {
    const int base=total/ranks,extra=total%ranks,begin=rank*base+std::min(rank,extra);
    return {begin,begin+base+(rank<extra)};
}
inline void RankSum(std::uint64_t*values,int count) {
    int rank=0,ranks=1;MPI_Comm_rank(MPI_COMM_WORLD,&rank);MPI_Comm_size(MPI_COMM_WORLD,&ranks);
    std::vector<int>counts(ranks,count),offsets(ranks);
    for(int r=0;r<ranks;++r)offsets[r]=r*count;
    std::vector<std::uint64_t>buffer(rank==0?std::size_t(count)*ranks:0);
    SyncBarrier("diagnostic barrier");
    MPICheck(MPI_Gatherv(values,count,MPI_UINT64_T,rank==0?buffer.data():nullptr,counts.data(),offsets.data(),MPI_UINT64_T,0,MPI_COMM_WORLD),"diagnostic Gatherv");
    if(rank==0)for(int k=0;k<count;++k){values[k]=0;for(int r=0;r<ranks;++r)values[k]+=buffer[r*count+k];}
    MPICheck(MPI_Bcast(values,count,MPI_UINT64_T,0,MPI_COMM_WORLD),"diagnostic Bcast");
}
class IndexedGather {
    int total_,local_,rank_=0,ranks_=1;
    const char*label_;
    std::vector<int>counts_,offsets_,order_;
    std::vector<double>buffer_;
  public:
    IndexedGather(int total,int begin,int end,const char*label):total_(total),local_(end-begin),label_(label) {
        MPI_Comm_rank(MPI_COMM_WORLD,&rank_);MPI_Comm_size(MPI_COMM_WORLD,&ranks_);
        if(rank_==0)Stage(label_,"COMM","validating owned observation indices...");
        counts_.resize(ranks_);offsets_.resize(ranks_);
        MPICheck(MPI_Gather(&local_,1,MPI_INT,counts_.data(),1,MPI_INT,0,MPI_COMM_WORLD),"owned-count Gather");
        int valid=total>=0;
        if(rank_==0){std::int64_t count=0;
            for(int r=0;r<ranks_;++r){if(count>std::numeric_limits<int>::max()||counts_[r]<0)valid=0;
                offsets_[r]=count<=std::numeric_limits<int>::max()?int(count):0;count+=counts_[r];}
            if(count!=total)valid=0;}
        MPICheck(MPI_Bcast(&valid,1,MPI_INT,0,MPI_COMM_WORLD),"owned-count validity Bcast");
        if(!valid)throw std::runtime_error("MPI ownership does not cover all observations");
        std::vector<int>rows(local_);for(int k=0;k<local_;++k)rows[k]=begin+k;
        if(rank_==0)order_.resize(total_);buffer_.resize(total_);
        int dummy=0;
        MPICheck(MPI_Gatherv(rows.empty()?&dummy:rows.data(),local_,MPI_INT,
            rank_==0&&!order_.empty()?order_.data():&dummy,counts_.data(),offsets_.data(),MPI_INT,0,MPI_COMM_WORLD),"owned-index Gatherv");
        if(rank_==0){std::vector<unsigned char>seen(total_,0);
            for(int row:order_){if(row<0||row>=total_||seen[row]){valid=0;break;}seen[row]=1;}}
        MPICheck(MPI_Bcast(&valid,1,MPI_INT,0,MPI_COMM_WORLD),"owned-index validity Bcast");
        if(!valid)throw std::runtime_error("MPI ownership contains missing or duplicate indices");
        if(rank_==0)Stage(label_,"COMM","ownership verified; indexed gather/broadcast ready");
    }
    void Collect(const std::vector<double>&local,double*global,int component,int fields) {
        if(local.size()!=std::size_t(local_)||(rank_==0&&total_&&!global))throw std::runtime_error("invalid MPI gather dimensions");
        if(rank_==0){std::printf("%s | GATHER component %d/%d BEGIN barrier/Gatherv/Bcast\n",label_,component,fields);std::fflush(stdout);}
        const double start=MPI_Wtime();double dummy=0;
        SyncBarrier("field gather barrier");
        MPICheck(MPI_Gatherv(local.empty()?&dummy:local.data(),local_,MPI_DOUBLE,
            rank_==0&&!buffer_.empty()?buffer_.data():&dummy,counts_.data(),offsets_.data(),MPI_DOUBLE,0,MPI_COMM_WORLD),"owned-field Gatherv");
        if(rank_==0)for(int k=0;k<total_;++k)global[order_[k]]=buffer_[k];
        MPICheck(MPI_Bcast(total_?(rank_==0?global:buffer_.data()):&dummy,total_,MPI_DOUBLE,0,MPI_COMM_WORLD),"owned-field Bcast");
        SyncBarrier("field gather completion barrier");
        if(rank_==0){std::printf("%s | GATHER component %d/%d COMPLETE elapsed_seconds=%.6f\n",label_,component,fields,MPI_Wtime()-start);std::fflush(stdout);}
    }
};
#endif

class Progress {
    using Clock=std::chrono::steady_clock;
    const char*label_;
    std::uint64_t total_;
    std::atomic<std::uint64_t>completed_{0};
    std::mutex mutex_;
    Clock::time_point started_=Clock::now(),displayed_=started_;
    int rank_=0;
    bool terminal_=false,distributed_=false;
#ifdef TSP_MPI
    int rounds_=0,tick_=0,interval_=1;
    std::uint64_t global_=0;
    std::vector<int>counts_,offsets_,gathered_;
    void Sync() {
        const int local=static_cast<int>(completed_.load(std::memory_order_relaxed));int global=0;
        // Mirrors R-SCHA v2.0mpi/MPIParallel.h, which deliberately avoids
        // MPI_Allreduce and uses barrier + Gatherv + root sum + Bcast instead.
        SyncBarrier("progress barrier");
        MPICheck(MPI_Gatherv(&local,1,MPI_INT,rank_==0?gathered_.data():nullptr,
            counts_.data(),offsets_.data(),MPI_INT,0,MPI_COMM_WORLD),"progress Gatherv");
        if(rank_==0){long long sum=0;for(int v:gathered_)sum+=v;global=static_cast<int>(sum);}
        MPICheck(MPI_Bcast(&global,1,MPI_INT,0,MPI_COMM_WORLD),"progress Bcast");
        global_=static_cast<std::uint64_t>(global);
    }
#endif
    static void Time(double seconds,char*text) {
        const auto t=static_cast<unsigned long long>(seconds);
        std::snprintf(text,32,"%02llu:%02llu:%02llu",t/3600,(t/60)%60,t%60);
    }
    void Print(Clock::time_point now,bool final=false,bool initial=false) {
        if(rank_!=0)return;
        if(!final&&!initial&&std::chrono::duration<double>(now-displayed_).count()<(terminal_?1:10))return;
        std::uint64_t done=completed_.load(std::memory_order_relaxed);
#ifdef TSP_MPI
        if(distributed_)done=global_;
#endif
        if(!final&&!initial&&done==total_)return;
        const double seconds=std::chrono::duration<double>(now-started_).count();
        char elapsed[32],remaining[32]="--:--:--";Time(seconds,elapsed);
        if(done&&(!final||done==total_))Time(std::ceil(seconds*(total_-done)/done),remaining);else if(!total_)Time(0,remaining);
        const double percent=total_?(done<total_?std::min(99.9,100.*done/total_):100.):100.;const int filled=static_cast<int>(percent/5);
        std::printf("%s%s [%.*s%.*s] %6.1f%% %llu/%llu elapsed %s ETA %s%s%s",
            terminal_?"\r":"",label_,filled,"####################",20-filled,"--------------------",percent,
            static_cast<unsigned long long>(done),static_cast<unsigned long long>(total_),elapsed,remaining,
            final&&done<total_?" stopped":"",terminal_?"   ":"");
        if(!terminal_||final)std::putchar('\n');std::fflush(stdout);displayed_=now;
    }
  public:
    Progress(const char*label,int total,bool distributed=false):label_(label),total_(total),distributed_(distributed) {
#ifdef TSP_MPI
        if(distributed_){int ranks=1;MPI_Comm_rank(MPI_COMM_WORLD,&rank_);MPI_Comm_size(MPI_COMM_WORLD,&ranks);
            rounds_=total/ranks+(total%ranks!=0);interval_=std::max(1,(rounds_-1)/25+1);
            counts_.assign(std::size_t(ranks),1);offsets_.resize(std::size_t(ranks));
            for(int i=0;i<ranks;++i)offsets_[std::size_t(i)]=i;
            if(rank_==0)gathered_.resize(std::size_t(ranks));}
#else
        (void)distributed_;
#endif
        if(rank_==0) {
#ifdef _WIN32
            terminal_=_isatty(_fileno(stdout))!=0;
#else
            terminal_=::isatty(::fileno(stdout))!=0;
#endif
            if(std::getenv("SLURM_JOB_ID"))terminal_=false;
            if(total_)Print(Clock::now(),false,true);
        }
    }
    void Step(bool completed=true) {
        if(completed)++completed_;
#ifdef TSP_MPI
        if(distributed_){++tick_;if(tick_%interval_==0||tick_==rounds_){Sync();Print(Clock::now());}return;}
#endif
        std::unique_lock<std::mutex>lock(mutex_,std::try_to_lock);if(lock.owns_lock())Print(Clock::now());
    }
    void Finish() {
        std::lock_guard<std::mutex>lock(mutex_);
#ifdef TSP_MPI
        if(distributed_)Sync();
#endif
        Print(Clock::now(),true);
        if(distributed_&&rank_==0)Stage(label_,"SYNC","all ranks completed observation rounds");
    }
};
}
#endif

#ifndef TSP_MMSI_INCLUDED
#define TSP_MMSI_INCLUDED
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace TSP_MMSI { namespace radial {
inline double Expm1(double x) {return std::expm1(x);}
inline double LogRatio(const double a, const double b) {
    const double r = (a - b) / b;
    if (std::abs(r) < 0.25) {
        return std::log1p(r);
    }
    return std::log(a / b);
}

inline double StablePowerIntegral(const int alpha, const double q1, const double q2) {
    if (q1 <= 0.0 || q2 <= 0.0) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    if (q1 == q2) {
        return 0.0;
    }

    const double logRatio = LogRatio(q2, q1);

    if (alpha == -1) {
        return logRatio;
    }

    const int betaInt = alpha + 1;
    const double beta = static_cast<double>(betaInt);
    return std::pow(q1, beta) * std::expm1(beta * logRatio) / beta;
}

inline double BaseMomentSmallB(const double u1, const double u2, const double b, const int n, const int k) {
    if (u1 * u2 <= 0.0) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    const double sigma = (u1 > 0.0) ? 1.0 : -1.0;
    const double q1 = std::min(std::abs(u1), std::abs(u2));
    const double q2 = std::max(std::abs(u1), std::abs(u2));

    const double eta=b/q1,logRatio=LogRatio(q2,q1);
    double coefficient=1,sum=0;
    for(int j=0;j<80;++j) {
        const int beta=k-n-2*j+1;
        const double integral=beta?std::expm1(beta*logRatio)/beta:logRatio;
        const double term=coefficient*integral;sum+=term;
        if(j>2&&std::abs(term)<2e-17*std::abs(sum))break;
        coefficient*=-(.5*n+j)/(j+1)*eta*eta;
    }
    return (k%2?sigma:1)*std::pow(q1,k-n+1)*sum;
}

inline double StableSqrtDifference(const double u1, const double u2, const double L1, const double L2) {
    return (u2 - u1) * (u2 + u1) / (L2 + L1);
}

inline double StableUOverPowerDifference(const double u1, const double u2, const double L1, const double L2,
                                         const int power, bool& reliable) {
    if (u1 == 0.0 || u2 == 0.0 || u1 * u2 < 0.0) {
        return u2 / std::pow(L2, power) - u1 / std::pow(L1, power);
    }

    const double logU = LogRatio(std::abs(u2), std::abs(u1));
    const double dL = StableSqrtDifference(u1, u2, L1, L2);
    const double logL = std::abs(dL)<0.25*L1?std::log1p(dL/L1):std::log(L2/L1);
    const double deltaLog = logU - static_cast<double>(power) * logL;
    const double logScale = std::abs(logU) + static_cast<double>(power) * std::abs(logL);
    if (std::abs(deltaLog) <= 1.0e-10 * logScale) {
        reliable = false;
    }

    const double value1 = u1 / std::pow(L1, power);
    return value1 * std::expm1(deltaLog);
}

inline double StableInversePowerDifference(const double u1, const double u2, const double L1, const double L2,
                                           const int power) {
    const double dL = StableSqrtDifference(u1, u2, L1, L2);
    const double logL = std::abs(dL)<0.25*L1?std::log1p(dL/L1):std::log(L2/L1);
    return std::pow(L1, -power) * std::expm1(-static_cast<double>(power) * logL);
}

inline void StableBaseMomentDiffTable(const double u1, const double u2, const double b, const int n, const int maxK,
                                      double* output, bool& reliable) {
    reliable = true;
    const double b2 = b * b;
    const double L1 = std::hypot(u1, b);
    const double L2 = std::hypot(u2, b);
    double D[4][7] = {};

    if (u1 > 0.0 && u2 > 0.0) {
        D[0][0] = LogRatio(u2 + L2, u1 + L1);
    } else if (u1 < 0.0 && u2 < 0.0) {
        D[0][0] = LogRatio(L1 - u1, L2 - u2);
    } else {
        D[0][0] = std::asinh(u2 / b) - std::asinh(u1 / b);
    }

    if (maxK >= 1) {
        D[0][1] = StableSqrtDifference(u1, u2, L1, L2);
    }
    if (maxK >= 2) {
        double deltaUL;
        if (u1 * u2 > 0.0) {
            const double denominator = u2 * L2 + u1 * L1;
            if (denominator == 0.0) {
                reliable = false;
                deltaUL = u2 * L2 - u1 * L1;
            } else {
                deltaUL = (u2 - u1) * (u2 + u1) * (u2 * u2 + u1 * u1 + b2) / denominator;
            }
        } else {
            deltaUL = u2 * L2 - u1 * L1;
        }
        D[0][2] = 0.5 * (deltaUL - b2 * D[0][0]);
    }

    const int targetIndex = (n - 1) / 2;
    for (int index = 1; index <= targetIndex; ++index) {
        const int currentN = 2 * index + 1;
        const int power = currentN - 2;
        const double inv = 1.0 / (static_cast<double>(power) * b2);
        bool endpointReliable = true;
        const double endpoint0 = StableUOverPowerDifference(u1, u2, L1, L2, power, endpointReliable);
        const double term0 = inv * endpoint0;
        const double term1 = static_cast<double>(currentN - 3) * inv * D[index - 1][0];
        D[index][0] = term0 + term1;
        reliable = reliable && endpointReliable;
        if (currentN > 3 && std::abs(D[index][0]) <= 1.0e-10 * (std::abs(term0) + std::abs(term1))) {
            reliable = false;
        }

        if (maxK >= 1) {
            D[index][1] = -StableInversePowerDifference(u1, u2, L1, L2, power) / static_cast<double>(power);
        }

        const int supportedMaxK = currentN == 3 ? 3 : (currentN == 5 ? 4 : 6);
        for (int currentK = 2; currentK <= std::min(maxK, supportedMaxK); ++currentK) {
            const double first = D[index - 1][currentK - 2];
            const double second = b2 * D[index][currentK - 2];
            D[index][currentK] = std::fma(-b2, D[index][currentK - 2], first);
            if (std::abs(D[index][currentK]) <= 1.0e-10 * (std::abs(first) + std::abs(second))) {
                reliable = false;
            }
        }
    }

    for (int currentK = 0; currentK <= maxK; ++currentK) {
        output[currentK] = D[targetIndex][currentK];
        reliable = reliable && std::isfinite(output[currentK]);
    }
}

inline void ExtendedMomentTable(const dd_real&a,const dd_real&b,const dd_real&h,dd_real output[4][7]) {
    for(int endpoint=0;endpoint<2;++endpoint) {
        const dd_real u=endpoint?b:a,h2=h*h,L=sqrt(u*u+h2),z=u/L,z2=z*z,z3=z*z2,z5=z3*z2;
        const dd_real A=u>=0?log(u+L)-log(h):log(h)-log(L-u);
        const dd_real inv=1/L,inv3=inv*inv*inv,inv5=inv3*inv*inv;
        dd_real v[4][7]{};
        v[0][0]=A;v[0][1]=L;v[0][2]=(u*L-h2*A)/2;
        v[1][0]=z/h2;v[1][1]=-inv;v[1][2]=A-z;v[1][3]=L+h2*inv;
        v[2][0]=(z-z3/3)/(h2*h2);v[2][1]=-inv3/3;v[2][2]=z3/(3*h2);
        v[2][3]=-inv+h2*inv3/3;v[2][4]=A-z-z3/3;
        v[3][0]=(z-2*z3/3+z5/5)/(h2*h2*h2);v[3][1]=-inv5/5;
        v[3][2]=(z3/3-z5/5)/(h2*h2);v[3][3]=-inv3/3+h2*inv5/5;
        v[3][4]=z5/(5*h2);v[3][5]=-inv+2*h2*inv3/3-h2*h2*inv5/5;
        v[3][6]=A-z-z3/3-z5/5;
        for(int n=0;n<4;++n)for(int k=0;k<7;++k)output[n][k]+=endpoint?v[n][k]:-v[n][k];
    }
}
inline dd_real LogRatio(const dd_real&a,const dd_real&b) {
    const dd_real r=(a-b)/b;
    if(abs(r)>=.25)return log(a/b);
    const dd_real t=r/(2+r),t2=t*t;dd_real sum=dd_real(1)/41;
    for(int k=19;k>=0;--k)sum=dd_real(1)/(2*k+1)+t2*sum;
    return 2*t*sum;
}
inline dd_real Expm1(const dd_real&x) {
    if(abs(x)>=.001)return exp(x)-1;
    dd_real sum=1,term=1;
    for(int k=2;k<=14;++k){term*=x/k;sum+=term;}
    return x*sum;
}
inline dd_real StablePowerIntegral(int alpha,const dd_real&a,const dd_real&b) {
    const dd_real logarithm=LogRatio(b,a);
    if(alpha==-1)return logarithm;
    return pow(a,alpha+1)*Expm1((alpha+1)*logarithm)/(alpha+1);
}
inline void BaseMomentDiffTable(const dd_real&a,const dd_real&b,const dd_real&h,int n,int maxK,dd_real*output) {
    if(a==b){std::fill(output,output+maxK+1,dd_real(0));return;}
    if(a*b>0) {
        const dd_real low=std::min(abs(a),abs(b)),high=std::max(abs(a),abs(b)),ratio=h/low;
        if(ratio<.5) {
            const dd_real logarithm=LogRatio(high,low);
            for(int k=0;k<=maxK;++k) {
                dd_real coefficient=1,sum=0;
                for(int j=0;j<80;++j) {
                    const int beta=k-n-2*j+1;
                    const dd_real term=coefficient*(beta?Expm1(beta*logarithm)/beta:logarithm);sum+=term;
                    if(j>2&&abs(term)<dd_real("1e-31")*abs(sum))break;
                    coefficient*=-(dd_real(n)/2+j)/(j+1)*ratio*ratio;
                }
                output[k]=(b>a?1:-1)*(k%2&&a<0?-1:1)*pow(low,k-n+1)*sum;
            }
            return;
        }
    }
    struct Cache {bool valid=false;dd_real a,b,h,values[4][7]{};};
    static thread_local Cache cache;
    if(!cache.valid||cache.a!=a||cache.b!=b||cache.h!=h) {
        cache=Cache{};cache.valid=true;cache.a=a;cache.b=b;cache.h=h;
        ExtendedMomentTable(a,b,h,cache.values);
    }
    std::copy_n(cache.values[(n-1)/2],maxK+1,output);
}

inline void BaseMomentDiffTable(const double u1, const double u2, const double bIn, const int n, const int maxK,
                                double* output) {
    if (u1 == u2) {
        std::fill(output, output + maxK + 1, 0.0);
        return;
    }

    double a = u1;
    double bnd = u2;
    double orderSign = 1.0;
    if (bnd < a) {
        std::swap(a, bnd);
        orderSign = -1.0;
    }

    double b = std::abs(bIn);
    if (b == 0.0) {
        b = std::numeric_limits<double>::min();
    }

    const bool sameSide = (a * bnd > 0.0);
    if (sameSide) {
        const double qmin = std::min(std::abs(a), std::abs(bnd));
        const double eta = b / std::max(qmin, std::numeric_limits<double>::min());
        if (eta < 0.5) {
            for (int k = 0; k <= maxK; ++k) {
                output[k] = orderSign * BaseMomentSmallB(a, bnd, b, n, k);
            }
            return;
        }
    }

    bool reliable = true;
    StableBaseMomentDiffTable(a, bnd, b, n, maxK, output, reliable);
    if (reliable) {
        for (int k = 0; k <= maxK; ++k) {
            output[k] *= orderSign;
        }
        return;
    }

    dd_real extended[7]{};
    BaseMomentDiffTable(dd_real(a),dd_real(bnd),dd_real(b),n,maxK,extended);
    for(int k=0;k<=maxK;++k)output[k]=orderSign*::to_double(extended[k]);
}

}}
#include <array>
#include <cstdint>
#include <stdexcept>
#include <type_traits>
#include <vector>

namespace TSP_MMSI {
constexpr double PI = 3.141592653589793238462643383279502884;
using Vec = std::array<double, 3>;
using Gravity = std::array<double, 20>;
using Magnetic = std::array<double, 10>;
namespace detail {
constexpr int Axes[20][3] = {
    {0,0,0},{1,0,0},{0,1,0},{0,0,1},{2,0,0},{1,1,0},{1,0,1},{0,2,0},{0,1,1},{0,0,2},
    {3,0,0},{2,1,0},{2,0,1},{1,2,0},{1,1,1},{1,0,2},{0,3,0},{0,2,1},{0,1,2},{0,0,3}};
constexpr int GravityOrder[20] = {0,1,2,3,4,5,7,6,8,9,10,11,12,14,13,16,17,15,18,19};
inline int Degree(int k) { return Axes[k][0]+Axes[k][1]+Axes[k][2]; }
inline int Component(int x,int y,int z) {
    static constexpr signed char index[4][4][4]={
        {{0,3,9,19},{2,8,18,-1},{7,17,-1,-1},{16,-1,-1,-1}},
        {{1,6,15,-1},{5,14,-1,-1},{13,-1,-1,-1},{-1,-1,-1,-1}},
        {{4,12,-1,-1},{11,-1,-1,-1},{-1,-1,-1,-1},{-1,-1,-1,-1}},
        {{10,-1,-1,-1},{-1,-1,-1,-1},{-1,-1,-1,-1},{-1,-1,-1,-1}}};
    if(unsigned(x)>3||unsigned(y)>3||unsigned(z)>3||index[x][y][z]<0)
        throw std::logic_error("invalid tensor component");
    return index[x][y][z];
}
template<class T> inline T Dot(const std::array<T,3>&a,const std::array<T,3>&b) { return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]; }
template<class T> inline std::array<T,3> Cross(const std::array<T,3>&a,const std::array<T,3>&b) {
    return {a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]};
}
inline Vec Add(Vec a,const Vec&b) { for(int i=0;i<3;++i)a[i]+=b[i];return a; }
inline Vec Subtract(Vec a,const Vec&b) { for(int i=0;i<3;++i)a[i]-=b[i];return a; }
inline Vec Scale(Vec a,double s) { for(auto&v:a)v*=s;return a; }
inline double Norm(const Vec&a) { const double scale=std::max({std::abs(a[0]),std::abs(a[1]),std::abs(a[2])});
    if(scale==0)return 0;
    const double x=a[0]/scale,y=a[1]/scale,z=a[2]/scale;
    return scale*std::sqrt(x*x+y*y+z*z); }
inline Vec Unit(const Vec&a) {
    const double r=Norm(a);if(!(r>0)||!std::isfinite(r))throw std::invalid_argument("invalid direction");return Scale(a,1/r);
}
inline Vec Direction(double lat,double lon) {
    if(!std::isfinite(lat)||!std::isfinite(lon)||std::abs(lat)>90)throw std::invalid_argument("invalid coordinates");
    const double p=lat*PI/180,l=lon*PI/180;
    return {std::cos(p)*std::cos(l),std::cos(p)*std::sin(l),std::sin(p)};
}
inline std::array<Vec,3> LocalAxes(double lat,double lon) {
    const double p=lat*PI/180,l=lon*PI/180;
    return {{{-std::sin(p)*std::cos(l),-std::sin(p)*std::sin(l),std::cos(p)},
             {-std::sin(l),std::cos(l),0},Direction(lat,lon)}};
}
inline void Accumulate(double value,double&sum,double&correction) {
    const double next=sum+value;
    correction+=std::abs(sum)>=std::abs(value)?(sum-next)+value:(value-next)+sum;sum=next;
}
inline double Area(const std::array<Vec,3>&v) {
    return 2*std::atan2(std::abs(Dot(v[0],Cross(v[1],v[2]))),
                       1+Dot(v[0],v[1])+Dot(v[1],v[2])+Dot(v[2],v[0]));
}
template<class T=double> struct Rule {
    std::vector<T> x,w;
    explicit Rule(int n):x(n),w(n) {
        using std::abs;
        for(int i=0;i<(n+1)/2;++i) {
            T z=std::cos(PI*(i+.75)/(n+.5)),dp=0;
            for(int j=0;j<12;++j) {
                T p0=1,p1=z;
                for(int k=2;k<=n;++k){const T next=((2*k-1)*z*p1-(k-1)*p0)/k;p0=p1;p1=next;}
                dp=n*(z*p1-p0)/(z*z-1);const T dz=p1/dp;z-=dz;
                if constexpr(std::is_same<T,double>::value){if(abs(dz)<3e-16)break;}
                else if(abs(dz)<T("1e-31"))break;
            }
            x[i]=(1-z)/2;x[n-1-i]=(1+z)/2;
            w[i]=w[n-1-i]=1/((1-z*z)*dp*dp);
        }
    }
};
template<class T> inline std::array<T,20> Rotate(const std::array<T,20>&a,const std::array<std::array<T,3>,3>&rotation,int components=20) {
    std::array<T,20>out{};out[0]=a[0];
    for(int k=1;k<components;++k) {
        int target[3],degree=0;
        for(int i=0;i<3;++i)for(int j=0;j<Axes[k][i];++j)target[degree++]=i;
        for(int i=0;i<3;++i)for(int j=0;j<(degree>1?3:1);++j)for(int h=0;h<(degree>2?3:1);++h) {
            int ex[3]={0,0,0};++ex[i];T weight=rotation[target[0]][i];
            if(degree>1){++ex[j];weight*=rotation[target[1]][j];}
            if(degree>2){++ex[h];weight*=rotation[target[2]][h];}
            out[k]+=weight*a[Component(ex[0],ex[1],ex[2])];
        }
    }
    return out;
}
inline Magnetic Rotate(const Magnetic&a,const std::array<Vec,3>&r,int components=10) {
    Magnetic out{};out[0]=a[0];
    if(components>1)for(int i=0;i<3;++i)out[i+1]=Dot(r[i],{a[1],a[2],a[3]});
    constexpr int ij[6][2]={{0,0},{0,1},{1,1},{2,0},{2,1},{2,2}};
    const double h[3][3]={{a[4],a[5],a[7]},{a[5],a[6],a[8]},{a[7],a[8],a[9]}};
    for(int k=0;k<components-4;++k)for(int i=0;i<3;++i)for(int j=0;j<3;++j)out[k+4]+=r[ij[k][0]][i]*r[ij[k][1]][j]*h[i][j];
    return out;
}

template<class T> inline std::array<T,20> RadialValues(T a,T b,const std::array<T,3>&u,const T&r,int power,int degree) {
    using std::pow;
    std::array<T,20> out{};if(a==b)return out;
    using std::sqrt;
    if(r==0) {
        if(a==0)throw std::domain_error("observer is inside the source");
        for(int k=0;k<20;++k) {
            const int d=Degree(k);if(d>degree)continue;
            int ax[3],at=0;for(int i=0;i<3;++i)for(int j=0;j<Axes[k][i];++j)ax[at++]=i;
            const T radial=radial::StablePowerIntegral(1+power-d,a,b);
            if(d==0)out[k]=radial;
            if(d==1)out[k]=u[ax[0]]*radial;
            if(d==2)out[k]=(3*u[ax[0]]*u[ax[1]]-(ax[0]==ax[1]))*radial;
            if(d==3){const int i=ax[0],j=ax[1],h=ax[2];out[k]=(15*u[i]*u[j]*u[h]
                -3*((i==j)*u[h]+(i==h)*u[j]+(j==h)*u[i]))*radial;}
        }
        return out;
    }
    const T c=std::min(std::max(u[2],T(-1)),T(1)),sine=sqrt(u[0]*u[0]+u[1]*u[1]);
    if(power==0||power==-3) {
        T base[4][7]{};constexpr int dn[4]={1,3,5,7};
        const T omc=c>=0?sine*sine/(1+c):1-c;
        const T zlo=power==0?(a-r)/r+omc:(r-b)/b+omc;
        const T zhi=zlo+(power==0?(b-a)/r:r*(b-a)/(a*b));
        for(int h=0;h<=degree;++h)radial::BaseMomentDiffTable(zlo,zhi,sine,dn[h],power==0?h+2:2*h,base[h]);
        constexpr double choose[6][6]={{1},{1,1},{1,2,1},{1,3,3,1},{1,4,6,4,1},{1,5,10,10,5,1}};
        T cp[7]={1},sp[4]={1};for(int i=1;i<7;++i)cp[i]=cp[i-1]*c;
        for(int i=1;i<4;++i)sp[i]=sp[i-1]*(-sine*sine);
        const auto m=[&](int h,int s,int q) {
            T result=0;
            if(power==-3) {
                const int exponent=dn[h]+1-s-q;
                for(int j=0;j<=exponent;++j)result+=choose[exponent][j]*cp[exponent-j]*base[h][j+q];
                return q%2?-result:result;
            }
            for(int j=0;j<=s;++j)for(int k=0;k<=q;++k)
                result+=choose[s][j]*choose[q][k]*cp[s-j+k]*sp[q-k]*base[h][j+k];
            return result;
        };
        out[0]=r*r*m(0,2,0);
        const T i32=degree>=1?m(1,2,0):0,i33=degree>=1?m(1,3,0):0,z321=degree>=1?m(1,2,1):0;
        const T i54=degree>=2?m(2,4,0):0,z531=degree>=2?m(2,3,1):0,z522=degree>=2?m(2,2,2):0;
        const T i53=degree>=3?m(2,3,0):0,z521=degree>=3?m(2,2,1):0;
        const T i75=degree>=3?m(3,5,0):0,z741=degree>=3?m(3,4,1):0;
        const T z732=degree>=3?m(3,3,2):0,z723=degree>=3?m(3,2,3):0;
        const T x=u[0],y=u[1],xx=x*x,xy=x*y,yy=y*y;
        if(degree>=1){out[1]=r*x*i33;out[2]=r*y*i33;out[3]=r*z321;}
        if(degree>=2){out[4]=3*xx*i54-i32;out[5]=3*xy*i54;out[6]=3*x*z531;
            out[7]=3*yy*i54-i32;out[8]=3*y*z531;out[9]=3*z522-i32;}
        if(degree>=3){out[10]=(15*xx*i75-9*i53)*x/r;out[11]=(15*xx*i75-3*i53)*y/r;
            out[12]=(15*xx*z741-3*z521)/r;out[13]=(15*yy*i75-3*i53)*x/r;out[14]=15*xy*z741/r;
            out[15]=(15*z732-3*i53)*x/r;out[16]=(15*yy*i75-9*i53)*y/r;out[17]=(15*yy*z741-3*z521)/r;
            out[18]=(15*z732-3*i53)*y/r;out[19]=(15*z723-9*z521)/r;}
        if(power==-3)for(auto&v:out)v/=r*r*r;
        return out;
    }
    throw std::invalid_argument("unsupported radial profile");
}
}
}
#include <numeric>
#include <queue>
#include <utility>

namespace TSP_MMSI {
struct SourceValue {
    double R1=0,R2=0;
    std::array<double,3> q{};
};
struct Triangle {
    std::array<Vec,3> vertex;
    std::array<SourceValue,3> value;
};
struct Options {
    double AbsoluteTolerance=1e-16,RelativeTolerance=1e-16;
    int LeafFaces=16;
    std::uint64_t MaxNodes=100000000;
    int DerivativeOrder=3;
};
struct Model {
    std::vector<Triangle> triangles;
    std::vector<Vec> vertices;
    std::vector<std::array<int,3>> faces;
    std::vector<double> data;
    int Fields=1;
    bool FullSphere=true;
    double MinimumRadius=0,MaximumRadius=0;
    std::size_t FaceCount()const{return faces.empty()?triangles.size():faces.size();}
    SourceValue State(std::size_t vertex)const {
        const auto at=vertex*(2+Fields);SourceValue value{data[at],data[at+1],{}};
        for(int f=0;f<Fields;++f)value.q[f]=data[at+2+f];return value;
    }
    void SetState(std::size_t vertex,const SourceValue&value) {
        const auto at=vertex*(2+Fields);data[at]=value.R1;data[at+1]=value.R2;
        for(int f=0;f<Fields;++f)data[at+2+f]=value.q[f];
    }
    Triangle GetTriangle(std::size_t face)const {
        if(faces.empty())return triangles[face];Triangle t;
        for(int k=0;k<3;++k){const auto id=faces[face][k];t.vertex[k]=vertices[id];t.value[k]=State(id);}return t;
    }
};
namespace detail {
using Values=std::array<double,60>;

}

class FastShellIntegrator;
}
#include <cstring>
#include <memory>
#include <unordered_map>
#ifdef TSP_MPI
#include <mpi.h>
#endif

namespace TSP_MMSI { namespace detail {
class IndexedSource {
    friend class ::TSP_MMSI::FastShellIntegrator;
    struct Node {Vec axis{};double cosine=-1;int begin=0,end=0,left=-1,right=-1;};
    struct Info {
        Options options;Vec maximum{};
        double length=0,inner=0,outer=0;int fields=1;bool constant=false;
        std::size_t vertices=0,faces=0,nodes=0;
    } info_;
    struct Storage {
        std::vector<Vec>vertices;std::vector<double>data;
        std::vector<std::array<int,3>>faces;std::vector<int>ids;std::vector<Node>nodes;
    };
    std::shared_ptr<Storage>owned_;
    Vec*vertices_=nullptr;double*data_=nullptr;std::array<int,3>*faces_=nullptr;int*ids_=nullptr;Node*nodes_=nullptr;
#ifdef TSP_MPI
    struct Shared {
        MPI_Comm local=MPI_COMM_NULL;MPI_Win window=MPI_WIN_NULL;unsigned char*base=nullptr;int ranks=1;
        Shared(MPI_Comm comm,std::size_t bytes):local(comm) {
            int rank;MPI_Comm_rank(local,&rank);MPI_Comm_size(local,&ranks);void*pointer=nullptr;
            MPI_Win_allocate_shared(rank==0?MPI_Aint(bytes):0,1,MPI_INFO_NULL,local,&pointer,&window);
            MPI_Aint size;int unit;MPI_Win_shared_query(window,0,&size,&unit,&pointer);
            base=static_cast<unsigned char*>(pointer);MPI_Win_lock_all(0,window);
        }
        ~Shared(){MPI_Win_unlock_all(window);MPI_Win_free(&window);MPI_Comm_free(&local);}
    };
    std::shared_ptr<Shared>shared_;
#endif
    struct Hash {
        std::size_t operator()(const std::array<double,8>&v)const {
            std::size_t h=0;for(double x:v)h^=std::hash<double>{}(x)+0x9e3779b9+(h<<6)+(h>>2);return h;
        }
    };
    static void Pack(Model&model) {
        if(model.triangles.empty())return;
        if(!model.faces.empty()||!model.vertices.empty()||!model.data.empty())throw std::invalid_argument("mixed source layouts");
        std::unordered_map<std::array<double,8>,int,Hash>vertices;
        model.faces.reserve(model.triangles.size());
        for(const auto&t:model.triangles){std::array<int,3>face;
            for(int k=0;k<3;++k){const auto&v=t.vertex[k];const auto&s=t.value[k];
                const std::array<double,8>key{v[0],v[1],v[2],s.R1,s.R2,s.q[0],model.Fields==3?s.q[1]:0,model.Fields==3?s.q[2]:0};
                auto found=vertices.find(key);
                if(found==vertices.end()){const int id=int(model.vertices.size());vertices.emplace(key,id);face[k]=id;model.vertices.push_back(v);
                    model.data.push_back(s.R1);model.data.push_back(s.R2);for(int f=0;f<model.Fields;++f)model.data.push_back(s.q[f]);}
                else face[k]=found->second;}
            model.faces.push_back(face);}
        std::vector<Triangle>().swap(model.triangles);
    }
    void BindOwned() {
        vertices_=owned_->vertices.data();data_=owned_->data.data();faces_=owned_->faces.data();ids_=owned_->ids.data();nodes_=owned_->nodes.data();
    }
    Vec Center(int face)const {
        const auto&v=faces_[face];return Scale(Add(Add(vertices_[v[0]],vertices_[v[1]]),vertices_[v[2]]),1.0/3);
    }
    int Build(int begin,int end,const std::vector<Vec>&centers) {
        Node node;node.begin=begin;node.end=end;Vec average{},lo{2,2,2},hi{-2,-2,-2};
        for(int j=begin;j<end;++j){const auto&c=centers[ids_[j]];average=Add(average,c);
            for(int k=0;k<3;++k){lo[k]=std::min(lo[k],c[k]);hi[k]=std::max(hi[k],c[k]);}}
        if(Norm(average)>1e-12){node.axis=Unit(average);node.cosine=1;
            for(int j=begin;j<end;++j)for(const auto&v:Vertices(ids_[j]))node.cosine=std::min(node.cosine,Dot(v,node.axis));
            if(node.cosine<0)node.cosine=-1;}
        const int id=int(owned_->nodes.size());owned_->nodes.push_back(node);
        if(end-begin>info_.options.LeafFaces){int axis=0;for(int k=1;k<3;++k)if(hi[k]-lo[k]>hi[axis]-lo[axis])axis=k;
            const int mid=(begin+end)/2;std::nth_element(ids_+begin,ids_+mid,ids_+end,[&](int a,int b){return centers[a][axis]<centers[b][axis];});
            const int left=Build(begin,mid,centers),right=Build(mid,end,centers);owned_->nodes[id].left=left;owned_->nodes[id].right=right;}
        return id;
    }
    template<class T>static T*Slice(unsigned char*base,std::size_t&offset,std::size_t count) {
        offset=(offset+alignof(T)-1)&~(alignof(T)-1);const auto at=offset;offset+=count*sizeof(T);return base?reinterpret_cast<T*>(base+at):nullptr;
    }
    std::size_t BindShared(unsigned char*base) {
        std::size_t at=0;vertices_=Slice<Vec>(base,at,info_.vertices);data_=Slice<double>(base,at,info_.vertices*(2+info_.fields));
        faces_=Slice<std::array<int,3>>(base,at,info_.faces);ids_=Slice<int>(base,at,info_.faces);nodes_=Slice<Node>(base,at,info_.nodes);return at;
    }
    IndexedSource()=default;
  public:
    explicit IndexedSource(Model model,Options options):owned_(std::make_shared<Storage>()) {
        if(!model.FullSphere)throw std::invalid_argument("fast MMSI requires a full-sphere model");
        if(model.Fields!=1&&model.Fields!=3)throw std::invalid_argument("source fields must be one or three");
        if(!(options.AbsoluteTolerance>0&&options.RelativeTolerance>0)||!std::isfinite(options.AbsoluteTolerance)||!std::isfinite(options.RelativeTolerance)||
           options.MaxNodes<80||options.LeafFaces<1||options.DerivativeOrder<0||options.DerivativeOrder>3)
            throw std::invalid_argument("invalid fast integration options");
        if(!(model.MinimumRadius>=0&&model.MaximumRadius>=0)||!std::isfinite(model.MinimumRadius)||!std::isfinite(model.MaximumRadius))
            throw std::invalid_argument("invalid radial bounds");
        Pack(model);if(model.faces.empty()||model.vertices.empty()||model.data.size()!=model.vertices.size()*(2+model.Fields))throw std::invalid_argument("invalid indexed source");
        info_.options=options;info_.fields=model.Fields;info_.constant=model.FullSphere;info_.inner=std::numeric_limits<double>::infinity();
        info_.vertices=model.vertices.size();info_.faces=model.faces.size();
        owned_->vertices=std::move(model.vertices);owned_->data=std::move(model.data);owned_->faces=std::move(model.faces);
        owned_->ids.resize(info_.faces);std::iota(owned_->ids.begin(),owned_->ids.end(),0);BindOwned();
        for(auto&v:owned_->vertices){const double length=Norm(v);if(!(length>0)||!std::isfinite(length))throw std::invalid_argument("invalid direction");
            if(std::abs(length-1)>8*std::numeric_limits<double>::epsilon())v=Scale(v,1/length);}
        for(const auto&face:owned_->faces)for(int id:face)if(id<0||std::size_t(id)>=info_.vertices)throw std::invalid_argument("triangle index out of range");
        const auto first=State(faces_[0][0]);double area=0;
        for(std::size_t i=0;i<info_.faces;++i){auto v=Vertices(int(i));double det=Dot(v[0],Cross(v[1],v[2]));
            if(det<0){std::swap(faces_[i][1],faces_[i][2]);v=Vertices(int(i));det=-det;}
            if(det<1e-15)throw std::invalid_argument("degenerate source triangle");area+=Area(v);
            for(int id:faces_[i]){const auto s=State(id);if(!(s.R1>=0&&s.R2>s.R1)||!std::isfinite(s.R2))throw std::invalid_argument("invalid source radii");
                info_.inner=std::min(info_.inner,s.R1);info_.length=std::max(info_.length,s.R2);info_.constant=info_.constant&&s.R1==first.R1&&s.R2==first.R2;
                for(int f=0;f<info_.fields;++f){if(!std::isfinite(s.q[f]))throw std::invalid_argument("nonfinite source");
                    info_.maximum[f]=std::max(info_.maximum[f],std::abs(s.q[f]));info_.constant=info_.constant&&s.q[f]==first.q[f];}}}
        if(model.FullSphere&&std::abs(area-4*PI)>1e-8)throw std::invalid_argument("source triangles do not cover the sphere");
        if(model.MinimumRadius>0)info_.inner=std::min(info_.inner,model.MinimumRadius);
        info_.length=std::max(info_.length,model.MaximumRadius);info_.outer=info_.length;
        int exponent;std::frexp(info_.length,&exponent);info_.length=std::ldexp(1.,exponent-1);
        info_.inner/=info_.length;info_.outer/=info_.length;
        if(!(info_.inner>0))throw std::invalid_argument("analytical near TSP requires positive inner radius");
        for(std::size_t i=0;i<info_.vertices;++i){data_[i*(2+info_.fields)]/=info_.length;data_[i*(2+info_.fields)+1]/=info_.length;}
        std::size_t leaves=1;while(info_.faces>std::size_t(options.LeafFaces)*leaves)leaves*=2;
        std::vector<Vec>centers(info_.faces);for(std::size_t i=0;i<info_.faces;++i)centers[i]=Center(int(i));
        owned_->nodes.reserve(2*leaves-1);Build(0,int(info_.faces),centers);info_.nodes=owned_->nodes.size();BindOwned();
    }
    std::array<Vec,3>Vertices(int face)const {const auto&v=faces_[face];return {vertices_[v[0]],vertices_[v[1]],vertices_[v[2]]};}
    SourceValue State(std::size_t id)const {
        const auto at=id*(2+info_.fields);SourceValue s{data_[at],data_[at+1],{}};
        for(int f=0;f<info_.fields;++f)s.q[f]=data_[at+2+f];return s;
    }
    Triangle GetTriangle(int face)const {
        Triangle t;t.vertex=Vertices(face);for(int k=0;k<3;++k)t.value[k]=State(faces_[face][k]);return t;
    }
    struct FaceGeometry {std::array<Vec,3>edge;Vec center,edgeNorm;};
    static FaceGeometry Geometry(const std::array<Vec,3>&v) {
        FaceGeometry g;g.center=Scale(Add(Add(v[0],v[1]),v[2]),1.0/3);
        for(int k=0;k<3;++k){g.edge[k]=Cross(v[(k+1)%3],v[(k+2)%3]);g.edgeNorm[k]=Norm(g.edge[k]);}return g;
    }
    bool Contains(int face,const Vec&u)const {
        const auto v=Vertices(face);for(int k=0;k<3;++k){const auto edge=Cross(v[(k+1)%3],v[(k+2)%3]);
            if(Dot(edge,u)<-32*std::numeric_limits<double>::epsilon()*Norm(edge))return false;}return true;
    }
    int Find(const Vec&u)const {
        int stack[64],size=1;stack[0]=0;
        while(size){const auto&node=nodes_[stack[--size]];
            if(node.cosine>-1&&Dot(node.axis,u)<node.cosine-4e-14)continue;
            if(node.left>=0){stack[size++]=node.left;stack[size++]=node.right;}
            else for(int j=node.begin;j<node.end;++j)if(Contains(ids_[j],u))return ids_[j];}
        // A direction coinciding with a mesh vertex falls through the gap between the incident
        // triangles; take the face that leaves it least outside. Keep failing if it is really outside.
        int best=-1;double bestMargin=-std::numeric_limits<double>::infinity();
        for(int id:Candidates(u,1e-6)) {const auto v=Vertices(id);double margin=std::numeric_limits<double>::infinity();
            for(int k=0;k<3;++k) {const auto edge=Cross(v[(k+1)%3],v[(k+2)%3]);margin=std::min(margin,Dot(edge,u)/Norm(edge));}
            if(margin>bestMargin){bestMargin=margin;best=id;}}
        return bestMargin<-1e-9?-1:best;
    }
    std::vector<int>Candidates(const Vec&u,double angle)const {
        std::vector<int>faces;int stack[64],size=1;stack[0]=0;
        while(size){const auto&node=nodes_[stack[--size]];
            if(node.cosine>-1&&Dot(node.axis,u)<std::cos(std::min(PI,angle+std::acos(std::min(std::max(node.cosine,-1.),1.))))-4e-14)continue;
            if(node.left>=0){stack[size++]=node.left;stack[size++]=node.right;}
            else for(int j=node.begin;j<node.end;++j)faces.push_back(ids_[j]);}
        std::sort(faces.begin(),faces.end());return faces;
    }
    Values Target(const Vec&point)const {
        const auto&o=info_.options;Values target{};const double r=Norm(point),mass=4*PI*((info_.outer-info_.inner)*
            (info_.outer*info_.outer+info_.outer*info_.inner+info_.inner*info_.inner)/3);
        for(int f=0;f<info_.fields;++f)for(int k=0;k<20;++k)if(Degree(k)<=o.DerivativeOrder)
            target[f*20+k]=o.AbsoluteTolerance*std::pow(info_.length,Degree(k)-2)/info_.fields+
                o.RelativeTolerance*mass*info_.maximum[f]/std::pow(std::max(r,info_.outer/2),Degree(k)+1);return target;
    }
    std::size_t Bytes()const {
        std::size_t at=0;Slice<Vec>(nullptr,at,info_.vertices);Slice<double>(nullptr,at,info_.vertices*(2+info_.fields));
        Slice<std::array<int,3>>(nullptr,at,info_.faces);Slice<int>(nullptr,at,info_.faces);Slice<Node>(nullptr,at,info_.nodes);return at;
    }
    int SharedRanks()const {
#ifdef TSP_MPI
        if(shared_)return shared_->ranks;
#endif
        return 1;
    }
#ifdef TSP_MPI
    static std::shared_ptr<const IndexedSource>PrepareMPI(Model model,Options options,MPI_Comm world) {
        int rank;MPI_Comm_rank(world,&rank);std::shared_ptr<IndexedSource>source;std::string error;
        if(rank==0)try{source=std::make_shared<IndexedSource>(std::move(model),options);}catch(const std::exception&e){error=e.what();}
        int length=int(error.size());MPI_Bcast(&length,1,MPI_INT,0,world);
        if(length){error.resize(length);MPI_Bcast(&error[0],length,MPI_CHAR,0,world);throw std::runtime_error(error);}
        Info info{};if(rank==0)info=source->info_;MPI_Bcast(&info,sizeof(info),MPI_BYTE,0,world);
        MPI_Comm local;MPI_Comm_split_type(world,MPI_COMM_TYPE_SHARED,rank,MPI_INFO_NULL,&local);int localRank;MPI_Comm_rank(local,&localRank);
        auto result=std::shared_ptr<IndexedSource>(new IndexedSource);result->info_=info;
        result->shared_=std::make_shared<Shared>(local,result->Bytes());result->BindShared(result->shared_->base);
        if(rank==0){const auto&s=*source->owned_;
            std::memcpy(result->vertices_,s.vertices.data(),s.vertices.size()*sizeof(Vec));
            std::memcpy(result->data_,s.data.data(),s.data.size()*sizeof(double));
            std::memcpy(result->faces_,s.faces.data(),s.faces.size()*sizeof(std::array<int,3>));
            std::memcpy(result->ids_,s.ids.data(),s.ids.size()*sizeof(int));std::memcpy(result->nodes_,s.nodes.data(),s.nodes.size()*sizeof(Node));source.reset();}
        MPI_Comm leaders;MPI_Comm_split(world,localRank==0?0:MPI_UNDEFINED,rank,&leaders);
        if(localRank==0){std::size_t remaining=result->Bytes();auto*data=result->shared_->base;
            while(remaining){const int count=int(std::min<std::size_t>(remaining,1u<<28));MPI_Bcast(data,count,MPI_BYTE,0,leaders);data+=count;remaining-=count;}}
        MPI_Win_sync(result->shared_->window);MPI_Barrier(local);MPI_Win_sync(result->shared_->window);
        if(leaders!=MPI_COMM_NULL)MPI_Comm_free(&leaders);return result;
    }
#endif
};
}}
#include <map>

namespace TSP_MMSI { namespace detail {
class NearBoundary {
    using WideVec=std::array<dd_real,3>;
    template<class R>using Frame=std::array<std::array<R,3>,3>;
    struct Arc {std::array<Vec,2>vertex;Vec normal;WideVec wideNormal;double c;dd_real wideC;};
    struct Sector {std::vector<int>ordinary,wide;};
    std::vector<Arc>arcs_;
    std::vector<Sector>sectors_;
    Frame<double>axes_;Frame<dd_real>wideAxes_;
    const dd_real pi_{"3.1415926535897932384626433832795028841971693993751"};
    template<class R>static std::array<R,3>Column(const Frame<R>&a,int k){return {a[0][k],a[1][k],a[2][k]};}
    template<class R>static std::array<R,3>Convert(const Vec&v){return {R(v[0]),R(v[1]),R(v[2])};}
    template<class R>const Frame<R>&Axes()const {if constexpr(std::is_same<R,double>::value)return axes_;else return wideAxes_;}
    template<class R>const std::array<R,3>&Normal(const Arc&a)const {if constexpr(std::is_same<R,double>::value)return a.normal;else return a.wideNormal;}
    template<class R>R C(const Arc&a)const {if constexpr(std::is_same<R,double>::value)return a.c;else return a.wideC;}
    template<class R>std::array<R,3>Direction(R cp,R sp)const {
        const auto&a=Axes<R>();std::array<R,3>d;for(int k=0;k<3;++k)d[k]=a[k][0]*cp+a[k][1]*sp;return d;
    }
    template<class R>std::vector<int>Active(R phi)const {
        using std::cos;using std::sin;const auto d=Direction<R>(cos(phi),sin(phi)),n=Column(Axes<R>(),2);std::vector<int>out;
        for(std::size_t i=0;i<arcs_.size();++i){const auto&a=arcs_[i];const R c=C<R>(a);if(c==0)continue;
            const R s=Dot(Normal<R>(a),d);std::array<R,3>u;
            for(int k=0;k<3;++k)u[k]=(c>0?-s:s)*n[k]+(c>0?c:-c)*d[k];
            if(Dot(Cross(Convert<R>(a.vertex[0]),u),Normal<R>(a))>=0&&
               Dot(Cross(u,Convert<R>(a.vertex[1])),Normal<R>(a))>=0)out.push_back(int(i));}
        return out;
    }
  public:
    std::vector<dd_real>breaks;
    NearBoundary(const std::map<std::array<Vec,2>,int>&edges,const Frame<double>&axes,const Frame<dd_real>&wideAxes):axes_(axes),wideAxes_(wideAxes) {
        breaks={dd_real(0),2*pi_};
        for(const auto&edge:edges)if(edge.second){
            if(std::abs(edge.second)!=1)throw std::runtime_error("nonmanifold near boundary");
            Arc a;a.vertex=edge.first;if(edge.second<0)std::swap_ranges(a.vertex[0].begin(),a.vertex[0].end(),a.vertex[1].begin());
            a.normal=Cross(a.vertex[0],a.vertex[1]);a.wideNormal=Cross(Convert<dd_real>(a.vertex[0]),Convert<dd_real>(a.vertex[1]));
            a.c=Dot(a.normal,Column(axes_,2));a.wideC=Dot(a.wideNormal,Column(wideAxes_,2));arcs_.push_back(a);
            for(const auto&v:edge.first){const auto w=Convert<dd_real>(v);dd_real p=atan2(Dot(w,Column(wideAxes_,1)),Dot(w,Column(wideAxes_,0)));
                if(p<0)p+=2*pi_;breaks.push_back(p);}}
        std::sort(breaks.begin(),breaks.end());breaks.erase(std::unique(breaks.begin(),breaks.end()),breaks.end());
        for(std::size_t i=1;i<breaks.size();++i){const dd_real mid=(breaks[i-1]+breaks[i])/2;
            sectors_.push_back({Active(::to_double(mid)),Active(mid)});}
    }
    template<class R>std::vector<std::pair<R,R>>Intervals(int sector,R cp,R sp)const {
        using std::atan2;const auto direction=Direction<R>(cp,sp);const auto&list=[&]()->const std::vector<int>&{
            if constexpr(std::is_same<R,double>::value)return sectors_[sector].ordinary;else return sectors_[sector].wide;}();
        struct Crossing{R first;int second;};std::vector<Crossing>crossings;crossings.reserve(list.size());
        for(int id:list){const auto&a=arcs_[id];const R c=C<R>(a),s=Dot(Normal<R>(a),direction);
            crossings.push_back({c>0?atan2(c,-s):atan2(-c,s),c>0?-1:1});}
        std::sort(crossings.begin(),crossings.end(),[](const Crossing&a,const Crossing&b){
            return a.first<b.first||(!(b.first<a.first)&&a.second<b.second);});
        std::vector<std::pair<R,R>>out;int inside=1;R begin=0;
        for(std::size_t i=0;i<crossings.size();){const R theta=crossings[i].first;int change=0;
            do{change+=crossings[i++].second;}while(i<crossings.size()&&crossings[i].first==theta);
            const int next=inside+change;if(inside>0&&next<=0)begin=theta;
            if(inside<=0&&next>0&&theta>begin)out.push_back({begin,theta});inside=next;}
        if(inside<=0){R pi;if constexpr(std::is_same<R,double>::value)pi=PI;else pi=pi_;out.push_back({begin,pi});}return out;
    }
};
}}
#include <map>

namespace TSP_MMSI { namespace detail {
class AnalyticalNear {
    using WideValues=std::array<dd_real,60>;
    using Weight=std::array<dd_real,3>;
    struct Group {Weight area{},density{};std::map<std::array<Vec,2>,Weight>edges;};
    std::map<std::array<double,2>,Group>groups_;
    Vec point_;TSP::Options options_;int fields_,components_;WideValues value_{};
  public:
    AnalyticalNear(const Vec&point,const Options&options,int fields,int components):point_(point),fields_(fields),components_(components) {
        options_.absolute_tolerance=options.AbsoluteTolerance;options_.relative_tolerance=options.RelativeTolerance;
    }
    bool Add(const std::array<Vec,3>&vertices,const SourceValue&state) {
        if(state.q==Vec{})return false;
        using namespace TSP::detail;
        const auto radius=Real(state.R2),inner=Real(state.R1)/radius;
        const auto p=rscale(convert(point_),1/radius);const Real R=rnorm(p);const double r=static_cast<double>(R);
        if(r<.03*static_cast<double>(inner)||r>20){
            const auto field=gravity_tensor(TSP::Prism{vertices,state.R1,state.R2,1},point_,options_,nullptr);
            for(int f=0;f<fields_;++f)for(int k=0;k<components_;++k)value_[f*20+k]+=state.q[f]*field[k].value;return true;}
        const Geometry geometry(TSP::Prism{vertices,state.R1,state.R2,1});
        Real fraction=0;const Real eps("2e-14");
        if(R>=inner-eps&&R<=1+eps){
            const RV n=rscale(p,1/R);int sides[3],count=0;fraction=1;
            for(int k=0;k<3;++k){const Real d=rdot(geometry.edges[k].N,n);
                if(d<-eps){fraction=0;break;}if(tsp_native_math::abs(d)<=eps)sides[count++]=k;}
            if(fraction!=0){
                if(tsp_native_math::abs(R-inner)<=eps||tsp_native_math::abs(R-1)<=eps)
                    throw std::domain_error("full third-order tensors are undefined on a source interface (2e-14 relative guard)");
                if(count==1)fraction=Real(.5);
                else if(count==2){const auto&a=geometry.edges[sides[0]].N;const auto&b=geometry.edges[sides[1]].N;
                    fraction=tsp_native_math::atan2(rnorm(rcross(a,b)),-rdot(a,b))/(2*tsp_math::PI);}
                else if(count>2)throw std::domain_error("unresolved source vertex");
            }
        }
        auto&group=groups_[{state.R1,state.R2}];const auto area=solid_angle(geometry);
        for(int f=0;f<fields_;++f){group.area[f]+=state.q[f]*area.value;group.density[f]+=state.q[f]*fraction.value;}
        for(int k=0;k<3;++k){auto a=vertices[k],b=vertices[(k+1)%3];const bool forward=a<b;if(!forward)std::swap_ranges(a.begin(),a.end(),b.begin());
            auto&weight=group.edges[{a,b}];for(int f=0;f<fields_;++f)weight[f]+=(forward?1:-1)*state.q[f];}
        return true;
    }
    WideValues Evaluate(std::size_t&edgeEvaluations) {
        using namespace TSP::detail;
        for(const auto&entry:groups_){const Real radius(entry.first[1]);
            const EllipticKernel kernel(Real(entry.first[0])/radius,rscale(convert(point_),1/radius));std::array<tsp_math::J,3>total;
            const auto poles=kernel.jm-kernel.jp;std::array<Real,4>scale;
            const int degree=components_==20?3:components_==10?2:components_==4?1:0;
            for(int d=0;d<=degree;++d)scale[d]=tsp_native_math::pow(radius,2-d);
            for(int f=0;f<fields_;++f)total[f]=Real(entry.second.area[f])*poles/2;
            const RV n{kernel.n[0].a[0],kernel.n[1].a[0],kernel.n[2].a[0]};
            for(const auto&arc:entry.second.edges){if(arc.second==Weight{})continue;
                const Edge edge(runit(convert(arc.first[0])),runit(convert(arc.first[1])));
                if(kernel.sign_change&&tsp_native_math::abs(rdot(edge.N,n))<Real("2e-14")){
                    const Real t=tsp_native_math::atan2(rdot(n,edge.b),rdot(n,edge.a));
                    if(t>=Real("-2e-14")&&t<=edge.gamma+Real("2e-14"))
                        throw std::domain_error("full third-order tensors are undefined on a noncancelled near-source boundary");}
                tsp_math::J field;kernel.AddEdge(field,edge);++edgeEvaluations;
                for(int f=0;f<fields_;++f)if(arc.second[f]!=0)total[f]+=Real(arc.second[f])*field;}
            for(int f=0;f<fields_;++f){const auto tensor=tensor_derivatives(total[f]*kernel.inverseR);
                auto check=[&](int i,int j,int k,const Real&target){
                    const Real budget=3*Real(options_.absolute_tolerance)+Real(options_.relative_tolerance)*
                        (tsp_native_math::abs(tensor[i])+tsp_native_math::abs(tensor[j])+tsp_native_math::abs(tensor[k]));
                    return tsp_native_math::abs(tensor[i]+tensor[j]+tensor[k]-target)<=budget;};
                if(((tsp_math::active_jet_mask&1023)==1023&&!check(4,7,9,-4*tsp_math::PI*Real(entry.second.density[f])))||
                   (tsp_math::active_jet_mask==tsp_math::JET_FULL_MASK&&(!check(10,13,15,Real(0))||!check(11,16,18,Real(0))||!check(12,17,19,Real(0)))))
                    throw std::runtime_error("near-source Poisson residual failed");
                for(int k=0;k<components_;++k)value_[f*20+k]+=(tensor[k]*scale[Degree(k)]).value;}}
        return value_;
    }
};
}}
#include <cstdlib>
#include <map>
#include <memory>
#include <string>

namespace TSP_MMSI {
struct FastDiagnostics {
    std::size_t AnalyticalTriangles=0,AnalyticalEvaluations=0,AnalyticalEdges=0,InitialPanels=0,Refinements=0;
    std::uint64_t KernelNodes=0;
    double ErrorRatio=0;
    bool QuadratureConverged=false,LimitReached=false,NearApproximation=false;
};
class FastShellIntegrator {
    std::shared_ptr<const detail::IndexedSource>source_;
    double angle_;
    using Values=detail::Values;
    using WideValues=std::array<dd_real,60>;
    struct Panel {
        double x0=0,x1=1,priority=0;int sector=0;
        dd_real p0=0,p1=0;
        WideValues value{};Values error{};
        bool operator<(const Panel&b)const{return priority<b.priority;}
    };
  public:
    using Settings=Options;
    using Diagnostic=FastDiagnostics;
    static constexpr std::size_t MaxRefinements=1000;
    explicit FastShellIntegrator(Model model,Options options={}):FastShellIntegrator(std::move(model),options,
        std::getenv("TSP_MMSI_NEAR_DEGREES")?std::stod(std::getenv("TSP_MMSI_NEAR_DEGREES")):1){}
    FastShellIntegrator(Model model,Options options,double nearDegrees):angle_(nearDegrees*PI/180) {
        if(!(nearDegrees>0&&nearDegrees<90))throw std::invalid_argument("fast MMSI requires 0 < near angle < 90");
        source_=std::make_shared<detail::IndexedSource>(std::move(model),options);
    }
    std::size_t SourceBytes()const{return source_->Bytes();}
    int SharedRanks()const{return source_->SharedRanks();}
#ifdef TSP_MPI
    static std::unique_ptr<FastShellIntegrator>PrepareMPI(Model model,Options options,MPI_Comm world) {
        int rank;MPI_Comm_rank(world,&rank);double angle=1;std::string error;
        if(rank==0)try {if(const char*p=std::getenv("TSP_MMSI_NEAR_DEGREES"))angle=std::stod(p);
            if(!(angle>0&&angle<90))throw std::invalid_argument("invalid near angle");
        }catch(const std::exception&e){error=e.what();}
        int length=int(error.size());MPI_Bcast(&length,1,MPI_INT,0,world);
        if(length){error.resize(length);MPI_Bcast(&error[0],length,MPI_CHAR,0,world);throw std::runtime_error(error);}
        MPI_Bcast(&angle,1,MPI_DOUBLE,0,world);
        auto source=detail::IndexedSource::PrepareMPI(std::move(model),options,world);
        return std::unique_ptr<FastShellIntegrator>(new FastShellIntegrator(std::move(source),angle));
    }
#endif
    Values Fields(const Vec&physicalPoint,FastDiagnostics*diagnostics=nullptr)const {
        using namespace detail;
        FastDiagnostics diag;const auto&source=*source_;const auto&options=source.info_.options;
        for(double v:physicalPoint)if(!std::isfinite(v))throw std::invalid_argument("nonfinite observation");
        const Vec point=Scale(physicalPoint,1/source.info_.length);const auto wp=TSP::detail::convert(point);
        const auto wr=TSP::detail::rnorm(wp);if(!(wr.value>0))throw std::invalid_argument("fast MMSI requires a noncentral observation");
        const auto wn=TSP::detail::runit(wp);
        const auto we1=TSP::detail::runit(TSP::detail::rcross(wn,TSP::detail::convert(std::abs(static_cast<double>(wn[2]))<.9?Vec{0,0,1}:Vec{0,1,0})));
        const auto we2=TSP::detail::rcross(wn,we1);Vec n;
        for(int k=0;k<3;++k)n[k]=static_cast<double>(wn[k]);
        auto axes=[&](auto tag){using R=decltype(tag);std::array<std::array<R,3>,3>a;
            for(int k=0;k<3;++k){if constexpr(std::is_same<R,double>::value)a[k]={static_cast<double>(we1[k]),static_cast<double>(we2[k]),n[k]};
                else a[k]={we1[k].value,we2[k].value,wn[k].value};}return a;};
        auto convert=[](const Vec&v,auto tag){using R=decltype(tag);return std::array<R,3>{R(v[0]),R(v[1]),R(v[2])};};
        auto scalar=[](const auto&v){if constexpr(std::is_same<typename std::decay<decltype(v)>::type,double>::value)return v;else return ::to_double(v);};
        const dd_real pi("3.1415926535897932384626433832795028841971693993751");
        const int components=options.DerivativeOrder==3?20:options.DerivativeOrder==2?10:options.DerivativeOrder==1?4:1;
        const tsp_math::JetScope scope((std::uint32_t(1)<<components)-1);
        const int containing=source.Find(n);if(containing<0)throw std::runtime_error("source mesh does not cover projection");
        auto stateAt=[&](const Triangle&t,const std::array<Vec,3>&edges,const Vec&u){Vec w;double sum=0;
            for(int k=0;k<3;++k){w[k]=Dot(edges[k],u);sum+=w[k];}for(auto&x:w)x/=sum;
            auto v=t.value[0];for(int k=1;k<3;++k){v.R1+=w[k]*(t.value[k].R1-t.value[0].R1);v.R2+=w[k]*(t.value[k].R2-t.value[0].R2);
                for(int f=0;f<source.info_.fields;++f)v.q[f]+=w[k]*(t.value[k].q[f]-t.value[0].q[f]);}return v;};
        const auto containingTriangle=source.GetTriangle(containing);const auto projection=stateAt(containingTriangle,source.Geometry(containingTriangle.vertex).edge,n);
        const dd_real guard=dd_real("2e-14")*projection.R2;
        if(abs(wr.value-projection.R1)<=guard||abs(wr.value-projection.R2)<=guard)
            throw std::domain_error("full third-order tensors are undefined on a source interface (2e-14 relative guard)");
        WideValues nearValue{};
        const auto target=source.Target(point);
        std::map<std::array<Vec,2>,int>edges;
        AnalyticalNear analyticalNear(point,options,source.info_.fields,components);
        for(int face:source.Candidates(n,angle_)){const auto t=source.GetTriangle(face);const auto geometry=source.Geometry(t.vertex);
            const Vec center=Unit(geometry.center);double cosine=1;for(const auto&v:t.vertex)cosine=std::min(cosine,Dot(center,v));
            if(Dot(center,n)<std::cos(std::min(PI,angle_+std::acos(std::min(std::max(cosine,-1.),1.)))))continue;
            Vec closest=t.vertex[0];double best=Dot(closest,n);const bool projected=source.Contains(face,n);
            if(projected){closest=n;best=1;}
            else for(int k=0;k<3;++k){if(Dot(t.vertex[k],n)>best){closest=t.vertex[k];best=Dot(closest,n);}
                const Vec normal=Scale(geometry.edge[k],1/geometry.edgeNorm[k]);const Vec projected=Subtract(n,Scale(normal,Dot(n,normal)));
                if(Norm(projected)==0)continue;const Vec foot=Unit(projected);
                if(Dot(Cross(t.vertex[(k+1)%3],foot),normal)>=0&&Dot(Cross(foot,t.vertex[(k+2)%3]),normal)>=0&&Dot(n,foot)>best){closest=foot;best=Dot(n,foot);}}
            if(best<std::cos(angle_))continue;++diag.AnalyticalTriangles;
            bool varying=false;for(int k=1;k<3;++k)varying=varying||t.value[k].R1!=t.value[0].R1||t.value[k].R2!=t.value[0].R2||t.value[k].q!=t.value[0].q;
            diag.NearApproximation=diag.NearApproximation||varying;
            auto state=varying?(projected?projection:stateAt(t,geometry.edge,closest)):t.value[0];
            diag.AnalyticalEvaluations+=analyticalNear.Add(t.vertex,state);
            for(int k=0;k<3;++k){auto a=t.vertex[k],b=t.vertex[(k+1)%3];const bool forward=a<b;if(!forward)std::swap_ranges(a.begin(),a.end(),b.begin());edges[{a,b}]+=forward?1:-1;}}
        const auto analyticalValue=analyticalNear.Evaluate(diag.AnalyticalEdges);
        for(int k=0;k<60;++k)nearValue[k]+=analyticalValue[k];
        auto farTarget=target;for(auto&x:farTarget)x*=.25;
        const auto wideAxes=axes(dd_real(0));
        const NearBoundary boundary(edges,axes(0.),wideAxes);const auto&breaks=boundary.breaks;
        auto quadTyped=[&](const Panel&panel,int nx,int ny,auto tag){using R=decltype(tag);using std::sin;using std::cos;using std::sqrt;using std::log;using std::exp;
            static const Rule<R>r8(8),r12(12);const auto&rx=nx==8?r8:r12;const auto&ry=ny==8?r8:r12;
            std::array<R,60>value{},correction{};const auto a=axes(R(0));R p0,dp,radius,p;
            if constexpr(std::is_same<R,double>::value){p0=scalar(panel.p0);dp=scalar(panel.p1-panel.p0);radius=scalar(wr.value);p=PI;}
            else {p0=panel.p0;dp=panel.p1-panel.p0;radius=wr.value;p=pi;}
            for(int j=0;j<ny;++j){const R phi=p0+dp*ry.x[j],cp=cos(phi),sp=sin(phi);
                for(auto range:boundary.Intervals(panel.sector,cp,sp)){const R sl=sin(range.first/2),su=sin(range.second/2),lower=2*sl*sl,upper=2*su*su;
                    if(!(lower>0))throw std::runtime_error("near union does not cover projection");const R span=log(upper/lower);const bool antipodal=range.second==p;int previous=-1;
                    for(int i=0;i<nx;++i){if(diag.KernelNodes>=options.MaxNodes)throw std::runtime_error("fast MMSI node budget exhausted");
                        const R x=panel.x0+(panel.x1-panel.x0)*rx.x[i];R d,jacobian=span,complement;
                        if(antipodal){const R tail=1-x,z=span*tail*tail;d=upper*exp(-z);complement=-2*radial::Expm1(-z);jacobian*=2*tail;}
                        else {d=lower*exp(span*x);complement=2-d;}
                        const R mu=std::min(std::max(R(1)-d,R(-1)),R(1)),sine=sqrt(std::max(R(0),d*complement));const std::array<R,3>local{sine*cp,sine*sp,mu};
                        auto t=containingTriangle;std::array<R,3>w{};
                        if(!source.info_.constant){std::array<R,3>u;Vec lookup;for(int k=0;k<3;++k){u[k]=a[k][0]*local[0]+a[k][1]*local[1]+a[k][2]*local[2];lookup[k]=scalar(u[k]);}
                            const int face=previous>=0&&source.Contains(previous,lookup)?previous:source.Find(lookup);previous=face;
                            if(face<0)throw std::runtime_error("source lookup failed");t=source.GetTriangle(face);R sum=0;
                            for(int k=0;k<3;++k){w[k]=Dot(Cross(convert(t.vertex[(k+1)%3],R(0)),convert(t.vertex[(k+2)%3],R(0))),u);sum+=w[k];}for(auto&v:w)v/=sum;}
                        R bottom=t.value[0].R1,top=t.value[0].R2;for(int k=1;k<3;++k){bottom+=w[k]*(t.value[k].R1-t.value[0].R1);top+=w[k]*(t.value[k].R2-t.value[0].R2);}
                        const auto kernel=RadialValues<R>(bottom,top,local,radius,0,options.DerivativeOrder);
                        const R weight=d*jacobian*(panel.x1-panel.x0)*dp*rx.w[i]*ry.w[j];
                        for(int f=0;f<source.info_.fields;++f){R density=t.value[0].q[f];for(int k=1;k<3;++k)density+=w[k]*(t.value[k].q[f]-t.value[0].q[f]);
                            for(int k=0;k<components;++k){const R term=weight*(density*kernel[k]);
                                if constexpr(std::is_same<R,double>::value)Accumulate(term,value[f*20+k],correction[f*20+k]);else value[f*20+k]+=term;}}
                        ++diag.KernelNodes;}}}
            WideValues out;for(int k=0;k<60;++k)out[k]=value[k]+correction[k];return out;};
        auto quad=[&](const Panel&p,int nx,int ny){return source.info_.constant?quadTyped(p,nx,ny,dd_real(0)):quadTyped(p,nx,ny,0.);};
        auto assess=[&](Panel&p){const auto low=quad(p,8,8);p.value=quad(p,12,12);p.priority=0;
            for(int f=0;f<source.info_.fields;++f)for(int k=0;k<components;++k){const int at=f*20+k;p.error[at]=2*::to_double(abs(p.value[at]-low[at]));p.priority=std::max(p.priority,p.error[at]/farTarget[at]);}};
        std::priority_queue<Panel>queue;WideValues total{},error{};
        auto update=[&](const Panel&p,double sign){for(int k=0;k<60;++k){total[k]+=sign*p.value[k];error[k]+=sign*p.error[k];}};
        for(std::size_t k=1;k<breaks.size();++k)if(breaks[k]>breaks[k-1]){Panel p;p.sector=int(k-1);p.p0=breaks[k-1];p.p1=breaks[k];assess(p);update(p,1);queue.push(std::move(p));}
        diag.InitialPanels=queue.size();for(;;){diag.ErrorRatio=0;for(int f=0;f<source.info_.fields;++f)for(int k=0;k<components;++k)
                diag.ErrorRatio=std::max(diag.ErrorRatio,std::max(0.,::to_double(error[f*20+k]))/farTarget[f*20+k]);
            if(diag.ErrorRatio<=1||diag.Refinements==MaxRefinements)break;
            const auto parent=queue.top();queue.pop();update(parent,-1);const auto xy=quad(parent,8,12),yx=quad(parent,12,8);double ex=0,ey=0;
            for(int f=0;f<source.info_.fields;++f)for(int k=0;k<components;++k){const int at=f*20+k;ex=std::max(ex,::to_double(abs(parent.value[at]-xy[at]))/farTarget[at]);ey=std::max(ey,::to_double(abs(parent.value[at]-yx[at]))/farTarget[at]);}
            for(int side=0;side<2;++side){auto child=parent;if(ex>=ey){const double mid=(parent.x0+parent.x1)/2;if(side)child.x0=mid;else child.x1=mid;}
                else {const dd_real mid=(parent.p0+parent.p1)/2;if(side)child.p0=mid;else child.p1=mid;}assess(child);update(child,1);queue.push(std::move(child));}++diag.Refinements;}
        diag.QuadratureConverged=diag.ErrorRatio<=1;
        diag.LimitReached=diag.ErrorRatio>1&&diag.Refinements==MaxRefinements;
        std::array<double,4>scale;for(int d=0;d<=options.DerivativeOrder;++d)scale[d]=std::pow(source.info_.length,2-d);
        Values value{};for(int f=0;f<source.info_.fields;++f){std::array<dd_real,20>local{};for(int k=0;k<components;++k)local[k]=total[f*20+k];
            const auto global=Rotate(local,wideAxes,components);for(int k=0;k<components;++k)value[f*20+k]=::to_double((nearValue[f*20+k]+global[k])*scale[Degree(k)]);}
        for(double v:value)if(!std::isfinite(v))throw std::runtime_error("nonfinite MMSI result");if(diagnostics)*diagnostics=diag;return value;
    }
    Gravity GravityCartesian(const Vec&point,FastDiagnostics*diagnostics=nullptr)const {
        if(source_->info_.fields!=1)throw std::invalid_argument("gravity requires a scalar source");const auto f=Fields(point,diagnostics);Gravity out{};std::copy_n(f.begin(),20,out.begin());return out;
    }
    Magnetic MagneticCartesian(const Vec&point,FastDiagnostics*diagnostics=nullptr)const {
        if(source_->info_.fields!=3||source_->info_.options.DerivativeOrder<1)throw std::invalid_argument("magnetics requires a vector source and first derivatives");const auto v=Fields(point,diagnostics);Magnetic out{};
        const int components=source_->info_.options.DerivativeOrder==3?10:source_->info_.options.DerivativeOrder==2?4:1;
        for(int k=0;k<components;++k)for(int f=0;f<3;++f){const auto&e=detail::Axes[detail::GravityOrder[k]];int a[3]={e[0],e[1],e[2]};++a[f];out[k]+=v[f*20+detail::Component(a[0],a[1],a[2])];}return out;
    }
  private:
    FastShellIntegrator(std::shared_ptr<const detail::IndexedSource>source,double angle):source_(std::move(source)),angle_(angle*PI/180){}
};
}
namespace TSP_MMSI {
inline void ValidateTriangularMesh(const Model&model) {
    std::vector<std::uint64_t>edges;if(model.FullSphere)edges.reserve(3*model.faces.size());
    for(const auto&face:model.faces) {
        for(int id:face)if(id<0||std::size_t(id)>=model.vertices.size())throw std::invalid_argument("triangle index out of range");
        if(model.FullSphere)for(int i=0;i<3;++i){const auto edge=std::minmax(face[i],face[(i+1)%3]);
            edges.push_back((std::uint64_t(edge.first)<<32)|std::uint32_t(edge.second));}
    }
    std::sort(edges.begin(),edges.end());
    for(std::size_t i=0;i<edges.size();){std::size_t end=i+1;while(end<edges.size()&&edges[end]==edges[i])++end;
        if(end-i!=2)throw std::invalid_argument("full-sphere mesh has an open/nonmanifold edge");i=end;}
}
inline Model MakeTriangularModel(const std::vector<Vec>&vertices,const std::vector<std::array<int,3>>&faces,
                                  const std::vector<SourceValue>&values,int fields=1,bool fullSphere=true) {
    if(vertices.size()!=values.size())throw std::invalid_argument("model vertex count mismatch");
    if(fields!=1&&fields!=3)throw std::invalid_argument("source fields must be one or three");
    Model model;model.Fields=fields;model.FullSphere=fullSphere;
    model.vertices=vertices;model.faces=faces;model.data.resize(vertices.size()*(2+fields));
    for(std::size_t i=0;i<vertices.size();++i)model.SetState(i,values[i]);
    ValidateTriangularMesh(model);
    return model;
}
inline Model MakeGridModel(double Lons,double Lone,double Lats,double Late,double dLon,double dLat,
                           const double*R2,const double*R1,const std::array<const double*,3>&q,int fields=1) {
    using namespace detail;
    if(fields!=1&&fields!=3)throw std::invalid_argument("source fields must be one or three");
    if(!(dLon>0&&dLat>0)||!std::isfinite(dLon)||!std::isfinite(dLat))
        throw std::invalid_argument("invalid grid steps");
    const int nlon=static_cast<int>(std::llround((Lone-Lons)/dLon))+1;
    const int nlat=static_cast<int>(std::llround((Late-Lats)/dLat))+1;
    if(nlon<3||nlat<2||std::abs(nlon*dLon-360)>1e-7||
       std::abs(Lats-dLat/2+90)>1e-7||std::abs(Late+dLat/2-90)>1e-7)
        throw std::invalid_argument("legacy adapter requires a global cell-centred geographic grid");
    if(!R1||!R2||!q[0]||(fields==3&&(!q[1]||!q[2])))throw std::invalid_argument("null source array");
    std::vector<Vec>vertices;std::vector<SourceValue>values;std::vector<std::array<int,3>>faces;
    vertices.reserve(nlon*nlat+2);values.reserve(nlon*nlat+2);
    for(int i=0;i<nlat;++i)for(int j=0;j<nlon;++j) {
        const int old=j*nlat+i;const double lat=Lats+i*dLat,lon=Lons+j*dLon;
        SourceValue value{R1[old],R2[old],{}};
        if(fields==1)value.q[0]=q[0][old];
        else {const auto axes=LocalAxes(lat,lon);for(int a=0;a<3;++a)for(int b=0;b<3;++b)value.q[b]+=q[a][old]*axes[a][b];}
        vertices.push_back(Direction(lat,lon));values.push_back(value);
    }
    for(int i=0;i<nlat-1;++i)for(int j=0;j<nlon;++j) {
        const int a=i*nlon+j,b=i*nlon+(j+1)%nlon,c=(i+1)*nlon+j,d=(i+1)*nlon+(j+1)%nlon;
        faces.push_back({a,b,d});faces.push_back({a,d,c});
    }
    for(int pole=0;pole<2;++pole) {
        const int row=pole?nlat-1:0;SourceValue average=values[row*nlon];const auto first=average;
        for(int j=1;j<nlon;++j) {
            const auto&v=values[row*nlon+j];average.R1+=(v.R1-first.R1)/nlon;average.R2+=(v.R2-first.R2)/nlon;
            for(int k=0;k<fields;++k)average.q[k]+=(v.q[k]-first.q[k])/nlon;
        }
        const int at=static_cast<int>(vertices.size());vertices.push_back({0,0,pole?1.0:-1.0});values.push_back(average);
        for(int j=0;j<nlon;++j)faces.push_back({at,row*nlon+j,row*nlon+(j+1)%nlon});
    }
    return MakeTriangularModel(vertices,faces,values,fields);
}
}
namespace TSP_MMSI { namespace work_detail {using Integrator=FastShellIntegrator;}}
#include <atomic>
#include <exception>
#include <mutex>
#include <memory>
#include <cstdio>
#ifdef _OPENMP
#include <omp.h>
#endif
#ifdef TSP_MPI
#include <mpi.h>
#endif

namespace TSP_MMSI { namespace work_detail {
using Settings=Integrator::Settings;
inline void CheckFailure(std::exception_ptr failure,int mode) {
#ifdef TSP_MPI
    if(mode==2) {
        std::uint64_t failed=failure?1:0;tsp_detail::RankSum(&failed,1);
        if(failed&&!failure)throw std::runtime_error("another MPI rank failed");
    }
#else
    (void)mode;
#endif
    if(failure)std::rethrow_exception(failure);
}
inline void Broadcast(void*data,std::size_t bytes) {
#ifdef TSP_MPI
    auto*p=static_cast<char*>(data);
    while(bytes){const int count=static_cast<int>(std::min<std::size_t>(bytes,1u<<28));
        MPI_Bcast(p,count,MPI_BYTE,0,MPI_COMM_WORLD);p+=count;bytes-=count;}
#else
    (void)data;(void)bytes;
#endif
}
inline int Run(Model model,int n,const double*lat,const double*lon,const double*radius,
                const std::array<double*,20>&output,Settings settings,int mode,std::uint32_t outputMask=tsp_math::JET_FULL_MASK) {
    int rank=0,ranks=1;
#ifdef TSP_MPI
    if(mode==2){MPI_Comm_rank(MPI_COMM_WORLD,&rank);MPI_Comm_size(MPI_COMM_WORLD,&ranks);}
#else
    if(mode==2)throw std::runtime_error("compile with TSP_MPI");
#endif
    if(mode==2){Broadcast(&n,sizeof(n));Broadcast(&settings,sizeof(settings));Broadcast(&outputMask,sizeof(outputMask));}
    if(n<0)throw std::invalid_argument("negative observation count");
    int active_workers=1;
    std::vector<double>latitude(n),longitude(n),radii(n);
    std::exception_ptr failure;
    if(rank==0&&n) {
        try {
            if(!lat||!lon||!radius)throw std::invalid_argument("null observation array");
            if(!(outputMask&((1u<<(model.Fields==1?20:10))-1)))throw std::invalid_argument("empty output selection");
            for(int k=0;k<(model.Fields==1?20:10);++k)if((outputMask>>k&1u)&&!output[k])throw std::invalid_argument("null output array");
            std::copy_n(lat,n,latitude.begin());std::copy_n(lon,n,longitude.begin());std::copy_n(radius,n,radii.begin());
        }catch(...){failure=std::current_exception();}
    }
    CheckFailure(failure,mode);
    if(mode==2) {
        Broadcast(latitude.data(),n*sizeof(double));Broadcast(longitude.data(),n*sizeof(double));Broadcast(radii.data(),n*sizeof(double));

        int sourceFields=model.Fields;Broadcast(&sourceFields,sizeof(sourceFields));model.Fields=sourceFields;
    }
    const int fields=model.Fields==1?20:10;outputMask&=(1u<<fields)-1;
    int selected=0;for(int k=0;k<fields;++k)selected+=int(outputMask>>k&1u);
    int begin=0,end=n;
#ifdef TSP_MPI
    if(mode==2){const auto owned=tsp_detail::OwnedRange(n,rank,ranks);begin=owned.first;end=owned.second;}
    std::unique_ptr<tsp_detail::IndexedGather>gather;
    if(mode==2)gather.reset(new tsp_detail::IndexedGather(n,begin,end,fields==20?"MMSI-G":"MMSI-M"));
#endif
    const int rounds=mode==2?n/ranks+(n%ranks!=0):end-begin;
    std::array<std::vector<double>,20>values;for(int k=0;k<fields;++k)if(outputMask>>k&1u)values[k].resize(end-begin);
    std::mutex mutex;std::atomic<bool>failed{false};
    using QueryDiagnostics=Integrator::Diagnostic;
    std::atomic<std::uint64_t>quadrature_unconverged{0},near_approximate{0};
    if(rank==0){std::printf("%s | PREPARE building source index...\n",fields==20?"MMSI-G":"MMSI-M");std::fflush(stdout);}
    std::unique_ptr<Integrator>prepared;
    try{
#ifdef TSP_MPI
        if(mode==2)prepared=Integrator::PrepareMPI(std::move(model),settings,MPI_COMM_WORLD);
        else
#endif
        prepared.reset(new Integrator(std::move(model),settings));
    }
    catch(...){failure=std::current_exception();}
    CheckFailure(failure,mode);
    if(rank==0)std::printf("%s | READY source_bytes=%llu (%.2f MiB) node_ranks=%d model_copies_per_node=1\n",
        fields==20?"MMSI-G":"MMSI-M",static_cast<unsigned long long>(prepared->SourceBytes()),prepared->SourceBytes()/1048576.,prepared->SharedRanks());
    tsp_detail::Progress progress(fields==20?"MMSI-G":"MMSI-M",n,mode==2);
#ifdef _OPENMP
    #pragma omp parallel if(mode==1)
#endif
    {
#ifdef _OPENMP
    #pragma omp master
    active_workers=omp_get_num_threads();
#endif
    std::unique_ptr<Integrator>engine;
    try{engine.reset(new Integrator(*prepared));}
    catch(...){std::lock_guard<std::mutex>guard(mutex);if(!failure)failure=std::current_exception();failed.store(true);}
#ifdef _OPENMP
    #pragma omp for schedule(static)
#endif
    for(int position=0;position<rounds;++position) {
        const int i=begin+position;bool completed=false;
        if(i<end&&!failed.load())try {
            if(!std::isfinite(radii[i])||radii[i]<0)throw std::invalid_argument("invalid observation radius");
            const auto axes=detail::LocalAxes(latitude[i],longitude[i]);
            const Vec point=detail::Scale(axes[2],radii[i]);
            QueryDiagnostics diagnostic;
            if(fields==20) {
                const int components=settings.DerivativeOrder==3?20:settings.DerivativeOrder==2?10:settings.DerivativeOrder==1?4:1;
                const auto value=detail::Rotate(engine->GravityCartesian(point,&diagnostic),axes,components);
                for(int k=0;k<20;++k)if(outputMask>>k&1u)values[k][i-begin]=value[detail::GravityOrder[k]];
            } else {
                const int components=settings.DerivativeOrder==3?10:settings.DerivativeOrder==2?4:1;
                const auto value=detail::Rotate(engine->MagneticCartesian(point,&diagnostic),axes,components);
                for(int k=0;k<10;++k)if(outputMask>>k&1u)values[k][i-begin]=value[k];
            }
            if(!diagnostic.QuadratureConverged)++quadrature_unconverged;
            if(diagnostic.NearApproximation)++near_approximate;
            completed=true;
        }catch(...){std::lock_guard<std::mutex>guard(mutex);if(!failure)failure=std::current_exception();failed.store(true);}
        progress.Step(completed);
    }
    }
    progress.Finish();
    if(rank==0)tsp_detail::Stage(fields==20?"MMSI-G":"MMSI-M","SYNC","checking worker status...");
    CheckFailure(failure,mode);
    if(rank==0)tsp_detail::Stage(fields==20?"MMSI-G":"MMSI-M","CHECKS","reducing diagnostic counts...");
    std::uint64_t status[2]={quadrature_unconverged.load(),near_approximate.load()};
#ifdef TSP_MPI
    if(mode==2)tsp_detail::RankSum(status,2);
#endif
    if(rank==0){std::printf("%s | CHECKS near_approximate=%llu quadrature_unconverged=%llu observations=%d\n",
        fields==20?"MMSI-G":"MMSI-M",static_cast<unsigned long long>(status[1]),static_cast<unsigned long long>(status[0]),n);std::fflush(stdout);}
    if(mode==2) {
#ifdef TSP_MPI
        int component=0;for(int k=0;k<fields;++k)if(outputMask>>k&1u)gather->Collect(values[k],output[k],++component,selected);
#endif
    } else for(int k=0;k<fields;++k)if(outputMask>>k&1u) {
        std::copy(values[k].begin(),values[k].end(),output[k]);
    }
    if(rank==0)tsp_detail::Stage(fields==20?"MMSI-G":"MMSI-M","RELEASE","releasing source storage...");
    prepared.reset();
    if(rank==0)tsp_detail::Stage(fields==20?"MMSI-G":"MMSI-M","RELEASE","source storage released");
    return active_workers;
}
inline void Legacy(double Lons,double Lone,double Lats,double Late,double dLon,double dLat,
                   const double*R2,const double*R1,std::array<const double*,3>q,int fields,int n,
                   const double*lat,const double*lon,const double*r,double absTol,double relTol,
                   const std::array<double*,20>&out,Settings options,int mode) {
    int rank=0;
#ifdef TSP_MPI
    if(mode==2)MPI_Comm_rank(MPI_COMM_WORLD,&rank);
#endif
    double header[8]={Lons,Lone,Lats,Late,dLon,dLat,absTol,relTol};
    if(mode==2)Broadcast(header,sizeof(header));
    options.AbsoluteTolerance=header[6];options.RelativeTolerance=header[7];
    Model model;std::exception_ptr failure;
    if(rank==0)try{model=MakeGridModel(header[0],header[1],header[2],header[3],header[4],header[5],R2,R1,q,fields);}
        catch(...){failure=std::current_exception();}
    CheckFailure(failure,mode);
    Run(std::move(model),n,lat,lon,r,out,options,mode);
}
}}

#define TSP_MMSI_MASS_ARGS \
    double Lons,double Lone,double Lats,double Late,double dLon,double dLat,const double*R2,const double*R1,const double*rho, \
    int NumPoint,const double*Latitude,const double*Longitude,const double*Radius,double AbsTol,double RelTol, \
    double*V,double*Vx,double*Vy,double*Vz,double*Vxx,double*Vxy,double*Vyy,double*Vzx,double*Vzy,double*Vzz, \
    double*Vxxx,double*Vxxy,double*Vxxz,double*Vxyz,double*Vyyx,double*Vyyy,double*Vyyz,double*Vzzx,double*Vzzy,double*Vzzz, \
    const TSP_MMSI::work_detail::Settings& options={}
#define TSP_MMSI_MAG_ARGS \
    double Lons,double Lone,double Lats,double Late,double dLon,double dLat,const double*R2,const double*R1, \
    const double*Mx,const double*My,const double*Mz,int NumPoint,const double*Latitude,const double*Longitude,const double*Radius, \
    double AbsTol,double RelTol,double*V,double*Vx,double*Vy,double*Vz,double*Vxx,double*Vxy,double*Vyy,double*Vzx,double*Vzy,double*Vzz, \
    const TSP_MMSI::work_detail::Settings& options={}
#define TSP_MMSI_MASS(NAME,MODE) inline void NAME(TSP_MMSI_MASS_ARGS) { \
    TSP_MMSI::work_detail::Legacy(Lons,Lone,Lats,Late,dLon,dLat,R2,R1,{rho,nullptr,nullptr},1,NumPoint,Latitude,Longitude,Radius,AbsTol,RelTol, \
        {V,Vx,Vy,Vz,Vxx,Vxy,Vyy,Vzx,Vzy,Vzz,Vxxx,Vxxy,Vxxz,Vxyz,Vyyx,Vyyy,Vyyz,Vzzx,Vzzy,Vzzz},options,MODE);}
#define TSP_MMSI_MAG(NAME,MODE) inline void NAME(TSP_MMSI_MAG_ARGS) { \
    TSP_MMSI::work_detail::Legacy(Lons,Lone,Lats,Late,dLon,dLat,R2,R1,{Mx,My,Mz},3,NumPoint,Latitude,Longitude,Radius,AbsTol,RelTol, \
        {V,Vx,Vy,Vz,Vxx,Vxy,Vyy,Vzx,Vzy,Vzz},options,MODE);}
namespace TSP {
TSP_MMSI_MASS(MassShellIntegration,0)
TSP_MMSI_MASS(MassShellIntegrationOpenMP,1)
TSP_MMSI_MAG(MagnetizationShellIntegration,0)
TSP_MMSI_MAG(MagnetizationShellIntegrationOpenMP,1)
#ifdef TSP_MPI
TSP_MMSI_MASS(MassShellIntegrationMPI,2)
TSP_MMSI_MAG(MagnetizationShellIntegrationMPI,2)
#endif
}
#undef TSP_MMSI_MASS
#undef TSP_MMSI_MAG
#undef TSP_MMSI_MASS_ARGS
#undef TSP_MMSI_MAG_ARGS
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <limits>
#include <string>
#include <chrono>

namespace TSP_MMSI { namespace work_detail {
struct Input {
    Model model;
    Settings options;
    int points=0;
    std::array<std::vector<double>,3>coordinates;
    template<class T> static void Read(std::ifstream&in,T*values,std::size_t n) {
        if(n>static_cast<std::size_t>(std::numeric_limits<std::streamsize>::max())/sizeof(T))
            throw std::runtime_error("parameter count overflow");
        in.read(reinterpret_cast<char*>(values),static_cast<std::streamsize>(n*sizeof(T)));
        if(!in)throw std::runtime_error("truncated parameter file");
    }
    explicit Input(bool magnetic) {
        model.Fields=magnetic?3:1;
        std::ifstream in(magnetic?"MSIM.ForPar":"MSIG.ForPar",std::ios::binary);
        if(!in)throw std::runtime_error("cannot open shell parameter file");
        char magic[8];Read(in,magic,8);
        if(std::string(magic,7)=="TSPMSI1") {
            double tol[2];std::int32_t header[4];Read(in,tol,2);Read(in,header,4);
            options.AbsoluteTolerance=tol[0];options.RelativeTolerance=tol[1];
            const int nv=header[0],nf=header[1],fields=header[2];
            if(header[3]!=0)throw std::runtime_error("MMSI requires radial power zero");
            if(nv<4||nf<4||fields!=(magnetic?3:1)||nv>100000000||nf>200000000)
                throw std::runtime_error("invalid native mesh dimensions");
            model.vertices.resize(nv);model.faces.resize(nf);model.data.resize(std::size_t(nv)*(2+fields));std::vector<double>array(nv);
            Read(in,model.vertices.data(),nv);
            Read(in,array.data(),nv);for(int i=0;i<nv;++i)model.data[std::size_t(i)*(2+fields)+1]=array[i];
            Read(in,array.data(),nv);for(int i=0;i<nv;++i)model.data[std::size_t(i)*(2+fields)]=array[i];
            for(int f=0;f<fields;++f){Read(in,array.data(),nv);for(int i=0;i<nv;++i)model.data[std::size_t(i)*(2+fields)+2+f]=array[i];}
            Read(in,model.faces.data(),nf);ValidateTriangularMesh(model);
        } else {
            in.seekg(0);double h[8];std::int32_t shape[2];Read(in,h,8);Read(in,shape,2);
            options.AbsoluteTolerance=h[0];options.RelativeTolerance=h[1];
            for(double value:h)if(!std::isfinite(value))throw std::runtime_error("nonfinite geographic header");
            if(h[6]<=0||h[7]<=0)throw std::runtime_error("invalid geographic grid steps");
            const std::int64_t count=std::int64_t(shape[0])*shape[1];
            if(shape[0]<3||shape[1]<2||count>100000000)throw std::runtime_error("invalid geographic grid dimensions");
            const int nlon=static_cast<int>(std::llround((h[3]-h[2])/h[6]))+1;
            const int nlat=static_cast<int>(std::llround((h[5]-h[4])/h[7]))+1;
            if(nlon!=shape[0]||nlat!=shape[1])throw std::runtime_error("geographic dimensions disagree with grid steps");
            std::vector<double>r2(count),r1(count);Read(in,r2.data(),count);Read(in,r1.data(),count);
            std::array<std::vector<double>,3>q;std::array<const double*,3>ptr{};
            for(int f=0;f<model.Fields;++f){q[f].resize(count);Read(in,q[f].data(),count);ptr[f]=q[f].data();}
            model=MakeGridModel(h[2],h[3],h[4],h[5],h[6],h[7],r2.data(),r1.data(),ptr,model.Fields);
        }
        Read(in,&points,1);if(points<0||points>100000000)throw std::runtime_error("invalid observation count");
        for(auto&c:coordinates){c.resize(points);Read(in,c.data(),points);}
        if(in.peek()!=std::char_traits<char>::eof())throw std::runtime_error("extra bytes in parameter file");
    }
};
inline void ReportParameters(const Input&in,const Settings&settings,int kind,int mode,int ranks) {
    const auto&model=in.model;const bool magnetic=kind==1||kind==3;
    tsp_detail::ParameterReport report(magnetic?"MMSI-M":"MMSI-G",kind,mode,ranks);
    report.Line("Method","analytical near | joint adaptive 2-D far");
    report.Line("Input / output","%s -> %s",magnetic?"MSIM.ForPar":"MSIG.ForPar",magnetic?"ResultM.dat":"ResultG.dat");
    report.Line("Source mesh","vertices=%llu | triangles=%llu",static_cast<unsigned long long>(model.vertices.size()),static_cast<unsigned long long>(model.FaceCount()));
    report.Line("Observations","%d",in.points);
    report.Line("Length convention","consistent input units | native MATLAB files: m");
    const std::size_t stride=2+model.Fields,count=model.vertices.size();
    report.Range("Inner radius R1",model.data.data(),count,stride);
    report.Range("Outer radius R2",model.data.data()+1,count,stride);
    if(magnetic) {
        report.Line("Magnetization","global Cartesian components");
        for(int k=0;k<3;++k)report.Range(k==0?"Mx":k==1?"My":"Mz",model.data.data()+2+k,count,stride);
    }else report.Range("Density",model.data.data()+2,count,stride);
    report.Range("Longitude (deg)",in.coordinates[0].data(),in.points);
    report.Range("Latitude (deg)",in.coordinates[1].data(),in.points);
    report.Range("Observer radius",in.coordinates[2].data(),in.points);
    if(count&&in.points) {
        const double top=model.data[1];bool uniform=true;
        for(std::size_t k=1;k<count;++k)if(model.data[k*stride+1]!=top){uniform=false;break;}
        if(uniform){const auto bounds=std::minmax_element(in.coordinates[2].begin(),in.coordinates[2].end());
            report.Line("Height above top","%.12g .. %.12g (input units)",*bounds.first-top,*bounds.second-top);}
    }
    report.Line("Requested abs/rel","%.3e / %.3e",settings.AbsoluteTolerance,settings.RelativeTolerance);
    report.Line("Near angle (deg)","%s",std::getenv("TSP_MMSI_NEAR_DEGREES")?std::getenv("TSP_MMSI_NEAR_DEGREES"):"1");
    report.Line("Far refine limit","%llu",static_cast<unsigned long long>(Integrator::MaxRefinements));
    report.Line("Kernel node budget","%llu",static_cast<unsigned long long>(settings.MaxNodes));
    report.Line("Normalization","%s",magnetic?"no permeability factor applied":"G = 1; density included");
    report.End();
}
template<int Kind,int Mode> int Main(int argc,char**argv) {
    constexpr bool magnetic=Kind==1||Kind==3;constexpr int fields=magnetic?10:20;
    int rank=0,ranks=1;
#ifdef TSP_MPI
    if constexpr(Mode==2){MPI_Init(&argc,&argv);MPI_Comm_rank(MPI_COMM_WORLD,&rank);MPI_Comm_size(MPI_COMM_WORLD,&ranks);
        tsp_detail::MPICheck(MPI_Comm_set_errhandler(MPI_COMM_WORLD,MPI_ERRORS_RETURN),"MPI error handler");}
#else
    (void)argc;(void)argv;
#endif
    int code=0;
    try {
        std::unique_ptr<Input>input;if(rank==0)input.reset(new Input(magnetic));
        int count=rank==0?input->points:0;
        Settings settings;if(rank==0)settings=input->options;
        if(rank==0) {
        settings.DerivativeOrder=Kind==2?1:Kind==3?2:3;
        if(const char*p=std::getenv("TSP_MMSI_MAX_NODES"))settings.MaxNodes=std::stoull(p);
        ReportParameters(*input,settings,Kind,Mode,ranks);

        }
        constexpr std::uint32_t outputMask=Kind==2?8u:Kind==3?14u:(1u<<fields)-1;
        std::array<std::vector<double>,20>values;std::array<double*,20>out{};
        if(rank==0)for(int k=0;k<fields;++k)if(outputMask>>k&1u){values[k].resize(count);out[k]=values[k].data();}
        Model model;if(rank==0)model=std::move(input->model);
        const auto start=std::chrono::steady_clock::now();
        const int workers=Run(std::move(model),count,rank==0?input->coordinates[1].data():nullptr,
            rank==0?input->coordinates[0].data():nullptr,rank==0?input->coordinates[2].data():nullptr,out,settings,Mode,outputMask);
        const double seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
        if(rank==0) {
            std::vector<double>empty;if constexpr(Kind==2||Kind==3)empty.assign(count,std::numeric_limits<double>::quiet_NaN());
            const char*name=magnetic?"ResultM.dat":"ResultG.dat";
            tsp_detail::Stage(magnetic?"MMSI-M":"MMSI-G","WRITE","writing result file...");
            std::ofstream result(name,std::ios::binary);
            result.write(reinterpret_cast<const char*>(&count),4);
            for(const auto&v:input->coordinates)result.write(reinterpret_cast<const char*>(v.data()),count*std::streamsize(8));
            for(int k=0;k<fields;++k)result.write(reinterpret_cast<const char*>(outputMask>>k&1u?values[k].data():empty.data()),count*std::streamsize(8));
            result.close();
            if(!result)throw std::runtime_error("cannot write shell result");
            std::printf("%s | DONE kind=%d mode=%d observations=%d threads=%d elapsed_seconds=%.9f output=%s\n",magnetic?"MMSI-M":"MMSI-G",Kind,Mode,count,workers,seconds,name);
        }
    }catch(const std::exception&e){std::fprintf(stderr,"Rank %d: %s\n",rank,e.what());code=1;
#ifdef TSP_MPI
        if constexpr(Mode==2)MPI_Abort(MPI_COMM_WORLD,1);
#endif
    }
#ifdef TSP_MPI
    if constexpr(Mode==2){
        if(rank==0)tsp_detail::Stage(magnetic?"MMSI-M":"MMSI-G","FINALIZE","closing MPI...");
        MPI_Finalize();
        if(rank==0)tsp_detail::Stage(magnetic?"MMSI-M":"MMSI-G","FINISHED","MPI closed; program exiting");
    }
#endif
    return code;
}
}}
#endif
