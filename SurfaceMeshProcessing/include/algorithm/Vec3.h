#ifndef VEC3_H
#define VEC3_H

#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif

// Let gcc optimize conditional branches a bit better...
#ifndef likely
#  if !defined(__GNUC__) || (__GNUC__ == 2 && __GNUC_MINOR__ < 96)
#    define likely(x) (x)
#    define unlikely(x) (x)
#  else
#    define likely(x)   (__builtin_expect((x), 1))
#    define unlikely(x) (__builtin_expect((x), 0))
#  endif
#endif

#include <cstddef>
#include <cmath>

template <typename T> class Array3
{
protected:
	T v[3];

public:
	Array3() { fill((T)0); }
	explicit Array3(const T &t) { fill(t); }
	Array3(const T &x, const T &y, const T &z) { setValue(x, y, z); }
	Array3(const Array3<T> &rh) : Array3<T>(rh[0], rh[1], rh[2]) {}

	void fill(const T &t) { v[0] = t; v[1] = t; v[2] = t; }
	void setValue(const T &x, const T &y, const T &z) { v[0] = x; v[1] = y; v[2] = z; }
	T& operator [](int i) { return v[i]; }
	const T& operator [](int i) const { return v[i]; }
	T* data() { return v; }
	const T* data() const { return v; }

	Array3<T>& operator =(const Array3<T> &rh) { v[0] = rh[0]; v[1] = rh[1]; v[2] = rh[2]; return *this; }
};

template <typename T> class Vec3
{
protected:
	T v[3];

public:
	Vec3() { fill((T)0); }
	explicit Vec3(const T &t) { fill(t); }
	Vec3(const T &x, const T &y, const T &z) { setValue(x, y, z); }
	Vec3(const Vec3<T> &rh) : Vec3<T>(rh[0], rh[1], rh[2]) {}

	void fill(const T &t) { v[0] = t; v[1] = t; v[2] = t; }
	void setValue(const T &x, const T &y, const T &z) { v[0] = x; v[1] = y; v[2] = z; }
	T& operator [](int i) { return v[i]; }
	const T& operator [](int i) const { return v[i]; }
	T* data() { return v; }
	const T* data() const { return v; }

	Vec3<T>& operator =(const Vec3<T> &rh) { v[0] = rh[0]; v[1] = rh[1]; v[2] = rh[2]; return *this; }
	Vec3<T>& operator +=(const Vec3<T> &rh) { v[0] += rh[0]; v[1] += rh[1]; v[2] += rh[2]; return *this; }
	Vec3<T>& operator -=(const Vec3<T> &rh) { v[0] -= rh[0]; v[1] -= rh[1]; v[2] -= rh[2]; return *this; }
	Vec3<T>& operator *=(T r) { v[0] *= r; v[1] *= r; v[2] *= r; return *this; }
	Vec3<T>& operator /=(T r) { v[0] /= r; v[1] /= r; v[2] /= r; return *this; }

	Vec3<T> operator +(const Vec3<T> &rh) const { Vec3<T> rt(*this); rt += rh; return rt; }
	Vec3<T> operator -(const Vec3<T> &rh) const { Vec3<T> rt(*this); rt -= rh; return rt; }
	Vec3<T> operator *(T r) const { Vec3<T> rt(*this); rt *= r; return rt; }
	Vec3<T> operator /(T r) const { Vec3<T> rt(*this); rt /= r; return rt; }

	T operator |(const Vec3<T> &rh) const { return v[0] * rh[0] + v[1] * rh[1] + v[2] * rh[2]; }
	Vec3<T> operator % (const Vec3<T> &rhs) const
	{
		return Vec3<T>(v[1] * rhs[2] - v[2] * rhs[1],
			v[2] * rhs[0] - v[0] * rhs[2], v[0] * rhs[1] - v[1] * rhs[0]);
	}

	T sqrnorm() const { return (*this) | (*this); }
	T norm() const { return std::sqrt(sqrnorm()); }
	T normalize() { T l = norm(); (*this) /= l; return l; }
	T normalize_cond() { T l = norm(); if (likely(l > T(0))) (*this) /= l; return l; }
	Vec3<T> normalized() const { Vec3<T> rt(*this); rt.normalize(); return rt; }
	Vec3<T> normalized_cond() const { Vec3<T> rt(*this); rt.normalize_cond(); return rt; }

	T maxv() const { T m = v[0]; if (m < v[1]) m = v[1]; if (m < v[2]) m = v[2]; return m; }
	T minv() const { T m = v[0]; if (m > v[1]) m = v[1]; if (m > v[2]) m = v[2]; return m; }
	void maxv(const Vec3<T> &rh)
	{
		if (v[0] < rh[0]) v[0] = rh[0];
		if (v[1] < rh[1]) v[1] = rh[1];
		if (v[2] < rh[2]) v[2] = rh[2];
	}
	void minv(const Vec3<T> &rh)
	{
		if (v[0] > rh[0]) v[0] = rh[0];
		if (v[1] > rh[1]) v[1] = rh[1];
		if (v[2] > rh[2]) v[2] = rh[2];
	}

};

template < typename T > static inline void add(const Vec3<T> &a, const Vec3<T> &b, Vec3<T> &c)
{
	c[0] = a[0] + b[0];
	c[1] = a[1] + b[1];
	c[2] = a[2] + b[2];
}
template < typename T > static inline void addOn(const Vec3<T> &a, const Vec3<T> &b, Vec3<T> &c)
{
	c[0] += a[0] + b[0];
	c[1] += a[1] + b[1];
	c[2] += a[2] + b[2];
}
template < typename T > static inline void sub(const Vec3<T> &a, const Vec3<T> &b, Vec3<T> &c)
{
	c[0] = a[0] - b[0];
	c[1] = a[1] - b[1];
	c[2] = a[2] - b[2];
}
template < typename T > static inline void subOn(const Vec3<T> &a, const Vec3<T> &b, Vec3<T> &c)
{
	c[0] += a[0] - b[0];
	c[1] += a[1] - b[1];
	c[2] += a[2] - b[2];
}
template < typename T > static inline void multi(const Vec3<T> &a, T t, Vec3<T> &b)
{
	b[0] = a[0] * t;
	b[1] = a[1] * t;
	b[2] = a[2] * t;
}
template < typename T > static inline void diff(const Vec3<T> &a, T t, Vec3<T> &b)
{
	b[0] = a[0] / t;
	b[1] = a[1] / t;
	b[2] = a[2] / t;
}

using Vec3i = Vec3<int>;
using Vec3f = Vec3<float>;
using Vec3d = Vec3<double>;

static inline void cot(const Vec3d &a, Vec3d &b)
{
	using namespace ::std;
	b[0] = 1.0 / tan(a[0]);
	b[1] = 1.0 / tan(a[1]);
	b[2] = 1.0 / tan(a[2]);
}

static inline void multiAdd(const Vec3d &a, double t, Vec3d &b)
{
	b[0] += t * a[0];
	b[1] += t * a[1];
	b[2] += t * a[2];
}

static inline void multiSub(const Vec3d &a, double t, Vec3d &b)
{
	b[0] -= t * a[0];
	b[1] -= t * a[1];
	b[2] -= t * a[2];
}

#endif // !VEC3_H

