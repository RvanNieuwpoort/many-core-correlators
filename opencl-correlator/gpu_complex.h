#ifndef _GPU_COMPLEX_H
#define _GPU_COMPLEX_H

#include <iostream>

// #include <complex> does not work in -deviceemu mode
template <typename T> struct complex
{
    complex<T>() {}
    complex<T>(T r, T i = 0) : real(r), imag(i) {}
    complex<T> operator ~() const { return complex<T>(real, -imag); }
    complex<T> operator *(const complex<T> r) const { return complex<T>(real * r.real - imag * r.imag, real * r.imag + imag * r.real); }
    complex<T> operator += (const complex<T> r) { real += r.real, imag += r.imag; return *this;}
    bool operator != (const complex<T> r) const { return real != r.real || imag != r.imag; }

    T real, imag;
};

template <typename T> std::ostream &operator << (std::ostream &str, complex<T> &v)
{
    return str << '(' << v.real << ',' << v.imag << ')';
}

#endif // _GPU_COMPLEX_H
