/*
Copyright (C) 2009 Rob van Nieuwpoort & John Romein
Astron
P.O.Box 2, 7990 AA Dwingeloo, The Netherlands, nieuwpoort@astron.nl

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/
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
