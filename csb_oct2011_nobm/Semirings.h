
#ifndef _SEMIRINGS_H_
#define _SEMIRINGS_H_

#include <utility>
#include <climits>
#include <cmath>
#include <tr1/array>
#include "promote.h"

template <typename T>
struct inf_plus{
  T operator()(const T& a, const T& b) const {
	T inf = std::numeric_limits<T>::max();
    	if (a == inf || b == inf){
      		return inf;
    	}
    	return a + b;
  }
};

// (+,*) on scalars
template <class T1, class T2>
struct PTSR
{
	typedef typename promote_trait<T1,T2>::T_promote T_promote;

	static T_promote add(const T1 & arg1, const T2 & arg2)
	{
		return (static_cast<T_promote>(arg1) +  
			static_cast<T_promote>(arg2) );
	}
	static T_promote multiply(const T1 & arg1, const T2 & arg2)
	{
		return (static_cast<T_promote>(arg1) * 
			static_cast<T_promote>(arg2) );
	}
	// y += ax overload with a=1
	static void axpy(const T2 & x, T_promote & y)
	{
		y += x;
	}
	
	static void axpy(T1 a, const T2 & x, T_promote & y)
	{
		y += a*x;
	}
};

// (+,*) on std:array's
template<class T1, class T2, unsigned D>
struct PTSRArray
{
	typedef typename promote_trait<T1,T2>::T_promote T_promote;

	// y <- a*x + y overload with a=1
	static void axpy(const tr1::array<T2, D> & b, tr1::array<T_promote, D> & c)
	{
		for(int i=0; i<D; ++i)
		{
			c[i] +=  b[i];
		}	
	}
	
	static void axpy(T1 a, const tr1::array<T2,D> & b, tr1::array<T_promote,D> & c)
	{
		for(int i=0; i<D; ++i)
		{
			c[i] +=  a* b[i];
		}	
	}
};

// (min,+) on scalars
template <class T1, class T2>
struct MPSR
{
	typedef typename promote_trait<T1,T2>::T_promote T_promote;

	static T_promote add(const T1 & arg1, const T2 & arg2)
	{
		return std::min<T_promote> 
		(static_cast<T_promote>(arg1), static_cast<T_promote>(arg2));
	}
	static T_promote multiply(const T1 & arg1, const T2 & arg2)
	{
		return inf_plus< T_promote > 
		(static_cast<T_promote>(arg1), static_cast<T_promote>(arg2));
	}
};


#endif
