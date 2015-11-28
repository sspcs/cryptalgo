#pragma once
#include "stdlib.h"
class Rho_pollard
{
private:
	static long long topow(long long a, long long b, long long n)// a to power b mod n
	{
		long long ans = 1;
		for (int i = 0; i < b; ++i)
		{
			ans = (ans*a) % n;
		}
		return ans;
	}
	bool is_gen(long long n);
	static inline long long s_indx(long long i, long long n) //returns index of the set s, to which i belongs
	{
		if (i >0 && i <= n / 3)
			return 1;
		if (i>n / 3 && i <= 2 * n / 3)
			return 2;
		return 3;
	}
	static long long gcd(long long a, long long b, long long & x, long long & y) {
		if (a == 0) {
			x = 0; y = 1;
			return b;
		}
		long long x1, y1;
		long long d = gcd(b%a, a, x1, y1);
		x = y1 - (b / a) * x1;
		y = x1;
		return d;
	}
	static bool rho_alg(long long a, long long b, long long c0, long long d0, long long n, long long &x);
	static long long pickrand(long long l, long long r)
	{
		return rand() % r + l;
	}
public:
	static long long dlog(long long a, long long b, long long n)
	{

		long long ans;
		if (rho_alg(a, b, 0, 0, n, ans))
		{
		}
		else
		{
			while (!rho_alg(a, b, pickrand(1, n - 1), pickrand(1, n - 1), n, ans))
			{
			}
		}
		return ans;
	};
};