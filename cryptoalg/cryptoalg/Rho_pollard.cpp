#include "rho_pollard.h"
bool Rho_pollard::rho_alg(long long a, long long b, long long c0, long long d0, long long n, long long &x) {
	long long x0, x1, c1, d1;
	long long x20, x21, c20, c21, d20, d21;
	long long start_x0;
	x0 = 1;
	for (int i = 0; i < c0; ++i)
	{
		x0 = (x0*a) % n;
	}
	for (int i = 0; i < d0; ++i)
	{
		x0 = (x0*b) % n;
	}
	start_x0 = x0;
	x20 = 1;
	c20 = c0;
	d20 = d0;
	do
	{
		if (s_indx(x0, n) == 1)
		{
			x1 = (b*x0) % (n);
			c1 = c0;
			d1 = (d0 + 1) % (n - 1);
		}
		if (s_indx(x0, n) == 2)
		{
			x1 = (x0*x0) % n;
			c1 = (2 * c0) % (n - 1);
			d1 = (2 * d0) % (n - 1);
		}
		if (s_indx(x0, n) == 3)
		{
			x1 = (x0*a) % n;
			c1 = (c0 + 1) % (n - 1);
			d1 = d0;
		}
		x0 = x1;
		c0 = c1;
		d0 = d1;
		for (long long i = 0; i<2; ++i)
		{

			if (s_indx(x20, n) == 1)
			{
				x21 = (b*x20) % (n);
				c21 = c20;
				d21 = (d20 + 1) % (n - 1);
			}
			if (s_indx(x20, n) == 2)
			{
				x21 = (x20*x20) % n;
				c21 = (2 * c20) % (n - 1);
				d21 = (2 * d20) % (n - 1);
			}
			if (s_indx(x20, n) == 3)
			{
				x21 = (x20*a) % n;
				c21 = (c20 + 1) % (n - 1);
				d21 = d20;
			}
			x20 = x21;
			c20 = c21;
			d20 = d21;
			//	cout << x20 << " " << c20 << " " << d20 << " " << " ---" << endl;
		}

	} while (x0 != x20);
	//cout << x20 << " " << d0 << " " << d20 << " " << " ---" << endl;
	long long r = (((d0 - d20) % (n - 1)) + n - 1) % (n - 1);
	long long rinv;
	long long temp;
	if (gcd(r, n - 1, rinv, temp) > 100)
		return false;
	long long d = gcd(r, n - 1, rinv, temp);
	//	cout << "d is" << d << endl;
	long long subc = ((c20 - c0) % (n - 1) + n - 1) % (n - 1);
	//	cout << "rinv " << rinv << endl;
	start_x0 = ((rinv*subc) / d) % (n - 1);
	if (d > 1)
	{
		long long m = 0;
		x = start_x0;
		while (topow(a, x, n) != b)
		{
			++m;
			x = start_x0 + m*(n - 1) / d;
		}
		return true;

	}
	rinv = ((rinv % (n - 1)) + n - 1) % (n - 1);
	x = ((rinv*(c20 - c0)) % (n - 1) + n - 1) % (n - 1);
	return true;

}