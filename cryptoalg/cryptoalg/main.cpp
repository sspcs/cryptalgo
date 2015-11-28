//#include "iprecision.h"
//#include "precisioncore.cpp"
#include "rho_pollard.h"
#include "stdlib.h"
#include <set>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <cmath>
#include <map>
#include <iterator>
#include "prime_tests.h"
#include "utilities.cpp"
using namespace std;
long long naive_discrete_log(long long b, long long g, long long n)
{
	long long i = 1;
	long long pow = g;
	while (i<n&&pow != b)
	{
		pow = (pow*g) % n;
		++i;
	}
	return i;
}
class Fermat_factor
{
public:
	static map<long long, long long> factor(long long n)
	{
		map<long long, long long> ans, a, b;
		if (PrimeTest::is_prime(n))
		{
			if (n>1)
				ans.insert(pair<int, int>(n, 1ll));
			return ans;
		}
		if (n % 2 == 0)
		{
			long long twopower = 0;
			while (n % 2ll == 0)
			{
				n /= 2ll;
				twopower++;
			}
			ans.insert(pair<int, int>(2ll, twopower));
		}
		long long x0, y0;
		x0 = round(sqrt(n));
		y0 = 0;
		int count = 0;
		while ((x0*x0 - y0*y0 - n) != 0)
		{
			count++;
			if (count>1000)
			{
				count = 0;
				//cout << x0 << " " << y0 << endl;

			}
			if (x0*x0 - y0*y0 - n>0)
				y0 = y0 + 1;
			if (x0*x0 - y0*y0 - n< 0)
				x0 = x0 + 1;
		}
	//	cout << x0 - y0 << " " << x0 + y0 << endl;
		a = factor(x0 - y0);
		b = factor(x0 + y0);
		merge_apply(a.begin(), a.end(), b.begin(), b.end(), inserter(ans, ans.begin()),
			compare_first<pair<long long, long long>>, sum_pairs<pair<long long, long long>>);
		return ans;
	}
};
class Polig_Hellman
{
public:
	static long long gcd(long long a, long long b, long long & x, long long & y) {
		if (a == 0ll) {
			x = 0ll; y = 1ll;
			return b;
		}
		long long x1, y1;
		long long d = gcd(b%a, a, x1, y1);
		x = y1 - (b / a) * x1;
		y = x1;
		return d;
	}
	static long long invmod(long long a, long long n)
	{
		long long ans,unused;
		gcd(a, n, ans, unused);
		cout << "gcd " << ans << endl;
		ans = (ans + n) % n;
		return ans;
	}
	static long long topow(long long a, long long b, long long n)
	{
		long long ans = 1;
		for (int i = 0; i < b; ++i)
		{
			ans = (ans*a) % n;
		}
		return ans;
	}
	void table(map<long long, long long> canon_decomp, map<long long, vector<long long> > tabl, long long n)
	{
		long long curprime;
		long long primepow;
		vector<long long> curvec;
		for (auto it = canon_decomp.begin(); it != canon_decomp.end(); ++it)
		{
			curvec.clear();
			curprime = (*it).first;
			primepow = (*it).second;
			for (int i = 0; i <= curprime - 1; ++i)
			{
				curvec.push_back(topow(curprime, i*(n - 1) / curprime,n));
			}
			tabl.insert(pair<long long, vector<long long> >(curprime, curvec));
		}
	}
	static long long logforpow(long long a, long long b, long long n, long long q,long long alpha)
	{
		vector<long long> curvec;
		long long curprime = q;
		long long primepow = alpha;
		for (int i = 0; i <= curprime - 1; ++i)
		{
			curvec.push_back(topow(a, i*(n - 1) / curprime, n));
		}
		for (int i = 0; i < curvec.size(); ++i)
		{
			cout << curvec[i] << " ";
		}
		cout << endl;
		long long xcur;
		vector<long long> vecofx;
		long long lastpow;
		long long apowbase=1;
		long long lastpowbase;
		lastpowbase = b;
		long long powdivider = q;
		long long powqinv;
		long long apow=1;
		lastpow =topow(b,(n-1)/q,n);
		cout << lastpow << endl;
		for (long long i = 0; i < alpha; ++i)
		{
			for (int j = 0; j < curvec.size(); ++j)
			{
				if ((lastpow+n)%n == (curvec[j]+n)%n)
				{
					xcur = j;
					vecofx.push_back(j);
				}
			}
			apow=(apow*(invmod(topow(a,vecofx[vecofx.size() - 1]*(powdivider/q), n),n)))%n;
			powdivider =powdivider* q;
			lastpow = (topow((b*apow)%n,(n-1)/powdivider,n))%n;
			cout << lastpow << endl;
		}
		cout << "IN" << invmod(4, 17);
		long long ans=0;
		long long qpow = 1;
		cout << "OLOLL" << endl;
		for (int i = 0; i < vecofx.size(); ++i)
		{
			cout << vecofx[i] << " ";
		}
		cout << endl;
		for (int i = 0; i < vecofx.size(); ++i)
		{
			ans += vecofx[i] * qpow;
			qpow *= q;
		}
		return ans;
	}
public:
	static long long dlog(long long a, long long b, long long n)
	{
		

	}
};
class Lenstra_factor
{
	static long long gcd(long long a, long long b, long long & x, long long & y) {
		if (a == 0ll) {
			x = 0ll; y = 1ll;
			return b;
		}
		long long x1, y1;
		long long d = gcd(b%a, a, x1, y1);
		x = y1 - (b / a) * x1;
		y = x1;
		return d;
	}
	static long long invmod(long long a, long long n)
	{
		long long ans, unused;
		gcd(a, n, ans, unused);
		//cout << "gcd " << ans << endl;
		ans = (ans + n) % n;
		return ans;
	}
	static long long lenalg(long long n, long long v, long long w, long long a, long long x, long long y)
	{
		long long k=1;
		//vector<long long> e;
		for (int r = 2; r <= w; ++r)
		{
			long long mr;
			long long rpow=1;
			mr = 0;
			for (; rpow <= v + 2 * round(sqrt(v)) + 1; ++mr)
			{
				rpow *= r;
			}
			if (PrimeTest::is_prime(r))
				k *= rpow;
		}
		long long unused1, unused2;
		long long xsum, ysum,x3,y3;
		long long lambda;
		xsum = x;
		ysum = y;
		long long d1,d,invxdif,nu;
		for (int i = 1; i < k; ++i)
		{
			d = gcd(xsum - x, n, invxdif, unused2);
			if (d > 1 && d < n)
				return d;
			if (d == 1)
			{
				lambda = ((ysum - y)*invxdif) % n;
				nu = (ysum - lambda*xsum) % n;
				x3 = ((-xsum - x + lambda*lambda) % n + n) % n;
				y3 = ((-lambda*x3 - nu) % n + n) % n;
				xsum = x3;
				ysum = y3;
			}
			else
			{
				d1 = gcd(ysum + y, n, invxdif, unused2);
				if (1 < d1&& d1 < n)
					return d1;
				if (d1 == n)
				{
					xsum = 0;
					ysum = 0;
				}
				if (d1 == 1)
				{
					lambda = (((3 * xsum*xsum + a)*invxdif) % n + n) % n;
					nu = ((ysum - lambda*xsum) % n + n) % n;
					x3 = ((-2 * xsum + lambda*lambda) % n + n) % n;
					y3 = ((-lambda*x3 - nu) % n + n) % n;
					xsum = x3;
					ysum = y3;

				}
			}
		}
		return -1;
	}
public:
	static long long factor(long long n)
	{
		if (n % 2 == 0)
			return 2;
		if (n % 3 == 0)
			return 3;
		long long x = rand() % n;
		long long y = rand() % n;
		long long a = rand() % n;
		long long v = 20;
		long long d;
		for(int w = 4; w < n; w+=10)
		{
			if (w % 10 == 0)
				v += 10;
			d = lenalg(n, v, w, a, x, y);
			if (d>1)
				return d;
			x = rand() % n;
			y = rand() % n;
			a = rand() % n;

		}
	}
};
class Index_calculus
{
	static vector<long long>  multrow(vector<long long> vec,long long alpha ,long long n)
	{
		for (int j = 0; j < vec.size(); ++j)
		{
			vec[j] = ((vec[j] * alpha) % n+n)%n;
		}
		return vec;
	}

	static vector<long long> subrows(vector<long long> a, vector <long long> b, long long n)
	{
		for (int i = 0; i < a.size(); ++i)
		{
			a[i] = ((a[i] - b[i]) % n + n) % n;
		}
		return a;
	}
	static long long topow(long long a, long long b, long long n)
	{
		long long ans = 1;
		for (int i = 0; i < b; ++i)
		{
			ans = (ans*a) % n;
		}
		return ans;
	}
	static vector<long long> basegen(long long r)
	{
		long long n;
		vector<long long> ans;
		vector<long long> fans;
		n = 10;
		do
		{
			n *= 2;
			ans.clear();
			vector<char> prime(n + 1, true);
			prime[0] = prime[1] = false;
			for (int i = 2; i <= n; ++i)
				if (prime[i])
					if (i * 1ll * i <= n)
						for (int j = i*i; j <= n; j += i)
							prime[j] = false;
			for (int i = 1; i <= n; ++i)
			{
				if (prime[i])
					ans.push_back(i);
			}
		} 
		while (ans.size() < r-1);
		for (int i = 0; i < r; ++i)
		{
			fans.push_back(ans[i]);
		}
		return fans;
	}
	static long long gcd(long long a, long long b, long long & x, long long & y) {
		if (a == 0ll) {
			x = 0ll; y = 1ll;
			return b;
		}
		long long x1, y1;
		long long d = gcd(b%a, a, x1, y1);
		x = y1 - (b / a) * x1;
		y = x1;
		return d;
	}
	public:
	static long long dlog(long long a, long long b, long long n)
	{
		long long r=3;
		long long c=1;// robust var
		vector<vector<long long> > relations;
		vector<long long> base= basegen(r);
		//base.insert(base.begin(), n - 1);
		//r++;
		vector<long long> fact;
		long long a0;
		long long temp;
		long long apow;
		long long inv;
		vector<long long> tempvec,difvec;
		
		long long unused;
		do
		{
			relations.clear();
			for (long long k = 1ll; relations.size() <r+c&&k<n - 1; ++k)
			{
				fact.clear();
				apow = topow(a, k, n);
				for (int i = 0; i < base.size(); ++i)
				{
					a0 = 0ll;
					while (apow%base[i] == 0)
					{
						apow /= base[i];
						a0++;
					}
					fact.push_back(a0);
				}
				if (apow == 1)
				{
					fact.push_back(k);
					for (int i = 0; i < relations.size(); ++i)
					{
						fact = subrows(fact, multrow(relations[i], fact[i], n - 1), n - 1);
					}
					if (gcd(fact[relations.size()], n - 1, inv, unused) == 1)
					{
						//fact.push_back(k);
						fact = multrow(fact, inv, n - 1);
						for (int i = 0; i < relations.size(); ++i)
						{
							relations[i] = subrows(relations[i], multrow(fact, relations[i][relations.size()], n - 1), n - 1);
						}
						relations.push_back(fact);
					}
				}
			}

		} while ((relations.size() < base.size()));
		vector<long long> constpowers;
		for (int i = 0; i < relations.size(); ++i)
		{
			constpowers.push_back(relations[i][relations[i].size()-1]);
		}
		long long k;
		long long ans;
		do
		{
			ans = 0;
			k = rand() % (n - 1);
			apow = topow(a, k, n);
			apow =(apow* b)%n;
			for (int i = 0; i < base.size(); ++i)
			{
				while (apow%base[i] == 0)
				{
					apow /= base[i];
					ans+=constpowers[i];
				}
			}
		} while (apow > 1);
		ans -= k;
		return (ans%(n-1)+n-1)%(n-1);
	}
};
class GroupGen
{
	static long long topow(long long a, long long b, long long n)
	{
		long long ans = 1;
		for (int i = 0; i < b; ++i)
		{
			ans = (ans*a) % n;
		}
		return ans;
	}
	static long long euler_func(long long n)
	{
		return n - 1;
	}
public:
	static long long findgen(long long n)//for prime n=order+1
	{
		map<long long, long long> canon_decomp;
		long long order = euler_func(n);
		canon_decomp = Fermat_factor::factor(order);
		long long b;
		long long g;
		b = 1;
		while(b==1)
		{
			g = rand() % n;
			for (auto it = canon_decomp.begin(); it != canon_decomp.end(); ++it)
			{
				b = topow(g, order / (*it).first, n);
				if (b == 1)
					break;
			}
		}
			return g;
	}
};
class Adleman_dl
{
	static vector<long long> basegen(long long boundary)
	{
		long long n;
		vector<long long> ans;
		vector<long long> fans;
		n = 10;
		do
		{
			n *= 2;
			ans.clear();
			vector<char> prime(n + 1, true);
			prime[0] = prime[1] = false;
			for (int i = 2; i <= n; ++i)
				if (prime[i])
					if (i * 1ll * i <= n)
						for (int j = i*i; j <= n; j += i)
							prime[j] = false;
			for (int i = 1; i <= n; ++i)
			{
				if (prime[i])
					ans.push_back(i);
			}
		} while (ans[ans.size()-1] < boundary);
		for (int i = 0; i < r; ++i)
		{
			if (ans[i] > boundary)
				break;
			fans.push_back(ans[i]);
		}
		return fans;
	}
public:
	static long long dlog(long long a, long long b, long long n)
	{
		double bound_const = 3;
		long long boundary = exp(bound_const*sqrt(log(n)*log(log(n))));
		vector<long long> base=basegen(boundary)
	}
};
int main()
{
	//cout << naive_discrete_log(15, 7, 41) << endl;;
	cout << Rho_pollard::dlog(10,3,47) << endl;
	//cout << PrimeTest::is_prime(2010) << endl;
	//map<long long,long long> ans=Fermat_factor::factor(63214);
	//cout << Polig_Hellman::logforpow(3, 11, 17, 2, 4);
	//cout << endl;
	//cout << Lenstra_factor::factor(11ll *7ll) << endl;
	//cout << Polig_Hellman::topow(10, 42, 47) << endl;
	//cout << GroupGen::findgen(139) << endl;
	cout<<Index_calculus::dlog(10, 17, 47);
	
	system("pause");
}
