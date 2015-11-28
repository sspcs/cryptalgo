#include "prime_tests.h"
#include "prime_tests.h"
#include <algorithm>
#include <cmath>
bool PrimeTest::is_prime(long long a)
{
	int sr = sqrt(a);
	for (int i = 2; i <= sr; ++i)
	{
		if (a%i == 0)
			return false;
	}
	return true;
}
