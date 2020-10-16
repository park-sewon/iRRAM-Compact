
// Euclidean space using alias declaration
template <int N>
using IR = std::array<REAL, N>;

// Hypercube as a (mental) subset
template <int N>
using HyperCube = IR<N>;

// just alias
template <int N>
using Point = IR<N>;


#define ZERO    RATIONAL(INTEGER(0),INTEGER(1))
#define ONE     RATIONAL(INTEGER(1),INTEGER(1))




// return 2^n
RATIONAL Exp(int n)
{
  if (n > 0)
    return (INTEGER(1) << n);
  else if (n == 0)
    return 1;
  else
    return (RATIONAL(INTEGER(1), INTEGER(1) << -n));
}

// Copied from iRRAM/src/stack.cc and modified a little.
// The original file is buggy in the case when the function is a constant function. (it never stops!)
// this should be reported sometime.
int module(std::function<REAL(REAL)> f, const REAL &x, int p)
{
  // Semantics: If m=module(f,x,p), then
  //    |x-z| <=2^m implies  |f(x)-f(z)| <= 2^p

  // If we are able to approximate f(x) by a DYADIC number d with an error of <=2^{p-1},
  // then for any z in the argument interval of x, f(z) must differ by at most
  //  2^{p-1} from d, hence |f(x)-f(z)|<=2^p

  int result;
  if ((ACTUAL_STACK.inlimit == 0) && iRRAM_thread_data_address->cache_i.get(result))
    return result;

  DYADIC d;
  REAL x_copy = x;

  sizetype argerror, testerror;

  x_copy.geterror(argerror);
  sizetype_set(testerror, 1, argerror.exponent);
  x_copy.adderror(testerror);
  {
    single_valued code;
    d = approx(f(x_copy), p - 1);
  }
  // At this line, we are sure that x_copy (and so also x) is precise enough to allow
  // the computation of f(x), even with a slightly increased error of the argument.

  // We now try to find the smallest p_arg such that the evaluation of f(x+- 2^p_arg)
  // is possible up to an error of at most  2^{p-1}

  ITERATION_STACK SAVED_STACK;

  // To do this, we start with p_arg=p.
  // If this is successfull, we increase the value of p_arg until the first failure
  // It it is not, then we decrease until the first success...
  int direction = 0, p_arg = p;
  bool try_it = true;

  while (try_it)
  {

    sizetype_set(testerror, 1, p_arg);
    x_copy.seterror(argerror);
    x_copy.adderror(testerror);
    bool fail = false;
    if (iRRAM_unlikely(iRRAM_debug > 0))
    {
      sizetype x_error;
      x_copy.geterror(x_error);
      iRRAM_DEBUG2(1, "Testing module: 1*2^%d + %d*2^%d\n", p_arg, argerror.mantissa, argerror.exponent);
      iRRAM_DEBUG2(1, "argument error: %d*2^%d\n", x_error.mantissa, x_error.exponent);
    }
    try
    {
      single_valued code;
      REAL z = f(x_copy);
      if (iRRAM_unlikely(iRRAM_debug > 0))
      {
        sizetype z_error;
        z.geterror(z_error);
        iRRAM_DEBUG2(1, "Module yields result %d*2^%d\n", z_error.mantissa, z_error.exponent);
      }
      d = approx(z, p - 1);
    }
    catch (Iteration it)
    {
      fail = true;
    }
    switch (direction)
    {
    case 0:;
      if (fail)
        direction = -1;
      else
        direction = 1;
      p_arg += direction;
      break;
    case 1:;
      if (fail || p_arg > 0)
      {
        try_it = false;
        p_arg -= direction;
      }
      else
      {
        p_arg += direction;
      };
      break;
    case -1:;
      if (fail)
      {
        p_arg += direction;
      }
      else
      {
        try_it = false;
      };
      break;
    }
  }
  iRRAM_DEBUG2(1, "Modules resulting in p_arg=%d\n", p_arg);

  sizetype_set(testerror, 1, p_arg);
  sizetype_inc(argerror, testerror);

  while (argerror.mantissa > 1)
  {
    argerror.mantissa = argerror.mantissa >> 1;
    argerror.exponent += 1;
  }

  result = argerror.exponent;
  if (ACTUAL_STACK.inlimit == 0)
    iRRAM_thread_data_address->cache_i.put(result);
  return result;
}


template <int N>
IR<N> IR_origin()
{
  std::array<REAL, N> res;
  for (int i = 0; i < N; i++)
    res = REAL(0);
  return res;
}

// metric
template <int N>
REAL IR_d(IR<N> x, IR<N> y)
{
  REAL sum = 0;
  REAL a, b;
  for (int i = 0; i < N; i++)
  {
    a = x[i];
    b = y[i];
    sum += (a - b) * (a - b);
  }
  return sqrt(sum);
}
