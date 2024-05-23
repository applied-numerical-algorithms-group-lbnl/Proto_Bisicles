
#define STATIC_assert(e) \
  typedef char (STATIC_assertion_failed) [(e)? 1 : -1]


/// Return the factorial n*(n-1)*...*(n-k+1)
template<int n, int k> struct Factorial
{
  STATIC_assert(k>0);
  enum{ val = n*Factorial<n-1,k-1>::val };
};

template<int n> struct Factorial<n,1>
{
  STATIC_assert(n>0);
  enum{ val = n };
};

template<int x, int p> struct Power
{
  STATIC_assert(p>0);
  enum{ val= x * Power<x,p-1>::val};
};
template<int x> struct Power<x,0>{enum { val=1 };};

template<int n, int k> struct BinomialCoefficient
{
  enum{ val = Factorial<n,k>::val/Factorial<k,k>::val };
};


/// Compute the number of coefficents in a polynomial at compile time 
template<int dim, int degree>
struct NumPolyCoefs
{
  /// sum up coefs of single degrees
  enum {
    val = NumPolyCoefs<dim,degree-1>::val
    + BinomialCoefficient<degree+dim-1, dim-1>::val
  };
};
    // the above calculation expands to the following
    //     c_numPolynomialCoefs  = 1 + SpaceDim
    //     + SpaceDim*(SpaceDim+1)/2 + SpaceDim*(SpaceDim+1)*(SpaceDim+2)/6
/// summation recursion terminator.
template<int dim> struct NumPolyCoefs<dim,0>{ enum {val=1}; };


template<int dim, int nCoef, int start=0>
struct HalfStencilWidth
{
  enum {
//     val = ((nCoef< Power<(1+2*start), dim>::val)? start :
//            HalfStencilWidth<dim,nCoef,start+1>::val)
// conditionals like the above will not work.
    val = (nCoef > Power<(1+2*start), dim>::val)
  };
};
