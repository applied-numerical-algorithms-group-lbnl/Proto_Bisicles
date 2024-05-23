#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _CFINTERPMATRIX_H
#define _CFINTERPMATRIX_H

#include "Combinatorics.H"
#include "IntVectSet.H"
#include "FArrayBox.H"
#include "Vector.H"

///
///  This class organizes the Coarse-Fine interpolation matrix
///   that maps the known values to unknowns.
///

class BaseCFInterpMatrix
{
public:

  enum PARAMETER {
    PolynomialDegree=0,
    Dimensionality,
    RefineRatio,
    Denominator,
    MatrixNorm,
    NumPolynmCoefs,
    nParameters
  };

  virtual ~BaseCFInterpMatrix() {;}

  static const int s_MaxD = 3;

  virtual int getParameter(const PARAMETER& p) const=0;

  virtual int const*const getIndexStar(void) const=0;

  virtual int getLattice(int i, int d) const=0;

  Real getInfNorm(void) const
  {
    return 1.0*getParameter(MatrixNorm)/getParameter(Denominator);
  }

  virtual void fillLattice(IntVectSet&, const IntVect& a_scaler) const=0;

  virtual int getMatrix(int i, int j) const=0;

};


///
/// This class encapsulates a CF interp matrix for a single case
/// Because of static variable initialization dependency problem
///  and the lack of a general singleton class in Chombo,
///  this is not designed as a singleton class.
///
template<int p, int D, int r, int iStar0=p, int iStar1=p, int iStar2=p>
class CFInterpMatrix : public BaseCFInterpMatrix
{
public:
  typedef BaseCFInterpMatrix Base;

  static const int s_lattice[][D];
  static const int s_matrix[][Power<r,D>::val];
  // i* characterizes the available source area, see my paper.
  static const int s_indexStar[Base::s_MaxD];
  static const int s_int[Base::nParameters];

  virtual ~CFInterpMatrix() {;}

  int getParameter(const Base::PARAMETER& par) const
  {
    return s_int[par];
  }

  int const*const getIndexStar(void) const
  {
    return s_indexStar;
  }

  int getLattice(int i, int d) const
  {
    return s_lattice[i][d];
  }

  /// Multiply each multi-index in the lattice by a_scaler 
  /// and put them into a_lattice.
  void fillLattice(IntVectSet& a_lattice, const IntVect& a_scaler) const
  {
    a_lattice.makeEmpty();
    for (int j=0; j<s_int[Base::NumPolynmCoefs]; j++)
      a_lattice |= IntVect(s_lattice[j])*a_scaler;
  }

  int getMatrix(int i, int j) const
  {
    return s_matrix[i][j];
  }

};


/// re-organizing iStar? into an integer array.
/// only the first D components are valid.
template<int p, int D, int r, int i0, int i1, int i2>
const int
CFInterpMatrix<p,D,r,i0,i1,i2>::s_indexStar[] = {i0, i1, i2};


#endif
