#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _PRINCIPALCFINTERPSTENCIL_H
#define _PRINCIPALCFINTERPSTENCIL_H

#include "FArrayBox.H"
#include "IntVectSet.H"
#include "Vector.H"
#include "NamespaceHeader.H"


class BaseCFInterpMatrix;

///
/// Spatial interpolator based on polynomial fitting
///  intended to be passed to CoarseFineInterp as template parameter Stencil.
/// This contrasts to LeastSquareInterpStencil in that 
///  the number of cells in the stencil equals the number of polynomial coefs.
//   Qinghai Zhang  August 13 2008

class PrincipalCFInterpStencil
{
public:
  /// Default constructor
  /// define() has to be called before fillFine().
  PrincipalCFInterpStencil();

  /// Destructor
  virtual ~PrincipalCFInterpStencil(){}

  /// Actual constructor.
  void define(/// shifted Coarse cells as possible fitting candidates
              ///  with zero indexing the cell to be interpolated.
              const IntVectSet& a_validArea,
              /// fine cells within this coarse cell to be interpoloated,
              const Box& a_fineInterpBox,
              /// refinement ratio
              const int a_refineCoarse,
              /// max degree of polynomials
              const int a_polyDegree,
              /// number of coefficients of the polynomial.
              const int a_numPolyCoefs);

  const IntVectSet& getStencil(void) const;

  /// Interpolate to all the fine cells within one coarse cell.
  /**
     Fill the fine cells inside a_coarseDataCell with interpolated data
     from a_coarseFab.
     The fine cells to be filled in are at
     m_baseFineBox + m_refineCoarse * (a_coarseDataCell + a_coarseToFineOffset).
     The coarse cells from which we draw data are
     {a_coarseDataCell + m_coarseBaseIndices[i*D+[0:2]]:  0 <= i < m_stencilSize}.
   */
  void fillFine(/// to fill at m_baseFineBox+m_refineCoarse*(a_coarseBase + a_coarseToFineOffset)
                FArrayBox&         a_fineFab,
                /// coarse data
                const FArrayBox&   a_coarseFab,
                /// coarse cell with fine subcells to fill in
                const IntVect&     a_coarseDataCell,
                /// offset from base coarse cell to coarsened fine cells
                const IntVect&     a_coarseToFineOffset = IntVect::Zero) const;

  /// output the stencil and its interpolation matrix.
//   friend std::ostream& operator<<(std::ostream& os, 
//                                   const PrincipalCFInterpStencil& s);

  /// converting a scalar-index to a multi-index for a fine cell
  /// r is the refinement ratio.
  static const IntVect fineIndexToVector(int k, int r);

  /// converting a multi-index to a scalar-index
  static const int vectorToFineIndex(const IntVect& iv, int r);

protected:

  /// try to initialize the stencil using a_mh, return true for a success.
  /// The 2^D-1 symmetries of a_mh are also tried.
  bool tryStencil(const IntVectSet& a_validArea,
                  Chombo::BaseCFInterpMatrix const*const a_mh);

  /// whether define() has been called
  bool m_defined;

  /// refinement ratio
  int m_refineCoarse;

  /// subset of 0:(m_refineCoarse-1) in every dimension
  /// the fine position to be filled by fillFine();
  Box m_fineInterpBox;

  IntVectSet m_stencil;

  /// m_stencil.size() coefficients of coarse cells
  ///   for each cell in [0,m_refineCoarse-1]^SpaceDim
  FArrayBox m_coarseToFineFab;
  /// index list of the coarse cells, size is SpaceDim*NumPolyCoefs
  Vector<int> m_coarseBaseIndices;

  /// for flipping the signs of a stencil.
  Vector<IntVect> m_signFlippers;
  /// The vector indices for the fine cells to be interpolated.
  Vector<IntVect> m_multiIndices;

private:

  // Disallowed for all the usual reasons
  void operator=(const PrincipalCFInterpStencil& a_input)
  {
    MayDay::Error("invalid operator");
  }

  // Disallowed for all the usual reasons
  PrincipalCFInterpStencil(const PrincipalCFInterpStencil& a_input)
  {
    MayDay::Error("invalid operator");
  }
};
#include "NamespaceFooter.H"
#endif
