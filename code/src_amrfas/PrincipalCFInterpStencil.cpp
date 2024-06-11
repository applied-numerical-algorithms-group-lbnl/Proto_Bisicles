#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "CFInterpMatrixLibrary.H"
#include "PrincipalCFInterpStencil.H"
#include "CFMatrixOldF_F.H"
#include "NamespaceHeader.H"
PrincipalCFInterpStencil::PrincipalCFInterpStencil()
{
  m_defined = false;
  // there are totally 2^D choices to flip the signs
  m_signFlippers.resize( 1<<SpaceDim );
  for (int i=0; i<m_signFlippers.size(); i++)
    {
      m_signFlippers[i] = fineIndexToVector(i,2);
      // map 0,1 to 1,-1, respectively.
      m_signFlippers[i] = m_signFlippers[i]*(-2)+1;
    }
}

void PrincipalCFInterpStencil::define(const IntVectSet& a_validArea,
                                      const Box& a_fineInterpBox,
                                      const int a_refineCoarse,
                                      const int a_polyDegree,
                                      const int a_numPolyCoefs)
{
  CH_TIME("PrincipalCFInterpStencil::define");
  CH_assert(a_validArea.numPts() > a_numPolyCoefs);
  CH_assert(a_validArea.contains(IntVect::Zero));

  m_fineInterpBox = a_fineInterpBox;
  m_refineCoarse = a_refineCoarse;
  // initialize stencil
  const Box baseFineBox(IntVect::Zero, (m_refineCoarse-1)*IntVect::Unit);
  m_coarseToFineFab.define(baseFineBox, a_numPolyCoefs);
  // the number of coarse indices is known.
  // pre-allocate the memory to avoid push_back.
  m_coarseBaseIndices.resize(SpaceDim*a_numPolyCoefs);
  int nFineCells = 1;
  for (int d=0; d<SpaceDim; d++)
    nFineCells *= m_refineCoarse;
  m_multiIndices.resize(nFineCells);
  for (int i=0; i<nFineCells; i++)
    m_multiIndices[i] = fineIndexToVector(i, m_refineCoarse);

  CFInterpMatrixLibrary::IMV quePtr = CFInterpMatrixLibrary::instance()
    ->getQueue(a_polyDegree, SpaceDim, a_refineCoarse);
  bool found = false;
  for (int i=0; i<quePtr->size(); i++)
    {
      found = tryStencil(a_validArea, quePtr->operator[](i));
      if (found) break;
    }
  if (!found)
    {
      pout()<< "No stencil fits into the following interpolation source:\n";
      a_validArea.printBoxes(pout());
      MayDay::Abort("\nAdd more stencils to the library!\n");
    }

  // Everything is defined now.
  m_defined = true;
}

const IntVect PrincipalCFInterpStencil::fineIndexToVector(int k, int r)
{
  IntVect mi;
  for (int d=0; d<SpaceDim; d++)
    {
      mi[d] = k%r;
      k = k/r;
    }
  return mi;
}

const int PrincipalCFInterpStencil::vectorToFineIndex(const IntVect& iv, int r)
{
  int id = 0;
  for (int d=SpaceDim-1; d>=0; d--)
    {
      id += iv[d];
      if (d>0) id *= r;
    }
  return id;
}

bool PrincipalCFInterpStencil::tryStencil
(const IntVectSet& a_validArea, Chombo::BaseCFInterpMatrix const*const a_mh)
{
  const int nP = a_mh->getParameter(Chombo::BaseCFInterpMatrix::NumPolynmCoefs);
  const IntVect shift = IntVect::Unit*static_cast<int>(m_refineCoarse/2);

  for (int i=0; i<m_signFlippers.size(); i++)
    {
      a_mh->fillLattice(m_stencil, m_signFlippers[i]);
      // don't bother if this stencil does not fit.
      if (!a_validArea.contains(m_stencil))
        continue;
      // now that a stencil is found, initialize members using it.
      for (int j=0; j<nP; j++)
        {
          for (int d=0; d<SpaceDim; d++)
            m_coarseBaseIndices[d+SpaceDim*j] =
              m_signFlippers[i][d]*a_mh->getLattice(j,d);
          for (int k=0; k<m_multiIndices.size(); k++)
            {
              const IntVect oldV = fineIndexToVector(k, m_refineCoarse);
              // In the d-th dimension, m_signFlippers[i][d] is 1 or -1;
              // For 1, newV[d] = oldV[d]
              // For -1, newV[d] = -oldV[d]+m_refineCoarse-1;
              const IntVect newV = oldV*m_signFlippers[i]
                + (m_signFlippers[i]-IntVect::Unit)*(1-m_refineCoarse)/2;
              m_coarseToFineFab(newV, j) =
                static_cast<double>(a_mh->getMatrix(j,k))/
                a_mh->getParameter(Chombo::BaseCFInterpMatrix::Denominator);
            }
        }
      return true;
    }
  // nothing fits.
  return false;
}

void PrincipalCFInterpStencil::fillFine(FArrayBox&       a_fineFab,
                                        const FArrayBox& a_coarseFab,
                                        const IntVect&   a_coarseDataCell,
                                        const IntVect&   a_coarseToFineOffset) const
{
  CH_TIME("PrincipalCFInterpStencil::fillFine");
  CH_assert(m_defined);
  IntVect fineBase = m_refineCoarse * a_coarseDataCell + a_coarseToFineOffset;
  FORT_APPLYCOARSEFINEINTERP2(CHF_FRA          (a_fineFab),
                             CHF_CONST_FRA    (a_coarseFab),
                             CHF_CONST_FRA    (m_coarseToFineFab),
                             CHF_CONST_INTVECT(fineBase),
                             CHF_CONST_INTVECT(a_coarseDataCell),
                             CHF_CONST_VI     (m_coarseBaseIndices),
                             CHF_BOX          (m_fineInterpBox));
}

const IntVectSet& PrincipalCFInterpStencil::getStencil(void) const
{
  return m_stencil;
}

// std::ostream& operator<<
// (std::ostream& os, const PrincipalCFInterpStencil& s)
// {
//   os << "Fitted Stencil = \n";
//   for (int k=0; k<15; k++)
//     {
//       for (int d=0; d<SpaceDim; d++)
//         os << s.m_coarseBaseIndices[k*SpaceDim+d] << ", ";
//       os << std::endl;
//     }
//   os << "Interp Matrix = \n";
//   for (int k=0; k<15; k++)
//     {
//       for (int fj=0; fj<2; fj++)
//         for (int fi=0; fi<2; fi++)
//           os << s.m_coarseToFineFab(IntVect(fi,fj),k) << ", ";
//       os << std::endl;
//     }
//   return os;
// }
#include "NamespaceFooter.H"
