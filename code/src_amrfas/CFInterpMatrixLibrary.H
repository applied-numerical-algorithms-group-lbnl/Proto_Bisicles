#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _CFINTERPMATRIXLIBRARY_H
#define _CFINTERPMATRIXLIBRARY_H

#include <queue>
#include <vector>

#include "CFInterpMatrix.H"
#include "RefCountedPtr.H"
#include "UsingNamespace.H"


/// This class organizes objects of CFInterpMatrix
///  using the `Singleton' design pattern.
class CFInterpMatrixLibrary
{
public:

  /// MatPtr is the pointer to CIM with a less-than symantics.
  struct MatPtr : public RefCountedPtr<BaseCFInterpMatrix>
  {
    MatPtr(BaseCFInterpMatrix* p=0) : RefCountedPtr<BaseCFInterpMatrix>(p)
    {}

    friend bool operator<(const MatPtr& m1, const MatPtr& m2)
    {
      if (m1.isNull() || m2.isNull())
        MayDay::Error("Comparing null pointers are illegal!");
      return (m1->getInfNorm() < m2->getInfNorm());
    }
  };
  typedef std::priority_queue<MatPtr >           Queue;
  typedef RefCountedPtr<Queue>                   QuePtr;
  typedef std::map<int, QuePtr >                 MatQueLib;
  typedef RefCountedPtr<std::vector<MatPtr > >   IMV;
  typedef std::map<int, IMV >                    MatVecLib;

  /// return the singleton object of this class.
  static CFInterpMatrixLibrary const*const instance()
  {
    if (!m_singleton)
      m_singleton = createObject();
    return m_singleton;
  }

  /// Given a three-tuple (p,D,r),
  ///  it returns a vector of CFInterpMatrix satisfying (p,D,r)
  const IMV getQueue(int p, int D, int r) const
  {
    MatVecLib::const_iterator iter = m_vecLib.find(hash(p,D,r));
    if ( iter ==m_vecLib.end() )
      MayDay::Abort("CFInterpMatrix not found in the library!");
    return iter->second;
  }

private:

  // Return the index of the three-tuple (p,D,r)
  // p and D must be less than 1024
  static int hash(int p, int D, int r)
  {
    return p + 1024*(D + 1024*r);
  }

  /// copy from queLib to vecLib with the matrices sorted.
  void sortMatrices(void)
  {
    if (!m_vecLib.empty())
      MayDay::Error("This routine has been called!");
    for (MatQueLib::iterator iter = m_queLib.begin();
         iter!=m_queLib.end(); iter++)
      {
        QuePtr& mque = iter->second;
        const int sz = (*mque).size();
        IMV mvec(new std::vector<MatPtr>(sz));
        for (int i=0; i<sz; i++)
          {
            (*mvec)[sz-1-i] = mque->top();
            //pout() << sz-1-i<<"-th norm: " << mque->top()->getInfNorm()<< std::endl;
            mque->pop();
          }
        m_vecLib[iter->first] = mvec;
        // debugging code
        // if (iter->first==hash(4,3,2))
        //   for (int i=0; i<sz; i++)
        //     {
        //       const MatPtr& tmp = mvec->operator[](i);
        //       pout() << "Adding lattice for (p=4,D=3,r=2) with i* = "
        //              << IntVect(tmp->getIndexStar())
        //              << "; norm = " << tmp->getInfNorm() << endl;
        //     }
      }
  }

  void add(BaseCFInterpMatrix* mp)
  {
    const int p = mp->getParameter(BaseCFInterpMatrix::PolynomialDegree);
    const int D = mp->getParameter(BaseCFInterpMatrix::Dimensionality);
    const int r = mp->getParameter(BaseCFInterpMatrix::RefineRatio);
    const int key = hash(p,D,r);
    MatQueLib::const_iterator iter = m_queLib.find(key);
    if ( iter ==m_queLib.end() )
      m_queLib[key] = QuePtr(new Queue);
    m_queLib[key]->push(MatPtr(mp));
  }

  static CFInterpMatrixLibrary* createObject();

  MatQueLib m_queLib;
  MatVecLib m_vecLib;

  static CFInterpMatrixLibrary* m_singleton;

  /// private constructors enforce the singleton pattern.
  CFInterpMatrixLibrary(){}

  CFInterpMatrixLibrary(CFInterpMatrixLibrary&);
};

#endif
