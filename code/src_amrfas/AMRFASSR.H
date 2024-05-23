#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _AMRFASSR_H_
#define _AMRFASSR_H_

#include "AMRFAS.H"

#include "NamespaceHeader.H"


///
/**
   function interface for setting a function on a patch in SR.  For RHS vectors.
*/
typedef void(*SRFunc)( FArrayBox&           a_state,
		       const Box&           a_valid,
		       Real                 a_dx );


//! \class SRHolder
//! This class is a catch-all that can use a function pointer or a SRFunction
//! object to set function (ie, RHS) for SR
class SRHolder
{
public:

  //! Creates a SRHolder using a function pointer.
  SRHolder( SRFunc funcptr ) : m_funcptr(funcptr)
  {
  }

  //! Set a (RHS) function in the given FArrayBox. 
  void operator()( FArrayBox&           a_state,
		   const Box&           a_valid,
		   Real                 a_dx )
  {
    if (m_funcptr)
      {
	m_funcptr( a_state, a_valid, a_dx );
      }
    else
      MayDay::Error("no function ptr");
  }

protected:
  SRFunc m_funcptr;
};

class AMRFASSR : public AMRFAS<LevelData<FArrayBox> >
{
public:
  AMRFASSR( SRUserFunctional &a_func, int a_refRat = 2 ) :
    AMRFAS<LevelData<FArrayBox> >(),
    m_numSRLevels(1),
    m_numBufferCells0(2), // need 4 to keep HO prolongation
    m_numGhostCells(1),  // 2nd order, but this is infered in define()
    m_refRatio(a_refRat),
    m_fmg_pre(2),
    m_rhsFunc(0),
    m_userFunctional(a_func)
  {
  }

  virtual ~AMRFASSR()
  {
  }

  // need to creat phi & RHS LDF for SR levels, and fill in RHS for all levels
  virtual Real solve( Vector<RefCountedPtr<LevelData<FArrayBox> > > &a_phi, 
		      Vector<RefCountedPtr<LevelData<FArrayBox> > > &a_rhs,		      
		      int *a_status = 0
		      );

  // need to run MG in a special way with SR, could put Vcycle in Op to avoid this.
  virtual void FMG( Vector<RefCountedPtr<LevelData<FArrayBox> > > & a_phi,
		    const Vector<RefCountedPtr<LevelData<FArrayBox> > > & a_rhs, 
		    int = -1 );
  
  virtual Real VCycle( Vector<RefCountedPtr<LevelData<FArrayBox> > >& a_phi,
		       const Vector<RefCountedPtr<LevelData<FArrayBox> > >& a_rhs,
		       int, int );

  virtual void SRVCycle( Vector<RefCountedPtr<LevelData<FArrayBox> > >& a_phi,
			 const Vector<RefCountedPtr<LevelData<FArrayBox> > >& a_rhs,
			 int a_lev
			 );

  // need to add SR levels below, before AMRFAS adds levels above and creates 
  virtual void define( const ProblemDomain&       a_coarseDomain,
		       const RealVect&            a_crsDx,
		       const Vector<DisjointBoxLayout>& a_grids,
		       const Vector<int>&               a_refRatios,
		       AMRFASOpFactory<LevelData<FArrayBox> > &a_factory,
		       int a_numSR = 0,
		       int a_num_levels = -1 
		       );

  void setRHSFunction( SRHolder a_func )  { m_rhsFunc = a_func; }
  void setNumSRLevels( int a_n ) { m_numSRLevels = a_n;  }
  void setNumBufferCells( int a_n )  { m_numBufferCells0 = a_n;  }
  int getNumSRLevels()const{ return m_numSRLevels; }
  int getNumBufferCells(int a_sr_ilev)const{ 
    return m_numBufferCells0 + m_numSRLevels - a_sr_ilev; 
  }
  int getNumGhostCells()const{ return m_numGhostCells; }
  void copyUFromCover( LevelData<FArrayBox> &R_u_f,
		       const RefCountedPtr<LevelData<FArrayBox> > a_crsCover );
  //
  virtual void setParameters( const char *name = "solver" );
  
private:
  int m_numSRLevels;
  // "ghost" region for cover is m_numBufferCells/refRatio + nesting <= m_numBufferCells
  // We need 2 nesting cells for high order interpolation (of ghost cells).  m_numBufferCells needs to be 4
  int m_numBufferCells0;  // num burfffer cells on finest grid
  int m_numGhostCells;
  const int m_refRatio;
  int m_fmg_pre;
  SRHolder                   m_rhsFunc;
  SRUserFunctional          &m_userFunctional;
  Vector<Copier>             m_inCopier;    // copy into SR cover, before prolongation
  Vector<Copier>             m_outCopier;   // copy out from SR cover data, after restriction
  Vector<RefCountedPtr<LevelData<FArrayBox> > > m_crsCover;
  Vector<DisjointBoxLayout>  m_SRGrids;
};

#include "NamespaceFooter.H"

#endif
