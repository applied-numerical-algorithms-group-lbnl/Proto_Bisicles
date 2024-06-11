#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _AMRFAS_H_
#define _AMRFAS_H_

#include "REAL.H"
#include "parstream.H"
#include "CH_Timer.H"
#include "Copier.H"
#include "LevelDataOps.H"
#include "CFRegion.H"
#include "BCFunc.H"
#include "LoadBalance.H"
#include "ProblemDomain.H"
#include "AMRIO.H"

#include "NamespaceHeader.H"

enum FAS_PROLONG_type {PROL_CONST_1,PROL_LINEAR_2,PROL_QUAD_3,PROL_CUBIC_4,PROL_QUART_5};
enum FAS_RESTRICT_type {REST_CONST_1};
enum FAS_SMOOTHER_type {FAS_GSRB=0,FAS_RICH=1};

//! \class SRUserFunctional
//! base class for usser functional to be evaluated in SR (4) kernel
class SRUserFunctional
{
public:
  //! Set a (RHS) function in the given FArrayBox. 
  virtual void operator()( const FArrayBox&           a_state,
			   const Box&           a_valid,
			   Real                 a_dx )
  {
    // default does nothing.
    pout() << "SRUserFunctional: box " << a_valid << endl; 
  }
};

///
/**
   Base class for FAS operator
 */
template <class T>
class AMRFASOp 
{
  template <typename TT> friend class AMRFAS;
public:
  //! Constructor.
  AMRFASOp( int a_order, const DisjointBoxLayout& a_grid ) : 
    m_order(a_order),
    m_smoother(FAS_GSRB), 
    m_FMGProlOrderP(-1),
    m_ProlOrderP(-1),
    m_grid(a_grid)
  {
    // m_dx = -1.0; m_dxCrse = -1.0;
  }
  //! Destructor.
  virtual ~AMRFASOp()
  {
  }

  virtual void exchangeState(RefCountedPtr<T> a_phi,
			     Copier &a_copier)   = 0;
  ///
  virtual void restrictState( RefCountedPtr<AMRFASOp<T> >,
			      Copier &a_copier ) = 0; 
  virtual bool computeState( RefCountedPtr<T> a_phi, 
			     const RefCountedPtr<T> a_CrsPhi,
			     const RefCountedPtr<T> a_FinePhi 
			     ) = 0; // default to noop

  /** full define function for AMRLevelOp with both coarser and finer levels */
  void define(const DisjointBoxLayout&   a_grids,
              const DisjointBoxLayout&   a_gridsFiner,
              const DisjointBoxLayout&   a_gridsCoarser,
              const RealVect&            a_dxLevel,
              int                        a_refRatio,
              int                        a_refRatioFiner,
              const ProblemDomain&       a_domain,
              BCHolder                   a_bc,
              const Copier&              a_exchange,
              const CFRegion&            a_cfregion,
	      int ncomp,
	      bool a_isSR 
	      );

  /** full define function for AMRLevelOp with finer levels, but no coarser - coarse grid */
  void define(const DisjointBoxLayout&   a_grids,
              const DisjointBoxLayout&   a_gridsFiner,
              const RealVect&            a_dxLevel,
              int                        a_refRatio, // dummy arg, send in 1
              int                        a_refRatioFiner,
              const ProblemDomain&       a_domain,
              BCHolder                   a_bc,
              const Copier&              a_exchange,
              const CFRegion&            a_cfregion, 
	      int ncomp 
	      );
  
  ///
  /**
     define function for AMRLevelOp for non-caarsest grid
  */
  virtual void define(const DisjointBoxLayout&   a_grids,
		      const DisjointBoxLayout&   a_coarse,
		      const RealVect&            a_dxLevel,
		      int                        a_refRatio,
		      const ProblemDomain&       a_domain,
		      BCHolder                   a_bc,
		      const Copier&              a_exchange,
		      const CFRegion&            a_cfregion,
		      bool a_isSR 
		      );

  ///
  /**
     define function for AMRLevelOp which has no finer or coarser AMR level
  */
  void define(const RealVect&            a_dx,
              const ProblemDomain&       a_domain,
              BCHolder                   a_bc,
              const Copier&              a_exchange,
              const CFRegion&            a_cfregion
	      );

  /**
     Apply the AMR operator, and store results in a_p.
  */
  virtual void apply( T& a_p,
		      const T& a_phi,
		      const T *a_phiCoarse = 0,
		      bool a_doExchange = true
		      ) = 0;
protected:
  virtual void applyLevel( T& a_LofPhi,
			   const T& a_phi
			   ) = 0;
public:
  virtual void residual( T& a_resid,
			 const T& a_phi,
			 const T& a_rhs,
			 const T *a_phiCoarse = 0,
			 T *a_phiFine = 0,
			 AMRFASOp<T> *a_finerOp = 0 );
  ///
  /**
      Smoother. a_phi updated in place.
  */
  virtual void smooth( RefCountedPtr<T> a_phi,
		       const RefCountedPtr<T> a_rhs,
		       const T* a_phiCoarse,
		       int a_num_it = 2,
		       bool a_print = false );
protected:
  virtual void smooth_private( RefCountedPtr<T> a_phi,
			       const RefCountedPtr<T> a_rhs,
			       bool do_exchange = true );

  virtual void levelGSRB( RefCountedPtr<T> a_phi,
			  const  RefCountedPtr<T> a_rhs
			  ) = 0;
  
  virtual void levelRich( RefCountedPtr<T> a_phi,
			  const RefCountedPtr<T> a_rhs
			  ) = 0;
  
  virtual void reflux( const T&        a_phiFine,
		       const T&        a_phi,
		       T&              a_residual,
		       AMRFASOp<T>*   a_finerOp ) = 0;

public:

  virtual void write(const T* a_data, const char* a_filename) = 0;

  ///
  /**
     basic vector ops
  */
  
  virtual void create( T& a_lhs,
		       const T& a_rhs) = 0;

  virtual void assign(T&       a_lhs,
                      const T& a_rhs) = 0;
  
  virtual void assignLocal(T&       a_lhs,
                           const T& a_rhs) = 0;
    
  virtual void assignCopier(T&       a_lhs,
                            const T& a_rhs,
                            const Copier&               a_copier) = 0;
  
  virtual void zeroCovered(T& a_lhs,
                           T& a_rhs,
                           const Copier&         a_copier) = 0;
  
  virtual Real dotProduct(const T& a_1,
                          const T& a_2) = 0;
  /* multiple dot products (for GMRES) */
  virtual void mDotProduct(const T& a_1,
                           const int a_sz,
                           const T a_2[],
                           Real a_mdots[]) = 0;
  
  virtual void incr(T&       a_lhs,
                    const T& a_x,
                    Real                        a_scale) = 0;
  
  virtual void axby(T&       a_lhs,
                    const T& a_x,
                    const T& a_y,
                    Real                        a_a,
                    Real                        a_b) = 0;
  
  virtual void mult( T& a_x,
		     const T& a_y ) = 0;
  
  virtual void scale(T& a_lhs,
                     const Real&           a_scale) = 0;
  
  virtual Real norm(const T& a_x,
                    int   a_ord) const = 0;

  virtual Real norm( const T& a_x,
		     int   a_ord,
		     int a_comp ) const = 0;
  
   
  virtual Real localMaxNorm( const T& a_x ) const = 0;
  
  virtual void setToZero( T& a_x)  = 0;

  virtual int refToCoarser() const {return m_refToCoarser;}
  virtual int refToFiner() const {return m_refToFiner;}

  virtual void CFInterp( T& a_phi,
			 const T& a_phiCoarse) = 0;
  // R and P 
  virtual void AMRProlong( T&       a_fineU,
			   const T& a_CrsU,
			   T&       a_temp,
			   RefCountedPtr<AMRFASOp<T> > a_crsOp ) = 0;

  virtual void AMRFMGProlong( T&       a_fineU,
			      const T& a_CrsU,
			      T&       a_temp,
			      RefCountedPtr<AMRFASOp<T> > a_crsOp ) = 0;
  
  virtual void AMRRestrict( T&       a_CrsU,     
			    const T& a_fineU,
			    T& a_crsCover,
			    const Copier &a_copier
			    ) const = 0;

  int getOrder() const;
  Real dx() const {return m_dx[0];}

  // data
  Copier                  m_HOCopier;
  RealVect                m_dx;
  RealVect                m_dxCrse;
  ProblemDomain           m_domain;
  BCHolder                m_bc;
  //protected:
  CFRegion                m_cfregion;
  Copier                  m_exchangeCopier;
  int                     m_refToCoarser;
  int                     m_refToFiner;
  const int               m_order;
  FAS_SMOOTHER_type       m_smoother;
public:
  Real                    m_smoothing_damping_factor;
  int                     m_FMGProlOrderP;
  int                     m_ProlOrderP;
  const DisjointBoxLayout m_grid;

private:
  // Forbidden copiers.
  AMRFASOp(const AMRFASOp<T>&);
  AMRFASOp& operator=(const AMRFASOp<T>&);
public:
  // SR kernels, default to noop
  // SR1: SR2 + FMG prolongation
  virtual void SR1_kernel( T& a_phi, const T& a_rhs, T& a_cover, int, int, int, int ) {}
  // SR2: smooth, residual, restrict u, form A( R(u) ) + R( f - A(u) )
  virtual void SR2_kernel( T& a_phi, const T& a_rhs, T& a_cover, int, int, int, int ) {}
  // SR3: subtract R(u)apply from u_H; prolongate & increment u, smooth
  virtual void SR3_kernel( T&, const T&, const T&, T&, int, int, int ) {}
  // SR4: SR3 + functional
  virtual void SR4_kernel( T&, const T&, const T&, T&, int, int, int, int, 
			   SRUserFunctional &a_func ) {}
};

///
/**
   LevelData<FArrayBox> derived class for FAS operator
*/
class AMRFAS_LDFOp : public AMRFASOp<LevelData<FArrayBox> >
{
public:  
  virtual void exchangeState(RefCountedPtr< LevelData<FArrayBox> > a_phi,
			     Copier &a_copier) 
  {
    a_phi->exchange(a_copier);
  }
  //! Constructor.
  AMRFAS_LDFOp( int a_order, const DisjointBoxLayout &a_grid ) : 
    AMRFASOp<LevelData<FArrayBox> >( a_order, a_grid )
  {
  }
  /// - default to noop
  virtual void restrictState( RefCountedPtr<AMRFASOp<LevelData<FArrayBox> > >,
			      Copier &a_copier 
			      ) 
  {
  }

  virtual bool computeState( RefCountedPtr<LevelData<FArrayBox> > a_phi, 
			     const RefCountedPtr<LevelData<FArrayBox> > a_CrsPhi,
			     const RefCountedPtr<LevelData<FArrayBox> > a_FinePhi 
			     )
  {
    return false; // did I do an exchange on phi (side effect)
  }

  virtual void apply( LevelData<FArrayBox>& a_p,
		      const LevelData<FArrayBox>& a_phi,
		      const LevelData<FArrayBox> *a_phiCoarse = 0,
		      bool a_doExchange = true );

  // general vector ops, with T=LevelData<FArrayBox>

  LevelDataOps<FArrayBox> m_levelOps;
  
  virtual void create(LevelData<FArrayBox>& a_lhs,
		      const LevelData<FArrayBox>& a_rhs);

  virtual void assign(LevelData<FArrayBox>&       a_lhs,
                      const LevelData<FArrayBox>& a_rhs);

  virtual void assignLocal(LevelData<FArrayBox>&       a_lhs,
                           const LevelData<FArrayBox>& a_rhs);

  virtual void assignCopier(LevelData<FArrayBox>&       a_lhs,
                            const LevelData<FArrayBox>& a_rhs,
                            const Copier&               a_copier);

  virtual void zeroCovered(LevelData<FArrayBox>& a_lhs,
                           LevelData<FArrayBox>& a_rhs,
                           const Copier&         a_copier);

  virtual Real dotProduct(const LevelData<FArrayBox>& a_1,
                          const LevelData<FArrayBox>& a_2);
  /* multiple dot products (for GMRES) */
  virtual void mDotProduct(const LevelData<FArrayBox>& a_1,
                           const int a_sz,
                           const LevelData<FArrayBox> a_2[],
                           Real a_mdots[]);

  virtual void incr(LevelData<FArrayBox>&       a_lhs,
                    const LevelData<FArrayBox>& a_x,
                    Real                        a_scale);

  virtual void axby(LevelData<FArrayBox>&       a_lhs,
                    const LevelData<FArrayBox>& a_x,
                    const LevelData<FArrayBox>& a_y,
                    Real                        a_a,
                    Real                        a_b);
  
  virtual void mult( LevelData<FArrayBox>& a_x,
		     const LevelData<FArrayBox>& a_y );
  
  virtual void scale(LevelData<FArrayBox>& a_lhs,
                     const Real&           a_scale);

  virtual Real norm(const LevelData<FArrayBox>& a_x,
                    int                         a_ord) const;

  virtual Real norm( const LevelData<FArrayBox>& a_x,
		     int   a_ord,
		     int a_comp ) const;

  virtual Real localMaxNorm(const LevelData<FArrayBox>& a_x) const ;

  virtual void setToZero( LevelData<FArrayBox>& a_x) ;

  void write( const LevelData<FArrayBox>* a_data,
	      const char*                 a_filename);
  
  // R and P, with explict types, protected, the real system implementations
  virtual void AMRProlong(LevelData<FArrayBox>&       a_fineU,
			  const LevelData<FArrayBox>& a_CrsU,
			  LevelData<FArrayBox>&       a_temp,
			  RefCountedPtr<AMRFASOp<LevelData<FArrayBox> > > a_crsOp,
			  FAS_PROLONG_type a_type );
  virtual void AMRRestrict( LevelData<FArrayBox>& a_CrsU,     
			    const LevelData<FArrayBox>& a_fineU,
			    FAS_RESTRICT_type a_type ) const;

  // R and P, without explict types, API, call real methods with default params
public:
  virtual void AMRProlong(LevelData<FArrayBox>&       a_fineU,
			  const LevelData<FArrayBox>& a_CrsU,
			  LevelData<FArrayBox>&       a_temp,
			  RefCountedPtr<AMRFASOp<LevelData<FArrayBox> > > a_crsOp );

  virtual void AMRFMGProlong(LevelData<FArrayBox>&       a_fineU,
			  const LevelData<FArrayBox>& a_CrsU,
			  LevelData<FArrayBox>&       a_temp,
			  RefCountedPtr<AMRFASOp<LevelData<FArrayBox> > > a_crsOp );

  
  virtual void AMRRestrict( LevelData<FArrayBox>& a_CrsU,     
			    const LevelData<FArrayBox>& a_fineU,
			    LevelData<FArrayBox>& a_crsCover,
			    const Copier &a_copier ) const;
};

///
/**
   Base class Factory to create AMRFASOp
 */
template <class T>
class AMRFASOpFactory
{
public:
  AMRFASOpFactory( int a_nc, int a_o = 2 ) : 
    m_ncomp( a_nc ),
    m_order( a_o ),
    m_FMGProlOrderP(-1),
    m_ProlOrderP(-1)
  {
  }

  virtual ~AMRFASOpFactory()
  {
  }

  ///
  /**
     return a new operator.  this is done with a new call.
     caller is responsible for deletion
  */
  virtual RefCountedPtr<AMRFASOp<T> > AMRNewOp( int, const DisjointBoxLayout&, bool = false ) = 0;

  ///
  /**
     return refinement ratio  to next finer level.
  */
  virtual int refToFiner( const ProblemDomain& a_indexSpace ) const;
  ///
  /**
     a_coarseDomain is the domain at the coarsest level.
     a_grids is the AMR  hierarchy.
     a_refRatios are the refinement ratios between levels.  The ratio lives
         with the coarser level so a_refRatios[ilev] is the ratio between
         ilev and ilev+1
     a_coarseDx is the grid spacing at the coarsest level.
     a_bc holds the boundary conditions.
  */
  void define( BCHolder a_bc ); // called by app.
  virtual void define( const ProblemDomain& a_coarseDomain,  // called by AMRFAS
		       const RealVect&      a_crsDx,
		       const Vector<DisjointBoxLayout>& a_grids,
		       const Vector<int>&   a_refRatios,
		       int a_nSRGrids = 0
		       );
  int nComp()const{return m_ncomp;}
  int getOrder()const{return m_order;}

  /// 
protected:
  Vector<ProblemDomain>     m_domains;
  Vector<RealVect> m_dx;
public:
  Vector<DisjointBoxLayout> m_grids;
  Vector<int>  m_refRatios; // refinement to next coarser level
  BCHolder m_bc;

  Vector<Copier>   m_exchangeCopiers;
  Vector<CFRegion> m_cfregion;
protected:
  const int m_ncomp; // these are baked in, order could be delegated
  const int m_order;
public:
  Real      m_smoothing_damping_factor;
  int       m_FMGProlOrderP;
  int       m_ProlOrderP;
};

///
/**
   Factory to create AMRFAS_LDFOps
*/
class AMRFAS_LDFOpFactory: public AMRFASOpFactory<LevelData<FArrayBox> >
{
public:
  AMRFAS_LDFOpFactory( int a_nc = 1, int a_order = 2 ) : 
    AMRFASOpFactory<LevelData<FArrayBox> >( a_nc, a_order )
  {
  }				 
  virtual ~AMRFAS_LDFOpFactory()
  {
  }

protected:
  ///
  void AMRNewOp( const int, 
		 RefCountedPtr<AMRFAS_LDFOp>,
		 bool a_isSR = false
		 ); 
};

///
/**
   Class to solve nonlinear equations L(phi) = rho using FAS scheme on an AMR grid.
 */

enum FASMG_type {FAS_FULL=0,FAS_VCYCLE=1,FAS_FCYCLE=2,FAS_SMOOTH=3};

template <class T>
class AMRFAS
{
public:
  AMRFAS();
  virtual ~AMRFAS();

  // 
  virtual void setParameters( const char *name = "solver" );

  ///
  /**
     Define the solver.
     a_coarseDomain is the index space on the coarsest AMR level.
     a_factory is the operator factory through which all special information is conveyed.
     a_numLevels is the number of AMR levels.
  */
  virtual void define( const ProblemDomain&             a_coarseDomain,
		       const RealVect&                  a_crsDx,
		       const Vector<DisjointBoxLayout>& a_grids,
		       const Vector<int>&               a_refRatios,
		       AMRFASOpFactory<T > &            a_factory,
		       int a_numSR = 0,
		       int a_num_levels = -1
		       );
protected:
  virtual void init( const Vector<RefCountedPtr<T> > &a_phi, const Vector<RefCountedPtr<T> > &a_rhs );
public:
  ///
  /**
     Solve L(phi) = rhs.
  */
  virtual Real solve( Vector<RefCountedPtr<T> > &a_phi, 
		      Vector<RefCountedPtr<T> > &a_rhs,		      
		      int *a_status = 0
		      );

  Real recoverInitialResidual() const
  {
    return m_initial_rnorm;
  }
protected:
  virtual Real solveNoInit( Vector<RefCountedPtr<T> >& a_phi,
			    const Vector<RefCountedPtr<T> >& a_rhs, 
			    int *a_status );

  ///
  /**
     Apply a single V-cycle. Returns the max norm of the residual.
  */
  virtual Real VCycle( Vector<RefCountedPtr<T> >& a_phi,
		       const Vector<RefCountedPtr<T> >& a_rhs,
		       int, int );

  virtual void FMG( Vector<RefCountedPtr<T> > & a_phi,
		    const Vector<RefCountedPtr<T> > & a_rhs, 
		    int = -1 );

 ///
  /**
     resid = rhs - L(phi) 
  */
  Real computeAMRResidual_private( Vector<RefCountedPtr<T> >&       a_resid,
				   Vector<RefCountedPtr<T> >&       a_phi,
				   const Vector<RefCountedPtr<T> >& a_rhs, 
				   const int a_base = 0 );
public:
  Real computeAMRResidual( Vector<RefCountedPtr<T> >&       a_resid,
			   Vector<RefCountedPtr<T> >&       a_phi,
			   const Vector<RefCountedPtr<T> >& a_rhs );

  void setCycleType( FASMG_type a_type )
  {
    m_type = a_type;
  }
  FASMG_type getCycleType() const
  {
    return m_type;
  }

  void setSmootherType( FAS_SMOOTHER_type a_type )
  {
    m_smoother = a_type;
  }
  void setAvoidNorms( bool b = true )
  {
    m_avoid_norms = b;
  }
  void setNumVcycles( int n )
  {
    m_num_vcycles = n;
  }
  int getNumVcycles() const
  {
    return m_num_vcycles;
  }
  void setSmoothingDampingFactor( Real n )
  {
    m_smoothing_damping_factor = n;
  }
  Real getSmoothingDampingFactor() const
  {
    return m_smoothing_damping_factor;
  }
  
  void setFMGProlOrderP( int n )
  {
    m_FMGProlOrderP = n;
  }
  int getFMGProlOrderP() const
  {
    return m_FMGProlOrderP;
  }
  void setProlOrderP( int n )
  {
    m_ProlOrderP = n;
  }
  int getProlOrderP() const
  {
    return m_ProlOrderP;
  }

  void setInternalCoarseningRate( int r );

  // get the operator at a level using the user's indexing
  RefCountedPtr<AMRFASOp<T> > getOp(int a_userLevel )
  {
    return m_op[a_userLevel + m_numInternalLevels];
  }
  // zero covered using ussers index
  void zeroCovered( Vector<RefCountedPtr<T> > &a_phi, int a_userLevel )
  {
    int myIndex = a_userLevel + m_numInternalLevels;
    m_op[myIndex]->zeroCovered(*a_phi[a_userLevel],*m_rhsC[myIndex+1],m_restrictCopier[myIndex+1]);
  }

protected:

  virtual void clear();
  void setFixedBoxSize( const DisjointBoxLayout &a_grid );

  // data
  Vector<RefCountedPtr<AMRFASOp<T> > > m_op;
  Real m_initial_rnorm; // default = 4586.  set during solve.   Used to get around solve interface
  Vector<RefCountedPtr<T> >  m_residual;
  Vector<RefCountedPtr<T> >  m_temp;
  Vector<RefCountedPtr<T> >  m_rhsC;
  Vector<RefCountedPtr<T> >  m_solC;

  Vector<Copier> m_restrictCopier;

  //Vector<RefCountedPtr<AMRFASOp<T> > > m_opInternal;
  ///
  /**
     m_opInternal[k+m_numInternalLevels] = m_opUser[k], 0 <= k < m_numLevels;
  */
  int m_numInternalLevels;

  FASMG_type    m_type;
  FAS_SMOOTHER_type m_smoother;
  Real          m_smoothing_damping_factor;
  int           m_FMGProlOrderP;
  int           m_ProlOrderP;
  bool          m_avoid_norms;  // flag to avoid norms and residuals (for convergence checking)
  int           m_num_vcycles;
  ProblemDomain m_coarseDomain; // need to cache this
  int           m_privateCoarseningRate;
public:
  int           m_max_iter, m_min_iter, m_verbosity;
  Real          m_rtol, m_stagnate, m_atol;
  int           m_pre, m_post, m_num_cycles, m_coarse_its;
  bool          m_init;
  bool          m_plot_residual;
  int           m_fixedBoxSize;
private:
  // Forbidden copiers.
  AMRFAS(const AMRFAS<T>&);
  AMRFAS& operator=(const AMRFAS<T>&);
};

#include "NamespaceFooter.H"

#ifndef CH_EXPLICIT_TEMPLATES
#include "AMRFASI.H"
#endif // CH_EXPLICIT_TEMPLATES

#endif
