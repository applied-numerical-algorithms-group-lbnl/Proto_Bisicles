#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _AMRFASPOISSONOP_H_
#define _AMRFASPOISSONOP_H_

#include "AMRFAS.H"
#include "LevelFluxRegister.H"
#include "CoarseFineInterp.H"
#include "NamespaceHeader.H"
#include "PrincipalCFInterpStencil.H"
///
/**
   Factory to create AMR FASPoissonOp(s) -- these two classes should be cloned to add operators
 */
class FASPoissonOpFactory: public AMRFAS_LDFOpFactory
{
public:
  FASPoissonOpFactory(int a_o) : AMRFAS_LDFOpFactory( 1, a_o )
  {
  }

  virtual ~FASPoissonOpFactory()
  {
  }

  ///
  /**
     a_bc holds the boundary conditions.
     a_alpha is the identity coefficient
     a_beta is the laplacian coefficient.
  */
  void define( BCHolder a_bc,
	       Real     a_alpha = 0.0,
	       Real     a_beta  = -1.0);

  ///
  virtual RefCountedPtr<AMRFASOp<LevelData<FArrayBox> > > AMRNewOp( int a_ilev, 
								    const DisjointBoxLayout &a_grid,
								    bool a_isSR = false
								    );  
  Real m_alpha;
  Real m_beta;
};

/**
   Poisson derived class for FAS operator
*/
class FASPoissonOp : public AMRFAS_LDFOp
{
  friend class FASPoissonOpFactory;
public:
  //! Constructor.
  FASPoissonOp( int a_o, const DisjointBoxLayout &a_grid );

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

protected:
  virtual void applyLevel( LevelData<FArrayBox>& a_LofPhi,
			   const LevelData<FArrayBox>& a_phi
			   );
public:

  virtual void reflux(const LevelData<FArrayBox>&        a_phiFine,
		      const LevelData<FArrayBox>&        a_phi,
                      LevelData<FArrayBox>&              a_residual,
		      AMRFASOp<LevelData<FArrayBox> >*   a_finerOp );

  // AMR stuff
  virtual void CFInterp( LevelData<FArrayBox>& a_phi,
			 const LevelData<FArrayBox>& a_phiCoarse);

protected:
  void getFlux(FArrayBox&       a_flux,
	       const FArrayBox& a_data,
	       int              a_dir,
	       const Box&       a_edgebox,
	       int              a_ref ) const;
  
  void getFlux( FArrayBox&       a_flux,
		const FArrayBox& a_data,
		int              a_dir,
		int              a_ref = 1 ) const;

  virtual void levelGSRB( RefCountedPtr<LevelData<FArrayBox> > a_phi,
			  const RefCountedPtr<LevelData<FArrayBox> > a_rhs
			  );
  
  virtual void levelRich( RefCountedPtr<LevelData<FArrayBox> >     a_phi,
			  const RefCountedPtr<LevelData<FArrayBox> > a_rhs
			  );
  
public:  
  Real m_alpha, m_beta;
  LevelFluxRegister       m_levfluxreg;
  CoarseFineInterp<PrincipalCFInterpStencil> m_interpWithCoarser;
  //QuadCFInterp m_interpWithCoarser;

  // SR stuff
  virtual void SR1_kernel( LevelData<FArrayBox>& a_phi, 
			   const LevelData<FArrayBox>& a_rhs, 
			   LevelData<FArrayBox>& a_cover,
			   int a_nsmooths,
			   int a_refrat,
			   int a_nGhosts,
			   int a_nBuffer );
  virtual void SR2_kernel( LevelData<FArrayBox>& a_phi, 
			   const LevelData<FArrayBox>& a_rhs, 
			   LevelData<FArrayBox>& a_cover,
			   int a_nsmooths,
			   int a_refrat,
			   int a_nGhosts,
			   int a_nBuffer );
  virtual void SR3_kernel( LevelData<FArrayBox>& a_phi, 
			   const LevelData<FArrayBox>& a_rhs, 
			   const LevelData<FArrayBox>& a_ru, 
			   LevelData<FArrayBox>& a_cover,
			   int a_nsmooths,
			   int a_refrat,
			   int a_nGhosts );
  virtual void SR4_kernel( LevelData<FArrayBox>& a_phi, 
			   const LevelData<FArrayBox>& a_rhs, 
			   const LevelData<FArrayBox>& a_ru, 
			   LevelData<FArrayBox>& a_cover,
			   int a_nsmooths,
			   int a_refrat,
			   int a_nGhosts,
			   int a_nBuffer, 
			   SRUserFunctional &a_func );
protected:
  virtual void smooth_sr( FArrayBox& phifab, const FArrayBox& rhsfab, Box phibox, int a_nsmooths );
};



#include "NamespaceFooter.H"

#endif
