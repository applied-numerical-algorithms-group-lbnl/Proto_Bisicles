#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// MFA, Sept 3, 2012


#include "FORT_PROTO.H"
#include "AMRFASPoissonOp.H"

#include "AMRFASOpF_F.H"
#include "AverageF_F.H"
#include "ParmParse.H"

#include "NamespaceHeader.H"

// ---------------------------------------------------------
// CFInterp
// ---------------------------------------------------------
void 
FASPoissonOp::CFInterp( LevelData<FArrayBox>& a_phi,
			const LevelData<FArrayBox>& a_phiCoarse )
{
  if( m_interpWithCoarser.isDefined() ) // SR levels do not have this
    {
      //m_interpWithCoarser.coarseFineInterp( a_phi, a_phiCoarse );
      m_interpWithCoarser.interpolate( a_phi, &a_phiCoarse );
    }
}

void 
FASPoissonOp::reflux(const LevelData<FArrayBox>&        a_phiFine, // fine grid to get flux from
		     const LevelData<FArrayBox>&        a_phiCrs,  // this grid function to use?
		     LevelData<FArrayBox>&              a_residual,// output -- why is this different from a_phiCrs???
		     AMRFASOp<LevelData<FArrayBox> >*   a_finerOp )
{
  CH_TIMERS("AMRFAS_LDFOp::reflux");

  m_levfluxreg.setToZero();
  Interval interv(0,a_phiCrs.nComp()-1);

  CH_TIMER("AMRFAS_LDFOp::reflux:incrementCoarse", t2);
  CH_START(t2);

  // get current flux at this level
  DataIterator dit = a_phiCrs.dataIterator();
  for (dit.reset(); dit.ok(); ++dit)
    {
      const FArrayBox& coarfab = a_phiCrs[dit];
      
      if (m_levfluxreg.hasCF(dit()))
        {
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              FArrayBox coarflux; // should take out of loop for eff?
              getFlux( coarflux, coarfab, idir );
	      
              Real scale = 1.0;
              m_levfluxreg.incrementCoarse( coarflux, scale, dit(),
					    interv, interv, idir);
            }
        }
    }
  CH_STOP(t2);

  // cast OK I guess because we're changing ghost cells only
  LevelData<FArrayBox>& phiFineRef = (LevelData<FArrayBox>&)a_phiFine;
  AMRFAS_LDFOp* finerAMRPOp = dynamic_cast<AMRFAS_LDFOp*>( a_finerOp );
  CH_assert(finerAMRPOp);
  // fill in fines ghostcells
  finerAMRPOp->CFInterp( phiFineRef, a_phiCrs );

  IntVect phiGhost = phiFineRef.ghostVect();
  int ncomps = a_phiFine.nComp();

  CH_TIMER("AMRFAS_LDFOp::reflux::incrementFine", t3);
  CH_START(t3);

  // get fine flux
  DataIterator ditf = a_phiFine.dataIterator();
  const  DisjointBoxLayout& dblFine = a_phiFine.disjointBoxLayout();
  for (ditf.reset(); ditf.ok(); ++ditf)
    {
      const FArrayBox& phifFab = a_phiFine[ditf];
      const Box& gridbox = dblFine[ditf];

      for (int idir = 0; idir < SpaceDim; idir++)
        {
          SideIterator sit;
          for (sit.begin(); sit.ok(); sit.next())
            {
              if (m_levfluxreg.hasCF(ditf(), sit()))
                {
                  Side::LoHiSide hiorlo = sit();
                  Box fluxBox = bdryBox(gridbox,idir,hiorlo,1);

                  FArrayBox fineflux(fluxBox,ncomps);
                  getFlux(fineflux, phifFab, idir, fluxBox, m_refToFiner);

                  Real scale = 1.0;
                  m_levfluxreg.incrementFine(fineflux, scale, ditf(),
                                             interv, interv, idir, hiorlo);
                }
            }
        }
    }

  CH_STOP(t3);

  // subtract off contents of flux registers: cellValue -= refluxScale*registerContents
  Real scale = 1.0/m_dx[0]; // assume isotropic mesh
  m_levfluxreg.reflux( a_residual, scale );
}


///
/**
   Poisson Op class -- these two classes should be cloned to add operators
*/
FASPoissonOp::FASPoissonOp( int a_o, const DisjointBoxLayout &a_grid ) : AMRFAS_LDFOp( a_o, a_grid )
{
  m_alpha = 0.0;
  m_beta  = -1.0;
}

void
FASPoissonOp::define( const DisjointBoxLayout& a_grids,
		      const DisjointBoxLayout& a_coarse,
		      const RealVect&          a_dxLevel,
		      int                      a_refRatioCrs,
		      const ProblemDomain&     a_domain,
		      BCHolder                 a_bc,
		      const Copier&            a_exchange,
		      const CFRegion&          a_cfregion,
		      bool a_isSR  )
{
  CH_TIME("FASPoissonOp::define");
  
  AMRFASOp<LevelData<FArrayBox> >::define( a_grids, a_coarse, a_dxLevel, a_refRatioCrs, a_domain,
					   a_bc, a_exchange, a_cfregion, a_isSR );

  if( !a_isSR )
    {
      // sets ghost cells with coarse grid data
      const int nGhost = m_order/2; // hack for number of ghosts
      m_interpWithCoarser.define( a_grids,   // fine (this) grid
				  a_coarse,  // coarse grid
				  m_domain,  // this domain
				  false,     // homogeneneous
				  m_order,   // degree of the fitting polynomial
				  a_refRatioCrs, // refinement ratio between this level and the coarser level
				  nGhost,     // number of layers of ghost cells to fill by interpolation
				  (m_order)/a_refRatioCrs+1, // proper nesting width
				  false,     // Should corner ghosts be interpolated?
				  true,      // if this is false, interpolate the whole ghosted fine patch
				  true       // if this is true, abort when there is not enough cells 
				  );

    }
}

// ---------------------------------------------------------
// applyLevel - apply operator on one level - do BC's but no exchange or C-F
// ---------------------------------------------------------
void FASPoissonOp::applyLevel( LevelData<FArrayBox>& a_lhs,
			       const LevelData<FArrayBox>& a_phi )
{
  CH_TIME("FASPoissonOp::applyLevel");

  LevelData<FArrayBox>& phi = (LevelData<FArrayBox>&)a_phi; // fortran chokes on const

  const DisjointBoxLayout& dbl = a_lhs.disjointBoxLayout();
  DataIterator dit = a_phi.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      FArrayBox& phiFab = phi[dit];
      const Box& valid = dbl[dit];
      FArrayBox& lhsFab = a_lhs[dit];

      m_bc( phiFab, valid, m_domain, m_dx[0], true );
      
      if( m_order == 2 )
	{	  
	  FORT_OPERATORLAP_FAS(CHF_FRA(lhsFab),
			       CHF_CONST_FRA(phiFab),
			       CHF_BOX(valid),
			       CHF_CONST_REAL(m_dx[0]),
			       CHF_CONST_REAL(m_alpha),
			       CHF_CONST_REAL(m_beta));
	}
      else  if( m_order == 4 )
	{
	  FORT_OPERATORLAP4(CHF_FRA(lhsFab),
			    CHF_CONST_FRA(phiFab),
			    CHF_BOX(valid),
			    CHF_CONST_REAL(m_dx[0]),
			    CHF_CONST_REAL(m_alpha),
			    CHF_CONST_REAL(m_beta));
	}
      else
	{
	  MayDay::Abort("FASPoissonOp::applyLevel bad order");
	}
    }
}

// ---------------------------------------------------------
void 
FASPoissonOp::getFlux( FArrayBox&       a_flux,
		       const FArrayBox& a_data,
		       int              a_dir,
		       int              a_ref /* =1 */) const
{
  CH_TIME("AMRFAS_LDFOpp::getFlux");
  
  CH_assert(a_dir >= 0);
  CH_assert(a_dir <  SpaceDim);
  CH_assert(!a_data.box().isEmpty());
  
  int nGhost = m_order/2;
  Box edgebox = surroundingNodes( a_data.box(), a_dir ); // increases upper corner of box 1 in a_dir direction
  edgebox.grow(a_dir, -nGhost);
  // if this fails, the data box was too small (one cell wide, in fact)
  CH_assert( !edgebox.isEmpty() );

  a_flux.resize( edgebox, a_data.nComp() );

  getFlux( a_flux, a_data, a_dir, edgebox, a_ref );
}

// ---------------------------------------------------------
// levelGSRB
// ---------------------------------------------------------
void FASPoissonOp::levelGSRB( RefCountedPtr<LevelData<FArrayBox> > a_phi,
			      const RefCountedPtr<LevelData<FArrayBox> > a_rhs )
{
  CH_TIME("FASPoissonOp::levelGSRB");
  
  CH_assert(a_phi->isDefined());
  CH_assert(a_rhs->isDefined());
  CH_assert(a_phi->ghostVect() >= IntVect::Unit);
  CH_assert(a_phi->nComp() == a_rhs->nComp());

  Real omega = m_smoothing_damping_factor;

  const DisjointBoxLayout& dbl = a_rhs->disjointBoxLayout();

  DataIterator dit = a_rhs->dataIterator();
  
  // do first red, then black passes
  bool do_exchange = false;
  for (int wpass = 0; wpass <= 1; wpass++ )
    {
      CH_TIME("FASPoissonOp::levelGSRB::Compute");
      
      if( do_exchange ){
        CH_TIME("FASPoissonOp::levelGSRB::exchange");
	a_phi->exchange( m_exchangeCopier );
      }
      do_exchange = true; // just cheat on first iteration

      for (dit.begin(); dit.ok(); ++dit)
        {
          const Box& valid = dbl[dit];
	  FArrayBox& phiFab = (*a_phi)[dit];

	  // BCs
	  m_bc( phiFab, valid, m_domain, m_dx[0], true );

	  if( m_order == 2 )
	    {
	      if (m_alpha == 0.0 && m_beta == -1.0 )
		{
		  FORT_GSRBLAPLACIAN_FAS(CHF_FRA(phiFab),
					 CHF_CONST_FRA((*a_rhs)[dit]),
					 CHF_BOX(valid),
					 CHF_CONST_REAL(m_dx[0]),
					 CHF_CONST_INT(wpass),
					 CHF_CONST_REAL(omega)
					 );
		}
	      else
		{
		  FORT_GSRBHELMHOLTZ_FAS(CHF_FRA(phiFab),
					 CHF_CONST_FRA((*a_rhs)[dit]),
					 CHF_BOX(valid),
					 CHF_CONST_REAL(m_dx[0]),
					 CHF_CONST_REAL(m_alpha),
					 CHF_CONST_REAL(m_beta),
					 CHF_CONST_INT(wpass),
					 CHF_CONST_REAL(omega)
					 );
		}
	    }
	  else if( m_order == 4 )
	    {
	      FArrayBox tmp(valid,1);
	      FORT_GSRBHELMHOLTZ4_FAS(CHF_FRA(phiFab),
				      CHF_CONST_FRA((*a_rhs)[dit]),
				      CHF_BOX(valid),
				      CHF_CONST_REAL(m_dx[0]),
				      CHF_FRA(phiFab),				      
				      CHF_CONST_REAL(m_alpha),
				      CHF_CONST_REAL(m_beta),
				      CHF_CONST_INT(wpass),
				      CHF_CONST_REAL(omega));
	    }
	  else
	    {
	      MayDay::Abort("FASPoissonOp::levelGSRB bad order");
	    }

        } // end loop through grids
    } // end loop through red-black

}

// ---------------------------------------------------------
// levelRich
// ---------------------------------------------------------
void FASPoissonOp::levelRich( RefCountedPtr<LevelData<FArrayBox> > a_phi,
			      const RefCountedPtr<LevelData<FArrayBox> > a_rhs )
{
  CH_TIME("FASPoissonOp::levelRich");
  
  CH_assert(a_phi->isDefined());
  CH_assert(a_rhs->isDefined());
  CH_assert(a_phi->ghostVect() >= IntVect::Unit);
  CH_assert(a_phi->nComp() == a_rhs->nComp());

  // BCs
  {
    Real diag = m_alpha - m_beta*2.*SpaceDim/(m_dx[0]*m_dx[0]);
    Real omega = m_smoothing_damping_factor;

    LevelData<FArrayBox> tmp;
    create( tmp, *a_rhs );

    apply( tmp, *a_phi, 0, false ); // applies BCs, excange done above

    // x = x + omega D^-1 (b - Ax) -- D == alaph - 2*beta*Dim/h^2
    axby( tmp, tmp, *a_rhs, -1.0, 1.0 );
    axby( *a_phi, *a_phi, tmp, 1.0, omega/diag );
  }

}

// ---------------------------------------------------------
void FASPoissonOp::getFlux(FArrayBox&       a_flux,
			   const FArrayBox& a_data,
			   int              a_dir,
			   const Box&       a_edgebox,
			   int              a_ref ) const
{
  // In this version of getFlux, the edgebox is passed in, and the flux array
  // is already defined.

  CH_TIME("FASPoissonOp::getFlux");
  CH_assert(a_dir >= 0 && a_dir <  SpaceDim);
  CH_assert(!a_data.box().isEmpty());
  // if this fails, the data box was too small (one cell wide, in fact)
  CH_assert(!a_edgebox.isEmpty());

  Real scale = m_beta * a_ref / m_dx[0];
  if( m_order == 2 )
    {
      FORT_GETFLUXLAP(CHF_FRA(a_flux),
		      CHF_CONST_FRA(a_data),
		      CHF_BOX(a_edgebox),
		      CHF_CONST_REAL(scale),
		      CHF_CONST_INT(a_dir));
    }
  else  if( m_order == 4 )
    {
      FORT_NEWGETFLUX4( CHF_FRA(a_flux),
			CHF_CONST_FRA(a_data),
			CHF_BOX(a_edgebox),
			CHF_CONST_REAL(scale),
			CHF_CONST_INT(a_dir));
    }
  else
    {
      MayDay::Abort("FASPoissonOp::getFlux bad order");
    }
}

// ---------------------------------------------------------
// smooth_sr
// ---------------------------------------------------------
void FASPoissonOp::smooth_sr( FArrayBox& phifab, const FArrayBox& rhsfab, Box phibox, int a_nsmooths )
{
  CH_TIME("private_sr");

  Real omega = m_smoothing_damping_factor, dx = m_dx[0];

  switch(m_smoother) 
    {
    case FAS_GSRB :
      for( int kk = 0 ; kk < a_nsmooths ; kk++ )
	{ 
	  for (int wpass = 0; wpass <= 1; wpass++ )
	    {
	      if (m_alpha == 0.0 && m_beta == -1.0 )
		{
		  FORT_GSRBLAPLACIAN_FAS( CHF_FRA(phifab),
					  CHF_CONST_FRA(rhsfab),
					  CHF_BOX(phibox),
					  CHF_CONST_REAL(dx),
					  CHF_CONST_INT(wpass),
					  CHF_CONST_REAL(omega)
					  );
		}
	      else
		{
		  FORT_GSRBHELMHOLTZ_FAS( CHF_FRA(phifab),
					  CHF_CONST_FRA(rhsfab),
					  CHF_BOX(phibox),
					  CHF_CONST_REAL(dx),
					  CHF_CONST_REAL(m_alpha),
					  CHF_CONST_REAL(m_beta),
					  CHF_CONST_INT(wpass),
					  CHF_CONST_REAL(omega));
		}
	    }
	}
      break;
    case FAS_RICH :       
      {
	Real diag = m_alpha - m_beta*2.*SpaceDim/(m_dx[0]*m_dx[0]);
	FArrayBox tmpFAB(phibox,1);
	// t <-- Ax
	for( int kk = 0 ; kk < a_nsmooths ; kk++ )
	  { 
	    
	    if( m_order == 2 )
	      {	  
		FORT_OPERATORLAP_FAS(CHF_FRA(tmpFAB),
				 CHF_CONST_FRA(phifab),
				 CHF_BOX(phibox),
				 CHF_CONST_REAL(dx),
				 CHF_CONST_REAL(m_alpha),
				 CHF_CONST_REAL(m_beta));
	      }
	    else  if( m_order == 4 )
	      {
		FORT_OPERATORLAP4(CHF_FRA(tmpFAB),
				  CHF_CONST_FRA(phifab),
				  CHF_BOX(phibox),
				  CHF_CONST_REAL(dx),
				  CHF_CONST_REAL(m_alpha),
				  CHF_CONST_REAL(m_beta));
	      }
	    else
	      {
		MayDay::Abort("FASPoissonOp::smooth_sr unnknown order");
	      }
	    
	    // x = x + omega D^-1 (b - Ax) -- D == alaph - 2*beta*Dim/h^2
	    tmpFAB.axby( tmpFAB, rhsfab, -1.0, 1.0 );
	    tmpFAB *= omega/diag;
	    phifab.plus( tmpFAB, phibox, 0, 0, 1 );
	  }
	break;
      }
    default:
      MayDay::Abort("FASPoissonOp::smooth_sr unknow smoother type");
    }
}


///
/**
   Poisson derived class for FAS operator factory
*/
// ---------------------------------------------------------
//  AMR Factory define function
void FASPoissonOpFactory::define(BCHolder                         a_bc,
				 Real                             a_alpha,
				 Real                             a_beta )
{
  CH_TIME("FASPoissonOpFactory::define");
  
  AMRFAS_LDFOpFactory::define( a_bc );  
  // 
  m_alpha = a_alpha;
  m_beta = a_beta;
}

// ---------------------------------------------------------
// FASPoissonOpFactory::AMRnewOp
RefCountedPtr<AMRFASOp<LevelData<FArrayBox> > > 
FASPoissonOpFactory::AMRNewOp( int a_ilev,
			       const DisjointBoxLayout &a_grid,
			       bool a_isSR )
{
  CH_TIME("FASPoissonOpFactory::AMRnewOp");

  RefCountedPtr<FASPoissonOp> Op_out = RefCountedPtr<FASPoissonOp>(new FASPoissonOp(m_order,a_grid));

  AMRFAS_LDFOpFactory::AMRNewOp( a_ilev, Op_out, a_isSR );

  Op_out->m_alpha = m_alpha;
  Op_out->m_beta  = m_beta;

  if( !a_isSR && a_ilev+1 < m_grids.size() )
    {
      Op_out->m_levfluxreg.define( m_grids[a_ilev+1],
				   a_grid,
				   m_domains[a_ilev+1],
				   m_refRatios[a_ilev],
				   1 );
    }

  return Op_out;
}

// ---------------------------------------------------------
// SR kernels -- too lazy to make a new op and factory so stuff in here
// ---------------------------------------------------------
//#define SR_PROL(f,c,b,r) FORT_PROLONGCONST(f,c,b,r)
//#define SR_PROL(f,c,b,r) FORT_PROLONGLINEAR(f,c,b,r)
//#define SR_PROL(f,c,b,r) FORT_PROLONGQUAD_V2(f,c,b,r)
//#define SR_PROL(f,c,b,r) FORT_PROLONGCUBIC(f,c,b,r)
//#define SR_PROL(f,c,b,r) FORT_PROLONGQUARTIC_R2(f,c,b)

Real sum( const LevelData<FArrayBox>& a_phi, int a_comp = 0 )
{
  DisjointBoxLayout dbl = a_phi.disjointBoxLayout();
  DataIterator dit = dbl.dataIterator();

  Real s = 0.0;
  for (dit.reset(); dit.ok(); ++dit)
    {
      const FArrayBox& phifab = a_phi[dit];
      Box phibox = dbl[dit];
      //phibox.grow(-1);
      FORT_SR_SUM( CHF_CONST_FRA(phifab),
		   CHF_BOX(phibox),
		   CHF_CONST_INT(a_comp),
		   CHF_REAL(s));
    }
  return s;
}

// ---------------------------------------------------------
// SR1_kernel: FMG prolongation + SR2:
//             smooth, restrict u, save R(u), form A(R(u)) + R(f-A(u))
//   'a_cover':  in: [u_H, ??]
//               output: [R(u), A(R(u)) + R(f-A(u))]
// ---------------------------------------------------------
void FASPoissonOp::SR1_kernel( LevelData<FArrayBox>& a_phi, // has buffer+ghost cells in LDF->ghost
			       const LevelData<FArrayBox>& a_rhs, 
			       LevelData<FArrayBox>& a_cover, // [u, r]
			       int a_nsmooths,
			       int a_refrat,
			       int a_nGhosts,
			       int a_nBuffer // num buffer cells on fine grid
			       )
{
  DisjointBoxLayout dbl = a_phi.disjointBoxLayout();
  DataIterator dit = dbl.dataIterator();

  for (dit.reset(); dit.ok(); ++dit)
    {
      FArrayBox& phifab = a_phi[dit];
      const FArrayBox& coverfab = a_cover[dit];
      Box phibox = phifab.box(); // prolongation to whole fine grid

      // H.O. FMG prolongation to whole fine grid, including ghosts
      // this only works because CH sees that coverfab has 2*nc vars and reads from 1st
      FORT_PROLONGQUAD_V2( CHF_FRA(phifab),
			   CHF_CONST_FRA(coverfab),
			   CHF_BOX(phibox),
			   CHF_CONST_INT(a_refrat) );
    }

  // just call SR2 for now, can fold later
  SR2_kernel( a_phi, a_rhs, a_cover, a_nsmooths, a_refrat, a_nGhosts, a_nBuffer );
}

// ---------------------------------------------------------
// SR2_kernel: smooth, restrict u, 
//             form A( R(u) ) + R( f - A(u) )
//   'a_cover':  in: [u_H, ??]
//               output: [R(u), A(R(u)) + R(f-A(u))] on reduced domain
//
//  Restriction
//    2 buffer, one ghost, update one buffer on coarse grid with buffer:
//              x oooo x ---- x ---- x ==== x ==== x ....
//                     |||||||||||||||||||||||||||||
//  oooo x ----------- x ----------- x =========== x ....
//
//   3 buffer, one ghost, update two buffer on coarse grid with ghost + buffer:
//       x oooo x ---- x ---- x ---- x ==== x ==== x ....
//       |||||||||||||||||||||||||||||||||||||||||||
//  ---- x ----------- x ----------- x =========== x ....
//
// ---------------------------------------------------------
void FASPoissonOp::SR2_kernel( LevelData<FArrayBox>& a_phi, 
			       const LevelData<FArrayBox>& a_rhs, 
			       LevelData<FArrayBox>& a_cover,
			       int a_nsmooths,
			       int a_refrat,
			       int a_nGhosts,
			       int a_nBuffer // num buffer cells on fine grid
			       ) 
{
  DisjointBoxLayout dbl = a_phi.disjointBoxLayout();
  Real dx = m_dx[0], crs_dx = a_refrat*dx;

  DataIterator dit = dbl.dataIterator();
  for (dit.reset(); dit.ok(); ++dit)
    {
      FArrayBox& phifab = a_phi[dit];
      const FArrayBox& rhsfab = a_rhs[dit];
      FArrayBox& coverfab = a_cover[dit];
      Box phibox = phifab.box(), cbox = phibox, cover_box = coverfab.box();
      phibox.grow( -a_nGhosts );          // only work in non-ghost region
      cbox.grow( -a_nGhosts-a_nBuffer );  // fine real part
      cbox.coarsen( a_refrat );           // real coarse grid
      cbox.grow( (a_nGhosts+a_nBuffer)/a_refrat ); // floor, coarse region that I can update
      CH_assert(cbox.longside()%2 == 0);

      // smooth
      smooth_sr( phifab, rhsfab, phibox, a_nsmooths );

      // restrict u to part of coverfab[0], outer buffer region is a bit stale now
      Box refbox( IntVect::Zero, (m_refToCoarser-1)*IntVect::Unit );
      FORT_AVERAGE_SR0( CHF_FRA(coverfab),       // out, will put in 1st comp(s)
			CHF_CONST_FRA(phifab),
			CHF_BOX(cbox),      // funny coarse box that includes only fine working region
			CHF_CONST_INT(a_refrat),
			CHF_BOX(refbox)
			);

      // special cover[1] = A(cover[0])
      FORT_OPERATORLAP_SR1( CHF_FRA(coverfab),
      			    CHF_BOX(cbox),            // we will do entire coarse grid, minus ghosts
      			    CHF_BOX(cover_box),       // box to zero out, for boundaries (debugging)
      			    CHF_CONST_REAL(crs_dx),
      			    CHF_CONST_REAL(m_alpha),
      			    CHF_CONST_REAL(m_beta));

      // special cover[1] += R( b - Au )
      FORT_RESTRICTRESADD_LAP_SR1( CHF_FRA(coverfab),
      				   CHF_CONST_FRA(phifab),
      				   CHF_CONST_FRA(rhsfab),
      				   CHF_CONST_REAL(m_alpha),
      				   CHF_CONST_REAL(m_beta),
      				   CHF_BOX(phibox),
      				   CHF_CONST_REAL(dx),
      				   CHF_CONST_INT(a_refrat)
      				   );
    }

}

//extern LevelData<FArrayBox> *pldf_exact_glob; // debug

// ---------------------------------------------------------
//  SR3_kernel: subtract init R(u) from coarse grid solution, prolongate & increment, smooth
// ---------------------------------------------------------
void FASPoissonOp::SR3_kernel( LevelData<FArrayBox>& a_phi,
			       const LevelData<FArrayBox>& a_rhs,
			       const LevelData<FArrayBox>& a_Ru,
			       LevelData<FArrayBox>& a_cover,
			       int a_nsmooths,
			       int a_refrat,
			       int a_nGhosts ) 
{
  const int nc = a_phi.nComp();
  DisjointBoxLayout dbl = a_phi.disjointBoxLayout();

  DataIterator dit = dbl.dataIterator();
  for (dit.reset(); dit.ok(); ++dit)
    {
      FArrayBox& phifab = a_phi[dit];
      const FArrayBox& rhsfab = a_rhs[dit];
      FArrayBox& coverfab = a_cover[dit];
      const FArrayBox& Rufab = a_Ru[dit];
      Box phibox = phifab.box(); 
      phibox.grow( -a_nGhosts ); // fine grid working region
      
      // subtract off initial solution to get increment
      coverfab.minus( Rufab, 0, 0, nc );

      // this only works because CH sees that coverfab has 2*nc vars and reads from 1st
      FORT_PROLONGLINEAR( CHF_FRA(phifab),
			  //FORT_PROLONGCONST( CHF_FRA(phifab),
			 CHF_CONST_FRA(coverfab),
			 CHF_BOX(phibox),
			 CHF_CONST_INT(a_refrat));
      //}
  



  // writeLevelname( &a_cover, "cover.hdf5" );
  // writeLevelname( &a_Ru, "Ru.hdf5" );
  // set error -- debug
  // if( pldf_exact_glob ) {
  //   LevelData<FArrayBox> err;
  //   err.define( a_phi.disjointBoxLayout(),1,a_phi.ghostVect());
  //   Real ldx = dx();
  //   for (DataIterator dit = a_phi.dataIterator(); dit.ok(); ++dit)
  //     {
  // 	const FArrayBox &exact = (*pldf_exact_glob)[dit];
  // 	FArrayBox &error = err[dit];
  // 	const FArrayBox &phi = a_phi[dit];
  // 	Box region = exact.box(); // not ghosted
  // 	for (BoxIterator bit(region); bit.ok(); ++bit)
  // 	  {
  // 	    const RealVect offset = bit(); //-domain.smallEnd();
  // 	    const RealVect x = ldx*(0.5+offset);
  // 	    error(bit()) = exact(bit()) - phi(bit());
  // 	  }
  // 	writeFABname( &error, "error2FAB.hdf5" );
  //     }
  //   writeLevelname( &err, "error2.hdf5" );
  // }



  // for (dit.reset(); dit.ok(); ++dit)
  //   {
      // FArrayBox& phifab = a_phi[dit];
      // const FArrayBox& rhsfab = a_rhs[dit];
      // FArrayBox& coverfab = a_cover[dit];
      // const FArrayBox& Rufab = a_Ru[dit];
      // Box phibox = phifab.box(); 
      // phibox.grow( -a_nGhosts ); // fine grid working region



      // smooth
      smooth_sr( phifab, rhsfab, phibox, a_nsmooths );
    }


  // set error -- debug
  // if( pldf_exact_glob ) {
  //   LevelData<FArrayBox> err;
  //   err.define( a_phi.disjointBoxLayout(),1,a_phi.ghostVect());
  //   Real ldx = dx();
  //   for (DataIterator dit = a_phi.dataIterator(); dit.ok(); ++dit)
  //     {
  // 	const FArrayBox &exact = (*pldf_exact_glob)[dit];
  // 	FArrayBox &error = err[dit];
  // 	const FArrayBox &phi = a_phi[dit];
  // 	Box region = exact.box(); // not ghosted
  // 	for (BoxIterator bit(region); bit.ok(); ++bit)
  // 	  {
  // 	    const RealVect offset = bit(); //-domain.smallEnd();
  // 	    const RealVect x = ldx*(0.5+offset);
  // 	    error(bit()) = exact(bit()) - phi(bit());
  // 	  }
  // 	writeFABname( &error, "error3FAB.hdf5" );
  //     }
  //   writeLevelname( &err, "error3.hdf5" );
  // }


}

// ---------------------------------------------------------
// SR4_kernel: SR3 + functional
// ---------------------------------------------------------
void FASPoissonOp::SR4_kernel( LevelData<FArrayBox>& a_phi,
			       const LevelData<FArrayBox>& a_rhs,
			       const LevelData<FArrayBox>& a_ru, 
			       LevelData<FArrayBox>& a_cover,
			       int a_nsmooths,
			       int a_refrat,
			       int a_nGhosts,
			       int a_nBuffer, 
			       SRUserFunctional &a_func ) 
{

  // SR3
  SR3_kernel(a_phi,a_rhs,a_ru,a_cover,a_nsmooths,a_refrat,a_nGhosts);

  DisjointBoxLayout dbl = a_phi.disjointBoxLayout();
  Real dx = m_dx[0];
  DataIterator dit = dbl.dataIterator();
  for (dit.reset(); dit.ok(); ++dit)
    {
      const FArrayBox& phifab = a_phi[dit];
      Box phibox = phifab.box();
      phibox.grow( -(a_nGhosts+a_nBuffer) ); 
      
      a_func(phifab,phibox,dx);
    }
}

#include "NamespaceFooter.H"
