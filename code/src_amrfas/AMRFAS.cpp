#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "FORT_PROTO.H"
#include "AMRFAS.H"

#include "AMRFASOpF_F.H"
#include "AverageF_F.H"
#include "ParmParse.H"

#include "NamespaceHeader.H"


///
/**
    LevelData<FArrayBox> derived class for FAS operator
*/

void 
AMRFAS_LDFOp::apply( LevelData<FArrayBox>& a_p,
		     const LevelData<FArrayBox>& a_phi,
		     const LevelData<FArrayBox> *a_phiCoarse,
		     bool a_doExchange )
{
  CH_TIME("AMRFAS_LDFOp::apply");
  CH_assert(a_phi.isDefined());

  LevelData<FArrayBox>& phi = (LevelData<FArrayBox>&)a_phi; // we update ghosts so descard const

  // fill in intersection of ghost cells and a_phi's boxes
  if ( a_phiCoarse && a_phiCoarse->isDefined() && a_doExchange )
    {
      CH_TIME("AMRFAS_LDFOp::apply:coarseFineInterp");
      CFInterp( phi, *a_phiCoarse );
    }
  else if( a_phiCoarse && !a_phiCoarse->isDefined() ) 
    MayDay::Abort("AMRFAS_LDFOp::apply !a_phiCoarse->isDefined()");

  if( a_doExchange )
    {
      phi.exchange( m_exchangeCopier );
    }

  // this applies BCs 
  applyLevel( a_p, a_phi );
}

// ---------------------------------------------------------
// AMRProlong
// ---------------------------------------------------------
void AMRFAS_LDFOp::AMRProlong( LevelData<FArrayBox>&       a_fu,
			       const LevelData<FArrayBox>& a_cu,
			       LevelData<FArrayBox>&       a_temp, // m_solC
			       RefCountedPtr<AMRFASOp<LevelData<FArrayBox> > > a_crsOp,
			       FAS_PROLONG_type a_type )
{ 
  CH_TIME("AMRFAS_LDFOp::AMRProlong");
  
  DisjointBoxLayout dbl = a_fu.disjointBoxLayout();
  DisjointBoxLayout cdbl = a_temp.disjointBoxLayout();
  
  a_cu.copyTo( a_temp.interval(), a_temp, a_temp.interval(), m_HOCopier );
  
  for ( DataIterator dit = a_fu.dataIterator(); dit.ok(); ++dit )
    {
      FArrayBox& phi =  a_fu[dit];
      FArrayBox& coarse = a_temp[dit];
      Box region = dbl[dit];
      const IntVect& iv = region.smallEnd();
      IntVect civ = coarsen(iv, m_refToCoarser);
 
      switch(a_type)
	{
	case PROL_CONST_1:
	  FORT_PROLONGCONST(CHF_FRA_SHIFT(phi, iv),
			    CHF_CONST_FRA_SHIFT(coarse, civ),
			    CHF_BOX_SHIFT(region, iv),
			    CHF_CONST_INT(m_refToCoarser)
			    );
	  break;
	default:

	  //a_crsOp->m_bc( coarse, cdbl[dit], a_crsOp->m_domain, a_crsOp->m_dx[0], true );
	  //a_temp.exchange( a_temp.interval(), m_HOCornerCopier ); // needed for AMR?
	  
	  switch(a_type)
	    {
	    case PROL_QUAD_3:
	      FORT_PROLONGQUAD_V2(CHF_FRA(phi),
				  CHF_CONST_FRA(coarse),
				  CHF_BOX(region),
				  CHF_CONST_INT(m_refToCoarser));
	      break;
	    case PROL_CUBIC_4:
	      FORT_PROLONGCUBIC(CHF_FRA(phi),
				CHF_CONST_FRA(coarse),
				CHF_BOX(region),
				CHF_CONST_INT(m_refToCoarser));
	      break;
	    case PROL_QUART_5:
	      
	      if( m_refToCoarser == 2 )
		{
		  FORT_PROLONGQUARTIC_R2(CHF_FRA(phi),
					 CHF_CONST_FRA(coarse),
					 CHF_BOX(region));
		}
	      else if( m_refToCoarser == 4 )
		{		  
		  FORT_PROLONGQUARTIC_R4(CHF_FRA(phi),
					 CHF_CONST_FRA(coarse),
					 CHF_BOX(region));
		}
	      else
		MayDay::Abort("AMRFAS_LDFOp::AMRProlong unsupported refinement ratio (2)");

	      break;

	    case PROL_LINEAR_2:
	      FORT_PROLONGLINEAR(CHF_FRA(phi),
				    CHF_CONST_FRA(coarse),
				    CHF_BOX(region),
				    CHF_CONST_INT(m_refToCoarser));
	      break;
	    default:
	      MayDay::Abort("AMRFAS_LDFOp::AMRProlong type undefined");
	    
	    }
	  break;
	}
    }  
}

// ---------------------------------------------------------
// AMRRestrict
// ---------------------------------------------------------
// ---------------------------------------------------------
void AMRFAS_LDFOp::AMRRestrict( LevelData<FArrayBox>&       a_cf,      // output
				const LevelData<FArrayBox>& a_ff,      // input
				FAS_RESTRICT_type a_type_dummy ) const // only one type implemented now
{
  CH_assert(a_cf.nComp() == a_ff.nComp());
  CH_TIME("AMRFAS_LDFOp::AMRRestrict");
  
  DisjointBoxLayout dblCoar = a_cf.disjointBoxLayout();
  
  DataIterator dit = a_ff.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      FArrayBox& coarse = a_cf[dit];
      const FArrayBox& fine = a_ff[dit];
      const Box& b = dblCoar[dit];
      Box refbox( IntVect::Zero, (m_refToCoarser-1)*IntVect::Unit);

      FORT_AVERAGE( CHF_FRA(coarse),
                    CHF_CONST_FRA(fine),
                    CHF_BOX(b),
                    CHF_CONST_INT(m_refToCoarser),
                    CHF_BOX(refbox)
                    );

    }
}

// ---------------------------------------------------------
//  basic vector ops 
/// ---------------------------------------------------------

// ---------------------------------------------------------
void AMRFAS_LDFOp::create( LevelData<FArrayBox>&       a_lhs,
			  const LevelData<FArrayBox>& a_rhs )
{
  CH_TIME("AMRFAS_LDFOp::create");

  m_levelOps.create( a_lhs, a_rhs );
}

void AMRFAS_LDFOp::assign( LevelData<FArrayBox>&       a_lhs,
                          const LevelData<FArrayBox>& a_rhs)
{
  CH_TIME("AMRFAS_LDFOp::assign");

  m_levelOps.assign(a_lhs, a_rhs);
}

// ---------------------------------------------------------
void AMRFAS_LDFOp::assignLocal(LevelData<FArrayBox>&        a_lhs,
                               const  LevelData<FArrayBox>& a_rhs)
{
  CH_TIME("AMRFAS_LDFOp::assignLocal");

  for (DataIterator dit= a_lhs.dataIterator(); dit.ok(); ++dit)
    {
      a_lhs[dit].copy(a_rhs[dit]);
    }
}

// ---------------------------------------------------------
void AMRFAS_LDFOp::assignCopier( LevelData<FArrayBox>&       a_lhs,
                                 const LevelData<FArrayBox>& a_rhs,
                                 const Copier&               a_copier)
{
  CH_TIME("AMRFAS_LDFOp::assignCopier"); // static method really
  a_rhs.copyTo(a_rhs.interval(), a_lhs, a_lhs.interval(), a_copier);
}

// ---------------------------------------------------------
void AMRFAS_LDFOp::zeroCovered( LevelData<FArrayBox>& a_lhs,
			       LevelData<FArrayBox>& a_rhs,
			       const Copier&         a_copier)
{
  CH_TIME("AMRFAS_LDFOp::zeroCovered");

  m_levelOps.copyToZero( a_lhs, a_copier );
}

// ---------------------------------------------------------
Real AMRFAS_LDFOp::dotProduct(const LevelData<FArrayBox>& a_1,
                              const LevelData<FArrayBox>& a_2)
{
  CH_TIME("AMRFAS_LDFOp::dotProduct");

  return m_levelOps.dotProduct(a_1, a_2);
}

// ---------------------------------------------------------
void AMRFAS_LDFOp::mDotProduct(const LevelData<FArrayBox>& a_1,
                               const int a_sz,
                               const LevelData<FArrayBox> a_2[],
                               Real a_mdots[])
{
  CH_TIME("AMRFAS_LDFOp::mDotProduct");

  m_levelOps.mDotProduct(a_1, a_sz, a_2, a_mdots);
}

// ---------------------------------------------------------
void AMRFAS_LDFOp::incr( LevelData<FArrayBox>&       a_lhs,
                         const LevelData<FArrayBox>& a_x,
                         Real                        a_scale)
{
  CH_TIME("AMRFAS_LDFOp::incr");

  m_levelOps.incr(a_lhs, a_x, a_scale);
}

// ---------------------------------------------------------
void AMRFAS_LDFOp::axby( LevelData<FArrayBox>&       a_lhs,
			const LevelData<FArrayBox>& a_x,
			const LevelData<FArrayBox>& a_y,
			Real                        a_a,
			Real                        a_b)
{
  CH_TIME("AMRFAS_LDFOp::axby");
  
  m_levelOps.axby(a_lhs, a_x, a_y, a_a, a_b);
}

// ---------------------------------------------------------
void AMRFAS_LDFOp::mult( LevelData<FArrayBox>&       a_lhs,
			 const LevelData<FArrayBox>& a_x )
{
  CH_TIME("AMRFAS_LDFOp::mult");
  
  m_levelOps.mult(a_lhs, a_x);
}

// ---------------------------------------------------------
void AMRFAS_LDFOp::scale(LevelData<FArrayBox>& a_lhs,
                         const Real&           a_scale)
{
  CH_TIME("AMRFAS_LDFOp::scale");

  m_levelOps.scale(a_lhs, a_scale);
}

// ---------------------------------------------------------
Real AMRFAS_LDFOp::norm(const LevelData<FArrayBox>& a_x,
			int                         a_ord) const
{
  CH_TIME("AMRFAS_LDFOp::norm");

  return CH_XD::norm(a_x, a_x.interval(), a_ord);
}

// ---------------------------------------------------------
Real AMRFAS_LDFOp::norm(const LevelData<FArrayBox>& a_x,
			int                         a_ord,
			int a_comp ) const
{
  CH_TIME("AMRFAS_LDFOp::norm");
  Interval interv(a_comp,a_comp);
  return CH_XD::norm(a_x, interv, a_ord);
}

// ---------------------------------------------------------
Real AMRFAS_LDFOp::localMaxNorm(const LevelData<FArrayBox>& a_x) const
{
  CH_TIME("AMRFAS_LDFOp::localMaxNorm");

  Real localMax = 0;
  int nComp=a_x.nComp();
  for (DataIterator dit=a_x.dataIterator(); dit.ok(); ++dit)
    {
      localMax = Max(localMax, a_x[dit].norm(a_x.box(dit()), 0, 0, nComp));
    }
  return localMax;
}

// ---------------------------------------------------------
void AMRFAS_LDFOp::setToZero(LevelData<FArrayBox>& a_lhs) 
{
  CH_TIME("AMRFAS_LDFOp::setToZero");

  m_levelOps.setToZero( a_lhs );
}
 
// ---------------------------------------------------------
void AMRFAS_LDFOp::write(const LevelData<FArrayBox>* a_data,
			const char*             a_filename)
{
#ifdef CH_USE_HDF5
  writeLevelname( a_data, a_filename );
#else
  MayDay::Warning("AMRFAS_LDFOp::write unimplemented since CH_USE_HDF5 undefined");
#endif
}

// ---------------------------------------------------------
//  default R and P 
void AMRFAS_LDFOp::AMRProlong( LevelData<FArrayBox>&       a_fineU,
			       const LevelData<FArrayBox>& a_CrsU,
			       LevelData<FArrayBox>&       a_temp,
			       RefCountedPtr<AMRFASOp<LevelData<FArrayBox> > > a_crsOp )
{
  switch(m_ProlOrderP)
    {
    case 0:
      AMRProlong( a_fineU, a_CrsU, a_temp, a_crsOp, PROL_CONST_1 );
      break;
    case 1:
      AMRProlong( a_fineU, a_CrsU, a_temp, a_crsOp, PROL_LINEAR_2 );
      break;
    case 2:
      AMRProlong( a_fineU, a_CrsU, a_temp, a_crsOp, PROL_QUAD_3 );
      break;
    case 3:
      AMRProlong( a_fineU, a_CrsU, a_temp, a_crsOp, PROL_CUBIC_4 );
      break;
    case 4:
      AMRProlong( a_fineU, a_CrsU, a_temp, a_crsOp, PROL_QUART_5 );
      break;
    default:
      MayDay::Abort("AMRFAS_LDFOp::AMRProlong unknow poly degree");
    }
}

void AMRFAS_LDFOp::AMRFMGProlong( LevelData<FArrayBox>&       a_fineU,
				  const LevelData<FArrayBox>& a_CrsU,
				  LevelData<FArrayBox>&       a_temp,
				  RefCountedPtr<AMRFASOp<LevelData<FArrayBox> > > a_crsOp )
{
  switch(m_FMGProlOrderP)
    {
    case 0:
      AMRProlong( a_fineU, a_CrsU, a_temp, a_crsOp, PROL_CONST_1 );
      break;
    case 1:
      AMRProlong( a_fineU, a_CrsU, a_temp, a_crsOp, PROL_LINEAR_2 );
      break;
    case 2:
      AMRProlong( a_fineU, a_CrsU, a_temp, a_crsOp, PROL_QUAD_3 );
      break;
    case 3:
      AMRProlong( a_fineU, a_CrsU, a_temp, a_crsOp, PROL_CUBIC_4 );
      break;
    case 4:
      AMRProlong( a_fineU, a_CrsU, a_temp, a_crsOp, PROL_QUART_5 );
      break;
    default:
      MayDay::Abort("AMRFAS_LDFOp::AMRFMGProlong unknow  poly degree");
    }
}

void AMRFAS_LDFOp::AMRRestrict( LevelData<FArrayBox>& a_CrsU,         // output
				const LevelData<FArrayBox>& a_fineU,  // input (this)
				LevelData<FArrayBox>& a_crsCover,
				const Copier &a_copier ) const
{
  AMRRestrict( a_crsCover, a_fineU, REST_CONST_1 ); // default
  
  // Overwrite R(u_f) on the valid region of the next coarser level a_phi[a_ilev-1]
  //assignCopier( a_CrsU, a_crsCover, a_copier );
  a_crsCover.copyTo(a_crsCover.interval(), a_CrsU, a_CrsU.interval(), a_copier);
}

///
/**
   LevelData<FArrayBox> > derived class for FAS operator factory
*/

// ---------------------------------------------------------
// AMRFAS_LDFOpFactory::AMRnewOp
void AMRFAS_LDFOpFactory::AMRNewOp( int a_ilev,
				    RefCountedPtr<AMRFAS_LDFOp> Op_out,
				    bool a_isSR )
{
  CH_TIME("AMRFAS_LDFOpFactory::AMRnewOp");

  CFRegion dummy_cf;
  const int finestLevel = m_domains.size() - 1;
  DisjointBoxLayout *this_dbl, *crs_dbl = 0; // flag for not the coarsest grid
  ProblemDomain crs_pdom;
  int refToCrs;

  RealVect dx = m_dx[a_ilev];

  if (a_ilev == 0) // coarsest AMR level
    {
      if ( a_ilev == finestLevel )
        {
          // no finer level -- one level solve
          Op_out->define( dx, m_domains[0], m_bc,
			  m_exchangeCopiers[0], 
			  m_cfregion[0] );
        }
      else 
        {
          // finer level exists but no coarser -- normal coarse grid 
          int dummyRat = 1;  // argument so compiler can find right function
          int refToFiner = m_refRatios[0]; // actual refinement ratio
          Op_out->define( m_grids[0], m_grids[1], dx,
			  dummyRat, refToFiner,
			  m_domains[0], m_bc,
			  m_exchangeCopiers[0],
			  m_cfregion[0],  
			  m_ncomp );	  
        }
    }
  else // not coarsest level
    {
      crs_pdom = m_domains[a_ilev];
      refToCrs = m_refRatios[a_ilev-1];
      this_dbl = &m_grids[a_ilev];
      crs_dbl = &m_grids[a_ilev-1];

      Copier copier = a_isSR ? Copier() :  m_exchangeCopiers[a_ilev];
      CFRegion cfr = a_isSR ? CFRegion() :  m_cfregion[a_ilev];
      if ( a_ilev == finestLevel )
	{
	  // finest level
	  Op_out->define( m_grids[a_ilev], 
			  m_grids[a_ilev-1], 
			  dx,
			  m_refRatios[a_ilev-1],
			  m_domains[a_ilev], 
			  m_bc,
			  copier,
			  cfr,
			  a_isSR
			  );
	}
      else 
	{
	  // intermediate user AMR level, full define
	  Op_out->define( m_grids[a_ilev], 
			  m_grids[a_ilev+1], 
			  m_grids[a_ilev-1], 
			  dx,
			  m_refRatios[a_ilev-1], 
			  m_refRatios[a_ilev],
			  m_domains[a_ilev], 
			  m_bc,
			  copier,
			  cfr,
			  m_ncomp,
			  a_isSR
			  );
	}
    }

  // pass down to op, coarse grid not damped???
  if( !crs_dbl && false ) Op_out->m_smoothing_damping_factor = 1.0;
  else Op_out->m_smoothing_damping_factor = m_smoothing_damping_factor;

  Op_out->m_FMGProlOrderP = m_FMGProlOrderP;
  Op_out->m_ProlOrderP = m_ProlOrderP;
}

#include "NamespaceFooter.H"
