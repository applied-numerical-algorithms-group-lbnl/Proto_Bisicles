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

#ifndef _AMRFASI_H_
#define _AMRFASI_H_

#include "FORT_PROTO.H"
#include "AMRFAS.H"
#include "ParmParse.H"
#include "BRMeshRefine.H"

#include "NamespaceHeader.H"

///
/**
   Base class for FAS operator
*/

static int s_fas_verbosity = 0;

// ---------------------------------------------------------
// residual
// ---------------------------------------------------------
template <class T>
void 
AMRFASOp<T>::residual( T& a_resid,
		       const T& a_phi,
		       const T& a_rhs,
		       const T *a_phiCoarse,
		       T *a_phiFine, // this has its ghost C-F points set with this level 
		       AMRFASOp<T> *a_finerOp )
{
  CH_TIME("AMRFASOp::AMRResidual");

  apply( a_resid, a_phi, a_phiCoarse ); // sets this C-F ghost point

  // fine grid could be SR and not defined refluxer yet (need to do for AMR)
  if ( a_phiFine && a_phiFine->isDefined() )
    {
      CH_assert(a_finerOp);
      reflux( *a_phiFine, a_phi, a_resid, a_finerOp ); // this sets fine ghosts with C-F interp with ops C-F interp
    }
  else if ( a_phiFine )
    {
      MayDay::Abort("AMRFASOp<T>::residual !a_phiFine->isDefined()");
    }

  // b - L u
  axby( a_resid, a_resid, a_rhs, -1.0, 1.0 );
}

// ---------------------------------------------------------
// smooth
// ---------------------------------------------------------
template <class T>
void
AMRFASOp<T>::smooth( RefCountedPtr<T> a_phi,
		     const  RefCountedPtr<T> a_rhs,
		     const T* a_phiCoarse,
		     int a_num_it,
		     bool a_print )
{
  CH_TIME("AMRFASOp::smooth");

  // compute state -- linear smoother (perhaps)
  bool did_exchange = false; 

  if ( a_phiCoarse && a_phiCoarse->isDefined() )
    {
      CH_TIME("AMRFASOp::smooth::coarseFineInterp");
      // fill in intersection of ghost cells and a_phi's boxes 
      CFInterp( *a_phi, *a_phiCoarse );
    }

  Real res, res0 = 0.; int i, interval = 10;
  for (i = 0; i < a_num_it; i++)
    {
      if( !did_exchange ){
	CH_TIME("smooth:exchange");
	a_phi->exchange( m_exchangeCopier );
      }
      did_exchange = false; // it will be dirty next

      if( a_print && i==0 )
	{
	  T tmp;
	  create( tmp, *a_rhs );
	  residual( tmp, *a_phi, *a_rhs, a_phiCoarse );
	  res0 = norm(tmp,0);
	  if( res0 == 0. )
	    {
	      if( s_fas_verbosity > 0 ) pout()<<"AMRFASOp<T>::smooth: inital residual is zero" << endl;
	      break;
	    }
	}

      smooth_private( a_phi, a_rhs, a_phiCoarse );

      if( a_print && i%interval == 0 )
	{
	  T tmp;
	  create( tmp, *a_rhs );
	  residual( tmp, *a_phi, *a_rhs, a_phiCoarse );
	  res = norm(tmp,0);
	  if( s_fas_verbosity > 6 ) pout()<<"AMRFASOp<T>::smooth: " << i << ") |res|=" << res << endl;
	  if( res/res0 < 1.e-10 ) break;
	}
    }

  // seeem to get stagnation at res ~ 1.e-10 in BISICLES
  if(a_print && i==a_num_it && (s_fas_verbosity>3 || (s_fas_verbosity>0 && res/res0 > 1.e-5 && res > 1.e-8 )))
    {
      pout()<<"AMRFASOp<T>::smooth WARNING coarse solver failed to converge: rel. res = " << res/res0 << ", " << i << " iterations" << ", res = " << res << endl;
      if(res > 1.e10)
	{
	  T tmp;
	  create( tmp, *a_phi );
	  residual( tmp, *a_phi, *a_rhs, a_phiCoarse );
	  //write(&tmp, "coarse_res.hdf5" );
	  //write(a_phi, "coarse_phi.hdf5" );
	}
    }
}

// ---------------------------------------------------------
/** full define function for AMRFASOp with both coarser and finer levels */
template <class T>
void 
AMRFASOp<T>::define(const DisjointBoxLayout& a_grids,
		    const DisjointBoxLayout& a_gridsFiner,
		    const DisjointBoxLayout& a_gridsCoarser,
		    const RealVect&          a_dxLevel,
		    int                      a_refRatioCrs,
		    int                      a_refRatioFiner,
		    const ProblemDomain&     a_domain,
		    BCHolder                 a_bc,
		    const Copier&            a_exchange,
		    const CFRegion&          a_cfregion,
		    int ncomp,
		    bool a_isSR  )
  
{
  CH_TIME("AMRFASOp::define1");

  // define a fine level
  define( a_grids, a_gridsCoarser, a_dxLevel, a_refRatioCrs, a_domain, a_bc,
	  a_exchange, a_cfregion, a_isSR );

  m_refToFiner = a_refRatioFiner;
}

// coarsest level ------------------------------------------------------
template <class T>
void 
AMRFASOp<T>::define(const DisjointBoxLayout& a_grids,
		    const DisjointBoxLayout& a_gridsFiner,
		    const RealVect&          a_dxLevel,
		    int                      a_refRatioCrs, // dummy arg
		    int                      a_refRatioFiner,
		    const ProblemDomain&     a_domain,
		    BCHolder                 a_bc,
		    const Copier&            a_exchange,
		    const CFRegion&          a_cfregion,
		    int ncomp )
{
  CH_TIME("AMRFASOp::define2");

  CH_assert(a_refRatioCrs == 1); // dummy arg

  // general level 
  define( a_dxLevel, a_domain, a_bc, a_exchange, a_cfregion );
  m_refToFiner = a_refRatioFiner;
  m_refToCoarser = 1;
}

// a fine level --------------------------------------------
template <class T>
void 
AMRFASOp<T>::define( const DisjointBoxLayout& a_grids,
		     const DisjointBoxLayout& a_coarse,
		     const RealVect&          a_dxLevel,
		     int                      a_refRatioCrs,
		     const ProblemDomain&     a_domain,
		     BCHolder                 a_bc,
		     const Copier&            a_exchange,
		     const CFRegion&          a_cfregion,
		     bool a_isSR  )
{
  CH_assert(m_order%2 == 0 && m_order/2 > 0);
  CH_TIME("AMRFASOp::define3");

  define( a_dxLevel, a_domain, a_bc, a_exchange, a_cfregion );
  m_refToCoarser = a_refRatioCrs;

  m_dxCrse = a_refRatioCrs*a_dxLevel;
  m_refToFiner = 1; 
}

// ---------------------------------------------------------
// base define
template <class T>
void 
AMRFASOp<T>::define( const RealVect&          a_dx,
		     const ProblemDomain&     a_domain,
		     BCHolder                 a_bc,
		     const Copier&            a_exchange,
		     const CFRegion&          a_cfregion )
{
  CH_TIME("AMRFASOp::define4");

  m_bc     = a_bc;
  m_domain = a_domain;
  m_dx     = a_dx;
  
  m_dxCrse = 2.*a_dx; // default ref of 2 -- this should always be overwriten

  m_exchangeCopier = a_exchange; // these can be dummies for SR grids
  m_cfregion = a_cfregion;
}

// ---------------------------------------------------------
// AMRFASOp::getOrder
template <class T>
int AMRFASOp<T>::getOrder() const
{ 
  return m_order; 
}

// ---------------------------------------------------------
// smooth_private
// ---------------------------------------------------------
template <class T>
void AMRFASOp<T>::smooth_private( RefCountedPtr<T> a_phi,
				  const  RefCountedPtr<T> a_rhs,
				  bool do_exchange 
				  )
{
  CH_TIME("AMRFASOp::smooth_private");

  switch(m_smoother) 
    {
    case FAS_GSRB :
      levelGSRB( a_phi, a_rhs );
      break;
    case FAS_RICH :
      levelRich( a_phi, a_rhs );
      break;
    default:
      MayDay::Abort("AMRFASOp::smooth_private unknow smoother type");
    }

}

///
/**
   base operator factory
*/

// ---------------------------------------------------------
//  AMR Factory define 
//    - set BC method, called in app
//
template <class T>
void AMRFASOpFactory<T>::define( BCHolder a_bc )
{
  m_bc = a_bc;
}
// ---------------------------------------------------------
//  AMR Factory define 
//    - make copiers, CFregions, domains, dx for all levels in 'a_grids'
//    - called by FASSolver with augmented grids
//
template <class T>
void AMRFASOpFactory<T>::define( const ProblemDomain& a_coarseDomain, 
				 const RealVect&      a_crsDx,
				 const Vector<DisjointBoxLayout>& a_grids,
				 const Vector<int>&   a_refRatios,
				 int a_nSRGrids
				 )
{
  CH_assert(m_order%2 == 0 && m_order/2 > 0);
  CH_TIME("AMRFASOpFactory::define");
  const int nGhost = m_order/2; // hack for number of ghosts!

  m_dx.resize( a_grids.size() ); // current coarse grid, before 
  m_domains.resize( a_grids.size() );

  m_dx[0] = a_crsDx;
  m_domains[0] = a_coarseDomain;

  m_cfregion.resize( a_grids.size() - a_nSRGrids );
  m_exchangeCopiers.resize( a_grids.size() - a_nSRGrids );

  m_grids = a_grids;
  m_refRatios = a_refRatios;

  m_exchangeCopiers[0].exchangeDefine( a_grids[0], nGhost*IntVect::Unit );
  m_exchangeCopiers[0].trimEdges(a_grids[0], nGhost*IntVect::Unit);
  m_cfregion[0].define( a_grids[0], a_coarseDomain );

  for (int i = 1; i < a_grids.size(); i++)
    {
      m_dx[i] = m_dx[i-1] / m_refRatios[i-1];

      m_domains[i] = m_domains[i-1];
      m_domains[i].refine(m_refRatios[i-1]);

      // skip construction of coperiers and CF region for SR grids
      if ( a_grids[i].isClosed() && i < a_grids.size()-a_nSRGrids )
        {
          m_exchangeCopiers[i].exchangeDefine(a_grids[i], nGhost*IntVect::Unit);
          m_exchangeCopiers[i].trimEdges(a_grids[i], nGhost*IntVect::Unit);
          m_cfregion[i].define(a_grids[i], m_domains[i]);
        }
      else if(i < a_grids.size()-a_nSRGrids)
	MayDay::Abort("!a_grids[i].isClosed() ???");
    }
}


// --------------------------------------------------------- 
template <class T>
int AMRFASOpFactory<T>::refToFiner( const ProblemDomain& a_domain ) const
{
  int retval = -1;
  bool found = false;

  for (int ilev = 0; ilev < m_domains.size(); ilev++)
    {
      if (m_domains[ilev].domainBox() == a_domain.domainBox())
        {
          retval = m_refRatios[ilev];
          found = true;
        }
    }

  if (!found)
    {
      MayDay::Abort("Domain not found in AMR hierarchy");
    }

  return retval;
}


//*******************************************************
// AMRFAS Implementation
//*******************************************************

////////////////////////////////////////////////////////////////////////
// AMRFAS<T>::AMRFAS
////////////////////////////////////////////////////////////////////////
template <class T>
AMRFAS<T>::AMRFAS() :
  m_numInternalLevels(0),
  m_type(FAS_FULL),
  m_smoother(FAS_GSRB),
  m_smoothing_damping_factor(1.0),
  m_FMGProlOrderP(-1),
  m_ProlOrderP(-1),
  m_avoid_norms(false),
  m_num_vcycles(1), // number of V cycles in FMG
  m_privateCoarseningRate(2),
  m_max_iter(20),
  m_min_iter(1),
  m_verbosity(0),
  m_rtol(1E-6),
  m_stagnate(1E-2),
  m_atol(1E-30),
  m_pre(2),
  m_post(2),
  m_num_cycles(1), // 1 for V-cycles, 2 for W-cycles
  m_coarse_its(50),
  m_init(false),
  m_plot_residual(false),
  m_fixedBoxSize(0)
{
}

////////////////////////////////////////////////////////////////////////
// AMRFAS<T>::~AMRFAS
////////////////////////////////////////////////////////////////////////
template <class T>
AMRFAS<T>::~AMRFAS()
{
  CH_TIME("~AMRFAS");
  clear();
}

////////////////////////////////////////////////////////////////////////
// AMRFAS<T>::clear
////////////////////////////////////////////////////////////////////////
template <class T>
void AMRFAS<T>::clear()
{
  m_init = false; // nothing to clean up from init() ???
}

////////////////////////////////////////////////////////////////////////
// AMRFAS<T>::setParameters
////////////////////////////////////////////////////////////////////////
template <class T>
void AMRFAS<T>::setParameters( const char *name )
{
  ParmParse ppSolver(name);

  // multigrid solver parameters2
  int numSmooth = 2, numMG = 1;
  ppSolver.query("num_smooth", numSmooth);
  m_pre = m_post = numSmooth;
  ppSolver.query("num_cycles", numMG);
  m_num_cycles = numMG; // 1 for V-cycles, 2 for W-cycles
  ppSolver.query("max_iterations", m_max_iter);
  ppSolver.query("min_iterations", m_min_iter);
  ppSolver.query("tolerance", m_rtol);
  int sm = (int)FAS_GSRB;
  ppSolver.query( "smoother", sm ); // FAS_GSRB=0,FAS_RICH=1
  setSmootherType( (FAS_SMOOTHER_type)sm );

  ppSolver.query("stagnate", m_stagnate);
  m_atol = 1.0e-30;
  ppSolver.query( "absolute_tolerance", m_atol ); 
  ppSolver.query("verbosity", m_verbosity);
  s_fas_verbosity = m_verbosity;
  // optional parameters
  ppSolver.query("num_pre", m_pre);
  ppSolver.query("num_post", m_post);
  // FAS stuff
  int type = (int)FAS_VCYCLE;
  ppSolver.query("cycle_type", type); // FASMG_type {FAS_FULL=0,FAS_VCYCLE=1,FAS_FCYCLE=2,FAS_SMOOTH=3};
  setCycleType( (FASMG_type)type );
  bool avoid_norms = false;
  ppSolver.query("avoid_norms", avoid_norms);
  setAvoidNorms(avoid_norms);
  int numv=1;
  ppSolver.query("num_v_cycles", numv);
  setNumVcycles(numv);
  ppSolver.query("plot_residual",m_plot_residual);
  ppSolver.query("smoothing_damping_factor",m_smoothing_damping_factor);

  ppSolver.query("prol_poly_order",m_ProlOrderP);
  ppSolver.query("fmg_prol_poly_order",m_FMGProlOrderP);
}

// hack to adjust this
template <class T>
void AMRFAS<T>::setInternalCoarseningRate( int r )
{
  m_privateCoarseningRate = r; // 2 or 4
}
////////////////////////////////////////////////////////////////////////
// AMRFAS<T>::setFixedBoxSize
////////////////////////////////////////////////////////////////////////
template <class T>
void AMRFAS<T>::setFixedBoxSize( const DisjointBoxLayout& a_grids )
{
  if(m_fixedBoxSize == 0) 
    {
      DataIterator dit = a_grids.dataIterator();
      for (dit.begin(); dit.ok(); ++dit)
	{
	  Box box = a_grids[dit];
	  int sz = box.longside();

	  if( box.shortside() != sz ) 
	    {
	      MayDay::Error("patches must be square");
	    }
	  
	  if( m_fixedBoxSize == 0 ) m_fixedBoxSize = sz;
	  else if( m_fixedBoxSize != sz ) 
	    {
	      MayDay::Error("all patches must be the same size");
	    }
	}

#ifdef CH_MPI
      {
	int rr; int ss = m_fixedBoxSize; 
	int result = MPI_Allreduce( &ss, &rr, 1, MPI_INT,
				    MPI_MAX, Chombo_MPI::comm );
	if (result != MPI_SUCCESS)
	  {
	    MayDay::Error("sorry, but I had a communcation error");
	  }
	if( m_fixedBoxSize ==  0 ) m_fixedBoxSize = rr; // an epty processor
	else if ( rr != m_fixedBoxSize )
	  {
	    MayDay::Error("all patches must be the same size");
	  }
      }
#endif
    }
}

////////////////////////////////////////////////////////////////////////
// AMRFAS<T>::define
//    - creates internal levels, keep track with m_numInternalLevels
//    - fully defines op factory
//    - allocaltes data
////////////////////////////////////////////////////////////////////////
template <class T>
void AMRFAS<T>::define( const ProblemDomain& a_coarseDomain,
			const RealVect&      a_crsDx,
			const Vector<DisjointBoxLayout>& a_grids,
			const Vector<int>&               a_refRatios,
			AMRFASOpFactory<T>&              a_factory,
			int a_numSRLevels /*=0*/, // I need to know about SR to construct comms.
			int a_num_levels /*=-1*/
			)
{
  CH_TIME("AMRFAS::define");
  
  const int num_lev_orig = (a_num_levels == -1) ? a_grids.size() : a_num_levels; // number of current (user) levels
  
  // get patch size from coarse grid
  setFixedBoxSize( a_grids[0] );
  
  // set prol order if not set yet, with m_order
  if( m_FMGProlOrderP < 0 ){
    m_FMGProlOrderP = a_factory.getOrder()/2 + 1;
  }
  if( m_ProlOrderP < 0 ){
    m_ProlOrderP = a_factory.getOrder()/2;
  }
  if (m_verbosity >= 2 )
    {
      pout() << "AMRFAS<T>::define: FMG Prol Order Poly = " << m_FMGProlOrderP << ",  Prol Order Poly = " << m_ProlOrderP << std::endl;
    }

  // passing down to the factory
  a_factory.m_smoothing_damping_factor = m_smoothing_damping_factor;
  a_factory.m_FMGProlOrderP = m_FMGProlOrderP; // not used for NEW prolongator
  a_factory.m_ProlOrderP = m_ProlOrderP;

  clear(); // noop
  
  if (m_verbosity >= 4 )
    {
      pout() << "AMRFAS<T>::define: user coarseDomain = " << a_coarseDomain << std::endl;
    }

  // count extra levels - set m_numInternalLevels
  const int top_grid_max_size = 8;
  ProblemDomain d = a_coarseDomain;
  m_numInternalLevels = 0;
  while( d.domainBox().shortside() >= top_grid_max_size*m_privateCoarseningRate )
    {
      d.coarsen(m_privateCoarseningRate);
      m_numInternalLevels++;
    }
  const int numExtraLev = m_numInternalLevels;
  const int numLevels = num_lev_orig + m_numInternalLevels;

  Vector<DisjointBoxLayout> new_grids(numLevels);
  Vector<int> new_ref_ratios(numLevels-1);

  // create coarse grids above users level zero, from fine to coarse
  m_coarseDomain = a_coarseDomain;
  RealVect crs_dx = a_crsDx;
  d = a_coarseDomain;
  m_numInternalLevels = 0;
  while( d.domainBox().shortside() >= top_grid_max_size*m_privateCoarseningRate )
    {
      Vector<Box>  newBoxes;
      Vector<int> procAssign;

      d.coarsen( m_privateCoarseningRate );
      m_coarseDomain = d;
      crs_dx *= (Real)m_privateCoarseningRate;

      domainSplit( d, newBoxes, m_fixedBoxSize ); // split it parallelize
      procAssign.resize( newBoxes.size() );
      LoadBalance( procAssign, newBoxes );
      if (m_verbosity >= 5 )
	{
	  pout() << "AMRFAS<T>::define: boxes:" << newBoxes << endl; 
	  pout() << "AMRFAS<T>::define: procs:" << procAssign << endl; 
	}

      DisjointBoxLayout dbl( newBoxes, procAssign, d );
      new_grids[numExtraLev-m_numInternalLevels-1] = dbl;
      new_ref_ratios[numExtraLev-m_numInternalLevels-1] = m_privateCoarseningRate;

      if (m_verbosity >= 5 )
  	{
  	  pout() << "\tAMRFAS::define: make new coarse domain = "<< m_coarseDomain <<std::endl;
  	}
      
      m_numInternalLevels++;
    }
  if (m_verbosity > 1 )
    {
      pout() << "AMRFAS::define: new coarse domain = " << m_coarseDomain << ", m_numInternalLevels = " << m_numInternalLevels << std::endl;
    }

  // copy original grids and refs
  for (int lev = 0; lev < num_lev_orig; lev++)
    {
      new_grids[numExtraLev+lev] = a_grids[lev];
    }
  for (int lev = 0; lev < num_lev_orig - 1; lev++)
    {
      new_ref_ratios[numExtraLev+lev] = a_refRatios[lev];
    }

  // fully define factory with these levels
  a_factory.define( m_coarseDomain, crs_dx, new_grids, new_ref_ratios, a_numSRLevels );

  // set number of iterations for coarse grid solver
  m_coarse_its = 20*m_coarseDomain.domainBox().longside();

  // create my data
  m_residual.  resize( numLevels );
  m_temp.      resize( numLevels );
  m_rhsC.      resize( numLevels - a_numSRLevels );
  m_solC.      resize( numLevels - a_numSRLevels ); 
  m_restrictCopier.resize( numLevels - a_numSRLevels ); // [0] is not used
  m_op.        resize( numLevels );

  for (int lev = 0; lev < numLevels; lev++)
    {
      bool isSR = (lev >= numLevels-a_numSRLevels);

      m_residual[lev] = RefCountedPtr<T>(new T());
      m_temp[lev]     = RefCountedPtr<T>(new T()); // only needed to compute residual
      if( !isSR )
	{
	  m_rhsC[lev]     = RefCountedPtr<T>(new T());
	  m_solC[lev]     = RefCountedPtr<T>(new T()); 
	}

      m_op[lev] = a_factory.AMRNewOp( lev, new_grids[lev], isSR );

      // ops have solvers so need to pass parameters down, 
      m_op[lev]->m_smoother = m_smoother; 
    }

}

////////////////////////////////////////////////////////////////////////
// AMRFAS<T>::init
////////////////////////////////////////////////////////////////////////
template <class T>
void AMRFAS<T>::init( const Vector<RefCountedPtr<T> > &a_phi, const Vector<RefCountedPtr<T> > &a_rhs )
{
  CH_TIME("AMRFAS::init");

  int l_base = 0, l_max = a_phi.size() - 1, numNonSRLevels = m_rhsC.size();

  // set new vectors, these are all ref counted so do not need to delete them in clear()
  for (int i = l_base; i <= l_max ; i++)
    {
      AMRFASOp<T>& op = *(m_op[i]);
      op.create( *m_residual[i], *a_rhs[i] ); // these are strange on SR levels
      op.create( *m_temp[i], *a_phi[i] );  // only used for residual
    }

  // create inter-level copiers
  for (int i = l_base+1; i < numNonSRLevels ; i++)
    {
      AMRFASOp<T>& op = *(m_op[i]);
      int r = op.refToCoarser();
      int nGhost = op.getOrder()/2;
      // restriction copier amd coarse cover, no ghosts, low order
      {
	int ncomp = a_rhs[i]->nComp();
	IntVect ghostVect = a_rhs[i]->ghostVect();
	const DisjointBoxLayout& dbl = a_rhs[i]->disjointBoxLayout();
	CH_assert(dbl.coarsenable(r));
	DisjointBoxLayout dblCoarsenedFine;
	coarsen( dblCoarsenedFine, dbl, r );
	m_rhsC[i]->define( dblCoarsenedFine, ncomp, ghostVect );
      }

      m_restrictCopier[i].define( m_rhsC[i]->disjointBoxLayout(), a_rhs[i-1]->disjointBoxLayout(), IntVect::Zero );
      // The number of coarse ghost required -- this is a hack, should query prolongators
      int maxpoly = (m_FMGProlOrderP>m_ProlOrderP) ? m_FMGProlOrderP :  m_ProlOrderP, n_cover_ghost = 1;
      switch(maxpoly)
	{
	case 0:
	  n_cover_ghost = 0; break;
	case 1:
	  n_cover_ghost = 1; break;
	case 2:
	  n_cover_ghost = 1; break;
	case 3:
	  n_cover_ghost = 2; break;
	case 4:
	  n_cover_ghost = 2; break;
	default:
	  MayDay::Abort("unknow poly undefined");
	}
      if (n_cover_ghost<nGhost) n_cover_ghost = nGhost;
      IntVect HOghostVect = n_cover_ghost*IntVect::Unit; 

      // prolongation copier and coarse cover
      const DisjointBoxLayout& dbl = a_phi[i]->disjointBoxLayout();
      DisjointBoxLayout dblCoarsenedFine;
      coarsen( dblCoarsenedFine, dbl, r );
      m_solC[i]->define( dblCoarsenedFine, a_phi[i]->nComp(), HOghostVect );
      op.m_HOCopier.define( a_phi[i-1]->disjointBoxLayout(),
			    m_solC[i]->disjointBoxLayout(),
			    HOghostVect
			    );
    }
}

////////////////////////////////////////////////////////////////////////
// AMRFAS<T>::solve
////////////////////////////////////////////////////////////////////////
template<class T>
Real 
AMRFAS<T>::solve( Vector<RefCountedPtr<T> >& a_phi, 
		  Vector<RefCountedPtr<T> >& a_rhs,
		  int *a_status
		  )
{
  CH_TIME("AMRFAS::solve");

  // convert to internal vectors: m_opInternal[k+m_numInternalLevels] = m_opUser[k], 0 <= k < m_numLevels;
  const int numUserLevels = a_phi.size(), nLevTot = numUserLevels+m_numInternalLevels;
  const int ncomp = a_phi[0]->nComp();
  Vector<RefCountedPtr<T> > phi(nLevTot);
  Vector<RefCountedPtr<T> > rhs(nLevTot);

  // copy users vectors in
  for( int ilev = m_numInternalLevels, userLev = 0 ; ilev < nLevTot ; ilev++, userLev++ )
    {
      phi[ilev] = a_phi[userLev];
      rhs[ilev] = a_rhs[userLev];
    }

  // create internal coarse levels
  for( int ilev = m_numInternalLevels-1 ; ilev >= 0 ; ilev-- )
    {
      phi[ilev] = RefCountedPtr<T>(new T());
      rhs[ilev] = RefCountedPtr<T>(new T());

      phi[ilev]->define( m_op[ilev]->m_grid, ncomp, a_phi[0]->ghostVect() ); 
      rhs[ilev]->define( m_op[ilev]->m_grid, ncomp, IntVect::Zero );
    }

  if( !m_init ) 
    {
      init( phi, rhs );
      m_init = true;
    }

  // restrict extra coarse grids rhs and states
  for( int ilev = m_numInternalLevels ; ilev > 0 ; ilev-- )
    {
      m_op[ilev]->AMRRestrict( *rhs[ilev-1],    // output
			       *rhs[ilev],      // input
			       *m_rhsC[ilev], 
			       m_restrictCopier[ilev] );

      m_op[ilev]->AMRRestrict( *phi[ilev-1],    // output
			       *phi[ilev],      // input
			       *m_rhsC[ilev], 
			       m_restrictCopier[ilev] );
      
      // this restricts (linear) state that is not set with recomputeState
      m_op[ilev-1]->restrictState( m_op[ilev], m_restrictCopier[ilev] );
    }

  return solveNoInit( phi, rhs, a_status );
}

////////////////////////////////////////////////////////////////////////
// AMRFAS<T>::solveNoInit
////////////////////////////////////////////////////////////////////////
template<class T>
Real AMRFAS<T>::solveNoInit( Vector<RefCountedPtr<T> >& a_phi,
			     const Vector<RefCountedPtr<T> >& a_rhs,
			     int *a_status )
{
  CH_assert(a_rhs.size() == a_phi.size());
  CH_TIMERS("AMRFAS::solveNoInit");
  CH_TIMER("AMRFAS::cycle", vtimer);

  //bool outputIntermediates = false;
  int l_max = a_phi.size() - 1;

  Real initial_norm = 0.;
  if( !m_avoid_norms )
    {
      CH_TIME("Initial AMR Residual"); // could just take norm of rhs
      initial_norm = computeAMRResidual_private( m_temp, a_phi, a_rhs );
      if(initial_norm!=initial_norm)
	{
	  MayDay::Error("Nan residual norm");
	}
    }

  Real rnorm = initial_norm;
  Real norm_last = 2*initial_norm; // keep this from triggering anything if no norms

  int iter=0;
  if (m_verbosity >= 2 && !m_avoid_norms )
    {
      pout() << "    AMRFAS:: iteration = " << iter << ", residual norm = " << rnorm << std::endl;
    }
  
  bool goNorm = rnorm > m_atol || m_avoid_norms;                   //iterate if norm is not small enough
  bool goRedu = rnorm > m_rtol*initial_norm || m_avoid_norms;     //iterate if initial norm is not reduced enough
  bool goIter = iter < m_max_iter;                                 //iterate if iter < max iteration count
  bool goHang = iter < m_min_iter || abs(rnorm-norm_last)/norm_last > m_stagnate; //iterate if we didn't hang
  bool goMin = iter < m_min_iter ; // iterate if iter < min
  while (goMin || (goIter && goRedu && goHang && goNorm))
    {
      norm_last = rnorm;

      //this generates a correction from the current residual
      CH_START(vtimer);
      if( m_type == FAS_FULL && iter == 0 )
	{
	  FMG( a_phi, a_rhs );
	}
      else if( m_type == FAS_VCYCLE || m_type == FAS_FULL )
	{
	  if( iter == 0 ) // fine grid does not get residual touched
	    {
	      // fine grid rhs should not be touched, so put it in once
	      m_op[l_max]->assignLocal( *m_residual[l_max], *a_rhs[l_max] ); 
	      // set fine grid C-F ghosts, coarse grids set in VCycle
	      if( l_max > 0 )
		{
		  m_op[l_max]->CFInterp( *a_phi[l_max], *a_phi[l_max-1] );
		}
	    }
	  VCycle( a_phi, a_rhs, l_max, l_max );
	}
      else if( m_type == FAS_SMOOTH )
	{
	  m_op[l_max]->assignLocal( *(m_residual[l_max]), *(a_rhs[l_max]) );
	  m_op[l_max]->smooth( a_phi[l_max], m_residual[l_max], 0, 4 ); //m_coarse_its, true );
	}
      else 
	{
	  MayDay::Error("unknown FAS type");
	}
      CH_STOP(vtimer);

      // For solvers with accuracy higher than 2nd order
      //  consistency between levels has to be explicitly enforced. -- does this go away!!!
      // if (m_op[0]->orderOfAccuracy()>2)
      //   {
      //     for (int ilev=l_max; ilev>l_base; ilev--)
      //       {
      //         m_op[ilev]->enforceCFConsistency(*a_phi[ilev-1], *a_phi[ilev]);
      //       }
      //   }

      // recompute residual
      iter++;
      goIter = iter < m_max_iter;                        //keep iterating if iter < max iteration count
      goMin = iter < m_min_iter ;                         // keep iterating if iter < min

      if( !m_avoid_norms )
	{
	  rnorm = computeAMRResidual_private( m_temp, a_phi, a_rhs );	  
	  goNorm = rnorm > m_atol;                        //keep iterating if norm is not small enough
	  goRedu = rnorm > m_rtol*initial_norm;          //keep iterating if initial norm is not reduced enough
	  goHang = iter < m_min_iter || abs(rnorm-norm_last)/norm_last > m_stagnate; //keep iterating if we didn't hang
	  if(rnorm!=rnorm)
	    {
	      MayDay::Error("Nan residual norm");
	    }
	}
      
      if( m_verbosity > 1 && !m_avoid_norms)
	{
	  pout() << "    AMRFAS:: iteration = " << iter << ", residual norm = " << rnorm;
	  if (rnorm > 0.0)
	    {
	      pout() << ", rate = " << rnorm/norm_last;
	    }
	  pout() << std::endl;
	}
    }

  int exitStatus = int(!goRedu) + int(!goIter)*2 + int(!goHang)*4 + int(!goNorm)*8;
  if (m_verbosity > 1 && !m_avoid_norms)
    {
      pout() << "    AMRFAS:: iteration = " << iter << ", residual norm = " << rnorm << std::endl;
    }
  if (m_verbosity > 1)
    {
      if (!goIter && goRedu && goNorm) // goRedu=T, goIter=F, goHang=?, goNorm=T
        { 
	  //exitStatus == 0 + 2 + 0|4 + 0 = 2|6
	  pout() << "    AMRFAS:: WARNING: Exit because max iteration count exceeded" << std::endl;
        }
      if (!goHang && goRedu && goNorm) // goRedu=T, goIter=?, goHang=F, goNorm=T
        { 
	  //exitStatus == 0 + 0|2 + 4 + 0 = 4|6
	  pout() << "    AMRFAS:: WARNING: Exit because of solver hang" << std::endl;
        }
      if (m_verbosity > 2)
        {
          pout() << "    AMRFAS:: exitStatus = " << exitStatus << std::endl;
        }
    }

  if( a_status ) *a_status = exitStatus;

  return rnorm;
}

////////////////////////////////////////////////////////////////////////
// AMRFAS<T>::FMG
////////////////////////////////////////////////////////////////////////
template<class T>
void 
AMRFAS<T>::FMG( Vector<RefCountedPtr<T> > & a_phi,
		const Vector<RefCountedPtr<T> > & a_rhs,
		int a_lmax /* = -1 */)
{
  int l_base = 0, l_max = (a_lmax == -1) ? a_phi.size() - 1 : a_lmax;

  // coarse grid
  m_op[l_base]->assignLocal( *m_residual[l_base], *a_rhs[l_base] ); // this is really the RHS
  VCycle( a_phi, a_rhs, l_base, l_max );
  
  for( int ilev = l_base+1 ; ilev <= l_max ; ilev++ )
    {
      if( m_verbosity > 5 )
	{
	  pout() << "AMRFAS::FMG visit level " << ilev << std::endl;
	}

      // FMG interpolation -- high order -- this can clobber a_phi[ilev]
      m_op[ilev]->AMRFMGProlong( *a_phi[ilev], *a_phi[ilev-1],
				 *m_solC[ilev], m_op[ilev-1] ); 
      
      // V-cycle
      m_op[ilev]->assignLocal( *(m_residual[ilev]), *(a_rhs[ilev]) ); // m_residual is really the RHS
      // set fine grid C-F ghosts, coarse grids set in VCycle
      m_op[ilev]->CFInterp( *a_phi[ilev], *a_phi[ilev-1] );

      for(int i=0;i<m_num_vcycles;i++)
	{
	  VCycle( a_phi, a_rhs, ilev, ilev );
	}
    }
}
#include <iomanip>
////////////////////////////////////////////////////////////////////////
// AMRFAS<T>::VCycle
//   'm_residual' is really the full RHS ala FAS.
//   'a_rhs' is only needed for coarse grids that do not cover the fine (ie, AMR) - fix this!
////////////////////////////////////////////////////////////////////////
template<class T>
Real AMRFAS<T>::VCycle( Vector<RefCountedPtr<T> >& a_phi,
			const Vector<RefCountedPtr<T> >& a_rhs,
			int a_ilev,
			int a_max
			)
{
  int l_base = 0;

  //pout() << "\t\t\t" << a_ilev << ")AMRFAS::VCycle |b| = " << m_op[a_ilev]->localMaxNorm( *m_residual[a_ilev] ) << std::endl;

  if (a_ilev == l_base )
    {
      CH_TIME("AMRFAS::VCycle:coarse_grid_solver");
      m_op[a_ilev]->computeState( a_phi[a_ilev], (RefCountedPtr<T>)0, (RefCountedPtr<T>)0 ); 
      // exact solver 
      m_op[a_ilev]->smooth( a_phi[a_ilev], m_residual[a_ilev], 0, m_coarse_its, true );
     }
  else
    {
      m_op[a_ilev]->computeState( a_phi[a_ilev], a_phi[a_ilev-1], (RefCountedPtr<T>)0 ); // redoing state at a high level

      //============= Downsweep ========================
      m_op[a_ilev]->smooth( a_phi[a_ilev], 
			    m_residual[a_ilev], 
			    0, // this sets this levels C-F ghosts on this level, done before Vcycle
			    m_pre ); 
      
      // Recompute residual on next coarser level for region NOT covered by this level.  Before restriction of u.
      m_op[a_ilev-1]->residual( *(m_residual[a_ilev-1]),  // set coarse grid residual, will soon be current resid
				*a_phi[a_ilev-1],
				*(a_rhs[a_ilev-1]), // real RHS for region not covered
				0, // (a_ilev-1 == l_base) ? 0 : &(*a_phi[a_ilev-2]), // set coarse grid C-F ghosts -- does not seem to make any difference
				&(*a_phi[a_ilev]),
				&(*m_op[a_ilev]) );
      
      // Compute the restriction of the solution to the coarser level a_phi[a_ilev-1] -- R(u_f) 
      m_op[a_ilev]->AMRRestrict( *a_phi[a_ilev-1],     // output
				 *a_phi[a_ilev],       // input
				 *m_rhsC[a_ilev], 
				 m_restrictCopier[a_ilev] );
      
      // Dedek's correction, add it to whole level at end
      m_op[a_ilev-1]->apply( *(m_temp[a_ilev-1]),  // coarse grid Au
			     *a_phi[a_ilev-1],
			     0 // this was iterpolated earlier so no need to do it again
			     );

      // Compute fine residual: temp <= f - L_f(u_f) 
      m_op[a_ilev]->residual( *(m_temp[a_ilev]), 
			      *a_phi[a_ilev],
			      *(m_residual[a_ilev]), // this is incoming RHS 'f'
			      0, // C-F interp is done
			      0, // not supposed to add refluxing here (Dedek)
			      0);

      // Compute the restriction of the _residual_ to the coarser level resC -- R(f - L_f(u_f))
      // This overwrites the supported coarse grid residual from above
      m_op[a_ilev]->AMRRestrict( *m_residual[a_ilev-1],      // output
				 *m_temp[a_ilev],            // input
				 *m_rhsC[a_ilev],            // coarse cover buffer 
				 m_restrictCopier[a_ilev]    // m_rhsC --> m_residual
				 );
      
      // add correction to whole level (same term in Dudek)
      m_op[a_ilev-1]->axby( *(m_residual[a_ilev-1]), // output: coarse grid residual
			    *(m_residual[a_ilev-1]), // coarse grid RHS
			    *(m_temp[a_ilev-1]), 
			    1.0, 1.0 );
      // end of coarse grid residual

      {
	// store correction R(u_f)
	T R_u_f;
	m_op[a_ilev-1]->create( R_u_f, *a_phi[a_ilev-1]);
	m_op[a_ilev-1]->assignLocal( R_u_f, *(a_phi[a_ilev-1]) ); 

	//============finish Compute residual for the next coarser level======
	for (int img = 0; img < m_num_cycles; img++)
	  {
	    VCycle( a_phi, a_rhs, a_ilev-1, a_max );
	  }

	// subtract off initial solution to get an increment (a_phi[a_ilev-1] is an increment!)
	m_op[a_ilev-1]->axby( *a_phi[a_ilev-1], *a_phi[a_ilev-1], R_u_f, 1.0, -1.0 );

	//================= Upsweep ======================
	// increment the correction with coarser version
	m_op[a_ilev]->AMRProlong( *(a_phi[a_ilev]), *(a_phi[a_ilev-1]), 
				  *(m_solC[a_ilev]), m_op[a_ilev-1] ); 

	// add initial solution back to get full solution
	m_op[a_ilev-1]->axby( *a_phi[a_ilev-1], *a_phi[a_ilev-1], R_u_f, 1.0, 1.0 );
      }

      m_op[a_ilev]->smooth( a_phi[a_ilev], 
			    m_residual[a_ilev], 
			    &(*a_phi[a_ilev-1]), // redo C-F ghost after V-cycle
			    m_post );
    }

  return 1.0; // norm needs to be done -- TODO
}

////////////////////////////////////////////////////////////////////////
// AMRFAS<T>::computeAMRResidual
////////////////////////////////////////////////////////////////////////
template <class T>
Real AMRFAS<T>::computeAMRResidual( Vector<RefCountedPtr<T> >&       a_resid,
				    Vector<RefCountedPtr<T> >&       a_phi,
				    const Vector<RefCountedPtr<T> >& a_rhs
				    )
{
  CH_TIME("AMRFAS::computeAMRResidual");

  return computeAMRResidual_private( a_resid, a_phi, a_rhs, m_numInternalLevels );
}
////////////////////////////////////////////////////////////////////////
// AMRFAS<T>::computeAMRResidual_private
////////////////////////////////////////////////////////////////////////
template <class T>
Real AMRFAS<T>::computeAMRResidual_private( Vector<RefCountedPtr<T> >&       a_resid,
					    Vector<RefCountedPtr<T> >&       a_phi,
					    const Vector<RefCountedPtr<T> >& a_rhs, 
					    const int a_base
					    )
{
  CH_TIME("AMRFAS::computeAMRResidual_private");
  Real rnorm = 0.;
  Real localNorm = 0.;
  int l_max = a_phi.size() - 1;  

  for (int ilev = 0 ; ilev <= l_max; ilev++)
    {
      // can recompute state & do exchange - this could be skipped in some cases ...
      m_op[ilev+a_base]->computeState( a_phi[ilev], 
				       (ilev == 0) ? (RefCountedPtr<T>)0 : a_phi[ilev-1],
				       (ilev == l_max) ? (RefCountedPtr<T>)0 : a_phi[ilev+1]
				       );

      m_op[ilev+a_base]->residual( *a_resid[ilev],
				   *a_phi[ilev],
				   *a_rhs[ilev],
				   0, // this does C-F interp which is done by coarser grid
				   (ilev == l_max) ? 0 : &(*a_phi[ilev+1]), // this does C-F interp on finer (next) grid
				   (ilev == l_max) ? 0 : &(*m_op[ilev+a_base+1])
				   );

      if (ilev == l_max)
	{
	  localNorm = m_op[ilev+a_base]->localMaxNorm(*a_resid[ilev]);
	}
      else
	{
	  if( ilev+1+a_base < m_rhsC.size() )
	    {
	      m_op[ilev+a_base]->zeroCovered(*a_resid[ilev],*m_rhsC[ilev+1+a_base],m_restrictCopier[ilev+1+a_base]);
	      localNorm = m_op[ilev+a_base]->localMaxNorm( *a_resid[ilev] );
	    }
	  else // SR 
	    {
	      localNorm = 0.0;
	    }
	}

      // const int normType = 0;
      // if (normType==2) localNorm = m_op[ilev+a_base]->norm(*a_resid[ilev],normType);
      rnorm = Max( localNorm, rnorm ); // takeing the max level norm
    }

#ifdef CH_MPI
  Real recv;
  MPI_Allreduce( &rnorm, &recv, 1, MPI_CH_REAL,MPI_MAX, Chombo_MPI::comm );
  rnorm = recv;
#endif
  
  return rnorm; // if a_computeNorm is false, then this just returns zero.
}

#include "NamespaceFooter.H"
#endif
