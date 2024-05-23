#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// MFA, December 2012

#include "AMRFASSR.H"
#include "AMRIO.H"

//*******************************************************
// AMRFASSR Implementation
//*******************************************************

////////////////////////////////////////////////////////////////////////
// AMRFASSR::setParameters
////////////////////////////////////////////////////////////////////////
void AMRFASSR::setParameters( const char *name )
{
  AMRFAS<LevelData<FArrayBox> >::setParameters( name );

  {
    ParmParse ppSolver( name );
    ppSolver.query( "num_buffer_cells", m_numBufferCells0 ); 
    ppSolver.query( "num_sr_levels", m_numSRLevels );
    m_fmg_pre = m_pre;
    ppSolver.query( "num_fmg_pre", m_fmg_pre ); 
    
    int num_coarse_update_fine = (getNumBufferCells(m_numSRLevels) + m_numGhostCells)/m_refRatio; // floor  
    pout() << "AMRFASSR::setParameters effective nesting region on coarse grid for finest = "<< (Real)getNumBufferCells(m_numSRLevels-1) - (Real)(getNumBufferCells(m_numSRLevels) + m_numGhostCells)/(Real)m_refRatio <<", updating "<< num_coarse_update_fine << " coarse cells from fine" << std::endl;
  }
}

////////////////////////////////////////////////////////////////////////
// AMRFASSR::copyUFromCover
////////////////////////////////////////////////////////////////////////
void AMRFASSR::copyUFromCover( LevelData<FArrayBox> &R_u_f,
			       const RefCountedPtr<LevelData<FArrayBox> > a_crsCover )
{
  int nc = a_crsCover->nComp();  CH_assert(nc%2==0);
  nc = nc / 2;
  
  R_u_f.define( a_crsCover->disjointBoxLayout(), 
		nc, // only Ru
		a_crsCover->ghostVect() );

  for (DataIterator dit= R_u_f.dataIterator(); dit.ok(); ++dit)
    {
      FArrayBox& rufab =  R_u_f[dit]; Box box = rufab.box();
      rufab *= 0.;
      rufab.plus( (*a_crsCover)[dit], box, 0, 0, nc );
    }
}

////////////////////////////////////////////////////////////////////////
// AMRFASSR::solve
//      need to create phi & RHS LDF for SR levels, and fill in RHS for all levels
//      this creates SR levels and returns them to the user (getting AMR into the solver)
////////////////////////////////////////////////////////////////////////
Real
AMRFASSR::solve( Vector<RefCountedPtr<LevelData<FArrayBox> > >& a_phi, 
		 Vector<RefCountedPtr<LevelData<FArrayBox> > >& a_rhs,
		 int *a_status
		 )
{
  CH_TIME("AMRFASSR::solve");

  // create phi & RHS LDF - add to a_phi, a_rhs
  int numCrsLevels = a_phi.size(), ii;
  RefCountedPtr<AMRFASOp<LevelData<FArrayBox> > > op = getOp(numCrsLevels-1); // get fine op with users index - fine grid 
  Real dx;
  const int ncomp = a_phi[0]->nComp();

  a_phi.resize(numCrsLevels+m_numSRLevels); // adding SR levels
  a_rhs.resize(numCrsLevels+m_numSRLevels);

  for (int i = 0, sr_ilev = 1; i < m_numSRLevels; i++)
    {
      IntVect ghostVect=(getNumBufferCells(sr_ilev)+m_numGhostCells)*IntVect::Unit;
      IntVect bufferVect = getNumBufferCells(sr_ilev)*IntVect::Unit; 

      DisjointBoxLayout &grid = m_SRGrids[i];
      RefCountedPtr<LevelData<FArrayBox> > srPhi = 
	RefCountedPtr<LevelData<FArrayBox> >(new LevelData<FArrayBox>(grid,a_phi[0]->nComp(),ghostVect*IntVect::Unit));
      a_phi[ numCrsLevels + i ] = srPhi;

      RefCountedPtr<LevelData<FArrayBox> > srRHS = 
	RefCountedPtr<LevelData<FArrayBox> >(new LevelData<FArrayBox>(grid,a_phi[0]->nComp(),bufferVect*IntVect::Unit));
      a_rhs[ numCrsLevels + i ] = srRHS;
    }
  
  // set RHS and phi=0 for all user and SR levels (but not solver coarse grids)
  for (ii = 0, dx = op->dx() ; ii < a_phi.size() ; ii++, dx /= m_refRatio )
    {
      const DisjointBoxLayout& dbl = a_phi[ii]->disjointBoxLayout();
      DataIterator dit = dbl.dataIterator();
      for (dit.begin(); dit.ok(); ++dit)
	{
	  FArrayBox& rhsFab = (*a_rhs[ii])[dit];
	  const Box fbox = rhsFab.box(); // fill entire box
	  m_rhsFunc( rhsFab, fbox, dx );
	  FArrayBox& phifab = (*a_phi[ii])[dit];
	  const Box pbox = phifab.box(); // fill entire box
	  phifab.setVal(0.0,pbox,0,ncomp);
	}
    }

  // solve
  return AMRFAS<LevelData<FArrayBox> >::solve( a_phi, a_rhs, a_status );
}

////////////////////////////////////////////////////////////////////////
// AMRFASSR::define
//      add SR levels to a_grids & a_refRatios
////////////////////////////////////////////////////////////////////////
void
AMRFASSR::define( const ProblemDomain&                   a_coarseDomain,
		  const RealVect&                        a_crsDx,
		  const Vector<DisjointBoxLayout>&       a_grids,
		  const Vector<int>&                     a_refRatios,
		  AMRFASOpFactory<LevelData<FArrayBox> > &a_factory,
		  int a_dummy,
		  int a_num_levels /*=-1*/ )
{
  CH_assert(a_num_levels==-1);
  CH_TIME("AMRFASSR::define");

  // get patch size from coarse grid
  setFixedBoxSize( a_grids[0] );

  const int sr_base_idx = a_grids.size();
  int idx, idx2,idx3, npe, finegridfactor = 1;                   // don't split first SR grid
  const int fixedBoxSize = m_refRatio*m_fixedBoxSize;              // don't split first level so increase box size
#ifdef CH_MPI
  MPI_Comm_size( MPI_COMM_WORLD, &npe );
#else
  npe = 1;
#endif

  m_numGhostCells = a_factory.getOrder()/2; // the hack for ghosts

  m_inCopier.resize( m_numSRLevels );
  m_outCopier.resize( m_numSRLevels );
  m_crsCover.resize( m_numSRLevels );
  m_SRGrids.resize( m_numSRLevels );

  // get coarse domain - 1
  ProblemDomain cdom = a_coarseDomain;
  for (int i = 1; i < sr_base_idx; i++)
    {
      cdom.refine(a_refRatios[i-1]);
    }
  ProblemDomain fdom = cdom;

  // need to append new vector of grids to send to parent
  Vector<DisjointBoxLayout> grids( a_grids );
  Vector<int> refRatios( a_refRatios );
  grids.resize( a_grids.size() + m_numSRLevels );
  refRatios.resize( a_grids.size() + m_numSRLevels - 1 );

  for (int clev = sr_base_idx-1, sr_idx = 0; sr_idx < m_numSRLevels ; clev++, sr_idx++ )
    {
      // refine domain
      cdom.refine(m_refRatio); 
      fdom.refine(m_refRatio); 

      refRatios[clev] = m_refRatio; // spoof parent define()

      // make a BoxLayout of SR level
      const DisjointBoxLayout& cdbl = a_grids[clev]; //
      Vector<int> coarseProcs = cdbl.procIDs();
      const int nCrsBoxes = coarseProcs.size(), nFineBoxes = nCrsBoxes*finegridfactor;
      Vector< Box > fineBoxes(finegridfactor*nCrsBoxes);
      Vector<int>   procAssign(finegridfactor*nCrsBoxes);

      LayoutIterator lit = cdbl.layoutIterator(); // we are iterating over the global grid -- would like to do this locally!
      for (lit.begin(),idx=0,idx2=0,idx3=0; lit.ok(); ++lit)
	{
	  Box fbox = cdbl[lit];
	  fbox.refine( m_refRatio );
	  Vector< Box > patch_fineboxes;
	  domainSplit( fbox, patch_fineboxes, fixedBoxSize, fixedBoxSize );
	  CH_assert(patch_fineboxes.size() == finegridfactor);

	  // I can decide to not refine here, see if cbox wants to be created, AMR here

	  if( nFineBoxes <= npe && npe%nFineBoxes == 0 )
	    {
	      int factor = npe/nFineBoxes;
	      int proc = coarseProcs[idx3++];
	      for(int i=0;i<patch_fineboxes.size();i++) procAssign[idx++] = proc + i*factor;
	    }
	  // make cover, copy boxes in
	  for(int i=0;i<finegridfactor;i++) 
	    {
	      Box box = patch_fineboxes[i];
	      fineBoxes[idx2++] = box;
	    }
	}
      CH_assert(fineBoxes.size()==procAssign.size());

      if( nFineBoxes <= npe && npe%nFineBoxes == 0 )
	{
	  pout() << "\tAMRFASSR::define: make new (perfect) fine grid. procs = "<< procAssign <<std::endl;
	}
      else
	{
	  LoadBalance( procAssign, fineBoxes );
	  pout() << "\tAMRFASSR::define: Warning, using LoadBalance"<< std::endl;
	}

      DisjointBoxLayout fgrid( fineBoxes, procAssign, fdom );
      m_SRGrids[sr_idx] = fgrid; // cache grid
      grids[clev+1] = fgrid;     // spoof parent define()

      // make cover patches
      DisjointBoxLayout coverGrid;
      coarsen( coverGrid, fgrid, m_refRatio );

      // "ghost" region for cover is m_numBufferCells + m_numGhostCells
      IntVect coverGhostVect = getNumBufferCells(sr_idx)*IntVect::Unit; 
      CH_assert(getNumBufferCells(sr_idx) >= (getNumBufferCells(sr_idx+1)+m_numGhostCells)/m_refRatio + m_numGhostCells ); // NB_c >= update_region + NG
      m_inCopier[sr_idx].define( cdbl, coverGrid, coverGhostVect );
      if( sr_idx == 0 )
	{
	  m_outCopier[sr_idx].define( coverGrid, cdbl, IntVect::Zero );
	}
      else
	{
	  IntVect crsRHSGhostVect = (getNumBufferCells(sr_idx+1)+m_numGhostCells)/m_refRatio*IntVect::Unit; 
	  m_outCopier[sr_idx].define( coverGrid, cdbl, crsRHSGhostVect ); // not correct!!!
	}

      RefCountedPtr<LevelData<FArrayBox> > crsLDF = 
	RefCountedPtr<LevelData<FArrayBox> >(new LevelData<FArrayBox>( coverGrid, 
								       2*a_factory.nComp(), // for u & resid
								       coverGhostVect) );
      m_crsCover[sr_idx] = crsLDF;

      // now lets split
      finegridfactor = 1;
      for(int i=0;i<SpaceDim;i++) finegridfactor *= m_refRatio; 
    }

  // call base class to create coarse grids and all internal levels stuff
  AMRFAS<LevelData<FArrayBox> >::define( a_coarseDomain,
					 a_crsDx,
					 grids,
					 refRatios,
					 a_factory,
					 m_numSRLevels
					 );
}

////////////////////////////////////////////////////////////////////////
// AMRFASSR::FMG
////////////////////////////////////////////////////////////////////////
void 
AMRFASSR::FMG( Vector<RefCountedPtr<LevelData<FArrayBox> > > & a_phi,
	       const Vector<RefCountedPtr<LevelData<FArrayBox> > > & a_rhs,
	       int a_dummy )
{
  CH_assert(a_dummy == -1); 
  const int l_max = a_phi.size() - 1, nc = a_phi[0]->nComp();
  const int crs_ilev = l_max - m_numSRLevels;
  Interval phiInterval = a_phi[0]->interval();
  Interval resInterval = Interval( nc, 2*nc - 1 );

  // coarse grid solve
  AMRFAS<LevelData<FArrayBox> >::FMG( a_phi, a_rhs, crs_ilev );

  // if not SR levels, done
  if( crs_ilev == l_max ) return;

  // bootstrap to a normal FMG loop  
  // copy-in: coarse grid to cover
  a_phi[crs_ilev]->copyTo( phiInterval, 
			   *m_crsCover[0], 
			   phiInterval,
			   m_inCopier[0] 
			   );

  // continue FMG in SR levels
  for( int ilev = crs_ilev + 1, sr_idx = 0 ; ilev <= l_max ; ilev++, sr_idx++ )
    {
      if( m_verbosity > 1 )
	{
	  pout() << "\t\t" << ilev << ") AMRFASSR::FMG visit SR level " << sr_idx << std::endl;
	}
      CH_assert(ilev == l_max); // only one SR now

      // SR1: FMG prolongation + SR2
      m_op[ilev]->SR1_kernel( *a_phi[ilev], // fine grid
			      *a_rhs[ilev], // would copy to residaul but no init guess here
			      *m_crsCover[sr_idx], // out: [R(u_f), A_H(R(u_f)) + R(f-A_h(u_h))]
			      m_fmg_pre,    // special smoothing count for first pre
			      m_refRatio,
			      m_numGhostCells,
			      getNumBufferCells(sr_idx+1)
			      );

      // save R(u)
      LevelData<FArrayBox> R_u_f;
      copyUFromCover( R_u_f, m_crsCover[sr_idx] );

      // coarse grid V-cycle, this must copy-out & copy-in
      SRVCycle( a_phi, a_rhs, ilev-1 );

      // m_residual is really the RHS, 
      m_op[ilev]->assignLocal( *(m_residual[ilev]), *(a_rhs[ilev]) ); 

      // post leg
      // SR3: subtract Ru from u_H, prlongate + increment u_h, smooth
      m_op[ilev]->SR3_kernel( *a_phi[ilev], 
			      *m_residual[ilev], 
			      R_u_f, // in
			      *m_crsCover[sr_idx], // out
			      m_post, 
			      m_refRatio,
			      m_numGhostCells );

      if( ilev != l_max )
	{
	  pout() << "\tNOT USED!!!!!----- AMRFASSR::FMG (last) COPY-IN: CLOBBER sr_idx = " << sr_idx << ", lev = " << ilev << std::endl;
	  // copy-in, an SR kernel always follows this, so get ready for it
	  a_phi[ilev]->copyTo( phiInterval, 
			       *m_crsCover[sr_idx], 
			       phiInterval, 
			       m_inCopier[sr_idx] 
			       );
	}
    }
 }

//extern LevelData<FArrayBox> *pldf_exact_glob; // debug

////////////////////////////////////////////////////////////////////////
// AMRFASSR::VCycle
////////////////////////////////////////////////////////////////////////
Real AMRFASSR::VCycle( Vector<RefCountedPtr<LevelData<FArrayBox> > >& a_phi,
		       const Vector<RefCountedPtr<LevelData<FArrayBox> > >& a_rhs,
		       int a_ilev,
		       int a_dummy // don't do C-F now so not needed
		       )
{
  const int l_max = a_phi.size() - 1, sr_idx = a_ilev - l_max + m_numSRLevels - 1;

  if( sr_idx < 0 )
    {
      // in coarse grid solver
      AMRFAS<LevelData<FArrayBox> >::VCycle( a_phi, a_rhs, a_ilev, a_ilev );
    }
  else if( a_ilev == l_max )
    {
      CH_assert(sr_idx==0); // deeper levels stay in SRVCycle


      // set error -- debug
      // if( pldf_exact_glob ) {
      // 	LevelData<FArrayBox> err;
      // 	err.define( a_phi[a_ilev]->disjointBoxLayout(), 
      // 		    1,(getNumBufferCells(sr_idx+1)+m_numGhostCells)*IntVect::Unit);
      // 	Real ldx = m_op[a_ilev]->dx();
      // 	for (DataIterator dit = a_phi[a_ilev]->dataIterator(); dit.ok(); ++dit)
      // 	  {
      // 	    const FArrayBox &exact = (*pldf_exact_glob)[dit];
      // 	    FArrayBox &error = err[dit];
      // 	    const FArrayBox &phi = (*a_phi[a_ilev])[dit];
      // 	    Box region = exact.box(); // ghosted
      // 	    for (BoxIterator bit(region); bit.ok(); ++bit)
      // 	      {
      // 		const RealVect offset = bit(); //-domain.smallEnd();
      // 		const RealVect x = ldx*(0.5+offset);
      // 		error(bit()) = exact(bit()) - phi(bit());
      // 	      }
      // 	    writeFABname( &error, "error0FAB.hdf5" );
      // 	  }
      // 	writeLevelname( &err, "error0.hdf5" );
      // }

      // SR2: smooth, residual, restrict u, form A( R(u) ) + R( f - A(u) )
      m_op[a_ilev]->SR2_kernel( *a_phi[a_ilev], 
				*m_residual[a_ilev], // rhs
				*m_crsCover[sr_idx], // out: [R(u_f), A_H(R(u_f)) + R(f-A_h(u_h))]
				m_pre,
				m_refRatio,
				m_numGhostCells,
				getNumBufferCells(sr_idx+1) 
				);

      // save R(u)
      LevelData<FArrayBox> R_u_f;
      copyUFromCover(R_u_f,m_crsCover[sr_idx]);
      
      // Vcycle, now in SR V cycle
      SRVCycle( a_phi, a_rhs, a_ilev-1 );
      
      // set error -- debug
      // if( pldf_exact_glob ) {
      // 	LevelData<FArrayBox> err;
      // 	err.define( a_phi[a_ilev]->disjointBoxLayout(), 
      // 		    1,(getNumBufferCells(sr_idx+1)+m_numGhostCells)*IntVect::Unit);
      // 	Real ldx = m_op[a_ilev]->dx();
      // 	for (DataIterator dit = a_phi[a_ilev]->dataIterator(); dit.ok(); ++dit)
      // 	  {
      // 	    const FArrayBox &exact = (*pldf_exact_glob)[dit];
      // 	    FArrayBox &error = err[dit];
      // 	    const FArrayBox &phi = (*a_phi[a_ilev])[dit];
      // 	    Box region = exact.box(); // ghosted
      // 	    for (BoxIterator bit(region); bit.ok(); ++bit)
      // 	      {
      // 		const RealVect offset = bit(); //-domain.smallEnd();
      // 		const RealVect x = ldx*(0.5+offset);
      // 		error(bit()) = exact(bit()) - phi(bit());
      // 	      }
      // 	    writeFABname( &error, "error1FAB.hdf5" );
      // 	  }
      // 	writeLevelname( &err, "error1.hdf5" );
      // }

      // SR3: post leg
      m_op[a_ilev]->SR3_kernel( *a_phi[a_ilev], 
				*m_residual[a_ilev],
				R_u_f,
				*m_crsCover[sr_idx], // out
				m_post,
				m_refRatio,
				m_numGhostCells 
				);

    }
  else // have SR & not fine grid
    {
      CH_assert(0);
      // deep in SR Vcycle, rhis should not happen
      SRVCycle( a_phi, a_rhs, a_ilev );
    }

  return 1.0; // norm 
}

////////////////////////////////////////////////////////////////////////
// AMRFASSR::SRVCycle
//   'm_residual' is really the full RHS ala FAS for coarse grid levels.
////////////////////////////////////////////////////////////////////////
void AMRFASSR::SRVCycle( Vector<RefCountedPtr<LevelData<FArrayBox> > >& a_phi,
			 const Vector<RefCountedPtr<LevelData<FArrayBox> > >& a_rhs,
			 int a_ilev 
			 )
{
  const int l_max = a_phi.size() - 1, nc = a_phi[0]->nComp();
  const int sr_idx = a_ilev - l_max + m_numSRLevels; CH_assert(sr_idx >= 0);
  Interval phiInterval = a_phi[0]->interval();
  Interval resInterval = Interval( nc, 2*nc - 1 );

  // copy-out, this always follows a SR kernel, so get data out

  m_crsCover[sr_idx]->copyTo( phiInterval, 
			      *a_phi[a_ilev], 
			      phiInterval, 
			      m_outCopier[sr_idx] 
			      );
  writeLevelname( a_phi[a_ilev], "copyoutphi.hdf5" );
  // copy out to m_residual, this is so AMRFAS coarse grid solve is correct
  m_crsCover[sr_idx]->copyTo( resInterval, 
			      *m_residual[a_ilev], // rhs set here
			      phiInterval, 
			      m_outCopier[sr_idx]
			      );
writeLevelname( m_residual[a_ilev], "copyoutres.hdf5" );
  if( sr_idx == 0 )
    {
      // coarse grid solve
      AMRFAS<LevelData<FArrayBox> >::VCycle( a_phi, a_rhs, a_ilev, a_ilev ); // no C-F here
    }
  else
    {
      CH_assert(0);
    }
  CH_assert(a_ilev != l_max);
  
  // copy-in, an SR kernel always follows this, so get ready for it
  a_phi[a_ilev]->copyTo( phiInterval, 
			 *m_crsCover[sr_idx], 
			 phiInterval, 
			 m_inCopier[sr_idx] 
			 );
}
