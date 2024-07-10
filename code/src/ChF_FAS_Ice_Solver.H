#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _ChF_FAS_Ice_SOLVER_H_
#define _ChF_FAS_Ice_SOLVER_H_

#include "IceVelocitySolver.H"
#include "ViscousTensorOp.H"
#include "Base_Simplified_AMRLevelOp.H"
#include "AMRFASMultiGrid.H"
///from Chombo/lib/src/AMRElliptic
#include "TensorCFInterp.H"
#include "ViscousTensorOpF_F.H"
#include "AverageF_F.H"
#include "InterpF_F.H"
#include <memory>  //good title for a book
#include "NamespaceHeader.H"

///
/**
   This is a pass-through class for ViscousTensorOp for the nonlinear ice problem.
   Original idea: Each operator owns a ViscousTensorOp.   It turns out this is a terrible idea.
   VTO precomputes tons of stuff.  
   I will poach heavily from  VTO.   I am also using the simpler Base_Simplified_AMRLevelOp
**/
class ChF_Ice_Viscous_Tensor_Op:  public Base_Simplified_AMRLevelOp< LevelData<FArrayBox> >
{
public:
  typedef                     ChF_Ice_Viscous_Tensor_Op  ch_ice_op;
  typedef                                     FArrayBox  ch_fab_cell;
  typedef                                       FluxBox  ch_fab_flux;
  typedef                             DisjointBoxLayout  ch_dbl;
  typedef                                 ProblemDomain  ch_dom;
  typedef                                  DataIterator  ch_dit;
  typedef                                        MayDay  ch_mayday;
  typedef                                           Box  ch_box;
  typedef                             LevelFluxRegister  ch_flux_reg;
  typedef                      LevelData< ch_fab_flux >  ch_ldf_flux;
  typedef                      LevelData< ch_fab_cell >  ch_ldf_cell;
  typedef     Base_Simplified_AMRLevelOp< ch_ldf_cell >  ch_amrlevelop;
  typedef                AMRFASMultiGrid< ch_ldf_cell >  ch_amrfasmg;
  typedef                      MGLevelOp< ch_ldf_cell >  ch_mglevelop;
  typedef                   LinearSolver< ch_ldf_cell >  ch_linear_solver;
  typedef                    RelaxSolver< ch_ldf_cell >  ch_relax_solver;
  typedef                 BiCGStabSolver< ch_ldf_cell >  ch_bicgstab_solver;
  typedef              AMRLevelOpFactory< ch_ldf_cell >  ch_amrlevelop_fact;

private:
  ///strong construction only
  ChF_Ice_Viscous_Tensor_Op();

  ///defining convenience data for this class
  void
  defineLocalData()
  {
    m_dom = m_coarseDomain;
    m_dx = m_dxCrse[0];
    for(int ilev = 0; ilev < m_start_level; ilev++)
    {
      m_dom.refine(m_ref_rat[ilev]);
      m_dx /= m_ref_rat[ilev];
    }
    m_dbl = m_grids[m_start_level];
    if(m_ref_ratioToAMR  > 1)
    {
      m_dom.coarsen(m_ref_ratioToAMR);
      m_dbl.coarsen(m_ref_ratioToAMR);
      m_dx *= m_ref_ratioToAMR;
    }

    m_grad.define(m_dbl, SpaceDim*SpaceDim, IntVect::Unit);
      
      
  }
        

  ///This is  the stuff that came in during construction
  ProblemDomain                m_coarseDomain       ;
  ConstitutiveRelation*        m_constRel           ;
  BasalFrictionRelation*       m_basalFrictionRel   ;
  Vector<DisjointBoxLayout>    m_grids              ;
  Vector<int>                  m_ref_rat             ;
  RealVect                     m_dxCrse             ;
  IceThicknessIBC*             m_bc                 ;
  int                          m_numLevels          ;
  int                          m_start_level         ;
  int                          m_ref_ratioToAMR      ;
  ///this is the local member data (defined at construction)
  /**
     I am trying to keep this at a minimum but more may creep in here for performance reasons.
  **/
  ch_dom                       m_dom;
  ch_dbl                       m_dbl;
  ch_ldf_cell                  m_grad;
  Real                         m_dx;
  LevelDataOps<FArrayBox>      m_levelOps;

      
public:
  ///destructor  is trivial but has to be virtual
  virtual ~ChF_Ice_Viscous_Tensor_Op(){;}

  ///
  /**
     compute coefficients needed for the operator given an input state and all that other stuff we saved.
  **/
  void  computeCoefficients( RefCountedPtr< ch_ldf_flux >       &      a_eta,
                             RefCountedPtr< ch_ldf_flux >       &      a_lambda,
                             RefCountedPtr< ch_ldf_cell >       &      a_acoef,
                             RefCountedPtr< ch_ldf_cell >       &      a_relco,
                             Real                               &      a_alpha,
                             Real                               &      a_beta,
                             const ch_ldf_cell                  &      a_phi)

  {
    ch_mayday::Error("ChF_Ice_VTO::computeCoefficients not implemented yet.");
  }
      
  
  /// 
  ChF_Ice_Viscous_Tensor_Op(const ch_dom                   &  a_coarseDomain       ,
                            ConstitutiveRelation           *  a_constRel           ,
                            BasalFrictionRelation          *  a_basalFrictionRel   ,
                            const Vector<DisjointBoxLayout>&  a_grids          ,
                            const Vector<int>              &  a_ref_rat       ,
                            const RealVect                 &  a_dxCrse             ,
                            IceThicknessIBC                *  a_bc                 ,
                            int                               a_numLevels          ,
                            int                               a_start_level         ,
                            int                               a_ref_ratioToAMR      )
  {

    m_coarseDomain       = a_coarseDomain     ;
    m_constRel           = a_constRel         ;
    m_basalFrictionRel   = a_basalFrictionRel ;
    m_grids              = a_grids            ;
    m_ref_rat             = a_ref_rat           ;
    m_dxCrse             = a_dxCrse           ;
    m_bc                 = a_bc               ;
    m_numLevels          = a_numLevels        ;
    m_start_level         = a_start_level       ;
    m_ref_ratioToAMR      = a_ref_ratioToAMR    ;
    defineLocalData();
  }

  ///mostly poached from lib/src/AMRElliptic/ViscousTensorOp.cpp
  virtual void
  prolongIncrement(ch_ldf_cell      & a_phiThisLevel,
                   const ch_ldf_cell& a_correctCoarse)
  {

    CH_TIME("ChF_ICE_VTO::prolongIncrement");
      
    DisjointBoxLayout dbl = a_phiThisLevel.disjointBoxLayout();
    int mgref = 2; //this is a multigrid func

    for (DataIterator dit = a_phiThisLevel.dataIterator(); dit.ok(); ++dit)
    {
      FArrayBox& phi =  a_phiThisLevel[dit];
      const FArrayBox& coarse = a_correctCoarse[dit];
      Box region = dbl.get(dit());
      Box cBox = coarsen(region, mgref);
      /**
         Standard VTO has a linear interpolation option
         This code matches the case in VTO where
         s_prolongType == piecewiseConstant
      **/
      FORT_PROLONGVTOP(CHF_FRA(phi),
                       CHF_CONST_FRA(coarse),
                       CHF_BOX(region),
                       CHF_CONST_INT(mgref));
    }
  } //end function prolongIncrement
    
  ///
  virtual ch_dom
  getDomainFromData(const ch_ldf_cell & a_data)
  {
    return a_data.disjointBoxLayout().physDomain();
  }
  ///
  /**
     I am indulging in the fantasy that I can finally write a bc function that makes more sense.  
     This fills corners outside the domain with something reasonable.    
     As the Boss said, love is bigger than a Honda.  It really is bigger than a Subaru.
  **/
  void reflectBCSide(ch_fab_cell      &      a_state,
                     const ch_box&      a_valid,
                     const ch_dom&      a_dom,
                     int                a_dir,
                     int                a_side,
                     int                a_comp,
                     int                a_nghost,
                     bool a_reflectOdd)
  {

    CH_TIME("CHF_ICE_VTO::reflectBCSide");
    ///only do stuff if we have ghost cells
    if(!a_dom.domainBox().contains(a_state.box()))
    {
      ///
      int isign = 1;
      if(a_side == 0)
      {
        isign = -1;
      }
      ch_box offdir_grown_dom_box = a_dom.domainBox();
      ch_box offdir_grown_valid_b = a_valid;
      for(int idir = 0; idir < SpaceDim; idir++)
      {
        if(idir != a_dir)
        {
          offdir_grown_dom_box.grow(idir, a_nghost);
          offdir_grown_valid_b.grow(idir, a_nghost);
        }
      }
      //only grow by one.  need the starting cell for a row of ghost
      ch_box start_box = offdir_grown_valid_b;
      start_box.grow(1);
      start_box.shift(a_dir, isign);
      start_box &= offdir_grown_dom_box;
      start_box.shift(a_dir, -isign);
      //this should cover (with something reasonable) corners outside the domain
      for (BoxIterator bit(start_box); bit.ok(); ++bit)
      {
        IntVect iv_first_ghost = bit();
        for(int ighost = 0; ighost < a_nghost; ighost++)
        {
          IntVect iv_src =  iv_first_ghost - (ighost+1)*isign*BASISV(a_dir);
          IntVect iv_dst =  iv_first_ghost + ( ighost )*isign*BASISV(a_dir);
          Real reflectSign  = 1.;
          if(a_reflectOdd)
          {
            reflectSign = -1;
          }
          a_state(iv_dst, a_comp) = reflectSign*a_state(iv_src, a_comp);
        }
      }
    }
  } // end function homogDiriSide


    ///
    /**
       From Base_Simplified_AMRLevelOP
    **/
  virtual void
  fillGhostDataOutsideDomain(ch_ldf_cell          & a_lhs,
                             const ch_dom    & a_domain) 
  {
    CH_TIME("CHF_ICE_VTO::fillGhostDataOutsideDomain");
    ch_dbl dbl = a_lhs.disjointBoxLayout();
    ch_dit dit = a_lhs.dataIterator();
    int nghost = a_lhs.ghostVect()[0];
    for(int ibox = 0; ibox < dit.size(); ibox++)
    {
      const ch_box& valid =   dbl[dit[ibox]];
      ch_fab_cell      & state = a_lhs[dit[ibox]];
      for(int iside = 0; iside < 1; iside++)
      {
        for(int idir = 0; idir < SpaceDim; idir++)
        {
          for(int icomp = 0; icomp < SpaceDim; icomp++)
          {
            //ye olde slipwall boundary condition
            bool is_diri = (idir == icomp);
            reflectBCSide(state,
                          valid,
                          m_dom,
                          idir,
                          iside,
                          icomp,
                          nghost,
                          is_diri);
          }
        }
      }
    }
  }

  ///
  virtual void
  fillGhostDataOverFineFineInterface(ch_ldf_cell &  a_phi) 
  {
    CH_TIME("CHF_ICE_VTO::fillDataOverFineFineInterface");
    a_phi.exchange();
  }

  ///
  virtual void
  copyTo(ch_ldf_cell       &  a_dst,
         const ch_ldf_cell &  a_src)
  {
    CH_TIME("CHF_ICE_VTO::copyTo");
    a_src.copyTo(a_src.interval(), a_dst, a_dst.interval());
  }
    

  ///if phi_coar == NULL, zeros on coarse-fine interface
  virtual void
  fillGhostDataOverCoarseFineInterface(ch_ldf_cell       &  a_phi_fine, 
                                       const ch_ldf_cell *  a_phi_coar_ptr)
  {
    CH_TIME("CHF_ICE_VTO::fillGhostDataOverCoarsefineInterface");
    if(m_start_level > 0)
    {
      Real    dx_fine = m_dx;
      Box     do_fine = m_dom.domainBox();
      int  refToCoar = m_ref_rat[m_start_level-1];
      setToZero(m_grad);
      int ncomp = SpaceDim;
      DisjointBoxLayout  grids_fine     = a_phi_fine.disjointBoxLayout();
      DisjointBoxLayout* grids_coar_ptr = &m_grids[m_start_level-1];
      TensorCFInterp    interp_with_coarser;  //_
      interp_with_coarser.define(grids_fine, grids_coar_ptr,
                                 dx_fine, refToCoar, ncomp,
                                 do_fine);
      if(a_phi_coar_ptr != NULL)
      {
        interp_with_coarser.coarseFineInterp(a_phi_fine, m_grad, *a_phi_coar_ptr);
      }
      else
      {
        interp_with_coarser.coarseFineInterpH(a_phi_fine, m_grad);
      }
    }//end if(we are on amrlevel > 0
  }

  ///mostly poached from lib/src/AMRElliptic/ViscousTensorOp.cpp
  void
  cellGrad(ch_fab_cell        &      a_gradPhi,
           const ch_fab_cell  &      a_phi,
           const ch_box  &      a_grid)
  {
    CH_assert(a_gradPhi.nComp() == SpaceDim*SpaceDim);
    CH_assert(a_phi.nComp() == SpaceDim);

    for (int derivDir = 0; derivDir < SpaceDim; derivDir++)
    {
      for (int phiDir = 0; phiDir < SpaceDim; phiDir++)
      {
        int gradcomp = TensorCFInterp::gradIndex(phiDir,derivDir);
        FORT_CELLGRADVTOP(CHF_FRA1(a_gradPhi, gradcomp),
                          CHF_CONST_FRA1(a_phi, phiDir),
                          CHF_BOX(a_grid),
                          CHF_CONST_REAL(m_dx),
                          CHF_CONST_INT(derivDir));

      }
    }
  }
  ///
  void
  fillGrad(const ch_ldf_cell&       a_phi)
  {
    CH_TIME("ViscousTensorOp::fillGrad");  
    LevelData<FArrayBox>& phi = (LevelData<FArrayBox>&)a_phi;
    phi.exchange();
    const DisjointBoxLayout& grids = a_phi.disjointBoxLayout();
    //compute gradient of phi for parts NOT in ghost cells
    for (DataIterator dit = grids.dataIterator(); dit.ok(); ++dit)
    {
      cellGrad(m_grad[dit()],
               a_phi [dit()],
               grids.get(dit()));
    }
    m_grad.exchange();
  }
  ///poached from Chombo/lib/src/AMRElliptic/ViscousTensorOp
  void
  loHiCenterFace(Box&                 a_loBox,
                 int&                 a_hasLo,
                 Box&                 a_hiBox,
                 int&                 a_hasHi,
                 Box&                 a_centerBox,
                 const ProblemDomain& a_eblg,
                 const Box&           a_inBox,
                 const int&           a_dir)
  {

    // if we're periodic in a_dir, then all boxes are
    // contained in the domain (DFM 2/5/10)
    if (a_eblg.isPeriodic(a_dir) )
    {
      a_centerBox = a_inBox;
      a_hasLo = 0;
      a_hasHi = 0;
    } // end if periodic in a_dir
    else
    {
      Box domainFaceBox = a_eblg.domainBox();
      domainFaceBox.surroundingNodes(a_dir);

      // Make a copy of the input box which can be modified
      Box inBox = a_inBox;

      inBox &= domainFaceBox;

      a_centerBox = inBox;
      a_centerBox.grow(a_dir, 1);
      a_centerBox &= domainFaceBox;
      a_centerBox.grow(a_dir,-1);

      // See if this chops off the high side of the input box
      Box tmp = inBox;
      tmp.shift(a_dir,1);
      tmp &= domainFaceBox;
      tmp.shift(a_dir,-1);

      // If so, set up the high, one-sided difference box, a_hiBox
      if (tmp != inBox)
      {
        a_hasHi = 1;
        a_hiBox = bdryHi(inBox,a_dir,1);
      }
      else
      {
        a_hasHi = 0;
      }

      // See if this chops off the low side of the input box
      tmp = inBox;
      tmp.shift(a_dir,-1);
      tmp &= domainFaceBox;
      tmp.shift(a_dir,1);

      // If so, set up the low, one-sided difference box, a_loBox
      if (tmp != inBox)
      {
        a_hasLo = 1;
        a_loBox = bdryLo(inBox,a_dir,1);
      }
      else
      {
        a_hasLo = 0;
      }
    } // end if not periodic in a_dir
  } //end function loHiCenterFace

  /// another refugee from AMRElliptic VTO
  void
  getFaceDivAndGrad(FArrayBox&             a_faceDiv,
                    FArrayBox&             a_faceGrad,
                    const FArrayBox&       a_data,
                    const FArrayBox&       a_gradData,
                    const ProblemDomain&   a_domain,
                    const Box&             a_faceBox,
                    const int&             a_faceDir,
                    const Real             a_dx)
  {
    //set divergence to zero everywhere
    a_faceDiv.setVal(0.);

    Box interiorFaceBox ;
    Box loBoxFace, hiBoxFace;
    int hasLo, hasHi;

    loHiCenterFace(loBoxFace, hasLo,
                   hiBoxFace, hasHi,
                   interiorFaceBox,
                   a_domain,
                   a_faceBox,
                   a_faceDir);

    for (int divDir = 0; divDir < SpaceDim; divDir++)
    {
      int gradComp = TensorCFInterp::gradIndex(divDir,divDir);
      FORT_FACEDIVINCRVTOP(CHF_FRA1(a_faceDiv, 0),
                           CHF_FRA(a_data),
                           CHF_FRA(a_gradData),
                           CHF_BOX(a_faceBox),
                           CHF_BOX(interiorFaceBox),
                           CHF_BOX(loBoxFace),
                           CHF_INT(hasLo),
                           CHF_BOX(hiBoxFace),
                           CHF_INT(hasHi),
                           CHF_REAL(a_dx),
                           CHF_INT(a_faceDir),
                           CHF_INT(divDir),
                           CHF_INT(gradComp));

      //now average cell-centered gradients to the face centers
      //use diffs in data if divDir == faceDir
    }

    for (int divDir = 0; divDir < SpaceDim; divDir++)
    {
      for (int velDir = 0; velDir < SpaceDim; velDir++)
      {
        int gradcomp = TensorCFInterp::gradIndex(velDir,divDir);
        FORT_GETFACEGRADVTOP(CHF_FRA1(a_faceGrad, gradcomp),
                             CHF_FRA1(a_gradData, gradcomp),
                             CHF_FRA1(a_data, velDir),
                             CHF_BOX(a_faceBox),
                             CHF_BOX(interiorFaceBox),
                             CHF_BOX(loBoxFace),
                             CHF_INT(hasLo),
                             CHF_BOX(hiBoxFace),
                             CHF_INT(hasHi),
                             CHF_REAL(a_dx),
                             CHF_INT(a_faceDir),
                             CHF_INT(divDir));

      }

    }
  } //end function getfacedivandgrad

  /**
     adapted from Chombo/lib/src/AMRElliptic/VTO
  **/
  void
  getFlux(FArrayBox&       a_flux,
          const FArrayBox& a_data,
          const FArrayBox& a_gradData,
          const FArrayBox& a_etaFace,
          const FArrayBox& a_lambdaFace,
          const Box&       a_faceBox,
          int              a_dir,
          Real             a_beta)
  {
    CH_TIME("ViscousTensorOp::getflux");  

    CH_assert(a_data.nComp() == SpaceDim);
    a_flux.resize(a_faceBox, a_data.nComp());
    FArrayBox  faceDiv(a_faceBox, 1);
    FArrayBox  faceGrad(a_faceBox, a_gradData.nComp());
    getFaceDivAndGrad(faceDiv, faceGrad, a_data, a_gradData, m_dom, a_faceBox, a_dir, m_dx);

    getFluxFromDivAndGrad(a_flux, faceDiv, faceGrad, a_etaFace, a_lambdaFace, a_faceBox, a_dir);
    a_flux *= a_beta;
  }
  
  /**
     From Ye Olde VTO.
  **/
  void
  getFluxFromDivAndGrad(FArrayBox&       a_flux,
                        const FArrayBox& a_faceDiv,
                        const FArrayBox& a_faceGrad,
                        const FArrayBox& a_etaFace,
                        const FArrayBox& a_lambdaFace,
                        const Box&       a_faceBox,
                        int              a_dir)

  {
    CH_TIME("ViscousTensorOp::getFluxFromDivAndGrad");  
    //copy the divergence into the diagonal component of the flux
    //diagonal component of the flux is the a_dir comp
    int dstcomp = a_dir;
    int srccomp = 0;
    int numcomp = 1;

    //this is important because we increment later
    a_flux.setVal(0.);

    //now make diag component of f = div v
    a_flux.copy(a_faceDiv, srccomp, dstcomp, numcomp);
    //the flux now contains div v * I

    //now flux gets multiplied by lambda
    a_flux.mult(a_lambdaFace, srccomp, dstcomp, 1);

    //the flux now contains div v * lambda * I
    //now add in eta*(grad v + (grad v)^T)
    for (int idir = 0; idir < SpaceDim; idir++)
    {
      int fluxComp = idir;
      int gradComp = TensorCFInterp::gradIndex(idir,a_dir);
      int tranComp = TensorCFInterp::gradIndex(a_dir, idir);
      FORT_ADDGRADTOFLUXVTOP(CHF_FRA(a_flux),
                             CHF_CONST_FRA1(a_etaFace,0),
                             CHF_CONST_INT(fluxComp),
                             CHF_CONST_FRA(a_faceGrad),
                             CHF_CONST_INT(gradComp),
                             CHF_CONST_INT(tranComp),
                             CHF_BOX(a_faceBox));
    }
  } //end function get flux from div and grad

  /// called by the coarser level
  virtual void
  refluxFineCorrection(ch_ldf_cell       &   a_lph_coar,
                       const ch_ldf_cell &   a_phi_fine, 
                       const ch_ldf_cell &   a_phi_coar,
                       AMRLevelOp<ch_ldf_cell>*   a_finerBase)
  {
    CH_TIME("CHF_ICE_VTO::refluxFineCorrection");
    Interval interv(0, SpaceDim-1);
    int ref_to_fine = m_ref_rat[m_start_level];
    ch_dbl dbl_fine = a_phi_fine.disjointBoxLayout();
    ch_dbl dbl_coar = a_phi_coar.disjointBoxLayout();
    ChF_Ice_Viscous_Tensor_Op*
      op_f = dynamic_cast<ChF_Ice_Viscous_Tensor_Op*>(a_finerBase);
    /**
       VTO skips this step (because of context, I guess).
       Ergo, the next two lines  may be unnecessary.
    **/
    this->fillGrad(a_phi_coar);
    ch_ldf_cell& phi_fine_cast = (ch_ldf_cell&)a_phi_fine;
    op_f->fillGhostDataOverCoarseFineInterface(phi_fine_cast, &a_phi_coar);
    op_f->fillGrad(a_phi_fine);
    ch_dom dom_fine = refine(m_dom, ref_to_fine);
      
    ch_flux_reg flux_reg(dbl_fine, dbl_coar, dom_fine, ref_to_fine, SpaceDim);
    flux_reg.setToZero();
    RefCountedPtr< ch_ldf_flux >          eta_fine,   eta_coar;
    RefCountedPtr< ch_ldf_flux >       lambda_fine,lambda_coar;
    RefCountedPtr< ch_ldf_cell >        acoef_fine, acoef_coar;
    RefCountedPtr< ch_ldf_cell >        relco_fine, relco_coar;
    Real                                alpha_fine, alpha_coar;
    Real                                 beta_fine,  beta_coar;
    this->computeCoefficients( eta_coar, lambda_coar, acoef_coar, relco_coar, alpha_coar, beta_coar, a_phi_coar);
    op_f->computeCoefficients( eta_fine, lambda_fine, acoef_fine, relco_fine, alpha_fine, beta_fine, a_phi_fine);
    CH_assert(std::abs(alpha_coar - alpha_fine) < 1.0e-10);
    CH_assert(std::abs( beta_coar -  beta_fine) < 1.0e-10);
    Real alpha = alpha_coar; Real beta = beta_coar;


    {//coarse bit
      CH_TIME("incrementCoarse");
      DataIterator dit = a_phi_coar.dataIterator();
      for (dit.reset(); dit.ok(); ++dit)
      {
        const FArrayBox& coarfab    =    a_phi_coar[dit];
        const FluxBox&   coareta    =    (*eta_coar)[dit];
        const FluxBox&   coarlambda = (*lambda_coar)[dit];
        const Box&       gridBox    =   a_phi_coar.getBoxes()[dit];

        for (int idir = 0; idir < SpaceDim; idir++)
        {
          FArrayBox coarflux;
          Box faceBox = surroundingNodes(gridBox, idir);
          //fillgrad was called before this when applyOp was called
          this->getFlux(coarflux, coarfab, m_grad[dit()], coareta[idir], coarlambda[idir], faceBox, idir, beta);

          Real scale = 1.0;
          flux_reg.incrementCoarse(coarflux, scale,dit(),
                                   interv,interv,idir);
        }
      }
    }//end coarse bit
    {//fine bit
      
      CH_TIME("incrementFine");
      DataIterator dit = a_phi_fine.dataIterator();
      for (dit.reset(); dit.ok(); ++dit)
      {
        const FArrayBox& finefab    =    a_phi_fine[dit];
        const FluxBox&   fineeta    =    (*eta_fine)[dit];
        const FluxBox&   finelambda = (*lambda_fine)[dit];
        const Box&       gridBox    =   a_phi_fine.getBoxes()[dit];

        for (int idir = 0; idir < SpaceDim; idir++)
        {
          FArrayBox fineflux;
          Box faceBox = surroundingNodes(gridBox, idir);
          //fillgrad was called before this when applyOp was called
          op_f->getFlux(fineflux, finefab, op_f->m_grad[dit()], fineeta[idir], finelambda[idir], faceBox, idir, beta);

          Real scale = 1.0;
          flux_reg.incrementFine(fineflux, scale,dit(),
                                 interv,interv,idir);
        }
      }
    } //end fine bit
    {
      ///and now for the big finish
      CH_TIME("reflux");
      Real scale =  beta/this->m_dx;
      flux_reg.reflux(a_lph_coar, scale);
    }
  }

  ///For this one, the derived class may assume all ghost cells have been set correctly
  virtual void 
  applyOperator(ch_ldf_cell        & a_lph,
                const ch_ldf_cell  & a_phi)
  {
    CH_TIME("CHF_ICE_VTO::applyOperator");
    RefCountedPtr< ch_ldf_flux >     eta_ptr, lambda_ptr;
    RefCountedPtr< ch_ldf_cell >   acoef_ptr,  relco_ptr;
    Real                          alpha, beta;
    ch_ldf_cell& phi_cast = (ch_ldf_cell &)(a_phi);
    this->computeCoefficients( eta_ptr, lambda_ptr, acoef_ptr, relco_ptr, alpha, beta, a_phi);

    /**
       already done 
       fillGhostDataOutsideDomain(        phi_cast);
       fillGhostDataOverFineFineInterface(phi_cast);
    **/
    fillGrad(a_phi);
    const DisjointBoxLayout dbl= a_phi.getBoxes();
    DataIterator dit = a_phi.dataIterator();
    //this makes the lhs = alpha*phi
    m_levelOps.setToZero(a_lph);
    incr(a_lph, a_phi, alpha);
    int sdim_comp = SpaceDim;
    for (dit.begin(); dit.ok(); ++dit)
    {
      for (int ivar = 0; ivar < SpaceDim; ivar++)
      {
        int src = 0; int dst = ivar; int one_comp = 1;
        a_lph[dit()].mult((*acoef_ptr)[dit()],src,dst,one_comp);
      }
      Box gridBox = dbl.get(dit());
      FluxBox flux(gridBox, SpaceDim);
      FluxBox& eta    = (   *eta_ptr)[dit];
      FluxBox& lambda = (*lambda_ptr)[dit];

      //operator is alpha I + (divF) = alpha I + beta div(eta(grad B + grad B^T) + lambda I div B )
      for (int idir = 0; idir < SpaceDim; idir++)
      {
        Box faceBox = surroundingNodes(gridBox, idir);
        getFlux(flux[idir], a_phi[dit()], m_grad[dit()], eta[idir], lambda[idir],faceBox, idir, beta);
        FArrayBox& lhsFAB = a_lph[dit()];
        const FArrayBox& fluxFAB = flux[idir];

        //beta and eta are part of the flux.
        FORT_ADDDIVFLUXDIRVTOP(CHF_FRA(lhsFAB),
                               CHF_CONST_FRA(fluxFAB),
                               CHF_BOX(gridBox),
                               CHF_REAL(m_dx),
                               CHF_CONST_INT(sdim_comp),
                               CHF_CONST_INT(idir));
      }
    } // end loop over boxes

  }

/// Restrict fine data to coarse data using volume-weighted averaging.   
  virtual void
  restriction(ch_ldf_cell & a_resCoar,
              ch_ldf_cell & a_resFine,
              int a_ref_rat)
  {
    CH_TIME("CHF_ICE_VTO::restriction");
    Box ref_box(IntVect::Zero, (a_ref_rat-1)*IntVect::Unit); //_
    DisjointBoxLayout  dbl_coar = a_resCoar.disjointBoxLayout();
    DataIterator dit = dbl_coar.dataIterator();
    for (dit.begin(); dit.ok(); ++dit)
    {
      FArrayBox&       res_coar = a_resCoar[dit()];
      const FArrayBox& res_fine = a_resFine[dit()];
      Box valid_coar            =  dbl_coar[dit()];
      FORT_AVERAGE( CHF_FRA(        res_coar),
                    CHF_CONST_FRA(  res_fine),
                    CHF_BOX(      valid_coar),
                    CHF_CONST_INT(a_ref_rat),
                    CHF_BOX(        ref_box));
    } 

  }// end function restriction

///from LinearOp
  virtual void
  create(ch_ldf_cell      & a_lhs,
         const ch_ldf_cell& a_rhs)
  {
    CH_TIME("CHF_ICE_VTO::create");
    m_levelOps.create(a_lhs, a_rhs);
  }

///from LinearOp
  virtual void
  assign(ch_ldf_cell      & a_lhs,
         const ch_ldf_cell& a_rhs) 
  {
    CH_TIME("CHF_ICE_VTO::assign");
    m_levelOps.assign(a_lhs, a_rhs);
  }
    
///from LinearOp
  virtual Real
  dotProduct(const ch_ldf_cell      & a_1,
             const ch_ldf_cell      & a_2)     
  {
    CH_TIME("CHF_ICE_VTO::dotProduct");
    Real retval = 4586;
    retval = m_levelOps.dotProduct(a_1, a_2);
    return retval;
  }


///from LinearOp
  virtual void
  axby(ch_ldf_cell      & a_lhs,
       const ch_ldf_cell& a_x,
       const ch_ldf_cell& a_y,
       Real          a_a,
       Real          a_b)
  {
    CH_TIME("CHF_ICE_VTO::axby");
    m_levelOps.axby(a_lhs, a_x, a_y, a_a, a_b);
  }

///from LinearOp
  virtual Real
  norm(const ch_ldf_cell & a_rhs, int a_ord)
  {
    CH_TIME("CHF_ICE_VTO::norm");
    Real retval = 4586;
    retval = CH_XD::norm(a_rhs, a_rhs.interval(), a_ord);

    return retval;
  }
    
///from LinearOp
  virtual void
  setToZero(ch_ldf_cell & a_lhs)
  {
    CH_TIME("CHF_ICE_VTO::setToZero");
    m_levelOps.setToZero(a_lhs);
  }

///from MGLevelOp
  virtual void
  createCoarser(ch_ldf_cell      & a_coarse,
                const ch_ldf_cell& a_fine,
                bool ghosted)
  {
    CH_TIME("CHF_ICE_VTO::createCoarser");

    IntVect ghost = a_fine.ghostVect();
    DisjointBoxLayout dbl;
    CH_assert(dbl.coarsenable(2));
    coarsen(dbl, a_fine.disjointBoxLayout(), 2); //multigrid, so coarsen by 2
    a_coarse.define(dbl, a_fine.nComp(), ghost);
  }
    
  void
  getMultiColors(Vector<IntVect>& a_colors) const
  {

#if CH_SPACEDIM==1
    a_colors.resize(2);
    a_colors[0] = IntVect::Zero;//(0,0)
    a_colors[1] = IntVect::Unit;//(1,1)
#elif CH_SPACEDIM==2
    a_colors.resize(4);
    a_colors[0] = IntVect::Zero;//(0,0)
    a_colors[1] = IntVect::Unit;//(1,1)
    a_colors[2] = IntVect::Zero + BASISV(1);//(0,1)
    a_colors[3] = IntVect::Zero + BASISV(0);//(1,0)
#elif CH_SPACEDIM==3
    a_colors.resize(8);
    a_colors[0] = IntVect::Zero;//(0,0,0)
    a_colors[1] = IntVect::Zero + BASISV(0) + BASISV(1);//(1,1,0)
    a_colors[2] = IntVect::Zero + BASISV(1) + BASISV(2);//(0,1,1)
    a_colors[3] = IntVect::Zero + BASISV(0) + BASISV(2);//(1,0,1)
    a_colors[4] = IntVect::Zero + BASISV(1);//(0,1,0)
    a_colors[5] = IntVect::Zero + BASISV(0);//(1,0,0)
    a_colors[6] = IntVect::Zero + BASISV(2);//(0,0,1)
    a_colors[7] = IntVect::Unit;//(1,1,1)
#endif
  }
///from MGLevelOp
  virtual void
  relax(ch_ldf_cell      & a_phi,
        const ch_ldf_cell& a_rhs,
        int a_iterations)
  {
    CH_TIME("CHF_ICE_VTO::relax");
    RefCountedPtr< ch_ldf_flux >    eta_ptr, lambda_ptr;
    RefCountedPtr< ch_ldf_cell >  acoef_ptr, relco_ptr;
    Real                          alpha, beta;
    ch_ldf_cell& phi_cast = (ch_ldf_cell &)(a_phi);
    this->computeCoefficients( eta_ptr, lambda_ptr, acoef_ptr, relco_ptr, alpha, beta, a_phi);

    LevelData<FArrayBox> lphi;
    m_levelOps.create(lphi, a_rhs);
    LevelData<FArrayBox> prevPhi;
    int var_comp = SpaceDim;
    m_levelOps.create(prevPhi, a_rhs);

    // do first red, then black passes
    const DisjointBoxLayout& dbl = a_phi.disjointBoxLayout();
    int whichIter = 0;
    bool done = false;

    Vector<IntVect> relax_colors;
    getMultiColors(relax_colors);
    for(int iiter = 0; iiter < a_iterations; iiter++)
    {
      for (int icolor = 0; icolor < relax_colors.size(); icolor++)
      {
        fillGhostData(a_phi, NULL);
        const IntVect& color= relax_colors[icolor];
        //after this lphi = L(phi)
        applyOperator(lphi, a_phi);

        for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit)
        {
          Box dblBox  = dbl.get(dit());

          IntVect loIV = dblBox.smallEnd();
          IntVect hiIV = dblBox.bigEnd();

          for (int idir = 0; idir < SpaceDim; idir++)
          {
            if (abs(loIV[idir]) % 2 != color[idir])
            {
              loIV[idir]++;
            }
          }

          if (loIV <= hiIV)
          {
            ch_fab_cell       & phi_fab =        a_phi[dit];
            const ch_fab_cell & lph_fab =         lphi[dit];
            const ch_fab_cell & rhs_fab =        a_rhs[dit];
            const ch_fab_cell & lam_fab = (*relco_ptr)[dit];
            Box colored_box(loIV, hiIV);
            FORT_GSRBVTOP(CHF_FRA(      phi_fab),
                          CHF_CONST_FRA(lph_fab),
                          CHF_CONST_FRA(rhs_fab),
                          CHF_CONST_FRA(lam_fab),
                          CHF_BOX(colored_box),
                          CHF_CONST_INT(var_comp));
          } // end weird loiv/hiiv check
        }   // end loop over grids
      }     // end loop through red-black
    }       // end loop over iterations
  }         // end function relax

  ///
  virtual int refToCoarser()
  {
    int retval = 4586;
    if(m_start_level > 0)
    {
      retval = m_ref_rat[m_start_level-1];
    }
    else
    {
      ch_mayday::Error("chf_fas::confusion: refToCoarser called when there is no coarser");
    }

    return retval;
  }
  
  ///borrowed from VTO
  virtual void
  AMRRestrict(ch_ldf_cell      & a_res_coar,
              const ch_ldf_cell& a_res_fine,
              const ch_ldf_cell& a_cor_fine,
              const ch_ldf_cell& a_cor_coar,
              bool a_skip_res)
  {
    CH_TIME("CHF_ICE_VTO::AMRRestrict");

    ch_ldf        resid_fine;
    create(       resid_fine, a_residual);
    ///a_resid_fine used as RHS here
    AMRResidualNF(resid_fine, a_correction, a_coarseCorrection, a_res_fine, true);
    ch_dbl  dbl_coar = a_res_coar.disjointBoxLayout();
    DataIterator dat_iter = a_res_fine.dataIterator();
    ///Think of ref_box as the box within the cell.  
    int  ref_to_coar = m_ref_rat[m_start_level-1];
    Box ref_box(IntVect::Zero, (m_refToCoar-1)*IntVect::Unit);
    for (int ibox = 0; ibox < dat_iter.size(); ibox++)
    {
      auto dit = dat_iter[ibox];
      ch_fab       & res_fab_coar = a_res_coar[dit];
      const ch_fab & res_fab_fine = resid_fine[dit];
      const Box    & valid        =    dblCoar[dit];
      FORT_AVERAGE( CHF_FRA(      res_fab_coar),
                    CHF_CONST_FRA(res_fab_fine),
                    CHF_BOX(valid),
                    CHF_CONST_INT(ref_to_coar),
                    CHF_BOX(ref_box) );
    }

  }

  ///borrowed from VTO
  /** a_correction += I[2h->h](a_coarseCorrection) */
  virtual void
  AMRProlong(ch_ldf_cell      & a_corr_fine,
             const ch_ldf_cell& a_corr_coar)
  {
    CH_TIME("CHF_ICE_VTO::AMRProlong");
    pr_dbl grid_fine = a_corr_fine.disjointBoxLayout();
    DisjointBoxLayout grid_cofi;
    int   ncomp = a_correction.nComp();
    Interval interv(0, ncomp-1);
    IntVect ivg = a_coarseCorrection.ghostVect();
    int  ref_to_coar = m_ref_rat[m_start_level-1];
    coarsen(grid_cofi,  grid_fine, ref_to_coar);
    ch_ldf corr_cofi(grids_cofi, ncomp , ivg);
    a_corr_coar.copyTo(interv, corr_cofi, interv);

    Real dx_coar = m_dx*ref_to_coar;
    DataIterator dit_fine = grid_fine.dataIterator();
    for(int ibox = 0; ibox < dit_fine.size(); ibox++)
    {
      auto dit = dit_fine[ibox];
      ch_ldf& corr_fab_fine =  a_corr_fine[dit];
      ch_ldf& corr_fab_cofi =    corr_cofi[dit];
      Box        valid_fine =   grid_fine[dit];

      // need this if we need to set physical BC's
      // on coarse-level correction.
      ProblemDomain crseDomain(m_domain);
      crseDomain.coarsen(m_refToCoar);

      FORT_PROLONGVTOP(CHF_FRA(phi),
                       CHF_CONST_FRA(coarse),
                       CHF_BOX(valid_fine),
                       CHF_CONST_INT(m_refToCoar));
  }

///
/** a_residual = a_residual - L(a_correction, a_coarseCorrection) */
  virtual void
  AMRUpdateResidual(ch_ldf_cell      & a_residual,
                    const ch_ldf_cell& a_correction,
                    const ch_ldf_cell& a_coarseCorrection)
  {
    CH_TIME("CHF_ICE_VTO::AMRUpdateResidual");
    ch_mayday::Error("not implmented");
  }

///
/**
 */
  virtual void createCoarsened(ch_ldf_cell&       a_lhs,
                               const ch_ldf_cell& a_rhs,
                               const int&     a_ref_rat)
  {
    CH_TIME("CHF_ICE_VTO::createCoarsened");
    ch_mayday::Error("not implmented");
  }
}; //end class ChF_Ice_Viscous_Tensor_Op
  
  
/// class to make ChF_Ice_VTO's.   
class ChF_Ice_VTO_Factory: public AMRLevelOpFactory< LevelData<FArrayBox> >
{
public:
  typedef                     ChF_Ice_Viscous_Tensor_Op  ch_ice_op;
  typedef                                     FArrayBox  ch_fab_cell;
  typedef                                       FluxBox  ch_fab_flux;
  typedef                             DisjointBoxLayout  ch_dbl;
  typedef                                 ProblemDomain  ch_dom;
  typedef                                  DataIterator  ch_dit;
  typedef                                        MayDay  ch_mayday;
  typedef                                           Box  ch_box;
  typedef                             LevelFluxRegister  ch_flux_reg;
  typedef                      LevelData< ch_fab_flux >  ch_ldf_flux;
  typedef                      LevelData< ch_fab_cell >  ch_ldf_cell;
  typedef     Base_Simplified_AMRLevelOp< ch_ldf_cell >  ch_amrlevelop;
  typedef                AMRFASMultiGrid< ch_ldf_cell >  ch_amrfasmg;
  typedef                      MGLevelOp< ch_ldf_cell >  ch_mglevelop;
  typedef                   LinearSolver< ch_ldf_cell >  ch_linear_solver;
  typedef                    RelaxSolver< ch_ldf_cell >  ch_relax_solver;
  typedef                 BiCGStabSolver< ch_ldf_cell >  ch_bicgstab_solver;
  typedef              AMRLevelOpFactory< ch_ldf_cell >  ch_amrlevelop_fact;
  ///
  /**
     The only member data for this class is the stuff that comes in through the constructor.   
  **/
  ProblemDomain              m_coarseDomain       ;
  ConstitutiveRelation*      m_constRel           ;
  BasalFrictionRelation*     m_basalFrictionRel   ;
  Vector<DisjointBoxLayout>  m_grids          ;
  Vector<int>                m_ref_rat       ;
  RealVect                   m_dxCrse             ;
  IceThicknessIBC*           m_bc                 ;
  int                        m_numLevels          ;


  ///
  ChF_Ice_VTO_Factory(const ProblemDomain&              a_coarseDomain       ,
                      ConstitutiveRelation*             a_constRel           ,
                      BasalFrictionRelation*            a_basalFrictionRel   ,
                      const Vector<DisjointBoxLayout>&  a_vectGrids          ,
                      const Vector<int>&                a_ref_rat       ,
                      const RealVect&                   a_dxCrse             ,
                      IceThicknessIBC*                  a_bc                 ,
                      int                               a_numLevels          )
  {
    m_coarseDomain       =  a_coarseDomain    ;
    m_constRel           =  a_constRel        ;
    m_basalFrictionRel   =  a_basalFrictionRel;
    m_grids          =  a_vectGrids       ;
    m_ref_rat       =  a_ref_rat    ;
    m_dxCrse             =  a_dxCrse          ;
    m_bc                 =  a_bc              ;
    m_numLevels          =  a_numLevels       ;
  }

  virtual int refToFiner(const ch_dom& a_indexSpace) const
  {
    int retval = 4586;
    bool found = false;
    ch_dom dom_lev = m_coarseDomain;
    for(int ilev = 0; ilev < m_numLevels; ilev++)
    {
      if(dom_lev == a_indexSpace)
      {
        found = true;
        retval = m_ref_rat[ilev];
      }
      ///Chombo changes in place here
      dom_lev.refine(m_ref_rat[ilev]);
    }
    if(!found)
    {
      ch_mayday::Error("ChF_Ice_VTO_Factory:domain not found");
    }
    return retval;
  }

  ///
  virtual ~ChF_Ice_VTO_Factory()
  {
  }

  /// 
  virtual ch_amrlevelop* AMRnewOp(const ch_dom& a_indexSpace)
  {
    int start_level = 4586;
    bool found = false;
    ch_dom dom_lev = m_coarseDomain;
    for(int ilev = 0; ilev < m_numLevels; ilev++)
    {
      if(a_indexSpace == dom_lev)
      {
        start_level = ilev;
        found = true;
      }
          
    }
      
    ch_amrlevelop* retval = NULL;
    if(found)
    {
      int refRatioToAMR = 1;
      ch_ice_op* ice_op_ptr
        = new ch_ice_op(  m_coarseDomain    ,
                          m_constRel        ,
                          m_basalFrictionRel,
                          m_grids           ,
                          m_ref_rat          ,
                          m_dxCrse          ,
                          m_bc              ,
                          m_numLevels       ,
                          start_level        ,
                          refRatioToAMR     );
      retval = static_cast<ch_amrlevelop*>(ice_op_ptr);
    }
    return retval;
  }

  ///
  virtual ch_mglevelop* MGnewOp(const ch_dom  & a_indexSpace,
                                int             a_depth,
                                bool            a_homoOnly = true)
  {
    int start_level = 4586;
    bool found = false;
    ch_dom dom_lev = m_coarseDomain;
    for(int ilev = 0; ilev < m_numLevels; ilev++)
    {
      if(a_indexSpace == dom_lev)
      {
        start_level = ilev;
        found = true;
      }
      dom_lev.refine(m_ref_rat[ilev]);
          
    }
    bool coarsenable = false;
    int refRatioToAMR = 4586;
    if(found)
    {
      const ch_dbl & startGrids = m_grids[start_level];
      int mg_ref_ratio = 1;  
      for(int img = 0; img < a_depth; img++)
      {
        mg_ref_ratio *= 2;
      }
      refRatioToAMR= mg_ref_ratio;
      //we want at least 2 points in each direction
      coarsenable = (startGrids.coarsenable(2*mg_ref_ratio));
    }
    ch_amrlevelop* retval = NULL;
    if(found && coarsenable)
    {
      int depth = 0;
      ch_ice_op* ice_op_ptr
        = new ch_ice_op(  m_coarseDomain       ,
                          m_constRel           ,
                          m_basalFrictionRel   ,
                          m_grids          ,
                          m_ref_rat       ,
                          m_dxCrse             ,
                          m_bc                 ,
                          m_numLevels          ,
                          start_level,  refRatioToAMR);
      retval = static_cast<ch_amrlevelop*>(ice_op_ptr);
    }
    return retval;
  }
    
private:
  ///The elusive quest for strong construction continues.
  ChF_Ice_VTO_Factory();
public:
    
};
  
typedef ChF_Ice_VTO_Factory ch_ice_fact;
///
/**
   This class is a pass-through to use AMRFASMultiGrid to solve the 
   non-linear viscous tensor equation for ice velocity.   
   The only member data of this class are just retentions of what comes in through the constructor.
**/
class ChF_FAS_Ice_Solver : public IceVelocitySolver
{
  typedef                               FArrayBox       ch_fab_cell;
  typedef                     LevelData< ch_fab_cell >  ch_ldf_cell;
  typedef                       DisjointBoxLayout       ch_dbl;
  typedef                           ProblemDomain       ch_dom;
  typedef                           DataIterator        ch_dit;
  typedef                                  MayDay       ch_mayday;
  typedef                                     Box       ch_box;
  typedef                       LevelFluxRegister       ch_flux_reg;
  typedef     Base_Simplified_AMRLevelOp< ch_ldf_cell > ch_amrlevelop;
  typedef                AMRFASMultiGrid< ch_ldf_cell > ch_amrfasmg;
  typedef                      MGLevelOp< ch_ldf_cell > ch_mglevelop;
  typedef                   LinearSolver< ch_ldf_cell > ch_linear_solver;
  typedef                    RelaxSolver< ch_ldf_cell > ch_relax_solver;
  typedef                 BiCGStabSolver< ch_ldf_cell > ch_bicgstab_solver;
  typedef              AMRLevelOpFactory< ch_ldf_cell > ch_amrlevelop_fact;
private:


  ///strong construction only 
  ChF_FAS_Ice_Solver();
  
public:
  ///
  /**
     The only member data for this class is the stuff that comes in through the constructor.   
  **/
  ProblemDomain              m_coarseDomain       ;
  ConstitutiveRelation*      m_constRel           ;
  BasalFrictionRelation*     m_basalFrictionRel   ;
  Vector<DisjointBoxLayout>  m_grids              ;
  Vector<int>                m_ref_rat             ;
  RealVect                   m_dxCrse             ;
  IceThicknessIBC*           m_bc                 ;
  int                        m_numLevels          ;


  ///strong construction only 
  ChF_FAS_Ice_Solver(const ProblemDomain&              a_coarseDomain    ,
                     ConstitutiveRelation*             a_constRel        ,
                     BasalFrictionRelation*            a_basalFrictionRel,
                     const Vector<DisjointBoxLayout>&  a_grids           ,
                     const Vector<int>&                a_ref_rat          ,
                     const RealVect&                   a_dxCrse          ,
                     IceThicknessIBC*                  a_bc              ,
                     int                               a_numLevels       )
  :IceVelocitySolver()
  {
    m_coarseDomain       =  a_coarseDomain    ;
    m_constRel           =  a_constRel        ;
    m_basalFrictionRel   =  a_basalFrictionRel;
    m_grids              =  a_grids       ;
    m_ref_rat             =  a_ref_rat    ;
    m_dxCrse             =  a_dxCrse          ;
    m_bc                 =  a_bc              ;
    m_numLevels          =  a_numLevels       ;
  }

  virtual ~ChF_FAS_Ice_Solver()
  {
  }
  
  /// full solve for non-isothermal ice
  virtual int solve(Vector<ch_ldf_cell* >& a_horizontalVel,
                    Vector<ch_ldf_cell* >& a_calvedIce, // not used
                    Vector<ch_ldf_cell* >& a_addedIce,  // not used
                    Vector<ch_ldf_cell* >& a_removedIce,//not used
                    Real& a_initialResidualNorm, //poached from AMRFASMultiGrid
                    Real& a_finalResidualNorm,   //poached from AMRFASMultiGrid
                    const Real a_convergenceMetric,
                    const Vector<ch_ldf_cell* >& a_rhs,      
                    const Vector<ch_ldf_cell* >& a_C,
                    const Vector<ch_ldf_cell* >& a_C0,
                    const Vector<ch_ldf_cell* >& a_A,
                    const Vector<ch_ldf_cell* >& a_muCoef,
                    Vector<RefCountedPtr<LevelSigmaCS > >& a_coordSys,
                    Real a_time,
                    int a_lbase, int a_maxLevel)
  {
    ch_ice_fact   ice_factory(m_coarseDomain       ,
                              m_constRel           ,
                              m_basalFrictionRel   ,
                              m_grids              ,
                              m_ref_rat            ,
                              m_dxCrse             ,
                              m_bc                 ,
                              m_numLevels          );
    ch_relax_solver relax_solver;

    ch_amrfasmg fas_solver;
    fas_solver.define(m_coarseDomain, ice_factory, &relax_solver, m_numLevels);

    bool zeroPhi = true; bool forceHomog = false;
    fas_solver.solve(a_horizontalVel, a_rhs, a_maxLevel, a_lbase, zeroPhi, forceHomog);
    fas_solver.recoverInitialAndFinalRNorms(a_initialResidualNorm, a_finalResidualNorm);
    /**
       AMRMultiGrid::solve has a void return value so I am tempted to return something whimsical.
       This value gets does checked  downstream, so I cannot return  anything interesting.
    **/                
    return 0; 
  } //end function solve

    ///
    /**
       This function is mandated by the IceVelocitySolver API, 
       so this function is needed to get stuff to compile.
       That said, I never call this particular function, so it is just a placeholder.
    **/
  virtual void define(const ProblemDomain& a_coarseDomain,
                      ConstitutiveRelation* a_constRel,
                      BasalFrictionRelation* a_FrictionRel,
                      const Vector<DisjointBoxLayout>& a_vectGrids,
                      const Vector<int>& a_ref_rat,
                      const RealVect& a_dxCrse,
                      IceThicknessIBC* a_bc,
                      int a_numLevels)
  {
    ch_mayday::Error("How did this get called?");
  }
}; //end  class fas ice solver 


#include "NamespaceFooter.H"
#endif
