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
  typedef                               FArrayBox  ch_fab;
  typedef                     LevelData< ch_fab >  ch_ldf;
  typedef                       DisjointBoxLayout  ch_dbl;
  typedef                           ProblemDomain  ch_dom;
  typedef                           DataIterator   ch_dit;
  typedef                                  MayDay  ch_mayday;
  typedef                                     Box  ch_box;
  typedef                       LevelFluxRegister  ch_flux_reg;
  typedef                AMRFASMultiGrid< ch_ldf > ch_amrfasmg;
  typedef                      MGLevelOp< ch_ldf > ch_mglevelop;
  typedef                   LinearSolver< ch_ldf > ch_linear_solver;
  typedef                    RelaxSolver< ch_ldf > ch_relax_solver;
  typedef                 BiCGStabSolver< ch_ldf > ch_bicgstab_solver;
private:
  ///strong construction only
  ChF_Ice_Viscous_Tensor_Op();

  ///defining convenience data for this class
  void
  defineLocalData()
  {
    m_dom = m_coarseDomain;
    m_dx = m_dxCrse[0];
    for(int ilev = 0; ilev < m_startLevel; ilev++)
    {
      m_dom.refine(m_refRat[ilev]);
      m_dx /= m_refRat[ilev];
    }
    m_dbl = m_grids[m_startLevel];
    if(m_refRatioToAMR  > 1)
    {
      m_dom.coarsen(m_refRatioToAMR);
      m_dbl.coarsen(m_refRatioToAMR);
      m_dx *= m_refRatioToAMR;
    }

    m_grad.define(m_dbl, SpaceDim*SpaceDim, IntVect::Unit);
      
      
  }
        

  ///This is  the stuff that came in during construction
  ProblemDomain                m_coarseDomain       ;
  ConstitutiveRelation*        m_constRel           ;
  BasalFrictionRelation*       m_basalFrictionRel   ;
  Vector<DisjointBoxLayout>    m_grids              ;
  Vector<int>                  m_refRat             ;
  RealVect                     m_dxCrse             ;
  IceThicknessIBC*             m_bc                 ;
  int                          m_numLevels          ;
  int                          m_startLevel         ;
  int                          m_refRatioToAMR      ;
  ///this is the local member data (defined at construction)
  /**
     I am trying to keep this at a minimum but more may creep in here for performance reasons.
  **/
  ch_dom                       m_dom;
  ch_dbl                       m_dbl;
  ch_ldf                       m_grad;
  Real                         m_dx;

      
public:
  ///destructor  is trivial but has to be virtual
  virtual ~ChF_Ice_Viscous_Tensor_Op(){;}

  /// 
  ChF_Ice_Viscous_Tensor_Op(const ch_dom                   &  a_coarseDomain       ,
                            ConstitutiveRelation           *  a_constRel           ,
                            BasalFrictionRelation          *  a_basalFrictionRel   ,
                            const Vector<DisjointBoxLayout>&  a_grids          ,
                            const Vector<int>              &  a_refRat       ,
                            const RealVect                 &  a_dxCrse             ,
                            IceThicknessIBC                *  a_bc                 ,
                            int                               a_numLevels          ,
                            int                               a_startLevel         ,
                            int                               a_refRatioToAMR      )
  {

    m_coarseDomain       = a_coarseDomain     ;
    m_constRel           = a_constRel         ;
    m_basalFrictionRel   = a_basalFrictionRel ;
    m_grids              = a_grids            ;
    m_refRat             = a_refRat           ;
    m_dxCrse             = a_dxCrse           ;
    m_bc                 = a_bc               ;
    m_numLevels          = a_numLevels        ;
    m_startLevel         = a_startLevel       ;
    m_refRatioToAMR      = a_refRatioToAMR    ;
    defineLocalData();
  }

  ///mostly poached from lib/src/AMRElliptic/ViscousTensorOp.cpp
  virtual void
  prolongIncrement(ch_ldf      & a_phiThisLevel,
                   const ch_ldf& a_correctCoarse)
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
  getDomainFromData(const ch_ldf & a_data)
  {
    return a_data.disjointBoxLayout().physDomain();
  }
  ///
  /**
     I am indulging in the fantasy that I can finally write a bc function that makes more sense.  
     This fills corners outside the domain with something reasonable.    
     As the Boss said, love is bigger than a Honda.  It really is bigger than a Subaru.
  **/
  void reflectBCSide(ch_fab      &      a_state,
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
  fillGhostDataOutsideDomain(ch_ldf          & a_lhs,
                             const ch_dom    & a_domain) 
  {
    CH_TIME("CHF_ICE_VTO::fillGhostDataOutsideDomain");
    ch_dbl dbl = a_lhs.disjointBoxLayout();
    ch_dit dit = a_lhs.dataIterator();
    int nghost = a_lhs.ghostVect()[0];
    for(int ibox = 0; ibox < dit.size(); ibox++)
    {
      const ch_box& valid =   dbl[dit[ibox]];
      ch_fab      & state = a_lhs[dit[ibox]];
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
  fillGhostDataOverFineFineInterface(ch_ldf &  a_phi) 
  {
    CH_TIME("CHF_ICE_VTO::fillDataOverFineFineInterface");
    a_phi.exchange();
  }

  ///
  virtual void
  copyTo(ch_ldf       &  a_dst,
         const ch_ldf &  a_src)
  {
    CH_TIME("CHF_ICE_VTO::copyTo");
    a_src.copyTo(a_src.interval(), a_dst, a_dst.interval());
  }
    

  ///mostly poached from lib/src/AMRElliptic/ViscousTensorOp.cpp
  virtual void
  fillGhostDataOverCoarseFineInterface(ch_ldf       &  a_phiFine, 
                                       const ch_ldf &  a_phiCoar)
  {
    CH_TIME("CHF_ICE_VTO::fillGhostDataOverCoarsefineInterface");
    Real    dx_fine = m_dx;
    Box     do_fine = m_dom.domainBox();
    int  refToCoar = m_refRat[m_startLevel-1];
    setToZero(m_grad);
    int ncomp = SpaceDim;
    DisjointBoxLayout grids_fine = a_phiFine.disjointBoxLayout();
    DisjointBoxLayout grids_coar = a_phiCoar.disjointBoxLayout();
    TensorCFInterp interpWithCoarser;
    interpWithCoarser.define(grids_fine, &grids_coar,
                             dx_fine, refToCoar, ncomp,
                             do_fine);
    interpWithCoarser.coarseFineInterp(a_phiFine, m_grad, a_phiCoar);
      
  }

  ///mostly poached from lib/src/AMRElliptic/ViscousTensorOp.cpp
  void
  cellGrad(ch_fab        &      a_gradPhi,
           const ch_fab  &      a_phi,
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
  fillGrad(const ch_ldf&       a_phi)
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
  /// called by the coarser level
  virtual void
  refluxFineCorrection(ch_ldf       &   a_lphCoar,
                       const ch_ldf &   a_phiFine, 
                       const ch_ldf &   a_phiCoar,
                       ch_amrlevelop*   a_finerBase)
  {
    CH_TIME("CHF_ICE_VTO::refluxFineCorrection");

    int refToFine = m_refRat[m_startLevel];
    ch_dbl dbl_fine = a_phiFine.disjointBoxLayout();
    ch_dbl dbl_coar = a_phiCoar.disjointBoxLayout();
    ChF_Ice_Viscous_Tensor_Op*
      op_f = dynamic_cast<ChF_Ice_Viscous_Tensor_Op*>(a_finerBase);
    /**
       VTO skips this step (because of context, I guess).
       Ergo, the next two lines  may be unnecessary.
    **/
    this->fillGrad(a_phiCoar);
    ch_ldf& phi_fine_cast = (ch_ldf&)a_phiFine;
    op_f->fillGhostDataOverCoarseFineInterface(phi_fine_cast, a_phiCoar);
    op_f->fillGrad(a_phiFine);
    ch_dom dom_fine = refine(m_dom, refToFine);
      
    ch_flux_reg flux_reg(dbl_fine, dbl_coar, dom_fine, refToFine, SpaceDim);
    flux_reg.setToZero();
    ///HERE I guess I need to figure out the coefficients now
#if 0      
    DataIterator dit = a_phiCoar.dataIterator();
    for (dit.reset(); dit.ok(); ++dit)
    {
      const FArrayBox& coarfab    =   a_phi    [dit];
      const FluxBox&   coareta    = (*m_eta)   [dit];
      const FluxBox&   coarlambda = (*m_lambda)[dit];
      const Box&       gridBox    =   a_phi.getBoxes()[dit];

      for (int idir = 0; idir < SpaceDim; idir++)
      {
        FArrayBox coarflux;
        Box faceBox = surroundingNodes(gridBox, idir);
        //fillgrad was called before this when applyOp was called
        getFlux(coarflux, coarfab, m_grad[dit()], coareta[idir], coarlambda[idir], faceBox, idir);

        Real scale = 1.0;
        levfluxreg.incrementCoarse(coarflux, scale,dit(),
                                   interv,interv,idir);
      }
    }
#endif    
  }
  ///For this one, the derived class may assume all ghost cells have been set correctly
  virtual void 
  applyOperator(ch_ldf        & a_lph,
                const ch_ldf  & a_phi)
  {
    CH_TIME("CHF_ICE_VTO::applyOperator");
    ch_mayday::Error("not implmented");
  }

  /// Restrict fine data to coarse data using volume-weighted averaging.   
  virtual void
  restriction(ch_ldf & a_resCoar,
              ch_ldf & a_resFine,
              int ref_rat)
  {
    CH_TIME("CHF_ICE_VTO::restriction");
    ch_mayday::Error("not implmented");
  }

  ///from LinearOp
  virtual void
  create(ch_ldf      & a_lhs,
         const ch_ldf& a_rhs)
  {
    CH_TIME("CHF_ICE_VTO::create");
    ch_mayday::Error("not implmented");
  }

  ///from LinearOp
  virtual void
  assign(ch_ldf      & a_lhs,
         const ch_ldf& a_rhs) 
  {
    CH_TIME("CHF_ICE_VTO::assign");
    ch_mayday::Error("not implmented");
  }
    
  ///from LinearOp
  virtual Real
  dotProduct(const ch_ldf      & a_1,
             const ch_ldf      & a_2)     
  {
    CH_TIME("CHF_ICE_VTO::dotProduct");
    Real retval = 4586;
    ch_mayday::Error("not implmented");
    return retval;
  }


  ///from LinearOp
  virtual void
  axby(ch_ldf      & a_lhs,
       const ch_ldf& a_x,
       const ch_ldf& a_y,
       Real          a_a,
       Real          a_b)
  {
    CH_TIME("CHF_ICE_VTO::axby");
    ch_mayday::Error("not implmented");
  }

  ///from LinearOp
  virtual Real
  norm(const ch_ldf & a_rhs, int a_ord)
  {
    CH_TIME("CHF_ICE_VTO::norm");
    Real retval = 4586;
    ch_mayday::Error("not implmented");
    return retval;
  }
    
  ///from LinearOp
  virtual void
  setToZero(ch_ldf & a_lhs)
  {
    CH_TIME("CHF_ICE_VTO::setToZero");
    ch_mayday::Error("not implmented");
  }

  ///from MGLevelOp
  virtual void
  createCoarser(ch_ldf      & a_coarse,
                const ch_ldf& a_fine,
                bool ghosted)
  {
    CH_TIME("CHF_ICE_VTO::createCoarser");
    ch_mayday::Error("not implmented");
  }
    
  ///from MGLevelOp
  virtual void
  relax(ch_ldf      & a_correction,
        const ch_ldf& a_residual,
        int a_iterations)
  {
    CH_TIME("CHF_ICE_VTO::relax");
    ch_mayday::Error("not implmented");
  }

  virtual int refToCoarser()
  {
    int retval = 4586;
    ch_mayday::Error("not implmented");
    return retval;
  }
    
  virtual void
  AMRRestrict(ch_ldf      & a_resCoarse,
              const ch_ldf& a_residual,
              const ch_ldf& a_correction,
              const ch_ldf& a_coarseCorrection,
              bool a_skip_res)
  {
    CH_TIME("CHF_ICE_VTO::AMRRestrict");
    ch_mayday::Error("not implmented");
  }

  ///
  /** a_correction += I[2h->h](a_coarseCorrection) */
  virtual void
  AMRProlong(ch_ldf      & a_correction,
             const ch_ldf& a_coarseCorrection)
  {
    CH_TIME("CHF_ICE_VTO::AMRProlong");
    ch_mayday::Error("not implmented");
  }

  ///
  /** a_residual = a_residual - L(a_correction, a_coarseCorrection) */
  virtual void
  AMRUpdateResidual(ch_ldf      & a_residual,
                    const ch_ldf& a_correction,
                    const ch_ldf& a_coarseCorrection)
  {
    CH_TIME("CHF_ICE_VTO::AMRUpdateResidual");
    ch_mayday::Error("not implmented");
  }

  ///
  /**
   */
  virtual void createCoarsened(ch_ldf&       a_lhs,
                               const ch_ldf& a_rhs,
                               const int&     a_refRat)
  {
    CH_TIME("CHF_ICE_VTO::createCoarsened");
    ch_mayday::Error("not implmented");
  }
}; //end class ChF_Ice_Viscous_Tensor_Op
  
  
/// class to make ChF_Ice_VTO's.   
class ChF_Ice_VTO_Factory: public AMRLevelOpFactory< LevelData<FArrayBox> >
{
public:
  typedef               ChF_Ice_Viscous_Tensor_Op  ch_ice_op;
  typedef                               FArrayBox  ch_fab;
  typedef                     LevelData< ch_fab >  ch_ldf;
  typedef                       DisjointBoxLayout  ch_dbl;
  typedef                           ProblemDomain  ch_dom;
  typedef                           DataIterator   ch_dit;
  typedef                                  MayDay  ch_mayday;
  typedef                                     Box  ch_box;
  typedef                       LevelFluxRegister  ch_flux_reg;
  typedef     Base_Simplified_AMRLevelOp< ch_ldf > ch_amrlevelop;
  typedef                AMRFASMultiGrid< ch_ldf > ch_amrfasmg;
  typedef                      MGLevelOp< ch_ldf > ch_mglevelop;
  typedef                   LinearSolver< ch_ldf > ch_linear_solver;
  typedef                    RelaxSolver< ch_ldf > ch_relax_solver;
  typedef                 BiCGStabSolver< ch_ldf > ch_bicgstab_solver;
  typedef              AMRLevelOpFactory< ch_ldf > ch_amrlevelop_fact;
  ///
  /**
     The only member data for this class is the stuff that comes in through the constructor.   
  **/
  ProblemDomain              m_coarseDomain       ;
  ConstitutiveRelation*      m_constRel           ;
  BasalFrictionRelation*     m_basalFrictionRel   ;
  Vector<DisjointBoxLayout>  m_grids          ;
  Vector<int>                m_refRat       ;
  RealVect                   m_dxCrse             ;
  IceThicknessIBC*           m_bc                 ;
  int                        m_numLevels          ;


  ///
  ChF_Ice_VTO_Factory(const ProblemDomain&              a_coarseDomain       ,
                      ConstitutiveRelation*             a_constRel           ,
                      BasalFrictionRelation*            a_basalFrictionRel   ,
                      const Vector<DisjointBoxLayout>&  a_vectGrids          ,
                      const Vector<int>&                a_refRat       ,
                      const RealVect&                   a_dxCrse             ,
                      IceThicknessIBC*                  a_bc                 ,
                      int                               a_numLevels          )
  {
    m_coarseDomain       =  a_coarseDomain    ;
    m_constRel           =  a_constRel        ;
    m_basalFrictionRel   =  a_basalFrictionRel;
    m_grids          =  a_vectGrids       ;
    m_refRat       =  a_refRat    ;
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
        retval = m_refRat[ilev];
      }
      ///Chombo changes in place here
      dom_lev.refine(m_refRat[ilev]);
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
    int startLevel = 4586;
    bool found = false;
    ch_dom dom_lev = m_coarseDomain;
    for(int ilev = 0; ilev < m_numLevels; ilev++)
    {
      if(a_indexSpace == dom_lev)
      {
        startLevel = ilev;
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
                          m_refRat          ,
                          m_dxCrse          ,
                          m_bc              ,
                          m_numLevels       ,
                          startLevel        ,
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
    int startLevel = 4586;
    bool found = false;
    ch_dom dom_lev = m_coarseDomain;
    for(int ilev = 0; ilev < m_numLevels; ilev++)
    {
      if(a_indexSpace == dom_lev)
      {
        startLevel = ilev;
        found = true;
      }
      dom_lev.refine(m_refRat[ilev]);
          
    }
    bool coarsenable = false;
    int refRatioToAMR = 4586;
    if(found)
    {
      const ch_dbl & startGrids = m_grids[startLevel];
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
                          m_refRat       ,
                          m_dxCrse             ,
                          m_bc                 ,
                          m_numLevels          ,
                          startLevel,  refRatioToAMR);
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
  typedef                               FArrayBox  ch_fab;
  typedef                     LevelData< ch_fab >  ch_ldf;
  typedef                       DisjointBoxLayout  ch_dbl;
  typedef                           ProblemDomain  ch_dom;
  typedef                           DataIterator   ch_dit;
  typedef                                  MayDay  ch_mayday;
  typedef                                     Box  ch_box;
  typedef                       LevelFluxRegister  ch_flux_reg;
  typedef     Base_Simplified_AMRLevelOp< ch_ldf > ch_amrlevelop;
  typedef                AMRFASMultiGrid< ch_ldf > ch_amrfasmg;
  typedef                      MGLevelOp< ch_ldf > ch_mglevelop;
  typedef                   LinearSolver< ch_ldf > ch_linear_solver;
  typedef                    RelaxSolver< ch_ldf > ch_relax_solver;
  typedef                 BiCGStabSolver< ch_ldf > ch_bicgstab_solver;
  typedef              AMRLevelOpFactory< ch_ldf > ch_amrlevelop_fact;
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
  Vector<int>                m_refRat             ;
  RealVect                   m_dxCrse             ;
  IceThicknessIBC*           m_bc                 ;
  int                        m_numLevels          ;


  ///strong construction only 
  ChF_FAS_Ice_Solver(const ProblemDomain&              a_coarseDomain    ,
                     ConstitutiveRelation*             a_constRel        ,
                     BasalFrictionRelation*            a_basalFrictionRel,
                     const Vector<DisjointBoxLayout>&  a_grids           ,
                     const Vector<int>&                a_refRat          ,
                     const RealVect&                   a_dxCrse          ,
                     IceThicknessIBC*                  a_bc              ,
                     int                               a_numLevels       )
  :IceVelocitySolver()
  {
    m_coarseDomain       =  a_coarseDomain    ;
    m_constRel           =  a_constRel        ;
    m_basalFrictionRel   =  a_basalFrictionRel;
    m_grids              =  a_grids       ;
    m_refRat             =  a_refRat    ;
    m_dxCrse             =  a_dxCrse          ;
    m_bc                 =  a_bc              ;
    m_numLevels          =  a_numLevels       ;
  }

  virtual ~ChF_FAS_Ice_Solver()
  {
  }
  
  /// full solve for non-isothermal ice
  virtual int solve(Vector<ch_ldf* >& a_horizontalVel,
                    Vector<ch_ldf* >& a_calvedIce, // not used
                    Vector<ch_ldf* >& a_addedIce,  // not used
                    Vector<ch_ldf* >& a_removedIce,//not used
                    Real& a_initialResidualNorm, //poached from AMRFASMultiGrid
                    Real& a_finalResidualNorm,   //poached from AMRFASMultiGrid
                    const Real a_convergenceMetric,
                    const Vector<ch_ldf* >& a_rhs,      
                    const Vector<ch_ldf* >& a_C,
                    const Vector<ch_ldf* >& a_C0,
                    const Vector<ch_ldf* >& a_A,
                    const Vector<ch_ldf* >& a_muCoef,
                    Vector<RefCountedPtr<LevelSigmaCS > >& a_coordSys,
                    Real a_time,
                    int a_lbase, int a_maxLevel)
  {
    ch_ice_fact   ice_factory(m_coarseDomain       ,
                              m_constRel           ,
                              m_basalFrictionRel   ,
                              m_grids              ,
                              m_refRat             ,
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
                      const Vector<int>& a_refRat,
                      const RealVect& a_dxCrse,
                      IceThicknessIBC* a_bc,
                      int a_numLevels)
  {
    ch_mayday::Error("How did this get called?");
  }
}; //end  class fas ice solver 


#include "NamespaceFooter.H"
#endif
