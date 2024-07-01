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
#include "ViscousTensorOpF_F.H"
#include <memory>  //good title for a book
///
namespace Chombo
{
  typedef           LevelData<FArrayBox>  ch_ldf;
  typedef              DisjointBoxLayout  ch_dbl;
  typedef     Base_Simplified_AMRLevelOp< ch_ldf > ch_amrlevelop;
  typedef                AMRFASMultiGrid< ch_ldf > ch_amrfasmg;
  typedef                      MGLevelOp< ch_ldf > ch_mglevelop;
  typedef                   LinearSolver< ch_ldf > ch_linear_solver;
  typedef                    RelaxSolver< ch_ldf > ch_relax_solver;
  typedef                 BiCGStabSolver< ch_ldf > ch_bicgstab_solver;
  typedef              AMRLevelOpFactory< ch_ldf > ch_amrlevelop_fact;
  typedef   ProblemDomain ch_dom;
  typedef          MayDay ch_mayday;

  ///
  /**
     This is a pass-through class for ViscousTensorOp for the nonlinear ice problem.
     Original idea: Each operator owns a ViscousTensorOp.   It turns out this is a terrible idea.
     VTO precomputes tons of stuff.  
     I will poach heavily from  VTO.   I am also using the simpler Base_Simplified_AMRLevelOp
  **/
  class ChF_Ice_Viscous_Tensor_Op:  public ch_amrlevelop
  {
  private:
    ///strong construction only
    ChF_Ice_Viscous_Tensor_Op();

    ///defining convenience data for this class
    void
    defineLocalData()
    {
      m_dom = m_coarseDomain;
      for(int ilev = 0; ilev < m_startLevel; ilev++)
      {
        m_dom.refine(m_refRat[ilev]);
      }
      m_dbl = m_grids[m_startLevel];
      if(m_refRatioToAMR  > 1)
      {
        m_dom.coarsen(m_refRatioToAMR);
        m_dbl.coarsen(m_refRatioToAMR);
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

    ///
    virtual void
    prolongIncrement(ch_ldf      & a_phiThisLevel,
                     const ch_ldf& a_correctCoarse)
    {

      CH_TIME("ChF_FAS_Ice_VTO::prolongIncrement");
      
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
       I am indulging in thinking I can finally write a bc function that makes more sense.  This fills corners outside the domain
       with something reasonable.    As the Boss said, love is bigger than a Honda.  It's bigger than a Subaru.
    **/
    void homogDiriBCSide(ch_fab      &      a_state,
                         const ch_box&      a_valid,
                         const ch_dom&      a_dom,
                         int                a_dir,
                         int                a_side,
                         int                a_comp,
                         int                a_nghost)
    {

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
        for(int idir = 0; idir < SpaceDim; idir++)
        {
          if(idir != a_dir)
          {
            offdir_grown_dom_box.grow(idir, a_nghost);
          }
        }
        //only grow by one.  need the starting cell for a row of ghost
        ch_box start_box = a_state.box();
        comp_box.grow(1);
        comp_box.shift(a_dir, isign);
        comp_box &= offdir_grown_dom_box;
        comp_box.shift(a_dir, -isign);
        //this should cover (with something reasonable) corners outside the domain
        for (BoxIterator bit(start_box); bit.ok(); ++bit)
        {
          IntVect iv_first_ghost = bit();
          for(int ighost = 0; ighost < a_nghost; ighost++)
          {
            IntVect iv_src =  iv_first_ghost - (ighost+1)*isign*BASISV(a_dir);
            IntVect iv_dst =  iv_first_ghost + ( ighost )*isign*BASISV(a_dir);
            a_state(iv_dst, a_comp) = -a_state(iv_src, a_comp);
          }
        }
      }
    } // end function homogDiriSide
    /// I know, I know.
    void homogNeumBCSide(ch_fab      &      a_state,
                         const ch_box&      a_valid,
                         const ch_dom&      a_dom,
                         int                a_dir,
                         int                a_side,
                         int                a_comp,
                         int                a_nghost)
    {

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
        for(int idir = 0; idir < SpaceDim; idir++)
        {
          if(idir != a_dir)
          {
            offdir_grown_dom_box.grow(idir, a_nghost);
          }
        }
        //only grow by one.  need the starting cell for a row of ghost
        ch_box start_box = a_state.box();
        comp_box.grow(1);
        comp_box.shift(a_dir, isign);
        comp_box &= offdir_grown_dom_box;
        comp_box.shift(a_dir, -isign);
        //this should cover (with something reasonable) corners outside the domain
        for (BoxIterator bit(start_box); bit.ok(); ++bit)
        {
          IntVect iv_first_ghost = bit();
          for(int ighost = 0; ighost < a_nghost; ighost++)
          {
            IntVect iv_src =  iv_first_ghost - (ighost+1)*isign*BASISV(a_dir);
            IntVect iv_dst =  iv_first_ghost + ( ighost )*isign*BASISV(a_dir);
            a_state(iv_dst, a_comp) = a_state(iv_src, a_comp);
          }
        }
      }
    } // end function homogNeumSide

    ///
    /**
       From Base_Simplified_AMRLevelOP
     **/
    virtual void
    fillGhostDataOutsideDomain(ch_ldf &          a_lhs,
                               const ch_dom    & a_domain) 
    {
      /**
         Buried in PetscCompGridVTO:
         public: ///thank god
         bool m_bcDiri[CH_SPACEDIM][2][CH_SPACEDIM]; // dim, side, Comp: false = Neumann, true: Diri
         It is all sorted out by the  bool isDiri function.
         PetscCompGrid does not use the LoHiSide class.    It always was an overdesign, I guess.
      **/
      DisjointBoxLayout dbl = a_lhs.disjointBoxLayout();
      DataIterator      dit = a_lhs.dataIterator();
      int nghost = a_lhs.ghostVect()[0];
      for(int ibox = 0; ibox < dit.size(); ibox++)
      {
        const ch_box& valid =   dbl[dit[ibox]];
        const ch_fab& state = a_lhs[dit[ibox]];
      for(int iside = 0; iside < 1; iside++)
      {
        for(int idir = 0; idir < SpaceDim; idir++)
        {
          for(int icomp = 0; icomp < SpaceDim; icomp++)
          {
            bool is_diri = m_bc->isDiri(idir, iside, icomp);
            if(is_diri)
            {
              void homogDiriBCSide(state,
                                   valid,
                                   a_dom,
                                   idir,
                                   iside,
                                   icomp,
                                   nghost);
            }
            else
            {
              void homogNeumBCSide(state,
                                   valid,
                                   a_dom,
                                   idir,
                                   iside,
                                   icomp,
                                   nghost);
            }
          }
        }
      }
    }

    ///
    virtual void
    fillGhostDataOverFineFineInterface(ch_ldf &  a_phi) 
    {
      a_phi.exchange();
    }

    ///
    virtual void
    copyTo(ch_ldf       &  a_dst,
           const ch_ldf &  a_src)
    {
      a_src.copyTo(a_src.interval(), a_dst, a_dst.interval());
    }
    
    /// 
    virtual void
    fillGhostDataOverCoarseFineInterface(ch_ldf       &  a_phiFine, 
                                         const ch_ldf &  a_phiCoar)
    {
      TensorCFInterp          m_interpWithCoarser;
      
    }
    /// amr thing
    virtual void
    refluxFineCorrection(ch_ldf       &  a_lphCoar,
                         const ch_ldf &  a_phiFine, 
                         const ch_ldf &  a_phiCoar,
                         ch_amrlevelop*   a_finerOp)
    {
      ch_mayday::Error("not implmented");
    }
    
    ///For this one, the derived class may assume all ghost cells have been set correctly
    virtual void 
    applyOperator(ch_ldf        & a_lph,
                  const ch_ldf  & a_phi)
    {
      ch_mayday::Error("not implmented");
    }

    /// Restrict fine data to coarse data using volume-weighted averaging.   
    /**
       Restrict is a c++ keyword, so we are restricted from naming a function restrict.
       To avoid this restriction, I named the function restriction instead of restrict.
    **/
    virtual void
    restriction(ch_ldf & a_resCoar,
                ch_ldf & a_resFine,
                int ref_rat)
    {
      ch_mayday::Error("not implmented");
    }




    ///from LinearOp
    virtual void
    create(ch_ldf      & a_lhs,
           const ch_ldf& a_rhs)
    {
      ch_mayday::Error("not implmented");
    }

    ///from LinearOp
    virtual void
    assign(ch_ldf      & a_lhs,
           const ch_ldf& a_rhs) 
    {
      ch_mayday::Error("not implmented");
    }
    
    ///from LinearOp
    virtual Real
    dotProduct(const ch_ldf      & a_1,
               const ch_ldf      & a_2)     
    {
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
      ch_mayday::Error("not implmented");
    }

    ///from LinearOp
    virtual Real
    norm(const ch_ldf & a_rhs, int a_ord)
    {
      Real retval = 4586;
      ch_mayday::Error("not implmented");
      return retval;
    }
    
    ///from LinearOp
    virtual void
    setToZero(ch_ldf & a_lhs)
    {
      ch_mayday::Error("not implmented");
    }

    ///from MGLevelOp
    virtual void
    createCoarser(ch_ldf      & a_coarse,
                  const ch_ldf& a_fine,
                  bool ghosted)
    {
      ch_mayday::Error("not implmented");
    }
    
    ///from MGLevelOp
    virtual void
    relax(ch_ldf      & a_correction,
          const ch_ldf& a_residual,
          int a_iterations)
    {
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
      ch_mayday::Error("not implmented");
    }

    ///
    /** a_correction += I[2h->h](a_coarseCorrection) */
    virtual void
    AMRProlong(ch_ldf      & a_correction,
               const ch_ldf& a_coarseCorrection)
    {
      ch_mayday::Error("not implmented");
    }

    ///
    /** a_residual = a_residual - L(a_correction, a_coarseCorrection) */
    virtual void
    AMRUpdateResidual(ch_ldf      & a_residual,
                      const ch_ldf& a_correction,
                      const ch_ldf& a_coarseCorrection)
    {
      ch_mayday::Error("not implmented");
    }

    ///
    /**
     */
    virtual void createCoarsened(ch_ldf&       a_lhs,
                                 const ch_ldf& a_rhs,
                                 const int&     a_refRat)
    {
      ch_mayday::Error("not implmented");
    }
  }; //end class ChF_Ice_Viscous_Tensor_Op
  
  typedef ChF_Ice_Viscous_Tensor_Op ch_ice_op;
  /// class to make ChF_Ice_VTO's.   
  class ChF_Ice_VTO_Factory: public ch_amrlevelop_fact
  {
  public:
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

    ///
    virtual ch_mglevelop* MGnewOp(const ch_dom       & a_indexSpace,
                                  int                  a_depth,
                                  bool                 a_homoOnly = true)
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
    Vector<DisjointBoxLayout>  m_grids          ;
    Vector<int>                m_refRat       ;
    RealVect                   m_dxCrse             ;
    IceThicknessIBC*           m_bc                 ;
    int                        m_numLevels          ;


    ///strong construction only 
    ChF_FAS_Ice_Solver(const ProblemDomain&              a_coarseDomain       ,
                       ConstitutiveRelation*             a_constRel           ,
                       BasalFrictionRelation*            a_basalFrictionRel   ,
                       const Vector<DisjointBoxLayout>&  a_grids          ,
                       const Vector<int>&                a_refRat       ,
                       const RealVect&                   a_dxCrse             ,
                       IceThicknessIBC*                  a_bc                 ,
                       int                               a_numLevels          )
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
       Mandated by the IceVelocitySolver API.
       This should not get called as this class uses strong construction
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
      ch_mayday::Error("How did this get called after all my fiendish plans?");
    }
  }; //end  class fas ice solver 

} // end namespace chombo

#endif
