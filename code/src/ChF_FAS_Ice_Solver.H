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
#include <memory>  //good title for a book
///
namespace Chombo
{
  typedef           LevelData<FArrayBox>  ch_ldf;
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

  public:
    ///This is  the stuff that came in during construction
    ProblemDomain                m_coarseDomain       ;
    ConstitutiveRelation*        m_constRel           ;
    BasalFrictionRelation*       m_basalFrictionRel   ;
    Vector<DisjointBoxLayout>    m_vectGrids          ;
    Vector<int>                  m_refRat       ;
    RealVect                     m_dxCrse             ;
    IceThicknessIBC*             m_bc                 ;
    int                          m_numLevels          ;
    ch_dom                       m_indexSpace         ;
    int                          m_depth              ;
    
    ///destructor  is trivial but has to be virtual
    virtual ~ChF_Ice_Viscous_Tensor_Op(){;}

    /// 
    ChF_Ice_Viscous_Tensor_Op(const ch_dom                   &  a_coarseDomain       ,
                              ConstitutiveRelation           *  a_constRel           ,
                              BasalFrictionRelation          *  a_basalFrictionRel   ,
                              const Vector<DisjointBoxLayout>&  a_vectGrids          ,
                              const Vector<int>              &  a_refRat       ,
                              const RealVect                 &  a_dxCrse             ,
                              IceThicknessIBC                *  a_bc                 ,
                              int                               a_numLevels          ,
                              ch_dom                            a_indexSpace         ,
                              int                               a_depth              )
    {

      m_coarseDomain       = m_coarseDomain     ;
      m_constRel           = m_constRel         ;
      m_basalFrictionRel   = m_basalFrictionRel ;
      m_vectGrids          = m_vectGrids        ;
      m_refRat       = m_refRat     ;
      m_dxCrse             = m_dxCrse           ;
      m_bc                 = m_bc               ;
      m_numLevels          = m_numLevels        ;
      m_indexSpace         = m_indexSpace       ;
      m_depth              = m_depth            ;
    }
    virtual void
    prolongIncrement(ch_ldf      & a_phiThisLevel,
                     const ch_ldf& a_correctCoarse)
    {
      ch_mayday::Error("not implmented");
    }
    
    ///
    virtual ch_dom
    getDomainFromData(const ch_ldf & a_data)
    {
      ch_mayday::Error("not implmented");
    }

    ///
    virtual void
    fillGhostDataOutsideDomain(ch_ldf & a_phi,
                               const ch_dom    & a_domain) 
    {
      ch_mayday::Error("not implmented");
    }

    ///
    virtual void
    fillGhostDataOverFineFineInterface(ch_ldf &  a_phiFine) 
    {
      ch_mayday::Error("not implmented");
    }

    ///
    virtual void
    copyTo(ch_ldf       &  a_dst,
           const ch_ldf &  a_src)
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


    /// amr thing
    virtual void
    fillGhostDataOverCoarseFineInterface(ch_ldf       &  a_phiFine, 
                                         const ch_ldf &  a_phiCoar)
    {
      ch_mayday::Error("not implmented");
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
    Vector<DisjointBoxLayout>  m_vectGrids          ;
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
      m_vectGrids          =  a_vectGrids       ;
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
      ch_ice_op * ice_op_ptr
        = new ch_ice_op(  m_coarseDomain       ,
                          m_constRel           ,
                          m_basalFrictionRel   ,
                          m_vectGrids          ,
                          m_refRat       ,
                          m_dxCrse             ,
                          m_bc                 ,
                          m_numLevels          ,
                          a_indexSpace,  0);
      ch_amrlevelop* retval = static_cast<ch_ice_op*>(ice_op_ptr);
      return retval;
    }

    ///
    virtual ch_mglevelop* MGnewOp(const ch_dom       & a_FineindexSpace,
                                  int                  a_depth,
                                  bool                 a_homoOnly = true)
    {
      ch_ice_op * ice_op_ptr
        = new ch_ice_op(  m_coarseDomain       ,
                          m_constRel           ,
                          m_basalFrictionRel   ,
                          m_vectGrids          ,
                          m_refRat       ,
                          m_dxCrse             ,
                          m_bc                 ,
                          m_numLevels          ,
                          a_FineindexSpace,
                          a_depth );
      return NULL;
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
    Vector<DisjointBoxLayout>  m_vectGrids          ;
    Vector<int>                m_refRat       ;
    RealVect                   m_dxCrse             ;
    IceThicknessIBC*           m_bc                 ;
    int                        m_numLevels          ;


    ///strong construction only 
    ChF_FAS_Ice_Solver(const ProblemDomain&              a_coarseDomain       ,
                       ConstitutiveRelation*             a_constRel           ,
                       BasalFrictionRelation*            a_basalFrictionRel   ,
                       const Vector<DisjointBoxLayout>&  a_vectGrids          ,
                       const Vector<int>&                a_refRat       ,
                       const RealVect&                   a_dxCrse             ,
                       IceThicknessIBC*                  a_bc                 ,
                       int                               a_numLevels          )
    :IceVelocitySolver()
    {
      m_coarseDomain       =  a_coarseDomain    ;
      m_constRel           =  a_constRel        ;
      m_basalFrictionRel   =  a_basalFrictionRel;
      m_vectGrids          =  a_vectGrids       ;
      m_refRat       =  a_refRat    ;
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
                                m_vectGrids          ,
                                m_refRat       ,
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
