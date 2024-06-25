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
#include "AMRFASMultiGrid.H"
#include <memory>  //good title for a book
///
namespace Chombo
{
  typedef     AMRFASMultiGrid< LevelData<FArrayBox> > ch_amrfasmg;
  typedef          AMRLevelOp< LevelData<FArrayBox> > ch_amrlevelop;
  typedef           MGLevelOp< LevelData<FArrayBox> > ch_mglevelop;
  typedef        LinearSolver< LevelData<FArrayBox> > ch_linear_solver;
  typedef         RelaxSolver< LevelData<FArrayBox> > ch_relax_solver;
  typedef      BiCGStabSolver< LevelData<FArrayBox> > ch_bicgstab_solver;
  typedef   AMRLevelOpFactory< LevelData<FArrayBox> > ch_amrlevelop_fact;
  typedef   ProblemDomain ch_dom;
  typedef          MayDay ch_mayday;
  

  class ChF_Ice_VTO_Factory: public ch_amrlevelop_fact
  {
  public:
    ///
    ChF_Ice_VTO_Factory(const ProblemDomain&              a_coarseDomain       ,
                        ConstitutiveRelation*             a_constRel           ,
                        BasalFrictionRelation*            a_basalFrictionRel   ,
                        const Vector<DisjointBoxLayout>&  a_vectGrids          ,
                        const Vector<int>&                a_vectRefRatio       ,
                        const RealVect&                   a_dxCrse             ,
                        IceThicknessIBC*                  a_bc                 ,
                        int                               a_numLevels          )
    {
      ch_mayday::Error("not implemented");
    }

    virtual int refToFiner(const ch_dom& a_indexSpace) const
    {
      ch_mayday::Error("not implemented");
      return 4586;
    }
    ///
    virtual ~ChF_Ice_VTO_Factory()
    {
    }

    /// 
    virtual ch_amrlevelop* AMRnewOp(const ch_dom& a_indexSpace)
    {
      ch_mayday::Error("not implemented");
      return NULL;
    }

    ///
    virtual ch_mglevelop* MGnewOp(const ch_dom       & a_FineindexSpace,
                                  int                  a_depth,
                                  bool                 a_homoOnly = true)
    {
      ch_mayday::Error("not implemented");
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
 **/
  class ChF_FAS_Ice_Solver : public IceVelocitySolver
  {
  private:


    ///strong construction only 
    ChF_FAS_Ice_Solver();
  
  public:
    ///sent in through constructor
    ProblemDomain              m_coarseDomain       ;
    ConstitutiveRelation*      m_constRel           ;
    BasalFrictionRelation*     m_basalFrictionRel   ;
    Vector<DisjointBoxLayout>  m_vectGrids          ;
    Vector<int>                m_vectRefRatio       ;
    RealVect                   m_dxCrse             ;
    IceThicknessIBC*           m_bc                 ;
    int                        m_numLevels          ;
    Real                       m_tol                ;
    int                        m_maxIter            ;
    int                        m_verbosity          ;


    ///strong construction only 
    ChF_FAS_Ice_Solver(const ProblemDomain&              a_coarseDomain       ,
                       ConstitutiveRelation*             a_constRel           ,
                       BasalFrictionRelation*            a_basalFrictionRel   ,
                       const Vector<DisjointBoxLayout>&  a_vectGrids          ,
                       const Vector<int>&                a_vectRefRatio       ,
                       const RealVect&                   a_dxCrse             ,
                       IceThicknessIBC*                  a_bc                 ,
                       int                               a_numLevels          )
    :IceVelocitySolver()
    {
      m_coarseDomain       =  a_coarseDomain    ;
      m_constRel           =  a_constRel        ;
      m_basalFrictionRel   =  a_basalFrictionRel;
      m_vectGrids          =  a_vectGrids       ;
      m_vectRefRatio       =  a_vectRefRatio    ;
      m_dxCrse             =  a_dxCrse          ;
      m_bc                 =  a_bc              ;
      m_numLevels          =  a_numLevels       ;
    }

    virtual ~ChF_FAS_Ice_Solver()
    {
    }
  
    /// full solve for non-isothermal ice
    virtual int solve(Vector<LevelData<FArrayBox>* >& a_horizontalVel,
                      Vector<LevelData<FArrayBox>* >& a_calvedIce, // not used
                      Vector<LevelData<FArrayBox>* >& a_addedIce,  // not used
                      Vector<LevelData<FArrayBox>* >& a_removedIce,//not used
                      Real& a_initialResidualNorm, //poached from AMRFASMultiGrid
                      Real& a_finalResidualNorm,   //poached from AMRFASMultiGrid
                      const Real a_convergenceMetric,
                      const Vector<LevelData<FArrayBox>* >& a_rhs,      
                      const Vector<LevelData<FArrayBox>* >& a_C,
                      const Vector<LevelData<FArrayBox>* >& a_C0,
                      const Vector<LevelData<FArrayBox>* >& a_A,
                      const Vector<LevelData<FArrayBox>* >& a_muCoef,
                      Vector<RefCountedPtr<LevelSigmaCS > >& a_coordSys,
                      Real a_time,
                      int a_lbase, int a_maxLevel)
    {
      ch_amrfasmg   fas_solver;
      ch_ice_fact   ice_factory(m_coarseDomain       ,
                                m_constRel           ,
                                m_basalFrictionRel   ,
                                m_vectGrids          ,
                                m_vectRefRatio       ,
                                m_dxCrse             ,
                                m_bc                 ,
                                m_numLevels          );
      ch_relax_solver relax_solver;

      ch_amrfasmg fas_solver;
      fas_solver.define(m_coarseDomain, ice_factory, relax_solver, m_numLevels);

      bool zeroPhi = true; bool forceHomog = false;
      fas_solver.solve(a_horizontalVel, a_rhs, a_maxLevel, a_lbase, zeroPhi, forceHomog);
      a_initialResidualNorm = m_solver->m_initial_rnom;
      a_finalResidualNorm   = m_solver->m_final_rnom;
                      
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
                        const Vector<int>& a_vectRefRatio,
                        const RealVect& a_dxCrse,
                        IceThicknessIBC* a_bc,
                        int a_numLevels)
    {
      ch_mayday::Error("How did this get called after all my fiendish plans?");
    }
  }; //end  class fas ice solver 

} // end namespace chombo

#endif
