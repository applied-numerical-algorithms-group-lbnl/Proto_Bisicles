
#include "AMRMultiGrid_Proto_Interface.H"
#include "ProtoInterface.H"
#include "PrCh_Utilities.H"
#include "PrCh_LoHiCenterStencil.H"
#include "Vector.H"

#ifndef __PrCh_IceVTO_H__
#define __PrCh_IceVTO_H__
/// The new namespace that all the cool operators want in on.   
namespace PrCh_AMR_Elliptic
{
  ///trying to stick with mostly Proto data holders here
  typedef Proto::Stencil<double>                                                              pr_sten;
  typedef Proto::InterpStencil<double>                                                        pr_interp_sten;
  typedef Proto::Point                                                                        pr_pt;
  typedef Proto::Box                                                                          pr_box;
  typedef Proto::BoxData<double, 1  >                                                         pr_box_data_sca;
  typedef Proto::BoxData<double, DIM>                                                         pr_box_data_vec;
  typedef Proto::BoxData<double, DIM*DIM>                                                     pr_box_data_ten;
  typedef Proto::DisjointBoxLayout                                                            pr_dbl;
  typedef Proto::Shift                                                                        pr_shift;
  typedef Proto::LevelBoxData<double,    1    , Proto::MEMTYPE_DEFAULT,  Proto::PR_CELL   >   pr_cell_data_sca;
  typedef Proto::LevelBoxData<double,    1    , Proto::MEMTYPE_DEFAULT,  Proto::PR_FACE_0 >   pr_xfac_data_sca;
  typedef Proto::LevelBoxData<double,    1    , Proto::MEMTYPE_DEFAULT,  Proto::PR_FACE_1 >   pr_yfac_data_sca;
  typedef Proto::LevelBoxData<double,    1    , Proto::MEMTYPE_DEFAULT,  Proto::PR_FACE_2 >   pr_zfac_data_sca;
  typedef Proto::LevelBoxData<double,  DIM    , Proto::MEMTYPE_DEFAULT,  Proto::PR_CELL   >   pr_cell_data_vec;
  typedef Proto::LevelBoxData<double,  DIM*DIM, Proto::MEMTYPE_DEFAULT,  Proto::PR_CELL   >   pr_cell_data_ten;
  typedef Proto::LevelBoxData<double,  DIM    , Proto::MEMTYPE_DEFAULT,  Proto::PR_FACE_0 >   pr_xfac_data_vec;
  typedef Proto::LevelBoxData<double,  DIM    , Proto::MEMTYPE_DEFAULT,  Proto::PR_FACE_1 >   pr_yfac_data_vec;
  typedef Proto::LevelBoxData<double,  DIM    , Proto::MEMTYPE_DEFAULT,  Proto::PR_FACE_2 >   pr_zfac_data_vec;
  ///AMRMultigrid forces me to use some of these.
  typedef Chombo::ProblemDomain                                                               ch_dom; // 
  typedef Chombo::MayDay                                                                      ch_mayday;
  typedef Chombo::AMRLevelOpFactory<pr_cell_data_vec>                                         ch_amrlevelop_fact_vec;
  typedef Chombo::AMRLevelOp<       pr_cell_data_vec>                                         ch_amrlevelop_vec;
  typedef Chombo::MGLevelOp<        pr_cell_data_vec>                                         ch_mglevelop_vec;

  typedef PrCh_Lo_Hi_Center::LHC_MetaData     lo_hi_center_t;
  typedef Chombo::AMRFASOp<pr_cell_data_vec> ch_amrfas_op;

  ///
  /**
     pr_cell_data_vec derived class for FAS operator
  **/
  class AMRFAS_Proto_Vec_Op : public ch_amrfas_op
  {
  public:  
    ///
    AMRFAS_Proto_Vec_Op( int a_order,const ch_dbl &a_grid )
      :AMRFASOp<pr_cell_data_vec >( a_order, a_grid )
    {
    }
    /// - default to noop
    virtual void restrictState( RefCountedPtr<AMRFASOp<pr_cell_data_vec > >,
                                Copier &a_copier 
      ) 
    {
    }

    ///return did I do an exchange on phi as side effect
    virtual bool computeState( RefCountedPtr<pr_cell_data_vec > a_phi, 
                               const RefCountedPtr<pr_cell_data_vec > a_CrsPhi,
                               const RefCountedPtr<pr_cell_data_vec > a_FinePhi)
    {
      return false; 
    }

    ///
    virtual void apply( pr_cell_data_vec& a_p,
                        const pr_cell_data_vec& a_phi,
                        const pr_cell_data_vec *a_phiCoarse = 0,
                        bool a_doExchange = true ) = 0;


    ///
    virtual void create(pr_cell_data_vec& a_lhs,
                        const pr_cell_data_vec& a_rhs)
    {
      ch_mayday::Error("not implemented");
    }

    ///
    virtual void assign(pr_cell_data_vec&       a_lhs,
                        const pr_cell_data_vec& a_rhs)
    {
      ch_mayday::Error("not implemented");
    }

    ///
    virtual void assignLocal(pr_cell_data_vec&       a_lhs,
                             const pr_cell_data_vec& a_rhs)
    {
      ch_mayday::Error("not implemented");
    }

    ///
    virtual void assignCopier(pr_cell_data_vec&       a_lhs,
                              const pr_cell_data_vec& a_rhs,
                              const Copier&               a_copier)
    {

      ch_mayday::Error("not implemented");
    }

    ///
    virtual void zeroCovered(pr_cell_data_vec& a_lhs,
                             pr_cell_data_vec& a_rhs,
                             const Copier&         a_copier)
    {
      ch_mayday::Error("not implemented");
    }
    

    ///
    virtual Real dotProduct(const pr_cell_data_vec& a_1,
                            const pr_cell_data_vec& a_2)
    {
      ch_mayday::Error("not implemented");
      return 4586.;
    }
    
    ///multiple dot products (for GMRES) 
    virtual void mDotProduct(const pr_cell_data_vec& a_1,
                             const int a_sz,
                             const pr_cell_data_vec a_2[],
                             Real a_mdots[])
    {
      ch_mayday::Error("not implemented");
    }

    ///
    virtual void incr(pr_cell_data_vec&       a_lhs,
                      const pr_cell_data_vec& a_x,
                      Real                        a_scale)
    {
      ch_mayday::Error("not implemented");
    }

    ///
    virtual void axby(pr_cell_data_vec&       a_lhs,
                      const pr_cell_data_vec& a_x,
                      const pr_cell_data_vec& a_y,
                      Real                        a_a,
                      Real                        a_b)
    {
      ch_mayday::Error("not implemented");
    }


    ///
    virtual void mult( pr_cell_data_vec& a_x,
                       const pr_cell_data_vec& a_y )
    {
      ch_mayday::Error("not implemented");
    }

    ///
    virtual void scale(pr_cell_data_vec& a_lhs,
                       const Real&           a_scale)
    {
      ch_mayday::Error("not implemented");
    }

    ///
    virtual Real norm(const pr_cell_data_vec& a_x,
                      int                         a_ord) const
    {
      ch_mayday::Error("not implemented");
      return 4586.;
    }

    ///
    virtual Real norm( const pr_cell_data_vec& a_x,
                       int   a_ord,
                       int a_comp ) const
    {
      ch_mayday::Error("not implemented");
      return 4586.;
    }

    ///
    virtual Real localMaxNorm(const pr_cell_data_vec& a_x) const
    {
      ch_mayday::Error("not implemented");
    }

    ///
    virtual void setToZero( pr_cell_data_vec& a_x)
    {
      ch_mayday::Error("not implemented");
    }

    ///
    void write( const pr_cell_data_vec* a_data,
                const char*                 a_filename)
    {
      ch_mayday::Error("not implemented");
    }
  
    /// R and P, with explict types, protected, the real system implementations
    /**
       If you say so. --dtg
    **/
    virtual void AMRProlong(pr_cell_data_vec&       a_fineU,
                            const pr_cell_data_vec& a_CrsU,
                            pr_cell_data_vec&       a_temp,
                            RefCountedPtr<AMRFASOp<pr_cell_data_vec > > a_crsOp,
                            FAS_PROLONG_type a_type )
    {
      ch_mayday::Error("not implemented");
    }

    ///
    virtual void AMRRestrict( pr_cell_data_vec& a_CrsU,     
                              const pr_cell_data_vec& a_fineU,
                              FAS_RESTRICT_type a_type ) const
    {
      ch_mayday::Error("not implemented");
    }

    /// R and P, without explict types, API, call real methods with default params
    /**
       If you say so. --dtg
    **/
    virtual void AMRProlong(pr_cell_data_vec&       a_fineU,
                            const pr_cell_data_vec& a_CrsU,
                            pr_cell_data_vec&       a_temp,
                            RefCountedPtr<AMRFASOp<pr_cell_data_vec > > a_crsOp )
    {
      ch_mayday::Error("not implemented");
    }

    ///
    virtual void AMRFMGProlong(pr_cell_data_vec&       a_fineU,
                               const pr_cell_data_vec& a_CrsU,
                               pr_cell_data_vec&       a_temp,
                               RefCountedPtr<AMRFASOp<pr_cell_data_vec > > a_crsOp )
    {
      ch_mayday::Error("not implemented");
    }

  
    ///
    virtual void AMRRestrict( pr_cell_data_vec& a_CrsU,     
                              const pr_cell_data_vec& a_fineU,
                              pr_cell_data_vec& a_crsCover,
                              const Copier &a_copier ) const
    {
      ch_mayday::Error("not implemented");
    }
  }; //end class AMRFAS_Proto_Vec_Op

  ///
  /**
   **/
  class PrCh_IceVTO: public   AMRFAS_Proto_Vec_Op
  {
  public:
    ///
    PrCh_IceVTO( int a_o, 
                 const ch_dbl &a_grid,
                 const ConstitutiveRelation*  a_constRelPtr,
                 const BasalFrictionRelation* a_basalFrictionRelPtr,
                 IceThicknessIBC* a_bc)
    {
      ch_mayday::Error("not implemented");
    }

    virtual void apply( pr_cell_data_vec& a_p,
                        const pr_cell_data_vec& a_phi,
                        const pr_cell_data_vec *a_phiCoarse = 0,
                        bool a_doExchange = true )
    {
      ch_mayday::Error("not implemented");
    }
    ///Because there a factory class, we can dispense with weak construction
    PrCh_IceVTO();


    virtual ~PrCh_IceVTO()
    {
    }

  }; //end class PrCh_IceVTO    


///
/**
   Factory to create AMRFAS_LDFOps
*/
  class PrCh_AMRFASOpVecFactory: public AMRFASOpFactory<pr_cell_data_vec >
  {
  public:
    PrCh_AMRFASOpVecFactory( int a_nc = 1, int a_order = 2 ) : 
      AMRFASOpFactory<pr_cell_data_vec >( a_nc, a_order )
    {
      ch_mayday::Error("not implemented");
    }				 
    virtual ~PrCh_AMRFASOpVecFactory()
    {
    }

    ///
    void AMRNewOp( const int, 
                   RefCountedPtr<AMRFAS_LDFOp>,
                   bool a_isSR = false
      )
    {
      ch_mayday::Error("not implemented");
    }				 
    
  };

 ///
  class PrCh_FASIceVTOFactory: public PrCh_AMRFASOpVecFactory
  {
  public:
    ///
    PrCh_FASIceVTOFactory( ConstitutiveRelation* a_constRelPtr,
                           BasalFrictionRelation* a_basalFrictionRelPtr,
                           IceThicknessIBC* a_bc
      )
    {
      ch_mayday::Error("not implemented");
    }

    ///
    virtual ~PrCh_FASIceVTOFactory()    
    {
    }

    ///
    virtual void define( const ProblemDomain& a_coarseDomain, 
                         const RealVect&      a_crsDx,
                         const Vector<ch_dbl>& a_grids,
                         const Vector<int>&   a_refRatios,
                         int a_nSRGrids = 0
      )
    {
      ch_mayday::Error("not implemented");
    }
    

    ///
    virtual RefCountedPtr<AMRFASOp<pr_cell_data_vec > > 
    AMRNewOp(int a_ilev, 
             const ch_dbl& a_grid,
             bool a_isSR = false )
    {
      ch_mayday::Error("not implemented");
    }

    const ConstitutiveRelation* m_constRelPtr;
    const BasalFrictionRelation* m_basalFrictionRelPtr;

    ViscousTensorOpFactory* m_VTOFactory;
    IceThicknessIBC* m_bc; // this masks m_bc in base class

  };

/// virtual class to manage the nonlinear solve for the ice-sheet
/// momentum 
  class FASIceSolver : public IceVelocitySolver, public AMRFAS<pr_cell_data_vec >
  {

  public:

    FASIceSolver(// int a_sizeA
      ) // : m_sizeA(a_sizeA)
    {
    }

    virtual ~FASIceSolver()
    {
    }

    virtual void define(const ProblemDomain& a_coarseDomain,
                        ConstitutiveRelation* a_constRel,
                        BasalFrictionRelation* a_basalFrictionRel,
                        const Vector<ch_dbl>& a_vectGrids,
                        const Vector<int>& a_vectRefRatio,
                        const RealVect& a_dxCrse,
                        IceThicknessIBC* a_bc,
                        int a_numLevels);
  
    virtual void setTolerance( Real a_tolerance ) { m_rtol = a_tolerance; }

    /// set "absolute tolerance"
    /** If norm(residual) drops below the absolute tolerance (regardless 
        of the initial residual), we declare ourselves "done" -- this is a 
        simple way to avoid solver hangs.
    */
    //virtual void setAbsoluteTolerance(Real a_tolerance) {m_atol = a_tolerance;}

    /// sets maximum number of FAS iterations
    virtual void setMaxIterations(int a_max_iter) { m_max_iter = a_max_iter; }
  
    virtual void setVerbosity(int a_verbosity) {
      IceVelocitySolver::m_verbosity = a_verbosity;
      AMRFAS<pr_cell_data_vec >::m_verbosity = a_verbosity;
    }
  
    /// full solve for non-isothermal ice
    virtual int solve(Vector<pr_cell_data_vec* >& a_horizontalVel,
                      Vector<pr_cell_data_vec* >& a_calvedIce,
                      Vector<pr_cell_data_vec* >& a_addedIce,
                      Vector<pr_cell_data_vec* >& a_removedIce,
                      Real& a_initialResidualNorm, Real& a_finalResidualNorm,
                      const Real a_convergenceMetric,
                      const Vector<pr_cell_data_vec* >& a_rhs,
                      const Vector<pr_cell_data_vec* >& a_C,
                      const Vector<pr_cell_data_vec* >& a_C0,
                      const Vector<pr_cell_data_vec* >& a_A,
                      const Vector<pr_cell_data_vec* >& a_muCoef,
                      Vector<RefCountedPtr<LevelSigmaCS > >& a_coordSys,
                      Real a_time,
                      int a_lbase, int a_maxLevel)
    {
      ch_mayday::Error("not implemented");
    }
  };
  
}    //end namespace Pr_Ch_AMR_AMRElliptic
#endif
