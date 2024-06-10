
#include "AMRMultiGrid_Proto_Interface.H"
#include "ProtoInterface.H"
#include "PrCh_Utilities.H"
#include "PrCh_LoHiCenterStencil.H"
#include "Vector.H"
#include "RefCountedPtr.H"
#include "Proto_FastVTO.H"

#ifndef __PrCh_IceVTO_H__
#define __PrCh_IceVTO_H__
/// The new namespace that all the cool operators want in on.   
namespace PrCh_AMR_Elliptic
{
  ///proto stuff
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
  
  ///chombo stuff
  typedef Chombo::ProblemDomain                                                               ch_dom; // 
  typedef Chombo::DisjointBoxLayout                                                           ch_dbl;
  typedef Chombo::MayDay                                                                      ch_mayday;

  typedef Chombo::LevelData<Chombo::FArrayBox>                                                ch_ldf;
  typedef Chombo::AMRFASOp<pr_cell_data_vec>                                                  ch_amrfas_op;
  typedef Chombo::AMRFASOpFactory<pr_cell_data_vec>                                           ch_amrfas_op_factory;
  typedef Chombo::Copier                                                                      ch_copier;
  typedef Chombo::FAS_PROLONG_type                                                            ch_fas_prolong_type;
  typedef Chombo::FAS_RESTRICT_type                                                           ch_fas_restrict_type;
  typedef Chombo::ConstitutiveRelation                                                        ch_consitutive_relation;
  typedef Chombo::BasalFrictionRelation                                                       ch_basal_friction_relation;
  typedef Chombo::IceThicknessIBC                                                             ch_ice_thickness_ibc;
  typedef Chombo::RealVect                                                                    ch_realvect;
  
  ///other stuff
  typedef PrCh_AMR_Elliptic::Proto_FastVTO_Factory                                            prch_fast_vto_factory;
  typedef PrCh_Lo_Hi_Center::LHC_MetaData                                                     lo_hi_center_t;
  using Chombo::RefCountedPtr;
  using Chombo::Vector;

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
                                ch_copier &a_copier 
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
                              const ch_copier&               a_copier)
    {

      ch_mayday::Error("not implemented");
    }

    ///
    virtual void zeroCovered(pr_cell_data_vec& a_lhs,
                             pr_cell_data_vec& a_rhs,
                             const ch_copier&         a_copier)
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
                            ch_fas_prolong_type a_type )
    {
      ch_mayday::Error("not implemented");
    }

    ///
    virtual void AMRRestrict( pr_cell_data_vec& a_CrsU,     
                              const pr_cell_data_vec& a_fineU,
                              ch_fas_restrict_type a_type ) const
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
                              const ch_copier &a_copier ) const
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
    PrCh_IceVTO( int a_order,  //called a_o in fortran version
                 const ch_dbl &a_grid,
                 const ch_consitutive_relation   *  a_constRelPtr,
                 const ch_basal_friction_relation*  a_basalFrictionRelPtr,
                 ch_ice_thickness_ibc*              a_bc)
      :AMRFAS_Proto_Vec_Op(a_order, a_grid)
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
  class PrCh_AMRFASOpVecFactory: public ch_amrfas_op_factory 
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
                   RefCountedPtr<Chombo::AMRFAS_LDFOp>,
                   bool a_isSR = false
      )
    {
      ch_mayday::Error("not implemented");
    }				 
    
  };

  ///because I loathe storing raw pointers.
  void null_deleter_void(void *) {}
  
  ///
  class PrCh_FASIceVTOFactory: public PrCh_AMRFASOpVecFactory
  {
  public:
    ///
    PrCh_FASIceVTOFactory( ch_consitutive_relation   * a_consitutive_relation_raw    ,
                           ch_basal_friction_relation* a_basal_friction_relation_raw ,
                           ch_ice_thickness_ibc      * a_ice_thickness_ibc_raw       )
    {
      m_consitutive_relation_ptr    =  std::shared_ptr< ch_consitutive_relation    >(a_consitutive_relation_raw    , &null_deleter_void);
      m_basal_friction_relation_ptr =  std::shared_ptr< ch_basal_friction_relation >(a_basal_friction_relation_raw , &null_deleter_void);
      m_ice_thickness_ibc_ptr       =  std::shared_ptr< ch_ice_thickness_ibc       >(a_ice_thickness_ibc_raw       , &null_deleter_void);  
      //m_fast_vto_factory_ptr        =  std::shared_ptr< prch_fast_vto_factory      >(new_vto_raw);                                           ;
      ch_mayday::Error("not implemented");
    }

    ///
    virtual ~PrCh_FASIceVTOFactory()    
    {
    }

    ///
    virtual void define( const ch_dom        &   a_coarseDomain, 
                         const ch_realvect   &   a_crsDx,
                         const Vector<ch_dbl>&   a_grids,
                         const Vector<int>   &   a_refRatios,
                         int a_nSRGrids = 0 )
    {
      ch_mayday::Error("not implemented");
    }
    

    ///
    virtual RefCountedPtr<Chombo::AMRFASOp<pr_cell_data_vec > > 
    AMRNewOp(int a_ilev, 
             const ch_dbl& a_grid,
             bool a_isSR = false )
    {
      ch_mayday::Error("not implemented");
    }

    std::shared_ptr< ch_consitutive_relation    > m_consitutive_relation_ptr;   
    std::shared_ptr< ch_basal_friction_relation > m_basal_friction_relation_ptr ;
    std::shared_ptr< ch_ice_thickness_ibc       > m_ice_thickness_ibc_ptr;
    std::shared_ptr< prch_fast_vto_factory      > m_fast_vto_factory_ptr;


  };

/// virtual class to manage the nonlinear solve for the ice-sheet
/// momentum 
  class Proto_FAS_IceSolver : public Chombo::IceVelocitySolver, public Chombo::AMRFAS<ch_ldf >
  {

  public:

    Proto_FAS_IceSolver(// int a_sizeA
      ) // : m_sizeA(a_sizeA)
    {
    }

    virtual ~Proto_FAS_IceSolver()
    {
    }

    virtual void define(const ch_dom& a_coarseDomain,
                        ch_consitutive_relation* a_constRel,
                        ch_basal_friction_relation* a_basalFrictionRel,
                        const Vector<ch_dbl>& a_vectGrids,
                        const Vector<int>& a_vectRefRatio,
                        const ch_realvect& a_dxCrse,
                        ch_ice_thickness_ibc* a_bc,
                        int a_numLevels)
    {
      ch_mayday::Error("not implemented");
    }
  
    //virtual void setTolerance( Real a_tolerance ) { m_rtol = a_tolerance; }

    /// set "absolute tolerance"
    /** If norm(residual) drops below the absolute tolerance (regardless 
        of the initial residual), we declare ourselves "done" -- this is a 
        simple way to avoid solver hangs.
    */
    //virtual void setAbsoluteTolerance(Real a_tolerance) {m_atol = a_tolerance;}

    /// sets maximum number of FAS iterations
   //virtual void setMaxIterations(int a_max_iter) { m_max_iter = a_max_iter; }
   //
   //virtual void setVerbosity(int a_verbosity)
   //{
   //  IceVelocitySolver::m_verbosity = a_verbosity;
   //  AMRFAS<pr_cell_data_vec >::m_verbosity = a_verbosity;
   //}
  
    /// full solve for non-isothermal ice
    virtual int solve(Vector<ch_ldf* >& a_horizontalVel,
                      Vector<ch_ldf* >& a_calvedIce,
                      Vector<ch_ldf* >& a_addedIce,
                      Vector<ch_ldf* >& a_removedIce,
                      Real& a_initialResidualNorm, Real& a_finalResidualNorm,
                      const Real a_convergenceMetric,
                      const Vector<ch_ldf* >& a_rhs,
                      const Vector<ch_ldf* >& a_C,
                      const Vector<ch_ldf* >& a_C0,
                      const Vector<ch_ldf* >& a_A,
                      const Vector<ch_ldf* >& a_muCoef,
                      Vector<RefCountedPtr<Chombo::LevelSigmaCS > >& a_coordSys,
                      Real a_time,
                      int a_lbase, int a_maxLevel)
    {
      ch_mayday::Error("not implemented");

      //copy data to device
      //solve on device
      //copy solution to host
    }
  };
  
}    //end namespace Pr_Ch_AMR_AMRElliptic
#endif
