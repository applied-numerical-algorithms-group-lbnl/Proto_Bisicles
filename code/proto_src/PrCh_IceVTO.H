
#include "AMRMultiGrid_Proto_Interface.H"
#include "PrChUtilities.H"
#include "ProtoInterface.H"
#include "Proto_LoHiCenterStencil.H"
#include "Vector.H"

#ifndef __PrCh_IceVTO_H__
#define __PrCh_IceVTO_H__
/****************/
PROTO_KERNEL_START 

unsigned int axbyIceVTOF(int                       a_pt[DIM],
                        Proto::Var<Real, DIM>     a_lph,
                        Proto::Var<Real, DIM>     a_x,
                        Proto::Var<Real, DIM>     a_y,
                        Real                      a_a,
                        Real                      a_b,
                        int                       a_iprint)
{
  for(int idir = 0; idir < DIM; idir++)
  {
    double xval = a_x(idir);
    double yval = a_y(idir);
    double lph  = a_a*xval + a_b*yval;
    a_lph(idir) = lph;
  }
  return 0;

}
PROTO_KERNEL_END(axbyIceVTOF, axbyIceVTO) 
///a_flux = eta( grad u + grad u T) + lambda I div U;
PROTO_KERNEL_START 
unsigned int  fluxFuncIceVTOF(int                     a_pt[DIM],
                             Proto::Var<Real, 1>     a_flux,
                             Proto::Var<Real, 1>     a_grad_u,
                             Proto::Var<Real, 1>     a_grad_u_t,
                             Proto::Var<Real, 1>     a_div_u,
                             Proto::Var<Real, 1>     a_eta,
                             Proto::Var<Real, 1>     a_lambda,
                             int                     a_face_dir,
                             int                     a_vel_comp)
{
  a_flux(0) = a_eta(0)*(a_grad_u(0) + a_grad_u_t(0));
  if(a_face_dir == a_vel_comp)
  {
    a_flux(0) += a_lambda(0)*a_div_u(0);
  }
  return 0;
}
PROTO_KERNEL_END(fluxFuncIceVTOF, fluxFuncIceVTO)


///a_lph = alpha*acoef*phi + beta(xincr + yincr + zincr)
PROTO_KERNEL_START 
unsigned int  operaIceVTOF(int                       a_pt[DIM],
                          Proto::Var<Real, DIM>     a_lph,
                          Proto::Var<Real, DIM>     a_phi,
                          Proto::Var<Real,  1 >     a_acoef,
                          Proto::Var<Real, DIM>     a_xincr,
                          Proto::Var<Real, DIM>     a_yincr,
                          Proto::Var<Real, DIM>     a_zincr,
                          Real                      a_alpha,
                          Real                      a_beta)
{
  //zincr is set to zero in 2d
  for(int idir = 0; idir < DIM; idir++)
  {
    double  divergence = a_xincr(idir) + a_yincr(idir) + a_zincr(idir);
    double  phival     =   a_phi(idir);  
    double  acoval     = a_acoef(0);   //same acoef for all variables
    double  lphval     = a_alpha*acoval*phival + a_beta*divergence;
    a_lph(idir) = lphval;
  }
  return 0;
}
PROTO_KERNEL_END(operaIceVTOF, operaIceVTO)

/// 
/**
   I am modelling this on lib/src/ViscousTensorOp but I will
   not account for differences in coefficients from one side of the cell
   to the other.   
   The VTO fortran amounts to 
   invlam(icomp) = acoef*a_alpha;
   for(int face_dir = 0; face_dir < SpaceDim; face_dir++)
   {
   double incr = 0;
   if(icomp == face_dir)
   {
   incr = beta*(2*eta(i+1/2) + lambda(i+1/2)
   +2*eta(i-1/2) + lambda(i-1/2))/(dx*dx);
   }
   else
   {
   incr = beta*(eta(i+1/2) + 2*eta(i-1/2));
   }
      
   inv_lam(idir) -= incr;
   }
   lam(idir) = -1/invlam(idir)

   Forall likes things to be at the same point so we will
   We will approximate this by;
   invlam(icomp) = acoef*a_alpha;
   for(int face_dir = 0; face_dir < SpaceDim; face_dir++)
   {
   double incr = 0;
   if(icomp == face_dir)
   {
   incr = beta*(4*eta(i-1/2) + 2*lambda(i-1/2))/(dx*dx);
                     
   }
   else
   {
   incr = 2*beta*eta(i-1/2);
   }
      
   inv_lam(idir) -= incr;
   }
   lam(idir) = 1/invlam(idir)
**/
PROTO_KERNEL_START 
unsigned int  relaxCoefIceVTOF(int                        a_pt[DIM],
                              Proto::Var<Real,  DIM>     a_relCoef,
                              Proto::Var<Real,  1  >     a_acoef,
                              Proto::Var<Real,  1  >     a_eta_x,
                              Proto::Var<Real,  1  >     a_eta_y,
                              Proto::Var<Real,  1  >     a_eta_z,
                              Proto::Var<Real,  1  >     a_lam_x,
                              Proto::Var<Real,  1  >     a_lam_y,
                              Proto::Var<Real,  1  >     a_lam_z,
                              double                     a_alpha,
                              double                     a_beta,
                              double                     a_dx)
{
  //this this is approximate, might have to have a real safety here
  double safety = 0.9; //same as ViscousTensorOpF.ChF
  double inv_lam[DIM];
  double eta_val[DIM];
  double lam_val[DIM];
  eta_val[0] = a_eta_x(0) ;
  lam_val[0] = a_lam_x(0) ;
  
  eta_val[1] = a_eta_y(0) ;
  lam_val[1] = a_lam_y(0) ;
#if DIM==3
  eta_val[2] = a_eta_z(0) ;
  lam_val[2] = a_lam_z(0) ;
#endif
  for(int icomp = 0; icomp < DIM; icomp++)
  {
    inv_lam[icomp] = a_alpha*a_acoef(0);
    for(int face_dir = 0; face_dir < DIM; face_dir++)
    {
      double incr = 0;
      if(icomp == face_dir)
      {
        incr = a_beta*(4.*eta_val[face_dir] + 2.*lam_val[face_dir])/(a_dx*a_dx);
      }
      else
      {
        incr = a_beta*(2.*eta_val[face_dir])/(a_dx*a_dx);
      }
      inv_lam[icomp] -= incr;
    }//end loop over face_dir
    a_relCoef(icomp) = safety/inv_lam[icomp];
  }

  return 0;
}
PROTO_KERNEL_END(relaxCoefIceVTOF, relaxCoefIceVTO)


///
/**
   phi = phi + lambda*(lphi - rho) 
   but only where 
   sum(i+j+k)%2== iredblack
**/
PROTO_KERNEL_START 
unsigned int  gsrbIceVTOF(int                       a_pt[DIM],
                         Proto::Var<Real, DIM>     a_phi,
                         Proto::Var<Real, DIM>     a_res,
                         Proto::Var<Real, DIM>     a_lam,
                         int                       a_iredBlack)
{
  int sumpt = 0;
  for(int idir = 0; idir < DIM; idir++)
  {
    sumpt += a_pt[idir];
  }
  if(sumpt%2 == a_iredBlack)
  {
    for(int icomp = 0; icomp < DIM; icomp++)
    {
      Real lambda = a_lam(icomp);
      Real phival = a_phi(icomp);
      Real resval = a_res(icomp);

      a_phi(icomp) = phival + lambda*resval;
    }
  }
  return 0;

}
PROTO_KERNEL_END(gsrbIceVTOF, gsrbIceVTO) 
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

///
/**
   pr_cell_data_vec derived class for FAS operator
*/
class AMRFAS_Proto_Vec_Op : public AMRFASOp<pr_cell_data_vec>
{
public:  
  //! Constructor.
  AMRFAS_Proto_Vec_Op( int a_order, const DisjointBoxLayout &a_grid ) : 
    AMRFASOp<pr_cell_data_vec >( a_order, a_grid )
  {
  }
  /// - default to noop
  virtual void restrictState( RefCountedPtr<AMRFASOp<pr_cell_data_vec > >,
			      Copier &a_copier 
			      ) 
  {
  }

  virtual bool computeState( RefCountedPtr<pr_cell_data_vec > a_phi, 
			     const RefCountedPtr<pr_cell_data_vec > a_CrsPhi,
			     const RefCountedPtr<pr_cell_data_vec > a_FinePhi 
			     )
  {
    return false; // did I do an exchange on phi (side effect)
  }

  virtual void apply( pr_cell_data_vec& a_p,
		      const pr_cell_data_vec& a_phi,
		      const pr_cell_data_vec *a_phiCoarse = 0,
		      bool a_doExchange = true ) = 0;

  // general vector ops, with T=pr_cell_data_vec

  LevelDataOps<FArrayBox> m_levelOps;
  
  virtual void create(pr_cell_data_vec& a_lhs,
		      const pr_cell_data_vec& a_rhs);

  virtual void assign(pr_cell_data_vec&       a_lhs,
                      const pr_cell_data_vec& a_rhs);

  virtual void assignLocal(pr_cell_data_vec&       a_lhs,
                           const pr_cell_data_vec& a_rhs);

  virtual void assignCopier(pr_cell_data_vec&       a_lhs,
                            const pr_cell_data_vec& a_rhs,
                            const Copier&               a_copier);

  virtual void zeroCovered(pr_cell_data_vec& a_lhs,
                           pr_cell_data_vec& a_rhs,
                           const Copier&         a_copier);

  virtual Real dotProduct(const pr_cell_data_vec& a_1,
                          const pr_cell_data_vec& a_2);
  /* multiple dot products (for GMRES) */
  virtual void mDotProduct(const pr_cell_data_vec& a_1,
                           const int a_sz,
                           const pr_cell_data_vec a_2[],
                           Real a_mdots[]);

  virtual void incr(pr_cell_data_vec&       a_lhs,
                    const pr_cell_data_vec& a_x,
                    Real                        a_scale);

  virtual void axby(pr_cell_data_vec&       a_lhs,
                    const pr_cell_data_vec& a_x,
                    const pr_cell_data_vec& a_y,
                    Real                        a_a,
                    Real                        a_b);
  
  virtual void mult( pr_cell_data_vec& a_x,
		     const pr_cell_data_vec& a_y );
  
  virtual void scale(pr_cell_data_vec& a_lhs,
                     const Real&           a_scale);

  virtual Real norm(const pr_cell_data_vec& a_x,
                    int                         a_ord) const;

  virtual Real norm( const pr_cell_data_vec& a_x,
		     int   a_ord,
		     int a_comp ) const;

  virtual Real localMaxNorm(const pr_cell_data_vec& a_x) const ;

  virtual void setToZero( pr_cell_data_vec& a_x) ;

  void write( const pr_cell_data_vec* a_data,
	      const char*                 a_filename);
  
  // R and P, with explict types, protected, the real system implementations
  virtual void AMRProlong(pr_cell_data_vec&       a_fineU,
			  const pr_cell_data_vec& a_CrsU,
			  pr_cell_data_vec&       a_temp,
			  RefCountedPtr<AMRFASOp<pr_cell_data_vec > > a_crsOp,
			  FAS_PROLONG_type a_type );
  virtual void AMRRestrict( pr_cell_data_vec& a_CrsU,     
			    const pr_cell_data_vec& a_fineU,
			    FAS_RESTRICT_type a_type ) const;

  // R and P, without explict types, API, call real methods with default params
public:
  virtual void AMRProlong(pr_cell_data_vec&       a_fineU,
			  const pr_cell_data_vec& a_CrsU,
			  pr_cell_data_vec&       a_temp,
			  RefCountedPtr<AMRFASOp<pr_cell_data_vec > > a_crsOp );

  virtual void AMRFMGProlong(pr_cell_data_vec&       a_fineU,
			  const pr_cell_data_vec& a_CrsU,
			  pr_cell_data_vec&       a_temp,
			  RefCountedPtr<AMRFASOp<pr_cell_data_vec > > a_crsOp );

  
  virtual void AMRRestrict( pr_cell_data_vec& a_CrsU,     
			    const pr_cell_data_vec& a_fineU,
			    pr_cell_data_vec& a_crsCover,
			    const Copier &a_copier ) const;
}; //end class AMRFAS_Proto_Vec_Op

  /// 
  class PrCh_IceVTO: public   AMRFAS_Proto_Vec_Op
  {
  public:


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
class AMRFAS_LDFOpFactory: public AMRFASOpFactory<pr_cell_data_vec >
{
public:
  AMRFAS_LDFOpFactory( int a_nc = 1, int a_order = 2 ) : 
    AMRFASOpFactory<pr_cell_data_vec >( a_nc, a_order )
  {
    ch_mayday::Error("not implemented");
  }				 
  virtual ~AMRFAS_LDFOpFactory()
  {
  }

protected:
  ///
  void AMRNewOp( const int, 
		 RefCountedPtr<AMRFAS_LDFOp>,
		 bool a_isSR = false
		 ); 
};
}    //end namespace Pr_Ch_AMR_AMRElliptic
#endif
