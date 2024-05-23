#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _COARSEFINEINTERP_H_
#define _COARSEFINEINTERP_H_

#include "AverageF_F.H"
#include "BoxIterator.H"
#include "DisjointBoxLayout.H"
#include "FArrayBox.H"
#include "IntVectSet.H"
#include "LevelData.H"
#include "ProblemDomain.H"
#include "RefCountedPtr.H"
#include "UsingNamespace.H"

///
/// Interpolation from coarse level to fine level for CELL-AVERAGED data.
/// Typical use case:
///  typedef CoarseFineInterp<> CFI;
///  std::Vector<int> t;
///  t[CFI::Homogeneous]      = 0;
///  t[CFI::PolynomialDegree] = 4;
///  t[CFI::RefineRatio]      = 2;
///   ... set all independent parameters this way.
///  CFI interpolator;
///  interpolator.define(...);
///  interpolator.interpolate(...);
///
/// Note for the special case that the fine box adjoins a non-periodic domain boundary:
///  the corner ghosts out of the physical domain will NOT be filled
///  even if FillCornerGhosts is set to true.



// Qinghai Zhang, August 13 2009
//  please report any bugs.

template<class Stencil>
class CoarseFineInterp
{
public:

  /// Parameters of this class
  /// The following is an iteration on all the parameters
  /// for (int i=0; i<nParameters; i++) {PARAMETER p = (PARAMETER) i; ...}
  enum PARAMETER {
    /// Is this interpolator homogeneous?
    Homogeneous = 0,
    /// The degree of the fitting polynomial
    PolynomialDegree,
    /// The number of coefficients in the polynomial
    NumPolynmCoefs,
    /// refinement ratio between this level and the next coarser
    RefineRatio,
    /// How many layers of fine ghosts to fill, = ghost width of LevelData.
    InterpRadius,
    /// CoarsenedInterpRadius = ceiling(InterpRadius/RefineRatio);
    CoarsenedInterpRadius,
    /// Proper nesting width.
    ///  it must satisfy (NestingWidth >= CoarsenedInterpRadius)
    ///  to guarantee that this is an interpolation.
    /// The stencil will be chosen within the box
    ///  grow(coarsenedFineBox,NestingWidth) & available data box.
    NestingWidth,
    /// Should corner ghosts be interpolated?
    /// In 3D, corner ghosts also include edge ghosts
    ///  that are not used in the standard 9-point Laplacian stencil
    FillCornerGhosts,
    /// In the non-homo case, interpolate ghost cells or the whole fine patch?
    InterpGhostsOnly,
    /// When there aren't sufficient cells for the specified order,
    /// this option determines whether to abort or try an lower order.
    InsistOnAccuracy,
    /// For homogeneous interpolation,
    ///  smoothers usually works fine with ghosts set to zero.
    /// This option means nothing for the non-homo case.
    UseZeroForHomo,
    /// Number of parameters for this class
    nParameters
  };

  /// output the name of a parameter.
  static char const*const toStr(const PARAMETER& p)
  {
    switch(p)
      {
      case Homogeneous           : return "Homogeneous          ";
      case PolynomialDegree      : return "PolynomialDegree     ";
      case NumPolynmCoefs        : return "NumPolynmCoefs       ";
      case RefineRatio           : return "RefineRatio          ";
      case InterpRadius          : return "InterpRadius         ";
      case CoarsenedInterpRadius : return "CoarsenedInterpRadius"; 
      case NestingWidth          : return "NestingWidth         ";
      case FillCornerGhosts      : return "FillCornerGhosts     ";
      case InterpGhostsOnly      : return "InterpGhostsOnly     ";
      case InsistOnAccuracy      : return "InsistOnAccuracy     ";
      case UseZeroForHomo        : return "UseZeroForHomo       ";
      default                    : return "";
      }
  }

  /// A dependent parameter can be derived from another,
  ///  thus not in the list of define().
  static bool dependent(const PARAMETER& p)
  {
    return ( (p==NumPolynmCoefs) || (p==CoarsenedInterpRadius) );
  }

  /// Return the value of a parameter from the enum list.
  int getParameter(const PARAMETER& p) const
  {
    return m_int[p];
  }

  /// Output all the parameters to a stream.
  friend std::ostream& operator<<(std::ostream& os,
                                  const CoarseFineInterp<Stencil>& cfi)
  {
    for (int i=0; i<nParameters; i++)
      {
        PARAMETER p = static_cast<PARAMETER>(i);
        os << "  " << toStr(p) << " = "
           << cfi.getParameter(p) << std::endl;
      }
    return os;
  }

  /// Default constructor.
  ///  define() has to be called before interpolate().
  CoarseFineInterp() : m_defined(false) {}

  /// Destructor
  virtual ~CoarseFineInterp(){}

  /// Actual constructor.
  /// The parameters passed here belong to LevelData's to be interpolated.
  /// Parameters starting from a_homogeneous are independent.
  void define(/// layout at this level
              const DisjointBoxLayout&  a_thisDisjointBoxLayout,
              /// layout at coarser level
              const DisjointBoxLayout&  a_coarserDisjointBoxLayout,
              /// problem domain on the finer level
              const ProblemDomain&      a_fineDomain,
              /// Is this interpolator homogeneous?
              const bool                a_homogeneous,
              /// The degree of the fitting polynomial
              const int                 a_polynomialDegree,
              /// refinement ratio between this level and the coarser level
              const int                 a_refineRatio,
              /// number of layers of ghost cells to fill by interpolation;
              /// it must be < the number of ghost layers of the fine data
              const int                 a_interpRadius,
              /// proper nesting width
              const int                 a_nestingWidth=3,
              /// Should corner ghosts be interpolated?
              const bool                a_cornersAreValid =true,
              /// if this is false, interpolate the whole ghosted fine patch.
              const bool                a_interpGhostsOnly=true,
              /// if this is true, abort when there is not enough cells
              const bool                a_insistOnAccuracy=false,
              /// set ghosts to zero for the homo case
              const bool                a_useZeroForHomo=true )
  {
    CH_TIME("CoarseFineInterp::define1");
    m_layout                  = a_thisDisjointBoxLayout;
    m_coarseLayout            = a_coarserDisjointBoxLayout;
    m_coarseDomain            = a_fineDomain;
    m_coarseDomain.coarsen(a_refineRatio);
    m_int[Homogeneous]        = a_homogeneous;
    m_int[PolynomialDegree]   = a_polynomialDegree;
    m_int[RefineRatio]        = a_refineRatio;
    m_int[InterpRadius]       = a_interpRadius;
    m_int[NestingWidth]       = a_nestingWidth;
    m_int[FillCornerGhosts]   = a_cornersAreValid;
    m_int[InterpGhostsOnly]   = a_interpGhostsOnly;
    m_int[InsistOnAccuracy]   = a_insistOnAccuracy;
    m_int[UseZeroForHomo]     = a_useZeroForHomo;

    initializeDependentMembers();
  }

  /// Syntactic sugar constructor .
  /// T must have an operator[] that maps the enum to values.
  ///  for example T can be Vector<int> or std::map<PARAMETER,int>.
  /// a_parameters must contain values of all the independent parameters,
  ///  i.e., the parameters in the other constructor.
  /// This is a good way to avoid bugs of wrong parameter order.
  template<class T>
  void define(const DisjointBoxLayout&  a_thisDisjointBoxLayout,
              const DisjointBoxLayout&  a_coarserDisjointBoxLayout,
              const ProblemDomain&      a_fineDomain,
              T&                        a_parameters )
  {
    CH_TIME("CoarseFineInterp::define2");
    // This is how enumerate helps to reduce code bloating.
    for (int i=0; i<nParameters; i++)
      {
        PARAMETER p = static_cast<PARAMETER>(i);
        if (!dependent(p))
          m_int[p] = a_parameters[p];
      }
    m_layout                = a_thisDisjointBoxLayout;
    m_coarseLayout          = a_coarserDisjointBoxLayout;
    m_coarseDomain          = coarsen(a_fineDomain, m_int[RefineRatio]);
    initializeDependentMembers();
  }

  ///
  /// Return a pointer to the stencil used at id for CoarseData,
  ///        0 if coarse cell id is not intended to be interpolated. 
  ///
  Stencil const*const getStencil(const IntVect& id) const
  {
    for (DataIterator dit=m_layout.dataIterator(); dit.ok(); ++dit)
      {
        if (m_coarsenedGhosts[dit].contains(id))
          return m_interpolators[dit](id);
      }
    // return 0 if no interpolator is found.
    return 0;
  }

  bool isDefined(void) const
  {
    return m_defined;
  }

  /// Interpolate at coarse-fine boundary
  /// If a_coarseData is NULL, do homogeneous interpolation,
  ///  otherwise, interpolate to a_fineData from a_coarseData.
  /// In the nonhomogeneous case,
  ///  the ghost cells of a_coarseData must be valid
  ///  since they might be used as interpolation sources.
  void interpolate(LevelData<FArrayBox>&            a_fineData,
                   LevelData<FArrayBox> const*const a_coarseData=0) const
  {
    CH_TIME("CoarseFineInterp::interpolate()");
    CH_assert(m_defined);
    if (m_int[Homogeneous]) CH_assert(!a_coarseData);

    IntVect gv = m_int[NestingWidth]*IntVect::Unit;
    for (int d=0; d<SpaceDim; d++)
      if (m_layout.physDomain().isPeriodic(d))
        gv[d] += m_int[NestingWidth];
    LevelData<FArrayBox> dataOnCoarsenedLayout
      (m_layoutCoarsened, a_fineData.nComp(), gv);

    // The coarse data comes from coarsening the fine data
    // and putting zeros to the boundary cells.
    if (m_int[Homogeneous])
      for (DataIterator dit = m_layout.dataIterator(); dit.ok(); ++dit)
        {
          const Box& bc = m_layoutCoarsened[dit];
          FArrayBox& fineFab = a_fineData[dit];
          FArrayBox& coarseFab = dataOnCoarsenedLayout[dit];
          coarseFab.setVal(0.0);
          // average the fine data to coarse data
          Box refbox(IntVect::Zero, (m_int[RefineRatio]-1)*IntVect::Unit);
          FORT_AVERAGE( CHF_FRA(coarseFab),
                        CHF_CONST_FRA(fineFab),
                        CHF_BOX(bc),
                        CHF_CONST_INT(m_int[RefineRatio]),
                        CHF_BOX(refbox) );
        }
    else
      // copy a_coarseData to fine layout with enough ghost layers.
      a_coarseData->copyTo(a_coarseData->interval(), dataOnCoarsenedLayout,
                           a_coarseData->interval());

    // try using zero to set the ghost cells in the homogeneous case.
    if (m_int[Homogeneous] && m_int[UseZeroForHomo])
      {
        LevelData<FArrayBox> fineCopy(a_fineData.disjointBoxLayout(),
                                      a_fineData.nComp(), IntVect::Zero);
        a_fineData.copyTo(fineCopy);
        for (DataIterator dit = m_layout.dataIterator(); dit.ok(); ++dit)
          {
            FArrayBox& fineFab = a_fineData[dit];
            fineFab.setVal(0.0);
          }
        fineCopy.copyTo(a_fineData);
      }
    else{

      for (DataIterator dit=m_layout.dataIterator(); dit.ok(); ++dit)
        // do the interpolation for each patch
        for (IVSIterator ivsit(m_coarsenedGhosts[dit]); ivsit.ok(); ++ivsit)
          m_interpolators[dit](ivsit())->fillFine
            (a_fineData[dit], dataOnCoarsenedLayout[dit], ivsit());
    }
  }


protected:

  /// define() has been called
  bool m_defined;

  /// value holder of the PARAMETER enum list
  int m_int[nParameters];

  /// Box layout for this level
  DisjointBoxLayout m_layout;     
  /// Coarsened layout of this level, the iterators are same for these two.
  DisjointBoxLayout m_layoutCoarsened;     
  /// Box layout for the coarse level
  DisjointBoxLayout m_coarseLayout;
  /// Box that contains the valid source data
  Box m_validSourceBox;
 
  /// problem domain at the coarser level
  ProblemDomain m_coarseDomain;

  /// coarsened ghost cells of fine patches
  LayoutData<IntVectSet> m_coarsenedGhosts;

  /// Indices of participating coarse cells uniquely determines the stencil.
  /// Using a map, the same interpolator coefs don't have to be recomputed.
  /// Reference counting is a must here.
  std::map<IntVectSet, RefCountedPtr<Stencil> >  m_interpMap;

  /// Associate each coarse cell with an interpolator in m_interpMap.
  LevelData<BaseFab<RefCountedPtr<Stencil> > >  m_interpolators;

private:

  /// assemble all the dependent members from independent members.
  void initializeDependentMembers()
  {
    CH_assert( m_int[RefineRatio]  >  1 );
    CH_assert( m_int[InterpRadius] >  0 );
    //    CH_assert( m_int[InterpRadius] <= m_int[NestingWidth]*m_int[RefineRatio] );
    CH_assert( m_int[PolynomialDegree] > 0 );

    m_int[NumPolynmCoefs] = numPolynmCoefs(m_int[PolynomialDegree]);
    m_int[CoarsenedInterpRadius] = m_int[InterpRadius]/m_int[RefineRatio];
    if (m_int[CoarsenedInterpRadius]*m_int[RefineRatio] < m_int[InterpRadius])
      m_int[CoarsenedInterpRadius]++;
    if(!m_int[Homogeneous] &&
       m_int[CoarsenedInterpRadius] > m_int[NestingWidth] )
      MayDay::Error("(NestingWidth >= CoarsenedInterpRadius) is a must!");

    // interpGhostOnly has to be true for the homogeneous case.
    if (m_int[Homogeneous])
      m_int[InterpGhostsOnly] = true;

    m_layoutCoarsened = DisjointBoxLayout();   // for re-define this object.
    coarsen(m_layoutCoarsened, m_layout, m_int[RefineRatio]);

    /// any cell inside the domain is a valid source
    /// extend the box to include NestingWidth layers of cell for periodic.
    m_validSourceBox = m_coarseDomain.domainBox();
    for (int d=0; d<SpaceDim; d++)
      if (m_layout.physDomain().isPeriodic(d))
        m_validSourceBox.grow(d, m_int[NestingWidth]);

    // initialization of more complicated members.
    initializeCoarsenedGhosts();
    initializeInterpolators();

    // Everything is defined now.
    m_defined = true;
  }


  /// Figure out the fine ghosts to be filled
  ///   and put the coarsened IntVect to m_coarsenedGhosts.
  /// Pay special attention to the following scenarios:
  /// (1) two fine boxes nested within the coarse box are touching each other;
  /// (2) on a periodic domain, one box is touching the low end of the domain
  ///      while another is touching the high end;
  /// (3) on a non-Periodic domain, a box is touching a domain boundary,
  ///      those adjCellBox's will be filled by BCHolder, however,
  ///      fill the corner ghosts outside domain when FillCornerGhosts=true;
  /// In summary, ghosts that should not be interpolated include
  ///  inter-Patch ghosts on same level and some outside-domain ghosts
  void initializeCoarsenedGhosts()
  {
    CH_TIME("CoarseFineInterp::initializeCoarsenedGhosts()");
    // the following is to initialize the set of coarsened Ghosts.
    m_coarsenedGhosts.define(m_layoutCoarsened);
    const Box domainBox = m_coarseDomain.domainBox();
    // periodicTestBox determines whether a box b adjoins a domain boundary
    //  by checking whether b is contained in it.
    // This is possible because of CoarsenedInterpRadius<NestingWidth
    Box periodicTestBox(domainBox);
    for (int d=0; d<SpaceDim; d++)
      periodicTestBox.grow(d, -m_int[CoarsenedInterpRadius]);

    // loop over each coarsened fine patch
    //  and determine the coarsened ghosts to be interpolated.
    for (DataIterator dit = m_layout.dataIterator(); dit.ok(); ++dit) {
      const Box& coarseBox = m_layoutCoarsened[dit];
      // Initially set coarsenedGhostsHere to the whole ghosted coarse box.
      const Box b = grow(coarseBox, m_int[CoarsenedInterpRadius]);
      m_coarsenedGhosts[dit] = IntVectSet(b);
      // remove corner ghosts of b if specified so
      if (!m_int[FillCornerGhosts])
        removeCornerGhosts(m_coarsenedGhosts[dit],
                           b, m_int[CoarsenedInterpRadius]);
      // We shouldn't worry about ghosts outside of a non-Periodic domain,
      // since they will be taken care of by the boundary conditions.
      for (int d=0; d<SpaceDim; d++)
        if (!m_coarseDomain.isPeriodic(d))
          {    
            for (SideIterator sit; sit.ok(); ++sit)
              if (coarseBox.sideEnd(sit())[d] == domainBox.sideEnd(sit())[d])
                m_coarsenedGhosts[dit] -=
                  adjCellBox(coarseBox, d, sit(), m_int[CoarsenedInterpRadius]);
          }

      // Now deal with inter-Patch ghosts.
      for (LayoutIterator lit=m_layoutCoarsened.layoutIterator();
           lit.ok(); ++lit)
        {
          const Box& coarseOtherBox = m_layoutCoarsened[lit];
          if (coarseBox==coarseOtherBox)
            {
              // if not the whole patch interp, subtract the interior
              if (m_int[InterpGhostsOnly])
                m_coarsenedGhosts[dit] -= coarseOtherBox;
            }
          // If the domain is periodic and 
          // coarseBox is touching the domain boundary,
          // some ghost cells might be set by others periodically.
          // subtract these possibilities off
          //  since they will be filled through exchange
          else if (m_coarseDomain.isPeriodic()
                   && !periodicTestBox.contains(coarseOtherBox)
                   && !periodicTestBox.contains(coarseBox))
            {
              ShiftIterator shiftIt = m_coarseDomain.shiftIterator();
              const IntVect shiftMult(m_coarseDomain.domainBox().size());
              Box shiftedBox(coarseOtherBox);
              for (shiftIt.begin(); shiftIt.ok(); ++shiftIt)
                {
                  const IntVect shiftVect = shiftMult*shiftIt();
                  shiftedBox.shift(shiftVect);
                  m_coarsenedGhosts[dit] -= shiftedBox;
                  shiftedBox.shift(-shiftVect);
                }
            }
        }
      m_coarsenedGhosts[dit] &= m_validSourceBox;
    }
  }

  /// assemble m_interpolators.
  void initializeInterpolators()
  {
    CH_TIME("CoarseFineInterp::initializeInterpolators()");
    m_interpMap.clear();
    const IntVect coarseGhostVect = m_int[CoarsenedInterpRadius] *IntVect::Unit;
    const IntVect fineGhostVect = m_int[InterpRadius] *IntVect::Unit;
    m_interpolators.define(m_layoutCoarsened, 1, coarseGhostVect);

    for (DataIterator dit = m_layoutCoarsened.dataIterator(); dit.ok(); ++dit)
      {
        const Box& smallCrseBox = m_layoutCoarsened[dit];
        const Box smallGhostedBox = grow(smallCrseBox, coarseGhostVect);
        const Box fineGhostedBox = grow(m_layout[dit], fineGhostVect);
        // Find the coarse level box that contains the coarsened fine box
        //  and serves as the valid source to choose stencil from.
        Box bigCrseBox = grow(smallCrseBox, m_int[NestingWidth]);
        bigCrseBox &= m_validSourceBox;
        const IntVectSet bigCrseArea(bigCrseBox);

        int degree = m_int[PolynomialDegree];
        int nCoefs = m_int[NumPolynmCoefs];
        const int nCells = bigCrseArea.numPts();
        // deal with the case that there aren't enough data for poly-fitting.
        if ( nCells<nCoefs )
          {
            if ( m_int[InsistOnAccuracy] )
              MayDay::Abort("Not enough cells for fitting polynomials");
            else
              pout() << "\n WARNING : Polynomial degree is reduced to ";
            // reduce the order
            while (nCells<nCoefs)
              nCoefs = numPolynmCoefs(--degree);
            pout() << degree << " \n  because there are only "
                   << bigCrseArea.numPts() << " valid cells in "
                   << bigCrseBox << std::endl;
          }
        //   pout()<< "ghosts to be interpolated:\n" << m_coarsenedGhosts[dit];

        // loop over the coarsened ghosts to cellwise initialize interpolators 
        for (IVSIterator ivsit(m_coarsenedGhosts[dit]); ivsit.ok(); ++ivsit)
          {
            const IntVect iv = ivsit();
            // figure out the fine box to be interpolated
            Box fineBase(iv, iv);
            fineBase.refine(m_int[RefineRatio]);
            fineBase &= fineGhostedBox;
            fineBase.shift( -scale(iv,m_int[RefineRatio]) );

            IntVectSet shiftedBigArea(bigCrseArea);
            // shift bigCrseArea so that the zero index points to ivsit()
            shiftedBigArea.shift(-iv);
            RefCountedPtr<Stencil> ps(new Stencil());
            ps->define(shiftedBigArea, fineBase, m_int[RefineRatio],
                       degree, nCoefs);
            const IntVectSet& stencil = ps->getStencil();
            //  pout()<< "Chosen stencil for "<< iv<< " is:\n"<< stencil<< endl;

            // initialize the associated interpolator
            if (m_interpMap.find(stencil)==m_interpMap.end())
              m_interpMap[stencil] = ps;
            // hook the pointer up
            m_interpolators[dit](ivsit()) = ps;
          }
      }
  }

  int factorial(const int n, const int k) const
  {
    int f = n;
    for (int i=1; i<k; i++)
      f *= n-i;
    return f;
  }

  int binomialCoefficient(const int n, int k) const
  {
    if (2*k>n) k = n-k;
    return factorial(n,k)/factorial(k,k);
  }

  int numPolynmCoefs(const int a_degreeOfPolynomial) const
  {
    // Putting PolynomialDegree indistinguishable balls
    //  into SpaceDim+1 distinguishable urns.
    return binomialCoefficient(a_degreeOfPolynomial+SpaceDim, SpaceDim);
  }

  /// Remove the ghosts cell in ghostedBox from ivSet
  /// ghostedBox must have at least 2*nGhosts for every dimension.
  void removeCornerGhosts(IntVectSet& ivSet, const Box& ghostedBox,
                          const int nGhosts, 
                          /// alsoRemoveEdges only makes a difference in 3D:
                          /// if true, the 12 edges of a cube are also removed.
                          const bool alsoRemoveEdges=true) const
  {
    if (alsoRemoveEdges)
      {
        const Box box = grow(ghostedBox, -nGhosts*IntVect::Unit);
        IntVectSet gCorners(ghostedBox);
        for (int i=0; i<SpaceDim; i++)
          gCorners -= grow(box, nGhosts*BASISV(i));
        ivSet -= gCorners;        
      }
    else
      {
        Box bLo(ghostedBox.smallEnd(), nGhosts-1+ ghostedBox.smallEnd());
        Box bHi(1-nGhosts+ ghostedBox.bigEnd(), ghostedBox.bigEnd());
        ivSet -= bLo;
        ivSet -= bHi;
        for (int i=0; i<SpaceDim; i++)
          {
            const int shift = ghostedBox.size(i) - nGhosts;
            bLo.shift(i, +shift);
            bHi.shift(i, -shift);
            ivSet -= bLo;
            ivSet -= bHi;
            bLo.shift(i, -shift);
            bHi.shift(i, +shift);
          }
      }
  }

  /// Disallowed for all the usual reasons
  void operator=(const CoarseFineInterp& a_input)
  {
    MayDay::Error("invalid operator");
  }

  /// Disallowed for all the usual reasons
  CoarseFineInterp(const CoarseFineInterp& a_input)
  {
    MayDay::Error("invalid operator");
  }
};

#endif
