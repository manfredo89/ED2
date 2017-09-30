module ed_nbg_init
  contains

!==========================================================================================!
!==========================================================================================!
!      This subroutine initializes a near-bare ground polygon.                             !
!------------------------------------------------------------------------------------------!
subroutine near_bare_ground_init(cgrid)
   use ed_state_vars  , only : edtype            & ! structure
                             , polygontype       & ! structure
                             , sitetype          & ! structure
                             , allocate_sitetype ! ! subroutine
   use ed_misc_coms   , only : ied_init_mode     ! ! intent(in)
   use physiology_coms, only : n_plant_lim       ! ! intent(in)
   use ed_type_init

   implicit none

   !----- Arguments. ----------------------------------------------------------------------!
   type(edtype)      , target  :: cgrid
   !----- Local variables. ----------------------------------------------------------------!
   type(polygontype) , pointer :: cpoly
   type(sitetype)    , pointer :: csite
   integer                     :: ipy
   integer                     :: isi
   !---------------------------------------------------------------------------------------!


   !----- Big loop ------------------------------------------------------------------------!
   do ipy=1,cgrid%npolygons
      cpoly => cgrid%polygon(ipy)
      do isi=1,cpoly%nsites
         csite => cpoly%site(isi)
         !----- We start with a single patch per site, always with primary vegetation. ----!
         csite%npatches = 1
         call allocate_sitetype(csite,1)
         csite%dist_type          (1) = 3
         csite%age                (1) = 0.0
         csite%area               (1) = 1.0


         !---------------------------------------------------------------------------------!
         !     Someone that uses the nitrogen model should check whether this is necessary !
         ! or not.  If nitrogen limitation is off, then we start all carbon and nitrogen   !
         ! pools with zeroes, otherwise we initialise with the former default values.      !
         !---------------------------------------------------------------------------------!
         select case (n_plant_lim)
         case (0)
            csite%fast_soil_C        (1) = 0.0
            csite%slow_soil_C        (1) = 0.0
            csite%structural_soil_C  (1) = 0.0
            csite%structural_soil_L  (1) = 0.0
            csite%mineralized_soil_N (1) = 0.0
            csite%fast_soil_N        (1) = 0.0

         case (1)
            csite%fast_soil_C        (1) = 0.2
            csite%slow_soil_C        (1) = 0.01
            csite%structural_soil_C  (1) = 10.0
            csite%structural_soil_L  (1) = csite%structural_soil_C (1)
            csite%mineralized_soil_N (1) = 1.0
            csite%fast_soil_N        (1) = 1.0

         end select
         !---------------------------------------------------------------------------------!

         csite%sum_dgd            (1) = 0.0
         csite%sum_chd            (1) = 0.0
         csite%plant_ag_biomass   (1) = 0.

         !----- We now populate the cohorts with near bare ground condition. --------------!
         select case (ied_init_mode)
         case (-1,0)
            call init_nbg_cohorts(csite,cpoly%lsl(isi),1,csite%npatches)
         end select

         !----- Initialise the patches now that cohorts are there. ------------------------!
         call init_ed_patch_vars(csite,1,csite%npatches,cpoly%lsl(isi))
      end do
      !----- Initialise some site-level variables. ----------------------------------------!
      call init_ed_site_vars(cpoly)
   end do

   !----- Last, but not the least, the polygons. ------------------------------------------!
   call init_ed_poly_vars(cgrid)

   return
end subroutine near_bare_ground_init
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!      This subroutine assigns a near-bare ground (NBG) state for some patches.            !
!------------------------------------------------------------------------------------------!
subroutine init_nbg_cohorts(csite,lsl,ipa_a,ipa_z)
   use ed_state_vars      , only : edtype             & ! structure
                                 , polygontype        & ! structure
                                 , sitetype           & ! structure
                                 , patchtype          & ! structure
                                 , allocate_sitetype  & ! subroutine
                                 , allocate_patchtype ! ! subroutine
   use ed_max_dims        , only : n_pft              ! ! intent(in)
   use ed_misc_coms       , only : ied_init_mode      ! ! intent(in)
   use pft_coms           , only : q                  & ! intent(in)
                                 , qsw                & ! intent(in)
                                 , sla                & ! intent(in)
                                 , hgt_min            & ! intent(in)
                                 , include_pft        & ! intent(in)
                                 , include_pft_ag     & ! intent(in)
                                 , include_pft_fp     & ! intent(in)
                                 , init_density       & ! intent(in)
                                 , agf_bs             ! ! intent(in)
   use consts_coms        , only : t3ple              & ! intent(in)
                                 , pio4               & ! intent(in)
                                 , kgom2_2_tonoha     & ! intent(in)
                                 , tonoha_2_kgom2     ! ! intent(in)
   use allometry          , only : h2dbh              & ! function
                                 , dbh2bd             & ! function
                                 , size2bl            & ! function
                                 , ed_biomass         & ! function
                                 , area_indices       ! ! subroutine
   use fuse_fiss_utils    , only : sort_cohorts    ! ! subroutine
   use ed_type_init

   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   type(sitetype)        , target     :: csite             ! Current site
   integer               , intent(in) :: lsl               ! Lowest soil level
   integer               , intent(in) :: ipa_a             ! 1st patch
   integer               , intent(in) :: ipa_z             ! Last patch
   !----- Local variables -----------------------------------------------------------------!
   type(patchtype)       , pointer    :: cpatch            ! Current patch
   integer                            :: ipa               ! Patch number
   integer                            :: ico               ! Cohort counter
   integer                            :: mypfts            ! Number of included PFTs
   integer                            :: ipft              ! PFT counter
   real                               :: salloc            ! balive/bleaf when on allom.
   real                               :: salloci           ! 1./salloc
   !---------------------------------------------------------------------------------------!

   !----- Patch loop. ---------------------------------------------------------------------!
   patchloop: do ipa=ipa_a,ipa_z
      cpatch => csite%patch(ipa)

      !----- Check which initialisation we are going to use. ------------------------------!
      select case (ied_init_mode)
      case (-1) !------ True bare ground simulation (absolute desert). --------------------!
         mypfts = 0
         return
      case ( 0) !------ Nearly bare ground simulation (start with a few seedlings). -------!

         !---------------------------------------------------------------------------------!
         !     Decide how many cohorts to allocate.                                        !
         !---------------------------------------------------------------------------------!
         select case (csite%dist_type(ipa))
         case (1)     !---- Agriculture. --------------------------------------------------!
            mypfts = count(include_pft_ag)
         case (2)     !---- Forest plantation. --------------------------------------------!
            mypfts = count(include_pft_fp)
         case default !---- Secondary or primary forest. ----------------------------------!
            mypfts = count(include_pft)
         end select
         !---------------------------------------------------------------------------------!
      end select

      !----- Perform cohort allocation. ---------------------------------------------------!
      call allocate_patchtype(cpatch,mypfts)

      !------------------------------------------------------------------------------------!
      !    Here we loop over PFTs rather than assigning the cohorts, so we ensure to only  !
      ! include the PFTs that should be  included (i.e., patches that should not happen in !
      ! agricultural patches will be skipped, and agricultural PFTs will not be included   !
      ! in the forest sites.                                                               !
      !------------------------------------------------------------------------------------!
      ico = 0
      pftloop: do ipft = 1,n_pft
         select case (csite%dist_type(ipa))
         case (1)
            !----- Agriculture, only the ones allowed are included. -----------------------!
            if (include_pft_ag(ipft)) then
               ico = ico + 1
            else
               cycle pftloop
            end if
         case (2)
            !----- Forest plantation, only the ones allowed are included. -----------------!
            if (include_pft_fp(ipft)) then
               ico = ico + 1
            else
               cycle pftloop
            end if
         case default
            !----- Primary or secondary vegetation, only the ones allowed are included. ---!
            if (include_pft(ipft)) then
               ico = ico + 1
            else
               cycle pftloop
            end if
         end select

         !----- The PFT is the plant functional type. -------------------------------------!
         cpatch%pft(ico)              = ipft
         !---------------------------------------------------------------------------------!
         !     Define the near-bare ground state using the standard minimum height and     !
         ! minimum plant density.  We assume all NBG PFTs to have leaves fully flushed,    !
         ! but with no storage biomass.  We then compute the other biomass quantities      !
         ! using the standard allometry for this PFT.                                      !
         !---------------------------------------------------------------------------------!
         cpatch%nplant(ico)           = init_density(ipft)
         cpatch%hite(ico)             = hgt_min(ipft)
         cpatch%phenology_status(ico) = 0
         cpatch%dbh(ico)              = h2dbh(cpatch%hite(ico),ipft)
         cpatch%bdead(ico)            = dbh2bd(cpatch%dbh(ico),ipft)
         cpatch%bleaf(ico)            = size2bl(cpatch%dbh(ico),cpatch%hite(ico),ipft)
         cpatch%sla(ico)              = sla(ipft)


         salloc                       = 1.0 + q(ipft) + qsw(ipft) * cpatch%hite(ico)
         salloci                      = 1. / salloc

         cpatch%balive(ico)           = cpatch%bleaf(ico) * salloc
         cpatch%broot(ico)            = q(ipft) * cpatch%balive(ico) * salloci
         cpatch%bsapwooda(ico)        = qsw(ipft) * cpatch%hite(ico) * cpatch%balive(ico)  &
                                      * salloci * agf_bs(ipft)
         cpatch%bsapwoodb(ico)        = qsw(ipft) * cpatch%hite(ico) * cpatch%balive(ico)  &
                                      * salloci * (1.-agf_bs(ipft))
         cpatch%bstorage(ico)         = 0.5 * ( cpatch%bleaf(ico)                          &
                                              + cpatch%broot(ico)                          &
                                              + cpatch%bsapwooda(ico)                      &
                                              + cpatch%bsapwoodb(ico))

         !----- Find the initial area indices (LAI, WAI, CAI). ----------------------------!
         call area_indices(cpatch, ico)

         !----- Find the above-ground biomass and basal area. -----------------------------!
         cpatch%agb(ico) = ed_biomass(cpatch, ico)
         cpatch%basarea(ico) = pio4 * cpatch%dbh(ico)*cpatch%dbh(ico)

         !----- Initialize other cohort-level variables. ----------------------------------!
         call init_ed_cohort_vars(cpatch,ico,lsl)

         !----- Update total patch-level above-ground biomass -----------------------------!
         csite%plant_ag_biomass(ipa) = csite%plant_ag_biomass(ipa)                         &
                                     + cpatch%nplant(ico) * cpatch%agb(ico)
      end do pftloop

      !------------------------------------------------------------------------------------!
      !     Since initial heights may not be constant, we must sort the cohorts.           !
      !------------------------------------------------------------------------------------!
      call sort_cohorts(cpatch)
   end do patchloop

   return
end subroutine init_nbg_cohorts
!==========================================================================================!
!==========================================================================================!

end module ed_nbg_init
