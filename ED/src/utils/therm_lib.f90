!==========================================================================================!
! BRAMS-4.0.6. File: therm_lib.f90                                                                            !
!                                                                                          !
!    This file contains most functions and subroutines that deal with several thermo-      !
! dynamic conversions. These procedures were built to avoid assumptions like hydrostatic   !
! and linearisation. Most equations could not be solved analytically, and the standard     !
! here was to use Newton's method as the default, always having bisection or, more often,  !
! the modified Regula Falsi (Illinois) method in case Newton's fails.                      !
!==========================================================================================!
!==========================================================================================!
module therm_lib
   implicit none

   !---------------------------------------------------------------------------------------!
   ! Constants that control the convergence for iterative methods                          !
   !---------------------------------------------------------------------------------------!
   real(kind=4), parameter ::   toler  = 10.* epsilon(1.) ! Relative tolerance for iter-
                                                          !    ative methods.  The smaller
                                                          !    the value, the more accurate
                                                          !    the result, but smaller
                                                          !    values will slow down the
                                                          !    run.
   integer     , parameter ::   maxfpo = 60               ! Maximum # of iterations before
                                                          !    crashing for false position
                                                          !    method.
   integer     , parameter ::   maxit  = 150              ! Maximum # of iterations before
                                                          !    crashing, for other methods.
   integer     , parameter ::   maxlev = 16               ! Maximum # of levels for adap-
                                                          !    tive quadrature methods.
   logical     , parameter ::   newthermo = .true.        ! Use new thermodynamics [T|F]
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !   This is the "level" variable, that used to be in micphys. Since it affects more the !
   ! thermodynamics choices than the microphysics, it was moved to here.                   !
   !---------------------------------------------------------------------------------------!
   integer, parameter ::   level = 3
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !    The following three variables are just the logical tests on variable "level",      !
   ! saved here to speed up checks for "li" functions.                                     !
   !---------------------------------------------------------------------------------------!
   logical, parameter ::   vapour_on = .true.
   logical, parameter ::   cloud_on  = .true.
   logical, parameter ::   bulk_on   = .true.
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     These constants came from the paper in which the saturation vapour pressure is    !
   ! based on:                                                                             !
   !                                                                                       !
   !  Murphy, D. M.; Koop, T., 2005: Review of the vapour pressures of ice and supercooled !
   !     water for atmospheric applications. Q. J. Royal Meteor. Soc., vol. 31, pp. 1539-  !
   !     1565 (hereafter MK05).                                                            !
   !                                                                                       !
   !  These equations give the triple point at t3ple, with vapour pressure being es3ple.   !
   !---------------------------------------------------------------------------------------!
   !----- Coefficients based on equation (7): ---------------------------------------------!
   real(kind=4), dimension(0:3), parameter :: iii_7   = (/  9.550426, -5723.265            &
                                                         ,  3.530680,    -0.00728332 /)
   !----- Coefficients based on equation (10), first fit ----------------------------------!
   real(kind=4), dimension(0:3), parameter :: l01_10  = (/ 54.842763, -6763.220            &
                                                         , -4.210   ,     0.000367   /)
   !----- Coefficients based on equation (10), second fit ---------------------------------!
   real(kind=4), dimension(0:3), parameter :: l02_10  = (/ 53.878   , -1331.22             &
                                                         , -9.44523 , 0.014025       /)
   !----- Coefficients based on the hyperbolic tangent ------------------------------------!
   real(kind=4), dimension(2)  , parameter :: ttt_10  = (/  0.0415  ,   218.80       /)
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     These constants came from the paper in which the saturation vapour pressure is    !
   ! based on:                                                                             !
   !                                                                                       !
   !  Flatau, P. J.; Walko, R. L.; Cotton, W. R., 1992: Polynomial fits to saturation      !
   !     vapor pressure. J. Appl. Meteor., vol. 31, pp. 1507-1513. (hereafter FWC92).      !
   !                                                                                       !
   !  These equations give the triple point at 273.004K.                                   !
   !  N.B.: The coefficients here don't seem to match those listed on FWC92, but that's    !
   !        what was on the original code...                                               !
   !---------------------------------------------------------------------------------------!
   !----- Coefficients for esat (liquid) --------------------------------------------------!
   real(kind=4), dimension(0:8), parameter :: cll = (/  .6105851e+03,  .4440316e+02        &
                                                     ,  .1430341e+01,  .2641412e-01        &
                                                     ,  .2995057e-03,  .2031998e-05        &
                                                     ,  .6936113e-08,  .2564861e-11        &
                                                     , -.3704404e-13                /)
   !----- Coefficients for esat (ice) -----------------------------------------------------!
   real(kind=4), dimension(0:8), parameter :: cii = (/  .6114327e+03,  .5027041e+02        &
                                                     ,  .1875982e+01,  .4158303e-01        &
                                                     ,  .5992408e-03,  .5743775e-05        &
                                                     ,  .3566847e-07,  .1306802e-09        &
                                                     ,  .2152144e-12                /)
   !----- Coefficients for d(esat)/dT (liquid) --------------------------------------------!
   real(kind=4), dimension(0:8), parameter :: dll = (/  .4443216e+02,  .2861503e+01        &
                                                     ,  .7943347e-01,  .1209650e-02        &
                                                     ,  .1036937e-04,  .4058663e-07        &
                                                     , -.5805342e-10, -.1159088e-11        &
                                                     , -.3189651e-14                /)
   !----- Coefficients for esat (ice) -----------------------------------------------------!
   real(kind=4), dimension(0:8), parameter :: dii = (/  .5036342e+02,  .3775758e+01        &
                                                     ,  .1269736e+00,  .2503052e-02        &
                                                     ,  .3163761e-04,  .2623881e-06        &
                                                     ,  .1392546e-08,  .4315126e-11        &
                                                     ,  .5961476e-14                /)
   !=======================================================================================!
   !=======================================================================================!


   contains



   !=======================================================================================!
   !=======================================================================================!
   !     This function calculates the liquid saturation vapour pressure as a function of   !
   ! Kelvin temperature. This expression came from MK05, equation (10).                    !
   !---------------------------------------------------------------------------------------!
   real(kind=4) function eslf(temp,l1funout,l2funout,ttfunout)
      use consts_coms, only : t00 ! ! intent(in)
      implicit none
      !----- Required arguments. ----------------------------------------------------------!
      real(kind=4), intent(in)            :: temp     ! Temperature                [     K]
      !----- Optional arguments. ----------------------------------------------------------!
      real(kind=4), intent(out), optional :: l1funout ! Function for high temperatures
      real(kind=4), intent(out), optional :: ttfunout ! Interpolation function
      real(kind=4), intent(out), optional :: l2funout ! Function for low temperatures
      !----- Local variables. -------------------------------------------------------------!
      real(kind=4)                        :: l1fun    !
      real(kind=4)                        :: ttfun    !
      real(kind=4)                        :: l2fun    !
      real(kind=4)                        :: x        !
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Choose between the old and the new thermodynamics.                             !
      !------------------------------------------------------------------------------------!
      if (newthermo) then
         !----- Updated method, using MK05 ------------------------------------------------!
         l1fun = l01_10(0) + l01_10(1)/temp + l01_10(2)*log(temp) + l01_10(3) * temp
         l2fun = l02_10(0) + l02_10(1)/temp + l02_10(2)*log(temp) + l02_10(3) * temp
         ttfun = tanh(ttt_10(1) * (temp - ttt_10(2)))
         eslf  = exp(l1fun + ttfun*l2fun)
         !---------------------------------------------------------------------------------!

         if (present(l1funout)) l1funout = l1fun
         if (present(l2funout)) l2funout = l2fun
         if (present(ttfunout)) ttfunout = ttfun
      else
         !----- Original method, using polynomial fit (FWC92) -----------------------------!
         x    = max(-80.,temp-t00)
         eslf = cll(0) + x * (cll(1) + x * (cll(2) + x * (cll(3) + x * (cll(4)             &
                       + x * (cll(5) + x * (cll(6) + x * (cll(7) + x * cll(8)) ) ) ) ) ) )
         !---------------------------------------------------------------------------------!

         if (present(l1funout)) l1funout = eslf
         if (present(l2funout)) l2funout = eslf
         if (present(ttfunout)) ttfunout = eslf
      end if

      return
   end function eslf
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This function calculates the ice saturation vapour pressure as a function of      !
   ! Kelvin temperature, based on MK05 equation (7).                                       !
   !---------------------------------------------------------------------------------------!
   real(kind=4) function esif(temp,iifunout)
      use consts_coms, only : t00 ! ! intent(in)
      implicit none
      !----- Required arguments. ----------------------------------------------------------!
      real(kind=4), intent(in)            :: temp     ! Temperature                 [    K]
      !----- Optional arguments. ----------------------------------------------------------!
      real(kind=4), intent(out), optional :: iifunout
      !----- Local variables. -------------------------------------------------------------!
      real(kind=4)                        :: iifun
      real(kind=4)                        :: x
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Choose between the old and the new thermodynamics.                             !
      !------------------------------------------------------------------------------------!
      if (newthermo) then

         !----- Updated method, using MK05 ------------------------------------------------!
         iifun = iii_7(0) + iii_7(1)/temp + iii_7(2) * log(temp) + iii_7(3) * temp
         esif  = exp(iifun)
         !---------------------------------------------------------------------------------!


         if (present(iifunout)) iifunout=iifun
      else
         !----- Original method, using polynomial fit (FWC92) -----------------------------!
         x=max(-80.,temp-t00)
         esif = cii(0) + x * (cii(1) + x * (cii(2) + x * (cii(3) + x * (cii(4)             &
                       + x * (cii(5) + x * (cii(6) + x * (cii(7) + x * cii(8)) ) ) ) ) ) )
         !---------------------------------------------------------------------------------!

         if (present(iifunout)) iifunout=esif
      end if
      !------------------------------------------------------------------------------------!

      return
   end function esif
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This function calculates the saturation vapour pressure as a function of  Kelvin  !
   ! temperature. It chooses which phase to look depending on whether the temperature is   !
   ! below or above the triple point.                                                      !
   !---------------------------------------------------------------------------------------!
   real(kind=4) function eslif(temp,useice)
      use consts_coms, only : t3ple ! ! intent(in)
      implicit none
      !----- Required arguments. ----------------------------------------------------------!
      real(kind=4), intent(in)           :: temp
      !----- Optional arguments. ----------------------------------------------------------!
      logical     , intent(in), optional :: useice
      !----- Local variables. -------------------------------------------------------------!
      logical                            :: frozen
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Decide which function to use (saturation for liquid water or ice).             !
      !------------------------------------------------------------------------------------!
      if (present(useice)) then
         frozen = useice  .and. temp < t3ple
      else
         frozen = bulk_on .and. temp < t3ple
      end if
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !      Call the appropriate function depending on the temperature and whether ice    !
      ! thermodynamics is to be used.                                                      !
      !------------------------------------------------------------------------------------!
      if (frozen) then
         !----- Saturation vapour pressure for ice. ---------------------------------------!
         eslif = esif(temp)
         !---------------------------------------------------------------------------------!
      else
         !----- Saturation vapour pressure for liquid. ------------------------------------!
         eslif = eslf(temp)
         !---------------------------------------------------------------------------------!
      end if
      !------------------------------------------------------------------------------------!

      return
   end function eslif
   !=======================================================================================!
   !=======================================================================================!




   !=======================================================================================!
   !=======================================================================================!
   !     This function calculates the saturation specific humidity, over liquid or ice     !
   ! depending on temperature, as a function of pressure and Kelvin temperature.           !
   !---------------------------------------------------------------------------------------!
   real(kind=4) function qslif(pres,temp,useice)
      use consts_coms, only : t3ple  & ! intent(in)
                            , ep     & ! intent(in)
                            , toodry ! ! intent(in)
      implicit none
      !----- Required arguments. ----------------------------------------------------------!
      real(kind=4), intent(in)           :: pres
      real(kind=4), intent(in)           :: temp
      !----- Optional arguments. ----------------------------------------------------------!
      logical     , intent(in), optional :: useice
      !----- Local variables. -------------------------------------------------------------!
      real(kind=4)                       :: esz
      logical                            :: frozen
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Decide which function to use (saturation for liquid water or ice).             !
      !------------------------------------------------------------------------------------!
      if (present(useice)) then
         frozen = useice  .and. temp < t3ple
      else
         frozen = bulk_on .and. temp < t3ple
      end if
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !      Call the appropriate function depending on the temperature and whether ice    !
      ! thermodynamics is to be used.                                                      !
      !------------------------------------------------------------------------------------!
      if (frozen) then
         !----- Saturation vapour pressure for ice. ---------------------------------------!
         esz = esif(temp)
         !---------------------------------------------------------------------------------!
      else
         !----- Saturation vapour pressure for liquid. ------------------------------------!
         esz = eslf(temp)
         !---------------------------------------------------------------------------------!
      end if
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Use the usual relation between pressure and vapour pressure to find the        !
      ! specific humidity.                                                                 !
      !------------------------------------------------------------------------------------!
      qslif = max(toodry, ep * esz/( pres - (1.0 - ep) * esz) )
      !------------------------------------------------------------------------------------!

      return
   end function qslif
   !=======================================================================================!
   !=======================================================================================!







   !=======================================================================================!
   !=======================================================================================!
   !     This function calculates the vapour-liquid equilibrium density for vapour, as a   !
   ! function of temperature in Kelvin.                                                    !
   !---------------------------------------------------------------------------------------!
   real(kind=4) function rhovsl(temp)
      use consts_coms, only : rh2o ! ! intent(in)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      real(kind=4), intent(in) :: temp ! Temperature                                [    K]
      !----- Local variables. -------------------------------------------------------------!
      real(kind=4)             :: eequ ! Saturation vapour pressure                 [   Pa]
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Find the equilibrium (saturation) vapour pressure.                             !
      !------------------------------------------------------------------------------------!
      eequ = eslf(temp)
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !    Find the saturation density.                                                    !
      !------------------------------------------------------------------------------------!
      rhovsl = eequ / (rh2o * temp)
      !------------------------------------------------------------------------------------!

      return
   end function rhovsl
   !=======================================================================================!
   !=======================================================================================!




   !=======================================================================================!
   !=======================================================================================!
   !     This function calculates the partial derivative of ice saturation vapour pressure !
   ! with respect to temperature as a function of Kelvin temperature.                      !
   !---------------------------------------------------------------------------------------!
   real(kind=4) function esifp(temp)
      use consts_coms, only : t00 ! ! intent(in)
      implicit none
      !------ Arguments. ------------------------------------------------------------------!
      real(kind=4), intent(in) :: temp
      !------ Local variables. ------------------------------------------------------------!
      real(kind=4)             :: esi
      real(kind=4)             :: iiprime
      real(kind=4)             :: x
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !      Decide which function to use, based on the thermodynamics.                    !
      !------------------------------------------------------------------------------------!
      if (newthermo) then
         !----- Updated method, using MK05 ------------------------------------------------!
         esi     = esif(temp)
         iiprime = -iii_7(1)/(temp*temp) + iii_7(2)/temp + iii_7(3)
         esifp   = esi * iiprime
      else
         !----- Original method, using polynomial fit (FWC92) -----------------------------!
         x=max(-80.,temp-t00)
         esifp = dii(0) + x * (dii(1) + x * (dii(2) + x * (dii(3) + x * (dii(4)            &
                        + x * (dii(5) + x * (dii(6) + x * (dii(7) + x * dii(8)) ) ) ) ) ) )
      end if
      !------------------------------------------------------------------------------------!

      return
   end function esifp
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This function calculates the partial derivative of liquid saturation vapour       !
   ! pressure with respect to temperature as a function of Kelvin temperature.             !
   !---------------------------------------------------------------------------------------!
   real(kind=4) function eslfp(temp)
      use consts_coms, only : t00 ! ! intent(in)
      implicit none
      !------ Arguments. ------------------------------------------------------------------!
      real(kind=4), intent(in) :: temp
      !------ Local variables. ------------------------------------------------------------!
      real(kind=4)             :: esl
      real(kind=4)             :: l2fun
      real(kind=4)             :: ttfun
      real(kind=4)             :: l1prime
      real(kind=4)             :: l2prime
      real(kind=4)             :: ttprime
      real(kind=4)             :: x
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !      Decide which function to use, based on the thermodynamics.                    !
      !------------------------------------------------------------------------------------!
      if (newthermo) then
         !----- Updated method, using MK05 ------------------------------------------------!
         esl     = eslf(temp,l2funout=l2fun,ttfunout=ttfun)
         l1prime = -l01_10(1)/(temp*temp) + l01_10(2)/temp + l01_10(3)
         l2prime = -l02_10(1)/(temp*temp) + l02_10(2)/temp + l02_10(3)
         ttprime =  ttt_10(1)*(1.-ttfun*ttfun)
         eslfp   = esl * (l1prime + l2prime*ttfun + l2fun*ttprime)
      else
         !----- Original method, using polynomial fit (FWC92) -----------------------------!
         x     = max(-80.,temp-t00)
         eslfp = dll(0) + x * (dll(1) + x * (dll(2) + x * (dll(3) + x * (dll(4)            &
                        + x * (dll(5) + x * (dll(6) + x * (dll(7) + x * dll(8)) ) ) ) ) ) )
      end if
      !------------------------------------------------------------------------------------!


      return
   end function eslfp
   !=======================================================================================!
   !=======================================================================================!





   !=======================================================================================!
   !=======================================================================================!
   !     This function calculates the partial derivative of saturation vapour pressure as  !
   ! a function of  Kelvin temperature. It chooses which phase to look depending on        !
   ! whether the temperature is below or above the triple point.                           !
   !---------------------------------------------------------------------------------------!
   real(kind=4) function eslifp(temp,useice)
      use consts_coms, only : t3ple ! ! intent(in)
      implicit none
      !------ Arguments. ------------------------------------------------------------------!
      real(kind=4), intent(in)           :: temp
      !------ Local variables. ------------------------------------------------------------!
      logical     , intent(in), optional :: useice
      logical                            :: frozen
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Decide which function to use (saturation for liquid water or ice).             !
      !------------------------------------------------------------------------------------!
      if (present(useice)) then
         frozen = useice  .and. temp < t3ple
      else
         frozen = bulk_on .and. temp < t3ple
      end if
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !      Call the appropriate function depending on the temperature and whether ice    !
      ! thermodynamics is to be used.                                                      !
      !------------------------------------------------------------------------------------!
      if (frozen) then
         !----- d(Saturation vapour pressure)/dT for ice. ---------------------------------!
         eslifp = esifp(temp)
         !---------------------------------------------------------------------------------!
      else
         !----- d(Saturation vapour pressure)/dT for liquid water. ------------------------!
         eslifp = eslfp(temp)
         !---------------------------------------------------------------------------------!
      end if
      !------------------------------------------------------------------------------------!

      return
   end function eslifp
   !=======================================================================================!
   !=======================================================================================!


   !=======================================================================================!
   !=======================================================================================!
   !     This function computes the vapour mixing ratio based (or specific humidity) based !
   ! on the pressure [Pa], temperature [K] and relative humidity [fraction].  It checks    !
   ! the temperature to decide between ice or liquid saturation.                           !
   !---------------------------------------------------------------------------------------!
   real(kind=4) function ptrh2rvapil(relh,pres,temp,out_shv,useice)
      use consts_coms, only : ep      & ! intent(in)
                            , toodry  & ! intent(in)
                            , t3ple   ! ! intent(in)

      implicit none
      !----- Required arguments. ----------------------------------------------------------!
      real(kind=4), intent(in)           :: relh    ! Relative humidity            [    --]
      real(kind=4), intent(in)           :: pres    ! Pressure                     [    Pa]
      real(kind=4), intent(in)           :: temp    ! Temperature                  [     K]
      logical     , intent(in)           :: out_shv ! Output is specific humidity  [   T|F]
      !----- Optional arguments. ----------------------------------------------------------!
      logical     , intent(in), optional :: useice  ! May use ice thermodynamics   [   T|F]
      !----- Local variables. -------------------------------------------------------------!
      real(kind=4)                       :: pvap    ! Vapour pressure              [    Pa]
      real(kind=4)                       :: relhh   ! Bounded relative humidity    [    --]
      logical                            :: frozen  ! Will use ice thermodynamics  [   T|F]
      !------------------------------------------------------------------------------------!


      !----- Check whether to use the user's or the default flag for ice saturation. ------!
      if (present(useice)) then
         frozen = useice  .and. temp < t3ple
      else
         frozen = bulk_on .and. temp < t3ple
      end if
      !------------------------------------------------------------------------------------!



      !---- Make sure relative humidity is bounded. ---------------------------------------!
      relhh = min(1.,max(0.,relh))
      !------------------------------------------------------------------------------------!


      !---- Find the vapour pressure (ice or liquid, depending on the value of frozen). ---!
      if (frozen) then
         pvap  = relhh * esif(temp)
      else
         pvap  = relhh * eslf(temp)
      end if
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Convert the output to the sought humidity variable.                            !
      !------------------------------------------------------------------------------------!
      if (out_shv) then
         !----- Specific humidity. --------------------------------------------------------!
         ptrh2rvapil = max(toodry, ep * pvap / (pres - (1.0 - ep) * pvap))
         !---------------------------------------------------------------------------------!
      else
         !----- Mixing ratio. -------------------------------------------------------------!
         ptrh2rvapil = max(toodry, ep * pvap / (pres - pvap))
         !---------------------------------------------------------------------------------!
      end if
      !------------------------------------------------------------------------------------!
      return
   end function ptrh2rvapil
   !=======================================================================================!
   !=======================================================================================!





   !=======================================================================================!
   !=======================================================================================!
   !     This function computes the relative humidity [fraction] based on pressure, tem-   !
   ! perature, and vapour mixing ratio (or specific humidity). Two important points:       !
   ! 1. It may consider whether the temperature is above or below the freezing point       !
   !    to choose which saturation to use. It is possible to explicitly force not to use   !
   !    ice in case level is 2 or if you have reasons not to use ice (e.g. reading data    !
   !    that did not consider ice).
   ! 2. IT DOESN'T PREVENT SUPERSATURATION TO OCCUR. This is because this subroutine is    !
   !    also used in the microphysics, where supersaturation does happen and needs to be   !
   !    accounted.                                                                         !
   !---------------------------------------------------------------------------------------!
   real(kind=4) function rehuil(pres,temp,humi,is_shv,useice)
      use consts_coms, only : t3ple  & ! intent(in)
                            , ep     & ! intent(in)
                            , toodry ! ! intent(in)

      implicit none
      !----- Required arguments. ----------------------------------------------------------!
      real(kind=4), intent(in)           :: pres    ! Air pressure                 [    Pa]
      real(kind=4), intent(in)           :: temp    ! Temperature                  [     K]
      real(kind=4), intent(in)           :: humi    ! Humidity                     [ kg/kg]
      logical     , intent(in)           :: is_shv  ! Input is specific humidity   [   T|F]
       !----- Optional arguments. ----------------------------------------------------------!
      logical     , intent(in), optional :: useice  ! May use ice thermodynamics   [   T|F]
      !----- Local variables --------------------------------------------------------------!
      real(kind=4)                       :: shv     ! Specific humidity            [ kg/kg]
      real(kind=4)                       :: pvap    ! Vapour pressure              [    Pa]
      real(kind=4)                       :: psat    ! Saturation vapour pressure   [    Pa]
      logical                            :: frozen  ! Will use ice saturation now  [   T|F]
      !------------------------------------------------------------------------------------!

      !------------------------------------------------------------------------------------!
      !    Check whether we should use ice or liquid saturation.                           !
      !------------------------------------------------------------------------------------!
      if (present(useice)) then
         frozen = useice  .and. temp < t3ple
      else
         frozen = bulk_on .and. temp < t3ple
      end if
      !------------------------------------------------------------------------------------!


      !---- Make sure that we have specific humidity. -------------------------------------!
      if (is_shv) then
         shv = max(toodry,humi)
      else
         shv = max(toodry,humi) / ( 1.0 + max(toodry,humi) )
      end if
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Find the vapour pressure and the saturation vapour pressure.                   !
      !------------------------------------------------------------------------------------!
      pvap = ( pres * shv ) / ( ep + (1.0 - ep) * shv )
      if (frozen) then
         psat = esif (temp)
      else
         psat = esif (temp)
      end if
      !------------------------------------------------------------------------------------!

      !------------------------------------------------------------------------------------!
      !     Find the relative humidity.                                                    !
      !------------------------------------------------------------------------------------!
      rehuil = max(0. ,pvap / psat)
      !------------------------------------------------------------------------------------!

      return
   end function rehuil
   !=======================================================================================!
   !=======================================================================================!




   !=======================================================================================!
   !=======================================================================================!
   !     This function computes the vapour pressure deficit based on pressure, temper-     !
   ! ature, and vapour mixing ratio (or specific humidity).                                !
   !                                                                                       !
   ! IMPORTANT: This fucntion may consider whether the temperature is above or below the   !
   !            freezing point to choose which saturation to use. It is possible to        !
   !            explicitly force not to use ice in case level is 2 or if you have reasons  !
   !            not to use ice (e.g. reading data that did not consider ice).              !
   !---------------------------------------------------------------------------------------!
   real(kind=4) function vpdefil(pres,temp,humi,is_shv,useice)
      use consts_coms, only : t3ple  & ! intent(in)
                            , ep     & ! intent(in)
                            , toodry ! ! intent(in)
      implicit none
      !----- Required arguments. ----------------------------------------------------------!
      real(kind=4), intent(in)           :: pres    ! Air pressure                 [    Pa]
      real(kind=4), intent(in)           :: temp    ! Temperature                  [     K]
      real(kind=4), intent(in)           :: humi    ! Humidity                     [ kg/kg]
      logical     , intent(in)           :: is_shv  ! Input is specific humidity   [   T|F]
       !----- Optional arguments. ----------------------------------------------------------!
      logical     , intent(in), optional :: useice  ! May use ice thermodynamics   [   T|F]
      !----- Local variables --------------------------------------------------------------!
      real(kind=4)                       :: shv     ! Specific humidity            [ kg/kg]
      real(kind=4)                       :: pvap    ! Vapour pressure              [    Pa]
      real(kind=4)                       :: psat    ! Saturation vapour pressure   [    Pa]
      logical                            :: frozen  ! Will use ice saturation now  [   T|F]
      !------------------------------------------------------------------------------------!

      !------------------------------------------------------------------------------------!
      !    Check whether we should use ice or liquid saturation.                           !
      !------------------------------------------------------------------------------------!
      if (present(useice)) then
         frozen = useice  .and. temp < t3ple
      else
         frozen = bulk_on .and. temp < t3ple
      end if
      !------------------------------------------------------------------------------------!


      !---- Make sure that we have specific humidity. -------------------------------------!
      if (is_shv) then
         shv = max(toodry,humi)
      else
         shv = max(toodry,humi) / ( 1.0 + max(toodry,humi) )
      end if
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Find the vapour pressure and the saturation vapour pressure.                   !
      !------------------------------------------------------------------------------------!
      pvap = ( pres * shv ) / ( ep + (1.0 - ep) * shv )
      if (frozen) then
         psat = esif(temp)
      else
         psat = esif(temp)
      end if
      !------------------------------------------------------------------------------------!

      !------------------------------------------------------------------------------------!
      !     Find the relative humidity.                                                    !
      !------------------------------------------------------------------------------------!
      vpdefil = max(0.0 , psat - pvap)
      !------------------------------------------------------------------------------------!

      return
   end function vpdefil
   !=======================================================================================!
   !=======================================================================================!



   !=======================================================================================!
   !=======================================================================================!
   !     This function finds the density based on the virtual temperature and the ideal    !
   ! gas law.  The only difference between this function and the one above is that here we !
   ! provide vapour and total specific mass (specific humidity) instead of mixing ratio.   !
   !---------------------------------------------------------------------------------------!
   real(kind=4) function idealdenssh(pres,temp,qvpr,qtot)
      use consts_coms, only : rdry & ! intent(in)
                            , epi  ! ! intent(in)
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real(kind=4), intent(in)           :: pres ! Pressure                        [    Pa]
      real(kind=4), intent(in)           :: temp ! Temperature                     [     K]
      real(kind=4), intent(in)           :: qvpr ! Vapour specific mass            [ kg/kg]
      real(kind=4), intent(in), optional :: qtot ! Total water specific mass       [ kg/kg]
      !----- Local variables. -------------------------------------------------------------!
      real(kind=4)                       :: qall ! Either qtot or qvpr...          [ kg/kg]
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !      Prefer using total specific humidity, but if it isn't provided, then use      !
      ! vapour phase as the total (no condensation).                                       !
      !------------------------------------------------------------------------------------!
      if (present(qtot)) then
        qall = qtot
      else
        qall = qvpr
      end if
      !------------------------------------------------------------------------------------!


      !----- Convert using a generalised function. ----------------------------------------!
      idealdenssh = pres / (rdry * temp * (1. - qall + epi * qvpr))
      !------------------------------------------------------------------------------------!

      return
   end function idealdenssh
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This function computes reduces the pressure from the reference height to the      !
   ! canopy height by assuming hydrostatic equilibrium.  For simplicity, we assume that    !
   ! R and cp are constants (in reality they are dependent on humidity).                   !
   !---------------------------------------------------------------------------------------!
   real(kind=4) function reducedpress(pres,thetaref,shvref,zref,thetacan,shvcan,zcan)
      use consts_coms, only : epim1    & ! intent(in)
                            , p00k     & ! intent(in)
                            , rocp     & ! intent(in)
                            , cpor     & ! intent(in)
                            , cpdry    & ! intent(in)
                            , grav     ! ! intent(in)
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real(kind=4), intent(in) :: pres     ! Pressure                            [      Pa]
      real(kind=4), intent(in) :: thetaref ! Potential temperature               [       K]
      real(kind=4), intent(in) :: shvref   ! Vapour specific mass                [   kg/kg]
      real(kind=4), intent(in) :: zref     ! Height at reference level           [       m]
      real(kind=4), intent(in) :: thetacan ! Potential temperature               [       K]
      real(kind=4), intent(in) :: shvcan   ! Vapour specific mass                [   kg/kg]
      real(kind=4), intent(in) :: zcan     ! Height at canopy level              [       m]
      !------Local variables. -------------------------------------------------------------!
      real(kind=4)             :: pinc     ! Pressure increment                  [ Pa^R/cp]
      real(kind=4)             :: thvbar   ! Average virtual pot. temperature    [       K]
      !------------------------------------------------------------------------------------!

      !------------------------------------------------------------------------------------!
      !      First we compute the average virtual potential temperature between the canopy !
      ! top and the reference level.                                                       !
      !------------------------------------------------------------------------------------!
      thvbar = 0.5 * (thetaref * (1. + epim1 * shvref) + thetacan * (1. + epim1 * shvcan))
      !------------------------------------------------------------------------------------!



      !----- Then, we find the pressure gradient scale. -----------------------------------!
      pinc = grav * p00k * (zref - zcan) / (cpdry * thvbar)
      !------------------------------------------------------------------------------------!



      !----- And we can find the reduced pressure. ----------------------------------------!
      reducedpress = (pres**rocp + pinc ) ** cpor
      !------------------------------------------------------------------------------------!

      return
   end function reducedpress
   !=======================================================================================!
   !=======================================================================================!







   !=======================================================================================!
   !=======================================================================================!
   !     This function computes the Exner function [J/kg/K], given the pressure.  It       !
   ! assumes for simplicity that R and Cp are constants and equal to the dry air values.   !
   !---------------------------------------------------------------------------------------!
   real(kind=4) function press2exner(pres)
      use consts_coms, only : p00i           & ! intent(in)
                            , cpdry          & ! intent(in)
                            , rocp           ! ! intent(in)

      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      real(kind=4), intent(in) :: pres   ! Pressure                               [     Pa]
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Find potential temperature.                                                    !
      !------------------------------------------------------------------------------------!
      press2exner = cpdry * ( pres * p00i ) ** rocp
      !------------------------------------------------------------------------------------!

      return
   end function press2exner
   !=======================================================================================!
   !=======================================================================================!


   !=======================================================================================!
   !=======================================================================================!
   !     This function computes the potential temperature [K], given the Exner function    !
   ! and temperature.  For simplicity we ignore the effects of humidity in R and cp and    !
   ! use the dry air values instead.                                                       !
   !---------------------------------------------------------------------------------------!
   real(kind=4) function extemp2theta(exner,temp)
      use consts_coms, only : cpdry          ! ! intent(in)

      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      real(kind=4), intent(in) :: exner  ! Exner function                         [ J/kg/K]
      real(kind=4), intent(in) :: temp   ! Temperature                            [      K]
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Find potential temperature.                                                    !
      !------------------------------------------------------------------------------------!
      extemp2theta = cpdry * temp / exner
      !------------------------------------------------------------------------------------!

      return
   end function extemp2theta
   !=======================================================================================!
   !=======================================================================================!







   !=======================================================================================!
   !=======================================================================================!
   !     This function computes the temperature [K], given the Exner function and          !
   ! potential temperature.  We simplify the equations by assuming that R and Cp are       !
   ! constants.                                                                            !
   !---------------------------------------------------------------------------------------!
   real(kind=4) function extheta2temp(exner,theta)
      use consts_coms, only : p00i           & ! intent(in)
                            , cpdryi         ! ! intent(in)

      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      real(kind=4), intent(in) :: exner  ! Exner function                         [ J/kg/K]
      real(kind=4), intent(in) :: theta  ! Potential temperature                  [      K]
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Find potential temperature.                                                    !
      !------------------------------------------------------------------------------------!
      extheta2temp = cpdryi * exner * theta
      !------------------------------------------------------------------------------------!

      return
   end function extheta2temp
   !=======================================================================================!
   !=======================================================================================!







   !=======================================================================================!
   !=======================================================================================!
   !     This function computes the specific (intensive) internal energy of water [J/kg],  !
   ! given the temperature and liquid fraction.                                            !
   !---------------------------------------------------------------------------------------!
   real(kind=4) function tl2uint(temp,fliq)
      use consts_coms, only : cice           & ! intent(in)
                            , cliq           & ! intent(in)
                            , tsupercool_liq ! ! intent(in)

      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real(kind=4), intent(in) :: temp  ! Temperature                              [     K]
      real(kind=4), intent(in) :: fliq  ! Fraction liquid water                    [ kg/kg]
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Internal energy is given by the sum of internal energies of ice and liquid     !
      ! phases.                                                                            !
      !------------------------------------------------------------------------------------!
      tl2uint = (1.0 - fliq) * cice * temp + fliq * cliq * (temp - tsupercool_liq)
      !------------------------------------------------------------------------------------!

      return
   end function tl2uint
   !=======================================================================================!
   !=======================================================================================!







   !=======================================================================================!
   !=======================================================================================!
   !     This function computes the extensive internal energy of water [J/m�] or [  J/m�], !
   ! given the temperature [K], the heat capacity of the "dry" part [J/m�/K] or [J/m�/K],  !
   ! water mass [ kg/m�] or [ kg/m�], and liquid fraction [---].                           !
   !---------------------------------------------------------------------------------------!
   real(kind=4) function cmtl2uext(dryhcap,wmass,temp,fliq)
      use consts_coms, only : cice           & ! intent(in)
                            , cliq           & ! intent(in)
                            , tsupercool_liq ! ! intent(in)

      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      real(kind=4), intent(in)  :: dryhcap ! Heat cap. of "dry" part   [J/m�/K] or [J/m�/K]
      real(kind=4), intent(in)  :: wmass   ! Water mass                [ kg/m�] or [ kg/m�]
      real(kind=4), intent(in)  :: temp    ! Temperature                           [     K]
      real(kind=4), intent(in)  :: fliq    ! Liquid fraction (0-1)                 [   ---]
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Internal energy is given by the sum of internal energies of dry part, plus the !
      ! contribution of ice and liquid phases.                                             !
      !------------------------------------------------------------------------------------!
      cmtl2uext = dryhcap * temp + wmass * ( (1.0 - fliq) * cice * temp                    &
                                           + fliq * cliq * (temp - tsupercool_liq) )
      !------------------------------------------------------------------------------------!

      return
   end function cmtl2uext
   !=======================================================================================!
   !=======================================================================================!







   !=======================================================================================!
   !=======================================================================================!
   !     This function computes the specific enthalpy [J/kg] given the temperature and     !
   ! humidity (either mixing ratio or specific humidity).  If we assume that latent heat   !
   ! of vaporisation is a linear function of temperature (equivalent to assume that        !
   ! specific heats are constants and that the thermal expansion of liquids and solids are !
   ! negligible), then the saturation disappears and the enthalpy becomes a straight-      !
   ! forward state function.  In case we are accounting for the water exchange only        !
   ! (latent heat), set the specific humidity to 1.0 and multiply the result by water mass !
   ! or water flux.                                                                        !
   !---------------------------------------------------------------------------------------!
   real(kind=4) function tq2enthalpy(temp,humi,is_shv)
      use consts_coms, only : cpdry          & ! intent(in)
                            , cph2o          & ! intent(in)
                            , tsupercool_vap ! ! intent(in)

      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      real(kind=4), intent(in) :: temp   ! Temperature                             [     K]
      real(kind=4), intent(in) :: humi   ! Humidity (spec. hum. or mixing ratio)   [ kg/kg]
      logical     , intent(in) :: is_shv ! Input humidity is specific humidity     [   T|F]
      !----- Local variables. -------------------------------------------------------------!
      real(kind=4)             :: shv    ! Specific humidity                       [ kg/kg]
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Copy specific humidity to shv.                                                 !
      !------------------------------------------------------------------------------------!
      if (is_shv) then
         shv = humi
      else
         shv = humi / (humi + 1.0)
      end if
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Enthalpy is the combination of dry and moist enthalpies, with the latter being !
      ! allowed to change phase.                                                           !
      !------------------------------------------------------------------------------------!
      tq2enthalpy = (1.0 - shv) * cpdry * temp + shv * cph2o * (temp - tsupercool_vap)
      !------------------------------------------------------------------------------------!

      return
   end function tq2enthalpy
   !=======================================================================================!
   !=======================================================================================!




   !=======================================================================================!
   !=======================================================================================!
   !     This function finds the latent heat of vaporisation for a given temperature.  If  !
   ! we use the definition of latent heat (difference in enthalpy between liquid and       !
   ! vapour phases), and assume that the specific heats are constants, latent heat becomes !
   ! a linear function of temperature.                                                     !
   !---------------------------------------------------------------------------------------!
   real(kind=4) function alvl(temp)
      use consts_coms, only : alvl3  & ! intent(in)
                            , dcpvl  & ! intent(in)
                            , t3ple  ! ! intent(in)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      real(kind=4), intent(in) :: temp
      !------------------------------------------------------------------------------------!


      !----- Linear function, using latent heat at the triple point as reference. ---------!
      alvl = alvl3 + dcpvl * (temp - t3ple)
      !------------------------------------------------------------------------------------!

      return
   end function alvl
   !=======================================================================================!
   !=======================================================================================!




   !=======================================================================================!
   !=======================================================================================!
   !     This function computes the ice-vapour equivalent potential temperature from       !
   ! theta_iland the total mixing ratio.  This is equivalent to the equivalent potential   !
   ! temperature considering also the effects of fusion/melting/sublimation.               !
   !     In case you want to find thetae (i.e. without ice) simply set the the logical     !
   ! useice to .false. .                                                                   !
   !---------------------------------------------------------------------------------------!
   real(kind=4) function thetaeiv(thil,pres,temp,rvap,rtot,useice)
      implicit none
      !----- Required arguments. ----------------------------------------------------------!
      real(kind=4), intent(in)           :: thil   ! Ice-liquid potential temp.    [     K]
      real(kind=4), intent(in)           :: pres   ! Pressure                      [    Pa]
      real(kind=4), intent(in)           :: temp   ! Temperature                   [     K]
      real(kind=4), intent(in)           :: rvap   ! Water vapour mixing ratio     [ kg/kg]
      real(kind=4), intent(in)           :: rtot   ! Total mixing ratio            [ kg/kg]
      !----- Optional arguments. ----------------------------------------------------------!
      logical     , intent(in), optional :: useice ! Should I use ice?             [   T|F]
      !----- Local variables. -------------------------------------------------------------!
      real(kind=4)                       :: tlcl   ! Internal LCL temperature      [     K]
      real(kind=4)                       :: plcl   ! Lifting condensation pressure [    Pa]
      real(kind=4)                       :: dzlcl  ! Thickness of lyr. beneath LCL [     m]
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Find the liquid condensation level (LCL).                                      !
      !------------------------------------------------------------------------------------!
      if (present(useice)) then
         call lcl_il(thil,pres,temp,rtot,rvap,tlcl,plcl,dzlcl,useice)
      else
         call lcl_il(thil,pres,temp,rtot,rvap,tlcl,plcl,dzlcl)
      end if
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     The definition of the thetae_iv is the thetae_ivs at the LCL. The LCL, in turn !
      ! is the point in which rtot = rvap = rsat, so at the LCL rliq = rice = 0.           !
      !------------------------------------------------------------------------------------!
      thetaeiv  = thetaeivs(thil,tlcl,rtot,0.,0.)
      !------------------------------------------------------------------------------------!

      return
   end function thetaeiv
   !=======================================================================================!
   !=======================================================================================!





   !=======================================================================================!
   !=======================================================================================!
   !    This function computes the saturation ice-vapour equivalent potential temperature  !
   ! from theta_il and the total mixing ratio (split into saturated vapour plus liquid and !
   ! ice. This is equivalent to the equivalent potential temperature considering also the  !
   ! effects of fusion/melting/sublimation, and it is done separatedly from the regular    !
   ! thetae_iv because it doesn't require iterations.                                      !
   !                                                                                       !
   !    References:                                                                        !
   !    Tripoli, J. T.; and Cotton, W.R., 1981: The use of ice-liquid water potential tem- !
   !        perature as a thermodynamic variable in deep atmospheric models. Mon. Wea.     !
   !        Rev., v. 109, 1094-1102. (TC81)                                                !
   !                                                                                       !
   !    Some algebra was needed to find this equation, essentially combining (TC81-26) and !
   ! (TC81-27), and the conservation of total water (TC81-16). It assumes that the divi-   !
   ! sion between the three phases is already taken care of.                               !
   !---------------------------------------------------------------------------------------!
   real(kind=4) function thetaeivs(thil,temp,rsat,rliq,rice)
      use consts_coms, only : cpdry    ! ! intent(in)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      real(kind=4), intent(in) :: thil  ! Theta_il, ice-liquid water pot. temp.    [     K]
      real(kind=4), intent(in) :: temp  ! Temperature                              [     K]
      real(kind=4), intent(in) :: rsat  ! Saturation water vapour mixing ratio     [ kg/kg]
      real(kind=4), intent(in) :: rliq  ! Liquid water mixing ratio                [ kg/kg]
      real(kind=4), intent(in) :: rice  ! Ice mixing ratio                         [ kg/kg]
      !----- Local variables --------------------------------------------------------------!
      real(kind=4)             :: rtots ! Saturated mixing ratio                   [     K]
      !------------------------------------------------------------------------------------!


      !------ Find the total saturation mixing ratio. -------------------------------------!
      rtots = rsat+rliq+rice
      !------------------------------------------------------------------------------------!


      !------ Find the saturation equivalent potential temperature. -----------------------!
      thetaeivs = thil * exp ( alvl(temp) * rtots / (cpdry * temp))
      !------------------------------------------------------------------------------------!

      return
   end function thetaeivs
   !=======================================================================================!
   !=======================================================================================!





   !=======================================================================================!
   !=======================================================================================!
   !    This subroutine finds the lifting condensation level given the ice-liquid          !
   ! potential temperature in Kelvin, temperature in Kelvin, the pressure in Pascal, and   !
   ! the mixing ratio in kg/kg. The output will give the LCL temperature and pressure, and !
   ! the thickness of the layer between the initial point and the LCL.                     !
   !                                                                                       !
   !    References:                                                                        !
   !    Tripoli, J. T.; and Cotton, W.R., 1981: The use of ice-liquid water potential      !
   !        temperature as a thermodynamic variable in deep atmospheric models. Mon. Wea.  !
   !        Rev., v. 109, 1094-1102. (TC81)                                                !
   !    Bolton, D., 1980: The computation of the equivalent potential temperature. Mon.    !
   !        Wea. Rev., v. 108, 1046-1053. (BO80)                                           !
   !                                                                                       !
   !    Some algebra was needed to find this equation, essentially combining (TC81-26) and !
   ! (TC81-27), and the conservation of total water (TC81-16). It assumes that the divi-   !
   ! sion between the three phases is already taken care of.                               !
   !    Iterative procedure is needed, and here we iterate looking for T(LCL). Theta_il    !
   ! can be rewritten in terms of T(LCL) only, and once we know this thetae_iv becomes     !
   ! straightforward. T(LCL) will be found using Newton's method, and in the unlikely      !
   ! event it fails,we will fall back to the modified regula falsi (Illinois method).      !
   !                                                                                       !
   ! Important remarks:                                                                    !
   ! 1. TLCL and PLCL are the actual TLCL and PLCL, so in case condensation exists, they   !
   !    will be larger than the actual temperature and pressure (because one would go down !
   !    to reach the equilibrium);                                                         !
   ! 2. DZLCL WILL BE SET TO ZERO in case the LCL is beneath the starting level. So in     !
   !    case you want to force TLCL <= TEMP and PLCL <= PRES, you can use this variable    !
   !    to run the saturation check afterwards. DON'T CHANGE PLCL and TLCL here, they will !
   !    be used for conversions between theta_il and thetae_iv as they are defined here.   !
   ! 3. In case you don't want ice, simply pass useice=.false.. Otherwise let the model    !
   !    decide by itself based on the LEVEL variable.                                      !
   !---------------------------------------------------------------------------------------!
   subroutine lcl_il(thil,pres,temp,rtot,rvap,tlcl,plcl,dzlcl,useice)
      use consts_coms, only : cpog     & ! intent(in)
                            , ep       & ! intent(in)
                            , p00      & ! intent(in)
                            , rocp     & ! intent(in)
                            , t3ple    & ! intent(in)
                            , t00      ! ! intent(in)
      implicit none
      !----- Required arguments. ----------------------------------------------------------!
      real(kind=4), intent(in)            :: thil      ! Ice liquid pot. temp. (*)[      K]
      real(kind=4), intent(in)            :: pres      ! Pressure                 [     Pa]
      real(kind=4), intent(in)            :: temp      ! Temperature              [      K]
      real(kind=4), intent(in)            :: rtot      ! Total mixing ratio       [  kg/kg]
      real(kind=4), intent(in)            :: rvap      ! Vapour mixing ratio      [  kg/kg]
      real(kind=4), intent(out)           :: tlcl      ! LCL temperature          [      K]
      real(kind=4), intent(out)           :: plcl      ! LCL pressure             [     Pa]
      real(kind=4), intent(out)           :: dzlcl     ! Sub-LCL layer thickness  [      m]
      !------------------------------------------------------------------------------------!
      ! (*) This is the most general variable. Thil is exactly theta for no condensation   !
      !     condition, and it is the liquid potential temperature if no ice is present.    !
      !------------------------------------------------------------------------------------!
      !----- Optional arguments. ----------------------------------------------------------!
      logical     , intent(in) , optional :: useice    ! May use ice thermodyn.?  [    T|F]
      !----- Local variables. -------------------------------------------------------------!
      real(kind=4)                        :: pvap      ! Sat. vapour pressure
      real(kind=4)                        :: deriv     ! Function derivative
      real(kind=4)                        :: funnow    ! Current function evaluation
      real(kind=4)                        :: funa      ! Smallest guess function
      real(kind=4)                        :: funz      ! Largest  guess function
      real(kind=4)                        :: tlcla     ! Smallest guess (Newton: previous)
      real(kind=4)                        :: tlclz     ! Largest guess (Newton: new)
      real(kind=4)                        :: es00      ! Defined as p00*rt/(epsilon + rt)
      real(kind=4)                        :: delta     ! Aux. variable for bisection
      integer                             :: itn       ! Iteration counter
      integer                             :: itb       ! Iteration counter
      logical                             :: converged ! Convergence flag
      logical                             :: zside     ! Flag to check sides
      logical                             :: frozen    ! Will use ice thermodyn.  [    T|F]
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Check whether ice thermodynamics is the way to go.                             !
      !------------------------------------------------------------------------------------!
      if (present(useice)) then
         frozen = useice
      else
         frozen = bulk_on
      end if
      !------------------------------------------------------------------------------------!

      !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
      !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
      !write (unit=21,fmt='(a)') '----------------------------------------------------------'
      !write (unit=21,fmt='(a,1x,i5,1x,5(a,1x,f11.4,1x))')                                  &
      !   'INPUT : it=',-1,'thil=',thil,'pres=',0.01*pres,'temp=',temp-t00                  &
      !        ,'rvap=',rvap*1000.,'rtot=',rtot*1000.
      !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
      !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!


      !----- Find es00, which is a constant. ----------------------------------------------!
      es00 = p00 * rtot / (ep+rtot)
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     The 1st. guess, use equation 21 from Bolton (1980). For this we'll need the    !
      ! vapour pressure.                                                                   !
      !------------------------------------------------------------------------------------!
      pvap      = pres * rvap / (ep + rvap)
      tlclz     = 55. + 2840. / (3.5 * log(temp) - log(0.01*pvap) - 4.805)
      pvap      = eslif(tlclz,frozen)
      funnow    = tlclz * (es00/pvap)**rocp - thil
      deriv     = (funnow+thil)*(1./tlclz - rocp*eslifp(tlclz,frozen)/pvap)
      !------------------------------------------------------------------------------------!

      !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
      !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
      !write (unit=21,fmt='(a,1x,i5,1x,2(a,1x,f11.4,1x),2(a,1x,es11.4,1x))')                &
      !   'NEWTON: it=',0,'tlclz=',tlclz-t00,'pvap=',0.01*pvap,'fun=',funnow,'deriv=',deriv
      !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
      !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
      tlcla     = tlclz
      funa      = funnow


      !------------------------------------------------------------------------------------!
      !      First loop: Newton's method.                                                  !
      !------------------------------------------------------------------------------------!
      newloop: do itn=1,maxfpo/6
         !----- If derivative is too small, skip Newton's and try bisection instead. ------!
         if (abs(deriv) < toler) exit newloop
         !---------------------------------------------------------------------------------!


         !----- Otherwise, update guesses. ------------------------------------------------!
         tlcla  = tlclz
         funa   = funnow
         tlclz  = tlcla - funnow/deriv
         pvap   = eslif(tlclz,frozen)
         funnow = tlclz * (es00/pvap)**rocp - thil
         deriv  = (funnow+thil)*(1./tlclz - rocp*eslifp(tlclz,frozen)/pvap)

         !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
         !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
         !write (unit=21,fmt='(a,1x,i5,1x,2(a,1x,f11.4,1x),2(a,1x,es11.4,1x))')             &
         !   'NEWTON: it=',itn,'tlclz=',tlclz-t00,'pvap=',0.01*pvap,'fun=',funnow           &
         !          ,'deriv=',deriv
         !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
         !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!

         !---------------------------------------------------------------------------------!
         !      Check for convergence.                                                     !
         !---------------------------------------------------------------------------------!
         converged = abs(tlcla-tlclz) < toler*tlclz
         if (converged) then
            !----- Guesses are almost identical, average them. ----------------------------!
            tlcl = 0.5*(tlcla+tlclz)
            funz = funnow
            exit newloop
            !------------------------------------------------------------------------------!
         elseif (funnow == 0.) then
            !----- We've hit the answer by luck, copy the answer. -------------------------!
            tlcl = tlclz
            funz = funnow
            converged = .true.
            exit newloop
            !------------------------------------------------------------------------------!
         end if
      end do newloop
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !      Check whether Newton's method has converged.                                  !
      !------------------------------------------------------------------------------------!
      if (.not. converged) then
         !---------------------------------------------------------------------------------!
         !     Newton's method has failed.  We use regula falsi instead.  First, we must   !
         ! find two guesses whose function evaluations have opposite signs.                !
         !---------------------------------------------------------------------------------!
         if (funa*funnow < 0. ) then
            !----- We already have two good guesses. --------------------------------------!
            funz  = funnow
            zside = .true.
            !------------------------------------------------------------------------------!
         else
            !------------------------------------------------------------------------------!
            !     We need to find another guess with opposite sign.                        !
            !------------------------------------------------------------------------------!

            !----- We fix funa, and try a funz that will work as 2nd guess ----------------!
            if (abs(funnow-funa) < toler*tlcla) then
               delta = 100.*toler*tlcla
            else
               delta = max(abs(funa*(tlclz-tlcla)/(funnow-funa)),100.*toler*tlcla)
            end if
            tlclz = tlcla + delta
            !------------------------------------------------------------------------------!

            !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
            !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
            !write (unit=21,fmt='(a,1x,i5,1x,3(a,1x,f11.4,1x),3(a,1x,es11.4,1x))')          &
            !   '2NGGSS: tt=',0,'tlclz=',tlclz-t00,'tlcla=',tlcla-t00,'pvap=',0.01*pvap     &
            !           ,'funa=',funa,'funz=',funnow,'delta=',delta
            !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
            !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
            zside = .false.
            zgssloop: do itb=1,maxfpo

               !----- So the first time tlclz = tlcla - 2*delta ---------------------------!
               tlclz = tlcla + real((-1)**itb * (itb+3)/2) * delta
               pvap  = eslif(tlclz,frozen)
               funz  = tlclz * (es00/pvap)**rocp - thil
               !---------------------------------------------------------------------------!


               !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
               !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
               !write (unit=21,fmt='(a,1x,i5,1x,2(a,1x,f11.4,1x),2(a,1x,es11.4,1x))')       &
               !   '2NGGSS: tt=',itb,'tlclz=',tlclz-t00,'pvap=',0.01*pvap,'fun=',funz       &
               !           ,'delta=',delta
               !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
               !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
               zside = funa*funz < 0.0
               if (zside) exit zgssloop
            end do zgssloop
            if (.not. zside) then
               write (unit=*,fmt='(a)') ' ====== No second guess for you... ======'
               write (unit=*,fmt='(a)') ' + INPUT variables: '
               write (unit=*,fmt='(a,1x,es14.7)') 'THIL =',thil
               write (unit=*,fmt='(a,1x,es14.7)') 'TEMP =',temp
               write (unit=*,fmt='(a,1x,es14.7)') 'PRES =',pres
               write (unit=*,fmt='(a,1x,es14.7)') 'RTOT =',rtot
               write (unit=*,fmt='(a,1x,es14.7)') 'RVAP =',rvap
               write (unit=*,fmt='(a)') ' ============ Failed guess... ==========='
               write (unit=*,fmt='(2(a,1x,es14.7))') 'TLCLA =',tlcla,'FUNA =',funa
               write (unit=*,fmt='(2(a,1x,es14.7))') 'TLCLZ =',tlclz,'FUNC =',funz
               write (unit=*,fmt='(2(a,1x,es14.7))') 'DELTA =',delta,'FUNN =',funnow
               call fatal_error('Failed finding the second guess for regula falsi'         &
                             ,'lcl_il','therm_lib.f90')
            end if
         end if
         !---------------------------------------------------------------------------------!


         !---------------------------------------------------------------------------------!
         !      We have the guesses, solve the regula falsi method.                        !
         !---------------------------------------------------------------------------------!
         fpoloop: do itb=itn+1,maxfpo
            !----- Update guess and function evaluation. ----------------------------------!
            tlcl   = (funz*tlcla-funa*tlclz)/(funz-funa)
            pvap   = eslif(tlcl,frozen)
            funnow = tlcl * (es00/pvap)**rocp - thil
            !------------------------------------------------------------------------------!

            !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
            !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
            !write (unit=21,fmt='(a,1x,i5,1x,2(a,1x,f11.4,1x),1(a,1x,es11.4,1x))')          &
            !   'REGFAL: it=',itb,'tlcl =',tlcl -t00,'pvap=',0.01*pvap,'fun=',funnow
            !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
            !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!

            !------------------------------------------------------------------------------!
            !    Check for convergence.  If it did, return, we found the solution.         !
            ! Otherwise, we update one of the guesses.                                     !
            !------------------------------------------------------------------------------!
            converged = abs(tlcl-tlcla) < toler*tlcl .and.  abs(tlcl-tlclz) < toler*tlcl
            if (funnow == 0. .or. converged) then
               converged = .true.
               exit fpoloop
            elseif (funnow*funa < 0.) then
               tlclz = tlcl
               funz  = funnow
               !----- If we are updating zside again, modify aside (Illinois method) ------!
               if (zside) funa=funa * 0.5
               !---------------------------------------------------------------------------!


               !----- We have just updated zside, sett zside to true. ---------------------!
               zside = .true.
               !---------------------------------------------------------------------------!
            else
               tlcla = tlcl
               funa  = funnow
               !----- If we are updating aside again, modify zside (Illinois method) ------!
               if (.not.zside) funz = funz * 0.5
               !---------------------------------------------------------------------------!


               !----- We have just updated aside, set zside to false. ---------------------!
               zside = .false.
               !---------------------------------------------------------------------------!
            end if
         end do fpoloop
      end if
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Check whether we have succeeded or not.                                        !
      !------------------------------------------------------------------------------------!
      if (converged) then
         !----- We have found a solution, find the remaining LCL properties. --------------!
         pvap  = eslif(tlcl,frozen)
         plcl  = (ep + rvap) * pvap / rvap
         dzlcl = max(cpog*(temp-tlcl),0.)
         !---------------------------------------------------------------------------------!


         !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
         !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
         !write (unit=21,fmt='(a,1x,i5,1x,3(a,1x,f11.4,1x),3(a,1x,es11.4,1x))')             &
         !   'ANSWER: itb=',itn,'tlcl=',tlcl-t00,'eslcl=',0.01*pvap                         &
         !        ,'dzlcl=',dzlcl,'plcl=',plcl*0.01,'funa=',funa,'funz=',funz
         !write (unit=21,fmt='(a)') '-------------------------------------------------------'
         !write (unit=21,fmt='(a)') ' '
         !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
         !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
      else
         write (unit=*,fmt='(a)') '-------------------------------------------------------'
         write (unit=*,fmt='(a)') ' LCL Temperature didn''t converge!!!'
         write (unit=*,fmt='(a,1x,i5,1x,a)') ' I gave up, after',maxfpo,'iterations...'
         write (unit=*,fmt='(a)') ' '
         write (unit=*,fmt='(a)') ' Input values.'
         write (unit=*,fmt='(a)') ' '
         write (unit=*,fmt='(a,1x,f12.4)' ) 'theta_il        [     K] =',thil
         write (unit=*,fmt='(a,1x,f12.4)' ) 'Pressure        [   hPa] =',0.01*pres
         write (unit=*,fmt='(a,1x,f12.4)' ) 'Temperature     [    �C] =',temp-t00
         write (unit=*,fmt='(a,1x,f12.4)' ) 'rtot            [  g/kg] =',1000.*rtot
         write (unit=*,fmt='(a,1x,f12.4)' ) 'rvap            [  g/kg] =',1000.*rvap
         write (unit=*,fmt='(a)') ' '
         write (unit=*,fmt='(a)') ' Last iteration outcome.'
         write (unit=*,fmt='(a,1x,f12.4)' ) 'tlcla           [    �C] =',tlcla-t00
         write (unit=*,fmt='(a,1x,f12.4)' ) 'tlclz           [    �C] =',tlclz-t00
         write (unit=*,fmt='(a,1x,f12.4)' ) 'fun             [     K] =',funnow
         write (unit=*,fmt='(a,1x,f12.4)' ) 'funa            [     K] =',funa
         write (unit=*,fmt='(a,1x,f12.4)' ) 'funz            [     K] =',funz
         write (unit=*,fmt='(a,1x,f12.4)' ) 'deriv           [  ----] =',deriv
         write (unit=*,fmt='(a,1x,es12.4)') 'toler           [  ----] =',toler
         write (unit=*,fmt='(a,1x,es12.4)') 'error           [  ----] ='                   &
                                                            ,abs(tlclz-tlcla)/tlclz
         write (unit=*,fmt='(a,1x,f12.4)' ) 'tlcl            [    �C] =',tlcl
         call fatal_error('TLCL didn''t converge, gave up!','lcl_il','therm_lib.f90')
      end if
      return
   end subroutine lcl_il
   !=======================================================================================!
   !=======================================================================================!




   !=======================================================================================!
   !=======================================================================================!
   !    This subroutine computes the temperature and fraction of liquid water from the     !
   ! intensive internal energy [J/kg].                                                     !
   !---------------------------------------------------------------------------------------!
   subroutine uint2tl(uint,temp,fliq)
      use consts_coms, only : cliqi          & ! intent(in)
                            , cicei          & ! intent(in)
                            , allii          & ! intent(in)
                            , t3ple          & ! intent(in)
                            , uiicet3        & ! intent(in)
                            , uiliqt3        & ! intent(in)
                            , tsupercool_liq ! ! intent(in)
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real(kind=4), intent(in)  :: uint     ! Internal energy                     [   J/kg]
      real(kind=4), intent(out) :: temp     ! Temperature                         [      K]
      real(kind=4), intent(out) :: fliq  ! Liquid Fraction (0-1)               [    ---]
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Compare the internal energy with the reference values to decide which phase    !
      ! the water is.                                                                      !
      !------------------------------------------------------------------------------------!
      if (uint <= uiicet3) then
         !----- Internal energy below qwfroz, all ice  ------------------------------------!
         fliq = 0.
         temp = uint * cicei
         !---------------------------------------------------------------------------------!
      elseif (uint >= uiliqt3) then
         !----- Internal energy, above qwmelt, all liquid ---------------------------------!
         fliq = 1.
         temp = uint * cliqi + tsupercool_liq
         !---------------------------------------------------------------------------------!
      else
         !----- Changing phase, it must be at freezing point ------------------------------!
         fliq = (uint - uiicet3) * allii
         temp = t3ple
         !---------------------------------------------------------------------------------!
      end if
      !------------------------------------------------------------------------------------!

      return
   end subroutine uint2tl
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !    This subroutine computes the temperature (Kelvin) and liquid fraction from         !
   ! extensive internal energy (J/m� or J/m�), water mass (kg/m� or kg/m�), and heat       !
   ! capacity (J/m�/K or J/m�/K).                                                          !
   !---------------------------------------------------------------------------------------!
   subroutine uextcm2tl(uext,wmass,dryhcap,temp,fliq)
      use consts_coms, only : cliqi          & ! intent(in)
                            , cliq           & ! intent(in)
                            , cicei          & ! intent(in)
                            , cice           & ! intent(in)
                            , allii          & ! intent(in)
                            , alli           & ! intent(in)
                            , t3ple          & ! intent(in)
                            , tsupercool_liq ! ! intent(in)
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real(kind=4), intent(in)  :: uext    ! Extensive internal energy [  J/m�] or [  J/m�]
      real(kind=4), intent(in)  :: wmass   ! Water mass                [ kg/m�] or [ kg/m�]
      real(kind=4), intent(in)  :: dryhcap ! Heat cap. of "dry" part   [J/m�/K] or [J/m�/K]
      real(kind=4), intent(out) :: temp    ! Temperature                           [     K]
      real(kind=4), intent(out) :: fliq    ! Liquid fraction (0-1)                 [   ---]
      !----- Local variable ---------------------------------------------------------------!
      real(kind=4)              :: uefroz  ! qw of ice at triple pt.   [  J/m�] or [  J/m�]
      real(kind=4)              :: uemelt  ! qw of liq. at triple pt.  [  J/m�] or [  J/m�]
      !------------------------------------------------------------------------------------!



      !----- Convert melting heat to J/m� or J/m� -----------------------------------------!
      uefroz = (dryhcap + wmass * cice) * t3ple
      uemelt = uefroz   + wmass * alli
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !    This is analogous to the uint2tl computation, we should analyse the magnitude   !
      ! of the internal energy to choose between liquid, ice, or both by comparing with    !
      ! the known boundaries.                                                              !
      !------------------------------------------------------------------------------------!
      if (uext < uefroz) then
         !----- Internal energy below qwfroz, all ice  ------------------------------------!
         fliq = 0.
         temp = uext  / (cice * wmass + dryhcap)
         !---------------------------------------------------------------------------------!
      elseif (uext > uemelt) then
         !----- Internal energy, above qwmelt, all liquid ---------------------------------!
         fliq = 1.
         temp = (uext + wmass * cliq * tsupercool_liq) / (dryhcap + wmass * cliq)
         !---------------------------------------------------------------------------------!
      elseif (uefroz == uemelt) then
         !---------------------------------------------------------------------------------!
         !    We are at the freezing point.  If water mass is so tiny that the internal    !
         ! energy of frozen and melted states are the same given the machine precision,    !
         ! then we assume that water content is negligible and we impose 50% frozen for    !
         ! simplicity.                                                                     !
         !---------------------------------------------------------------------------------!
         fliq = 0.5
         temp = t3ple
         !---------------------------------------------------------------------------------!
      else
         !---------------------------------------------------------------------------------!
         !    Changing phase, it must be at freezing point.  The max and min are here just !
         ! to avoid tiny deviations beyond 0. and 1. due to floating point arithmetics.    !
         !---------------------------------------------------------------------------------!
         fliq = min(1.,max(0.,(uext - uefroz) * allii / wmass))
         temp = t3ple
         !---------------------------------------------------------------------------------!
      end if
      !------------------------------------------------------------------------------------!

      return
   end subroutine uextcm2tl
   !=======================================================================================!
   !=======================================================================================!
end module therm_lib
!==========================================================================================!
!==========================================================================================!
