!==========================================================================================!
!. File: therm_lib8.f90                                                                    !
!                                                                                          !
!  Based on BRAMS-4.0.6   This file contains most functions and subroutines that deal with !
! several thermodynamic conversions that are needed in double precision.  Most of them     !
! have the equivalent in single precision in therm_lib.  These procedures were built to    !
! avoid assumptions like hydrostatic and linearisation.  Most equations could not be       !
! solved analytically, and the standard here was to use Newton's method as the default,    !
! always having bisection or, more often, the modified Regula Falsi (Illinois) method in   !
! case Newton's fails.                                                                     !
!==========================================================================================!
!==========================================================================================!
module therm_lib8
   use therm_lib, only : toler4     => toler     & ! intent(in)
                       , maxfpo4    => maxfpo    & ! intent(in)
                       , maxit4     => maxit     & ! intent(in)
                       , maxlev4    => maxlev    & ! intent(in)
                       , newthermo4 => newthermo & ! intent(in)
                       , level4     => level     & ! intent(in)
                       , vapour_on4 => vapour_on & ! intent(in)
                       , cloud_on4  => cloud_on  & ! intent(in)
                       , bulk_on4   => bulk_on   ! ! intent(in)

   !---------------------------------------------------------------------------------------!
   !     Relative tolerance for iterative methods. The smaller the value, the more         !
   ! accurate the result, but it will slow down the run.  Notice that we are using the     !
   ! tolerance that is based on the single precision...                                    !
   !---------------------------------------------------------------------------------------!
   real(kind=8), parameter ::   toler8  = dble(toler4) ! Relative tolerance for iterative
                                                       !    methods.  The smaller the
                                                       !    value, the more accurate the
                                                       !    result, but smaller values will
                                                       !    slow down the run.
   integer     , parameter ::   maxfpo = maxfpo4       ! Maximum # of iterations before
                                                       !    crashing for false position
                                                       !    method.
   integer     , parameter ::   maxit  = maxit4        ! Maximum # of iterations before
                                                       !    crashing, for other methods.
   integer     , parameter ::   maxlev = maxlev4       ! Maximum # of levels for adaptive
                                                       !    quadrature methods.
   logical     , parameter ::   newthermo = newthermo4 ! Use new thermodynamics [T|F]
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !   This is the "level" variable, that used to be in micphys. Since it affects more the !
   ! thermodynamics choices than the microphysics, it was moved to here.                   !
   !---------------------------------------------------------------------------------------!
   integer, parameter ::   level = level4
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !    The following three variables are just the logical tests on variable "level",      !
   ! saved here to speed up checks for "li" functions.                                     !
   !---------------------------------------------------------------------------------------!
   logical, parameter ::   vapour_on = vapour_on4
   logical, parameter ::   cloud_on  = cloud_on4
   logical, parameter ::   bulk_on   = bulk_on4
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
   real(kind=8), dimension(0:3), parameter :: iii_78 = (/ 9.550426d0, -5.723265d3          &
                                                        , 3.53068d0,  -7.28332d-3 /)
   !----- Coefficients based on equation (10), first fit ----------------------------------!
   real(kind=8), dimension(0:3), parameter :: l01_108 = (/ 5.4842763d1,-6.76322d3          &
                                                         ,-4.210d0    , 3.67d-4 /)
   !----- Coefficients based on equation (10), second fit ---------------------------------!
   real(kind=8), dimension(0:3), parameter :: l02_108 = (/ 5.3878d1   ,-1.33122d3          &
                                                         ,-9.44523d0  , 1.4025d-2 /)
   !----- Coefficients based on the hyperbolic tangent ------------------------------------!
   real(kind=8), dimension(2)  , parameter :: ttt_108 = (/4.15d-2     , 2.188d2 /)
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
   real(kind=8), dimension(0:8), parameter :: cll8 = (/  .6105851d+03,  .4440316d+02       &
                                                      ,  .1430341d+01,  .2641412d-01       &
                                                      ,  .2995057d-03,  .2031998d-05       &
                                                      ,  .6936113d-08,  .2564861d-11       &
                                                      , -.3704404d-13                /)
   !----- Coefficients for esat (ice) -----------------------------------------------------!
   real(kind=8), dimension(0:8), parameter :: cii8 = (/  .6114327d+03,  .5027041d+02       &
                                                      ,  .1875982d+01,  .4158303d-01       &
                                                      ,  .5992408d-03,  .5743775d-05       &
                                                      ,  .3566847d-07,  .1306802d-09       &
                                                      ,  .2152144d-12                /)
   !----- Coefficients for d(esat)/dT (liquid) --------------------------------------------!
   real(kind=8), dimension(0:8), parameter :: dll8 = (/  .4443216d+02,  .2861503d+01       &
                                                      ,  .7943347d-01,  .1209650d-02       &
                                                      ,  .1036937d-04,  .4058663d-07       &
                                                      , -.5805342d-10, -.1159088d-11       &
                                                      , -.3189651d-14                /)
   !----- Coefficients for esat (ice) -----------------------------------------------------!
   real(kind=8), dimension(0:8), parameter :: dii8 = (/  .5036342d+02,  .3775758d+01       &
                                                      ,  .1269736d+00,  .2503052d-02       &
                                                      ,  .3163761d-04,  .2623881d-06       &
                                                      ,  .1392546d-08,  .4315126d-11       &
                                                      ,  .5961476d-14                /)
   !---------------------------------------------------------------------------------------!
   !=======================================================================================!
   !=======================================================================================!


   contains



   !=======================================================================================!
   !=======================================================================================!
   !     This function calculates the liquid saturation vapour pressure as a function of   !
   ! Kelvin temperature. This expression came from MK05, equation (10).                    !
   !---------------------------------------------------------------------------------------!
   real(kind=8) function eslf8(temp,l1funout,l2funout,ttfunout)
      use consts_coms, only : t008 ! ! intent(in)
      implicit none
      !----- Required arguments. ----------------------------------------------------------!
      real(kind=8), intent(in)            :: temp     ! Temperature                [     K]
      !----- Optional arguments. ----------------------------------------------------------!
      real(kind=8), intent(out), optional :: l1funout ! Function for high temperatures
      real(kind=8), intent(out), optional :: ttfunout ! Interpolation function
      real(kind=8), intent(out), optional :: l2funout ! Function for low temperatures
      !----- Local variables. -------------------------------------------------------------!
      real(kind=8)                        :: l1fun    !
      real(kind=8)                        :: ttfun    !
      real(kind=8)                        :: l2fun    !
      real(kind=8)                        :: x        !
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Choose between the old and the new thermodynamics.                             !
      !------------------------------------------------------------------------------------!
      if (newthermo) then
         !----- Updated method, using MK05 ------------------------------------------------!
         l1fun = l01_108(0) + l01_108(1)/temp + l01_108(2)*log(temp) + l01_108(3) * temp
         l2fun = l02_108(0) + l02_108(1)/temp + l02_108(2)*log(temp) + l02_108(3) * temp
         ttfun = tanh(ttt_108(1) * (temp - ttt_108(2)))
         eslf8 = exp(l1fun + ttfun*l2fun)
         !---------------------------------------------------------------------------------!

         if (present(l1funout)) l1funout = l1fun
         if (present(l2funout)) l2funout = l2fun
         if (present(ttfunout)) ttfunout = ttfun
      else
         !----- Original method, using polynomial fit (FWC92) -----------------------------!
         x     = max(-8.0d1,temp-t008)
         eslf8 = cll8(0) + x * (cll8(1) + x * (cll8(2) + x * (cll8(3) + x * (cll8(4)       &
                         + x * (cll8(5) + x * (cll8(6) + x * (cll8(7) + x * cll8(8)) ))))))
         !---------------------------------------------------------------------------------!

         if (present(l1funout)) l1funout = eslf8
         if (present(l2funout)) l2funout = eslf8
         if (present(ttfunout)) ttfunout = eslf8
      end if

      return
   end function eslf8
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This function calculates the ice saturation vapour pressure as a function of      !
   ! Kelvin temperature, based on MK05 equation (7).                                       !
   !---------------------------------------------------------------------------------------!
   real(kind=8) function esif8(temp,iifunout)
      use consts_coms, only : t008 ! ! intent(in)
      implicit none
      !----- Required arguments. ----------------------------------------------------------!
      real(kind=8), intent(in)            :: temp     ! Temperature                 [    K]
      !----- Optional arguments. ----------------------------------------------------------!
      real(kind=8), intent(out), optional :: iifunout
      !----- Local variables. -------------------------------------------------------------!
      real(kind=8)                        :: iifun
      real(kind=8)                        :: x
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Choose between the old and the new thermodynamics.                             !
      !------------------------------------------------------------------------------------!
      if (newthermo) then

         !----- Updated method, using MK05 ------------------------------------------------!
         iifun = iii_78(0) + iii_78(1)/temp + iii_78(2) * log(temp) + iii_78(3) * temp
         esif8 = exp(iifun)
         !---------------------------------------------------------------------------------!


         if (present(iifunout)) iifunout=iifun
      else
         !----- Original method, using polynomial fit (FWC92) -----------------------------!
         x     = max(-8.d1,temp-t008)
         esif8 = cii8(0) + x * (cii8(1) + x * (cii8(2) + x * (cii8(3) + x * (cii8(4)       &
                        + x * (cii8(5) + x * (cii8(6) + x * (cii8(7) + x * cii8(8))))))))
         !---------------------------------------------------------------------------------!

         if (present(iifunout)) iifunout=esif8
      end if
      !------------------------------------------------------------------------------------!

      return
   end function esif8
   !=======================================================================================!
   !=======================================================================================!





   !=======================================================================================!
   !=======================================================================================!
   !     This function calculates the saturation specific humidity, over liquid or ice     !
   ! depending on temperature, as a function of pressure and Kelvin temperature.           !
   !---------------------------------------------------------------------------------------!
   real(kind=8) function qslif8(pres,temp,useice)
      use consts_coms, only : t3ple8  & ! intent(in)
                            , ep8     & ! intent(in)
                            , toodry8 ! ! intent(in)
      implicit none
      !----- Required arguments. ----------------------------------------------------------!
      real(kind=8), intent(in)           :: pres
      real(kind=8), intent(in)           :: temp
      !----- Optional arguments. ----------------------------------------------------------!
      logical     , intent(in), optional :: useice
      !----- Local variables. -------------------------------------------------------------!
      real(kind=8)                       :: esz
      logical                            :: frozen
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Decide which function to use (saturation for liquid water or ice).             !
      !------------------------------------------------------------------------------------!
      if (present(useice)) then
         frozen = useice  .and. temp < t3ple8
      else
         frozen = bulk_on .and. temp < t3ple8
      end if
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !      Call the appropriate function depending on the temperature and whether ice    !
      ! thermodynamics is to be used.                                                      !
      !------------------------------------------------------------------------------------!
      if (frozen) then
         !----- Saturation vapour pressure for ice. ---------------------------------------!
         esz = esif8(temp)
         !---------------------------------------------------------------------------------!
      else
         !----- Saturation vapour pressure for liquid. ------------------------------------!
         esz = eslf8(temp)
         !---------------------------------------------------------------------------------!
      end if
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Use the usual relation between pressure and vapour pressure to find the        !
      ! specific humidity.                                                                 !
      !------------------------------------------------------------------------------------!
      qslif8 = max(toodry8, ep8 * esz/( pres - (1.d0 - ep8) * esz) )
      !------------------------------------------------------------------------------------!

      return
   end function qslif8
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
   real(kind=8) function rehuil8(pres,temp,humi,is_shv,useice)
      use consts_coms, only : t3ple8  & ! intent(in)
                            , ep8     & ! intent(in)
                            , toodry8 ! ! intent(in)
      implicit none
      !----- Required arguments. ----------------------------------------------------------!
      real(kind=8), intent(in)           :: pres    ! Air pressure                 [    Pa]
      real(kind=8), intent(in)           :: temp    ! Temperature                  [     K]
      real(kind=8), intent(in)           :: humi    ! Humidity                     [ kg/kg]
      logical     , intent(in)           :: is_shv  ! Input is specific humidity   [   T|F]
       !----- Optional arguments. ----------------------------------------------------------!
      logical     , intent(in), optional :: useice  ! May use ice thermodynamics   [   T|F]
      !----- Local variables --------------------------------------------------------------!
      real(kind=8)                       :: shv     ! Specific humidity            [ kg/kg]
      real(kind=8)                       :: pvap    ! Vapour pressure              [    Pa]
      real(kind=8)                       :: psat    ! Saturation vapour pressure   [    Pa]
      logical                            :: frozen  ! Will use ice saturation now  [   T|F]
      !------------------------------------------------------------------------------------!

      !------------------------------------------------------------------------------------!
      !    Check whether we should use ice or liquid saturation.                           !
      !------------------------------------------------------------------------------------!
      if (present(useice)) then
         frozen = useice  .and. temp < t3ple8
      else
         frozen = bulk_on .and. temp < t3ple8
      end if
      !------------------------------------------------------------------------------------!


      !---- Make sure that we have specific humidity. -------------------------------------!
      if (is_shv) then
         shv = max(toodry8,humi)
      else
         shv = max(toodry8,humi) / ( 1.d0 + max(toodry8,humi) )
      end if
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Find the vapour pressure and the saturation vapour pressure.                   !
      !------------------------------------------------------------------------------------!
      pvap = ( pres * shv ) / ( ep8 + (1.d0 - ep8) * shv )
      if (frozen) then
         psat = esif8(temp)
      else
         psat = esif8(temp)
      end if
      !------------------------------------------------------------------------------------!

      !------------------------------------------------------------------------------------!
      !     Find the relative humidity.                                                    !
      !------------------------------------------------------------------------------------!
      rehuil8 = max(0.d0 ,pvap / psat)
      !------------------------------------------------------------------------------------!

      return
   end function rehuil8
   !=======================================================================================!
   !=======================================================================================!



   !=======================================================================================!
   !=======================================================================================!
   !     This function finds the density based on the virtual temperature and the ideal    !
   ! gas law.  The only difference between this function and the one above is that here we !
   ! provide vapour and total specific mass (specific humidity) instead of mixing ratio.   !
   !---------------------------------------------------------------------------------------!
   real(kind=8) function idealdenssh8(pres,temp,qvpr,qtot)
      use consts_coms, only : rdry8 & ! intent(in)
                            , epi8  ! ! intent(in)
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real(kind=8), intent(in)           :: pres ! Pressure                        [    Pa]
      real(kind=8), intent(in)           :: temp ! Temperature                     [     K]
      real(kind=8), intent(in)           :: qvpr ! Vapour specific mass            [ kg/kg]
      real(kind=8), intent(in), optional :: qtot ! Total water specific mass       [ kg/kg]
      !----- Local variables. -------------------------------------------------------------!
      real(kind=8)                       :: qall ! Either qtot or qvpr...          [ kg/kg]
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
      idealdenssh8 = pres / (rdry8 * temp * (1.d0 - qall + epi8 * qvpr))
      !------------------------------------------------------------------------------------!

      return
   end function idealdenssh8
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This function computes reduces the pressure from the reference height to the      !
   ! canopy height by assuming hydrostatic equilibrium.  For simplicity, we assume that    !
   ! R and cp are constants (in reality they are dependent on humidity).                   !
   !---------------------------------------------------------------------------------------!
   real(kind=8) function reducedpress8(pres,thetaref,shvref,zref,thetacan,shvcan,zcan)
      use consts_coms, only : epim18    & ! intent(in)
                            , p00k8     & ! intent(in)
                            , rocp8     & ! intent(in)
                            , cpor8     & ! intent(in)
                            , cpdry8    & ! intent(in)
                            , grav8     ! ! intent(in)
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real(kind=8), intent(in) :: pres     ! Pressure                            [      Pa]
      real(kind=8), intent(in) :: thetaref ! Potential temperature               [       K]
      real(kind=8), intent(in) :: shvref   ! Vapour specific mass                [   kg/kg]
      real(kind=8), intent(in) :: zref     ! Height at reference level           [       m]
      real(kind=8), intent(in) :: thetacan ! Potential temperature               [       K]
      real(kind=8), intent(in) :: shvcan   ! Vapour specific mass                [   kg/kg]
      real(kind=8), intent(in) :: zcan     ! Height at canopy level              [       m]
      !------Local variables. -------------------------------------------------------------!
      real(kind=8)             :: pinc     ! Pressure increment                  [ Pa^R/cp]
      real(kind=8)             :: thvbar   ! Average virtual pot. temperature    [       K]
      !------------------------------------------------------------------------------------!

      !------------------------------------------------------------------------------------!
      !      First we compute the average virtual potential temperature between the canopy !
      ! top and the reference level.                                                       !
      !------------------------------------------------------------------------------------!
      thvbar = 5.d-1 * ( thetaref * (1.d0 + epim18 * shvref)                               &
                       + thetacan * (1.d0 + epim18 * shvcan) )
      !------------------------------------------------------------------------------------!



      !----- Then, we find the pressure gradient scale. -----------------------------------!
      pinc   = grav8 * p00k8 * (zref - zcan) / (cpdry8 * thvbar)
      !------------------------------------------------------------------------------------!


      !----- And we can find the reduced pressure. ----------------------------------------!
      reducedpress8 = (pres**rocp8 + pinc ) ** cpor8
      !------------------------------------------------------------------------------------!

      return
   end function reducedpress8
   !=======================================================================================!
   !=======================================================================================!







   !=======================================================================================!
   !=======================================================================================!
   !     This function computes the Exner function [J/kg/K], given the pressure.  It       !
   ! assumes for simplicity that R and Cp are constants and equal to the dry air values.   !
   !---------------------------------------------------------------------------------------!
   real(kind=8) function press2exner8(pres)
      use consts_coms, only : p00i8           & ! intent(in)
                            , cpdry8          & ! intent(in)
                            , rocp8           ! ! intent(in)

      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      real(kind=8), intent(in) :: pres   ! Pressure                               [     Pa]
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Find potential temperature.                                                    !
      !------------------------------------------------------------------------------------!
      press2exner8 = cpdry8 * ( pres * p00i8 ) ** rocp8
      !------------------------------------------------------------------------------------!

      return
   end function press2exner8
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This function computes the potential temperature [K], given the Exner function    !
   ! and temperature.  For simplicity we ignore the effects of humidity in R and cp and    !
   ! use the dry air values instead.                                                       !
   !---------------------------------------------------------------------------------------!
   real(kind=8) function extemp2theta8(exner,temp)
      use consts_coms, only : cpdry8          ! ! intent(in)

      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      real(kind=8), intent(in) :: exner  ! Exner function                         [ J/kg/K]
      real(kind=8), intent(in) :: temp   ! Temperature                            [      K]
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Find potential temperature.                                                    !
      !------------------------------------------------------------------------------------!
      extemp2theta8 = cpdry8 * temp / exner
      !------------------------------------------------------------------------------------!

      return
   end function extemp2theta8
   !=======================================================================================!
   !=======================================================================================!







   !=======================================================================================!
   !=======================================================================================!
   !     This function computes the temperature [K], given the Exner function and          !
   ! potential temperature.  We simplify the equations by assuming that R and Cp are       !
   ! constants.                                                                            !
   !---------------------------------------------------------------------------------------!
   real(kind=8) function extheta2temp8(exner,theta)
      use consts_coms, only : p00i8           & ! intent(in)
                            , cpdryi8         ! ! intent(in)

      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      real(kind=8), intent(in) :: exner  ! Exner function                         [ J/kg/K]
      real(kind=8), intent(in) :: theta  ! Potential temperature                  [      K]
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Find potential temperature.                                                    !
      !------------------------------------------------------------------------------------!
      extheta2temp8 = cpdryi8 * exner * theta
      !------------------------------------------------------------------------------------!

      return
   end function extheta2temp8
   !=======================================================================================!
   !=======================================================================================!







   !=======================================================================================!
   !=======================================================================================!
   !     This function computes the specific internal energy of water [J/kg], given the    !
   ! temperature and liquid fraction.                                                      !
   !---------------------------------------------------------------------------------------!
   real(kind=8) function tl2uint8(temp,fliq)
      use consts_coms, only : cice8           & ! intent(in)
                            , cliq8           & ! intent(in)
                            , tsupercool_liq8 ! ! intent(in)

      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real(kind=8), intent(in) :: temp  ! Temperature                              [     K]
      real(kind=8), intent(in) :: fliq  ! Fraction liquid water                    [ kg/kg]
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Internal energy is given by the sum of internal energies of ice and liquid     !
      ! phases.                                                                            !
      !------------------------------------------------------------------------------------!
      tl2uint8 = (1.d0 - fliq) * cice8 * temp + fliq * cliq8 * (temp - tsupercool_liq8)
      !------------------------------------------------------------------------------------!

      return
   end function tl2uint8
   !=======================================================================================!
   !=======================================================================================!







   !=======================================================================================!
   !=======================================================================================!
   !     This function computes the internal energy of water [J/m�] or [  J/m�], given the !
   ! temperature [K], the heat capacity of the "dry" part [J/m�/K] or [J/m�/K], water mass !
   ! [ kg/m�] or [ kg/m�], and liquid fraction [---].                                      !
   !---------------------------------------------------------------------------------------!
   real(kind=8) function cmtl2uext8(dryhcap,wmass,temp,fliq)
      use consts_coms, only : cice8           & ! intent(in)
                            , cliq8           & ! intent(in)
                            , tsupercool_liq8 ! ! intent(in)

      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      real(kind=8), intent(in)  :: dryhcap ! Heat cap. of "dry" part   [J/m�/K] or [J/m�/K]
      real(kind=8), intent(in)  :: wmass   ! Mass                      [ kg/m�] or [ kg/m�]
      real(kind=8), intent(in)  :: temp    ! Temperature                           [     K]
      real(kind=8), intent(in)  :: fliq    ! Liquid fraction (0-1)                 [   ---]
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Internal energy is given by the sum of internal energies of dry part, plus the !
      ! contribution of ice and liquid phases.                                             !
      !------------------------------------------------------------------------------------!
      cmtl2uext8 = dryhcap * temp + wmass * ( (1.d0 - fliq) * cice8 * temp                 &
                                            + fliq * cliq8 * (temp - tsupercool_liq8) )
      !------------------------------------------------------------------------------------!

      return
   end function cmtl2uext8
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
   real(kind=8) function tq2enthalpy8(temp,humi,is_shv)
      use consts_coms, only : cpdry8          & ! intent(in)
                            , cph2o8          & ! intent(in)
                            , tsupercool_vap8 ! ! intent(in)

      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      real(kind=8), intent(in) :: temp   ! Temperature                             [     K]
      real(kind=8), intent(in) :: humi   ! Humidity (spec. hum. or mixing ratio)   [ kg/kg]
      logical     , intent(in) :: is_shv ! Input humidity is specific humidity     [   T|F]
      !----- Local variables. -------------------------------------------------------------!
      real(kind=8)             :: shv    ! Specific humidity                       [ kg/kg]
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Copy specific humidity to shv.                                                 !
      !------------------------------------------------------------------------------------!
      if (is_shv) then
         shv = humi
      else
         shv = humi / (humi + 1.d0)
      end if
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Enthalpy is the combination of dry and moist enthalpies, with the latter being !
      ! allowed to change phase.                                                           !
      !------------------------------------------------------------------------------------!
      tq2enthalpy8 = (1.d0 - shv) * cpdry8 * temp + shv * cph2o8 * (temp - tsupercool_vap8)
      !------------------------------------------------------------------------------------!

      return
   end function tq2enthalpy8
   !=======================================================================================!
   !=======================================================================================!







   !=======================================================================================!
   !=======================================================================================!
   !     This function computes the temperature [K] given the specific enthalpy and        !
   ! humidity.  If we assume that latent heat of vaporisation is a linear function of      !
   ! temperature (equivalent to assume that specific heats are constants and that the      !
   ! thermal expansion of liquid and water are negligible), then the saturation disappears !
   ! and the enthalpy becomes a straightforward state function.  In case you are looking   !
   ! at water exchange only, set the specific humidity to 1.0 and multiply the result by   !
   ! the water mass or water flux.                                                         !
   !---------------------------------------------------------------------------------------!
   real(kind=8) function hq2temp8(enthalpy,humi,is_shv)
      use consts_coms, only : cpdry8          & ! intent(in)
                            , cph2o8          & ! intent(in)
                            , tsupercool_vap8 ! ! intent(in)

      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real(kind=8), intent(in) :: enthalpy ! Specific enthalpy                     [  J/kg]
      real(kind=8), intent(in) :: humi     ! Humidity (spec. hum. or mixing ratio) [ kg/kg]
      logical     , intent(in) :: is_shv   ! Input humidity is specific humidity   [   T|F]
      !----- Local variables. -------------------------------------------------------------!
      real(kind=8)             :: shv      ! Specific humidity                     [ kg/kg]
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Copy specific humidity to shv.                                                 !
      !------------------------------------------------------------------------------------!
      if (is_shv) then
         shv = humi
      else
         shv = humi / (humi + 1.d0)
      end if
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Enthalpy is the combination of dry and moist enthalpies, with the latter being !
      ! allowed to change phase.                                                           !
      !------------------------------------------------------------------------------------!
      hq2temp8 = ( enthalpy + shv * cph2o8 * tsupercool_vap8 )                             &
               / ( (1.d0 - shv) * cpdry8 + shv * cph2o8 )
      !------------------------------------------------------------------------------------!

      return
   end function hq2temp8
   !=======================================================================================!
   !=======================================================================================!





   !=======================================================================================!
   !=======================================================================================!
   !     This function finds the latent heat of vaporisation for a given temperature.  If  !
   ! we use the definition of latent heat (difference in enthalpy between liquid and       !
   ! vapour phases), and assume that the specific heats are constants, latent heat becomes !
   ! a linear function of temperature.                                                     !
   !---------------------------------------------------------------------------------------!
   real(kind=8) function alvl8(temp)
      use consts_coms, only : alvl38  & ! intent(in)
                            , dcpvl8  & ! intent(in)
                            , t3ple8  ! ! intent(in)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      real(kind=8), intent(in) :: temp
      !------------------------------------------------------------------------------------!


      !----- Linear function, using latent heat at the triple point as reference. ---------!
      alvl8 = alvl38 + dcpvl8 * (temp - t3ple8)
      !------------------------------------------------------------------------------------!

      return
   end function alvl8
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This function finds the latent heat of sublimation for a given temperature.  If   !
   ! we use the definition of latent heat (difference in enthalpy between ice and vapour   !
   ! phases), and assume that the specific heats are constants, latent heat becomes a      !
   ! linear function of temperature.                                                       !
   !---------------------------------------------------------------------------------------!
   real(kind=8) function alvi8(temp)
      use consts_coms, only : alvi38  & ! intent(in)
                            , dcpvi8  & ! intent(in)
                            , t3ple8  ! ! intent(in)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      real(kind=8), intent(in) :: temp
      !------------------------------------------------------------------------------------!


      !----- Linear function, using latent heat at the triple point as reference. ---------!
      alvi8 = alvi38 + dcpvi8 * (temp - t3ple8)
      !------------------------------------------------------------------------------------!

      return
   end function alvi8
   !=======================================================================================!
   !=======================================================================================!



   !=======================================================================================!
   !=======================================================================================!
   !    This subroutine computes the temperature and fraction of liquid water from the     !
   ! intensive internal energy [J/kg].                                                     !
   !---------------------------------------------------------------------------------------!
   subroutine uint2tl8(uint,temp,fliq)
      use consts_coms, only : cliqi8          & ! intent(in)
                            , cicei8          & ! intent(in)
                            , allii8          & ! intent(in)
                            , t3ple8          & ! intent(in)
                            , uiicet38        & ! intent(in)
                            , uiliqt38        & ! intent(in)
                            , tsupercool_liq8 ! ! intent(in)
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real(kind=8), intent(in)  :: uint     ! Internal energy                     [   J/kg]
      real(kind=8), intent(out) :: temp     ! Temperature                         [      K]
      real(kind=8), intent(out) :: fliq     ! Liquid Fraction (0-1)               [    ---]
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Compare the internal energy with the reference values to decide which phase    !
      ! the water is.                                                                      !
      !------------------------------------------------------------------------------------!
      if (uint <= uiicet38) then
         !----- Internal energy below qwfroz, all ice  ------------------------------------!
         fliq = 0.d0
         temp    = uint * cicei8
         !---------------------------------------------------------------------------------!
      elseif (uint >= uiliqt38) then
         !----- Internal energy, above qwmelt, all liquid ---------------------------------!
         fliq = 1.d0
         temp    = uint * cliqi8 + tsupercool_liq8
         !---------------------------------------------------------------------------------!
      else
         !----- Changing phase, it must be at freezing point ------------------------------!
         fliq = (uint - uiicet38) * allii8
         temp    = t3ple8
         !---------------------------------------------------------------------------------!
      end if
      !------------------------------------------------------------------------------------!

      return
   end subroutine uint2tl8
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !    This subroutine computes the temperature (Kelvin) and liquid fraction from         !
   ! extensive internal energy (J/m� or J/m�), water mass (kg/m� or kg/m�), and heat       !
   ! capacity (J/m�/K or J/m�/K).                                                          !
   !---------------------------------------------------------------------------------------!
   subroutine uextcm2tl8(uext,wmass,dryhcap,temp,fliq)
      use consts_coms, only : cliqi8          & ! intent(in)
                            , cliq8           & ! intent(in)
                            , cicei8          & ! intent(in)
                            , cice8           & ! intent(in)
                            , allii8          & ! intent(in)
                            , alli8           & ! intent(in)
                            , t3ple8          & ! intent(in)
                            , tsupercool_liq8 ! ! intent(in)
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real(kind=8), intent(in)  :: uext    ! Extensive internal energy [  J/m�] or [  J/m�]
      real(kind=8), intent(in)  :: wmass   ! Water mass                [ kg/m�] or [ kg/m�]
      real(kind=8), intent(in)  :: dryhcap ! Heat cap. of "dry" part   [J/m�/K] or [J/m�/K]
      real(kind=8), intent(out) :: temp    ! Temperature                           [     K]
      real(kind=8), intent(out) :: fliq    ! Liquid fraction (0-1)                 [   ---]
      !----- Local variable ---------------------------------------------------------------!
      real(kind=8)              :: uefroz  ! qw of ice at triple pt.   [  J/m�] or [  J/m�]
      real(kind=8)              :: uemelt  ! qw of liq. at triple pt.  [  J/m�] or [  J/m�]
      !------------------------------------------------------------------------------------!



      !----- Convert melting heat to J/m� or J/m� -----------------------------------------!
      uefroz = (dryhcap + wmass * cice8) * t3ple8
      uemelt = uefroz   + wmass * alli8
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !    This is analogous to the uint2tl8 computation, we should analyse the magnitude  !
      ! of the internal energy to choose between liquid, ice, or both by comparing with    !
      ! the known boundaries.                                                              !
      !------------------------------------------------------------------------------------!
      if (uext < uefroz) then
         !----- Internal energy below qwfroz, all ice  ------------------------------------!
         fliq = 0.d0
         temp = uext  / (cice8 * wmass + dryhcap)
         !---------------------------------------------------------------------------------!
      elseif (uext > uemelt) then
         !----- Internal energy, above qwmelt, all liquid ---------------------------------!
         fliq = 1.d0
         temp = (uext + wmass * cliq8 * tsupercool_liq8) / (dryhcap + wmass * cliq8)
         !---------------------------------------------------------------------------------!
      elseif (uefroz == uemelt) then
         !---------------------------------------------------------------------------------!
         !    We are at the freezing point.  If water mass is so tiny that the internal    !
         ! energy of frozen and melted states are the same given the machine precision,    !
         ! then we assume that water content is negligible and we impose 50% frozen for    !
         ! simplicity.                                                                     !
         !---------------------------------------------------------------------------------!
         fliq = 5.d-1
         temp = t3ple8
         !---------------------------------------------------------------------------------!
      else
         !---------------------------------------------------------------------------------!
         !    Changing phase, it must be at freezing point.  The max and min are here just !
         ! to avoid tiny deviations beyond 0. and 1. due to floating point arithmetics.    !
         !---------------------------------------------------------------------------------!
         fliq = min(1.d0,max(0.d0,(uext - uefroz) * allii8 / wmass))
         temp = t3ple8
         !---------------------------------------------------------------------------------!
      end if
      !------------------------------------------------------------------------------------!

      return
   end subroutine uextcm2tl8
   !=======================================================================================!
   !=======================================================================================!
end module therm_lib8
!==========================================================================================!
!==========================================================================================!

