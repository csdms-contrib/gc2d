#!/usr/bin/env python

import sys
import getopt
import numpy
import time
import scipy
import logging
from scipy import interpolate
from scipy import signal
from scipy.io.numpyio import fwrite

# Available Mass Balance
class MassBalance:
   ( BAD_VAL , 
     ZERO_BALANCE ,
     CONSTANT_ELA ,
     ELA_LOWERING ,
     ELA_TIME_SERIES ,
     EXTERNAL_FUNC ,
     ELA_LOWERING2 ,
     BALANCE_FILE ,
     D180_TIME_SERIES ) = range( 9 )
        
class BoundaryCond:
   ( BAD_VAL ,
     ICE_FREE_BOUND ,
     ZERO_FLUX_BOUND ,
     CONST_FLUX_BOUND ,
     SURF_ELEV_BOUND ,
     SURF_SLOPE_BOUND ) = range( 6 )

class Parameters:
   g        = numpy.longdouble(9.81)                    # gravitional acceleration
   rhoI     = numpy.longdouble(917)                     # density of ice
   rhoW     = numpy.longdouble(1000)                    # density of water
   glensA   = numpy.longdouble( (6.8e-15)*3.15e7/(1e9) )   # Patterson, 1994; MacGregor, 2000
   day      = numpy.longdouble(0.00274)                 # length of a day in years

   # Time
   t         = numpy.longdouble(0)                  # set time to zero
   tMax      = numpy.longdouble(100000)             # maximum simulation time in years
   dtMax     = numpy.longdouble(0.4 * 365*day)      # maximum timestep in years
   dtDefault = numpy.longdouble(0.4 * 365*day)      # timestep if VARIABLE_DT_TOGGLE==0

   # Attractor Sliding -- only if ICESLIDE_TOGGLE==1 (generally used)
   UsChar   = numpy.longdouble(10)
   taubChar = numpy.longdouble(100000)

   # Glacier Properties
   MinGlacThick = numpy.longdouble(1)

   WEST_BC_TOGGLE  = BoundaryCond.ICE_FREE_BOUND
   EAST_BC_TOGGLE  = BoundaryCond.ICE_FREE_BOUND
   NORTH_BC_TOGGLE = BoundaryCond.ICE_FREE_BOUND
   SOUTH_BC_TOGGLE = BoundaryCond.ICE_FREE_BOUND

   MASS_BALANCE_TOGGLE = MassBalance.ELA_LOWERING    # select climate scenerio   (off|on|select)

   initELA         = numpy.longdouble(3000)
   gradBz          = numpy.longdouble(0.01)
   maxBz           = numpy.longdouble(2)
   ELAStepSize     = numpy.longdouble(-50)
   ELAStepInterval = numpy.longdouble(500)
   tmin            = numpy.longdouble(200)          # Years, spin-up time

   # Standard Sliding -- used if ICESLIDE_TOGGLE==2 (generally not used)
   B                 = numpy.longdouble(0.0012)     # m/(Pa*yr) -- MacGregor, 2000
   DepthToWaterTable = numpy.longdouble(20)         # distance from ice surface to water table
   MaxFloatFraction  = numpy.longdouble(80)         # limits water level in ice
   Hpeff             = numpy.longdouble(20)         # effective pressure (meters of water)

   # Avalanching
   angleOfRepose = numpy.longdouble(30)
   avalanchFreq  = numpy.longdouble(3)              # average number per year

   # Calving
   seaLevel    = numpy.longdouble(-100)             # meters
   calvingCoef = numpy.longdouble(2)                # year^-1

   # Thermal
   c      = numpy.longdouble(2060)                  # specific heat capacity (J/(kg*K))
   Qg     = numpy.longdouble(0.05*3.15e7)           # Geothermal heat flux (W/m^2)*seconds/year = (J/year)/(m^2)
   gradTz = numpy.longdouble(-0.0255)               # Geothermal Gradient


# Available Boundary Conditions
ICE_FREE_BOUND      = 1   # Ice Free Boundary
ZERO_FLUX_BOUND     = 2   # Zero Ice Flux
CONST_FLUX_BOUND    = 3   # Constant Ice Flux
SURF_ELEV_BOUND     = 4   # Constant Surface Elevation
SURF_SLOPE_BOUND    = 5   # Continuous Ice Surface Slope

#g    = numpy.longdouble(9.81)                    # gravitional acceleration
#rhoI = numpy.longdouble(917)                     # density of ice

#glensA = numpy.longdouble( (6.8e-15)*3.15e7/(1e9) )   # Patterson, 1994; MacGregor, 2000

# Attractor Sliding -- only if ICESLIDE_TOGGLE==1 (generally used)
#UsChar   = numpy.longdouble(10)
#taubChar = numpy.longdouble(100000)

# Glacier Properties
#MinGlacThick = numpy.longdouble(1)

#WEST_BC_TOGGLE  = ICE_FREE_BOUND
#EAST_BC_TOGGLE  = ICE_FREE_BOUND
#NORTH_BC_TOGGLE = ICE_FREE_BOUND
#SOUTH_BC_TOGGLE = ICE_FREE_BOUND

# Available Mass Balance
ZERO_BALANCE        = 1   # Constant Ice Flux
CONSTANT_ELA        = 2   # Ice Free Boundary
ELA_LOWERING        = 3   # Zero Ice Flux
ELA_TIME_SERIES     = 4   # Continuous Ice Surface Slope
EXTERNAL_FUNC       = 5   # Constant Surface Elevation
ELA_LOWERING2       = 6   # Zero Ice Flux
BALANCE_FILE        = 7   # Zero Ice Flux
D18O_TIME_SERIES    = 8   # Load d18O record and convert to ELA history
        
#MASS_BALANCE_TOGGLE = ELA_LOWERING    # select climate scenerio   (off|on|select)
        
#initELA         = numpy.longdouble(3000)
#gradBz          = numpy.longdouble(0.01)
#maxBz           = numpy.longdouble(2)
#ELAStepSize     = numpy.longdouble(-50)
#ELAStepInterval = numpy.longdouble(500)
#tmin            = numpy.longdouble(200)          # Years, spin-up time

def compress_grid( H , Zb , COMPRESS_TOGGLE=False , RESTART_TOGGLE=0 ):

   # COMPRESS - ONLY SIMULATE SUB-RECTANGLE THAT CONTAINS ICE
   if COMPRESS_TOGGLE and H.max() > 1 and RESTART_TOGGLE != 2:
      H_FullSpace  = H.copy()
      Zb_FullSpace = Zb.copy()
      if THERMAL_TOGGLE:
         Ts_FullSpace = Ts.copy()
         Tb_FullSpace = Tb.copy()
         Tm_FullSpace = Tm.copy()

      #[indrw,indcl] = find(H ~= 0);
      indrw,indcl   = numpy.where( H!=0 )

      mxrw,mxcl = Zb.shape

      mnrw = max( 0    , min(indrw) - 2 )
      mxrw = min( mxrw , max(indrw) + 2 )
      mncl = max( 0    , min(indcl) - 2 )
      mxcl = min( mxcl , max(indcl) + 2 )

      H  = H [ mnrw:mxrw , mncl:mxcl ]
      Zb = Zb[ mnrw:mxrw , mncl:mxcl ]
      Zi = Zb + numpy.choose( H<0 , (H,0) )
      #Zi = Zb + numpy.choose( numpy.less(H,0) , (H,0) )
      #Zi = Zb + max( H, 0 ) ;

      rws,cls = H.shape

      if THERMAL_TOGGLE:
         Ts = Ts[ mnrw:mxrw , mncl:mxcl ]
         Tb = Tb[ mnrw:mxrw , mncl:mxcl ]
         Tm = Tm[ mnrw:mxrw , mncl:mxcl ]

      mxrws,mxcls       = Zb_FullSpace.shape
      rws,cls           = Zb.shape
      compression_ratio = (mxcls*mxrws)/(cls*rws)
      COMPRESSED_FLAG   = 1
   else:
      #Zi = Zb + max( H, 0 ) # included for restarts
      Zi = Zb + numpy.choose( H<0 , (H,0) )
      compression_ratio = 1.
      COMPRESSED_FLAG   = 0

   return ( Zi , compression_ratio , COMPRESSED_FLAG )

def add_halo( x ):

   x_ext = numpy.concatenate( ( x[:,0,numpy.newaxis] , x     , x[:,-1,numpy.newaxis] ) , axis=1 )
   x_ext = numpy.concatenate( ( [x_ext[0,:]]         , x_ext , [x_ext[-1,:]]         ) )

   return x_ext

def set_bc( H , Zb , Zi , THERMAL_TOGGLE=False , WEST_BC_TOGGLE=ICE_FREE_BOUND , EAST_BC_TOGGLE=ICE_FREE_BOUND , SOUTH_BC_TOGGLE=ICE_FREE_BOUND , NORTH_BC_TOGGLE=ICE_FREE_BOUND ):

   ######
   ### MODIFY BOUNDARY CELLS TO ENFORCE BOUNDARY CONDITIONS

   # DEFAULT BOUNDARY CONDITION IS ZERO FLUX
   H_ext  = add_halo( H )
   Zb_ext = add_halo( Zb )
   Zi_ext = add_halo( Zi )
        
   if THERMAL_TOGGLE:
      Ts_ext = add_halo( Ts )
      Tb_ext = add_halo( Tb )
      Tm_ext = add_halo( Tm )

   # WESTERN BOUNDARY CONDTION
   if WEST_BC_TOGGLE == SURF_ELEV_BOUND:          # Constant Ice Surface Height
      ZiBound    = mean(Zb[:,0]) + Hbound
      H_ext[:,0] = ZiBound - Zb_ext[:,0]
   elif WEST_BC_TOGGLE == CONST_FLUX_BOUND:       # Constant Ice Flux B.C.
      pass
   elif WEST_BC_TOGGLE == SURF_SLOPE_BOUND:       # Constant Ice Surface Slope
      Zi_ext[:,0] = 2*Zi_ext[:,1] - Zi_ext[:,2]
      H_ext [:,0] = Zi_ext[:,0] - Zb_ext[:,0]
      H_ext [:,0] = numpy.choose( H_ext[:,0]<0 , (H_ext[:,0],0) )
   elif WEST_BC_TOGGLE == ICE_FREE_BOUND:         # Ice Free Boundary
      H_ext[:,0] = 0

   # EASTERN BOUNDARY CONDTION
   if EAST_BC_TOGGLE == SURF_ELEV_BOUND:          # Constant Ice Surface Height
      ZiBound     = mean(Zb[:,-1]) + Hbound
      H_ext[:,-1] = ZiBound - Zb_ext[:,-1]
   elif EAST_BC_TOGGLE == CONST_FLUX_BOUND:       # Constant Ice Flux B.C.
      pass
   elif EAST_BC_TOGGLE == SURF_SLOPE_BOUND:       # Constant Ice Surface Slope
      Zi_ext[:,-1] = 2*Zi_ext[:,-2] - Zi_ext[:,-3]
      H_ext [:,-1] = Zi_ext[:,-1] - Zb_ext[:,-1]
      H_ext [:,-1] = numpy.choose( H_ext[:,-1]<0 , (H_ext[:,-1],0) )
   elif EAST_BC_TOGGLE == ICE_FREE_BOUND:         # Ice Free Boundary
      H_ext[:,-1] = 0
            
   # SOUTHERN BOUNDARY CONDTION
   if SOUTH_BC_TOGGLE == SURF_ELEV_BOUND:         # Constant Ice Surface Height
      ZiBound    = mean(Zb[0,:]) + Hbound
      H_ext[0,:] = ZiBound - Zb_ext[0,:]
   elif SOUTH_BC_TOGGLE == CONST_FLUX_BOUND:      # Constant Ice Flux B.C.
      pass
   elif SOUTH_BC_TOGGLE == SURF_SLOPE_BOUND:      # Constant Ice Surface Slope
      Zi_ext[0,:] = 2*Zi_ext[1,:] - Zi_ext[2,:]
      H_ext [0,:] = Zi_ext[0,:] - Zb_ext[0,:]
      H_ext [0,:] = numpy.choose( H_ext[0,:]<0 , (H_ext[0,:],0) )
   elif SOUTH_BC_TOGGLE == ICE_FREE_BOUND:        # Ice Free Boundary
      H_ext[0,:] = 0
            
   # NORTHERN BOUNDARY CONDTION
   if NORTH_BC_TOGGLE == SURF_ELEV_BOUND:         # Constant Ice Surface Height
      ZiBound     = mean(Zb[-1,:]) + Hbound
      H_ext[-1,:] = ZiBound - Zb_ext[-1,:]
   elif NORTH_BC_TOGGLE == CONST_FLUX_BOUND:      # Constant Ice Flux B.C.
      pass
   elif NORTH_BC_TOGGLE == SURF_SLOPE_BOUND:      # Constant Ice Surface Slope
      Zi_ext[-1,:] = 2*Zi_ext[-2,:] - Zi_ext[-3,:]
      H_ext [-1,:] = Zi_ext[-1,:] - Zb_ext[-1,:]
      H_ext [-1,:] = numpy.choose( H_ext[-1,:]<0 , (H_ext[-1,:],0) )
   elif NORTH_BC_TOGGLE == ICE_FREE_BOUND:        # Ice Free Boundary
      H_ext[-1,:] = 0
        
   Zi_ext = Zb_ext + H_ext

   return ( H_ext , Zb_ext , Zi_ext )

def difference_grid( A , dx , dy ):

   dAdx_ext   = ( A[:,1:] - A[:,:-1] ) / dx
   dAdy_ext   = ( A[1:,:] - A[:-1,:] ) / dy
   dAdx       = dAdx_ext[1:-1,:]
   dAdy       = dAdy_ext[:,1:-1]

   return ( dAdx , dAdy )

def basal_shear_stress( H_ext , Zi_ext , dx=1. , dy=1. , g=Parameters.g , rhoI=Parameters.rhoI ):

   ######
   ### CALCULATE THE BASAL SHEAR STRESS
            
   # forward differences
   dZidxX_ext = ( Zi_ext[:,1:] - Zi_ext[:,:-1] ) / dx
   dZidyY_ext = ( Zi_ext[1:,:] - Zi_ext[:-1,:] ) / dy

   dZidxX     = dZidxX_ext[1:-1,:]
   dZidyY     = dZidyY_ext[:,1:-1]

   HX_ext     = ( H_ext[:,1:] + H_ext[:,:-1] ) / 2.
   HY_ext     = ( H_ext[1:,:] + H_ext[:-1,:] ) / 2.
   HX         = HX_ext[1:-1,:]
   HY         = HY_ext[:,1:-1]
            
   taubxX_ext = -rhoI * g * HX_ext * dZidxX_ext
   taubyY_ext = -rhoI * g * HY_ext * dZidyY_ext

   taubxX     = taubxX_ext[1:-1,:]
   taubyY     = taubyY_ext[:,1:-1]

   taubxY = ( taubxX_ext[:-1,:-1] + taubxX_ext[:-1,1:] + 
              taubxX_ext[1: ,:-1] + taubxX_ext[1: ,1:] ) / 4.
            
   taubyX = ( taubyY_ext[:-1,:-1] + taubyY_ext[:-1,1:] +
              taubyY_ext[1: ,:-1] + taubyY_ext[1: ,1:] ) / 4.
            
   taubX  = numpy.sqrt( taubxX**2 + taubyX**2 )
   taubY  = numpy.sqrt( taubxY**2 + taubyY**2 )

   taubX  = numpy.choose( HX>0 , (0,taubX) )
   taubY  = numpy.choose( HY>0 , (0,taubY) )

   # fill in zero values with 1 for use in division
   xcmpnt = numpy.choose( numpy.abs(taubX)<1e-5 , ( taubxX / taubX , 0. ) )
   ycmpnt = numpy.choose( numpy.abs(taubY)<1e-5 , ( taubyY / taubY , 0. ) )

   return ( ( xcmpnt , ycmpnt ) , ( taubX , taubY ) , ( HX , HY ) )

def iceflow( taubX , taubY , HX , HY , xcmpnt , ycmpnt , THERMAL_TOGGLE=False , MinGlacThick=1. , glensA=Parameters.glensA ):

   ######
   ### CALCULATE ICE VELOCITY DUE TO DEFORMATION
        
   if THERMAL_TOGGLE:
      A_ext = numpy.zeros(H_ext.shape , dtype=numpy.longdouble )
      ind = nonzero( ravel(H_ext) >= MinGlacThick )

      Ts_ext = To + lapseRate*( Zi_ext - Elev0 )

      #A_ext(ind) = interp3( eHs, eTs, eTm, eA, H_ext(ind), Ts_ext(ind), Tm_ext(ind) ) ;
      try:
         put( A_ext , ind , interpolate.interp3d( eHs , eTs , eTm )( take(H_ext,ind) , take(Ts_ext,ind) , take(Tm_ext,ind) ) )
      except:
         logging.error( "NaN in A, likely H_node exceeds H_glens limits" )
         return -1
            
      AX = ( A_ext[1:-1, :-1] + A_ext[1:-1,1:  ] ) / 2.
      AY = ( A_ext[ :-1,1:-1] + A_ext[1:  ,1:-1] ) / 2.

   else:
      AX = glensA
      AY = glensA

   # here's the guts of calculating the depth averaged velocity
   UdxX = numpy.abs( .4 * AX * taubX*taubX*taubX * HX ) * xcmpnt
   UdyY = numpy.abs( .4 * AY * taubY*taubY*taubY * HY ) * ycmpnt

   #UdxX = numpy.fix(UdxX*1e6)*1e-6
   #UdyY = numpy.fix(UdyY*1e6)*1e-6

   return ( UdxX , UdyY )

def ice_sliding( taubX , taubY , xcmpnt , ycmpnt , THERMAL_TOGGLE=False , FREEZEON_TOGGLE=0 , UsChar=Parameters.UsChar , taubChar=Parameters.taubChar ):

   ######
   ### CALCULATE SLIDING VELOCITY
        
   # here's the guts of calculating the sliding velocity 
   UsxX = numpy.choose( numpy.abs(taubX)<1e-5 , ( UsChar * numpy.exp(1 - taubChar / taubX) * xcmpnt ,
                                                  UsChar * numpy.exp(1 - taubChar        ) * xcmpnt ) )
   UsyY = numpy.choose( numpy.abs(taubY)<1e-5 , ( UsChar * numpy.exp(1 - taubChar / taubY) * ycmpnt , 
                                                  UsChar * numpy.exp(1 - taubChar        ) * ycmpnt ) )

   if THERMAL_TOGGLE and FREEZEON_TOGGLE:
      notFrozen  = Tb_ext > -.5 or Zb_ext < seaLevel
      notFrozenX = ( notFrozen[1:-1, :-1] + notFrozen[1:-1,1:  ] ) / 2.
      notFrozenY = ( notFrozen[ :-1,1:-1] + notFrozen[1:  ,1:-1] ) / 2.

      UsxX *= notFrozenX
      UsyY *= notFrozenY

   return ( UsxX , UsyY )

def sum_ice_motion( UdxX , UdyY , UsxX , UsyY ):

   UxX = UdxX + UsxX
   UyY = UdyY + UsyY

   return ( UxX , UyY )

def mass_conservation( H_ext , UxX , UyY , HX , HY , dZidxX , dZidyY , dx=1. , dy=1. , MinGlacThick=1. , WEST_BC_TOGGLE=ICE_FREE_BOUND , EAST_BC_TOGGLE=ICE_FREE_BOUND , SOUTH_BC_TOGGLE=ICE_FREE_BOUND , NORTH_BC_TOGGLE=ICE_FREE_BOUND ):

   ######
   ### MASS CONSERVATION -- CONTINUITY

   # ensure that no ice is drawn from the rock
   #CLASS = H_ext >= MinGlacThick
   CLASS = numpy.choose( H_ext>=MinGlacThick , (0.,1.) )
            
   DCLASSx = ( CLASS[1:-1,1:  ] - CLASS[1:-1, :-1] ) * numpy.sign( dZidxX )
   DCLASSy = ( CLASS[1:  ,1:-1] - CLASS[ :-1,1:-1] ) * numpy.sign( dZidyY )
            
   UxX = numpy.choose( numpy.abs(DCLASSx+1)<1e-5 , (UxX,0.) )
   UyY = numpy.choose( numpy.abs(DCLASSy+1)<1e-5 , (UyY,0.) )

   # calculate both components of the ice flux
   qxX = UxX * HX
   qyY = UyY * HY
            
   if WEST_BC_TOGGLE  == CONST_FLUX_BOUND: qxX[: , 0] = BoundaryFlux
   if EAST_BC_TOGGLE  == CONST_FLUX_BOUND: qxX[: ,-1] = BoundaryFlux
   if SOUTH_BC_TOGGLE == CONST_FLUX_BOUND: qyY[0 , :] = BoundaryFlux
   if NORTH_BC_TOGGLE == CONST_FLUX_BOUND: qyY[-1, :] = BoundaryFlux
            
   # here's the guts of the continuity equation
   dqdxX = ( qxX[ :,1:] - qxX[:  ,:-1] ) / dx
   dqdyY = ( qyY[1:, :] - qyY[:-1,:  ] ) / dy
   dHdt  = -dqdxX - dqdyY

   return ( dHdt , ( qxX , qyY ) )

def mass_balance( Zi , t , MASS_BALANCE_TOGGLE=MassBalance.ELA_LOWERING , initELA=Parameters.initELA , tmin=Parameters.tmin , ELAStepSize=Parameters.ELAStepSize , ELAStepInterval=Parameters.ELAStepInterval , gradBz=Parameters.gradBz , maxBz=Parameters.maxBz ):

   ######
   ### CALCULATE MASS BALANCE
            
   # the imposed mass balance is the imposed climate
   # there are many possibilities, here are only a few
   # all must populate the 2D matrix Bxy of size = size(Zb)
   # with values of net precip/melt rate in m/yr
   # define the scalar, ELA (m), for plotting

   if MASS_BALANCE_TOGGLE == CONSTANT_ELA:
      # Simple ELA, maxBz, gradBz

      ELA = initELA
      #Bxy = min( maxBz , gradBz * ( Zi - ELA ) )
      Bxy = gradBz * ( Zi - ELA )
      Bxy = numpy.choose( Bxy>maxBz , (Bxy,maxBz) )
            
   elif MASS_BALANCE_TOGGLE == ELA_LOWERING:
      # ELA changing with time experiment
            
      # ELAStepSize = -10 ;       # positive/negative values raise/lower ELA
      # ELAStepInterval = 500 ;
                
      ELA = initELA + ELAStepSize * max( 0 , numpy.floor( (t-tmin)/ELAStepInterval ) )
      Bxy = gradBz * ( Zi - ELA )
      Bxy = numpy.choose( Bxy>maxBz , (Bxy,maxBz) )

   elif MASS_BALANCE_TOGGLE == ELA_LOWERING2:
      # ELA changing with time experiment
            
      tau      = numpy.longdouble(25)          # intrinsic timescale of ice dynamics 
      tmin     = numpy.longdouble(0)           # time to begin ELA modification
      initELA  = numpy.longdouble(4200)        # initial ELA
      stepSize = numpy.longdouble(-10)         # positive/negative values raise/lower ELA
      dELAdt   = numpy.longdouble(-0.1) 
                
      ELA = initELA + stepSize * max( 0, numpy.floor( (t-tmin) / (8*tau) ) )
      Bxy = gradBz * ( Zi - ELA )
      Bxy = numpy.choose( Bxy>maxBz , (Bxy,maxBz) )
            
   elif MASS_BALANCE_TOGGLE == EXTERNAL_FUNC:
      # external mass balance function
      try: Bxy
      except NameError:
         # Mass Balance 2D Must Return Bxy (2d Matrix)
         Bxy = mass_balance_gc2d( t , cellsize , Zi )
         nextGetBxy = t + getBxyInterval
      else:
         if t >= nextGetBxy:
            Bxy = mass_balance_gc2d( t , cellsize , Zi )
            nextGetBxy = t + getBxyInterval

   elif MASS_BALANCE_TOGGLE == ELA_TIME_SERIES or MASS_BALANCE_TOGGLE == D18O_TIME_SERIES:
      # ELA time series
      ELA = interpolate.interp1d( trecord , ELArecord )( t )
      Bxy = gradBz * ( Zi - ELA )
      Bxy = numpy.choose( Bxy>maxBz , (Bxy,maxBz) )

   elif MASS_BALANCE_TOGGLE == BALANCE_FILE:
      # external mass balance file
      Bxy = load_dem_var( filenameDEM, 'Bxy' )
      ind = nonzero( ravel(abs(Bxy)==min(abs(Bxy))) )
      ELA = mean( take( ravel(Zi) , ind ) )

   elif MASS_BALANCE_TOGGLE == ZERO_BALANCE:
      ELA = 0
      Bxy = numpy.zeros( Zb.shape , dtype=numpy.longdouble )
                
   else:
      logging.error( "Unrecognized Mass Balance" )
      return -1

   return ( Bxy , ELA )
            
def get_timestep( H , Zi_ext , Zi , dHdt , Bxy , dtMax=None , dtDefault=None ):

   #######
   ### CALCULATE TIMESTEP
            
   # now that we know the rate of change in ice surface heights due to  
   # ice motion and due to precipitation or melt we need to know over 
   # what period of time we can project forward with these rates and 
   # maintain stability of the ice surface.  The basic idea here is that
   # we don't want to take a timestep any longer then it would take to 
   # reverse the ice surface slope between two cells, such that ice 
   # should be flowing in the other direction.  In fact, let's make our 
   # timestep much less then that.
            
   # this calculation sets the timestep such that the change
   # in ice surface elevation nowhere exceeds a set fraction
   # of the local standard deviation in ice surface elevations
            
   # include ice changes by precip and melt
   dHdtTot = dHdt + Bxy
   adHdt   = numpy.abs(dHdtTot)
            
   # something like standard deviation of 3x3 cell areas around each cell
   filt   = numpy.ones( (3,3) , dtype=numpy.longdouble ) / 9.
   ZiMean = filter2d( filt , Zi_ext , 'valid' )
   dHmax  = numpy.sqrt( filter2d( filt, (ZiMean - Zi)**2 ) )
            
   # only consider cells with ice thickness > 10 m
   isGlac = H>10.
            
   # find limiting timestep for each considered cell
   ind = ( numpy.logical_and( numpy.logical_and( adHdt!=0 , dHmax!=0 ) , isGlac!=0 ) ).flatten().nonzero()

   if ind[0].size>0:
      dtLimits = dHmax.flatten()[ ind ] / adHdt.flatten()[ ind ]
      dt       = dtLimits.min()
      idt      = ( dtLimits==dt ).nonzero()

      #ind = find( adHdt~=0 & dHmax~=0 & isGlac~=0 ) ;    
      #dtLimits = dHmax(ind)./adHdt(ind) ;
      #[dt, idt] = min( dtLimits ) ;
            
      # locate the x and y position of limiting cell for plotting
      #[rwDT,clDT] = ind2sub( size(adHdt), ind(idt) ) ; 
            
      # limit timestep to dtMax or some fraction of the calculated timestep
      if dtMax is not None :
         dt = min( dtMax, dt/2. )
            
   else:
      # catch an error, (e.g. if H<10 in all cells )
      #if dt.size==0:
      dt = dtDefault

   #dt = numpy.fix(dt*1e6)*1e-6

   return dt

def update_vars( H , Zb , Zi , Bxy , qxX , qyY , dHdt , t , dt , conserveIce , dx=1. , dy=1. ):

   t = t + dt

   # numTimeSteps = numTimeSteps + 1 ;
   # timeSteps(numTimeSteps) = dt ;
            
   # increase in ice thicknesses due to precip
   Bxy_pos  = numpy.choose( Bxy>0 , (0,Bxy) )
   H       += Bxy_pos * dt
            
   # change ice thicknesses due to ice motion
   H       += dHdt*dt
            
   # decrease in ice thicknesses due to melt
   Bxy_neg  =   numpy.choose( Bxy<0 , (0,Bxy) )
   Bxy_neg  = - numpy.choose( H<-Bxy_neg , (-Bxy_neg,H) )
   H       += Bxy_neg * dt
            
   # record ice addition or removal by climate
   snowFall    = ( Bxy_neg + Bxy_pos ) * dt
   conserveIce = conserveIce + snowFall.sum(axis=0).sum()
            
   # record ice flux through boundaries
   qbound = qyY[0,:].sum(axis=0).sum() - qyY[-1,:].sum(axis=0).sum() + qxX[:,0].sum(axis=0).sum() - qxX[:,-1].sum(axis=0).sum()
   conserveIce = conserveIce + dt * qbound / dx
            
   Zi = Zb + numpy.choose( H<0 , (H,0) )
        
   if numpy.isnan(Zi).any():
      #save workspacedump
      logging.error( "NaN in ice thickness" )
      return -1

   return ( t , H , Zi , conserveIce )
            
def avalanche( H , angleOfRepose=30. ):

   ######
   ### AVALANCHE SNOW OFF OF STEEP SURFACES

   # move ice downslope until the ice surface is everywhere
   # less then or near the angle of repose
            
   rws,cls  = Zb.shape
   dHRepose = dx*numpy.tan(angleOfRepose*numpy.pi/180.)
   Ho       = numpy.choose( H<0 , (H,0) )
            
   while True:
      dZidx_down        = numpy.zeros( (rws,cls) , dtype=numpy.longdouble )
      dZidx_up          = numpy.zeros( (rws,cls) , dtype=numpy.longdouble )
      dZidx_down[:,1:]  = numpy.choose( Zi[:,1:]  < Zi[:,:-1] , ( Zi[:,1:]  - Zi[:,:-1] , 0 ) )
      dZidx_up  [:,:-1] = numpy.choose( Zi[:,:-1] < Zi[:,1:]  , ( Zi[:,:-1] - Zi[:,1:]  , 0 ) )
      dZidx             = numpy.choose( dZidx_up > dZidx_down , ( dZidx_down , dZidx_up ) )

      dZidy_left         = numpy.zeros( (rws,cls) , dtype=numpy.longdouble )
      dZidy_right        = numpy.zeros( (rws,cls) , dtype=numpy.longdouble )
      dZidy_left [1:,:]  = numpy.choose( Zi[1:,:] < Zi[:-1,:] , ( Zi[1:,:] - Zi[:-1,:] , 0 ) )
      dZidy_right[:-1,:] = numpy.choose( Zi[:-1,:] < Zi[1:,:] , ( Zi[:-1,:] - Zi[1:,:] , 0 ) )
      dZidy              = numpy.choose( dZidy_left > dZidy_right , ( dZidy_right , dZidy_left ) )

      grad  = numpy.sqrt( dZidx**2 + dZidy**2 )
      gradT = dZidy_left + dZidy_right + dZidx_down + dZidx_up
      gradT = numpy.choose( gradT==0   , (gradT,1) )
      grad  = numpy.choose( Ho    <0.1 , (grad ,0) )

      mxGrad = grad.max()
        
      if mxGrad <= 1.1*dHRepose:
         break

      delH = numpy.choose( grad<dHRepose , ( ( grad-dHRepose)/3. , 0 ) )
        
      Htmp = Ho.copy()
      Ho   = numpy.choose( Htmp<delH , ( Htmp-delH , 0 ) )
      delH = Htmp - Ho
        
      delHdn = numpy.zeros( (rws,cls) , dtype=numpy.longdouble )
      delHup = numpy.zeros( (rws,cls) , dtype=numpy.longdouble )
      delHlt = numpy.zeros( (rws,cls) , dtype=numpy.longdouble )
      delHrt = numpy.zeros( (rws,cls) , dtype=numpy.longdouble )
        
      delHup[:,1:  ] = delH[:, :-1] * dZidx_up  [:, :-1] / gradT[:, :-1]
      delHdn[:, :-1] = delH[:,1:  ] * dZidx_down[:,1:  ] / gradT[:,1:  ]
      delHrt[1:  ,:] = delH[ :-1,:] * dZidy_right[ :-1,:] / gradT[ :-1,:]
      delHlt[ :-1,:] = delH[1:  ,:] * dZidy_left [1:  ,:] / gradT[1:  ,:]
        
      Ho = Ho + delHdn + delHup + delHlt + delHrt
      Ho = numpy.choose( Ho<0 , (Ho,0) )
        
      Zi = Zb + Ho
            
   #H = Ho + (H<0).*H ;
   H = Ho + numpy.choose( H<0 , (0,H) )

   return H

def calve( H , dt , CALVING_TOGGLE=True ):

      ######
      ### CALVING GLACIER FRONT
        
      if CALVING_TOGGLE:
         # one reason this is difficult is that the height of ice in the cell
         # is really just recording the volume of ice, the position of the 
         # margin in the cell not the actual ice height.  Here floation
         # height is assumed (or higher if necessary to account for ice volume)
            
         Hold      = H.copy()
         calvedIce = 0
            
         # count time backwards with a sshorted timestep until the whole 
         # timestep used during this itteration has been simulated
            
         dtTot = dt
         while dtTot > 0:
            # find the calving front, aka the wet glacier margin
            G = H > 1
            W = numpy.logical_and( G==0 , Zb <= seaLevel )
            filt = numpy.array( [[0,1,0],[1,1,1],[0,1,0]] , dtype=numpy.longdouble )
            Wfilt = filter2d( filt , W )
            Wfilt[:,(0,-1)] = Wfilt[:,(2,-3)]
            Wfilt[(0,-1),:] = Wfilt[(2,-3),:]
            wetGmargin = Gi * Wfilt > 0
            indWGM = wetGmargin.ravel().nonzero()
                
            # if calving front exists, find water depth, ensure it's positive
            if indWGM.size>0:
               WDmarg = seaLevel - Zb.flatten()[indWGM]
               WDmarg = numpy.choose( WDmarg<0 , (WDmarg,0) )
               ind    = (WDmarg!=0).nonzero()
               indWGM = take( indWGM , ind )
               WDmarg = take( WDmarg , ind )

               #WDmarg = max( 0, seaLevel - Zb(indWGM) ) ;
               #ind = find( WDmarg == 0 ) ;
               #indWGM(ind) = [] ;
               #WDmarg(ind) = [] ;
                
            # if calving front exists, remove some ice
            if indWGM.size>0:
               # ice thickness in calving cells
               Hmarg = H.flatten()[indWGM]
               Hmarg = numpy.choose( Hmarg<WDmarg/0.917 , (Hmarg,WDmarg/0.917) )
                    
               # a new timestep is calculated such that the calving rate times the 
               # timesstep does not exceed the total contents of any calving cell.
                    
               dLinCalvdt   = calvingCoef * WDmarg                         # front migration rate
               dVolCalvdt   = dx * dLinCalvdt * Hmarg                      # rate of volume calved
               volAvailMarg = dx * dx * H.flatten()[indWGM]                # ice volume available
               calvDt = min( dtTot, ( volAvailMarg / dVolCalvdt ).min() )  # calving timestep

               # remove this calving timestep from total time to calve
               dtTot = dtTot - calvDt
                    
               # convert the volume calved to ice thickness and remove
               calve = dVolCalvdt * calvDt / ( dx * dx )
               H[indWGM] = H[indWGM] - calve
                
               # record total volume calved for posterity
               calvedIce = calvedIce + calve.sum(asis=0).sum() * dx * dx
                    
            else:
               dtTot = 0
                
         # record ice removal by calving for conservation test
         conserveIce = conserveIce + ( H - Hold ).sum(axis=0).sum()

def print_watch_point( fd , x ):
   y = numpy.double( x )
   fwrite( fd , y.size , y )
   fd.flush()

def filter2d( b , x , shape='same' ):
   y = scipy.signal.convolve( b ,x , mode=shape )
   return y

def load_dem( file ):

   vars = scipy.io.loadmat( file )

   cellsize = numpy.longdouble(vars['cellsize'])
   easting  = numpy.longdouble(vars['easting'])
   northing = numpy.longdouble(vars['northing'])
   topo     = numpy.longdouble(vars['topo'])

   n_rows , n_cols = topo.shape

   logging.info( 'Shape of topo is %d by %d' , n_rows , n_cols )
   logging.info( 'Shape of easting is %d'    , easting.size )
   logging.info( 'Shape of northing is %d'   , northing.size )

   if easting.size != n_cols:
      sys.exit( 'Easting does not match dimension of topo (%d != %d)' % (easting.size,n_cols) )
   if northing.size != n_rows:
      sys.exit( 'Northing does not match dimension of topo (%d != %d)' % (northing.size,n_rows) )

   return ( topo , easting , northing , cellsize )

def load_dem_var( file , val_s ):
   vars = scipy.io.loadmat( file )

   if vars.has_key( val_s ):
      var = vars[val_s]
   else:
      var = None

   return var

def load_input_args( ):
   CLEAR_FIGURE          = 1
   CONTOUR_INTERVAL      = 50.
   DEBUG_TOGGLE          = 0
   DT_LIMIT              = 0
   ELA_CONTOUR           = 1.
   ICE_CONTOUR           = 1.
   NEW_FIGURE            = 0
   QUIVER_VECS           = 0
   RECONSTRUCT           = 0
   SUBFIGURE             = 0
   THERMAL_CONTOUR       = 0

   return 1

class Usage( Exception ):
   def __init__( self , msg ):
      self.msg = msg

def old_gc2d( argv=None , inputFile='Animas_200.mat' ):

   if argv is None:
      argv = sys.argv
   try:
      try:
         opts, args = getopt.getopt( argv[1:] , "h" , ["help"] )
      except getopt.error, msg:
         raise Usage(msg)
   except Usage, err:
      print >> sys.stderr, err.msg
      print >> sys.stderr, "for help use --help"
      return 2

   RESTART_TOGGLE = 0

   ######
   ### Load a saved state

   if RESTART_TOGGLE == 0  or  RESTART_TOGGLE == 3: # LOAD A SAVED STATE

      # CODE BEHAVIOR TOGGLES
      # toggles turn on/off segments of the code or select 
      # between multiple possibilities for a given process
      # values can be reset in INIT_COND segment
        
      GUISTART_TOGGLE     = 0   # started simulation with the gui   (off|on)
        
      SAVE_TOGGLE         = 1   # saving                            (off|on)
      PLOT_TOGGLE         = 1   # plotting                          (off|on)
      REPORT_TOGGLE       = 1   # reporting                         (off|on)
        
      COMPRESS_TOGGLE     = 0   # only simulate area with ice       (off|on)
      VARIABLE_DT_TOGGLE  = 1   # state dependent time step         (off|on)

      INIT_COND_TOGGLE    = 1   # load DEM and climate              (synth|valley|sheet)
      GENERIC_ICE_TOGGLE  = 0   # start with generic ice surface    (off|on)       
        
      ICEFLOW_TOGGLE      = 1   # ice motion by deformation         (off|on)
      ICESLIDE_TOGGLE     = 0   # ice motion by sliding             (off|on|select)        
        
      THERMAL_TOGGLE      = 0   # temp dependance of flow           (off|on)
      FREEZEON_TOGGLE     = 0   # basal ice freeze to bed           (off|on)

      AVALANCH_TOGGLE     = 0   # avalanch off steep surfaces       (off|on)
      ERODE_TOGGLE        = 0   # erode the bed                     (off|on|select)
      CALVING_TOGGLE      = 0   # calving front                     (off|on)
        
      CRN_TOGGLE          = 0   # CRN accumulation                  (off|on)
        
        
      # Available Mass Balance
      ZERO_BALANCE        = 1   # Constant Ice Flux
      CONSTANT_ELA        = 2   # Ice Free Boundary
      ELA_LOWERING        = 3   # Zero Ice Flux
      ELA_TIME_SERIES     = 4   # Continuous Ice Surface Slope
      EXTERNAL_FUNC       = 5   # Constant Surface Elevation
      ELA_LOWERING2       = 6   # Zero Ice Flux
      BALANCE_FILE        = 7   # Zero Ice Flux
      D18O_TIME_SERIES    = 8   # Load d18O record and convert to ELA history
        
      MASS_BALANCE_TOGGLE = ELA_LOWERING    # select climate scenerio   (off|on|select)
        
      # Available Boundary Conditions
      ICE_FREE_BOUND      = 1   # Ice Free Boundary
      ZERO_FLUX_BOUND     = 2   # Zero Ice Flux
      CONST_FLUX_BOUND    = 3   # Constant Ice Flux
      SURF_ELEV_BOUND     = 4   # Constant Surface Elevation
      SURF_SLOPE_BOUND    = 5   # Continuous Ice Surface Slope
        
      WEST_BC_TOGGLE = ICE_FREE_BOUND   # boundary condition    (no ice|reflect|no flow)
      EAST_BC_TOGGLE = ICE_FREE_BOUND   # boundary condition    (no ice|reflect|no flow)
      SOUTH_BC_TOGGLE = ICE_FREE_BOUND  # boundary condition    (no ice|reflect|no flow)
      NORTH_BC_TOGGLE = ICE_FREE_BOUND  # boundary condition    (no ice|reflect|no flow)
        
        
      # OUTPUT BEHAVIOR
      plotInterval = 60 * 120       # seconds
      saveInterval = 100            # whole years
      reportInterval = 30           # seconds
        
      nextPlot = 0                  # initialize to plot on first timestep
      nextSave = 0                  # initialize to save on first timestep
      nextReport = 0                # initialize to report on first timestep
        
      outputFile = 'savetmp'

      ######
      ### Set numerical and physical constants
      # Constants
      g    = numpy.longdouble(9.81)                    # gravitional acceleration
      rhoI = numpy.longdouble(917)                     # density of ice
      rhoW = numpy.longdouble(1000)                    # density of water
      day  = numpy.longdouble(0.00274)                 # length of a day in years

      # Time
      t         = numpy.longdouble(0)                  # set time to zero
      tMax      = numpy.longdouble(100000)             # maximum simulation time in years
      dtMax     = numpy.longdouble(0.4 * 365*day)      # maximum timestep in years
      dtDefault = numpy.longdouble(0.4 * 365*day)      # timestep if VARIABLE_DT_TOGGLE==0

      # Glacier Properties
      MinGlacThick = numpy.longdouble(1)

      # Ice Deformation
      glensA = numpy.longdouble((6.8e-15)*3.15e7/(1e9))    # Patterson, 1994; MacGregor, 2000


      # Attractor Sliding -- only if ICESLIDE_TOGGLE==1 (generally used)
      UsChar   = numpy.longdouble(10)
      taubChar = numpy.longdouble(100000)

      # Standard Sliding -- used if ICESLIDE_TOGGLE==2 (generally not used)
      B                 = numpy.longdouble(0.0012)     # m/(Pa*yr) -- MacGregor, 2000
      DepthToWaterTable = numpy.longdouble(20)         # distance from ice surface to water table
      MaxFloatFraction  = numpy.longdouble(80)         # limits water level in ice
      Hpeff             = numpy.longdouble(20)         # effective pressure (meters of water)

      # Avalanching
      angleOfRepose = numpy.longdouble(30)
      avalanchFreq  = numpy.longdouble(3)              # average number per year

      # Calving
      seaLevel    = numpy.longdouble(-100)             # meters
      calvingCoef = numpy.longdouble(2)                # year^-1


      # Thermal
      c      = numpy.longdouble(2060)                  # specific heat capacity (J/(kg*K))
      Qg     = numpy.longdouble(0.05*3.15e7)           # Geothermal heat flux (W/m^2)*seconds/year = (J/year)/(m^2)
      gradTz = numpy.longdouble(-0.0255)               # Geothermal Gradient


      # Mass Balance
      initELA         = numpy.longdouble(4500)
      initELA         = numpy.longdouble(3000)
      gradBz          = numpy.longdouble(0.01)
      maxBz           = numpy.longdouble(2)
      ELAStepSize     = numpy.longdouble(-50)
      ELAStepInterval = numpy.longdouble(500)
      tmin            = numpy.longdouble(200)         # Years, spin-up time

      ######
      ### RELOAD INPUT ARGUMENTS
    
      #load inputArgs
      inputArgs = load_input_args

      #if ( GUISTART_TOGGLE & exist('guiSimParams.mat','file') )
      #   load guiSimParams
      #   delete guiSimParams.mat
      #   clear newInitFile
      #elseif ( ~GUISTART_TOGGLE & exist( './guiPlotParams.mat', 'file' ) )
      #   delete guiPlotParams.mat
      #end
        
        
      ######
      ### INITIALIZE COUNTERS
    
      # numTimeSteps = 0 ;
      # timeSteps = zeros(1000000,1) ;
        
    
      ######
      ### INITIALIZE BED and ICE TOPOGRAPHY, and CLIMATE VARIABLES
        
      # Must define topo, cellsize, dx, and dy
        
      if INIT_COND_TOGGLE:
        
         ### .mat file contains: 'topo' = matrix of bed elevations and 'cellsize', 
         ### both in meters. 'easting' and 'northing' are included for plotting
            
         if INIT_COND_TOGGLE == 1:    # Valley glaciers
                
#                 filenameDEM = 'Yosemite200_rot35_400x650' ;
#                 filenameDEM = 'Nederland100' ;
#                 filenameDEM = 'KingsCanyon200Rot256x256shift' ;
#                 filenameDEM = 'sample200' ;
#                 filenameDEM = 'animas_200' ;
#                 filenameDEM = '4J_newDEM_200' ;
#                 filenameDEM = 'reproj4j_200' ;
            filenameDEM = inputFile
            filenameDEM = 'Animas_200.mat'

            #load( filenameDEM ) ;
            ( topo , easting , northing , cellsize ) = load_dem( filenameDEM )

            dx = numpy.longdouble(200)  # set a new dx
            dy = numpy.longdouble(dx)
                
            # AAR and eroded volume watershed mask
            mask_file = 'watershed_mask'
                
            try:
               #load( mask_file );
               watershed_mask = load_mask( mask_file )
            except:
               watershed_mask = numpy.ones( topo.shape , dtype=numpy.longdouble )  # Use the whole grid if no watershed mask is available
               logging.warning( 'No watershed mask found; using the whole grid for AAR and eroded flux calculations.' )
                
            # Mass Balance

            try: initELA
            except NameError:
               initELA = numpy.longdouble(3350)
               maxBz   = numpy.longdouble(2)
               gradBz  = numpy.longdouble(1./100.)
            else:
               if INIT_COND_TOGGLE == 2:    # Ice sheets
                  filenameDEM = 'Baffin200d'
                  filenameDEM = 'ValleyNonFjordTopo'
                
                  #load( filenameDEM ) ;
                  ( topo , easting , northing ) = load_dem( filenameDEM )
                            
                  dx       = numpy.longdouble(2000)  # set a new dx
                  dy       = dx 
                
                  UsChar   = numpy.longdouble(100)
                  taubChar = numpy.longdouble(50000)
                
                  #load( filenameDEM, 'Bxy' ) ;
                  Bxy = load_dem_var( filenameDEM , 'Bxy' )
                
                  # Mass Balance
                  initELA   = numpy.longdouble(3500)
                  maxBz     = numpy.longdouble(0)
                  gradBz    = numpy.longdouble(1./100)
                
                  Hbound    = numpy.longdouble(2000)
                
                  Elev0     = numpy.longdouble(0)           # reference elevation
                  To        = numpy.longdouble(-30)         # temperature at Elev0
                  lapseRate = numpy.longdouble(-0.0065)     # degrees per meter
                   
                  COMPRESS_TOGGLE         = 0
                  GENERIC_ICE_TOGGLE      = 0
                  MASS_BALANCE_TOGGLE     = ELA_TIME_SERIES
                  CALVING_TOGGLE          = 1
                  ERODE_TOGGLE            = 0
        
                  THERMAL_TOGGLE          = 0
                  FREEZEON_TOGGLE         = 0
                  HORZTL_ADVECT_TOGGLE    = 0
                  GEOTHERMAL_HEAT_TOGGLE  = 0
                  STRAIN_HEAT_TOGGLE      = 0
                  SLIDING_HEAT_TOGGLE     = 0
                  SURFACE_HEAT_FLUX_TOGGLE= 0
                  THERMAL_3D_TOGGLE       = 0
        
                  WEST_BC_TOGGLE      = ZERO_FLUX_BOUND
                  EAST_BC_TOGGLE      = ZERO_FLUX_BOUND
                  SOUTH_BC_TOGGLE     = ZERO_FLUX_BOUND
                  NORTH_BC_TOGGLE     = ZERO_FLUX_BOUND
                

               elif INIT_COND_TOGGLE == 3:    # gui_start
                  #load( filenameDEM ) ;
                  ( topo , easting , northing ) = load_dem( filenameDEM )
                  dy = dx
             
            rws,cls = topo.shape
            #if !exist('easting') : easting  = numpy.arange( cls )
            #if !exist('northing'): northing = numpy.arange( rws )

            try:              easting
            except NameError: easting  = numpy.arange( cls )
            try:              northing
            except NameError: northing = numpy.arange( rws )
         
                    
            # resample DEM at new node spacing
            if cellsize != dx:
                
               rws,cls = topo.shape
               xOld = numpy.arange(cls-1)*cellsize
               yOld = numpy.arange(rws-1)*cellsize

               #xOld = (0:cls-1)*cellsize ;
               #yOld = (0:rws-1)*cellsize ;

               XOld,YOld = numpy.meshgrid( xOld , yOld )
                
               #if rem(max(xOld),dx) == 0 and rem(max(yOld),dy) == 0:
               if max(xOld) % dx == 0 and max(yOld) % dy == 0:
                  clsNew = max(xOld)/dx + 1
                  rwsNew = max(yOld)/dy + 1
               else:
                  clsNew = numpy.ceil( xOld[-1] / dx )
                  rwsNew = numpy.ceil( yOld[-1] / dy )
                    
               x = numpy.arange(clsNew)*dx
               y = numpy.arange(rwsNew)*dy

               X,Y = numpy.meshgrid( x , y )
               
               topo     = interpolate.interp2d( XOld , YOld , topo , kind='linear' )( X , Y )
               #topo     = interpolate.interp2d( XOld , YOld , topo, X, Y ) ;

               easting  = interpolate.interp1d( xOld , easting  , kind='linear' )( x )
               northing = interpolate.interp1d( yOld , northing , kind='linear' )( y )
               cellsize = dx

            # Set the bed elevation to 'topo'
            Zb     = topo.copy()
            initZb = Zb.copy()

            #if !exist('H'): H = numpy.zeros(Zb.shape)
            try: H
            except NameError: H = numpy.zeros( Zb.shape , dtype=numpy.longdouble )

            Zi     = H + Zb
            #clear topo
            
            rws,cls = Zb.shape
            x   = numpy.arange( cls )*dx
            y   = numpy.arange( rws )*dy
            X,Y = numpy.meshgrid( x , y )

            # Create a generic ice surface
            if GENERIC_ICE_TOGGLE:
                
               # This code segment rotates the topo such that the 
               # ice boundary is on the left side of the simulation
               # need to check code; better to rotate DEM prior to use
                
               ZiBound = numpy.mean(Zb[:,0]) + Hbound
               taub    = 200000
               H       = numpy.zeros(Zb.shape, dtype=numpy.longdouble )
               rws,cls = Zb.shape
               beta    = taub/(rhoI*g)
               jtermlast = cls-2
               icefree = 0
                
               # for each row, find the cell for which the ice surface 
               # height at the left boundary would be ZiBound if the 
               # terminus that starts in that cell
               #for i =1:rws
               for i in range(rws):
            
                  mZb   = Zb[i,:]
                  slope = -numpy.diff(mZb)/dx
                    
                  # search starts in front of the terminus
                  # of the adjacent row that was just found  
                  jterm = min( jtermlast+1, cls-2 )
                  while jterm > 0:
                
                     # backwater calculation
                     mH = numpy.zeros(mZb.shape, dtype=numpy.longdouble )
                     for j in range(jterm-1,-1,-1):
            
                        term1  = ( -slope[j]/2. - (mH[j+1]/dx) )**2
                        term2  = -(2./dx) * ( slope[j] * mH[j+1] - beta )
                        deltaH = -slope[j]*dx/2. - mH[j+1] + dx * numpy.sqrt(term1+term2)
                        mH[j]  = mH[j+1] + deltaH
                        
                     # the following ensures that the search for
                     # the terminus was started beyond the terminus
                     mZi = mZb + mH
                     if mZi[0] > ZiBound:
                        icefree = 1
                     elif icefree and mZi[0] < ZiBound:
                        H[i,:] = mH
                        jtermlast = jterm
                        icefree = 0
                        break
                     else:
                        jterm = jterm + 2
                        if jterm >= cls-1:
                           logging.error( "Generic ice overruns boundary" )
                           return -1
                        
                     jterm = jterm - 1
            
               Zi = Zb + H
                           
               rws,cls = Zb.shape
               filt    = numpy.ones( (3,3) , dtype=numpy.longdouble ) / 9
               ZiBig   = numpy.zeros( (rws+2,cls+2) , dtype=numpy.longdouble )
               ZiBig[1:-1,1:-1] = Zi
            
               for i in range(10):
                  ZiBig[(0,-1),:]  = ZiBig[(1,-2),:]
                  ZiBig[:,(0,-1)]  = ZiBig[:,(1,-2)]
                  ZiBig            = filter2d( filt , ZiBig )
            
               Zi = ZiBig[1:-2,1:-2]
            
            ind = H == 0
            Zi[ind] = Zb[ind]
            
            conserveIce = H.sum(axis=0).sum()
            iceVolumeLast = conserveIce*dx*dy
            
         else:  # SYNTHETIC BEDROCK TOPOGRAPHY
            logging.error( "Must code synthetic initial condition" )
            return -1
        
      ### INIT_COND_TOGGLE
      ######

   ### Load a saved state
   ######

   # Initialize matrices
   #n_rows = 100
   #n_cols = 200
   #H  = numpy.ones( ( n_rows , n_cols ) )*100
   #Zb = numpy.ones( ( n_rows , n_cols ) )
   #Tb = numpy.ones( ( n_rows , n_cols ) )
   #Tm = numpy.ones( ( n_rows , n_cols ) )
   #Ts = numpy.ones( ( n_rows , n_cols ) )
   #
   #COMPRESS_TOGGLE = True
   #THERMAL_TOGGLE  = True
   #
   #RESTART_TOGGLE = 1

   # Start the time loop
   fd_watch = {}
   fd_watch['thick']  = open( 'thickness_py.bin' , 'wb' )
   fd_watch['taubxx'] = open( 'taubxX_py.bin' , 'wb' )
   fd_watch['taubyy'] = open( 'taubyY_py.bin' , 'wb' )
   fd_watch['taubx'] = open( 'taubX_py.bin' , 'wb' )
   fd_watch['tauby'] = open( 'taubY_py.bin' , 'wb' )
   fd_watch['xcmpnt'] = open( 'xcmpnt_py.bin' , 'wb' )
   fd_watch['ycmpnt'] = open( 'ycmpnt_py.bin' , 'wb' )
   fd_watch['udxx'] = open( 'UdxX_py.bin' , 'wb' )
   fd_watch['udyy'] = open( 'UdyY_py.bin' , 'wb' )
   fd_watch['usxx'] = open( 'UsxX_py.bin' , 'wb' )
   fd_watch['usyy'] = open( 'UsyY_py.bin' , 'wb' )

   fd_csv = open( 'dt.csv' , 'w' )

   ( H , Zb , dx , dy ) = load_state( inputFile )

   run_for( t , tMax , H , Zb , dx , dy )

   return

   counter = 0
   tic = time.time()
   while t<tMax or RESTART_TOGGLE==2:
      
      # COMPRESS - ONLY SIMULATE SUB-RECTANGLE THAT CONTAINS ICE
      ( Zi , compression_ratio , COMPRESSED_FLAG ) = compress_grid( H , Zb , COMPRESS_TOGGLE=False )

      ######
      ### MODIFY BOUNDARY CELLS TO ENFORCE BOUNDARY CONDITIONS

      ( H_ext , Zb_ext , Zi_ext ) = set_bc( H , Zb , Zi )

      ######
      ### CALCULATE THE BASAL SHEAR STRESS

      # forward differences
      #dZidxX_ext = ( Zi_ext[:,1:] - Zi_ext[:,:-1] ) / dx
      #dZidyY_ext = ( Zi_ext[1:,:] - Zi_ext[:-1,:] ) / dy
      #dZidxX     = dZidxX_ext[1:-1,:]
      #dZidyY     = dZidyY_ext[:,1:-1]

      ( dZidxX , dZidyY ) = difference_grid( Zi_ext , dx , dy )
            
      ( ( xcmpnt , ycmpnt ) , ( taubX , taubY ) , ( HX , HY ) ) = basal_shear_stress( H_ext , Zi_ext , dx=dx , dy=dy )

      ######
      ### CALCULATE ICE VELOCITY DUE TO DEFORMATION
        
      if ICEFLOW_TOGGLE:
         ( UdxX , UdyY ) = iceflow( taubX , taubY , HX , HY , xcmpnt , ycmpnt )
      else:
         UdxX    = numpy.zeros( xcmpnt.shape , dtype=numpy.longdouble )
         UdyY    = numpy.zeros( ycmpnt.shape , dtype=numpy.longdouble )

      ######
      ### CALCULATE SLIDING VELOCITY
        
      if ICESLIDE_TOGGLE:
         ( UsxX , UsyY ) = ice_sliding( taubX , taubY , xcmpnt , ycmpnt )
      else:
         UsxX = numpy.zeros( xcmpnt.shape , dtype=numpy.longdouble )
         UsyY = numpy.zeros( ycmpnt.shape , dtype=numpy.longdouble )

      # sum all contributions to ice motion
      ( UxX , UyY ) = sum_ice_motion( UdxX , UdyY , UsxX , UsyY )

      ######
      ### MASS CONSERVATION -- CONTINUITY

      ( dHdt , ( qxX , qyY ) ) = mass_conservation( H_ext , UxX , UyY , HX , HY , dZidxX , dZidyY , dx=dx , dy=dy );

      ######
      ### CALCULATE MASS BALANCE

      ( Bxy , ELA ) = mass_balance( Zi , t )
            
      #######
      ### CALCULATE TIMESTEP

      if VARIABLE_DT_TOGGLE:
         dt = get_timestep( H , Zi_ext , Zi , dHdt , Bxy , dtMax=dtMax , dtDefault=dtDefault )
      else:
         dt = dtDefault
        
      ######
      ### UPDATE the TIME and ICE THICKNESS
        
      ( t , H , Zi , conserveIce ) = update_vars( H , Zb , Zi , Bxy , qxX , qyY , dHdt , t , dt , conserveIce , dx=dx , dy=dy )

      fd_csv.write( '%f\n' % t )
      fd_csv.flush()
            
      # Calculate AAR
            
      # AccumGrid = (Zi > ELA) .* (H > 0);
      IndGlacier = numpy.choose( H  >0 , (0,watershed_mask) )
      AccumGrid  = numpy.choose( Bxy>0 , (0,IndGlacier    ) )
      AccumArea  = AccumGrid.sum(axis=0).sum()
      TotArea    = IndGlacier.sum(axis=0).sum()
      AAR        = AccumArea / TotArea

      ######
      ###  CALCULATION OF ICE TEMPERATURES
        
      if THERMAL_TOGGLE == 0:
         pass
      elif THERMAL_TOGGLE == 1:
         Ts = To + lapseRate*( Zi - Elev0 )
         Tb = Ts - gradTz * H
         Tb = numpy.choose( Tb>0 , (Tb,0) )
         Tm = Ts.copy()

         Htemp = Ts / gradTz

         ind = nonzero( H.flatten() <= Htemp )
         put( Tm , ind , ( Ts.flatten()[ind] + Tb.flatten()[ind] ) * .5 )

         ind = nonzero( H.flatten() > Htemp )
         put( Tm , ind , Ts.flatten()[ind] * (1. - Htemp.flatten()[ind] / ( 2.*H.flatten()[ind] ) ) )
            
      elif THERMAL_TOGGLE == 2:
         thermal_gc2d
        
      ######
      ### COMPRESS - ONLY SIMULATE SUB-RECTANGLE THAT CONTAINS ICE
        
      #if COMPRESS_TOGGLE and H.max() > 1 and RESTART_TOGGLE != 2:

      #   disp( 'Error!!!' )
      #   H_FullSpace  = H.copy()
      #   Zb_FullSpace = Zb.copy()
      #   if THERMAL_TOGGLE:
      #      Ts_FullSpace = Ts.copy()
      #      Tb_FullSpace = Tb.copy()
      #      Tm_FullSpace = Tm.copy()
            
      #   indrw,indcl = (H!=0).nonzero()
      #   mxrw ,mxcl  = Zb.shape
            
      #   mnrw = max( 0    , indrw.min() - 2 )
      #   mxrw = min( mxrw , indrw.max() + 2 )
      #   mncl = max( 0    , indcl.min() - 2 )
      #   mxcl = min( mxcl , indcl.max() + 2 )
            
      #   H  = H [mnrw:mxrw,mncl:mxcl]
      #   Zb = Zb[mnrw:mxrw,mncl:mxcl]
      #   Zi = Zb + numpy.choose( H<0 , (H,0) )
            
      #   rws,cls = H.shape
        
      #   if THERMAL_TOGGLE:
      #      Ts = Ts[mnrw:mxrw,mncl:mxcl]
      #      Tb = Tb[mnrw:mxrw,mncl:mxcl]
      #      Tm = Tm[mnrw:mxrw,mncl:mxcl]
            
      #   mxrws,mxcls       = Zb_FullSpace.shape
      #   rws  ,cls         = Zb.shape
      #   compression_ratio = (mxcls*mxrws)/(cls*rws)
      #   COMPRESSED_FLAG   = 1
      #else:
      #   Zi = Zb + numpy.choose( H<0 , (H,0) ) # included for restarts
      #   compression_ratio = 1
      #   COMPRESSED_FLAG   = 0
        

      # THIS IS THE END OF THE CONTINUUM CALCULATIONS 
      # NOW SIMULATE PROCESSES FOR WHICH WE HAVE NO EQUATIONS

      ######
      ### AVALANCHE SNOW OFF OF STEEP SURFACES
        
      if AVALANCH_TOGGLE and ( numpy.random.uniform() < dt*avalanchFreq ):
         avalanche( H )

      ######
      ### CALVING GLACIER FRONT
        
      if CALVING_TOGGLE:
         calve( H , dt )

      if counter%1==0:
         print_watch_point( fd_watch['thick']  , H            )
      #print_watch_point( fd_watch['taubxx'] , taubxX[:,1:] )
      #print_watch_point( fd_watch['taubyy'] , taubyY[1:,:] )
      #print_watch_point( fd_watch['taubx']  , taubX [:,1:] )
      #print_watch_point( fd_watch['tauby']  , taubY [1:,:] )
      #print_watch_point( fd_watch['xcmpnt'] , xcmpnt[:,1:] )
      #print_watch_point( fd_watch['ycmpnt'] , ycmpnt[1:,:] )
      #print_watch_point( fd_watch['udxx']   , UdxX  [:,1:] )
      #print_watch_point( fd_watch['udyy']   , UdyY  [1:,:] )
      #print_watch_point( fd_watch['usxx']   , UsxX  [:,1:] )
      #print_watch_point( fd_watch['usyy']   , UsyY  [1:,:] )

      counter += 1
      if counter > 3000: return

      ######
      ### ERODE THE BED and TRACK CRN INVENTORY

      if CRN_TOGGLE:
         CRN_gc2d             # Call the CRN module

      ######
      ### ERODE THE BED - now handled in CRN module
# 
#         if ERODE_TOGGLE:
#             erode_gc2d
#             
        

      ######
      ### REPORT SOME STUFF
      toc = time.time()
        
      if REPORT_TOGGLE and  toc >= nextReport:
         logging.info( 'elapsed time: %1.2f seconds' , (toc-tic) )
         logging.info( 'simulation time: %1.2f yr'   , t         )
         logging.info( 'timestep: %1.2e yr'          , dt        )
         logging.info( 'ELA: %1.0f m'                , ELA       )
         logging.info( 'AAR: %1.2f'                  , AAR       )
#        print 'Erosion mass flux: %1.1e kg/yr' % eroded_mass_flux
            
         # fractional ice conservation
         iceVolume = numpy.choose( H<0 , (H,0) ).sum(axis=0).sum()*dx*dy
         logging.info( 'total ice: %1.2e km^3'       , (iceVolume*1e-9)                )
         logging.info( 'excess ice: %1.2f m^3'       , (iceVolume - conserveIce*dx*dy) )
         logging.info( 'ice change: %f m^3'          , (iceVolume - iceVolumeLast)     )
         logging.info( 'max ice thickness: %1.2e km' , (H.max()/1000.)                 )
         if iceVolume != 0:
            logging.info( 'ice conservation (%%): %1.15f' , (100 - 100*( iceVolume - conserveIce*dx*dy ) / iceVolume) )

         iceVolumeLast = iceVolume
            
         if CALVING_TOGGLE:
            logging.info( 'calved ice volume: %1.2e m^3' , calvedIce )
            
         if COMPRESS_TOGGLE:
            logging.info( 'compression ratio = %f' , compression_ratio )
            
         nextReport = toc + reportInterval

   fd_watch.close()
   logging.info( "Finished!" )

   return 0

def run_for( t , t_max , H , Zb , dx , dy , ICEFLOW_TOGGLE=True , ICESLIDE_TOGGLE=False , VARIABLE_DT_TOGGLE=True , dtDefault=Parameters.dtDefault , dtMax=Parameters.dtMax):

   fd_watch = {}
   fd_watch['thick']  = open( 'thickness_py.bin' , 'wb' )

   conserveIce = numpy.longdouble(0.)
   counter = 0
   tic     = time.time()

   while t<t_max:

      # COMPRESS - ONLY SIMULATE SUB-RECTANGLE THAT CONTAINS ICE
      ( Zi , compression_ratio , COMPRESSED_FLAG ) = compress_grid( H , Zb , COMPRESS_TOGGLE=False )

      ######
      ### MODIFY BOUNDARY CELLS TO ENFORCE BOUNDARY CONDITIONS

      ( H_ext , Zb_ext , Zi_ext ) = set_bc( H , Zb , Zi )

      ( dZidxX , dZidyY ) = difference_grid( Zi_ext , dx , dy )

      ######
      ### CALCULATE THE BASAL SHEAR STRESS

      ( ( xcmpnt , ycmpnt ) , ( taubX , taubY ) , ( HX , HY ) ) = basal_shear_stress( H_ext , Zi_ext , dx=dx , dy=dy )

      ######
      ### CALCULATE ICE VELOCITY DUE TO DEFORMATION
        
      if ICEFLOW_TOGGLE:
         ( UdxX , UdyY ) = iceflow( taubX , taubY , HX , HY , xcmpnt , ycmpnt )
      else:
         UdxX    = numpy.zeros( xcmpnt.shape , dtype=numpy.longdouble )
         UdyY    = numpy.zeros( ycmpnt.shape , dtype=numpy.longdouble )

      ######
      ### CALCULATE SLIDING VELOCITY
        
      if ICESLIDE_TOGGLE:
         ( UsxX , UsyY ) = ice_sliding( taubX , taubY , xcmpnt , ycmpnt )
      else:
         UsxX = numpy.zeros( xcmpnt.shape , dtype=numpy.longdouble )
         UsyY = numpy.zeros( ycmpnt.shape , dtype=numpy.longdouble )

      # sum all contributions to ice motion
      ( UxX , UyY ) = sum_ice_motion( UdxX , UdyY , UsxX , UsyY )

      ######
      ### MASS CONSERVATION -- CONTINUITY

      ( dHdt , ( qxX , qyY ) ) = mass_conservation( H_ext , UxX , UyY , HX , HY , dZidxX , dZidyY , dx=dx , dy=dy );

      ######
      ### CALCULATE MASS BALANCE

      ( Bxy , ELA ) = mass_balance( Zi , t )
            
      #######
      ### CALCULATE TIMESTEP

      if VARIABLE_DT_TOGGLE:
         dt = get_timestep( H , Zi_ext , Zi , dHdt , Bxy , dtMax=dtMax , dtDefault=dtDefault )
      else:
         dt = dtDefault
        
      ######
      ### UPDATE the TIME and ICE THICKNESS
        
      ( t , H , Zi , conserveIce ) = update_vars( H , Zb , Zi , Bxy , qxX , qyY , dHdt , t , dt , conserveIce , dx=dx , dy=dy )

      if counter%1==0:
         print_watch_point( fd_watch['thick']  , H            )
      counter = counter + 1


class Toggles:
   # CODE BEHAVIOR TOGGLES
   # toggles turn on/off segments of the code or select 
   # between multiple possibilities for a given process
   # values can be reset in INIT_COND segment
        
   GUISTART_TOGGLE     = 0   # started simulation with the gui   (off|on)
        
   SAVE_TOGGLE         = 1   # saving                            (off|on)
   PLOT_TOGGLE         = 1   # plotting                          (off|on)
   REPORT_TOGGLE       = 1   # reporting                         (off|on)
        
   COMPRESS_TOGGLE     = 0   # only simulate area with ice       (off|on)
   VARIABLE_DT_TOGGLE  = 1   # state dependent time step         (off|on)

   INIT_COND_TOGGLE    = 1   # load DEM and climate              (synth|valley|sheet)
   GENERIC_ICE_TOGGLE  = 0   # start with generic ice surface    (off|on)       
        
   ICEFLOW_TOGGLE      = 1   # ice motion by deformation         (off|on)
   ICESLIDE_TOGGLE     = 0   # ice motion by sliding             (off|on|select)        
        
   THERMAL_TOGGLE      = 0   # temp dependance of flow           (off|on)
   FREEZEON_TOGGLE     = 0   # basal ice freeze to bed           (off|on)

   AVALANCH_TOGGLE     = 0   # avalanch off steep surfaces       (off|on)
   ERODE_TOGGLE        = 0   # erode the bed                     (off|on|select)
   CALVING_TOGGLE      = 0   # calving front                     (off|on)
        
   CRN_TOGGLE          = 0   # CRN accumulation                  (off|on)

   MASS_BALANCE_TOGGLE = MassBalance.ELA_LOWERING # select climate scenerio   (off|on|select)

   WEST_BC_TOGGLE      = BoundaryCond.ICE_FREE_BOUND  # boundary condition    (no ice|reflect|no flow)
   EAST_BC_TOGGLE      = BoundaryCond.ICE_FREE_BOUND  # boundary condition    (no ice|reflect|no flow)
   SOUTH_BC_TOGGLE     = BoundaryCond.ICE_FREE_BOUND  # boundary condition    (no ice|reflect|no flow)
   NORTH_BC_TOGGLE     = BoundaryCond.ICE_FREE_BOUND  # boundary condition    (no ice|reflect|no flow)
        

def init_valley_glacier( file='Animas_200.mat' ):

#  filenameDEM = 'Yosemite200_rot35_400x650' ;
#  filenameDEM = 'Nederland100' ;
#  filenameDEM = 'KingsCanyon200Rot256x256shift' ;
#  filenameDEM = 'sample200' ;
#  filenameDEM = 'animas_200' ;
#  filenameDEM = '4J_newDEM_200' ;
#  filenameDEM = 'reproj4j_200' ;

   ( topo , easting , northing , cellsize ) = load_dem( file )

   dx = numpy.longdouble(200)  # set a new dx
   dy = numpy.longdouble(dx)
                
   # AAR and eroded volume watershed mask
   mask_file = 'watershed_mask'
                
   try:
      #load( mask_file );
      watershed_mask = load_mask( mask_file )
   except:
      watershed_mask = numpy.ones( topo.shape , dtype=numpy.longdouble )  # Use the whole grid if no watershed mask is available
      logging.warning( 'No watershed mask found; using the whole grid for AAR and eroded flux calculations.' )
                
   # Mass Balance

   try: initELA
   except NameError:
      initELA = numpy.longdouble(3350)
      maxBz   = numpy.longdouble(2)
      gradBz  = numpy.longdouble(1./100.)

   return ( topo , easting , northing , cellsize )

def init_ice_sheet( file ):

   file = 'Baffin200d'
   file = 'ValleyNonFjordTopo'
                
   #load( filenameDEM ) ;
   ( topo , easting , northing ) = load_dem( file )
                            
   dx       = numpy.longdouble(2000)  # set a new dx
   dy       = dx 
                
   UsChar   = numpy.longdouble(100)
   taubChar = numpy.longdouble(50000)
                
   #load( filenameDEM, 'Bxy' ) ;
   Bxy = load_dem_var( filenameDEM , 'Bxy' )
                
   # Mass Balance
   initELA   = numpy.longdouble(3500)
   maxBz     = numpy.longdouble(0)
   gradBz    = numpy.longdouble(1./100)
                
   Hbound    = numpy.longdouble(2000)
                
   Elev0     = numpy.longdouble(0)           # reference elevation
   To        = numpy.longdouble(-30)         # temperature at Elev0
   lapseRate = numpy.longdouble(-0.0065)     # degrees per meter
                   
   COMPRESS_TOGGLE         = 0
   GENERIC_ICE_TOGGLE      = 0
   MASS_BALANCE_TOGGLE     = ELA_TIME_SERIES
   CALVING_TOGGLE          = 1
   ERODE_TOGGLE            = 0
        
   THERMAL_TOGGLE          = 0
   FREEZEON_TOGGLE         = 0
   HORZTL_ADVECT_TOGGLE    = 0
   GEOTHERMAL_HEAT_TOGGLE  = 0
   STRAIN_HEAT_TOGGLE      = 0
   SLIDING_HEAT_TOGGLE     = 0
   SURFACE_HEAT_FLUX_TOGGLE= 0
   THERMAL_3D_TOGGLE       = 0
        
   WEST_BC_TOGGLE      = ZERO_FLUX_BOUND
   EAST_BC_TOGGLE      = ZERO_FLUX_BOUND
   SOUTH_BC_TOGGLE     = ZERO_FLUX_BOUND
   NORTH_BC_TOGGLE     = ZERO_FLUX_BOUND

   return ( topo , easting , northing )
                
def load_state( file , RESTART_TOGGLE=0 , INIT_COND_TOGGLE=True , GENERIC_ICE_TOGGLE=False ):

   ######
   ### Load a saved state

   if RESTART_TOGGLE == 0 or RESTART_TOGGLE == 3: # LOAD A SAVED STATE

      # CODE BEHAVIOR TOGGLES
      # toggles turn on/off segments of the code or select 
      # between multiple possibilities for a given process
      # values can be reset in INIT_COND segment

      toggles = Toggles
        
      # OUTPUT BEHAVIOR
      plotInterval = 60 * 120       # seconds
      saveInterval = 100            # whole years
      reportInterval = 30           # seconds
        
      nextPlot = 0                  # initialize to plot on first timestep
      nextSave = 0                  # initialize to save on first timestep
      nextReport = 0                # initialize to report on first timestep
        
      outputFile = 'savetmp'

      ######
      ### Set numerical and physical constants
      params = Parameters

      # Constants
      g    = numpy.longdouble(9.81)                    # gravitional acceleration
      rhoI = numpy.longdouble(917)                     # density of ice
      rhoW = numpy.longdouble(1000)                    # density of water
      day  = numpy.longdouble(0.00274)                 # length of a day in years

      # Time
      t         = numpy.longdouble(0)                  # set time to zero
      tMax      = numpy.longdouble(100000)             # maximum simulation time in years
      dtMax     = numpy.longdouble(0.4 * 365*day)      # maximum timestep in years
      dtDefault = numpy.longdouble(0.4 * 365*day)      # timestep if VARIABLE_DT_TOGGLE==0

      # Glacier Properties
      MinGlacThick = numpy.longdouble(1)

      # Ice Deformation
      glensA = numpy.longdouble((6.8e-15)*3.15e7/(1e9))    # Patterson, 1994; MacGregor, 2000


      # Attractor Sliding -- only if ICESLIDE_TOGGLE==1 (generally used)
      UsChar   = numpy.longdouble(10)
      taubChar = numpy.longdouble(100000)

      # Standard Sliding -- used if ICESLIDE_TOGGLE==2 (generally not used)
      B                 = numpy.longdouble(0.0012)     # m/(Pa*yr) -- MacGregor, 2000
      DepthToWaterTable = numpy.longdouble(20)         # distance from ice surface to water table
      MaxFloatFraction  = numpy.longdouble(80)         # limits water level in ice
      Hpeff             = numpy.longdouble(20)         # effective pressure (meters of water)

      # Avalanching
      angleOfRepose = numpy.longdouble(30)
      avalanchFreq  = numpy.longdouble(3)              # average number per year

      # Calving
      seaLevel    = numpy.longdouble(-100)             # meters
      calvingCoef = numpy.longdouble(2)                # year^-1


      # Thermal
      c      = numpy.longdouble(2060)                  # specific heat capacity (J/(kg*K))
      Qg     = numpy.longdouble(0.05*3.15e7)           # Geothermal heat flux (W/m^2)*seconds/year = (J/year)/(m^2)
      gradTz = numpy.longdouble(-0.0255)               # Geothermal Gradient


      # Mass Balance
      initELA         = numpy.longdouble(4500)
      initELA         = numpy.longdouble(3000)
      gradBz          = numpy.longdouble(0.01)
      maxBz           = numpy.longdouble(2)
      ELAStepSize     = numpy.longdouble(-50)
      ELAStepInterval = numpy.longdouble(500)
      tmin            = numpy.longdouble(200)         # Years, spin-up time

      ######
      ### RELOAD INPUT ARGUMENTS
    
      #load inputArgs
      inputArgs = load_input_args

      #if ( GUISTART_TOGGLE & exist('guiSimParams.mat','file') )
      #   load guiSimParams
      #   delete guiSimParams.mat
      #   clear newInitFile
      #elseif ( ~GUISTART_TOGGLE & exist( './guiPlotParams.mat', 'file' ) )
      #   delete guiPlotParams.mat
      #end
        
        
      ######
      ### INITIALIZE COUNTERS
    
      # numTimeSteps = 0 ;
      # timeSteps = zeros(1000000,1) ;
        
    
      ######
      ### INITIALIZE BED and ICE TOPOGRAPHY, and CLIMATE VARIABLES
        
      # Must define topo, cellsize, dx, and dy
        
      if INIT_COND_TOGGLE:
        
         ### .mat file contains: 'topo' = matrix of bed elevations and 'cellsize', 
         ### both in meters. 'easting' and 'northing' are included for plotting
            
         if INIT_COND_TOGGLE == 1:    # Valley glaciers
                
#                 filenameDEM = 'Yosemite200_rot35_400x650' ;
#                 filenameDEM = 'Nederland100' ;
#                 filenameDEM = 'KingsCanyon200Rot256x256shift' ;
#                 filenameDEM = 'sample200' ;
#                 filenameDEM = 'animas_200' ;
#                 filenameDEM = '4J_newDEM_200' ;
#                 filenameDEM = 'reproj4j_200' ;
            filenameDEM = file
            filenameDEM = 'Animas_200.mat'

            #load( filenameDEM ) ;
            ( topo , easting , northing , cellsize ) = load_dem( filenameDEM )

            dx = numpy.longdouble(200)  # set a new dx
            dy = numpy.longdouble(dx)
                
            # AAR and eroded volume watershed mask
            mask_file = 'watershed_mask'
                
            try:
               #load( mask_file );
               watershed_mask = load_mask( mask_file )
            except:
               watershed_mask = numpy.ones( topo.shape , dtype=numpy.longdouble )  # Use the whole grid if no watershed mask is available
               logging.warning( 'No watershed mask found; using the whole grid for AAR and eroded flux calculations.' )
                
            # Mass Balance

            try: initELA
            except NameError:
               initELA = numpy.longdouble(3350)
               maxBz   = numpy.longdouble(2)
               gradBz  = numpy.longdouble(1./100.)

         elif INIT_COND_TOGGLE==2:    # Ice sheets

            filenameDEM = 'Baffin200d'
            filenameDEM = 'ValleyNonFjordTopo'
                
            #load( filenameDEM ) ;
            ( topo , easting , northing ) = load_dem( filenameDEM )
                            
            dx       = numpy.longdouble(2000)  # set a new dx
            dy       = dx 
                
            UsChar   = numpy.longdouble(100)
            taubChar = numpy.longdouble(50000)
                
            #load( filenameDEM, 'Bxy' ) ;
            Bxy = load_dem_var( filenameDEM , 'Bxy' )
                
            # Mass Balance
            initELA   = numpy.longdouble(3500)
            maxBz     = numpy.longdouble(0)
            gradBz    = numpy.longdouble(1./100)
                
            Hbound    = numpy.longdouble(2000)
                
            Elev0     = numpy.longdouble(0)           # reference elevation
            To        = numpy.longdouble(-30)         # temperature at Elev0
            lapseRate = numpy.longdouble(-0.0065)     # degrees per meter
                   
            COMPRESS_TOGGLE         = 0
            GENERIC_ICE_TOGGLE      = 0
            MASS_BALANCE_TOGGLE     = ELA_TIME_SERIES
            CALVING_TOGGLE          = 1
            ERODE_TOGGLE            = 0
        
            THERMAL_TOGGLE          = 0
            FREEZEON_TOGGLE         = 0
            HORZTL_ADVECT_TOGGLE    = 0
            GEOTHERMAL_HEAT_TOGGLE  = 0
            STRAIN_HEAT_TOGGLE      = 0
            SLIDING_HEAT_TOGGLE     = 0
            SURFACE_HEAT_FLUX_TOGGLE= 0
            THERMAL_3D_TOGGLE       = 0
        
            WEST_BC_TOGGLE      = ZERO_FLUX_BOUND
            EAST_BC_TOGGLE      = ZERO_FLUX_BOUND
            SOUTH_BC_TOGGLE     = ZERO_FLUX_BOUND
            NORTH_BC_TOGGLE     = ZERO_FLUX_BOUND

         elif INIT_COND_TOGGLE == 3:    # gui_start
            #load( filenameDEM ) ;
            ( topo , easting , northing ) = load_dem( filenameDEM )
            dy = dx
             
         rws,cls = topo.shape
         #if !exist('easting') : easting  = numpy.arange( cls )
         #if !exist('northing'): northing = numpy.arange( rws )

         try:              easting
         except NameError: easting  = numpy.arange( cls )
         try:              northing
         except NameError: northing = numpy.arange( rws )
         
                    
         # resample DEM at new node spacing
         if cellsize != dx:
                
            rws,cls = topo.shape
            xOld = numpy.arange(cls-1)*cellsize
            yOld = numpy.arange(rws-1)*cellsize

            #xOld = (0:cls-1)*cellsize ;
            #yOld = (0:rws-1)*cellsize ;

            XOld,YOld = numpy.meshgrid( xOld , yOld )
                
            #if rem(max(xOld),dx) == 0 and rem(max(yOld),dy) == 0:
            if max(xOld) % dx == 0 and max(yOld) % dy == 0:
               clsNew = max(xOld)/dx + 1
               rwsNew = max(yOld)/dy + 1
            else:
               clsNew = numpy.ceil( xOld[-1] / dx )
               rwsNew = numpy.ceil( yOld[-1] / dy )
                    
            x = numpy.arange(clsNew)*dx
            y = numpy.arange(rwsNew)*dy

            X,Y = numpy.meshgrid( x , y )
               
            topo     = interpolate.interp2d( XOld , YOld , topo , kind='linear' )( X , Y )
            #topo     = interpolate.interp2d( XOld , YOld , topo, X, Y ) ;

            easting  = interpolate.interp1d( xOld , easting  , kind='linear' )( x )
            northing = interpolate.interp1d( yOld , northing , kind='linear' )( y )
            cellsize = dx

         # Set the bed elevation to 'topo'
         Zb     = topo.copy()
         initZb = Zb.copy()

         #if !exist('H'): H = numpy.zeros(Zb.shape)
         try: H
         except NameError: H = numpy.zeros( Zb.shape , dtype=numpy.longdouble )

         Zi     = H + Zb
         #clear topo
            
         rws,cls = Zb.shape
         x   = numpy.arange( cls )*dx
         y   = numpy.arange( rws )*dy
         X,Y = numpy.meshgrid( x , y )

         # Create a generic ice surface
         if GENERIC_ICE_TOGGLE:
                
            # This code segment rotates the topo such that the 
            # ice boundary is on the left side of the simulation
            # need to check code; better to rotate DEM prior to use
                
            ZiBound = numpy.mean(Zb[:,0]) + Hbound
            taub    = 200000
            H       = numpy.zeros(Zb.shape, dtype=numpy.longdouble )
            rws,cls = Zb.shape
            beta    = taub/(rhoI*g)
            jtermlast = cls-2
            icefree = 0
                
            # for each row, find the cell for which the ice surface 
            # height at the left boundary would be ZiBound if the 
            # terminus that starts in that cell
            #for i =1:rws
            for i in range(rws):
            
               mZb   = Zb[i,:]
               slope = -numpy.diff(mZb)/dx
                    
               # search starts in front of the terminus
               # of the adjacent row that was just found  
               jterm = min( jtermlast+1, cls-2 )
               while jterm > 0:
                
                  # backwater calculation
                  mH = numpy.zeros(mZb.shape, dtype=numpy.longdouble )
                  for j in range(jterm-1,-1,-1):
            
                     term1  = ( -slope[j]/2. - (mH[j+1]/dx) )**2
                     term2  = -(2./dx) * ( slope[j] * mH[j+1] - beta )
                     deltaH = -slope[j]*dx/2. - mH[j+1] + dx * numpy.sqrt(term1+term2)
                     mH[j]  = mH[j+1] + deltaH
                        
                  # the following ensures that the search for
                  # the terminus was started beyond the terminus
                  mZi = mZb + mH
                  if mZi[0] > ZiBound:
                     icefree = 1
                  elif icefree and mZi[0] < ZiBound:
                     H[i,:] = mH
                     jtermlast = jterm
                     icefree = 0
                     break
                  else:
                     jterm = jterm + 2
                     if jterm >= cls-1:
                        logging.error( "Generic ice overruns boundary" )
                        return -1
                        
                  jterm = jterm - 1
            
            Zi = Zb + H
                           
            rws,cls = Zb.shape
            filt    = numpy.ones( (3,3) , dtype=numpy.longdouble ) / 9
            ZiBig   = numpy.zeros( (rws+2,cls+2) , dtype=numpy.longdouble )
            ZiBig[1:-1,1:-1] = Zi
            
            for i in range(10):
               ZiBig[(0,-1),:]  = ZiBig[(1,-2),:]
               ZiBig[:,(0,-1)]  = ZiBig[:,(1,-2)]
               ZiBig            = filter2d( filt , ZiBig )
            
            Zi = ZiBig[1:-2,1:-2]
            
         ind = H == 0
         Zi[ind] = Zb[ind]
            
         conserveIce = H.sum(axis=0).sum()
         iceVolumeLast = conserveIce*dx*dy
            
      else:  # SYNTHETIC BEDROCK TOPOGRAPHY
         logging.error( "Must code synthetic initial condition" )
         return -1
        
      ### INIT_COND_TOGGLE
      ######

   ### Load a saved state
   ######

   return ( H , Zb , dx , dy )

def gc2d( argv=None , inputFile='Animas_200.mat' ):

#   if argv is None:
#      argv = sys.argv
#   try:
#      try:
#         opts, args = getopt.getopt( argv[1:] , "h" , ["help"] )
#      except getopt.error, msg:
#         raise Usage(msg)
#   except Usage, err:
#      print >> sys.stderr, err.msg
#      print >> sys.stderr, "for help use --help"
#      return 2

   ( H , Zb , dx , dy ) = load_state( inputFile )

   run_for( Parameters.t , Parameters.tMax , H , Zb , dx , dy )

   return 0

from csdms import Component

class gc2d( Component ):
   def __init__( self ):
      self._name = 'GC2D'
      self._vars = {}
      self.set_var( 'H'  , None )
      self.set_var( 'Zb' , None )
      self.set_var( 'dx' , None )
      self.set_var( 'dy' , None )

   def setup( self , file=None ):
      Component( self._name ).setup()

      ( H , Zb , dx , dy ) = load_state( file )
      self.set_var( 'H'  , H )
      self.set_var( 'Zb' , Zb )
      self.set_var( 'dx' , dx )
      self.set_var( 'dy' , dy )

   def run_for( self , duration , start=0. ):

      #if type( duration ) == unum.Unum:
      #   duration_in_y = duration.convert( YR ).asNumber()
      #else:
      duration_in_y = duration

      #if type( start ) == unum.Unum:
      #   start_in_y = start.convert( YR ).asNumber()
      #else:
      start_in_y = start

      Component( self._name ).run()

      H  = self.get_var( 'H'  )
      Zb = self.get_var( 'Zb' )
      dx = self.get_var( 'dx' )
      dy = self.get_var( 'dy' )

      run_for( start_in_y , start_in_y+duration_in_y , H , Zb , dx , dy )

   def teardown( self ):
      Component( self._name ).teardown()


if __name__ == "__main__":
   logging.basicConfig( level=logging.INFO ,
                        format='%(asctime)s %(levelname)-8s %(message)s' ,
                        datefmt='%a, %d %b %Y %H:%M:%S' ,
                        filename='gc2d.log' ,
                        filemode='w' )
   sys.exit( gc2d() )

