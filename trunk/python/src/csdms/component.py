#!/usr/bin/env python

#import unum
#from unum.units import *

#MIN = unum.Unum.unit( 'minute' , 60.0 * unum.units.S )
#HR  = unum.Unum.unit( 'hour'   , 60.0 * MIN          )
#DAY = unum.Unum.unit( 'day'    , 24.0 * HR           )
#YR  = unum.Unum.unit( 'year'   , 365.0 * DAY         )

command_pre = "\x1B[00;34m"
command_post = "\x1B[00m"

class TimeUnit:
   ( Seconds ,
     Days ,
     Years ) = range( 3 )

class Component:
   def __init__( self , name ):
      self._name = name
      self._vars = {}

   def name( self ):
      return self._name

   def setup( self ):
      print "   %s%s/%s%s" % (command_pre,self._name,'Setup',command_post)

   def run( self ):
      print "   %s%s/%s%s" % (command_pre,self._name,'Run',command_post)

   def run_for( self , duration , start=0. , units=TimeUnit.Seconds ):
      print "   %s%s/%s%s" % (command_pre,self._name,'Run',command_post)

   def teardown( self ):
      print "   %s%s/%s%s" % (command_pre,self._name,'Teardown',command_post)

   def set_var( self , key , value ):
      self._vars[key] = value

   def get_var( self , key , units=None ):
      if self._vars.has_key( key ): return self._vars[key]
      else:                         return None

