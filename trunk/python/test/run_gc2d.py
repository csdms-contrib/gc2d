#!/usr/bin/env python

import models

glacier_model = models.gc2d()

glacier_model.setup( 'Animas_200.mat' )
#gc2d.run_for( 1. * YR )
#gc2d.run_for( 25. * DAY )
glacier_model.run_for( 1. )
glacier_model.run_for( .25 )

t = glacier_model.get_var( 'Thickness' )

glacier_model.teardown()

