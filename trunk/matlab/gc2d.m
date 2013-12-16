%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%  BEGINNING OF gc2d.m  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%  version = 2007.11.29
%%  Based on MK version 2006.04.15
%%  Adjusted by Dylan Ward to incorporate CRN accumulation/erosion layers,
%%  calculate AAR given a watershed mask.
%%
%%  function   gc2d( RESTART_TOGGLE, inputFile, varargin )
%%
%%    RESTART_TOGGLE = (0|1|2|3) = ( start | restart | reconstruct | initialize )
%%
%%    <varargin> is a list of property names each followed by a value
%%         this allows parameters to be changed from the command line
%%         when restarting a simulation ( i.e. RESTART_TOGGLE ~= 0 ).
%%
%%         example (and only for example):
%%            gc2d( 1, 'year1000', 'tMax', 2000 ) ;
%%
%%         would restart with data in year1000.mat and continue to year 2000
%%
%%    If <inputFile> is not specified, it defaults to 'savetmp.mat'.  Thus, gc2d(1)
%%         will restart from the last plotted state.
%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function    gc2d( RESTART_TOGGLE, inputFile, varargin )

if ( nargin == 0 )
	RESTART_TOGGLE = 0 ;
elseif ( nargin > 2 )
	%% Extract the input arguments from 'varargin' in the form
	%% 'name' 'value', then sets variable 'name' equal to 'value'
    
	for n=1:2:length(varargin)-1
        property = cat(1,varargin{n}) ;
        eval( sprintf( '%s = cat(1,varargin{n+1}) ;', property ) ) ;
        disp( sprintf( '%s = %s ;', property, num2str(cat(1,varargin{n+1}))) ) ;
	end

end
save inputArgs
    
        
if ( RESTART_TOGGLE == 0  |  RESTART_TOGGLE == 3 ) %%% LOAD A SAVED STATE


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% CODE BEHAVIOR TOGGLES

        %% toggles turn on/off segments of the code or select 
        %% between multiple possibilities for a given process
        %% values can be reset in INIT_COND segment
        
        GUISTART_TOGGLE     = 0 ;   % started simulation with the gui   (off|on)
        
        SAVE_TOGGLE         = 1 ;   % saving                            (off|on)
        PLOT_TOGGLE         = 1 ;   % plotting                          (off|on)
        REPORT_TOGGLE       = 1 ;   % reporting                         (off|on)
        
        COMPRESS_TOGGLE     = 0 ;   % only simulate area with ice       (off|on)
        VARIABLE_DT_TOGGLE  = 1 ;   % state dependent time step         (off|on)

        INIT_COND_TOGGLE    = 1 ;   % load DEM and climate              (synth|valley|sheet)
        GENERIC_ICE_TOGGLE  = 0 ;   % start with generic ice surface    (off|on)       
        
        ICEFLOW_TOGGLE      = 1 ;   % ice motion by deformation         (off|on)
        ICESLIDE_TOGGLE     = 1 ;   % ice motion by sliding             (off|on|select)        
        
        THERMAL_TOGGLE      = 0 ;   % temp dependance of flow           (off|on)
        FREEZEON_TOGGLE     = 0 ;   % basal ice freeze to bed           (off|on)

        AVALANCH_TOGGLE     = 1 ;   % avalanch off steep surfaces       (off|on)
        ERODE_TOGGLE        = 0 ;   % erode the bed                     (off|on|select)
        CALVING_TOGGLE      = 0 ;   % calving front                     (off|on)
        
        CRN_TOGGLE          = 0 ;   % CRN accumulation                  (off|on)
        
        
        %%% Available Mass Balance
        ZERO_BALANCE        = 1 ;   % Constant Ice Flux
        CONSTANT_ELA        = 2 ;   % Ice Free Boundary
        ELA_LOWERING        = 3 ;   % Zero Ice Flux
        ELA_TIME_SERIES     = 4 ;   % Continuous Ice Surface Slope
        EXTERNAL_FUNC       = 5 ;   % Constant Surface Elevation
        ELA_LOWERING2       = 6 ;   % Zero Ice Flux
        BALANCE_FILE        = 7 ;   % Zero Ice Flux
        D18O_TIME_SERIES    = 8 ;   % Load d18O record and convert to ELA history
        
        MASS_BALANCE_TOGGLE = ELA_LOWERING ;    % select climate scenerio   (off|on|select)
        
        
        %%% Available Boundary Conditions
        ICE_FREE_BOUND      = 1 ;   % Ice Free Boundary
        ZERO_FLUX_BOUND     = 2 ;   % Zero Ice Flux
        CONST_FLUX_BOUND    = 3 ;   % Constant Ice Flux
        SURF_ELEV_BOUND     = 4 ;   % Constant Surface Elevation
        SURF_SLOPE_BOUND    = 5 ;   % Continuous Ice Surface Slope
        
        WEST_BC_TOGGLE = ICE_FREE_BOUND ;   % boundary condition    (no ice|reflect|no flow)
        EAST_BC_TOGGLE = ICE_FREE_BOUND ;   % boundary condition    (no ice|reflect|no flow)
        SOUTH_BC_TOGGLE = ICE_FREE_BOUND ;  % boundary condition    (no ice|reflect|no flow)
        NORTH_BC_TOGGLE = ICE_FREE_BOUND ;  % boundary condition    (no ice|reflect|no flow)
        
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% OUTPUT BEHAVIOR
    
        plotInterval = 60 * 120;            % seconds
        saveInterval = 100 ;             % whole years
        reportInterval = 30 ;           % seconds
        
        nextPlot = 0 ;                  % initialize to plot on first timestep
        nextSave = 0 ;                  % initialize to save on first timestep
        nextReport = 0 ;                % initialize to report on first timestep
        
        outputFile = 'savetmp' ;
        
        
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% NUMERICAL and PHYSICAL CONSTANTS
        
        %%% Constants
        g = 9.81;                       % gravitional acceleration
        rhoI = 917;                     % density of ice
        rhoW = 1000;                    % density of water
        day = 0.00274 ;                 % length of a day in years
        

        %%% Time
        t = 0 ;                         % set time to zero
        tMax = 100000 ;                 % maximum simulation time in years
        dtMax = 0.4 * 365*day ;         % maximum timestep in years
        dtDefault = 0.4 * 365*day ;     % timestep if VARIABLE_DT_TOGGLE==0
        
        
        %%% Glacier Properties
        MinGlacThick = 1 ;
        
        
        %%% Ice Deformation
        glensA = (6.8e-15)*3.15e7/(1e9) ;    % Patterson, 1994; MacGregor, 2000
        
        
        %%% Attractor Sliding -- only if ICESLIDE_TOGGLE==1 (generally used)
        UsChar = 10 ;
        taubChar = 100000 ;
        
        
        %%% Standard Sliding -- used if ICESLIDE_TOGGLE==2 (generally not used)
        B = 0.0012 ;                    % m/(Pa*yr) -- MacGregor, 2000
        DepthToWaterTable = 20 ;        % distance from ice surface to water table
        MaxFloatFraction = 80 ;         % limits water level in ice
        Hpeff = 20 ;                    % effective pressure (meters of water)
         
        
        %%% Avalanching
        angleOfRepose = 30 ;
        avalanchFreq = 3 ;              % average number per year
        
        
        %%% Calving
        seaLevel = -100 ;               % meters
        calvingCoef = 2 ;               % year^-1
        
        
        %%% Thermal
        c = 2060 ;                      % specific heat capacity (J/(kg*K))
        Qg = 0.05*3.15e7 ;              % Geothermal heat flux (W/m^2)*seconds/year = (J/year)/(m^2)
        gradTz = -0.0255 ;              % Geothermal Gradient
        
        
        %%% Mass Balance
        initELA = 4500 ;
        gradBz = 0.01 ;
        maxBz = 2 ;
        ELAStepSize = -50 ;
        ELAStepInterval = 500 ;
        tmin = 200;                % Years, spin-up time
        
        
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% RELOAD INPUT ARGUMENTS
    
        load inputArgs
        if ( GUISTART_TOGGLE & exist('guiSimParams.mat','file') )
            load guiSimParams
            delete guiSimParams.mat
            clear newInitFile
        elseif ( ~GUISTART_TOGGLE & exist( './guiPlotParams.mat', 'file' ) )
            delete guiPlotParams.mat
        end
        
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% INITIALIZE COUNTERS
    
         % numTimeSteps = 0 ;
         % timeSteps = zeros(1000000,1) ;
        
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% INITIALIZE BED and ICE TOPOGRAPHY, and CLIMATE VARIABLES
        
        %%% must define topo, cellsize, dx, and dy
        
        if ( INIT_COND_TOGGLE )
        
            %%% .mat file contains: 'topo' = matrix of bed elevations and 'cellsize', 
            %%% both in meters. 'easting' and 'northing' are included for plotting
            
            if ( INIT_COND_TOGGLE == 1 )    %% valley glaciers
                
%                 filenameDEM = 'Yosemite200_rot35_400x650' ;
%                 filenameDEM = 'Nederland100' ;
%                 filenameDEM = 'KingsCanyon200Rot256x256shift' ;
%                 filenameDEM = 'sample200' ;
%                 filenameDEM = 'animas_200' ;
%                 filenameDEM = '4J_newDEM_200' ;
%                 filenameDEM = 'reproj4j_200' ;
                filenameDEM = inputFile ;

                load( filenameDEM ) ;
                
                dx = 200 ;  % set a new dx
                dy = dx ;
                
                % AAR and eroded volume watershed mask
                mask_file = 'watershed_mask';
                
                try
                    load( mask_file );
                catch
                    watershed_mask = ones(size(topo));  % Use the whole grid if no watershed mask is available
                    disp('No watershed mask found; using the whole grid for AAR and eroded flux calculations.')
                end
                
                %%% Mass Balance
                
                if exist('initELA') ~= 1
                initELA = 3350 ;
                maxBz = 2 ;
                gradBz = 1/100 ;
                end
                
            elseif ( INIT_COND_TOGGLE == 2 )    %% ice sheets
            
                filenameDEM = 'Baffin200d' ;
                filenameDEM = 'ValleyNonFjordTopo' ;
                
                load( filenameDEM ) ;
                            
                dx = 2000 ;  % set a new dx
                dy = dx ;
                
                UsChar = 100 ;
                taubChar = 50000 ;
                
                load( filenameDEM, 'Bxy' ) ;
                
                %%% Mass Balance
                initELA = 3500 ;
                maxBz= 0 ;
                gradBz = 1/100 ;
                
                Hbound = 2000 ;
                
                Elev0 = 0 ;             % reference elevation
                To = -30 ;              % temperature at Elev0
                lapseRate = -0.0065 ;   % degrees per meter
                   
                COMPRESS_TOGGLE         = 0 ;
                GENERIC_ICE_TOGGLE      = 0 ;
                MASS_BALANCE_TOGGLE     = ELA_TIME_SERIES ;
                CALVING_TOGGLE          = 1 ;
                ERODE_TOGGLE            = 0 ;
        
                THERMAL_TOGGLE          = 0 ;
                FREEZEON_TOGGLE         = 0 ;
                HORZTL_ADVECT_TOGGLE    = 0 ;
                GEOTHERMAL_HEAT_TOGGLE  = 0 ;
                STRAIN_HEAT_TOGGLE      = 0 ;
                SLIDING_HEAT_TOGGLE     = 0 ;
                SURFACE_HEAT_FLUX_TOGGLE= 0 ;
                THERMAL_3D_TOGGLE       = 0 ;
        
                WEST_BC_TOGGLE      = ZERO_FLUX_BOUND ;
                EAST_BC_TOGGLE      = ZERO_FLUX_BOUND ;
                SOUTH_BC_TOGGLE     = ZERO_FLUX_BOUND ;
                NORTH_BC_TOGGLE     = ZERO_FLUX_BOUND ;
                
            elseif ( INIT_COND_TOGGLE == 3 )    %% gui_start
                
                load( filenameDEM ) ;
                dy = dx ;
                
            end 
     
             
            [rws,cls] = size(topo) ;
            if ( ~exist('easting') )
                easting = 1:cls ;
            end
            if ( ~exist('northing') )
                northing = 1:rws ;
            end
            
                    
            %% resample DEM at new node spacing
            if ( cellsize ~= dx )
                
                [rws,cls] = size(topo) ;
                xOld = (0:cls-1)*cellsize ;
                yOld = (0:rws-1)*cellsize ;
                [XOld,YOld] = meshgrid(xOld,yOld);
                
                if ( rem(max(xOld),dx) == 0 & rem(max(yOld),dy) == 0 )
                    clsNew = max(xOld)/dx + 1 ;
                    rwsNew = max(yOld)/dy + 1 ;
                else
                    clsNew = ceil( xOld(end) / dx ) ;
                    rwsNew = ceil( yOld(end) / dy ) ;
                end
                    
                x = (0:clsNew-1)*dx ;
                y = (0:rwsNew-1)*dy ;
                [X, Y] = meshgrid(x,y);
            
                topo = interp2( XOld, YOld, topo, X, Y ) ;
                easting = interp1( xOld, easting, x ) ;
                northing = interp1( yOld, northing, y ) ;
                cellsize = dx ;
                
            end
            
            % Set the bed elevation to 'topo'
            Zb = topo ;
            initZb = Zb;
            if( ~exist('H') ), H = zeros(size(Zb)) ; end
            Zi = H + Zb ;
            clear topo
            
            [rws,cls] = size(Zb) ;
            x = 0:dx:(cls-1)*dx ;
            y = 0:dy:(rws-1)*dy ;
            [X,Y] = meshgrid(x,y);
            
                
            %%% Create a generic ice surface
            if ( GENERIC_ICE_TOGGLE )
                
                %% This code segment rotates the topo such that the 
                %% ice boundary is on the left side of the simulation
                %% need to check code; better to rotate DEM prior to use
                
                ZiBound = mean(Zb(:,1)) + Hbound ;
                taub = 200000 ;
                H = zeros(size(Zb)) ;
                [rws,cls] = size(Zb) ;
                beta = taub/(rhoI*g) ;
                jtermlast = cls-1 ;
                icefree = 0 ;
                
                %% for each row, find the cell for which the ice surface 
                %% height at the left boundary would be ZiBound if the 
                %% terminus that starts in that cell
                for i =1:rws
            
                    mZb = Zb(i,:) ;
                    slope = -diff(mZb)/dx ;
                    
                    % search starts in front of the terminus
                    % of the adjacent row that was just found  
                    jterm = min( jtermlast+1, cls-1 ) ;
                    while jterm > 1 
                
                        %% backwater calculation
                        mH = zeros(size(mZb)) ;
                        for j = jterm:-1:1
            
                            term1 = ( -slope(j)/2 - (mH(j+1)/dx) ).^2 ;
                            term2 = -(2/dx) * ( slope(j) * mH(j+1) - beta ) ;
                            deltaH = -slope(j)*dx/2 - mH(j+1) + dx * sqrt(term1+term2) ;
                            mH(j) = mH(j+1) + deltaH;
            
                        end
                        
                        % the following ensures that the search for
                        % the terminus was started beyond the terminus
                        mZi = mZb + mH ;
                        if ( mZi(1) > ZiBound )
                            icefree = 1 ;
                        elseif ( icefree & mZi(1) < ZiBound )
                            H(i,:) = mH ;
                            jtermlast = jterm ;
                            icefree = 0 ;
                            break ;
                        else
                            jterm = jterm + 2 ;
                            if( jterm >= cls )
                                error('generic ice overruns boundary') ;
                            end
                        end
                        
                        jterm = jterm - 1 ;
                        
                    end
                    
                end
            
                Zi = Zb + H ;
                           
                [rws,cls] = size(Zb) ;
                filt = ones(3)/9 ;
                ZiBig = zeros(rws+2,cls+2) ;
                ZiBig(2:end-1,2:end-1) = Zi ;
            
                for i=1:10
                    ZiBig([1 end],:) = ZiBig([2 end-1],:) ;
                    ZiBig(:,[1 end]) = ZiBig(:,[2 end-1]) ;
                    ZiBig = filter2(filt,ZiBig) ;
                end
            
                Zi = ZiBig(2:end-1,2:end-1) ;
       
            end
            
            ind = find( H == 0 ) ;
            Zi(ind) = Zb(ind) ;
            
            conserveIce = sum(sum(H)) ;
            iceVolumeLast = conserveIce*dx*dy ;
            
        else  %%% SYNTHETIC BEDROCK TOPOGRAPHY
        
            disp('must code synthetic initial condition') ;
            return
        
        end    %- INIT_COND_TOGGLE -%
        
         
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% CALCULATE MASS BALANCE TIME SERIES
    
        %% load in Delta O18 record and create ELA time series
        if ( MASS_BALANCE_TOGGLE == D18O_TIME_SERIES )
            
            maxELA = 3825 ;
            minELA = 3250 ;
            load delO18record
            
            dO18_max = 0;       % Filter threshold for wacko values
            
            dO18 = delO18record(:,2) ;
            
            dO18_valid = find( dO18 <= dO18_max );  % Find non-wacko values
            
            dO18 = dO18( dO18_valid );
            yrBP = delO18record(dO18_valid,3) + 2006 - 1949 ;
           
            n = 5 ;
            dO18f = filter2( ones(n,1)/n, dO18, 'valid' ) ;
            ELArecord = flipud(minELA + (maxELA-minELA) * ...
                    ((dO18f-min(dO18f))/max(dO18f-min(dO18f)))) ;
            trecord = flipud(yrBP(end) - yrBP) ;
            trecord = trecord(3:end-2) ;
            trecord = trecord - min(trecord) ;

            % spin-up model to steady state
            tmin = 1000 ;
            trecord = [0 trecord' + tmin]' ;
            ELArecord = [ELArecord(1) ELArecord']' ;                
        end

        %% Load in pre-made ELA time series (C1: Years BP  C2: ELA (m))
        if ( MASS_BALANCE_TOGGLE == ELA_TIME_SERIES )
            
            load NewELArecord
            ELAseries = NewELArecord(:,2) ;
            tSeries = NewELArecord(:,1) ;
            
            ELArecord = flipud(ELAseries) ;

            % spin-up model to steady state
            tmin = 1000 ;
            trecord = flipud((tMax-tmin) - tSeries) ;
            trecord = [0 trecord' + tmin]' ;
            ELArecord = [ELArecord(1) ELArecord']' ;             
            
        end
     
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% SET INITIAL THERMAL STATE
    
        if ( THERMAL_TOGGLE )
        
            % thermal_gc2d.m creates this file
            try
                load GlensFlowCoef
            catch
                getGlensA_gc2d
                load GlensFlowCoef
            end
            
            Ts = To + lapseRate*( Zi - Elev0 ) ;
            Tb = min( 0, Ts - gradTz * H ) ;
            Tm = Ts ;

            Htemp = Ts/gradTz ;

            ind = find( H <= Htemp ) ;
            Tm(ind) = ( Ts(ind) + Tb(ind) ) / 2 ;

            ind = find( H > Htemp ) ;
            Tm(ind) = Ts(ind) .* ( 1 - Htemp(ind)./(2*H(ind)) ) ;

            if ( THERMAL_3D_TOGGLE )
                thermal_layers = 10 ;
                [rws,cls] = size(Zb) ;
                T3D = zeros(rws,cls,thermal_layers) ;
                T3D = getTemp3D_gc2d( H, Ts, Tm, T3D ) ;
            end
        
        end
        

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% INITIALIZE RANDOM NUMBER GENERATOR
    
        rand('state',1);


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% SAVE INITIAL SIMULATION STATE
    
        rndstate = rand('state') ;
        initFilename = 'yearInit0' ;
        save( initFilename ) ;
        if ( RESTART_TOGGLE == 3 )
            return
        end
        
        
elseif ( RESTART_TOGGLE == 1  |  RESTART_TOGGLE == 2 ) %%% LOAD A SAVED STATE
    
    disp(sprintf('\nRECONSTRUCTING WORKSPACE')) ;
    
    %% Set the default input file to 'savetmp' if 'inputFile' is not specified
    
        if ( ~exist('inputFile') | isempty(inputFile) )
            inputFile = 'savetmp' ;
            save inputArgs inputFile
        end 
        disp( ['inputFile = '  inputFile] ) ;
        
    
    %% Initialize simulation and load desired state
        
        try
            try
                load( inputFile, 'initFilename' )     % load input initFilename
                load( initFilename )         % initialize simulation
            catch
                gc2d(3) ;                   % create 'yearInit0'
                load yearInit0              % initialize simulation
            end
            
            load inputArgs inputFile        % load input filename
            load( inputFile )               % load specified file
        catch
            disp( lasterr )
            return
        end
        
        
    %% Changes to toggles, parameters or variables
        
        nextPlot = 0 ;      % initialize to plot on next timestep
        nextReport = 0 ;    % initialize to report on next timestep
        
        load inputArgs      % reload the input args in case they were
                            % overwritten when inputFile was loaded
    
        if ( GUISTART_TOGGLE & exist('guiSimParams.mat','file') )
            load guiSimParams
            delete guiSimParams.mat
            clear newInitFile
        end
        
    
    %% Set random number generator to previous position
    
        rand( 'state', rndstate );
        
       
    %% Save simulation state   
    
        initFilename = sprintf( 'yearInit%d', round(t) ) ;
        save( initFilename ) ;
        
        
else

    error('Unrecognized value for RESTART_TOGGLE') ;
        
end %- RESTART_TOGGLE -%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%              START THE TIME LOOP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%try     % this try/catch statement saves the state if the simulation dumps
    
    diary_name = sprintf('%s_log', date ) ;
    diary( diary_name ) ;   
    tic
    while ( t < tMax | RESTART_TOGGLE==2 )
            
            
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% INTERACT WITH GUI
        
        if ( RESTART_TOGGLE ~= 2 & GUISTART_TOGGLE )
            
            if( exist( './guiSimParams.mat', 'file' ) )
                
                load guiSimParams
                
                disp('Loaded changed parameters:')
                params = who( '-file', 'guiSimParams.mat' ) ;
                for p=1:length(params)
                    eval( params{p} )
                end
                
                delete guiSimParams.mat
                
                if ( exist( 'stopSimulation' ) )
                    rndstate = rand('state') ;
                    clear stopSimulation
                    save savetmp
                    return
                end
                
                if ( exist('newInitFile') )
                    initFilename = sprintf( 'yearInit%d', round(t) ) ;
                    save( initFilename ) ;
                    clear newInitFile
                end
                
            else
                drawnow
            end
            
        end
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% COMPRESS - ONLY SIMULATE SUB-RECTANGLE THAT CONTAINS ICE
        
        if ( COMPRESS_TOGGLE & max(max(H)) > 1 & RESTART_TOGGLE ~= 2 )
            
            H_FullSpace = H ;
            Zb_FullSpace = Zb ;
            if ( THERMAL_TOGGLE )
                Ts_FullSpace = Ts ;
                Tb_FullSpace = Tb ;
                Tm_FullSpace = Tm ;
            end
            
            [indrw,indcl] = find(H ~= 0) ;
            [mxrw,mxcl] = size(Zb) ;
            
            mnrw = max( 1, min(indrw) - 2 ) ;
            mxrw = min( mxrw, max(indrw) + 2 ) ;
            mncl = max( 1, min(indcl) - 2 ) ;
            mxcl = min( mxcl, max(indcl) + 2 ) ;
            
            H = H(mnrw:mxrw,mncl:mxcl) ;
            Zb = Zb(mnrw:mxrw,mncl:mxcl) ;
            Zi = Zb + max( H, 0 ) ;
            
            [rws,cls] = size(H) ;
        
            if ( THERMAL_TOGGLE )
                Ts = Ts(mnrw:mxrw,mncl:mxcl) ;
                Tb = Tb(mnrw:mxrw,mncl:mxcl) ;
                Tm = Tm(mnrw:mxrw,mncl:mxcl) ;
            end
            
            [mxrws,mxcls] = size(Zb_FullSpace) ;
            [rws,cls] = size(Zb) ;
            compression_ratio = (mxcls*mxrws)/(cls*rws) ;
            COMPRESSED_FLAG = 1 ;
            
        else
        
            Zi = Zb + max( H, 0 ) ; %% included for restarts
            compression_ratio = 1 ;
            COMPRESSED_FLAG = 0 ;
            
        end
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% MODIFY BOUNDARY CELLS TO ENFORCE BOUNDARY CONDITIONS
        
            %%% DEFAULT BOUNDARY CONDITION IS ZERO FLUX
            H_ext = [H(1,1) H(1,:) H(1,end)
                    H(:,1) H H(:,end)
                    H(end,1) H(end,:) H(end,end)] ;
                 
            Zb_ext = [Zb(1,1) Zb(1,:) Zb(1,end)
                    Zb(:,1) Zb Zb(:,end)
                    Zb(end,1) Zb(end,:) Zb(end,end)] ;
        
            Zi_ext = [Zi(1,1) Zi(1,:) Zi(1,end)
                    Zi(:,1) Zi Zi(:,end)
                    Zi(end,1) Zi(end,:) Zi(end,end)] ;
        
            if ( THERMAL_TOGGLE )
        
                Ts_ext = [Ts(1,1) Ts(1,:) Ts(1,end)
                        Ts(:,1) Ts Ts(:,end)
                        Ts(end,1) Ts(end,:) Ts(end,end)] ;
        
                Tb_ext = [Tb(1,1) Tb(1,:) Tb(1,end)
                        Tb(:,1) Tb Tb(:,end)
                        Tb(end,1) Tb(end,:) Tb(end,end)] ;
                 
                Tm_ext = [Tm(1,1) Tm(1,:) Tm(1,end)
                        Tm(:,1) Tm Tm(:,end)
                        Tm(end,1) Tm(end,:) Tm(end,end)] ;
            end
        
            %%% WESTERN BOUNDARY CONDTION
            if ( WEST_BC_TOGGLE == SURF_ELEV_BOUND )            % Constant Ice Surface Height
                ZiBound = mean(Zb(:,1)) + Hbound ;
                H_ext(:,1) = ZiBound - Zb_ext(:,1)  ;
            elseif ( WEST_BC_TOGGLE == CONST_FLUX_BOUND )         % Constant Ice Flux B.C.
            elseif ( WEST_BC_TOGGLE == SURF_SLOPE_BOUND )       % Constant Ice Surface Slope
                Zi_ext(:,1) = 2*Zi_ext(:,2) - Zi_ext(:,3) ;
                H_ext(:,1) = Zi_ext(:,1) - Zb_ext(:,1) ;
                H_ext(:,1) = max( 0, H_ext(:,1) ) ;
            elseif( WEST_BC_TOGGLE == ICE_FREE_BOUND )          % Ice Free Boundary
                H_ext(:,1) = 0 ;
            end

            %%% EASTERN BOUNDARY CONDTION
            if ( EAST_BC_TOGGLE == SURF_ELEV_BOUND )            % Constant Ice Surface Height
                ZiBound = mean(Zb(:,end)) + Hbound ;
                H_ext(:,end) = ZiBound - Zb_ext(:,end)  ;
            elseif ( EAST_BC_TOGGLE == CONST_FLUX_BOUND )         % Constant Ice Flux B.C.
            elseif ( EAST_BC_TOGGLE == SURF_SLOPE_BOUND )       % Constant Ice Surface Slope
                Zi_ext(:,end) = 2*Zi_ext(:,end-1) - Zi_ext(:,end-2) ;
                H_ext(:,end) = Zi_ext(:,end) - Zb_ext(:,end) ;
                H_ext(:,end) = max( 0, H_ext(:,end) ) ;
            elseif( EAST_BC_TOGGLE == ICE_FREE_BOUND )          % Ice Free Boundary
                H_ext(:,end) = 0 ;
            end
            
            %%% SOUTHERN BOUNDARY CONDTION
            if ( SOUTH_BC_TOGGLE == SURF_ELEV_BOUND )           % Constant Ice Surface Height
                ZiBound = mean(Zb(1,:)) + Hbound ;
                H_ext(1,:) = ZiBound - Zb_ext(1,:)  ;
            elseif ( SOUTH_BC_TOGGLE == CONST_FLUX_BOUND )        % Constant Ice Flux B.C.
            elseif ( SOUTH_BC_TOGGLE == SURF_SLOPE_BOUND )      % Constant Ice Surface Slope
                Zi_ext(1,:) = 2*Zi_ext(2,:) - Zi_ext(3,:) ;
                H_ext(1,:) = Zi_ext(1,:) - Zb_ext(1,:) ;
                H_ext(1,:) = max( 0, H_ext(1,:) ) ;
            elseif( SOUTH_BC_TOGGLE == ICE_FREE_BOUND )          % Ice Free Boundary
                H_ext(1,:) = 0 ;
            end
            
            %%% NORTHERN BOUNDARY CONDTION
            if ( NORTH_BC_TOGGLE == SURF_ELEV_BOUND )           % Constant Ice Surface Height
                ZiBound = mean(Zb(end,:)) + Hbound ;
                H_ext(end,:) = ZiBound - Zb_ext(end,:)  ;
            elseif ( NORTH_BC_TOGGLE == CONST_FLUX_BOUND )        % Constant Ice Flux B.C.
            elseif ( NORTH_BC_TOGGLE == SURF_SLOPE_BOUND )      % Constant Ice Surface Slope
                Zi_ext(end,:) = 2*Zi_ext(end-1,:) - Zi_ext(end-2,:) ;
                H_ext(end,:) = Zi_ext(end,:) - Zb_ext(end,:) ;
                H_ext(end,:) = max( 0, H_ext(end,:) ) ;
            elseif( NORTH_BC_TOGGLE == ICE_FREE_BOUND )         % Ice Free Boundary
                H_ext(end,:) = 0 ;
            end
        
            Zi_ext = Zb_ext + H_ext ;
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% CALCULATE THE BASAL SHEAR STRESS
            
            % forward differences
            dZidxX_ext = ( Zi_ext(:,2:end) - Zi_ext(:,1:end-1) ) / dx ;
            dZidyY_ext = ( Zi_ext(2:end,:) - Zi_ext(1:end-1,:) ) / dy ;
            dZidxX = dZidxX_ext(2:end-1,:) ;
            dZidyY = dZidyY_ext(:,2:end-1) ;

            HX_ext = ( H_ext(:,2:end) + H_ext(:,1:end-1) ) / 2 ;
            HY_ext = ( H_ext(2:end,:) + H_ext(1:end-1,:) ) / 2 ;
            HX = HX_ext(2:end-1,:) ;
            HY = HY_ext(:,2:end-1) ;
            
            taubxX_ext = -rhoI * g * HX_ext .* dZidxX_ext ;
            taubyY_ext = -rhoI * g * HY_ext .* dZidyY_ext ;
            
            taubxX = taubxX_ext(2:end-1,:) ;
            taubyY = taubyY_ext(:,2:end-1) ;
            
            taubxY = ( taubxX_ext(1:end-1,1:end-1) + taubxX_ext(1:end-1,2:end) + ...
                       taubxX_ext(2:end,1:end-1) + taubxX_ext(2:end,2:end) ) / 4;
            
            taubyX = ( taubyY_ext(1:end-1,1:end-1) + taubyY_ext(1:end-1,2:end) + ...
                       taubyY_ext(2:end,1:end-1) + taubyY_ext(2:end,2:end) ) / 4;
            
            taubX = sqrt( taubxX.^2 + taubyX.^2 ) ;
            taubY = sqrt( taubxY.^2 + taubyY.^2 ) ;
       
            taubX = taubX .* ( HX > 0 ) ;
            taubY = taubY .* ( HY > 0 ) ;
            %taubX = numpy.choose( HX<0 , (HX,0) )
            %taubY = numpy.choose( HY<0 , (HY,0) )
            
            % fill in zero values with 1 for use in division
            taubOnesX    = taubX
            taubNot_indX = find( taubX == 0 ) ;
            taubOnesX(taubNot_indX) = 1 ;
            taubOnesY = taubY ;
            taubNot_indY = find( taubY == 0 ) ;
            taubOnesY(taubNot_indY) = 1 ;
            
            xcmpnt = taubxX./taubOnesX ;
            xcmpnt(taubNot_indX) = 0 ;
            
            ycmpnt = (taubyY./taubOnesY) ;
            ycmpnt(taubNot_indY) = 0 ;
            
            
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% CALCULATE ICE VELOCITY DUE TO DEFORMATION
        
        if ( ICEFLOW_TOGGLE )
            
            if ( THERMAL_TOGGLE )

                A_ext = zeros(size(H_ext)) ;
                ind = find( H_ext >= MinGlacThick ) ;
                Ts_ext = To + lapseRate*( Zi_ext - Elev0 ) ;

                A_ext(ind) = interp3( eHs, eTs, eTm, eA, H_ext(ind), Ts_ext(ind), Tm_ext(ind) ) ;
                                
                if( length( find( isnan(A_ext) ) ) > 0 )
                    save workspacedump
                    error('NaN in A, likely H_node exceeds H_glens limits') ;
                end
            
                AX = ( A_ext(2:end-1,1:end-1) + A_ext(2:end-1,2:end) )/2 ;
                AY = ( A_ext(1:end-1,2:end-1) + A_ext(2:end,2:end-1) )/2 ;
                
            else
            
                AX = glensA ;
                AY = glensA ;
                
            end

            %% here's the guts of calculating the depth averaged velocity
            UdxX = abs( (2/5) * AX .* (taubX.*taubX.*taubX) .* HX ) .* xcmpnt ;
            UdyY = abs( (2/5) * AY .* (taubY.*taubY.*taubY) .* HY ) .* ycmpnt ;
                    
        else  % need variables filled with zero values
        
            [rws,cls] = size(Zb) ;
            UdxX = zeros(rws,cls+1) ;
            UdyY = zeros(rws+1,cls) ;
            
        end %- ICEFLOW_TOGGLE -%
            
            
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% CALCULATE SLIDING VELOCITY
        
        if ( ICESLIDE_TOGGLE == 1  )        %%% ATTRACTOR FORMULATION
       
            if ( THERMAL_TOGGLE & FREEZEON_TOGGLE )
            
                notFrozen = ( Tb_ext > -0.5 | Zb_ext < seaLevel );
                notFrozenX = ( notFrozen(2:end-1,1:end-1) + notFrozen(2:end-1,2:end) )/2 ;
                notFrozenY = ( notFrozen(1:end-1,2:end-1) + notFrozen(2:end,2:end-1) )/2 ;
               
                %% here's the guts of calculating the sliding velocity 
                UsxX = notFrozenX .* UsChar .* exp(1 - taubChar./taubOnesX) .* xcmpnt ;
                UsyY = notFrozenY .* UsChar .* exp(1 - taubChar./taubOnesY) .* ycmpnt ;
                
            else
            
                %% here's the guts of calculating the sliding velocity
                UsxX = UsChar .* exp(1 - taubChar./taubOnesX) .* xcmpnt ; 
                UsyY = UsChar .* exp(1 - taubChar./taubOnesY) .* ycmpnt ; 
                                                  
            end             
                        
        else  % need variables filled with zero values
        
            [rws,cls] = size(Zb) ;
            UsxX = zeros(rws,cls+1) ;
            UsyY = zeros(rws+1,cls) ;

        end %- ICESLIDE_TOGGLE -%
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% MASS CONSERVATION -- CONTINUITY

            % ensure that no ice is drawn from the rock
            CLASS = ( H_ext >= MinGlacThick ) ;
            
            DCLASSx = ( CLASS(2:end-1,2:end) - CLASS(2:end-1,1:end-1) ) .* ...
                        sign( dZidxX ) ;
                     
            DCLASSy = ( CLASS(2:end,2:end-1) - CLASS(1:end-1,2:end-1) ) .* ...
                        sign( dZidyY ) ;
            
            % sum all contributions to ice motion
            UxX = UdxX + UsxX ;
            UyY = UdyY + UsyY ;
                        
            ind = find( DCLASSx == -1 ) ;
            UxX(ind) = 0;
                        
            ind = find( DCLASSy == -1 ) ;
            UyY(ind) = 0;
            
            % calculate both components of the ice flux
            qxX = UxX .* HX ;
            qyY = UyY .* HY ;
            
            if ( WEST_BC_TOGGLE == CONST_FLUX_BOUND )
				qxX(:,1) = BoundaryFlux ;
            end
            
            if ( EAST_BC_TOGGLE == CONST_FLUX_BOUND )
				qxX(:,end) = BoundaryFlux ;
            end
            
            if ( SOUTH_BC_TOGGLE == CONST_FLUX_BOUND )
				qyY(1,:) = BoundaryFlux ;
            end
            
            if ( NORTH_BC_TOGGLE == CONST_FLUX_BOUND )
				qyY(end,:) = BoundaryFlux ;
            end
            
            % here's the guts of the continuity equation
            dqdxX = ( qxX(:,2:end) - qxX(:,1:end-1) ) / dx ;
            dqdyY = ( qyY(2:end,:) - qyY(1:end-1,:) ) / dy ;
            dHdt = -dqdxX -dqdyY ;
            
            
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% CALCULATE MASS BALANCE
            
            %% the imposed mass balance is the imposed climate
            %% there are many possibilities, here are only a few
            %% all must populate the 2D matrix Bxy of size = size(Zb)
            %% with values of net precip/melt rate in m/yr
            %% define the scalar, ELA (m), for plotting
            
            if ( MASS_BALANCE_TOGGLE == CONSTANT_ELA )
            % Simple ELA, maxBz, gradBz

                ELA = initELA ;
                Bxy = min( maxBz, gradBz * ( Zi - ELA ) ) ;
            
            elseif ( MASS_BALANCE_TOGGLE == ELA_LOWERING )
            % ELA changing with time experiment
            
                % ELAStepSize = -10 ;       % positive/negative values raise/lower ELA
                % ELAStepInterval = 500 ;
                
                ELA = initELA + ELAStepSize * max( 0, floor( (t-tmin)/ELAStepInterval ) ) ;
                Bxy = min( maxBz, gradBz * ( Zi - ELA ) ) ;
            
            elseif ( MASS_BALANCE_TOGGLE == ELA_LOWERING2 )
            % ELA changing with time experiment
            
                tau = 25 ;              % intrinsic timescale of ice dynamics 
                tmin = 0 ;              % time to begin ELA modification
                initELA = 4200 ;        % initial ELA
                stepSize = -10 ;        % positive/negative values raise/lower ELA
                dELAdt = -0.1 ;
                
                ELA = initELA + stepSize * max( 0, floor( (t-tmin) / (8*tau) ) ) ;
                Bxy = min( maxBz, gradBz * ( Zi - ELA ) ) ;
            
            elseif ( MASS_BALANCE_TOGGLE == EXTERNAL_FUNC )
            % external mass balance function
                
                if ( ~exist('Bxy') | t >= nextGetBxy )
            
                    % Mass Balance 2D Must Return Bxy (2d Matrix)
                    Bxy = mass_balance_gc2d( t, cellsize, Zi ) ;
                    nextGetBxy = t + getBxyInterval ;
                    
                end
            
            elseif ( MASS_BALANCE_TOGGLE == ELA_TIME_SERIES | MASS_BALANCE_TOGGLE == D18O_TIME_SERIES)
            % ELA time series

                ELA = interp1( trecord, ELArecord, t ) ;
                Bxy = min( maxBz, gradBz * ( Zi - ELA ) ) ;
            
            elseif ( MASS_BALANCE_TOGGLE == BALANCE_FILE )
            % external mass balance file
                
                load( filenameDEM, 'Bxy' )
                ind = find( abs(Bxy) == min(min(abs(Bxy))) ) ;
                ELA = mean( Zi(ind) ) ;
        
            elseif ( MASS_BALANCE_TOGGLE == ZERO_BALANCE )
            
                ELA = 0 ;
                Bxy = zeros( size(Zb) ) ;
                
            else
            
                error('Unrecognized Mass Balance')
                
            end    %- MASS_BALANCE_TOGGLE -%
            
            
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% CALCULATE TIMESTEP
        
        if ( VARIABLE_DT_TOGGLE )
            
            %% now that we know the rate of change in ice surface heights due to  
            %% ice motion and due to precipitation or melt we need to know over 
            %% what period of time we can project forward with these rates and 
            %% maintain stability of the ice surface.  The basic idea here is that
            %% we don't want to take a timestep any longer then it would take to 
            %% reverse the ice surface slope between two cells, such that ice 
            %% should be flowing in the other direction.  In fact, let's make our 
            %% timestep much less then that.
            
            %% this calculation sets the timestep such that the change
            %% in ice surface elevation nowhere exceeds a set fraction
            %% of the local standard deviation in ice surface elevations
            
            % include ice changes by precip and melt
            dHdtTot = dHdt + Bxy ;
            adHdt = abs(dHdtTot) ;
            
            % something like standard deviation of 3x3 cell areas around each cell
            filt = [1 1 1;1 1 1;1 1 1] / 9 ;
            ZiMean = filter2( filt, Zi_ext, 'valid' ) ;
            dHmax = sqrt( filter2( filt, (ZiMean - Zi).^2 ) ) ;
            
            % only consider cells with ice thickness > 10 m
            isGlac = H > 10 ;
            
            % find limiting timestep for each considered cell
            ind = find( adHdt~=0 & dHmax~=0 & isGlac~=0 ) ;    
            dtLimits = dHmax(ind)./adHdt(ind) ;
            [dt, idt] = min( dtLimits ) ;
            
            % locate the x and y position of limiting cell for plotting
            [rwDT,clDT] = ind2sub( size(adHdt), ind(idt) ) ; 
            
            % limit timestep to dtMax or some fraction of the calculated timestep
            dt = min( dtMax, dt/2 ) ;
            
            % catch an error, (e.g. if H<10 in all cells )
            if isempty(dt)  
                dt = dtDefault ;
            end
            
        else

            dt = dtDefault ;
                        
        end %- VARIABLE_DT_TOGGLE -%
        
        
            
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% UPDATE the TIME and ICE THICKNESS
        
            %% time update
            t = t + dt ;
            % numTimeSteps = numTimeSteps + 1 ;
            % timeSteps(numTimeSteps) = dt ;
            
            %% increase in ice thicknesses due to precip
            Bxy_pos = Bxy .* ( Bxy > 0 ) ;
            H = H + Bxy_pos * dt ;
            
            %% change ice thicknesses due to ice motion
            H = H + dHdt*dt ;
            
            %% decrease in ice thicknesses due to melt
            Bxy_neg = Bxy .* ( Bxy < 0 ) ;
            Bxy_neg = -min( H, -Bxy_neg ) ;
            H = H + Bxy_neg * dt ;
            
            %% record ice addition or removal by climate
            snowFall = ( Bxy_neg + Bxy_pos ) * dt ;
            conserveIce = conserveIce + sum(sum(snowFall));
            
            %% record ice flux through boundaries
            qbound = sum(qyY(1,:)) - sum(qyY(end,:)) + sum(qxX(:,1)) - sum(qxX(:,end)) ;
            conserveIce = conserveIce + dt * qbound / dx ;
            
            Zi = Zb + max( H, 0 );
        
            if ( isnan(Zi) )
                save workspacedump
                error('NaN in ice thickness') ;
            end
            
            % Calculate AAR
            
            % AccumGrid = (Zi > ELA) .* (H > 0);
            IndGlacier = (H > 0) .* watershed_mask;
            AccumGrid = (Bxy > 0) .* IndGlacier;
            AccumArea = sum(sum(AccumGrid));
            TotArea = sum(sum(IndGlacier));
            AAR = AccumArea / TotArea;
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%  CALCULATION OF ICE TEMPERATURES
        
        if( THERMAL_TOGGLE == 0 )
        
        elseif( THERMAL_TOGGLE == 1 )

            Ts = To + lapseRate*( Zi - Elev0 ) ;
            Tb = min( 0, Ts - gradTz * H ) ;
            Tm = Ts ;

            Htemp = Ts/gradTz ;

            ind = find( H <= Htemp ) ;
            Tm(ind) = ( Ts(ind) + Tb(ind) ) / 2 ;

            ind = find( H > Htemp ) ;
            Tm(ind) = Ts(ind) .* ( 1 - Htemp(ind)./(2*H(ind)) ) ;
            
        elseif( THERMAL_TOGGLE == 2 )
        
            thermal_gc2d

        end
        
    
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% UNCOMPRESS - FILL SPACE WITH SUB-RECTANGLE THAT CONTAINS ICE
        
        if ( COMPRESSED_FLAG )
            
            H_FullSpace(mnrw:mxrw,mncl:mxcl) = H ;
            Zb_FullSpace(mnrw:mxrw,mncl:mxcl) = Zb ;
        
            if ( THERMAL_TOGGLE )
                Ts_FullSpace(mnrw:mxrw,mncl:mxcl) = Ts ;
                Tb_FullSpace(mnrw:mxrw,mncl:mxcl) = Tb ;
                Tm_FullSpace(mnrw:mxrw,mncl:mxcl) = Tm ;
            end
            
            H = H_FullSpace ;
            Zb = Zb_FullSpace ;
            Zi = Zb + max( H, 0 ) ;
            
            if ( THERMAL_TOGGLE )
                Ts = Ts_FullSpace ;
                Tb = Tb_FullSpace ;
                Tm = Tm_FullSpace ;
            end
        
            COMPRESSED_FLAG = 0 ;
        
        end



    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% THIS IS THE END OF THE CONTINUUM CALCULATIONS 
    %%% NOW SIMULATE PROCESSES FOR WHICH WE HAVE NO EQUATIONS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% AVALANCH SNOW OFF OF STEEP SURFACES
        
        if ( AVALANCH_TOGGLE & ( rand < dt*avalanchFreq ) )
            
            %% move ice downslope until the ice surface is everywhere
            %% less then or near the angle of repose
            
            [rws,cls] = size(Zb) ;
            dHRepose = dx*tan(angleOfRepose*pi/180) ;
            Ho = max( H, 0 ) ;
            
            while ( 1 )

                dZidx_down = zeros(rws,cls);
                dZidx_down(:,2:end) = max( 0, ( Zi(:,2:end) - Zi(:,1:end-1) ) ) ;
                dZidx_up = zeros(rws,cls);
                dZidx_up(:,1:end-1) = max( 0, ( Zi(:,1:end-1) - Zi(:,2:end) ) ) ;
                dZidx = max( dZidx_up, dZidx_down ) ;

                dZidy_left = zeros(rws,cls);
                dZidy_left(2:end,:) = max( 0, ( Zi(2:end,:) - Zi(1:end-1,:) ) ) ;
                dZidy_right = zeros(rws,cls);
                dZidy_right(1:end-1,:) = max( 0, ( Zi(1:end-1,:) - Zi(2:end,:) ) ) ;
                dZidy = max( dZidy_left, dZidy_right ) ;

                grad = sqrt( dZidx.^2 + dZidy.^2 );
                gradT =  dZidy_left + dZidy_right + dZidx_down + dZidx_up ;
                gradT(find(gradT==0)) = 1;
                grad( find(Ho < 0.1) ) = 0 ;

                mxGrad = max(max( grad ) ) ;
        
                if ( mxGrad <= 1.1*dHRepose )
                    break ;
                end

                delH = max( 0, ( grad - dHRepose ) / 3 ) ;
        
                Htmp = Ho ;        
                Ho = max( 0, Htmp - delH );
                delH = Htmp - Ho ;
        
                delHdn = zeros(rws,cls) ; delHup = zeros(rws,cls) ;
                delHlt = zeros(rws,cls) ; delHrt = zeros(rws,cls) ;
        
                delHup(:,2:end) = delH(:,1:end-1) .* dZidx_up(:,1:end-1)./gradT(:,1:end-1) ;
                delHdn(:,1:end-1) = delH(:,2:end) .* dZidx_down(:,2:end)./gradT(:,2:end) ;
                delHrt(2:end,:) = delH(1:end-1,:) .* dZidy_right(1:end-1,:)./gradT(1:end-1,:) ;
                delHlt(1:end-1,:) = delH(2:end,:) .* dZidy_left(2:end,:)./gradT(2:end,:) ;
        
                Ho = max( 0, Ho + delHdn + delHup + delHlt + delHrt ) ;
        
                Zi = Zb + Ho ;
        
            end
            
            H = Ho + (H<0).*H ;
            
        end %- AVALANCH_TOGGLE -%
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% CALVING GLACIER FRONT
        
        if ( CALVING_TOGGLE )
            
            %% one reason this is difficult is that the height of ice in the cell
            %% is really just recording the volume of ice, the position of the 
            %% margin in the cell not the actual ice height.  Here floation
            %% height is assumed (or higher if necessary to account for ice volume)
            
            Hold = H ;
            calvedIce = 0 ;
            
            % count time backwards with a sshorted timestep until the whole 
            % timestep used during this itteration has been simulated
            
            dtTot = dt ;
            while ( dtTot > 0 )

                % find the calving front, aka the wet glacier margin
                G = H > 1 ;                
                W = ( G==0 & Zb <= seaLevel ) ;
                filt = [0 1 0; 1 1 1; 0 1 0] ;
                Wfilt = filter2( filt, W ) ;
                Wfilt(:,[1 end]) = Wfilt(:,[3 end-2]) ;
                Wfilt([1 end],:) = Wfilt([3 end-2],:) ;
                wetGmargin = G.*Wfilt > 0 ;
                indWGM = find( wetGmargin ) ;
                
                % if calving front exists, find water depth, ensure it's positive
                if ( ~isempty(indWGM) )
                    WDmarg = max( 0, seaLevel - Zb(indWGM) ) ;
                    ind = find( WDmarg == 0 ) ;
                    indWGM(ind) = [] ;
                    WDmarg(ind) = [] ;
                end
                
                % if calving front exists, remove some ice
                if ( ~isempty(indWGM) )
                
                    % ice thickness in calving cells
                    Hmarg = max( H(indWGM), WDmarg/0.917 ) ;
                    
                    % a new timestep is calculated such that the calving rate times the 
                    % timesstep does not exceed the total contents of any calving cell.
                    
                    dLinCalvdt = calvingCoef .* WDmarg ;                        % front migration rate
                    dVolCalvdt = dx * dLinCalvdt .* Hmarg ;                     % rate of volume calved
                    volAvailMarg = dx * dx * H(indWGM) ;                        % ice volume available
                    calvDt = min( dtTot, min( volAvailMarg ./ dVolCalvdt ) ) ;  % calving timestep

                    % remove this calving timestep from total time to calve
                    dtTot = dtTot - calvDt ;
                    
                    % convert the volume calved to ice thickness and remove
                    calve = dVolCalvdt * calvDt / ( dx * dx ) ;
                    H(indWGM) = H(indWGM) - calve ;
                
                    % record total volume calved for posterity
                    calvedIce = calvedIce + sum(calve) * dx * dx;
                    
                else
                    dtTot = 0 ;
                end
                
            end
                
            %% record ice removal by calving for conservation test
            conserveIce = conserveIce + sum(sum( H - Hold ) );
            
        end %- CALVING_TOGGLE -%
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% ERODE THE BED and TRACK CRN INVENTORY

        if ( CRN_TOGGLE )
                    
            CRN_gc2d             % Call the CRN module
            
        end %- ERODE_TOGGLE -%
    
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% ERODE THE BED - now handled in CRN module
% 
%         if ( ERODE_TOGGLE )
%                     
%             erode_gc2d
%             
%         end %- ERODE_TOGGLE -%
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% REPORT SOME STUFF
        
        if ( REPORT_TOGGLE &  toc >= nextReport )
                
            disp( sprintf('\nelapsed time: %1.2f s', toc ) ) ;
            disp( sprintf('simulation time: %1.2f yr', t) ) ;
            disp( sprintf('timestep: %1.2e yr', dt) ) ;
            disp( sprintf('ELA: %1.0f m', ELA) ) ;
            disp( sprintf('AAR: %1.2f', AAR) ) ;
%            disp( sprintf('Erosion mass flux: %1.1e kg/yr', eroded_mass_flux) ) ;
            
            %% fractional ice conservation
            iceVolume = sum(sum(H.*(H>0)))*dx*dy ;
            disp( sprintf('total ice: %1.2e km^3', iceVolume/1e9 ) ) ;
            disp( sprintf('excess ice: %1.2f m^3', (( iceVolume - conserveIce*dx*dy )) ) ) ;
            disp( sprintf('ice change: %f m^3', iceVolume - iceVolumeLast ) ) ;
            disp( sprintf('max ice thickness: %1.2e km', max(max(H))/1000 ) ) ;
            if ( iceVolume ~= 0 )
                disp( sprintf('ice conservation (%%): %1.15f', 100 - 100*( iceVolume - conserveIce*dx*dy ) / iceVolume ) ) ;
            end
            iceVolumeLast = iceVolume ;
            
            if( CALVING_TOGGLE )
                disp( sprintf('calved ice volume: %1.2e m^3', calvedIce ) ) ;
            end
            
            if( COMPRESS_TOGGLE )
                disp( sprintf('compression ratio = %f', compression_ratio ) ) ;
            end
            
            nextReport = toc + reportInterval ;
            
        end %- REPORT -%
        
    
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% SAVE RECONSTRUCTED WORKSPACE and EXIT
    
        if  ( RESTART_TOGGLE == 2 )
        
            save( outputFile )
            return
            
        end
        
            
            
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% PLOT SOME STUFF
        
        if ( PLOT_TOGGLE & toc >= nextPlot & RESTART_TOGGLE ~= 2 )
            
            rndstate = rand('state') ;
            save( outputFile )
            plot_gc2d( outputFile )
            
            nextPlot = toc + plotInterval ;
            
        end %- PLOT -%
        
            
            
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% SAVE SOME STUFF    
                       
        if ( SAVE_TOGGLE & t >= nextSave )
        
            nextSave = nextSave + saveInterval ;
            rndstate = rand('state') ;
            elapsedtime = toc ;
            savefile = sprintf( 'year%d', round(t) ) ;
            
            if ( THERMAL_TOGGLE & CRN_TOGGLE )
                eval( sprintf('save %s initFilename ELA AAR BeGrid P0Be P0Al Tm H Zb cumE eroded_mass_flux tExposed t dt conserveIce rndstate nextSave', savefile ) ) ;
            elseif ( THERMAL_TOGGLE )
                eval( sprintf('save %s initFilename ELA AAR Tm H Zb t dt conserveIce rndstate nextSave', savefile ) ) ;
            elseif ( CRN_TOGGLE )
                eval( sprintf('save %s initFilename ELA AAR BeGrid P0Be P0Al H Zb cumE eroded_mass_flux tExposed t dt conserveIce rndstate nextSave', savefile ) ) ;
            else
                eval( sprintf('save %s initFilename ELA AAR H Zb t dt conserveIce rndstate nextSave', savefile ) ) ;
            end
            
        end %- SAVE -%
            
    end %- TIME LOOP -%
    diary off
    
% catch
%   
%     err = lasterror ;
%     errdepth = length(err.stack) ;
%   
%     disp( sprintf('\n%s\n', err.message) ) ;
%     disp( sprintf('error at line %d in %s', err.stack(1).line, err.stack(1).file ) ) ;
%     for i=2:errdepth
%         disp( sprintf('called at line %d in %s', err.stack(i).line, err.stack(i).file ) ) ;
%     end
%     
%     rndstate = rand('state') ;
%     save workspacedump
%     disp(sprintf('\nworkspace saved in workspacedump.mat\n')) ;
%     diary off
%     return
%     
% end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% ONE LAST REPORT AND PLOT BEFORE EXIT

if ( REPORT_TOGGLE )
    disp( sprintf('\nelapsed time: %1.2f s', toc ) ) ;
    disp( sprintf('simulation time: %1.2f yr', t) ) ;
    disp( sprintf('timestep: %1.2e yr', dt) ) ;
    disp( sprintf('ELA: %1.0f m', ELA) ) ;
    disp( sprintf('AAR: %1.2f', AAR) ) ;
%    disp( sprintf('Erosion mass flux: %1.1e kg/yr', eroded_mass_flux) ) ;    

    %% fractional ice conservation
    iceVolume = sum(sum(H.*(H>0)))*dx*dy ;
    disp( sprintf('total ice: %1.2e km^3', iceVolume/1e9 ) ) ;
    disp( sprintf('excess ice: %1.2f m^3', (( iceVolume - conserveIce*dx*dy )) ) ) ;
    disp( sprintf('ice conservation (%%): %1.15f', 100 - 100*( iceVolume - conserveIce*dx*dy ) / iceVolume ) ) ;
    disp( sprintf('max ice thickness: %1.2e km', max(max(H))/1000 ) ) ;        
end

if ( PLOT_TOGGLE )
    rndstate = rand('state') ;
    save( outputFile )
    plot_gc2d( outputFile )
end

diary off
 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       END OF gc2d.m     %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

