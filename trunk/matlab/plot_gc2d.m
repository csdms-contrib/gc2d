%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%  BEGINNING OF plot_gc2d.m  %%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%  version = 2006.04.15
%%
%%  function  plot_gc2d( inputFile, varargin )
%%
%%      <varargin> is a list of property names each followed by a value
%%
%%          example: 
%%              plot_gc2d( 'year1000', 'azView', 30, 'elView', 60 ) ;
%%
%%          would plot the data in year1000.mat with a 3D perspective view
%%          from an azimuth of 30 degees and an elevation of 60 degrees.
%%
%%          Variables that can be set and their default values:
%%
%%          DEBUG_TOGGLE  = 0 ;
%%          azView = 0 ;
%%          elView = 90 ;
%%          minGlacThick = 10 ;
%%          zScale = 1 ;
%%          plotContent = 'Extent' ;
%%          figureNumber = 0 ;
%%          NEW_FIGURE = 0 ;
%%          CLEAR_FIGURE = 1 ;
%%          QUIVER_VECS = 0 ;
%%          ELA_CONTOUR = 0 ;
%%          ICE_CONTOUR = 1 ;
%%          THERMAL_CONTOUR = 0 ;
%%          DT_LIMIT = 0 ;
%%          SUBFIGURE = 0 ;
%%          CONTOUR_INTERVAL = 100 ;
%%          RECONSTRUCT = 0 ;
%%          
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function    plot_gc2d( inputFile, varargin )
        
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% SET DEFAULT VALUES

        DEBUG_TOGGLE  = 0 ;
        [azView, elView] = view ;
        minGlacThick = 10 ;
        zScale = 1 ;
        plotContent = 'Extent' ;
        figureNumber = 0 ;
        NEW_FIGURE = 0 ;
        CLEAR_FIGURE = 1 ;
        QUIVER_VECS = 0 ;
        ELA_CONTOUR = 1 ;
        ICE_CONTOUR = 1 ;
        THERMAL_CONTOUR = 0 ;
        DT_LIMIT = 0 ;
        SUBFIGURE = 0 ;
        CONTOUR_INTERVAL = 50 ;
        RECONSTRUCT = 0 ;

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%  Extract the input arguments from 'varargin' in the form
    %%%  'variable' 'value', then set 'variable' equal to 'value'
  
        if ( exist('guiPlotParams.mat', 'file' ) )
            load guiPlotParams
        end
        
        if ( ~exist('inputFile') ), inputFile  = 'savetmp' ; end
        disp( ['inputFile = '  inputFile] ) ;
        

%         for n=1:2:length(varargin)-1
%             property = num2str(cat(1,varargin{n})) ;
%                 eval( sprintf( '%s = cat(1,varargin{n+1}) ;', property ) ) ;
%                 disp( sprintf( '%s = %s ;', property, num2str(cat(1,varargin{n+1}))) ) ;
%         end
        
        save inputArgs
        
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% LOAD FILE TO PLOT
    
        if ( RECONSTRUCT )
            gc2d( 2, inputFile, 'outputFile', 'recontmp', 'dtMax', 0 ) ;
            load recontmp
        else
            load( inputFile )
        end
 
        load inputArgs
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% SETUP FIGURE FOR PLOTTING

        if ( figureNumber == 0 )
            figureNumber = get(0,'CurrentFigure') ;
            if( isempty(figureNumber) )
                figureNumber = 1 ;
                NEW_FIGURE = 1 ;
            end
        elseif ( NEW_FIGURE == 0 )
            try
                h = get( figureNumber ) ;
            catch
                NEW_FIGURE = 1 ;
            end
        end     
   
        if ( NEW_FIGURE )
    
            h = figure( figureNumber ) ; clf
            s = get(0,'ScreenSize') ;
            figsize = 400 ;
            [rws,cls] = size(H) ;
            if( cls > rws )
                scale = cls/rws ;
                set( figureNumber, 'Position', [10 s(4)-(figsize+85) figsize*scale figsize] )
            else
                scale = 0.9*(rws/cls) ;
                set( figureNumber, 'Position', [10 s(4)-(figsize*scale+85) figsize figsize*scale] )
            end
    
        elseif ( CLEAR_FIGURE )
    
            figure(figureNumber) ;
            [azView,elView] = view ;
            if ( SUBFIGURE )
                AxisPosition = get( gca, 'Position') ;
                cla( figureNumber );
            else
                clf( figureNumber );
            end
        end
        
    
% try      
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% INITIALIZE SOME PLOTTING VARIABLES
        
        [rws,cls] = size(Zb) ;
        x = 0:dx:(cls-1)*dx ;
        y = 0:dy:(rws-1)*dy ;
        [X,Y] = meshgrid(x,y);
    
        x_ext = -dx:dx:cls*dx ;
        y_ext = -dy:dy:rws*dy ;
        [X_ext,Y_ext] = meshgrid(x_ext,y_ext);
    
        Xin = X(2:end-1,2:end-1) ;
        Yin = Y(2:end-1,2:end-1) ;
        
    
        if ( zScale ~= 0 )
            Zb = zScale * Zb ;
            Zi = zScale * Zi ;
            H = zScale * H ;
            minGlacThick = zScale * minGlacThick ;
        end
          
            
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% PLOTTING
    
        if ( ~exist('landsat') )
            cmap = colormap('gray') ;
            minZb = min(min(Zb)) ;
            maxZb = max(max(Zb)) ;
            textureZb = 1 + 50*(Zb - minZb)/(maxZb - minZb) ;
            landsat = uint8( round( 255*interp1(1:64, cmap, textureZb )) );
        end
        
        [rws,cls,lrs] = size(landsat) ;
        xintrp = linspace( x(1), x(end), cls ) ;
        yintrp = linspace( y(1), y(end), rws ) ;
        [Xintrp,Yintrp] = meshgrid(xintrp,yintrp) ;
        Hintrp = interp2( X, Y, H, Xintrp, Yintrp ) ;
        Ziintrp = interp2( X, Y, Zi, Xintrp, Yintrp ) ;
        Zbintrp = interp2( X, Y, Zb, Xintrp, Yintrp ) ;
        
        landsatsum = double(landsat(:,:,1)) + double(landsat(:,:,2)) ...
                     + double(landsat(:,:,3)) ;
        landsatsum3 = ones(size(landsat)) ;
        landsatsum3(:,:,1) = landsatsum ;
        landsatsum3(:,:,2) = landsatsum ;
        landsatsum3(:,:,3) = landsatsum ;
        
        ind = find( landsatsum3 > 700 ) ;
        landsat(ind) = 210 ;
            
        if ( strcmp( plotContent, 'Extent' ) )
        
            if (0)
        
                ind = find(Hintrp>minGlacThick) ;
                landsat([ind rws*cls+ind 2*rws*cls+ind]) = 255 ;
            
            elseif (1)
        
                ind = find(Hintrp>minGlacThick) ;
                bw = max( 0, min( 255, round(landsatsum(ind)/3) ) ) ;
                landsat([ind]) = bw ;
                landsat([rws*cls+ind]) = bw ;
                landsat([2*rws*cls+ind]) = bw ;
            
            elseif (1)
        
                ind = find(Hintrp>minGlacThick & Ziintrp>ELA) ;
                landsat([ind rws*cls+ind 2*rws*cls+ind]) = 255 ;
            
                ind = find(Hintrp>minGlacThick & Ziintrp<=ELA) ;
                landsat([ind]) = 245 ;
                landsat([rws*cls+ind]) = 245 ;
                landsat([2*rws*cls+ind]) = 255 ;
            
            end
            
            titleText = 'Ice Extent' ; 
            
            h = surf(X/1000,Y/1000,Zi/1000,landsat, ... 
                'EdgeColor','none','FaceColor','texturemap');
            set(gca,'Visible','off')
            h = title(titleText,'fontname','times','fontsize',16) ;
            set( h, 'visible', 'on' ) ;
            colorbar('SouthOutside','off' )
        
        else
            
            ind = find(Hintrp>minGlacThick) ;
            
            if ( strcmp( plotContent, 'Discharge' ) )
    
                Uy = ( UyY(1:end-1,:) + UyY(2:end,:) )/ 2 ;
                Ux = ( UxX(:,1:end-1) + UxX(:,2:end) )/ 2 ;
                U = sqrt(Ux.^2+Uy.^2) ;
                Q = U.*H ;
                P = Q ;
            
                titleText = 'Ice Discharge' ;
                colorbarText = 'discharge (m^2/yr)' ;
        
            elseif ( strcmp( plotContent, 'Velocity' ) )
    
                Uy = ( UyY(1:end-1,:) + UyY(2:end,:) )/ 2 ;
                Ux = ( UxX(:,1:end-1) + UxX(:,2:end) )/ 2 ;
                U = sqrt(Ux.^2+Uy.^2) ;
                P = U ;
            
                titleText = 'Total Ice Velocity' ;
                colorbarText = 'velocity (m/yr)' ;
            
            elseif ( strcmp( plotContent, 'Slide' ) )
    
                Usy = ( UsyY(1:end-1,:) + UsyY(2:end,:) )/ 2 ;
                Usx = ( UsxX(:,1:end-1) + UsxX(:,2:end) )/ 2 ;
                Us = sqrt(Usx.^2+Usy.^2) ;
                P = Us ;
            
                titleText = 'Ice Sliding Velocity' ;
                colorbarText = 'velocity (m/yr)' ;
            
            elseif ( strcmp( plotContent, 'Deform' ) )
    
                Udy = ( UdyY(1:end-1,:) + UdyY(2:end,:) )/ 2 ;
                Udx = ( UdxX(:,1:end-1) + UdxX(:,2:end) )/ 2 ;
                Ud = sqrt(Udx.^2+Udy.^2) ;
                P = Ud ;
            
                titleText = 'Depth Averaged Ice Deformation Velocity' ;
                colorbarText = 'velocity (m/yr)' ;
            
            elseif ( strcmp( plotContent, 'TauB' ) )
            
                tauby = ( taubyY(1:end-1,:) + taubyY(2:end,:) )/ 2 ;
                taubx = ( taubxX(:,1:end-1) + taubxX(:,2:end) )/ 2 ;
                taub = sqrt(taubx.^2+tauby.^2) ;
                P = taub/1000 ;
            
                titleText = 'Basal Shear Stress' ;
                colorbarText = 'shear stress (kPa)' ;
            
            elseif ( strcmp( plotContent, 'Balance' ) )
    
                P = Bxy ;
                titleText = 'Net Mass Balance' ;
                colorbarText = 'balance (m/yr)' ;
            
            elseif ( strcmp( plotContent, 'CRN_Be' ) )
    
                P = BeGrid ;
                titleText = '^1^0Be concentration' ;
                colorbarText = '(at/g-qtz)' ;
                ind = find( Hintrp < inf ) ;

            elseif ( strcmp( plotContent, 'CRN_Age_Be' ) )

                if ~exist('P0Be')         % Recalculate production grid if necessary

                    % Load the calculated topographic shielding grid if it exists

                    TopoShieldFile = 'FourthOfJulyShielding.mat';       % File extension required! (i.e., '.mat')

                    if exist(TopoShieldFile)
                        load(TopoShieldFile);
                    else
                        corrGrid = topo_corr;   % else use the default shielding value
                        disp('Topographic shielding grid not found. Default shielding factor used.')
                    end

                        P0Be = production(P_Be_SLHL,rho,zmax,lat,Zb,corrGrid);    
                        P0Al = 6.1 .* P0Be;  
                        % P0Cl = 0;
                        % P0He = 0;

                end
                
                P = BeGrid ./ P0Be;
                titleText = 'Apparent ^1^0Be age' ;
                colorbarText = '(years)' ;
                ind = find( Hintrp < inf ) ;
                
            elseif ( strcmp( plotContent, 'Erosion' ) )
    
                P = cumE ;
                titleText = 'Erosion (m)' ;
                colorbarText = 'Meters' ;
                ind = find( Hintrp < inf ) ;
                
            elseif ( strcmp( plotContent, 'Exposure_Time' ) )
    
                P = tExposed ;
                titleText = 'Exposure time (yr)' ;
                colorbarText = 'Years' ;
                ind = find( Hintrp < inf ) ;
                
                            
            elseif ( strcmp( plotContent, 'Exposure_Time_frac' ) )
    
                P = tExposed ./ tMax ;
                titleText = 'Fractional exposure time' ;
                colorbarText = 't/t_m_a_x' ;
                ind = find( Hintrp < inf ) ;
                
            end
            
            Pmax = max(max(P(1:end-10,11:end))) ;
            Pmin = min(min(P(1:end-10,11:end))) ;
            Pintrp = interp2( X, Y, P, Xintrp, Yintrp ) ;
            Pintrp = min( Pmax, Pintrp );
            
            clr = max( 1, min( 64, ceil(64*Pintrp(ind)/Pmax) ) ) ;
            cmap = colormap('jet') ;
            landsat([ind; rws*cls+ind; 2*rws*cls+ind]) = ceil(255*cmap(clr,:)) ;    
            
            h = surf(X/1000,Y/1000,Zi/1000,landsat, ... 
                'EdgeColor','none','FaceColor','texturemap');
            set(gca,'Visible','off')
           
            hc = colorbar( 'SouthOutside' ) ;
            set(hc, 'XTick', linspace(0,1,5) ) ;
            xticks = linspace(Pmin,Pmax,5) ;
            eval( sprintf('xtlbls = [%1.1f %1.1f %1.1f %1.1f %1.1f] ;',...
                xticks(1), xticks(2), xticks(3), xticks(4), xticks(5) ) ) ;
            set(hc, 'XTickLabel', xtlbls ) ;
            
            hct = get( hc,'title' ) ;
            set( hct, 'String', colorbarText, 'FontSize', 12 ) ;
            set( hct, 'Position', [0.5 -4 1.00005] ) ;

        end
        
        hold on
        
        
        
        if ( ICE_CONTOUR )
            
            V = min(min(Zb)):CONTOUR_INTERVAL*zScale:max(max(Zb)) ;
            c = contourc( X(1,:)/1000, Y(:,1)/1000, Zi/1000, V/1000 ) ;
            mask = H > minGlacThick/2 ;
            i = 1 ;
            while i < length(c)
            
                cz = c(1,i) ;
                num = c(2,i) ;
                cnext = i+num+1;
                
                glac = interp2( X/1000, Y/1000, mask, c(1,i+1:i+num), c(2,i+1:i+num) ) ;
                if( sum(glac) >=2 )
                    for j = 1:num-1
                        i = i+1 ;
                        if ( glac(j)+glac(j+1) == 2 )
	                        plot3( c(1,i:i+1), c(2,i:i+1), [cz cz], 'k', 'linewidth', 1 ) ;
                        end
                    end
                end
                
                i = cnext ;
            end
            
        end
        
        if ( ELA_CONTOUR )
        
            if ( ~exist('ELA') )
                if ( exist('initELA') )
                    ELA = initELA ;
                elseif ( exist('Bxy') )
                    ind = find( abs(Bxy) == min(min(abs(Bxy))) ) ;
                    ELA = mean( Zi(ind) ) ;
                else
                    ELA = 0 ;
                end 
            end
        
            if ( ELA ~= 0 )
                V = [ELA ELA]/1000  ;
                [c,h] = contour3(Xin/1000, Yin/1000, Zi(2:end-1,2:end-1)/1000, V ) ;
                set(h, 'LineWidth', 2 );
            end
            
        end
        
        if ( 0 & exist('seaLevel') )
        
            V = 1 + [seaLevel seaLevel]/1000  ;
            [c,h] = contour3(X(2:end-1,2:end-1)/1000, Y(2:end-1,2:end-1)/1000, 1 + Zi(2:end-1,2:end-1)/1000, V ) ;
            set(h, 'LineWidth', 2, 'EdgeColor', 'k' );
        
        end
        
        if ( THERMAL_TOGGLE )
          
            if( 0 )
            
                % changing this only illustrates what the thermal
                % cold/warm boundary would approximately look like
                To = -35 ;                                 
                lapseRate = -0.0065 ;                      
                gradTz = -0.0255 ;                         
                Ts = To + lapseRate*( Zi - Elev0 ) ;       
                Tb = min( 0, Ts - gradTz * H ) ;           
                Tb_ext = [Tb(1,1) Tb(1,:) Tb(1,end)        
                        Tb(:,1) Tb Tb(:,end)               
                        Tb(end,1) Tb(end,:) Tb(end,end)] ;
                
                notFrozen = ( Tb_ext > -0.01 | Zb_ext < seaLevel ) ;

            end
            
            nF = notFrozen.*(H_ext>minGlacThick) ;
            S = 9999*(2*nF-1) ;
            Ss = filter2(ones(3)/9, S, 'valid') ;
            
            V = [5 5]  ;
            [c,h] = contour3(X/1000, Y/1000, Ss, V ) ;
            set(h, 'LineWidth', 2 );
           
        end
        
        if ( ~exist('landsat') )
            colormap( cmap );
            shading interp;
        end
        
        cpatch = [0,0,0] ;
        
        xpatch = [round(X(1,1)/1000) X(1,:)/1000 round(X(1,end)/1000)] ;
        ypatch = [round(Y(1,1)/1000) Y(1,:)/1000 round(Y(1,end)/1000)] ;
        zmin = min(min(Zb))/1000 ;
        zpatch = [ zmin Zb(1,:)/1000 zmin ] ;
        fill3( xpatch, ypatch, zpatch, cpatch ) ;
        
        xpatch = [round(X(end,1)/1000) X(end,:)/1000 round(X(end,end)/1000)] ;
        ypatch = [round(Y(end,1)/1000) Y(end,:)/1000 round(Y(end,end)/1000)] ;
        zmin = min(min(Zb))/1000 ;
        zpatch = [ zmin Zb(end,:)/1000 zmin ] ;
        fill3( xpatch, ypatch, zpatch, cpatch ) ;
    
        xpatch = [round(X(1,1)/1000); X(:,1)/1000; round(X(end,1)/1000)] ;
        ypatch = [round(Y(1,1)/1000); Y(:,1)/1000; round(Y(end,1)/1000)] ;
        zmin = min(min(Zb))/1000 ;
        zpatch = [ zmin; Zb(:,1)/1000; zmin ] ;
        fill3( xpatch, ypatch, zpatch, cpatch ) ;
        
        xpatch = [round(X(1,end)/1000); X(:,end)/1000; round(X(end,end)/1000)] ;
        ypatch = [round(Y(1,end)/1000); Y(:,end)/1000; round(Y(end,end)/1000)] ;
        zmin = min(min(Zb))/1000 ;
        zpatch = [ zmin; Zb(:,end)/1000; zmin ] ;
        fill3( xpatch, ypatch, zpatch, cpatch ) ;
        
        minE = min(easting) ;
        minN = min(northing) ;
    
        xterm = round( (344754-minE)/dx ) ;
        yterm = round( (4074019-minN)/dx ) ;
        
        
        if ( 0 & CALVING_TOGGLE )
        
            G = H > 10 ;                
            W = ( G==0 & Zb <= seaLevel ) ;
            filt = [0 1 0; 1 1 1; 0 1 0] ;
            Wfilt = filter2( filt,     W ) ;
            Wfilt(:,[1 end]) = Wfilt(:,[3 end-2]) ;
            Wfilt([1 end],:) = Wfilt([3 end-2],:) ;
            
            wetGmargin = G.*Wfilt > 0 ;
            indWGM = find( wetGmargin ) ;
            
            if ( ~isempty(indWGM) )
                plot3( X(indWGM)/1000, Y(indWGM)/1000, Zi(indWGM)/1000, 'ko' ) ;
            end
            
        end
            
        if ( QUIVER_VECS )
    
            Usx = (UsxX(:,1:end-1)+UsxX(:,2:end))/2 ;
            Usy = (UsyY(1:end-1,:)+UsyY(2:end,:))/2 ;
            Udx = (UdxX(:,1:end-1)+UdxX(:,2:end))/2 ;
            Udy = (UdyY(1:end-1,:)+UdyY(2:end,:))/2 ;
            
            Ux = Usx + Udx ;
            Uy = Usy + Udy ;
            U = sqrt( Ux.*Ux + Uy.*Uy ) ;
            
            n = 3 ;
            Xq = X(1:n:end,1:n:end)/1000 ;
            Yq = Y(1:n:end,1:n:end)/1000 ;
            Zq = Zi(1:n:end,1:n:end)/1000 ;
            Uxq = Ux(1:n:end,1:n:end) ;
            Uyq = Uy(1:n:end,1:n:end) ;
            Usxq = Usx(1:n:end,1:n:end) ;
            Usyq = Usy(1:n:end,1:n:end) ;
            Udxq = Udx(1:n:end,1:n:end) ;
            Udyq = Udy(1:n:end,1:n:end) ;
            Uq = U(1:n:end,1:n:end) ;
            Hq = H(1:n:end,1:n:end) ;
            
            ind = find( Hq < minGlacThick ) ;
            Xq(ind) = [] ; Yq(ind) = [] ; Zq(ind) = [] ;
            Uxq(ind) = [] ; Uyq(ind) = [] ;
            Usxq(ind) = [] ; Usyq(ind) = [] ;
            Udxq(ind) = [] ; Udyq(ind) = [] ;
            Uq(ind) = [] ;
                
            % for n=1:10
            %     ind = find( Uq == max(max(Uq)) ) ;
            %     Uq(ind) = 0 ;
            %     Uxq(ind) = 0 ;
            %     Uyq(ind) = 0 ;
            % end
            
            if ( QUIVER_VECS == 1 )
                quiver3( Xq, Yq, Zq, Uxq, Uyq, 0, 5, 'k', 'linewidth', 1 ) ;
            elseif ( QUIVER_VECS == 2 )
                quiver3( Xq, Yq, Zq, Usxq, Usyq, 0, 5, 'k', 'linewidth', 1 ) ;
            elseif ( QUIVER_VECS == 3 )
                quiver3( Xq, Yq, Zq, Udxq, Udyq, 0, 5, 'k', 'linewidth', 1 ) ;
            end
            
        end

            
        if ( DT_LIMIT & exist('rwDT') )
    
            if ( ~isempty(rwDT) )
    
                mxZi = max(max(Zi)) / 1000 ;
    
                plot3(X(rwDT,clDT)/1000, Y(rwDT,clDT)/1000, mxZi, 'm*', 'MarkerSize', 8, 'linewidth', 2 ) ;
                plot3(X(rwDT,clDT)/1000, Y(rwDT,clDT)/1000, mxZi, 'mo', 'MarkerSize', 10, 'linewidth', 2 ) ;
                plot3(X(rwDT,clDT)/1000, Y(rwDT,clDT)/1000, mxZi, 'mo', 'MarkerSize', 15, 'linewidth', 2 ) ;
                
            end
        
        end
  
        timeTitleText = sprintf('%s, Year=%d', titleText, round(t) );
        h = title(timeTitleText,'fontname','times','fontsize',16) ;
        set( h, 'visible', 'on' );
        ylabel('(km)','fontname','times','fontsize',16)
        xlabel('(km)','fontname','times','fontsize',16)
        axis([min(min(X)) max(max(X)) min(min(Y)) max(max(Y))]/1000) ;
        view(azView,elView)

    if ( zScale ~= 0 )
        axis equal
    end
    
    if ( exist('AxisPosition') )
        set( gca, 'Position', AxisPosition ) ;
    end
    
    drawnow

% catch
%         
%     if ( DEBUG_TOGGLE )
%         err = lasterror ;
%         errdepth = length(err.stack) ;
%   
%         disp( sprintf('\n%s\n', err.message) ) ;
%         disp( sprintf('error at line %d in %s', err.stack(1).line, err.stack(1).file ) ) ;
%         for i=2:errdepth
%             disp( sprintf('called at line %d in %s', err.stack(i).line, err.stack(i).file ) ) ;
%         end
%     end
% 
%     if ( RECONSTRUCT == 0 )
%     
%         cmd = sprintf( '%s( inputFile', mfilename ) ;
%         if ( nargin > 1 )
%             for n=1:2:length(varargin)-1
%                 cmd = sprintf( '%s, cat(1,varargin{%d}), cat(1,varargin{%d})', cmd, n, n+1 ) ;
%             end  
%         end
%         
%         property = 'RECONSTRUCT' ;
%         value = 1 ;
%         cmd = sprintf( '%s, property, %d)', cmd, value ) ;
%         disp(sprintf('\nFULL WORKSPACE MUST BE RECONSTRUCTED')) ;
%         eval(cmd) ;
%         
%     else
%         
%         err = lasterror ;
%         errdepth = length(err.stack) ;
%   
%         disp( sprintf('\n%s\n', err.message) ) ;
%         disp( sprintf('error at line %d in %s', err.stack(1).line, err.stack(1).file ) ) ;
%         for i=2:errdepth
%             disp( sprintf('called at line %d in %s', err.stack(i).line, err.stack(i).file ) ) ;
%         end
% 
%         error('unrecoverable error') ;  
%     end
%     
% end
                

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    END OF plot_gc2d.m     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

