%function gis_2km_jigsaw()
% 1km gis mesh

    name = 'gis1km';

    %addpath('dual-mesh');

%------------------------------------ setup files for JIGSAW

    opts.geom_file = ...                % GEOM file
    [name '-GEOM.msh'];
    
    opts.jcfg_file = ...                % JCFG file
    [name '.jig'];
    
    opts.mesh_file = ...                % MESH file
    [name '-MESH.msh'];
    
    opts.hfun_file = ...                % HFUN file
    [name '-HFUN.msh'];
    
%------------------------------------ define JIGSAW geometry

    geom.mshID = 'EUCLIDEAN-MESH';

    x0=-750000.0; x1=1000000.0;
    y0=-3400000.0; y1=-600000.0;
    
    buf=50000.0
    x0=-750000.0+buf; x1=1000000.0-buf;
    y0=-3400000.0+buf; y1=-600000.0-buf;

    
    geom.point.coord = [    % list of xy "node" coordinates
        x0, y0, 0
        x1, y0, 0
        x1, y1, 0
        x0, y1, 0 ] ;
    
    geom.edge2.index = [    % list of "edges" between nodes
        1, 2, 0
        2, 3, 0
        3, 4, 0
        4, 1, 0  ] ;    
        
    savemsh(opts.geom_file, geom) ;
    
%------------------------------------ compute HFUN over GEOM
        
%%  read density function from a file 
% densFile='density.nc';
% xpos=ncread(densFile, 'x');
% ypos=ncread(densFile, 'y');
% dens = ncread(densFile, 'density')';
% hfun = dens.^-0.25*1000.0 * 10.0;
% 
% 
% 
% 
%    [XPOS,YPOS] = meshgrid(xpos,ypos) ;
%     
% 
%     hmat.mshID = 'EUCLIDEAN-GRID' ;
%     hmat.point.coord{1} = xpos ;
%     hmat.point.coord{2} = ypos ;
%     hmat.value = hfun ;
% 
%         
%     savemsh(opts.hfun_file,hmat) ;

    
%% calculate density based on new criteria related to geometry

%% load GIS data file

gisFile = '/Users/mhoffman/Documents/greenland_geometry/jkennedy_20190227/greenland_1km_2019_02_11.epsg3413.nc';
x=ncread(gisFile, 'x1');
y=ncread(gisFile, 'y1');
thk = ncread(gisFile, 'thk')';
topg = ncread(gisFile, 'topg')';
vx = ncread(gisFile, 'vx')';
vy = ncread(gisFile, 'vy')';

dx = x(2)-x(1);
nx = length(x);
ny = length(y);


%% Decimate fields
% mesh decimation. Keep a point every nsteps (every nsteps indices)
% if the initial mesh is 1km resol., then the resolution of the uniform mesh will be nsteps. 
nsteps = 16;

x = x(1:nsteps:end);
y = y(1:nsteps:end);
vx = vx(1:nsteps:end, 1:nsteps:end, 1);
vy = vy(1:nsteps:end, 1:nsteps:end, 1);
thk = thk(1:nsteps:end, 1:nsteps:end, 1);
topg = topg(1:nsteps:end, 1:nsteps:end, 1);

dx = x(2)-x(1);
nx = length(x);
ny = length(y);


[YPOS,XPOS] = meshgrid(x,y);

maskSize = size(thk);




%%  Calculate masks


ice_tol = 1.0
icemask = int32(thk>ice_tol);

% floating mask - defined as floatation criterion + bed below sea level + presence of ice
rhoi = 910.0;
rhoo = 1027.0; 
floatmask = (rhoi * thk) ./ (rhoo * -1.0 * topg) < 1.0 & topg <0.0 & icemask;

neighbors=[[1,0]; [-1,0]; [0,1]; [0,-1]]';
% ice margin mask
marginMask = icemask*0;
iceMask = thk>0;
for n=neighbors;
   marginMask = marginMask | ~(circshift(iceMask, n));
end
marginMask = marginMask & iceMask;  % where ice exists and neighbors non-ice locations



% plot masks
figure(10); clf; 
subplot(1,3,1); hold all; axis equal
imagesc(icemask); colorbar; title('ice mask')
subplot(1,3,2); hold all; axis equal
imagesc(floatmask); colorbar; title('float mask')
subplot(1,3,3); hold all; axis equal
imagesc(marginMask); colorbar; title('margin mask')

%% calculate distance to margin
distToMargin = thk*0.0+5e4;

% -- KEY PARAMETER: how big of a search 'box' (one-directional) to use.
% Bigger number makes search slower, but if too small, the transition zone
% could get truncated.
% (could automatically set this from maxDist variables used in next section.)
windowSize = 200.0e3;
% ---

d = int32(ceil(windowSize / dx));
%d=80;
rng = [-1*d:d];
maxdist = double(d) * dx * 1;


for i=d+1:nx-d-1
    for j=d+1:ny-d-1

    irng = i+rng;
    jrng = j+rng;

   % irng = irng(find(irng>0 & irng < nx));
   % jrng = jrng(find(jrng>0 & jrng < ny));

    dist2Here = ((XPOS(jrng,irng)-y(j)).^2 + (YPOS(jrng,irng)-x(i)).^2).^0.5;
    dist2Here(marginMask(jrng,irng)==0) = maxdist;  % put large value everywhere that isn't a margin
    distToMargin(j,i) = min(min(dist2Here(:)));
    end
end


figure (33); clf;
pcolor(distToMargin), shading flat; colorbar, axis equal

%% Calculate velo gradient

% get mesh info
dx = x(2)-x(1);
dy = y(2)-y(1);

% calculate gradient
[v1x,v1y] = gradient(vx, dx, dy); 
[v2x,v2y] = gradient(vy, dx, dy); 
dv2 = (v1x.^2 + v1y.^2 + v2x.^2 + v2y.^2);

dv2 = del2( (vx.^2+vy.^2).^0.5 , dx);
dv2 = (vx.^2+vy.^2).^0.5;
dv2(icemask==0)=4e2;
% plot the gradient
figure(20); clf; 
subplot(1,4,1);
pcolor((dv2)); shading flat; colorbar; title('(dv2)');  axis equal
caxis([5e-1,10^2])

subplot(1,4,2);
pcolor(log10(dv2)); shading flat; colorbar; title('log10(dv2)');  axis equal
caxis([-0,3])

% print max/min values
max(max(dv2))
min(min(dv2))


%%  QAQC - fill in missing data with something reasonable for now
spacing = log10(dv2);
figure(20);  
subplot(1,4,3); axis equal
pcolor(spacing); shading flat; colorbar; title('spacing');  axis equal


minvalue = 0.5;
maxvalue = 2.5;

dmin = 1.0;
dmax=8.0;

spacing(spacing<minvalue)=minvalue;
spacing(spacing>maxvalue)=maxvalue;
spacing(icemask==0)=maxvalue - (maxvalue-minvalue)*0.3;


%spacingClean = interp1([minvalue, maxvalue], [dmin, dmax], spacing);
spacingClean = interp1([maxvalue, minvalue], [dmin, dmax], spacing);

rimDist = 50.0e3;
ind = distToMargin<rimDist;
newvals = dmin*(dmax-dmin)*0.2 + (dmax-dmin)*0.3 * (max(0.0, distToMargin(ind)-rimDist*0.25)/rimDist);  % put min value near margin with fade
spacingClean(ind) = min(newvals , spacingClean(ind))

figure(20);  
subplot(1,4,4); axis equal
pcolor(spacingClean); shading flat; colorbar; title('spacingClean');  axis equal


hfun = spacingClean *1000.0 * 4.0;   
    
figure(96); clf; hold all
pcolor(hfun/1000.0)
shading flat
colorbar
axis equal


%% save to jigsaw format        

    hmat.mshID = 'EUCLIDEAN-GRID' ;
    hmat.point.coord{1} = x;
    hmat.point.coord{2} = y ;
    hmat.value = hfun ;


    savemsh(opts.hfun_file,hmat) ;


%% ------------------------------------ build mesh via JIGSAW! 
  
    opts.hfun_scal = 'absolute';
    opts.hfun_hmax = +inf ;             % null HFUN limits
    opts.hfun_hmin = 0.00 ;
  
    opts.mesh_dims = +2 ;               % 2-dim. simplexes
    
    opts.optm_qlim = 0.9375 ;
   
    opts.mesh_top1 = true ;             % for sharp feat's
    opts.geom_feat = true ;
    
    mesh = jigsaw  (opts) ;
 
%------------------------------------ draw mesh/cost outputs

    ang2 = triang2( ...                 % calc. tri-angles
        mesh.point.coord(:,1:2), ...
        mesh.tria3.index(:,1:3)) ;
            
    t_90 = max(ang2,[],2) > 90.0 ;
    t_95 = max(ang2,[],2) > 95.0 ;
    
    figure(1); clf;
    patch ('faces',geom.edge2.index(:,1:2), ...
        'vertices',geom.point.coord(:,1:2), ...
        'facecolor','w', ...
        'edgecolor',[.1,.1,.1], ...
        'linewidth',1.5) ;
    hold on; axis image;
    title('JIGSAW GEOM data') ;
%%
    figure(2); clf; hold all;
    pcolor(XPOS,YPOS,(hmat.value));
    axis equal; 
    shading interp ;
    title('JIGSAW HFUN data') ;
        patch ('faces',geom.edge2.index(:,1:2), ...
        'vertices',geom.point.coord(:,1:2), ...
        'facecolor','w', ...
        'edgecolor',[1,.1,.8], ...
        'linewidth',1.5) ;
    colorbar();
%%

    figure(3); clf;
    patch ('faces',mesh.tria3.index(:,1:3), ...
        'vertices',mesh.point.coord(:,1:2), ...
        'facecolor','w', ...
        'edgecolor',[.2,.2,.2]) ;

    hold on; axis image;
    patch ('faces',mesh.tria3.index(t_90,1:3), ...
        'vertices',mesh.point.coord(:,1:2), ...
        'facecolor','y', ...
        'edgecolor',[.2,.2,.2]) ;
    patch ('faces',mesh.tria3.index(t_95,1:3), ...
        'vertices',mesh.point.coord(:,1:2), ...
        'facecolor','r', ...
        'edgecolor',[.2,.2,.2]) ;
    patch ('faces',mesh.edge2.index(:,1:2), ...
        'vertices',mesh.point.coord(:,1:2), ...
        'facecolor','w', ...
        'edgecolor',[.1,.1,.1], ...
        'linewidth',1.5) ;
    axis equal
    title('JIGSAW TRIA mesh') ;
%%
    drawscr(mesh.point.coord (:,1:2), ...
            mesh.edge2.index (:,1:2), ...
            mesh.tria3.index (:,1:3)) ;
    
    drawnow ;        
%     set(figure(1),'units','normalized', ...
%         'position',[.05,.55,.30,.35]) ;
%     set(figure(2),'units','normalized', ...
%         'position',[.35,.55,.30,.35]) ;
%     set(figure(3),'units','normalized', ...
%         'position',[.35,.10,.30,.35]) ;
%     set(figure(4),'units','normalized', ...
%         'position',[.05,.10,.30,.35]) ;
    drawnow ;

%end
