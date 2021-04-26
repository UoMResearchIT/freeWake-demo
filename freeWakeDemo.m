close all;
clear
clc;
format short g;

% define simulation parameters
b    = 3;  % number of blades (typically 3-5)
Nr   = 25; % number of panels per blade (typically 15-30)
Npsi = 36; % number of time steps per revolution (typically 24 - 72)
revs = 01; % number of revolutions of wake to simulate (typically 3 - 10)

% wake nodes per blade
wnpb = (Nr+1) * (Npsi*revs + 1);
% total number of wake nodes
wn   = b*wnpb;

% create a random wake geometry
rows = Nr + 1;
cols = b;
pags = Npsi*revs + 1;
x = rand(rows,cols,pags);
y = rand(rows,cols,pags);
z = rand(rows,cols,pags);

% preallocate the induction at each of the wake nodes
[u, v, w] = deal( nan(size(x)) );

tic

% loop over each blade
for ii = 1:b
    
    % create a column vector of coordinates of the wake nodes of the ii-th
    % blades
    xp = reshape(x(:,ii,:), [], 1);
    yp = reshape(y(:,ii,:), [], 1);
    zp = reshape(z(:,ii,:), [], 1);
    
    % loop over the wake from each blade
    [ublade, vblade, wblade] = deal( zeros(length(xp),b) );
    for jj = 1:b
        
        % trailing filament induction
        [x1,y1,z1,x2,y2,z2] = trailed_nodes(x(:,jj,:),y(:,jj,:),z(:,jj,:));
        Gamma = rand(size(x1));
        [utrailed, vtrailed, wtrailed] = vortexLine(x1, y1, z1, x2, y2, z2, Gamma, xp, yp, zp);
        
        % shed filament induction
        [x1,y1,z1,x2,y2,z2] = shed_nodes(x(:,jj,:),y(:,jj,:),z(:,jj,:));
        Gamma = rand(size(x1));
        [ushed, vshed, wshed] = vortexLine(x1, y1, z1, x2, y2, z2, Gamma, xp, yp, zp);
        
        % sum the trailed and shed inductions
        ublade(:,jj) = sum(utrailed,2) + sum(ushed, 2);
        vblade(:,jj) = sum(vtrailed,2) + sum(vshed, 2);
        wblade(:,jj) = sum(wtrailed,2) + sum(wshed, 2);
        
    end
    
    % sum the induction from all the blades on the ii-th blade and reshape
    u(:,ii,:) = reshape( sum(ublade,2), [rows 1 pags] );
    v(:,ii,:) = reshape( sum(vblade,2), [rows 1 pags] );
    w(:,ii,:) = reshape( sum(wblade,2), [rows 1 pags] );
    
end

toc

function [x1,y1,z1,x2,y2,z2] = trailed_nodes(x,y,z)
    
    % start nodes
    x1 = reshape( x(:,:,1:end-1), 1, []);
    y1 = reshape( y(:,:,1:end-1), 1, []);
    z1 = reshape( z(:,:,1:end-1), 1, []);
    
    % end nodes
    x2 = reshape( x(:,:,2:end), 1, []);
    y2 = reshape( y(:,:,2:end), 1, []);
    z2 = reshape( z(:,:,2:end), 1, []);
    
end

function [x1,y1,z1,x2,y2,z2] = shed_nodes(x,y,z)
    
    % start nodes
    x1 = reshape( x(1:end-1,:,:), 1, []);
    y1 = reshape( y(1:end-1,:,:), 1, []);
    z1 = reshape( z(1:end-1,:,:), 1, []);
    
    % end nodes
    x2 = reshape( x(2:end,:,:), 1, []);
    y2 = reshape( y(2:end,:,:), 1, []);
    z2 = reshape( z(2:end,:,:), 1, []);
    
end


