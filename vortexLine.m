function [u, v, w] = vortex_line(x1, y1, z1, x2, y2, z2, Gamma, xp, yp, zp)
    
    % r1 vector (from x1 to xp)
    r1x = xp - x1;
    r1y = yp - y1;
    r1z = zp - z1;
    r1  = sqrt(r1x.^2 + r1y.^2 + r1z.^2);
    
    % r2 vector (from x2 to xp)
    r2x = xp - x2;
    r2y = yp - y2;
    r2z = zp - z2;
    r2  = sqrt(r2x.^2 + r2y.^2 + r2z.^2);
    
    % r0 vector (from x1 to x2)
    r0x = r1x - r2x;
    r0y = r1y - r2y;
    r0z = r1z - r2z;
    r0  = sqrt(r0x.^2 + r0y.^2 + r0z.^2);
    
    % cross product: r3 = r1 x r2
    r3x = r1y.*r2z - r1z.*r2y;
    r3y = r1z.*r2x - r1x.*r2z;
    r3z = r1x.*r2y - r1y.*r2x;
    r3  = sqrt(r3x.^2 + r3y.^2 + r3z.^2);
    
    % perpendicular distance
    h = r3./r0;
    
    % cosine terms
    cos1 = (r0x.*r1x + r0y.*r1y + r0z.*r1z )./(r0.*r1);
    cos2 = (r0x.*r2x + r0y.*r2y + r0z.*r2z )./(r0.*r2);
    
    % Cartesian unit vectors of the swirl profile
    ux = r3x ./ r3;
    uy = r3y ./ r3;
    uz = r3z ./ r3;
    
    % velocity components at p induced by each filament
    B = Gamma.*(cos1 - cos2)./(4*pi*h);
    u = B .* ux;
    v = B .* uy;
    w = B .* uz;
    
    % set undefined induced velocities to zero
    u(~isfinite(u)) = 0;
    v(~isfinite(v)) = 0;
    w(~isfinite(w)) = 0;
    
end