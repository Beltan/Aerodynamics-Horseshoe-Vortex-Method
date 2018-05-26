function [ind] = induccion(x1,y1,z1,x2,y2,z2,xc,yc,zc)
  r0=[x2-x1,y2-y1,z2-z1];
  r1=[xc-x1,yc-y1,zc-z1];
  r2=[xc-x2,yc-y2,zc-z2];
  ind=(1/(4*pi))*dot(((1/norm(r1))*r1-(1/norm(r2))*r2),r0)*(1/((norm(cross(r1,r2)))^2))*cross(r1,r2);
end
