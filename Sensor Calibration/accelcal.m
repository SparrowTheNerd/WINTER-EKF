% Assumes sensor axes are orthogonal bc I don't wanna do matrix soft calibration

clear; clc;
accelDat = readtable('Accel.csv');
accelArray = table2array(accelDat);

[Center_LSE,Radius_LSE] = sphereFit(accelArray);
% [X, Y, Z] = sphere(50);

x = accelDat.X; y = accelDat.Y; z = accelDat.Z;

% do the fitting (From https://www.mathworks.com/matlabcentral/fileexchange/24693-ellipsoid-fit)
[ center, radii, evecs, v, chi2 ] = ellipsoid_fit( [ x y z ], '0' );
fprintf( 'Sphere center: %.5g %.5g %.5g\n', Center_LSE);
fprintf( 'Sphere radius: %.5g\n============\n', Radius_LSE);
fprintf( 'Ellipsoid center: %.5g %.5g %.5g\n', center );
fprintf( 'Ellipsoid radii: %.5g %.5g %.5g\n', radii );
% fprintf( 'Ellipsoid evecs:\n' );
% fprintf( '%.5g %.5g %.5g\n%.5g %.5g %.5g\n%.5g %.5g %.5g\n', ...
%     evecs(1), evecs(2), evecs(3), evecs(4), evecs(5), evecs(6), evecs(7), evecs(8), evecs(9) );
% fprintf( 'Algebraic form:\n' );
% fprintf( '%.5g ', v );
fprintf( 'Average deviation of the fit: %.5f\n', sqrt( chi2 / size( x, 1 ) ) );
fprintf( '\n' );

% calibration parameters
scale = 9.80665./radii;
ofst = center;

[centerC, radiusC] = sphereFit((accelArray-ofst').*scale');
[ centerEC, radiiC, evecsC, ~ , chi2C ] = ellipsoid_fit( accelArray.*scale'-ofst', '0');
fprintf( 'Corrected sphere center: %.5g %.5g %.5g\n', centerC);
fprintf( 'Corrected sphere radius: %.5g\n============\n', radiusC);
fprintf( 'Corrected ellipsoid center: %.5g %.5g %.5g\n', centerEC );
fprintf( 'Corrected ellipsoid radii: %.5g %.5g %.5g\n', radiiC );
fprintf( 'Average deviation of the corrected fit: %.5f\n\n=============\n\n', sqrt( chi2C / size( x, 1 ) ) );

fprintf( 'Correction Parameters A*s-o:\nXo: %.5g  Yo: %.5g  Zo: %.6g\nXs: %.5g  Ys: %.5g  Zs: %.5g\n',ofst(1),ofst(2),ofst(3),scale(1),scale(2),scale(3));

% draw  data
figure,
plot3( x, y, z, '.r' );
hold on;
%draw fit
mind = min( [ x y z ] );
maxd = max( [ x y z ] );
nsteps = 50;
step = ( maxd - mind ) / nsteps;
[ x, y, z ] = meshgrid( linspace( mind(1) - step(1), maxd(1) + step(1), nsteps ), linspace( mind(2) - step(2), maxd(2) + step(2), nsteps ), linspace( mind(3) - step(3), maxd(3) + step(3), nsteps ) );
Ellipsoid = v(1) *x.*x +   v(2) * y.*y + v(3) * z.*z + ...
          2*v(4) *x.*y + 2*v(5)*x.*z + 2*v(6) * y.*z + ...
          2*v(7) *x    + 2*v(8)*y    + 2*v(9) * z;
p = patch( isosurface( x, y, z, Ellipsoid, -v(10) ) );
hold off;
set( p, 'FaceColor', 'g', 'EdgeColor', 'none' );
axis equal;
grid on;
camlight;
lighting phong;
xlabel('X'); ylabel('Y'); zlabel('Z');

% % draw data
% plot3(accelDat.X,accelDat.Y,accelDat.Z,'.r');
% hold on;
% % surf(X*Radius_LSE+Center_LSE(1),Y*Radius_LSE+Center_LSE(2),Z*Radius_LSE+Center_LSE(3),'faceAlpha',0.3,'Facecolor','b','EdgeAlpha',0.2);
% axis square;
% grid on;
