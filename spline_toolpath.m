function spline_toolpath

% Christopher Mannes
% Student ID: 20743505

% This function employs a cubic spline interpolation technique to generate
% a fan-shaped toolpath, then subsequently calls the function cont_traj to
% obtain arc length and feedrate profiles to employ natural and 1st order
% Taylor series trajectory interpolation.

% xp and yp give the coordinates of the waypoints for a fan-shape.
xp = 1.0*[  -10.7125  -0.0375   5.2475   9.0075  10.0325   9.0575   6.0750   5.4750   5.6100   6.1350   7.0275   8.0550   9.0550  10.1100  12.1300  17.0600  21.0975  22.8925  25.0675  27.5550 ...
             30.0850  33.0800  37.9900  38.1050  36.1225  32.0925  26.1375  20.1075  10.1025   7.9100   6.9450   6.0200   5.5975   5.5625   5.6875   6.0650   7.0975  12.2125  18.2025  20.1600 ...
             21.6450  22.3925  22.1600  20.1425  10.7125   0.0375  -5.2475  -9.0075 -10.0325  -9.0575  -6.0750  -5.4750  -5.6100  -6.1350  -7.0275  -8.0550  -9.0550 -10.1100 -12.1300 -17.0600 ...
            -21.0975 -22.8925 -25.0675 -27.5550 -30.0850 -33.0800 -37.9900 -38.1050 -36.1225 -32.0925 -26.1375 -20.1075 -10.1025  -7.9100  -6.9450  -6.0200  -5.5975  -5.5625  -5.6875  -6.0650 ...
             -7.0975 -12.2125 -18.2025 -20.1600 -21.6450 -22.3925 -22.1600 -20.1425 -10.7125 ]';

yp = 1.0*[  -37.9900 -38.1050 -36.1225 -32.0925 -26.1375 -20.1075 -10.1025  -7.9100  -6.9450  -6.0200  -5.5975  -5.5625  -5.6875  -6.0650  -7.0975 -12.2125 -18.2025 -20.1600 -21.6450 -22.3925 ...
            -22.1600 -20.1425 -10.7125  -0.0375   5.2475   9.0075  10.0325   9.0575   6.0750   5.4750   5.6100   6.1350   7.0275   8.0550   9.0550  10.1100  12.1300  17.0600  21.0975  22.8925 ...
             25.0675  27.5550  30.0850  33.0800  37.9900  38.1050  36.1225  32.0925  26.1375  20.1075  10.1025   7.9100   6.9450   6.0200   5.5975   5.5625   5.6875   6.0650   7.0975  12.2125 ...
             18.2025  20.1600  21.6450  22.3925  22.1600  20.1425  10.7125   0.0375  -5.2475  -9.0075 -10.0325  -9.0575  -6.0750  -5.4750  -5.6100  -6.1350  -7.0275  -8.0550  -9.0550 -10.1100 ...
            -12.1300 -17.0600 -21.0975 -22.8925 -25.0675 -27.5550 -30.0850 -33.0800 -37.9900 ]';
        

% Values for the desired feed profile.
v_0 = 0;
v_f = 100;
v_e = 0;
A = 1000;
D = 1000;
J = 30000;
T_dwell = 0;
T_s = 0.001;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Number of segments.
N = length(xp)-1;

% Cubic spline interpolation by chord length parametrization.

% The spline curve employed here is parametized with respect to the chord
% length (l_k) st. the spline parameter is bounded by u = 0 and u = l_k in
% each segment.

% Generation of a vector containing the chord lengths for each segment on
% the spline curve.
l = zeros(1,N);
for k = 1:N
    l(k) = sqrt( ( xp(k+1) - xp(k) )^2 + ( yp(k+1) - yp(k) )^2 );
end

% Renamed variable for natural interpolation.
L_k = l;
lk = l;

% Stacking of the position Boundary conditions (BCs) gives a 2Nx4N diagonal 
% matrix in which each element is 2x4 matrix.
counter1 = 0;
counter2 = 0;
L0 = zeros(2*N,4*N);
for i = 1:N
    L1 = [ 0 0 0 1; ( l(i) )^3 ( l(i) )^2 ( l(i) ) 1 ];
    L2 = vertcat(L1, zeros( (2*N) - 2 , 4 ));
    L3 = [L2 zeros( 2*N , (4*N) - 4 ) ];
    L4 = circshift(L3,[counter1,counter2]);
    L0 = L0 + L4;
    counter1 = counter1 + 2;
    counter2 = counter2 + 4;
end

% Stacking of the C^1 and C^2 continuity conditions to obtain a 2(N-1)x4N 
% matrix in which each element is a 2x8 matrix horizontally offset by 4 to 
% the right with respect to the preceeding matrix.
counter1 = 0;
counter2 = 0;
L12 = zeros(2*(N-1),4*N);
for i = 1:N-1
     L1 = [ 3*(( l(i) )^2)  2*( l(i) ) 1 0 0 0 -1 0; 6*( l(i) ) 2 0 0 0 -2 0 0 ];
     L2 = vertcat( L1,zeros( 2*(N-2), 8) );
     L3 = [ L2 zeros( 2*(N-1) , (4*N) - 8 ) ];
     L4 = circshift(L3,[counter1,counter2]);
     L12 = L12 + L4;
     counter1 = counter1 + 2;
     counter2 = counter2 + 4;
end

% Construction of a 2x4N matrix for the continuity of the first and second
% derivatives across the first/last control point in order to give a
% closed spline curve.
L1 = [ 0 0 -1 0; 0 -2 0 0];
LN = [ 3*(( l(N) )^2) 2*( l(N) ) 1 0; 6*( ( l(N) ) ) 2 0 0 ];
L_center = zeros( 2,( (4*N) - 8) );
L1N = [L1 L_center LN];

% Formation of the 4Nx4N matrix of equation relating the BCs to the 4
% unknowns of a cubic polynomial in each segment of the spline curve.
L = vertcat(L0,L12,L1N);

% Formation of the position solution vectors for x and y coordinates. Since
% the position solution vectors are the position boundary conditions in
% which the intermediate points are double counted to satify the segment of
% each interval. A vector is created with the positions separated by zero 
% elements, then a second identical vector is shifted by one element and 
% added to form the final vector.

% x-component.
zeta_x01 = zeros(2*N+2,1); 
for i = 2:2:2*N+2
    zeta_x01(i) = xp(i/2);
end

% Formation of the two appropriate staggered vectors.
zeta_x02 = zeta_x01(2:end-1,1);
zeta_x03 = vertcat( zeta_x01(3:end-2,1), 0,xp(end) );

% y-component.
zeta_y01 = zeros(2*N+2,1); 
for i = 2:2:2*N+2
    zeta_y01(i) = yp(i/2);
end

% Formation of the two appropriate staggered vectors.
zeta_y02 = zeta_y01(2:end-1,1);
zeta_y03 = vertcat( zeta_y01(3:end-2,1), 0,yp(end) )  ;

% Vector of solutions for the x component for  C^0, C^1, and C^2 boundary
% conditions.
zeta_x = [ ( zeta_x02 + zeta_x03 )' zeros(1,2*N) ]';

% Vector of solutions for the y component for  C^0, C^1, and C^2 boundary
% conditions.
zeta_y = [ ( zeta_y02 + zeta_y03 )' zeros(1,2*N) ]';

% 4Nx1 vector of coefficients for the x component. 
theta_x = L\zeta_x;

% Reshape the vector theta_x into a matrix in which the number of each row
% gives the coefficients corresponding to that segment and A is the first
% column, B is the second column, C in the third column, and D is the
% fourth column.
Cx = (reshape(theta_x,[4,N]))';

% 4Nx1 vector of coefficients for the y component. 
theta_y = L\zeta_y;

% Reshape the vector theta_x into a matrix in which the number of each row
% gives the coefficients corresponding to that segment and A is the first
% column, B is the second column, C in the third column, and D is the
% fourth column.
Cy = (reshape(theta_y,[4,N]))';

% Generation of a N(500)x1 vector for the x and y components for u with a
% domain between 0 and l_k for each k segment using the A, B, C, and D
% coefficents in Cx and Cy, respectively.
counter = 0;
M_k = 500;
x_spline = zeros( N*(M_k),1);
y_spline = zeros( N*(M_k),1);
for i = 1:N
    u = linspace( 0,l(i),M_k );
    for j = 1:length(u)
        counter = counter + 1;
        x_spline(counter,1) = Cx(i,1)*(( u(j) )^3) + Cx(i,2)*(( u(j) )^2) + Cx(i,3)*( u(j) ) + Cx(i,4);
        y_spline(counter,1) = Cy(i,1)*(( u(j) )^3) + Cy(i,2)*(( u(j) )^2) + Cy(i,3)*( u(j) ) + Cy(i,4);
    end
end

% Conversion of x_spline from a single column vector into a matrix with N
% rows and M_k columns. 
x_spline_matrix = ( reshape(x_spline,[M_k,N]) )';

% Conversion of y_spline from a single column vector into a matrix with N
% rows and M_k columns. 
y_spline_matrix = ( reshape(y_spline,[M_k,N]) )';

% The arc length of each segment is obtained by taking the square root of
% the sum of the squares for the corresponding increments in x and y.
s_k = zeros(N,1);
for l = 1:N
    for m = 1:M_k-1
    s_k(l) = s_k(l) + sqrt( ( x_spline_matrix(l,m+1) - x_spline_matrix(l,m) )^2 + ( y_spline_matrix(l,m+1) - y_spline_matrix(l,m) )^2 );
    end
end

% Total arc length distance.
s_L = sum(s_k);

% Feed profile for the arc length.
[s_toolpath, s_toolpath_dot] = cont_traj(s_L,v_0,v_f,v_e,A,D,J,T_dwell);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Natural cubic spline interpolation.

% Linear mapping by natural interpoalation between u and s calculated for 
% the first interval. Each value of s_toolpath that occurs in the first
% interval is scaled by the ratio of the chord length over the arc length.
% The resulting x and y components are then calculated using the cubic
% polynomial.
counter = 0;
for i = 1:length(s_toolpath)
    if s_toolpath(i) < s_k(1)
        counter = counter + 1;
        U(counter) = ( lk(1)/s_k(1) )*s_toolpath(i);
        u = ( lk(1)/s_k(1) )*s_toolpath(i);
        x_natural(counter,1) = Cx(1,1)*(( u )^3) + Cx(1,2)*(( u )^2) + Cx(1,3)*( u ) + Cx(1,4);
        y_natural(counter,1) = Cy(1,1)*(( u )^3) + Cy(1,2)*(( u )^2) + Cy(1,3)*( u ) + Cy(1,4);
    else
        break
    end
end

% Linear mapping by natural interpolation of the remaining N-1 intervals.
% The interval for each element is determined by iterating the sum of the
% arc lengths in reverse, that is total arc length to the arc length of the
% first interval and checking when the toolpath distance is greater.
for i = 183:length(s_toolpath)
    for j = length(s_k):-1:2
        if s_toolpath(i) < sum( s_k(1:j-1) )
            continue
        else
            counter = counter + 1;
            s_prime = s_toolpath(i) - sum(s_k(1:j-1));
            U(counter) = ( lk(j)/s_k(j) )*s_prime;
            u = ( lk(j)/s_k(j) )*s_prime;
            x_natural(counter,1) = Cx(j,1)*(( u )^3) + Cx(j,2)*(( u )^2) + Cx(j,3)*( u ) + Cx(j,4);
            y_natural(counter,1) = Cy(j,1)*(( u )^3) + Cy(j,2)*(( u )^2) + Cy(j,3)*( u ) + Cy(j,4);
            break
        end
    end
end


% Numerical differentiation wrt. time.

% Numerically differentiated x position.
x_dot_natural = zeros(1,length(x_natural)-2);
for i = 1:length(x_natural)-2
    x_dot_natural(i) = ( x_natural(i+2) - x_natural(i) )/(2*T_s);
end

% Numerically differentiated y position.
y_dot_natural = zeros(1,length(x_natural)-2);
for i = 1:length(x_natural)-2
    y_dot_natural(i) = ( y_natural(i+2) - y_natural(i) )/(2*T_s);
end

% Feedrate obtained by natural interpolation.
f_natural = zeros(1,length(x_natural)-2);
for i = 1:length(x_natural)-2
    f_natural(i) = sqrt( x_dot_natural(i)^2 + y_dot_natural(i)^2 );
end

% Numerically differentiation of the x component of velocity.
x_ddot_natural = zeros(1,length(x_natural)-4);
for i = 1:length(x_natural)-4
    x_ddot_natural(i) = ( x_dot_natural(i+2) - x_dot_natural(i) )/(2*T_s);
end

% Numerically differentiation of the y component of velocity.
y_ddot_natural = zeros(1,length(x_natural)-4);
for i = 1:length(x_natural)-4
    y_ddot_natural(i) = ( y_dot_natural(i+2) - y_dot_natural(i) )/(2*T_s);
end

% Acceleration obtained from natural Interpolation.
A_natural = zeros(1,length(x_natural)-4);
for i = 1:length(x_natural)-4
    A_natural(i) = sqrt( x_ddot_natural(i)^2 + y_ddot_natural(i)^2 );
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% 1st Order Taylor Series Interpolation.

% The value of the toolpath distance is incrementally increased and saved
% in s_prime. The value of s_prime is compared with the arc length of the
% corresponding interval. If the toolpath distance in less than the arc
% length interval then the toolpath distance is mapped to the spline
% parameter u by the same equation as natural interpolation and
% subsequently used to calculate dx/du and dy/du. The resulting derivative
% values are used to obtain u_dot for the previous time step and added to 
% the value of u at the previous time step. The resulting u is then used to
% obtain x_taylor and y_taylor. Lastly, u_km1 is reset to the current u
% value.

% Initialize parameters.
j = 1; 
u_km1 = 0;     
s_prime = 0;
counter = 0;
x_taylor(1) = xp(1);
y_taylor(1) = yp(1);
for i = 2:length(s_toolpath)-34
    s_prime =  s_prime + (s_toolpath(i) - s_toolpath(i-1));
    if s_prime > s_k(j)
        s_prime = s_prime - s_k(j);
        j = j + 1;
        u_km1 = u_km1 - lk(j-1);
    end
    counter = counter + 1;
    u_t = ( lk(j)/s_k(j) )*s_prime;
    xu_km1 = 3*Cx(j,1)*u_t^2 + 2*Cx(j,2)*u_t + Cx(j,3);
    yu_km1 = 3*Cy(j,1)*u_t^2 + 2*Cy(j,2)*u_t + Cy(j,3);
    u_t = u_km1 + ( s_toolpath_dot(i)*T_s )/( sqrt( xu_km1^2 + yu_km1^2 ) );
    x_taylor(counter,1) = Cx(j,1)*(( u_t )^3) + Cx(j,2)*(( u_t )^2) + Cx(j,3)*( u_t ) + Cx(j,4);
    y_taylor(counter,1) = Cy(j,1)*(( u_t )^3) + Cy(j,2)*(( u_t )^2) + Cy(j,3)*( u_t ) + Cy(j,4);
    u_km1 = u_t;
end

% Due to numerical issues, the last 34 points were obtained by a separate
% For loop.
for i = length(s_toolpath)-34:length(s_toolpath)
    counter = counter + 1;
    s_prime =  s_prime + (s_toolpath(i) - s_toolpath(i-1));
    u_t = ( lk(N)/s_k(N) )*s_prime;
    xu_km1 = 3*Cx(j,1)*u_t^2 + 2*Cx(j,2)*u_t + Cx(j,3);
    yu_km1 = 3*Cy(j,1)*u_t^2 + 2*Cy(j,2)*u_t + Cy(j,3);
    u_t = u_km1 + ( s_toolpath_dot(i)*T_s )/( sqrt( xu_km1^2 + yu_km1^2 ) );
    x_taylor(counter,1) = Cx(N,1)*(( u_t )^3) + Cx(N,2)*(( u_t )^2) + Cx(N,3)*( u_t ) + Cx(N,4);
    y_taylor(counter,1) = Cy(N,1)*(( u_t )^3) + Cy(N,2)*(( u_t )^2) + Cy(N,3)*( u_t ) + Cy(N,4);
    u_km1 = u_t;
end


% Numerical differentiation wrt. time.

% Numerically differentiation of x.
x_dot_taylor = zeros(1,length(s_toolpath)-2);
for i = 1:length(x_taylor)-2
    x_dot_taylor(i) = ( x_taylor(i+2) - x_taylor(i) )/(2*T_s);
end

% Numerically differentiation of y.
y_dot_taylor = zeros(1,length(x_taylor)-2);
for i = 1:length(x_taylor)-2
    y_dot_taylor(i) = ( y_taylor(i+2) - y_taylor(i) )/(2*T_s);
end

% Feedrate obtained from 1st Order Taylor Series Interpolation.
f_taylor = zeros(1,length(x_taylor)-2);
for i = 1:length(x_taylor)-2
    f_taylor(i) = sqrt( x_dot_taylor(i)^2 + y_dot_taylor(i)^2 );
end

% Numerically differentiation of the x component of velocity.
x_ddot_taylor = zeros(1,length(x_taylor)-4);
for i = 1:length(x_taylor)-4
    x_ddot_taylor(i) = ( x_dot_taylor(i+2) - x_dot_taylor(i) )/(2*T_s);
end

% Numerically differentiation of the x component of velocity.
y_ddot_taylor = zeros(1,length(x_taylor)-4);
for i = 1:length(x_taylor)-4
    y_ddot_taylor(i) = ( y_dot_taylor(i+2) - y_dot_taylor(i) )/(2*T_s);
end

% Acceleration obtained from 1st Order Taylor Series Interoolation.
A_taylor = zeros(1,length(x_taylor)-4);
for i = 1:length(x_taylor)-4
    A_taylor(i) = sqrt( x_ddot_taylor(i)^2 + y_ddot_taylor(i)^2 );
end

% Jerk obtained from 1st Order Taylor Series Interoolation.
df = zeros(1,length(x_taylor)-4);
for i = 1:length(x_taylor)-4
    df(i) = ( f_taylor(i+2) - f_taylor(i) )/(2*T_s);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Plotting

T = T_s:T_s:T_s*length(s_toolpath);

% Fan-shaped trajectory
figure(2)
set(gcf,'color','white')
box on

hold on
plot(xp,yp,'.k','MarkerSize',25)
plot(x_spline,y_spline,'-y','Linewidth',6)
plot(x_natural,y_natural,'-.b','Linewidth',2)
plot(x_taylor,y_taylor,'--m','Linewidth',2)
plot(xp,yp,'.k','MarkerSize',15)
% title('Fan-Shaped Cubic Spline Trajectory')
xlabel('x-Coordinate (mm)','fontsize',16)
ylabel('y-Coordinate (mm)','fontsize',16)
legend('Control Points','Cubic Spline Toolpath','Natural Interpolation','1^{st} Order Taylor Interpolation','location','Northwest')
hold off

% x component of position. 
figure(3)
set(gcf,'color','white')
box on

plot(T,x_natural,'-b',T,x_taylor,'--m','Linewidth',2)
% title('x-Coordinate of Position Profile')
xlabel('Time (ms)','fontsize',20)
ylabel('x-Coordinate (mm)','fontsize',20)
legend('Natural Interpolation','1^{st} Order Taylor Series Interpolation')

% y component of position.  
figure(4)
set(gcf,'color','white')
box on

plot(T,y_natural,'-b',T,y_taylor,'--m','Linewidth',2)
% title('y-Coordinate of Position Profile')
xlabel('Time (ms)')
ylabel('y-Coordinate (mm)')
legend('Natural Interpolation','1^{st} Order Taylor Series Interpolation','location','Northeast')

% x component of velocity. 
figure(5)
set(gcf,'color','white')
box on

plot(T(2:end-1),x_dot_natural,'-b',T(2:end-1),x_dot_taylor,'--m','Linewidth',2)
% title('x Component of Velocity Profile')
xlabel('Time (ms)')
ylabel('x Component of Velocity (mm/s)')
legend('Natural Interpolation','1^{st} Order Taylor Series Interpolation','location','Northeast')

% y component of velocity. 
figure(6)
set(gcf,'color','white')
box on

plot(T(2:end-1),y_dot_natural,'-b',T(2:end-1),y_dot_taylor,'--m','Linewidth',2)
% title('y Component of Velocity Profile')
xlabel('Time (ms)')
ylabel('y Component of Velocity (mm/s)')
legend('Natural Interpolation','1^{st} Order Taylor Series Interpolation','location','Northeast')

% Feedrate. 
figure(7)
set(gcf,'color','white')
box on

plot(T(2:end-1),f_natural,'-b',T(2:end-1),f_taylor,'--m','Linewidth',2)
% title('Feedrate Profile of a Fan-Shaped Trajectory')
xlabel('Time (ms)')
ylabel('Tangential Velocity (mm/s)')
legend('Natural Interpolation','1^{st} Order Taylor Series Interpolation','location','Northeast')

% x component of acceleration. 
figure(8)
set(gcf,'color','white')
box on

plot(T(3:end-2),x_ddot_natural,'-b',T(3:end-2),x_ddot_taylor,'--m','Linewidth',2)
% title('x Component of Acceleration Profile')
xlabel('Time (ms)')
ylabel('x Component of Acceleration (mm/s^3)')
legend('Natural Interpolation','1^{st} Order Taylor Series Interpolation','location','Northeast')

% y Component of acceleration. 
figure(9)
set(gcf,'color','white')
box on

plot(T(3:end-2),y_ddot_natural,'-b',T(3:end-2),y_ddot_taylor,'--m','Linewidth',2)
% title('y Component of Acceleration Profile')
xlabel('Time (ms)')
ylabel('y Component of Acceleration (mm/s^3)')
legend('Natural Interpolation','1^{st} Order Taylor Series Interpolation','location','Northeast')

% Acceleration. 
figure(10)
set(gcf,'color','white')
box on

plot(T(3:end-2),A_natural,'-b',T(3:end-2),A_taylor,'--m','Linewidth',2)
% title('Jerk Profile of a Fan-Shaped Trajectory')
xlabel('Time (ms)')
ylabel('Acceleration (mm/s^3)')
legend('Natural Interpolation','1^{st} Order Taylor Series Interpolation','location','Northeast')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Subplot format for report.
figure(11)
set(gcf,'color','white')

sp(1) = subplot(2,2,1);
plot(T,x_natural,'-b',T,x_taylor,'--m','Linewidth',2)
xlabel('Time (ms)','fontsize',15)
ylabel('x-Coordinate (mm)','fontsize',15)

sp(2) = subplot(2,2,2);
plot(T,y_natural,'-b',T,y_taylor,'--m','Linewidth',2)
xlabel('Time (ms)','fontsize',15)
ylabel('y-Coordinate (mm)','fontsize',15)

sp(3) = subplot(2,2,3);
plot(T(2:end-1),x_dot_natural,'-b',T(2:end-1),x_dot_taylor,'--m','Linewidth',2)
xlabel('Time (ms)','fontsize',15)
ylabel('x Velocity (mm/s)','fontsize',15)

sp(4) = subplot(2,2,4);
plot(T(2:end-1),y_dot_natural,'-b',T(2:end-1),y_dot_taylor,'--m','Linewidth',2)
xlabel('Time (ms)','fontsize',15)
ylabel('y Velocity (mm/s)','fontsize',15)
legend('Natural Interpolation','1^{st} Order Taylor Series','location','Southoutside')

figure(12)
set(gcf,'color','white')

sp(5) = subplot(2,1,1);
plot(T(2:end-1),f_natural,'-b',T(2:end-1),f_taylor,'--m','Linewidth',2)
xlabel('Time (ms)','fontsize',15)
ylabel('Feedrate (mm/s^2)','fontsize',15)
legend('Natural Interpolation','1^{st} Order Taylor Series','location','Southoutside')

sp(6) = subplot(2,1,2);
plot(T(2:end-1),f_natural,'-b',T(2:end-1),f_taylor,'--m','Linewidth',2)
xlabel('Time (ms)','fontsize',15)
ylabel('Feedrate (mm/s^2)','fontsize',15)
axis([ 0.2 3.665 97.5 102.5 ])

figure(13)
set(gcf,'color','white')

sp(7) = subplot(3,1,1);
plot(T(3:end-2),x_ddot_natural,'-b',T(3:end-2),x_ddot_taylor,'--m','Linewidth',2)
xlabel('Time (ms)','fontsize',15)
ylabel('x Acceleration (mm/s^2)','fontsize',15)

sp(8) = subplot(3,1,2);
plot(T(3:end-2),y_ddot_natural,'-b',T(3:end-2),y_ddot_taylor,'--m','Linewidth',2)
xlabel('Time (ms)','fontsize',15)
ylabel('y Acceleration (mm/s^2)','fontsize',15)

sp(9) = subplot(3,1,3);
plot(T(3:end-2),A_natural,'-b',T(3:end-2),A_taylor,'--m','Linewidth',2)
xlabel('Time (ms)','fontsize',15)
ylabel('Accelerartion (mm/s^2)','fontsize',15)
legend('Natural Interpolation','1^{st} Order Taylor Series','location','Southoutside')




end

