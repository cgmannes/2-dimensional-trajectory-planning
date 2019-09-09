function [ s,s_dot,s_ddot,s_tdot ] = cont_traj(d,v_0,v_f,v_e,A,D,J,T_dwell)

% Christopher Mannes
% Student ID: 20743505

% This function calculates and plots the distance, velocity, acceleration,
% and jerk profiles for the given inputs.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Inputs: 
%
% d is the length of the commanded position.
% v_0 is the initial velocity.
% v_f is the final velocity.
% v_e is the exit velocity.
% A is the acceleration magnitude.
% D is the decceleration magnitude.
% J is the jerk magnitude.
% T_dwell is the time interval in the position d=L is constant.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 7
    error('This function requires a minimum of 7 inputs.')
else if nargin == 7
        T_dwell = 0;
    else if nargin >8
        error('This function accepts a maximum of 8 inputs.')
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Determination of the jerk, acceleration, and velocity time intervals.

% Time step.
T_s = 0.001;

% Jerk time interval.
TJ = (3*A)/(2*J);

% Constant acceleration time interval.
TA = ( (v_f - v_0)/A ) - TJ;

% Constant velocity time interval.
Tv = (1/v_f)*( abs(d) - (1/2)*(v_0 + v_f)*(2*TJ + TA) - (1/2)*(v_f + v_e)*(2*TJ + TA) );

% Constant decceleration time interval.
TD = ( (v_f - v_e)/D ) - TJ;

% Verification that commanded position is sufficiently long to acheive the
% desired maximum feedrate (tangential velocity).
l_123 = ( (v_0 + v_f)/2 )*(TJ + TA + TJ);
l_567 = ( (v_f + v_e)/2 )*(TJ + TD + TJ);

% Implementation of a triangular velocity profile for commanded positions
% that are insufficiently long to achieve the final velocity for a
% specified acceleration.
% if abs(d) < A*( ( (v_f - v_0)/A )^2 );
if (l_123 + l_567) > abs(d);
    Tv = 0;
    TA = sqrt(abs(d)/A) - TJ;
    TD = sqrt(abs(d)/D) - TJ;
end

% Total time required to move a distance L.
TL = 4*TJ + TA + Tv + TD;

% Verification that the time required to move a distance, L is an integer
% multiple of the sampling period, T_s by using the modulo operation.
N = TL/T_s;

if mod(N,1) == 0
    T_J = TJ;
    T_A = TA;
    T_v = Tv;
    T_D = TD;
    T_L = TL;
else
    [T_J,T_A,T_v,T_D,T_L] = segments(TJ,TA,Tv,TD,T_s);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Preallocate memory for tangential jerk in T1
s_tdot1 = zeros(1,round((T_J/T_s),0));
for i = 1:length(s_tdot1)
    s_tdot1(i) = ( -1*6*A/(T_J^3) )*( (i*T_s)^2 ) + ( 6*A/(T_J^2) )*(i*T_s);
end

% Preallocate memory for tangential acceleration in T2
s_tdot2 = zeros(1,round((T_A/T_s),0));

% The tangential jerk in T3 is the negative of T1
s_tdot3 = -1*s_tdot1;

% Preallocate memory for tangential velocity in T4
s_tdot4 = zeros(1,round((T_v/T_s),0));

% The tangential jerk in T5 is equal the tangential jerk in T3
s_tdot5 = s_tdot3;

% Preallocate memory for tangential decceleration in T6
s_tdot6 = zeros(1,round((T_D/T_s),0));

% The tangential jerk in T7 is equal the tangential jerk in T1
s_tdot7 = s_tdot1;

% % % The tangential jerk for T_L is obtained by concatenation.
if T_v > 0
    s_tdot = [s_tdot1 s_tdot2 s_tdot3 s_tdot4 s_tdot5 s_tdot6 s_tdot7];
else T_v = 0;
    s_tdot = [s_tdot1 s_tdot2 s_tdot3 s_tdot5 s_tdot6 s_tdot7];
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Preallocate memory for tangential acceleration in T1
s_ddot1 = zeros(1,round((T_J/T_s),0));

for i = 1:length(s_ddot1)
    s_ddot1(i) = ( -1*2*A/(T_J^3) )*( (i*T_s)^3 ) + ( 3*A/(T_J^2) )*((i*T_s)^2);
end

% Constant acceleration in T2 and T6.
s_ddot2 = A*ones(1,round((T_A/T_s),0));
s_ddot6 = -A*ones(1,round((T_A/T_s),0));

% The tangential acceleration in T3 is the negative of T1 + A
s_ddot3 = -1*s_ddot1 + A;

% The tangential acceleration in T5 is equal to the tangential acceleration
% in T3.
s_ddot5 = -1*s_ddot1;

% The tangential acceleration in T7 is equal to the tangential acceleration
% in T1.
s_ddot7 = s_ddot1 - A;

% % % % The tangential acceleration is obtained by concatenation.
if T_v > 0
    s_ddot = [s_ddot1 s_ddot2 s_ddot3 s_tdot4 s_ddot5 s_ddot6 s_ddot7];
else T_v = 0;
    s_ddot = [s_ddot1 s_ddot2 s_ddot3 s_ddot5 s_ddot6 s_ddot7];
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Preallocate memory for feedrate (tangential velocity) in T1.
s_dot1 = zeros(1,round((T_J/T_s),0));

for i = 1:length(s_dot1)
    s_dot1(i) = v_0 + ( -1*A/(2*(T_J^3)) )*( (i*T_s)^4 ) + ( A/(T_J^2) )*((i*T_s)^3);
end

% Preallocate memory for feedrate (tangential velocity) in T2.
s_dot2 = zeros(1,round((T_A/T_s),0));

% Feedrate in T2.
for i = 1:length(s_dot2)
    s_dot2(i) = s_dot1(end) + A*(i*T_s);
end

% The feed rate in T3 is the s_dot2(end) + negative of the feedrate in T1
s_dot3 = zeros(1,round((T_J/T_s),0));
for i = 1:length(s_dot3)
    s_dot3(i) = s_dot2(end) + -1*s_dot1(i) + A*(i*T_s);
end

% The feedrate in T4 is constant.
if T_v > 0
    s_dot4 = s_dot3(end)*ones(1,round((T_v/T_s),0));
else T_v = 0;
    s_dot4 = s_dot3(end);
end

% The feedrate in T5 is the s_dot4(end) + negative of the feedrate in T3.
for i = 1:length(s_dot1)
    s_dot5(i) = s_dot4(end) - s_dot1(i);
end

% The feedrate in T6 is the s_dot5(end) + negative of the feedrate in T2.
s_dot6 = zeros(1,round((T_A/T_s),0));
for i = 1:length(s_dot6)
    s_dot6(i) = s_dot5(end) - A*(i*T_s);
end

% The feed rate in T7 is the s_dot6(end) + the feedrate in T1 tranlated by
% A*T_s.
s_dot7 = zeros(1,round((T_J/T_s),0));
for i = 1:length(s_dot1)
    s_dot7(i) = s_dot6(end) + s_dot1(i) - A*(i*T_s);
end

% % % % The feedrate is obtained by concatenation.
if T_v > 0
    s_dot = [s_dot1 s_dot2 s_dot3 s_dot4 s_dot5 s_dot6 s_dot7];
else T_v = 0;
    s_dot = [s_dot1 s_dot2 s_dot3 s_dot5 s_dot6 s_dot7];
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Preallocate memory distance travelled in T1.
s1 = zeros(1,round((T_J/T_s),0));
for i = 1:length(s1)
    s1(i) = ( -1*A/(10*(T_J^3)) )*((i*T_s)^5) + ( A/(4*(T_J^2)) )*((i*T_s)^4) + v_0*(i*T_s);
end

% Distance travelled in T2.
s2 = zeros(1,round((T_A/T_s),0));
for i = 1:(T_A/T_s)
    s2(i) = s1(end) + (A/2)*((i*T_s)^2) + s_dot1(end)*(i*T_s);
end

% The distance travelled in T3 is the distance at the end of T2 + the
% distance travelled in T3.
s3 = zeros(1,round((T_J/T_s),0));
for i = 1:length(s3)
     s3(i) = s2(end) + s_dot2(end)*(i*T_s) + -1*( s1(i) ) + (A/2)*(i*T_s)^2;
end

% The distance travelled in T4 is the distance at the end of T3 + the
% distance travelled in T4.
if T_v > 0
    s4 = zeros(1,round((T_v/T_s),0));
    for i = 1:length(s4)
        s4(i) = s3(end) + v_f*(i*T_s);
    end
else T_v = 0;
    s4 = 0;
end

% The distance travelled in T5 is the distance at the end of T4 + the
% distance travelled in T5.
s5 = zeros(1,round((T_J/T_s),0));
if T_v > 0
    for i = 1:length(s5)
        s5(i) = s4(end) + -1*s1(i) + s_dot4(end)*(i*T_s);
    end
else T_v = 0;
    for i = 1:length(s5)
        s5(i) = s3(end) + -1*s1(i) + s_dot3(end)*(i*T_s);
    end
end

% The distance travelled in T6 is the distance at the end of T5 + the
% distance travelled in T6.
s6 = zeros(1,round((T_A/T_s),0));
for i = 1:length(s6)
    s6(i) = s5(end) - (A/2)*(i*T_s)^2 + s_dot5(end)*(i*T_s);
end

% The distance travelled in T7 is the distance at the end of T6 + the
% distance travelled in T7.
s7 = zeros(1,round((T_J/T_s),0));
for i = 1:length(s7)
    s7(i) = s6(end) + s1(i) + s_dot6(end)*(i*T_s) - (A/2)*(i*T_s)^2;
end

% The distance travelled is obtained by concatenation.
if T_v > 0
    s = [s1 s2 s3 s4 s5 s6 s7];
else T_v = 0;
    s = [s1 s2 s3 s5 s6 s7];
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Generates a positive displacement vector if the change in distance is
% positive and generates a negative displacement vector if the change is 
% distance is negative.
if d > 0
    s = s;
else d < 0;
    s = -1*s;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% % Generation of vectors for the dwell time.
% sdwell = s7(end)*ones(1,(T_dwell/T_s));
% sdwell_dot = s_dot7(end)*ones(1,(T_dwell/T_s));
% sdwell_ddot = s_ddot7(end)*ones(1,(T_dwell/T_s));
% sdwell_tdot = s_tdot7(end)*ones(1,(T_dwell/T_s));
% 
% % Generation of vectors for return displacement.
% sreturn = -1*s + s7(end);
% sreturn_dot = -1*s_dot + s_dot7(end);
% sreturn_ddot = -1*s_ddot + s_ddot7(end);
% sreturn_tdot = -1*s_tdot + s_tdot7(end);
% 
% % Generation of vectors for forward and return displacement.
% scomplete = [s sdwell sreturn];
% scomplete_dot = [s_dot sdwell_dot sreturn_dot];
% scomplete_ddot = [s_ddot sdwell_ddot sreturn_ddot];
% scomplete_tdot = [s_tdot sdwell_tdot sreturn_tdot];
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Verification by numerical differentiation for s.

% The distance travelled is given by s, therefore the feedrate is
% obtained by the differene between proceeding value and the preceeding
% value at each value divided by twice the sampling period.
ds = zeros(1,length(s)-2);
for i = 1:length(s)-2
    ds(i) = ( s(i+2) - s(i) )/(2*T_s);
end

dds = zeros(1,length(s)-2);
for i = 1:length(s)-2
    dds(i) = ( s_dot(i+2) - s_dot(i) )/(2*T_s);
end

ddds = zeros(1,length(s)-2);
for i = 1:length(s)-2
    ddds(i) = ( s_ddot(i+2) - s_ddot(i) )/(2*T_s);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% T = T_s:T_s:T_s*length(scomplete);
% 
% figure(1)
% set(gcf,'color','white')
% 
% sp(1) = subplot(4,1,1);
% plot(T,scomplete,'k-')
% title('Feed Profile for a Jerk-continous Trajectory')
% xlabel('Time (ms)')
% ylabel('Displacement (mm)')
% 
% sp(2) = subplot(4,1,2);
% plot(T(2:end-1),scomplete_dot(2:end-1),'k-',T(2:end-1),ds,'--m')
% xlabel('Time (ms)')
% ylabel('Feedrate (mm/s)')
% 
% sp(3) = subplot(4,1,3);
% plot(T(2:end-1),scomplete_ddot(2:end-1),'k-',T(2:end-1),dds,'--m')
% xlabel('Time (ms)')
% ylabel({'Tangential','Acceleration','(mm/s^2)'})
% 
% sp(4) = subplot(4,1,4);
% plot(T(2:end-1),scomplete_tdot(2:end-1),'k-',T(2:end-1),ddds,'--m')
% xlabel('Time (ms)')
% ylabel({'Tangential','Jerk','(mm/s^3)'})
% legend('Analytical Solution','Numerical Derivative Solution','location','Southoutside')
% 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

T = T_s:T_s:T_s*length(s);

figure(1)
set(gcf,'color','white')

sp(1) = subplot(4,1,1);
plot(T,s,'k-')
% title('Feed Profile for the Total Arc Length of a Fan-Shaped Toolpath')
xlabel('Time (ms)','fontsize',15)
ylabel({'Displacement','(mm)'},'fontsize',15)
% legend('Analytical Solution','location','Northwest')

sp(2) = subplot(4,1,2);
plot(T(2:end-1),s_dot(2:end-1),'-k',T(2:end-1),ds,'--m')
xlabel('Time (ms)','fontsize',15)
ylabel({'Feedrate','(mm/s)'},'fontsize',15)

sp(3) = subplot(4,1,3);
plot(T(2:end-1),s_ddot(2:end-1),'k-',T(2:end-1),dds,'--m')
xlabel('Time (ms)','fontsize',15)
ylabel({'Tangential','Acceleration','(mm/s^2)'},'fontsize',15)

sp(4) = subplot(4,1,4);
plot(T(2:end-1),s_tdot(2:end-1),'k-',T(2:end-1),ddds,'--m')
xlabel('Time (ms)','fontsize',15)
ylabel({'Tangential','Jerk','(mm/s^3)'},'fontsize',15)
legend('Analytical Solution','Numerical Derivative Solution','location','Southoutside')


end

