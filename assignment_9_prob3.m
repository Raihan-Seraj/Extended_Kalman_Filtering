close all
clear all
clc
% Initialize simulation variables
SigmaW = [0 0; 0 0.05];% Process noise covariance
muW = [0 0]
SigmaV = 1; % Sensor noise covariance
delta_t = 0.05
g  = 9.81

maxIter = 1000;

xtrue =[0;0]+ transpose(mvnrnd(muW,SigmaW,1)) % Initialize true system initial state
ztrue = 0 +randn(1)
xhat = [0;0]; % Initialize Kalman filter initial estimate
SigmaX = [1 1; 1 1]; % Initialize Kalman filter covariance
u = 0; % Unknown initial driving input: assume zero

% Reserve storage for variables we might want to plot/evaluate
xstore = zeros(maxIter+1,1); xstore(1,:) = xtrue(1);
xhatstore = zeros(maxIter,1);
zstore = zeros(maxIter+1,length(xtrue)); zstore(1,:)=ztrue;
SigmaXstore = zeros(maxIter,length(xhat)^2);

for k = 1:maxIter
% EKF Step 0: Compute Ahat, Bhat
% Note: For this example, x(k+1) = x_k^3/(1+x_k^2) + w(k)
Ahat = [1 delta_t;-g*cos(xhat(1))*delta_t 1]; Bhat=[1 1];
% EKF Step 1a: State estimate time update
% Note: You need to insert your system's f(...) equation here
xhat = [xhat(1)+xhat(2)*delta_t; xhat(2)-g*sin(xhat(1))*delta_t];
% EKF Step 1b: Error covariance time update
SigmaX = Ahat*SigmaX*Ahat' + Bhat*SigmaW*Bhat';
% [Implied operation of system in background, with
% input signal u, and output signal z]
w = transpose(mvnrnd(muW,SigmaW,1));
v = chol(SigmaV)'*randn(1);

ztrue = sin(xtrue(1)) + v; % z is based on present x and u
xtrue = [xtrue(1)+xtrue(2)*delta_t; xtrue(2)-g*sin(xtrue(1))*delta_t] + w; % future x is based on present u
xtrue = sin(xtrue)
xtrue = asin(xtrue)
% EKF Step 1c: Estimate system output
% Note: You need to insert your system's h(...) equation here
Chat = [cos(xhat(1)) 0];Dhat=[1 1]
zhat = sin(xhat(1));
% EKF Step 2a: Compute Kalman gain matrix
L = SigmaX*Chat'/(Chat*SigmaX*Chat' + Dhat*SigmaV*Dhat');
% EKF Step 2b: State estimate measurement update
xhat = xhat + L*(ztrue - zhat);
xhat = sin(xhat)
xhat  = asin(xhat)
SigmaX = SigmaX - L*Chat*SigmaX;
% [Store information for evaluation/plotting purposes]
xstore(k+1,:) = xtrue(1); xhatstore(k,:) = xhat(1);
zstore(k+1,:) = ztrue; 
SigmaXstore(k,:) = SigmaX(:);
end
SigmaXstore = SigmaXstore(:,4)
% SigmaXstore = sin(SigmaXstore)
% SigmaXstore = asin(SigmaXstore)
figure(1); clf;
% xaxis = 0:delta_t:10
xaxis = 0:delta_t:50
plot(xaxis(1:maxIter),xstore(1:maxIter),'k-', 'LineWidth',2)
hold on
scatter(xaxis(1:maxIter),zstore(1:maxIter),'filled','r')
hold on
plot(xaxis(1:maxIter),xhatstore,'b-','LineWidth',2)
hold on 
plot(xaxis(1:maxIter),xhatstore+2*sqrt(SigmaXstore),'m-.',...
xaxis(1:maxIter),xhatstore-2*sqrt(SigmaXstore),'m-.'); grid;
legend('true','observation','estimate','bounds'); xlabel('Time'); ylabel('State');
title('Extended Kalman filter in action');
figure(2); clf;
plot(xaxis(1:maxIter),xstore(1:maxIter)-xhatstore,'b-',0:maxIter-1, ...
2*sqrt(SigmaXstore),'m--',0:maxIter-1,-2*sqrt(SigmaXstore),'m--');
grid; 
legend('Error','bounds');
title('EKF Error with bounds');
xlabel('Time'); ylabel('Estimation Error');