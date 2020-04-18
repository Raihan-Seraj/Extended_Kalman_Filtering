close all
clear all
clc
% Initialize simulation variables
SigmaW = 1; % Process noise covariance
SigmaV = 1; % Sensor noise covariance
maxIter = 50;

xtrue =1+ randn(1); % Initialize true system initial state
ztrue = 1 +randn(1)
xhat = 1; % Initialize Kalman filter initial estimate
SigmaX = 1; % Initialize Kalman filter covariance
u = 0; % Unknown initial driving input: assume zero

% Reserve storage for variables we might want to plot/evaluate
xstore = zeros(maxIter+1,length(xtrue)); xstore(1,:) = xtrue;
xhatstore = zeros(maxIter,length(xhat));
zstore = zeros(maxIter+1,length(xtrue)); zstore(1,:)=ztrue;
SigmaXstore = zeros(maxIter,length(xhat)^2);
SigmaXstoreplus = zeros(maxIter, length(xhat));

for k = 1:maxIter,
% EKF Step 0: Compute Ahat, Bhat
% Note: For this example, x(k+1) = x_k^3/(1+x_k^2) + w(k)
Ahat = (3*xhat^2-xhat^4)/(1+xhat^2)^2; Bhat=1;
% EKF Step 1a: State estimate time update
% Note: You need to insert your system's f(...) equation here
xhat = xhat^3/(1+xhat^2)

% EKF Step 1b: Error covariance time update
SigmaX = Ahat*SigmaX*Ahat' + Bhat*SigmaW*Bhat';
SigmaXstoreplus(k,:)=SigmaX(:)
% [Implied operation of system in background, with
% input signal u, and output signal z]
w = chol(SigmaW)'*randn(1);
v = chol(SigmaV)'*randn(1);

ztrue = xtrue^5 + v; % z is based on present x and u
xtrue = xtrue^3/(1+xtrue^2) + w; % future x is based on present u

% EKF Step 1c: Estimate system output
% Note: You need to insert your system's h(...) equation here
Chat = 5*xhat^4; Dhat = 1;
zhat = xhat^5;
% EKF Step 2a: Compute Kalman gain matrix
L = SigmaX*Chat'/(Chat*SigmaX*Chat' + Dhat*SigmaV*Dhat');
% EKF Step 2b: State estimate measurement update
xhat = xhat + L*(ztrue - zhat);

SigmaX = SigmaX - L*Chat*SigmaX;
% [Store information for evaluation/plotting purposes]
xstore(k+1,:) = xtrue; xhatstore(k,:) = xhat;
zstore(k+1,:) = ztrue; 
SigmaXstore(k,:) = SigmaX(:);
end

figure(1); clf;
plot(0:maxIter-1,xstore(1:maxIter),'k-',0:maxIter-1,zstore(1:maxIter),'r-',0:maxIter-1,xhatstore,'b--', ...
0:maxIter-1,xhatstore+2*sqrt(SigmaXstore),'m-.',...
0:maxIter-1,xhatstore-2*sqrt(SigmaXstore),'m-.'); grid;
legend('true','observation','estimate','bounds'); xlabel('Iteration'); ylabel('State');
title('Extended Kalman filter in action');
figure(2); clf;
plot(0:maxIter-1,xstore(1:maxIter)-xhatstore,'b-',0:maxIter-1, ...
2*sqrt(SigmaXstore),'m--',0:maxIter-1,-2*sqrt(SigmaXstore),'m--');
grid; 
legend('Error','bounds');
title('EKF Error with bounds');
xlabel('Iteration'); ylabel('Estimation Error');
figure(3);clf
plot(0:maxIter-1,SigmaXstore,'r','LineWidth',2)
hold on 
plot(0:maxIter-1,SigmaXstoreplus,'b','LineWidth',2)
legend('V_{k|k}','V_{k+1|k}')
xlabel('Iteration'); ylabel('Ricatti Euqation Trajectories');