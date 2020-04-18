close all
clear all
clc

maxIter = 50;
samplesize = 20
samples_xstore = zeros(20,maxIter+1);
samples_xhatstore = zeros(20,maxIter);
samples_zstore = zeros(20,maxIter+1);
samples_SigmaXstore = zeros(20,maxIter);
%generating samples for the state

SigmaW=1
xtrue= 1+randn(1)
xstore = zeros(maxIter+1,length(xtrue));
xstore(1,:)=xtrue;
for states = 1:maxIter
    w = chol(SigmaW)'*randn(1);
    xtrue = -(1/3)*(xtrue^3-1) -(1/3)* w;
    xstore(states+1,:)=xtrue
end




for samples=1:samplesize
% Initialize simulation variables
SigmaV = 1; % Sensor noise covariance
ztrue = 1 +randn(1)
xhat = 1; % Initialize Kalman filter initial estimate
SigmaX = 1; % Initialize Kalman filter covariance
u = 0; % Unknown initial driving input: assume zero

% Reserve storage for variables we might want to plot/evaluate
xhatstore = zeros(maxIter,length(xhat));
zstore = zeros(maxIter+1,length(xtrue)); zstore(1,:)=ztrue;
SigmaXstore = zeros(maxIter,length(xhat)^2);

for k = 1:maxIter
% EKF Step 0: Compute Ahat, Bhat
% Note: For this example, x(k+1) = x_k^3/(1+x_k^2) + w(k)
Ahat = -xhat^2; Bhat=1;
% EKF Step 1a: State estimate time update
% Note: You need to insert your system's f(...) equation here
xhat =  -(1/3)*(xhat^3-1)

% EKF Step 1b: Error covariance time update
SigmaX = Ahat*SigmaX*Ahat' + Bhat*SigmaW*Bhat';
% [Implied operation of system in background, with
% input signal u, and output signal z]
w = chol(SigmaW)'*randn(1);
v = chol(SigmaV)'*randn(1);

ztrue = xstore(k+1,:)

% xtrue = xtrue^3/(1+xtrue^2) + w; % future x is based on present u

% EKF Step 1c: Estimate system output
% Note: You need to insert your system's h(...) equation here
Chat = 1; Dhat = 1;
zhat = xhat;
% EKF Step 2a: Compute Kalman gain matrix
L = SigmaX*Chat'/(Chat*SigmaX*Chat' + Dhat*SigmaV*Dhat');
% EKF Step 2b: State estimate measurement update
xhat = xhat + L*(ztrue - zhat)

SigmaX = SigmaX - L*Chat*SigmaX;
% [Store information for evaluation/plotting purposes]
xhatstore(k,:)=xhat;
zstore(k+1,:) = ztrue; 
SigmaXstore(k,:) = SigmaX(:);
end
samples_xstore(samples,:)= transpose(xstore);
samples_xhatstore(samples,:) = transpose(xhatstore);
samples_zstore(samples,:) = transpose(zstore);
samples_SigmaXstore(samples,:) = transpose(SigmaXstore);

end 

average_SigmaXstore = mean(samples_SigmaXstore,1);
average_xhatstore = mean(samples_xhatstore,1);
average_zstore = mean(samples_zstore,1);

figure(1); clf;
plot(0:maxIter-1,xstore(1:maxIter),'k-','LineWidth',2)
hold on
plot(0:maxIter-1,average_zstore(1:maxIter),'r-',0:maxIter-1,average_xhatstore,'b--', ...
0:maxIter-1,average_xhatstore+2*sqrt(average_SigmaXstore),'m-.',...
0:maxIter-1,average_xhatstore-2*sqrt(average_SigmaXstore),'m-.'); grid;
legend('true','observation','estimate','bounds'); xlabel('Iteration'); ylabel('State');
title('Extended Kalman filter in action');
figure(2); clf;
plot(0:maxIter-1,xstore(1:maxIter)-average_xhatstore,'b-',0:maxIter-1, ...
2*sqrt(average_SigmaXstore),'m--',0:maxIter-1,-2*sqrt(average_SigmaXstore),'m--');
grid; 
legend('Error','bounds');
title('EKF Error with bounds');
xlabel('Iteration'); ylabel('Estimation Error');