% written by Shahrokh Akhlaghi
%    Date: 01/05/2015
clc
clear all
% close all


%% define our meta-variables (i.e. how long and often we will sample)
duration = 20;       % how long the sample
dt       = 0.05;    % The observation continuously looks for data,

%% Define update equations (Coefficent matrices):
A  = 0.5;      %[1 dt; 0 1] ; % state transition matrix
C1 = 0.5;


Cseq       = 0 : 0.5 : 10;
for cIndex = 1:length(Cseq)
    C      = Cseq(cIndex);

    randn('state', 40);
        
    %% define main variables
    u                 = 1;                       % define acceleration magnitude
    X                 = 0.01;                       % initized state
    X_estimate        = X;                       % x_estimate of initial location estimation
    State_noise_mag   = 0.1 ;%* (randn);         % w: process noise
    Measurement_Noise = 0.010 ;%* (randn);         % v: measurement noise
    
    EQ                = State_noise_mag;       % convert the   process   noise (stdv) into covariance matrix
    ER                = Measurement_Noise;     % convert the measurement noise (stdv) into covariance matrix
    P                 = EQ;                      % (covariance matrix)
    
    %% initize result variables
    % Initialize for speed
    X_loc      = [];
    X_loc_meas = [];
    
    
    %% simulate over time
    A0 = 0 : dt: duration;
    
    for t = 0 : dt: duration
        
        State_noise       = State_noise_mag  * (randn);     
        X                 = A * X + State_noise; %+ 0.01*(randn); % Q = A*x
        
        NinjaVision_noise = Measurement_Noise;
        y                 = C1 * X + C * X.^3 + NinjaVision_noise; %+ (randn);
        X_loc             = [X_loc; X];
        X_loc_meas        = [X_loc_meas; y];
        
    end
   
    if 1
        figure(11); clf;
        hold on
        plot(0:dt:t, X_loc, 'ro-', 'linewidth', 1);
        plot(0:dt:t, X_loc_meas, '-k.');
        
        legend('State', 'Meas');
        xlabel('Time (s)');
        ylabel('Amp');
    end
    
    
    %% ======================= ** kalman filtering ** =========================
    %initize estimation variables
    
    X_loc_estimate = [];
    X_loc_preEst   = [];
    vel_estimate   = [];
    X              = 0;  
    P_estimate     = P;
    P_mag_estimate = [];
    predic_state   = [];
    predic_var     = [];
    
    
    for t = 1:length(X_loc)
          P_mag_estimate = [P_mag_estimate; P];
% %         
        
        fstate                  = @(XXX) A * XXX;
        hmeas                   = @(XXX) C1* XXX + C * XXX.^3;
        X_loc_preEst            = [X_loc_preEst; X_estimate];
        [X_estimate,predic_var] = ekf(fstate,X_estimate,P,hmeas,X_loc_meas(t),EQ,ER);
        X_loc_estimate          = [X_loc_estimate; X_estimate];
        P_mag_estimate          = [P_mag_estimate; predic_var];
        
       
    end
    
    MyColor  = {'-*b','-+r','-om','-*c','-+k','-og','-*y'};
    %% Plot the results
    figure(2);clf;
    tt = 0 : dt : duration;
    plot(tt,X_loc,'-r.','linewidth',2); hold on;
    plot(tt,X_loc_estimate,'-bo');
    legend('True States', 'Estimated');
    xlabel('Time (s)'); ylabel('Amp');
    ErrAbs = abs(X_loc_estimate-X_loc);
    ErrRel = sum(ErrAbs)/sum(abs(X_loc));
    sTmp   = sprintf('Non-linear Index C = %4.3f (Relative Error = %4.1f%%)', C, ErrRel*100);
    title(sTmp);

    
    if 0
        figure(5);
        plot(tt,delta_x,'-r.');
        ylabel('delta\_x')
        xlabel('Time (s)');
        title('state perturbation')
        
        figure(6);
        subplot(2,1,1)
        plot(tt,Epsilon_S,'-r.');
        ylabel('Epsilon State')
        subplot(2,1,2)
        plot(tt,Epsilon_M,'-r.');
        ylabel('Epsilon Measurement')
    
    
    figure(3); clf;
    subplot(2,1,1)
    plot(tt,abs(Quasi_verify_S),'-r.',tt,.1,'--');
    legend('State Non-linearity','Threshold');
    ylabel('State Non-linearity')
%     axis([0.0 7 -.1 .5])
    subplot(2,1,2)
    plot(tt,abs(Quasi_verify_M),'-r.',tt,.1,'--');
    legend('Meas Non-linearity','Threshold');
    ylabel('Meas Non-linearity')
    
%     axis([0.0 7 -.1 .5])
    xlabel('time (s)');
    end
%     pause
end



