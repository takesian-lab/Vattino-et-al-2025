function [Rin, Rs, Cm, error] = calcRs_CC(data,dx,stimTable)
% function for calculating Rs for current clamp experiments
% Input:
%   1. data: single data trace
%   2. acq+freq: acquisition frequence (dx),
%   3. stimTable: RC stim settings for a single input -- could be already loaded in the prot folder with data or read from the raw data trace's user data (need to be loaded before)  a
% Output:
%   1. Rin
%   2. Rs
%   3. Cm
%   4. error

% turn off warnings
w1 = 'curvefit:fit:noStartPoint';
warning('off',w1);

w2 = 'MATLAB:colon:nonIntegerIndex';
warning('off',w2);


%% define stim variables from the stim table for this acq
pulseAmp = stimTable.amplitude;
pulseStart = stimTable.delay;
pulseStop = stimTable.delay + stimTable.pulseWidth;
baselineStart = 1;
baselineEnd = stimTable.delay-1;
%%
Rin=NaN;
Rs=NaN;
Cm=NaN;
error=NaN;

% Time variables in sampling points:
databaselineStart = baselineStart/dx;
databaselineEnd = baselineEnd/dx;
datapulseStart = pulseStart/dx;
datapulseStart = datapulseStart+1;
datapulseStop = pulseStop/dx;


%% Confirm that stimTable is correct (we do have recorded RC data)
% sometimes (for unknown reason, even when checking the CC_RC box during
% acquisition, and the stim data records RC data being recorded, the data
% clearly do not have CC_RC drop

if std(data(datapulseStart:datapulseStop))< 0.1
    error = 1;
else % proceed as normal
    %% calculate Rs, Rm, Cm and tau
    I = abs(pulseAmp); %in pA
    I = I*1e-12; %in A
    
    %V drop in linear phase:
    v1 = abs(data(round(datapulseStart+datapulseStart*0.002))-data(datapulseStart)); %in mV
    v1 = v1*1e-3; %in V
    
    Rs = v1/I; %in V/A = Ohm
    
    
    %V drop in steady state:
    v2 = abs(mean(data(round(datapulseStop-datapulseStop*0.005:datapulseStop)))-mean(data(databaselineStart:databaselineEnd)));
    v2 = v2*1e-3; %in V
    
    Rm = v2/I-Rs; %in V/A = Ohm
    
    V = data(datapulseStart:datapulseStop-datapulseStop*0.005);
    x = 1:length(V);
    
    % now try to fit data to exponential function and calculate tau
    fo = fitoptions('Method', 'NonlinearLeastSquares');
    % Fit exponential decay curve to data
    ft = fittype('m*exp(-x/tau)+b');
    % Try to fit the curve and catch any errors
    done = false; % flag to indicate if fit is successful
    counter = 1;
      
    
    while ~done
        try
            [f, gof] = fit(x', V, ft, fo);
            counter = counter;
            done = true; % fit was successful, set flag to true
        catch
            fo = fitoptions('Method', 'NonlinearLeastSquares');
            counter = counter + 1;
        end
    end
    % Plot the fit
    %     figure
    %     plot(f, x, V)
    %     title('Exponential Decay Fit')
    %     xlabel('Time(t)')
    %     ylabel('Voltage(V)')
    %     legend('Data', 'Exponential Decay Fit')
    %
    
    tau_ms = f.tau*dx; %in ms
    tau_sec = tau_ms/1e+3;
    Cm = tau_sec/ Rm; %s/ohm = farad
    
    Cm = Cm / 1e-12; %in pF
    Rs = Rs/1e+6; %in MOhm
    Rm = Rm/1e+6; %in MOhm
    
    if any(~isreal([tau_ms Rin, Rs, Cm]))
        disp('calcRs: complex number returned.  No fit possible')
        Rin= [];
        Rs= [];
        Cm= [];
        error=1;
        return
    end
    
end
end
