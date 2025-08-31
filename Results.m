%% NIKUNJ BHEDA
clc 
clear all 
close all
format long
u = xlsread("E:\19042025\DM with Updated Mirror\Dig_42.5_4\9_3\combined_output.xlsx");
n = size(u);
X = u(:,1);
Y = u(:,2);
% % scatter(X,Y)
% title('2D Surface of mirror')
rad = sqrt(X.^2+Y.^2);
R= max(rad);

%%
for i=1:1:n
    x= X(i);
    y=Y(i);
    r=rad(i);
    theta = atan2(y,x);
    r = r/R;
    a = theta;

%% Computation of Zernike fringe polynomails for each point
    c(i,1)=1;
    c(i,2)= r*cos(a);  %X tilt
    c(i,3)= r*sin(a);  %Y tilt
    c(i,4)=(2*(r^2))-1 ; % defocus
    c(i,5)=(r^2)*cos(2*a); % Astig0
    c(i,6)= (r^2)*sin(2*a);
    c(i,7)=((3*(r^3))-(2*r))*cos(a); %coma
    c(i,8)= ((3*(r^3))-(2*r))*sin(a); %coma
    c(i,9)=((6*(r^4))-(6*(r^2)))+1; %Spherical
    c(i,10)=(r^3)*cos(3*a); %trefoil X
    c(i,11)= (r^3)*sin(3*a); %trefoil Y
    c(i,12) = ((4*(r^4))-(3*(r^2)))*cos(2*a);
    c(i,13) = ((4*(r^4))-(3*(r^2)))*sin(2*a);
    c(i,14) = ((10*(r^5))-(12*(r^3)) + (3*r))*cos(a);
    c(i,15) = ((10*(r^5))-(12*(r^3)) + (3*r))*sin(a);
    c(i,16) = (20*(r^6))-(30*(r^4))+(12*(r^2))-1;
    c(i,17) = (r^4)*cos(4*a);                                                   
    c(i,18) = (r^4)*sin(4*a);
    c(i,19) = ((5*(r^5))-(4*(r^3)))*cos(3*a);
    c(i,20) = ((5*(r^5))-(4*(r^3)))*sin(3*a);
    c(i,21) = ((15*(r^6))-(20*(r^4))+(6*(r^2)))*cos(2*a);
    c(i,22) = ((15*(r^6))-(20*(r^4))+(6*(r^2)))*sin(2*a);
    c(i,23) = ((35*(r^7))-(60*(r^5))+(30*(r^3))-(4*r))*cos(a);
    c(i,24) = ((35*(r^7))-(60*(r^5))+(30*(r^3))-(4*r))*sin(a);
    c(i,25) = ((70*(r^8))-(140*(r^6))+(90*(r^4))-(20*(r^2)))+1;
    c(i,26) = (r^5)*cos(5*a);
    c(i,27) = (r^5)*sin(5*a);
    c(i,28) = ((6*(r^6))-(5*(r^4)))*cos(4*a);
    c(i,29) = ((6*(r^6))-(5*(r^4)))*sin(4*a);
    c(i,30) = ((21*(r^7))-(30*(r^5))+(10*(r^3)))*cos(3*a);
    c(i,31) = ((21*(r^7))-(30*(r^5))+(10*(r^3)))*sin(3*a);
    c(i,32) = ((56*(r^8))-(105*(r^6))+(60*(r^4))-(10*(r^2)))*cos(2*a);
    c(i,33) = ((56*(r^8))-(105*(r^6))+(60*(r^4))-(10*(r^2)))*sin(2*a);
    c(i,34) = ((126*(r^9))-(280*(r^7))+(210*(r^5))-(60*(r^3))+(5*r))*cos(a);
    c(i,35) = ((126*(r^9))-(280*(r^7))+(210*(r^5))-(60*(r^3))+(5*r))*sin(a);
    c(i,36) = ((252*(r^10))-(630*(r^8))+(560*(r^6))-(210*(r^4))+(30*(r^2)))-1;
    c(i,37) = (924*(r^12))-(2772*(r^10))+(3150*(r^8))-(1680*(r^6))+(420*(r^4))-(42*(r^2))+1;

end
%% Initialize base z vector
z_base = zeros(37,1);

% Define your z values as index-value pairs (index, value)
z_values = [
    4, -4.0e-05;
    5, 3.5e-04;
    6, 3.5e-04;
    7, 5.2e-04;
    8, 5.2e-01;
    9, 1.0e-04;
    10, 0.8e-04;
    11, 0.8e-04;
    12, 0.85e-04;
    13, 0.86e-04;
    17, 0.86e-04;
    18, 0.86e-04;
    19, 0.9e-04;
    20, 0.9e-04;
    26, 0.92e-04;
    27, 0.92e-04
];

% Initialize results matrix
results = zeros(size(z_values,1), 3);

% Loop through each z value individually
for k = 1:size(z_values,1)
    
    % Reset z to zeros
    z = z_base;
    
    % Set only current z value
    idx = z_values(k,1);
    val = z_values(k,2);
    z(idx,1) = val;
    
    %% ======== Computation block ========
    for i=1:n
        dz(i,1) = c(i,:)*z;  
    end
    
    % Remove piston, tilt X & Y, focus
    dzf = dz;
    dzfm = mean(dzf);
    
    % Surface RMS
    dzfsq = (dzf - dzfm).^2;
    dzrms = sqrt(sum(dzfsq)./n)*1e+06;  % nm
    
    %% Actuator commands using least fit method
    B = u(:,3:19);
    Alpha = (inv(B'*B))*B'*dzf;
    
    % Compensated wavefront
    beta = B*Alpha;
    betam = mean(beta);
    
    betaSq = (beta - betam).^2;
    Crrms = sqrt(sum(betaSq)./n)*1e+06;
    
    % Residual wavefront
    RWF = dzf - beta;
    RWFm = mean(RWF);
    
    RWFsq = (RWF - RWFm).^2;
    RWFrms = sqrt(sum(RWFsq)./n)*1e+06;  
    
    %% Correction Accuracy
    CA = 100*(1-(RWFrms/dzrms));
    accuracies(k) = round(CA,2);
    % Store results
    results(k,:) = [idx, val, CA];
end

%% Display all results
fprintf('Index\tValue\t\tAccuracy(%%)\n');
for k = 1:size(results,1)
    fprintf('%d\t%.2e\t%.2f\n', results(k,1), results(k,2), results(k,3));
end
%% Save accuracies as a single row in CSV
csv_filename = 'accuracies_row.csv';
writematrix(accuracies, csv_filename);

fprintf('\nOnly accuracy values saved as a single row to %s\n', csv_filename);