% % ELEC 4700 Assignment 1: Monte-Carlo Modelling of Electrons
% % Mike Crawford
% % 100952432
% % Question1
%%Would like the one on one evaluation please
% This code is designed to simulate a specified number of electrons movement 
% in silicon using the Monte-Carlo method. The effective mass is known and 
% each electron is given a random direction. Boundaries were put on the box 
% to allow electrons to pass through the other side of the y plane and 
% reflect off the x plane. The thermal velocity is calculated by the equation
% 
% Vth = sqrt((C.k*T)/(mn)
% 
% The parameters are labeled below. The mean free path is can be found with
% the thermal velocity and the mean time between collisions, which is given
% as 0.2 ps. The simulation plots these electrons during the specified 
% timeframe as well as the average semiconductor temperature over time.

%reset
clearvars
clearvars -GLOBAL
close all
format shortE
global C
global Ecount
global Vx Vy Vtotal x y
Ecount =5000;   %setting amount of electrons used
C.mo = 9.10938215e-31;  %Electron mass
C.k = 1.3806504e-23;    %Boltzmann constant

T =300;     %Temperature (K)
mn = 0.26*C.mo;
L = 200e-9; %Length of Border
W = 100e-9; %Width of Border
Vth = sqrt(2*(C.k*T)/mn); %Calculation of Thermal Velocity
dt = 10e-15; %timestep
Stop =500*dt; %timeframe
x = zeros(Ecount, 2);   %array of x positions
y = zeros(Ecount, 2);   %array of y positions

Temperature = zeros(1,2);   %array of temperatures
Time = 0;
VisibleEcount = 10; %number of visible electrons

for i = 1:Ecount %setting random locations for eact electron
    x(i,1) = rand()*200e-9; 
    y(i,1) = rand()*100e-9;
end

Vx(1:Ecount) = Vth * cos(2*pi*randn(1, Ecount)); %setting velocities
Vy(1:Ecount) = Vth * sin(2*pi*randn(1, Ecount));
for i = 1:Ecount
    Vtotal(i) = sqrt(Vx(i)^2 + Vy(i)^2); %summing x and y velocities
end


figure(1)   %plotting electrons
subplot(2,1,1);
axis([0 L 0 W]);
title('Monte-Carlo model of electron movement in Silicon');
xlabel('X');
ylabel('Y');
hold on;

subplot(2,1,2); %plotting avg temperature
axis([0 Stop 0 500]);
title('Semiconductor Temperature');
xlabel('Time (s)');
ylabel('Temperature (K)');
hold on;

for i = 1:Ecount %Sum of all initial temperatures
   Temperature(1,2) = Temperature(1,2) + (mn*Vtotal(i)^2)/(2*C.k);
end
%intitial average temperature
AvgTemperature = Temperature(1,2)/Ecount;
TemperaturePlot = [300 AvgTemperature]; %setting temperature line to plot
TimePlot = [0 Time];    %setting time line to plot
plot(Time, AvgTemperature);

%reseting values before loop
VTotal = 0;
Temperature(1,2) = 0;
AvgTemperature = 0;

%looping through timesteps
while Time < Stop
    subplot(2,1,1)
    for j = 1:Ecount %looping through electrons
        x(j,2) = x(j,1); %updating previous positions
        y(j,2) = y(j,1);
        x(j,1) = x(j,1) + (dt * Vx(j));%updating new positions
        y(j,1) = y(j,1) + (dt * Vy(j));
        if x(j,1) > L %setting right side boundary condition
            x(j,2) = 0;
            x(j,1) = dt * Vx(j);
        end
        if x(j,1) < 0 %setting left side boundary condition
            x(j,2) = L;
            x(j,1) = x(j,2) + (dt * Vx(j));
        end 
        %setting roof/floor boundary condition
        if y(j,1) > W || y(j,1) < 0 
            Vy(j) = -Vy(j);
        end
        %setting up position line vectors to plot
        XPlot = [x(j,2) x(j,1)];
        YPlot = [y(j,2) y(j,1)];
        if j < VisibleEcount
            %plot visible position line vectors
        plot(XPlot,YPlot);
        end
        
        %sum of x and y  velocities
       VTotal = sqrt(Vx(j)^2 + Vy(j)^2);
       
        %sum of all temperatures
       Temperature(1,2) = Temperature(1,2) + (mn*Vtotal(j)^2)/(2*C.k);
         
    end
    %averaging temperature sum and plotting against each timestep
    AvgTemperature = Temperature(1,2)/Ecount;
    TemperaturePlot = [Temperature(1,1) AvgTemperature];
    TimePlot = [(Time - dt) Time];
    subplot(2,1,2);
    plot(TimePlot, TemperaturePlot);
    Temperature(1,1) = AvgTemperature;
    AvgTemperature = 0;
    Temperature(1,2) = 0;
    pause(1e-19)
    Time = Time + dt;
end