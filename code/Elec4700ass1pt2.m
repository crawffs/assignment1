% % ELEC 4700 Assignment 1: Monte-Carlo Modelling of Electrons
% % Mike Crawford
% % 100952432
% % Question2
% This next part is similar to the previous, with the same model of electrons
% and average tmeperatures. The difference is it adds a collission factor 
% with a scattering probability as welll as a random velocity with each
% electron. A histogram showing the distribution of diferent velocities is
% also given. This part is using the same code skelaton as the first part
%so much of the code will look the same

%reset
clearvars
clearvars -GLOBAL
close all
format shortE
global C
global Ecount
global Vx Vy Vtotal x y
Ecount =1000;   %setting amount of electrons used
C.mo = 9.10938215e-31;%Electron mass
C.k = 1.3806504e-23; %Boltzmann constant

T =300; %Temperature (K)
mn = 0.26*C.mo;
L = 200e-9;  %Length of Border
W = 100e-9; %Width of Border
Vth = sqrt((2*C.k*T)/mn); %Calculation of Thermal Velocity
dt = 10e-15; %timestep
Stop = 100*dt; %timeframe
x = zeros(Ecount, 2); %array of x positions
y = zeros(Ecount, 2); %array of y positions

Temperature = zeros(1,2); %array of temperatures
Time = 0;
VisibleEcount = 50; %visibble electrons
tmn = 0.2e-12; %mean time between collisions
PScat = 1 - exp(-dt/tmn); %scattering probability equation
VTotalHistogram = zeros(Ecount, 1);%array of thermal velocities
for i = 1:Ecount
    x(i,1) = rand()*200e-9;%setting random electron positions
    y(i,1) = rand()*100e-9;
end
%setting random velocities
Vx(1:Ecount) = Vth*rand * cos(2*pi*randn(1, Ecount));
Vy(1:Ecount) = Vth*rand * sin(2*pi*randn(1, Ecount));
%sum of x and y velocities
for i = 1:Ecount
    Vtotal(i) = sqrt(Vx(i)^2 + Vy(i)^2);
    VTotalHistogram(i) = sqrt(Vx(i)^2 + Vy(i)^2);

end

%plot of electrons
figure(1)
subplot(2,1,1);
axis([0 L 0 W]);
title('Monte-Carlo model of electron movement in Silicon');
xlabel('X');
ylabel('Y');
hold on;
%plot of average temperatures
subplot(2,1,2);
axis([0 Stop 0 500]);
title('Semiconductor Temperature');
xlabel('Time (s)');
ylabel('Temperature (K)');
hold on;
%initial temperature sum
for i = 1:Ecount
   Temperature(1,2) = Temperature(1,2) + (mn*Vtotal(i)^2)/(2*C.k);
end
%initial average temperatue
AvgTemperature = Temperature(1,2)/Ecount;
TemperaturePlot = [300 AvgTemperature];
TimePlot = [0 Time];
plot(Time, AvgTemperature);
%resetting before loop
VTotal = 0;
Temperature(1,2) = 0;
AvgTemperature = 0;
%looping timestep
while Time < Stop
    subplot(2,1,1)
    %looping electrons
    for j = 1:Ecount
        %Scattering probability conditioned with random value to simulate
        %random scattering at  set rate
        if PScat> rand
                Vx(j) = Vth * randn;
                Vy(j) = Vth * randn;
        end
        %update previous/new x and y positions
        x(j,2) = x(j,1);
        y(j,2) = y(j,1);
        x(j,1) = x(j,1) + (dt * Vx(j));
        y(j,1) = y(j,1) + (dt * Vy(j));
        %check right wall border
        if x(j,1) > L
            x(j,2) = 0;
            x(j,1) = dt * Vx(j);
        end
        %check left wall border
        if x(j,1) < 0
            x(j,2) = L;
            x(j,1) = x(j,2) + (dt * Vx(j));
        end
        %check roof/floor border
        if y(j,1) > W || y(j,1) < 0
            Vy(j) = -Vy(j);
        end
        %set line vectores for x and y positions
        XPlot = [x(j,2) x(j,1)];
        YPlot = [y(j,2) y(j,1)];
        %plot visible x and y line vectors
        if j < VisibleEcount
        plot(XPlot,YPlot);
        end
        
        %sum of themal velocities and temperatures
       VTotal = sqrt(Vx(j)^2 + Vy(j)^2);
       Temperature(1,2) = Temperature(1,2) + (mn*Vtotal(j)^2)/(2*C.k);
         
    end
    %plotting average temperature and ressetting values for next timestep
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
%histogram of thermal velocities
figure(2)
histogram(VTotalHistogram)
title('Histogram of Electron Voltage Distribution')
xlabel('Velocity (m/s)')
ylabel('Count')

