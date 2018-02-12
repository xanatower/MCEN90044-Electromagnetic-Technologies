%% Core Loss Calculation Script

% David Meeker
% dmeeker@ieee.org
% 17Oct2017
clear all
clc 
close all
addpath('C:\femm42\mfiles');
savepath;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Start User-Defined Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Model Name
MyModel = 'untitled_part2.fem';

% base frequency in Hz
wbase=50; %(4000 rev/minute)*(minute/(60*seconds))

% range of speeds over which to evaluate losses
SpeedMin = 100; % in RPM
SpeedMax = 8000; % in RPM
SpeedStep = 100; % in RPM

% Winding properties
MyIdCurrent = 0; % direct current in phase current amplitude scaling
MyIqCurrent = 13.9; % quadrature current phase current amplitude scaling
MyLowestHarmonic = 2; % lowest numbered harmonic present in the stator winding
AWG=30;	% Magnet wire gauge used in winding
WindingFill=0.867;%from calculation in Q7 part 1
PhaseResistance = 0.125; % from Q1 part 2
TemperatureRise = 100; % NOT SURE

% Magnet properties
RotorMagnets = 14;  %Number of poles
omag = 0.556*10^6;	% conductivity of sintered NdFeB in S/m

% Core properties ce and ch cant be found from the official data sheet
% therefore assume silicon steel M19 is used
ce = 0.530; % Eddy current coefficient in (Watt/(meter^3 * T^2 * Hz^2)
ch = 143.; % Hysteresis coefficient in (Watts/(meter^3 * T^2 * Hz)
cs = 0.97;  % Lamination stacking factor (nondimensional)
             %http://www.femm.info/wiki/StackingFactor

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% End User-Defined Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% helpful unit definitions
PI=pi; Pi=pi;
deg=Pi/180.;

% angle through which to spin the rotor in degrees
n = 360/MyLowestHarmonic;

% angle increment in degrees
dk = 1;

% A similar loss/volume expression can be derived to compute proximity
% effect losses in the windings.  Use the low frequency approximiation
% from the paper, since in this case, wire size is a lot smaller than skin
% depth at the frequencies of interest.

% Get parameters for proximity effect loss computation for phase windings
dwire=0.324861*0.0254*exp(-0.115942*AWG); % wire diameter in meters as a function of AWG
owire = (58.*10^6)/(1+TemperatureRise*0.004); % conductivity of the wire in S/m at prescribed deltaT
cePhase = (Pi^2/8)*dwire^2*WindingFill*owire;

%% Perform a series of finite element analyses

openfemm;
opendocument(MyModel);
smartmesh(1);
mi_saveas('temp.fem');

% Run an analysis through an entire spin of the rotor and record
% element centroid flux density and vector potential (using the mesh from the first iteration)
% at every step.  This information will then be used to estimate
% core losses
for k = 0:dk:(n-1)
	starttime=clock;
    kk = k/dk + 1;
    
    % make sure that the current is set to the appropriate value for this iteration.  
    tta = (RotorMagnets/2)*k*deg;
	Id = [cos(tta), cos(tta-2*pi/3), cos(tta+2*pi/3)];
    Iq =-[sin(tta), sin(tta-2*pi/3), sin(tta+2*pi/3)];
	Itot =  MyIdCurrent*Id + MyIqCurrent*Iq;
    mi_setcurrent('A', Itot(1));
    mi_setcurrent('B', Itot(2));
    mi_setcurrent('C', Itot(3));
    
    mi_analyze(1);
    mi_loadsolution;
    if (k == 0)
        % Record the initial mesh elements if the first time through the loop
        nn = mo_numelements;
        b = zeros(floor(n/dk),nn); % matrix that will hold the flux density info
		A = zeros(floor(n/dk),nn); % matrix that will hold the vector potential info
        z = zeros(nn,1);
        a = zeros(nn,1);
        g = zeros(nn,1);
        tq = zeros(floor(n/dk),1);
        for m = 1:nn
            elm = mo_getelement(m);
            % z is a vector of complex numbers that represents the location of
            % the centroid of each element.  The real part is x, the
            % imaginary part is y.  The purpose of representing the
            % location in this way is that it is easy to rotate the
            % location of the centroid (by multiplying by a complex number)
            % to find the point that corresponds to the centroid when the
            % rotor is rotated.
            z(m) = elm(4) + 1j*elm(5);
            % element area in the length units used to draw the geometry
            a(m) = elm(6);
            % group number associated with the element
            g(m) = elm(7);
        end
    end
    
    % Store element flux densities *)
    u=exp(1j*k*pi/180.);
    for m = 1:nn
        % Element is on the rotor.  Elements on the rotor have been
        % assigned to groups 10 and higher (by assigning the block label that marks them
        % to group 10, 11, 12, etc.) so that the program can tell if an element is part of the rotor. 
        if(g(m)>=10)
            % the location in the original mesh is rotated so that the
            % flux density at the same point with the rotor rotated through
            % angle k can be evaluated. Multiplying by the complex number u
            % does the rotation.
            p = z(m)*u;
            % Flux densities bx and by are evaluated and rolled into a complex number.
            % Dividing by the complex number u rotates the flux density
			pv = mo_getpointvalues(real(p),imag(p));
            % back into a rotor-fixed reference frame.
			if (g(m)==10)
				% store flux density for elements in rotor core
				b(kk,m) = (pv(2)+1j*pv(3))/u;
			else
				% store vector potential for elements that are in PMs
				A(kk,m) = pv(1);
			end 
        elseif (g(m) > 0) % element is on the stator
            % since the stator doesn't move, no games have to be played
            % with rotations.
            p = z(m);
            b(kk,m) = (mo_getb(real(p),imag(p))*[1;1j]);
        end
    end

    % mo_getprobleminfo returns, among other things, the depth of the
    % machine in the into-the-page direction and the length units used to
    % draw the geometry. Both of these pieces of information will be needed
    % to integrate the losses over the volume of the machine.
    probinfo=mo_getprobleminfo;
    
	% select all blocks on the rotor and compute the torque
	for m=10:(10+RotorMagnets)
        mo_groupselectblock(m);
    end
    tq(kk)=mo_blockintegral(22);  
    mo_close;
    
    % rotate the rotor to the position for the next iteration.
    for m=10:(10+RotorMagnets)
        mi_selectgroup(m);
    end
    mi_moverotate(0, 0, dk); 
	mi_clearselected();
    
    fprintf('% i of % i :: %f seconds ::  %f N*m \n',k,n,etime(clock,starttime),tq(kk));
end

% clean up after finite element runs are finished
closefemm;
delete('temp.fem');
delete('temp.ans');

%% Add Up Core Losses

% Compute the square of the amplitude of each harmonic at the centroid of
% each element in the mesh. Matlab's built-in FFT function makes this easy.
ns=n/dk;
bxfft=abs(fft(real(b)))*(2/ns);
byfft=abs(fft(imag(b)))*(2/ns);
bsq=(bxfft.*bxfft) + (byfft.*byfft);

% Compute the volume of each element in units of meter^3
h = probinfo(3);            % Length of the machine in the into-the-page direction
lengthunits = probinfo(4);  % Length of drawing unit in meters
v = a*h*lengthunits^2;

% compute fft of A at the center of each element
Jm=fft(A)*(2/ns);
for k=1:RotorMagnets
	g3=(g==(10+k));
	% total volume of the magnet under consideration;
	vmag=v'*g3;
	% average current in the magnet for each harmonic
	Jo=(Jm*(v.*g3))/vmag;    
	% subtract averages off of each each element in the magnet
	Jm = Jm - Jo*g3';
end

Iphase=sqrt(MyIdCurrent^2+MyIqCurrent^2)/sqrt(2);
PhaseOhmic = 3*(PhaseResistance*(1+TemperatureRise*0.004))*Iphase^2;
	
results=[];

for thisSpeed=SpeedMin:SpeedStep:SpeedMax

	thisFrequency = thisSpeed/60; % mechanical speed in Hz
	
	% Make a vector representing the frequency associated with each harmonic
	% The last half of the entries are zeroed out so that we don't count each
	% harmonic twice--the upper half of the FFT a mirror of the lower half
	w=0:(ns-1);
	w=MyLowestHarmonic*thisFrequency*w.*(w<(ns/2));  

	% Now, total core loss can be computed in one fell swoop...
	% Dividing the result by cs corrects for the lamination stacking factor
	g1=(g==10);
	rotor_loss = ((ch*w+ce*w.*w)*bsq*(v.*g1))/cs;

	g2=(g==1);
	stator_loss = ((ch*w+ce*w.*w)*bsq*(v.*g2))/cs;

	% and prox losses can be totalled up in a similar way
	g4=(g==2);
	prox_loss = ((cePhase*w.*w)*bsq*(v.*g4));

	% Add up eddy current losses in the magnets
	magnet_loss = (1/2)*((omag*(2*pi*w).^2)*(abs(Jm).^2)*v);

	total_loss = rotor_loss + stator_loss + prox_loss + PhaseOhmic + magnet_loss;

	results = [results; thisSpeed, rotor_loss, stator_loss, magnet_loss, PhaseOhmic, prox_loss, total_loss];
end

% save('c:\\temp\\myLossData.m','results','-ASCII');

%% Loss plot
plot(results(:,1),results(:,2:7));
xlabel('Speed, RPM');
ylabel('Total Losses, Watts');
title('Loss versus Speed');
legend('Rotor Core','Stator Core','Magnets','Coil Ohmic','Coil Proximity','Total Loss','Location','northwest');

% %% Torque plot
% plot(0:dk:(n-1),tq);
% xlabel('Rotor Angle, Degrees');
% ylabel('Torque, N*m');
% title('Torque on Rotor vs. Angle');

% %% Plot of heating
%
% wbase=4000/60; %(4000 rev/minute)*(minute/(60*seconds))
% w=0:(ns-1);
% w=MyLowestHarmonic*wbase*w.*(w<(ns/2)); 
%
% % Rotor loss
% g1=(g==10);
% ptloss = (transpose ((ch*w+ce*w.*w)*bsq).*g1)/cs;
%
% % Stator loss
% g2=(g==1);
% ptloss = ptloss + (transpose ((ch*w+ce*w.*w)*bsq).*g2)/cs;
% 
% % Prox loss
% g4=(g==2);
% ptloss = ptloss + (transpose ((cePhase*w.*w)*bsq).*g4);
% 
% % PM contribution
% ptloss = ptloss + ((1/2)*(omag*(2*pi*w).^2)*(abs(Jm).^2))';
% 
% %% point location in z
% heating = [real(z),imag(z),ptloss];
% save('c:\\temp\\myLossData','heating','-ASCII');

