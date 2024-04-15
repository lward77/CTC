% program to simulate new two-signal Kang model
% There are two cortical layers and two regions within each layer
% all modeled by EI pairs oscillating (quasi-cycle) at 40 Hz
% Cortical regions in layer 2 receive inpts from the same region in layer 1
% with a delay of 10 ms - constant phase offset of about \pi?
% Cortical regions in the two layers receive 10 Hz inputs at a variable 
% phase offset from "pulvinar" - this phase offset ideally is adjustable by top
% down signals so as to create a larger input response fromn either of the
% two cortical regions (tuned to one of the two inputs - or signals). This
% was not implemented in the final sims for the PLoS CB paper
% L. Ward 2021-2024

% Define parameters of subpopulations of Excitatory neurons and Inhibitory neurons
clear all;
rng='shuffle';
N=1;           %number of subpopulations/signals
nits=10;        %number of instantiations
runtime = 100000;    % run simulation for runtime moments delt
delt=0.00005;               %time step size
t=0:delt:((runtime-1).*delt);            %time stepsfor sinusoids

S1EE=1.4.*ones(N,1);      %synaptic efficacies - cortical 40 Hz -layer 1
S1EI=0.78.*ones(N,1);      %I on E or E from I
S1IE=2.1105.*ones(N,1);
S1II=0.1.*ones(N,1);
S2EE=1.4.*ones(N,1);      %synaptic efficacies - cortical 40 Hz - layer 2
S2EI=0.78.*ones(N,1);      %I on E or E from I
S2IE=2.1105.*ones(N,1);
S2II=0.1.*ones(N,1);
SPEE1=1.5.*ones(N,1);      %synaptic efficacies -pulvinar 10 Hz - to layer 1
SPEI1=0.6.*ones(N,1);      %I on E
SPIE1=1.0372.*ones(N,1);
SPII1=0.1.*ones(N,1);
SPEE2=1.5.*ones(N,1);      %synaptic efficacies - pulvinar 10 hz - to layer 2
SPEI2=0.6.*ones(N,1);      %I on E
SPIE2=1.0372.*ones(N,1);
SPII2=0.1.*ones(N,1);

S1EPE1=1.5.*ones(N,1);    %input from VPE1/alpha1 to V1E top-down control
S2EPE2=1.5.*ones(N,1);   %input from VPE2/alpha2 to V2E via TRN

tauE=0.003;             %synaptic time constants
tauI=0.006;

sig1E=0.15.*ones(N,1);       %cortical noise intensities
sig1I=0.15.*ones(N,1);
sig2E=0.15.*ones(N,1);
sig2I=0.15.*ones(N,1);

Pnse=0.08;                    %standard noise for pulvinar pairs
sigPE1=Pnse.*ones(N,1);       %pulvinar noise intensities
sigPI1=Pnse.*ones(N,1);
sigPE2=Pnse.*ones(N,1);       %pulvinar noise intensities
sigPI2=Pnse.*ones(N,1);

%voltage variables cortical layer 1 H=horizontal line tuning, V = vertical
V1E=zeros(N,runtime);       %voltage variables
V1I=zeros(N,runtime);       %V1,V2 run at 40 Hz, VP1, VP2 at 10 Hz when no interaction
%V1EV=zeros(N,runtime);       %voltage variables
%V1IV=zeros(N,runtime);

%voltage variables cortical layer 2
V2E=zeros(N,runtime);       %see efficacies above - horizontal line tuning
V2I=zeros(N,runtime);
%V2EV=zeros(N,runtime);       %see efficacies above - vertical line tuning
%V2IV=zeros(N,runtime);

%voltage variables pulvinar to cortical layer 1
VPE1=zeros(N,runtime);       %see efficacies above
VPI1=zeros(N,runtime);
%VPE1V=zeros(N,runtime);       %see efficacies above
%VPI1V=zeros(N,runtime);

%voltage variable pulvinar to cortical layer 2
VPE2=zeros(N,runtime);       %see efficacies above
VPI2=zeros(N,runtime);
%VPE2V=zeros(N,runtime);       %see efficacies above
%VPI2V=zeros(N,runtime);

%differentials of voltage variables cortical layer 1
dV1E=zeros(N,runtime);      
dV1I=zeros(N,runtime);       %V1,V2 run at 40 Hz, VP1, VP2 at 10 Hz when no interaction
%dV1EV=zeros(N,runtime);      
%dV1IV=zeros(N,runtime); 

%differentials of voltage variables cortical layer 2
dV2E=zeros(N,runtime);       %see efficacies above
dV2I=zeros(N,runtime);
%dV2EV=zeros(N,runtime);       %see efficacies above
%dV2IV=zeros(N,runtime);

%differentials of voltage variables pulvinbar to cortical layer 1
dVPE1=zeros(N,runtime);       %see efficacies above
dVPI1=zeros(N,runtime);
%dVPE1V=zeros(N,runtime);       %see efficacies above
%dVPI1V=zeros(N,runtime);

%differentials of voltage variables pulvinar to cortical layer c2
dVPE2=zeros(N,runtime);       %see efficacies above
dVPI1V=zeros(N,runtime);
%dVPE2V=zeros(N,runtime);       %see efficacies above
%dVPI2V=zeros(N,runtime);

V1E(:,1)=rand(N,1);        %start at random positions
V1I(:,1)=rand(N,1);
V2E(:,1)=rand(N,1);        %start at random positions
V2I(:,1)=rand(N,1);
VPE1(:,1)=rand(N,1);        %start at random positions
VPI1(:,1)=rand(N,1);
VPE2(:,1)=rand(N,1);        %start at random positions
VPI2(:,1)=rand(N,1);
%V1EV(:,1)=rand(N,1);        %start at random positions
%V1IV(:,1)=rand(N,1);
%V2EV(:,1)=rand(N,1);        %start at random positions
%V2IV(:,1)=rand(N,1);
%VPE1V(:,1)=rand(N,1);        %start at random positions
%VPI1V(:,1)=rand(N,1);
%VPE2V(:,1)=rand(N,1);        %start at random positions
%VPI2V(:,1)=rand(N,1);

phdiffa=zeros(nits,runtime);
phdiffatopi=zeros(nits,runtime);

twopi=2.*pi;
sqrtdelt=delt^0.5;
nbins=59.0;                 %actually nbins+1 = number of phase bins
binpcoh=zeros(nits,nbins);
pcoh=zeros(nits,1);
V1Earp=zeros(nits,nbins);     %raw amp of V1E as function of alpha phase
V2Earp=zeros(nits,nbins);
V1Iarp=zeros(nits,nbins);
nV1Iarp=zeros(nits,nbins);
alpha1arp=zeros(nits,nbins);
nV1Earp=zeros(nits,nbins);      %numnber of occurrences in alpha phase bins
nV2Earp=zeros(nits,nbins);
aV2Earp=zeros(nits,nbins);
VPE1arp=zeros(nits,nbins);      %raw amplitude of VPE1 as function of alpha phase
aV1Earp=zeros(nits,nbins);     %Hilbert amp of V1E as function of alpha phase
nVPE1arp=zeros(nits,nbins);
alphafact1=3;                 %amount of pulvinar alpha from VPE1 or alpha1 to cortex  
alphafact2=3;             %amount of pulvinar alpha from VPE2 or alpha2 to cortex  
ialpha=1;            %top-down alpha trigger=0 if no alpha to cortex
%alagindex=123;          %>lag of second alpha generator
%alag=120;              %lag of second alpha in msec=alag*delt
tdtrn=-0.5;              %input to TRN that decreases PI processes (in pulvinar)
P2IP1I=0.0;             %input to VP2I from VP1I
nsesig=0.0;
signal(1,1:50000)=0;
signal(1,50001:runtime)=20+nsesig.*randn(N,50000);       %signal - can be noisy
S1E2E=0.5.*ones(N,1);   %input from V1E to V2E Quax


mxV1E=zeros(nits,1);
mnV1E=zeros(nits,1);
mxV2E=zeros(nits,1);
mnV2E=zeros(nits,1);

determ=1;                   % determ=1 for noisy 10 Hz sine waves
                            % determ=0 for EI pair 10 Hz quasi-cycles
amag=2;                     % initial amplitude of 10 Hz sine waves

V2Eoff=125;             % delay to adjust phase offset between cortical 
                        % layers 1 and 2 - same for both input signals
                        % 125 timepoints = 6.25 msec = -\pi/2 radians at 40 Hz
sigon=50000;            % time at which signal input turned on
sigoff=100000;           % time at which signal input turned off

offset=0;           %alpha offset for H-tuned cortical region
%offsetV=pi;             %alpha offset for V-tuned cortical region
noV2=0;             %=1 if no V2 processes - only noise
nV2nse=0.02;         %noise for noV2 processes
    
    for jj=1:nits          % iteration loop
      
      if determ==1
         alpha1=amag.*sin(twopi.*10.*t+(0.4.*randn(1,runtime)));        %alpha inputs to gamma processes
         alpha2=amag.*sin((twopi.*10.*t)+offset+(0.4.*randn(1,runtime)));
         %alpha1V=amag.*sin(twopi.*10.*t+(0.4.*randn(1,runtime+1)));        %alpha inputs to gamma processes
         %alpha2V=amag.*sin((twopi.*10.*t)+offsetV+(0.4.*randn(1,runtime+1)));
      else
         alag=round((offset./(2.*pi)).*2000);     %>lag of second alpha generator
         alagindex=alag+3;          
      end
        
      if ialpha == 0
            S1IP1I=0.0.*ones(N,1);    %input from VPE1/alpha1 to V1E/V1I top-down control
            S2IP2I=0.0.*ones(N,1);    %input from VPE2/alpha2 to V2E/V2I
        else     
            S1IP1I=alphafact1.*ones(N,1);    %input from VPE1/alpha1 to V1E top-down control
            S2IP2I=alphafact2.*ones(N,1);
        end
        
        %S1IE(:,1)=S1IEr;%(S1IPE.*alpha1(:,k))+
        %S2IE(:,1)=S2IEr;%(S2IPI.*alpha2(:,k))+

    noiseE1=0.2*randn(N,runtime);      %independent noises - same for H and V
    noiseI1=0.2*randn(N,runtime);      % layer 1
    noiseE2=0.2*randn(N,runtime);      %independent noises - same for H and V
    noiseI2=0.2*randn(N,runtime);      % layer 2

%1/f noise
%Nx = round(2^16.6097);  % number of samples to synthesize
%B = [0.049922035 -0.095993537 0.050612699 -0.004408786];
%A = [1 -2.494956002   2.017265875  -0.522189400];
%nT60 = round(log(1000)/(1-max(abs(roots(A))))); % T60 est.
%v = randn(1,Nx+nT60); % Gaussian white noise: N(0,1)
%noiseE1 = filter(B,A,v);    % Apply 1/f roll-off to PSD
%noiseE1 = 0.05.*noiseE1(nT60+1:end);    % Skip transient response

%[pxx,fx]=pmtm(noiseE1,[],length(noiseE1),Nx);
%plot(log(pxx),fx);
       
for k=2:runtime      %basic loop
        
    %noiseE1V=0.2.*randn(1,1);      %independent noises - same for H and V
    %noiseI1V=0.2.*randn(1,1);      % layer 1
    %noiseE2V=0.2.*randn(1,1);      %independent noises - same for H and V
    %noiseI2V=0.2.*randn(1,1);      % layer 2
    
    %if k<100000
     %   S1IP1I=alphafact1.*ones(N,1);    %input from VPE1/alpha1 to V1E top-down control
     %   S2IP2I=alphafact2.*ones(N,1);
    %else
      %  S1IP1I=0.0.*ones(N,1);    %input from VPE1/alpha1 to V1E/V1I top-down control
      %  S2IP2I=0.0.*ones(N,1);    %input from VPE2/alpha2 to V2E/V2I
    %end


    if determ==0                   %10 Hz from EI pairs
        noiseEP1H=0.2.*randn(1,1);      %independent noises
        noiseIP1H=0.2.*randn(1,1);
        noiseEP2H=0.2.*randn(1,1);      %independent noises
        noiseIP2H=0.2.*randn(1,1);
    
        %noiseEP1V=0.2.*randn(1,1);      %independent noises
        %noiseIP1V=0.2.*randn(1,1);
        %noiseEP2V=0.2.*randn(1,1);      %independent noises
        %noiseIP2V=0.2.*randn(1,1);
  
    %now compute voltages - solve SDE using Euler-Maruyama
        %EI pulvinar
        if k>sigon && k<sigoff
            dV1E(:,k)=((delt.*(-V1E(:,k-1)+(S1EE(:,1).*V1E(:,k-1))-(S1EI(:,1).*V1I(:,k-1))+signal))+(sig1E(:,1).*sqrtdelt.*noiseE1(:,k)))./tauE;
        else
            dV1E(:,k)=((delt.*(-V1E(:,k-1)+(S1EE(:,1).*V1E(:,k-1))-(S1EI(:,1).*V1I(:,k-1))))+(sig1E(:,1).*sqrtdelt.*noiseE1(:,k)))./tauE;
        end
        dV1I(:,k)=((delt.*(-V1I(:,k-1)-(S1II(:,1).*V1I(:,k-1))+(S1IE(:,1).*V1E(:,k-1))+(S1IP1I(:,1).*VPE1(:,k-1))))+(sig1I(:,1).*sqrtdelt.*noiseI1(:,k)))./tauI;
        if k<=V2Eoff
            dV2E(:,k)=0.0001.*rand(N,1);   %((delt.*(-V2E(:,k-1)+(S2EE(:,1).*V2E(:,k-1))-(S2EI(:,1).*V2I(:,k-1))+(S1E2E.*V1E(:,1))))+(sig2E(:,1).*sqrtdelt.*noiseE2))./tauE;
        else if k>V2Eoff
                if k==V2Eoff+1
                    dV2E(:,k)=dV1E(:,k);
                else
                    dV2E(:,k)=((delt.*(-V2E(:,k-1)+(S2EE(:,1).*V2E(:,k-1))-(S2EI(:,1).*V2I(:,k-1))+(S1E2E.*V1E(:,k-V2Eoff))))+(sig2E(:,1).*sqrtdelt.*noiseE2(:,k)))./tauE;   
                end
            end
        end
        dV2I(:,k)=((delt.*(-V2I(:,k-1)-(S2II(:,1).*V2I(:,k-1))+(S2IE(:,1).*V2E(:,k-1))+(S2IP2I(:,1).*VPE2(:,k-1))))+(sig2I(:,1).*sqrtdelt.*noiseI2(:,k)))./tauI;
        %if alag<=0
            dVPE1(:,k)=((delt.*(-VPE1(:,k-1)+(SPEE1(:,1).*VPE1(:,k-1))-(SPEI1(:,1).*VPI1(:,k-1))))+(sigPE1(:,1).*sqrtdelt.*noiseEP1(1,1)))./tauE;
        %    if k<abs(alag)
         %       dVPE2(:,k)=0.0001.*rand(N,1);
         %   else
             dVPE2(:,k)=((delt.*(-VPE2(:,k-1)+(SPEE2(:,1).*VPE2(:,k-1))-(SPEI2(:,1).*VPI2(:,k-1))))+(sigPE2(:,1).*sqrtdelt.*noiseEP2(1,1)))./tauE;
         %   end
        %else
          %  if k<alag
          %      VPE1(:,k)=0.0001.*rand(N,1);
          %  else
          %      dVPE1(:,k)=((delt.*(-VPE1(:,k-1)+(SPEE1(:,1).*VPE1(:,k-1))-(SPEI1(:,1).*VPI1(:,k-1))))+(sigPE1(:,1).*sqrtdelt.*noiseEP1(1,1)))./tauE;
           % end
            %dVPE2(:,k)=((delt.*(-VPE2(:,k-1)+(SPEE2(:,1).*VPE2(:,k-1))-(SPEI2(:,1).*VPI2(:,k-1))))+(sigPE2(:,1).*sqrtdelt.*noiseEP2(1,1)))./tauE;
         
        dVPI1(:,k)=((delt.*(-VPI1(:,k-1)-(SPII1(:,1).*VPI1(:,k-1))+(SPIE1(:,1).*VPE1(:,k-1))))+(sigPI1(:,1).*sqrtdelt.*noiseIP1(1,1)))./tauI;   
        dVPI2(:,k)=((delt.*(-VPI2(:,k-1)-(SPII2(:,1).*VPI2(:,k-1))+(SPIE2(:,1).*VPE2(:,k-1))))+(sigPI2(:,1).*sqrtdelt.*noiseIP2(1,1)))./tauI;       
    
    else
 
        %noisy alpha sin waves
        if k>=sigon && k<=sigoff   %signal
            dV1E(:,k)=((delt.*(-V1E(:,k-1)+(S1EE(:,1).*V1E(:,k-1))-(S1EI(:,1).*V1I(:,k-1))+signal(:,k-1)))+(sig1E(:,1).*sqrtdelt.*noiseE1(:,k)))./tauE;
        else
            dV1E(:,k)=((delt.*(-V1E(:,k-1)+(S1EE(:,1).*V1E(:,k-1))-(S1EI(:,1).*V1I(:,k-1))))+(sig1E(:,1).*sqrtdelt.*noiseE1(:,k)))./tauE;
        end
        dV1I(:,k)=((delt.*(-V1I(:,k-1)-(S1II(:,1).*V1I(:,k-1))+(S1IE(:,1).*V1E(:,k-1))+(S1IP1I.*alpha1(:,k-1))))+(sig1I(:,1).*sqrtdelt.*noiseI1(:,k)))./tauI;
        
        if noV2==0

        if k<=V2Eoff
            dV2E(:,k)=((delt.*(-V2E(:,k-1)+(S2EE(:,1).*V2E(:,k-1))-(S2EI(:,1).*V2I(:,k-1))))+(sig2E(:,1).*sqrtdelt.*noiseE2(:,k)))./tauE;
        else if k>V2Eoff
            dV2E(:,k)=((delt.*(-V2E(:,k-1)+(S2EE(:,1).*V2E(:,k-1))-(S2EI(:,1).*V2I(:,k-1))+(S1E2E.*V1E(:,k-V2Eoff))))+(sig2E(:,1).*sqrtdelt.*noiseE2(:,k)))./tauE;   
            end
        end
        
        dV2I(:,k)=((delt.*(-V2I(:,k-1)-(S2II(:,1).*V2I(:,k-1))+(S2IE(:,1).*V2E(:,k-1))+(S2IP2I.*alpha2(:,k-1))))+(sig2I(:,1).*sqrtdelt.*noiseI2(:,k)))./tauI;
        else
            if k<=V2Eoff
                dV2E(:,k)=nV2nse.*randn(1,1)+(delt.*(-V2E(:,k-1)+(S1E2E.*V1E(:,k))));
            else
                dV2E(:,k)=nV2nse.*randn(1,1)+(delt.*(-V2E(:,k-1)+S1E2E.*V1E(:,k-V2Eoff)));
            end
            dV2I(:,k)=nV2nse.*randn(1,1);
        
        %dVPE(:,k)=((delt.*(-VPE(:,k-1)+(SPEE(:,1).*VPE(:,k-1))-(SPEI(:,1).*VPI(:,k-1))))+(sigPE(:,1).*sqrtdelt.*noiseEP))./tauE;
        %dVPI(:,k)=((delt.*(-VPI(:,k-1)-(SPII(:,1).*VPI(:,k-1))+(SPIE(:,1).*VPE(:,k-1))))+(sigPI(:,1).*sqrtdelt.*noiseIP))./tauI;
        end
    end
        V1E(:,k)=V1E(:,k-1)+dV1E(:,k);
        V1I(:,k)=V1I(:,k-1)+dV1I(:,k);
        V2E(:,k)=V2E(:,k-1)+dV2E(:,k);
        V2I(:,k)=V2I(:,k-1)+dV2I(:,k);
        if determ==0
            if alag<=0
                VPE2(:,k)=VPE2(:,k-1)+dVPE2(:,k);
                VPI2(:,k)=VPI2(:,k-1)+dVPI2(:,k);
                if k<=abs(alag)
                    VPE1(:,k)=0.0001.*rand(N,1);
                    VPI1(:,k)=0.0001.*rand(N,1);
                else
                    VPE1(:,k)=VPE2(:,k+alag)+0.1.*randn(N,1);
                    VPI1(:,k)=VPI2(:,k+alag)+0.1.*randn(N,1);
                end
            else         %alag>0
          
                VPE1(:,k)=VPE1(:,k-1)+dVPE1(:,k);
                VPI1(:,k)=VPI1(:,k-1)+dVPI1(:,k);
                if k<=abs(alag)
                    VPE2(:,k)=0.0001.*rand(N,1);
                    VPI2(:,k)=0.0001.*rand(N,1);
                else
                    VPE2(:,k)=VPE1(:,k-alag)+0.1.*randn(N,1);
                    VPI2(:,k)=VPI1(:,k-alag)+0.1.*randn(N,1);
                end     
            end
          
            VPE1(:,k)=VPE1(:,k-1)+dVPE1(:,k);
            VPI1(:,k)=VPI1(:,k-1)+dVPI1(:,k);
            VPE2(:,k)=VPE2(:,k-1)+dVPE2(:,k);
            VPI2(:,k)=VPI2(:,k-1)+dVPI2(:,k);
        end
               
    end           %basic runtime loop index=k
    
indx=0;
for g=1:100:runtime
    indx=indx+1;
    V2Edown(1,indx)=V2E(1,g);
    V1Edown(1,indx)=V1E(1,g);
end


  %filter V1E, V2E at gamma and alpha-note only 10 Hz labeled after this
bpdownfil=fir1(100,[0.1,0.2]);
bpdown10fil=fir1(100,[0.01,0.011]);
%bp40fil=fir1(50,0.0035,'high');
%bp10fil=fir1(50,0.0015);
%V1Efilt=filter(bp40fil,1,V1E);
%V1E10filt=filter(bp10fil,1,V1E);
%V2Efilt=filter(bp40fil,1,V2E);
%V2E10filt=filter(bp10fil,1,V2E);
V2Edownfilt=filter(bpdownfil,1,V2Edown);
V2Edown10filt=filter(bpdown10fil,1,V2Edown);
V1Edownfilt=filter(bpdownfil,1,V1Edown);
V1Edown10filt=filter(bpdown10fil,1,V1Edown);

hV1E=hilbert(V1E);     
hV2E=hilbert(V2E);
hV1I=hilbert(V1I);
hV2I=hilbert(V2I);
if determ==0
    hVPE1=hilbert(VPE1);
    hVPE2=hilbert(VPE2);
    hVPI1=hilbert(VPI1);
    hVPI2=hilbert(VPI2);
else
    halpha1=hilbert(alpha1);
    halpha2=hilbert(alpha2);
end

%halpha1=hilbert(alpha1);
ampV1E=abs(hV1E);
phaseV1E=angle(hV1E);
ampV1I=abs(hV1I);
phaseV1I=angle(hV1I);
ampV2E=abs(hV2E);
phaseV2E=angle(hV2E);
ampV2I=abs(hV2I);
phaseV2I=angle(hV2I);

EI1phdiff=wrapToPi(phaseV1E-phaseV1I);
EI2phdiff=wrapToPi(phaseV2E-phaseV2I);

if determ==0
    ampVPE1=abs(hVPE1);
    phaseVPE1nat=angle(hVPE1);
    phaseVPE1=wrapToPi(angle(hVPE1));
    ampVPE2=abs(hVPE2);
    phaseVPE2=angle(hVPE2);
    ampVPI1=abs(hVPI1);
    phaseVPI1=angle(hVPI1);
    ampVPI2=abs(hVPI2);
    phaseVPI2=angle(hVPI2);
else
    ampalpha1=abs(halpha1);
    phasealpha1nat=wrapToPi(angle(halpha1));
    phasealpha1=wrapToPi(angle(halpha1));
    ampalpha2=abs(halpha2);
    phasealpha2=wrapToPi(angle(halpha2));    
end

hV1Efilt=hilbert(V1Edownfilt);   %
hV2Efilt=hilbert(V2Edownfilt);
ampV1Efilt=abs(hV1Efilt);
phaseV1Efilt=angle(hV1Efilt);
ampV2Efilt=abs(hV2Efilt);
phaseV2Efilt=angle(hV2Efilt);

hV1E10filt=hilbert(V1Edown10filt);
hV2E10filt=hilbert(V2Edown10filt);
ampV1E10filt=abs(hV1E10filt);
phaseV1E10filt=angle(hV1E10filt);
ampV2E10filt=abs(hV2E10filt);
phaseV2E10filt=angle(hV2E10filt);

V1V2diff=zeros(2000,nbins);
nV1V2diff=zeros(1,nbins);
pinc=6.28318./nbins;
for j=1:runtime
    lowp=-3.14159;
    hip=lowp+pinc;
    for i=1:nbins
      if determ==0
        if (phaseVPE1nat(j) >= lowp) && (phaseVPE1nat(j) < hip)
            V1Earp(jj,i)=V1Earp(jj,i)+V1E(j);         
            nV1V2diff(1,i)=nV1V2diff(1,i)+1;
            V1V2diff(nV1V2diff(1,i),i)=phaseV1E(1,j)-phaseV2E(1,j);            
            VPE1arp(jj,i)=VPE1arp(jj,i)+VPE1(j);
            aV1Earp(jj,i)=aV1Earp(jj,i)+ampV1E(j);
            nV1Earp(jj,i)=nV1Earp(jj,i)+1;
            V1Iarp(jj,i)=V1Iarp(jj,i)+V1I(j);           
            nV1Iarp(jj,i)=nV1Iarp(jj,i)+1;
            nVPE1arp(jj,i)=nVPE1arp(jj,i)+1;           
        end
      else
         if (phasealpha1(j) >= lowp) && (phasealpha1(j) < hip)
            V1Earp(jj,i)=V1Earp(jj,i)+V1E(j);          
            nV1V2diff(1,i)=nV1V2diff(1,i)+1;
            V1V2diff(nV1V2diff(1,i),i)=phaseV1E(1,j)-phaseV2E(1,j);            
            %VPE1arp(jj,i)=VPE1arp(jj,i)+VPE1(j);
            aV1Earp(jj,i)=aV1Earp(jj,i)+ampV1E(j);
            nV1Earp(jj,i)=nV1Earp(jj,i)+1;
            V1Iarp(jj,i)=V1Iarp(jj,i)+V1I(j);            
            nV1Iarp(jj,i)=nV1Iarp(jj,i)+1; 
            %nVPE1arp(jj,i)=nVPE1arp(jj,i)+1;
         end
         
      end
    
        if determ==0
           if (phaseVPE2(j) >= lowp) && (phaseVPE2(j) < hip)
            V2Earp(jj,i)=V2Earp(jj,i)+V2E(j);           
            nV2Earp(jj,i)=nV2Earp(jj,i)+1;                   
            end
        else
           if (phasealpha2(j) >= lowp) && (phasealpha2(j) < hip)            
            V2Earp(jj,i)=V2Earp(jj,i)+V2E(j);           
            nV2Earp(jj,i)=nV2Earp(jj,i)+1;           
           end
        end   
     
        lowp=lowp+pinc;
        hip=hip+pinc;
    end   %end of nbins loop index i

end                 %end of runtime loop index j

    lowp=-3.14159;
    hip=lowp+pinc;
    for i=1:nbins
        midp(i)=(lowp+hip)./2;
        lowp=lowp+pinc;
        hip=hip+pinc;
    end
    
%compute phase differences and phase coherences
wphaseV1E=wrapToPi(phaseV1E);
wphaseV2E=wrapToPi(phaseV2E);
wphaseV1Efilt=wrapToPi(phaseV1Efilt);
wphaseV2Efilt=wrapToPi(phaseV2Efilt);
wphaseV1E10filt=wrapToPi(phaseV1E10filt);
wphaseV2E10filt=wrapToPi(phaseV2E10filt);

wwphdiff=wrapToPi(wphaseV2E-wphaseV1E);
wwphdifffilt=wrapToPi(wphaseV1Efilt-wphaseV2Efilt);
wwphdiff10filt=wrapToPi(wphaseV1E10filt-wphaseV2E10filt);
phdiff=phaseV2E-phaseV1E;
wphdiff=wrapToPi(phdiff);
filtphdiff=phaseV1Efilt-phaseV2Efilt;
wfiltphdiff=wrapToPi(filtphdiff);
filt10phdiff=phaseV1E10filt-phaseV2E10filt;
wfilt10phdiff=wrapToPi(filt10phdiff);
pcoh(jj,1)=(1./runtime).*abs(sum(exp(1i*phdiff)));
filtpcoh(jj,1)=(1./runtime).*abs(sum(exp(1i*filtphdiff)));
wfiltpcoh(jj,1)=(1./runtime).*abs(sum(exp(1i*wfiltphdiff)));
filt10pcoh(jj,1)=(1./runtime).*abs(sum(exp(1i*filt10phdiff)));
wfilt10pcoh(jj,1)=(1./runtime).*abs(sum(exp(1i*wfilt10phdiff)));
wpcoh(jj,1)=(1./runtime).*abs(sum(exp(1i*wphdiff)));

%for k=1:runtime
   % if phdiffa(1,k)<0
   %     phdiffa0topi(jj,k)=phdiffa(1,k)+twopi;
   % else
   %     phdiffa0topi(jj,k)=phdiffa(1,k);
   % end
%end
if determ==0
    phdiffa=phaseVPI1-phaseVPI2;
else
    phdiffa=phasealpha1-phasealpha2;
end
    pcoha(jj,1)=(1./runtime).*abs(sum(exp(1i*phdiffa)));


%compute phase coherence by pulv phase bin
for i=1:nbins
    binpcoh(jj,i)=(1./nV1V2diff(1,i)).*abs(sum(exp(1i*V1V2diff(:,i)')));
end

%corrcoef(binpcoh',V1Earp')

for i=1:nbins
    if nV1Earp(jj,i) >0
        V1Earp(jj,i)=V1Earp(jj,i)./nV1Earp(jj,i);
        aV1Earp(jj,i)=aV1Earp(jj,i)./nV1Earp(jj,i);        
    else 
        V1Earp(jj,i)=0;
        aV1Earp(jj,i)=0;
    end
    if nV1Iarp(jj,i) >0
        V1Iarp(jj,i)=V1Iarp(jj,i)./nV1Iarp(jj,i);  
    else 
        V1Iarp(jj,i)=0;
    end
    if nV2Earp(jj,i) >0
        V2Earp(jj,i)=V2Earp(jj,i)./nV2Earp(jj,i);
        aV2Earp(jj,i)=aV2Earp(jj,i)./nV2Earp(jj,i);        
    else 
        V2Earp(jj,i)=0;
        aV2Earp(jj,i)=0;
    end
end  %end of nbins loop index i

%compute max and mean signal effects - signal on at 50000
[mxV1E(jj,1),me1]=max(subplus(V1E(1,50000:55000)));
mnV1E(jj,1)=mean(subplus(V1E(1,me1+50000-1000:me1+50000+1000)));
[mxV2E(jj,1),me2]=max(subplus(V2E(1,50000:55000)));
mnV2E(jj,1)=mean(subplus(V2E(1,me2+50000-1000:me2+50000+1000)));
%mnalphoff(jj,1)=mean(wrapToPi(phaseVPI1-phaseVPI2));

rawsigeffmax(jj,1)=mxV2E(jj,1)-mxV1E(jj,1);   %sig onset effects by iteration
rawsigeffmn(jj,1)=mnV2E(jj,1)-mnV1E(jj,1);

%compute mutual information    
mutualinfpre(1,jj)=mi(V1E(1,1:50000)',V2E(1,1:50000)');
mutualinfpost(1,jj)=mi(V1E(1,50001:100000)',V2E(1,50001:100000)');
mutualamppre(1,jj)=mi(ampV1E(1,1:50000)',ampV2E(1,1:50000)');
mutualamppost(1,jj)=mi(ampV1E(1,50001:100000)',ampV2E(1,50001:100000)');
mutualphpre(1,jj)=mi(phaseV1E(1,1:50000)',phaseV2E(1,1:50000)');
mutualphpost(1,jj)=mi(phaseV1E(1,50001:100000)',phaseV2E(1,50001:100000)');
jj
  end         %end of iteration loop jj index

 
sdrawmaxsoe=std(rawsigeffmax);         %mean and sd of signal onset effects
sdrawmnsoe=std(rawsigeffmn);
meanmaxsoe=mean(rawsigeffmax);
meanmnsoe=mean(rawsigeffmn);
  

mutualinfpremn=mean(mutualinfpre);
    mutualinfpostmn=mean(mutualinfpost);
    mutualamppremn=mean(mutualamppre);
    mutualamppostmn=mean(mutualamppost);
    mutualphpremn=mean(mutualphpre);
    mutualphpostmn=mean(mutualphpost);
mutualinfprese=std(mutualinfpre)/sqrt(nits);
    mutualinfpostse=std(mutualinfpost)/sqrt(nits);
    mutualampprese=std(mutualamppre)/sqrt(nits);
    mutualamppostse=std(mutualamppost)/sqrt(nits);
    mutualphprese=std(mutualphpre)/sqrt(nits);
    mutualphpostse=std(mutualphpost)/sqrt(nits);

  
  
  if jj>1
    V1Earpave=mean(V1Earp);
    V2Earpave=mean(V2Earp);
    V1Iarpave=mean(V1Iarp);
    binpcohave=mean(binpcoh);
    V1Earpsd=std(V1Earp);
    V2Earpsd=std(V2Earp);
    V1Iarpsd=std(V1Iarp);
    binpcohsd=std(binpcoh);
  else 
    V1Earpave=V1Earp;
    V2Earpave=V1Earp;
    V1Iarpave=V1Iarp;
    binpcohave=binpcoh;
    V1Earpsd=0;
    V2Earpsd=0;
    V1Iarpsd=0;
    binpcohsd=binpcoh;
  end
  
V1Earpaveo(:,:)=V1Earpave(:,:)./max(V1Earpave(:,:));
V2Earpaveo(:,:)=V2Earpave(:,:)./max(V2Earpave(:,:));
V1Iarpaveo(:,:)=V1Iarpave(:,:)./max(V1Iarpave(:,:));


Fs = 20000;                  % Sampling frequency in time: samples/sec
L = runtime;                    % Length of signal in time
NFFT = 2.^nextpow2(L);       % Next power of 2 from length of y
Y = fft(V1E,NFFT)/L;
%do FFT on V1E
f = Fs/2*linspace(0,1,NFFT/2+1);    %get X-axis values
YY=2.*abs(Y(1:NFFT/2+1));
Y2 = fft(V2E,NFFT)/L;
%do FFT on V2E
YY2=2.*abs(Y2(1:NFFT/2+1));
Y5=fft(V1I,NFFT)/L;
YY5=2.*abs(Y5(1:NFFT/2+1));

[pxx,fx]=pmtm(V1E,[],length(V1E),Fs);
[pxy,fy]=pmtm(V2E,[],length(V2E),Fs);
%[pxx10,fx10]=pmtm(V1Edown10filt,[],length(V1Edown10filt),Fs);
%[pxy10,fy10]=pmtm(V2Edown10filt,[],length(V2Edown10filt),Fs);
[pxydown,fydown]=pmtm(V2Edownfilt,[],length(V2Edownfilt),1000);
[pxydown10,fydown10]=pmtm(V2Edown10filt,[],length(V2Edown10filt),1000);

f1=figure; hold on
plot(fx(20:400),pxx(20:400));   %V1E all
plot(fy(20:400),pxy(20:400));    %V2E  all
%f2=figure;
%plot(fx10(20:400),pxx10(20:400));  %V1E 10Hz
%f3=figure;
%plot(fy10(20:400),pxy10(20:400));   %V2E 10 Hz
%f4=figure;
%plot(fydown(20:400),pxydown(20:400));
%f5=figure;
%plot(fydown10(20:400),pxydown10(20:400));

%time-frequency analysis - difficult to get time-freq tradeoff right
%fs = 20000;
%t = 0:1/fs:2;
%AA=pspectrum(V1Erect,fs,'FrequencyResolution',10);

%pspectrum(V1Erect,fs,'spectrogram','TimeResolution',0.1)


%[cxy,fc] = mscohere(V1E,V2E,hamming(512),500,2048);
%plot(fc/pi,cxy);

Bmeanpcohgamma=mean(pcoh)
wBmeanpcohgamma=mean(wpcoh)
sdpcohgamma=std(pcoh)
BBmeanpcohalpha=mean(pcoha)
%BBBmeanfiltpcohgamma=mean(filtpcoh)
C=mean(binpcohave)
if determ==0
    A=phaseVPE1-phaseVPE2;
else
    A=phasealpha1-phasealpha2;
end

[ft,gof]=fit(midp',V1Earpaveo','smoothingspline');
[ft1,gof]=fit(midp',V1Iarpaveo','smoothingspline');
f6=figure;
plot(ft,midp,V1Earpaveo);hold on;
plot(ft1,midp,V1Iarpaveo);hold off;
f7=figure;
[ft2,gof2]=fit(midp',V2Earpaveo','smoothingspline');
plot(ft2,midp,V2Earpaveo);

%plot(midp,VPEarp./max(VPEarp));hold on
%f4=figure;
%plot(midp,binpcohave);
%plot(midp,binpcohave+binpcohsd);
%plot(midp,binpcohave-binpcohsd);hold off;
%f3=figure;
%plot(V1E);hold on;
%if determ==0
%    plot(VPE1);hold off;
%else
%    plot(alpha1);hold off;
%end
%f4=figure;
%plot(V2E);hold on;
%if determ==0
%    plot(VPE2);hold off;
%else
%    plot(alpha2);hold off;
%end
%f5=figure;
%plot(f(1,1:1000),YY(1,1:1000));hold on;
%plot(f(1,1:1000),YY2(1,1:1000));hold off;

%f6=figure;
%plot(f(1,1:1000),YY4(1,1:1000));

D=V1E(1,5000:10000);
DD=V1I(1,5000:10000);
if determ==0
    E=VPE1(1,5000:10000);
else
    E=alpha1(1,5000:10000);
end
%f7=figure;
%plot(D); hold on
%plot(DD);
%plot(E);

F=V2E(1,5000:10000);
FF=V2I(1,5000:10000);
if determ==0
    G=VPE2(1,5000:10000);
else
    G=alpha2(1,5000:10000);
end
f8=figure;
plot(F); hold on
plot(FF);
plot(G);

f9=figure;
histogram(wphdiff);

f10=figure; hold on;
plot(V1E(1,:),'b');
plot(V2E(1,:),'r');
axis([0 10e4 -60 80]);
title(['Signal onset effects']);
xlabel('Timepoints');
ylabel('V1E and V2E');
zlabel('Time');
legend('V1E','V2E');

