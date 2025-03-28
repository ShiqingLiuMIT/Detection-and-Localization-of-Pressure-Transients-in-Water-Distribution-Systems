% simulation result analysis. Find the accurate/inaccurate detection range of events


%% add noise to signal by function awgn

% original noise free data from simulation, the result has N cells, each cell is pressure at nodes of a simulation
data=ky2_pd5m; 

numsource=numel(data); % number of sources in the simulation results

data_awgn10=cell(1,numsource); % create the noisy data variable

for i=1:numsource
    data_awgn10{i}=awgn(data{i},10); % add white gaussian noise to the noise-free data
end 
ky2_pd5m_awgn10=data_awgn10;% noisy data variable

% note: Due to the randomness of the 'awgn' function, the generated noisy signals are not exactly the same each time

%% 
% find the detected arrival times of noisy signals. 
% compare the arrival times with theoretical shortest times from sources
% decide if the arrival times are detected, and if yes, if the times is accurate/inaccurate 

simdataset=ky2_pd5m_awgn10; % input dataset

clear    ST_adj ST_det      ismatch_temp
numsource=810;% for ky2, its 810. ky15=577, ky2=810, ky18 763, ky10 884
numrec=810;% for ky2, its 810. %ky15=661, ky2=810, ky18 763 ky10 922
ShortestT=ShortestT2; %shortest time in network KY2
goodsource=goodsource2_pd5m; % only calculate the good source (source with correct signal)
simdur=10; % simulation duration in seconds.

% a variable to record if the detected time match the theoretical time
ismatch_temp=zeros(numsource,numrec);

tic

%parameter for modcusum
drift_mul=0.0143; % drift of the modcusum algorithm
threshold_mul=0.82; 
denoise=2.21; 

for sourcenode=goodsource % for each source node
    sourcenode
    simudata_sig=simdataset{sourcenode}'; % signal matrix
    simudata_sig=round(simudata_sig,3);
    timestep=0.005;            
    % for each node, find its transient arrival, transient end... from
    % simulation data
    ST_det=zeros(1,numrec); % detection result vector
    for k=1:numrec
        sig=simudata_sig(:,k); % signal at sensor k

        % find the arrival time
        [~,~,ta]=modcusum(sig,denoise,drift_mul,threshold_mul,200,0);
        %adjust shortest path time to 10 if too far away
        if ShortestT(sourcenode,k)>simdur-0.1 % far away nodes
            ismatch_temp(sourcenode,k)=3;
            ST_adj(k)=simdur; %shortest time, adjusted for t>10 nodes
            ST_det(k)=simdur+0.1; % Detected arrival time
        else
            ST_adj(k)=ShortestT(sourcenode,k); % theoretical shortest time 
            ST_det(k)=ta*timestep+0.5*timestep; % arrival time from detection result
            % check the accuracy of detected arrival time, the first 0.1 is the transient start time, the second 0.1 is
            % accuracy threshold
            if abs(ST_det(k)-ST_adj(k)-0.1)<0.1 % if accurate
                ismatch_temp(sourcenode,k)=1; % mark the accuracy result as 1
            elseif ST_det(k)<simdur %if not accurate, but detected
                ismatch_temp(sourcenode,k)=2; % mark the accuracy result as 2
            else % if not accurate, not detected at all
                ismatch_temp(sourcenode,k)=4; % mark the accuracy result as 4
            end
        end
    end
end

toc



disp(["Number of Accurate Detected Signals" numel(find(ismatch_temp==1))])
disp(["Number of Inaccurate Detected Signals" numel(find(ismatch_temp==2))])

ismatch_ky2_pd5m_10=ismatch_temp;