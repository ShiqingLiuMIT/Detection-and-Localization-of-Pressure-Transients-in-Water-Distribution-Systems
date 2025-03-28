% function of modified CUSUM, for anomaly detection in time-series signal
function [ta,tai,tai1]=modcusum(x,denoise,drift_p,threshold_p,feq,plott)

% x is signal, 1xN vector
% denoise is to decide if denoise the signal first
% denoise=0, do not denoise
% denoise=1, use function wdenoise
% denoise=2.x, use moving average, with a window of x
% e.g., denoise=2.3, moving average denoise with a window size of 3,
% e.g., denoise=2.21, moving average denoise with a window size of 21.
% drift_p and threshold_p are parameters
% feq is the frequency of the sample, only apply when plott=1
% plott, if plot or not, =1 is plot.


if denoise==0
    xde=x;
end
if denoise==1
    xde=wdenoise(x);
end
if denoise>2 && denoise <3
    decimal_part = denoise - floor(denoise);
    mmwindow = str2double(regexprep(num2str(decimal_part), '0\.', ''));
    % disp(mmwindow);
    xde=movmean(x,mmwindow);
end


drift=drift_p;
threshold=threshold_p;


N=numel(xde);
gp=zeros(N,1);gn=zeros(N,1);
gp_real=zeros(N,1);gn_real=zeros(N,1);
ta=[];
tai=[];
taf=[];
tai_o=[];
tap=1;
tan=1;


for i=2:N
   s(i)=xde(i)-xde(i-1);
   gp(i)=gp(i-1)+s(i)-drift;
   gp_real(i)=gp_real(i-1)+s(i);
   gn(i)=gn(i-1)-s(i)-drift;
   gn_real(i)=gn_real(i-1)-s(i);
   
   if gp(i)<0
       gp(i)=0;gp_real(i)=0;tap=i;
   end
   if gn(i)<0
       gn(i)=0;gn_real(i)=0;tan=i;
   end
   if gp_real(i)>threshold
       ta=[ta;i];
       tai=[tai;tap];
       gp_real(i)=0;
   else
       if gn_real(i)>threshold
           ta=[ta;i];
           tai=[tai;tan];
           gn_real(i)=0;
       end
   end
end

if isempty(tai)
   tai1=numel(xde);
else
   tai1=tai(1);
end

timestep=1/feq;

if plott==1
    figure
    p=plot(linspace(0,timestep*(numel(x)-1),numel(x)),x);
    hold on
    % p1=plot(linspace(0,timestep*(numel(xde)-1),numel(xde)),xde);
    p2=plot((tai1-1)*timestep,xde(tai1),'rd',"Markersize",3,'MarkerFaceColor','red');
    p3=plot((tai1-1)*timestep+10*timestep,xde(tai1),'rd',"Markersize",3,'MarkerFaceColor','red');
    p4=plot((tai1-1)*timestep-10*timestep,xde(tai1),'rd',"Markersize",3,'MarkerFaceColor','red');

    text((tai1-1)*timestep,xde(tai1),{' Wave Front Start'},'FontSize',6,'Color','r');
        xlabel("Time(s)")
        ylabel("Pressure(m)")

    
end


end