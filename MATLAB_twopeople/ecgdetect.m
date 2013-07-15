function [beats] = ecgdetect(x1,fs,threshmultiple,test)

fa = 1;                 % Bandpass filter
fb = 25;                
fc = 1.5;
fd = 0.5;

wa = fa/(fs/2);
wb = fb/(fs/2);
%[B,A] = butter(2,[wa wb]);     %bandpass option
[B,A] = butter(2,wb,'low');     %Low pass option
x1 = filtfilt(B,A,x1);
x2 = x1;
x1 = diff(x1);
x1 = x1.^2;
wc = fc/(fs/2);
[B,A] = butter(2,wc,'low');
x1 = filtfilt(B,A,x1);
wd = fd/(fs/2);
[B,A] = butter(2,wd,'low');
xthresh = filtfilt(B,A,x1);
%find maximum of X1

thresh = max(x1);

j = 1;
for i=2:(length(x1)-1) 
        if (x1(i-1) < x1(i)) && (x1(i) > x1(i+1) && (xthresh(i) < thresh) && (x1(i) >= xthresh(i)))
               mxt(j) = i;
               mx2(j) = x1(i);
           if j > 1
                  mx2mean = mean(mx2);
                  thresh = threshmultiple*mx2mean;
           end    
           j = j + 1;                   
      end    
end 

%Find local maximum or minimum of filter eg (x2)

j = 1;

for i=1:length(mxt)
   
    while(1)
            
        if ( x2(mxt(i)+1) < x2(mxt(i)))
           mxt(i) = mxt(i) - 1;
        else
            mxt(i) = mxt(i) + 1;
        end
        
        if (x2(mxt(i)-1) <  x2(mxt(i))) && (x2(mxt(i)) >  x2(mxt(i) + 1)) 
            break;
        end
        
    end
 
end

mx = x2(mxt);

tau = diff(mxt);
taut = mxt(1:(length(mxt)-1));

%filter non phiosiologically possible values

lb = int32(fs/2); %0.5s => 120bpm
ub = int32(fs*2); %2s => 30bpm

j = 1;

for i=1:length(tau)
   
    if (tau(i) <= ub) && (tau(i) >= lb)
       R_R(j) = tau(i);
       R_t(j) = taut(i);
       j = j + 1;
    end
end

hbtau = x2(R_t);

if test == 1
     
    figure;
    subplot(3,1,1);
    hold on;
    plot(x2);
    plot(R_t,hbtau,'+k');
    hold off;
    subplot(3,1,2);
    hold on;
    plot(x1);
    plot(xthresh,':r');
    plot(mxt,mx2,'+k');
    plot([1 length(x1)],[thresh thresh],'-g');
    hold off;
    subplot(3,1,3);
    plot(R_t,R_R,'-');
    
    
end

beats = [R_t', R_R']/fs;

end