%Creates subdata files from the matlab data file 

clear all;
close all;

%Load global variables
filename = ['..\..\Data\global'];
load(filename);

%Load data
load(filenameFout);

n = find(value1 == - 8388608,1,'first') - 1;

for i=1:(n - 1) 
   if time(i+1) < time(i)  
       time(i+1) =  8388608 + time(i+1); 
   end
end

dtime = diff(time);
locs = find(dtime > 1);

%Prefilter to remove 50Hz

f1 = 20;
w1 = f1/(fs/2);

[B,A] = butter(4,w1,'low');
out1 = filtfilt(B,A,value1(1:n));
out2 = filtfilt(B,A,value2(1:n));

%Calibration data
load(filenameCalib);

F1 = (out1 - c1).*m1;
F2 = (out2 - c2).*m2;

clear out1;
clear out2;
clear value1;
clear value2;

save(filenameF,'F1','F2');


%BCG signal

wlo = 1.5/(fs/2);
whi = 8/(fs/2);
[B,A] = butter(2,[wlo whi]);
Fbcg1 = filtfilt(B,A,F1);
Fbcg2 = filtfilt(B,A,F2);

save(filenameBCG,'Fbcg1','Fbcg2','time');
clear Fbcg1;
clear Fbcg2;

%Breathing signal
wlo = 0.05/(fs/2);
whi = 0.5/(fs/2);
[B,A] = butter(2,[wlo whi]);
Fbr1 = filtfilt(B,A,F1);
Fbr2 = filtfilt(B,A,F2);

save(filenameBR,'Fbr1','Fbr2','time');

clear Fbr1;
clear Fbr2;

%DC signal for mass
f1 = 0.2;
w1 = f1/(fs/2);

[B,A] = butter(2,w1,'low');
Fdc1 = filtfilt(B,A,F1);
Fdc2 = filtfilt(B,A,F2);

save(filenameDC,'Fdc1','Fdc2','time');

clear Fdc1;
clear Fdc2;









