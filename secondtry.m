%{
first read from folder each file y
read data from each file 
call the find function for all
    in this do filtering 
    do further based on power
    find the note and place it in an array
    put that array in the second cell of matrix
check with function return typeshttps://github.com/charithaskt/guitar-note-recoznition-data-ground-truth.git
print the cell array
%}
%%
clear, clc, close all
abc=input('enter 1 for folder and 0 for a single file');
if abc==1
myFolder=input('Enter folder path:','s');
%myFolder = '~/matlabfiles/project/GroundTruth/';
%myFolder = '~/matlabfiles/project/Guitar_Only/a';
%myFolder = '~/matlabfiles/project/try';
if ~isfolder(myFolder)
  errorMessage = sprintf('Error: The following folder does not exist:\n%s', myFolder);
  uiwait(warndlg(errorMessage));
  return;
end
count=0;
filePattern = fullfile(myFolder, '*.wav');
theFiles = dir(filePattern);
matrix=cell(length(theFiles),2);
for jk = 1 : length(theFiles)
  baseFileName = theFiles(jk).name;
  matrix{jk,1}=baseFileName;
  fullFileName = fullfile(myFolder, baseFileName);
  matrix{jk,2}= find(fullFileName,abc);
end
%celldisp(matrix)
%matrix{2,2}{2}
for ij=1:length(matrix)
    fprintf("%40s      -",matrix{ij,1});
    if length(matrix{ij,2})==1
        count=count+1;
    end
    for jk=1:length(matrix{ij,2})
        fprintf(" %s  ",matrix{ij,2}{jk});
    end
    fprintf("\n");
end
fprintf("No.of files is");
disp(length(matrix));
count=(count/length(matrix))*100;
fprintf("%2.3d percent data showed perfect outputs\n",double(count));
elseif abc==0
        file=input('Enter file complete path:','s');
        if ~isfile(file)
            errorMessage = sprintf('Error: The following file does not exist:\n%s', file);
            uiwait(warndlg(errorMessage));
            return;
        else
            out=find(file,abc);
            fprintf("notes for the file %s are:\n",file);
            for jk=1:lengdisp(snr(x1));
th(out)
                fprintf("%s ",out{jk});
            end
            fprintf("\n");
        end    
else
    fprintf("please rerun code as you have entered neither 0 or 1\n");
    return;
end
%%
%[x1,Fs] = audioread('./80/07-magg.wav');
%[x1,Fs]=audioread('~/Downloads/Titanic.wav');
%[x1,Fs] = audioread('./tones/Gsharp-28.wav');
%info=audioinfo('./beatles_14.wav');
%x1 = x1(:,1);
%% making a bandpass filter
%{
n=7;
low_freq=15.1/(Fs/2);
highfreq=700/(Fs/2);
[z,p,k]=butter(n,[low_freq,highfreq],'bandpass');
[sos,glp] = zp2sos(z,p,k)  ;       % Obtain second-order sections form
%fvtool(sos) ;
xout=filtfilt(sos,glp,x1);
%x1out=bpfilt(x1,15,700,Fs,0);      %not used
%}
%%
% [S,f,TT,POW,F,TF]=spectrogram(x1out,[],[],[],Fs);%floor(info.Duration)*2
% %f=zeros(floor(info.Duration)*2);
% POW
% frequency=F(F<700);
% fre=zeros(length(S(1,:)));
% for i=1:length(fre)
% p=max(POW(:,i));
% fre(i)=F(POW(:,i)==p);
% end
%%
%{
[POWER,F]=pwelch(xout,[],[],[],Fs);%computes power density vs frequency power density=power/frequency
p1=POWER;
pavg1=(max(p1)-min(p1))/8;
%pav=(max(p1)-max(p1(~(p1==max(p1)))))
%pavg1=max(p1)-6*pav;
pavg=sum(p1(p1>=pavg1))/length(p1(p1>=pavg1));
final_fr=round((F(p1>=pavg1)));
%final_fr=floor(F(p1>pavg/1.45));
final_freq=final_fr(final_fr<700)
j=findforsonglist(final_freq);
%}
%used
%%
% Nsamps = length(xout);% Layer 2
% t = (1/Fs)*(1:Nsamps);          %Prepare time data for plot
% 
% %Do Fourier Transform
% y_fft = abs(fft(xout));            %Retain Magnitude
% y_fft = y_fft(1:Nsamps/2);      %Discard Half of Points
% ff = Fs*(0:Nsamps/2-1)/Nsamps;
% ff1=ff(ff<700);
% p2=y_fft;
% max(p2)
% max(p2(~(p2==ma% Nsamps = length(xout);% Layer 2
% t = (1/Fs)*(1:Nsamps);          %Prepare time data for plot
% 
% %Do Fourier Transform
% y_fft = abs(fft(xout));            %Retain Magnitude
% y_fft = y_fft(1:Nsamps/2);      %Discard Half of Points
% ff = Fs*(0:Nsamps/2-1)/Nsamps;
% ff1=ff(ff<700);
% pav2=(max(p2)-max(p2(~(p2==max(p2)))))
% pavg2=max(p2)-6*pav2;
% pavg22=sum(p2(p2>=pavg2))/length(p2(p2>=pavg2));
% final_fr2=floor((F(p2>=pavg2)));
% %final_fr=floor(F(p1>pavg/1.45));
% final_freq2=final_fr2(final_fr2<700);
%pavg2=(max(p2)+min(p2))/10;
%pavg1=(max(p1)-min(p1))/16;
%pavg22=sum(p2(p2>=pavg2))/length(p2(p2>=pavg2));
%final_fr2=floor((ff(p2>=pavg22)));
%final_fr=floor(F(p1>pavg/1.45));
%final_freq2=final_fr2(final_fr2<700);
% final_fff2=zeros(ceil(length(y_fft)/length(ff(ff<=5))),1);
% j=1;
% k=length(ff<=5);
% y_fft(j:k);
% for i=1:length(final_fff2)
%     final_fff2(i)=ff(y_fft==max(y_fft(j:k)));
%     j=j+length(ff<=5);
%     k=k+length(ff<=5);
%  end
%  final_ff2=final_fff2(final_fff2<700);
% j=findforsonglist(final_freq2);
%Prepare freq data for plot
%Plot Sound File in Time Domain
% figure
% plot(t, xout);
% xlabel('Time (s)')
% ylabel('Amplitude')
% title('Tuning Fork A4 in Time Domain')
% 
% %Plot Sound File in Frequency Domain
%  figure
%  plot(ff, y_fft);
%  p2=y_fft;
% [pks,locs]=findpeaks(y_fft,Fs*5);
% pavg2=(max(p2)+min(p2))/10;
% %pavg1=(max(p1)-min(p1))/16;
% pavg22=sum(p2(p2>=pavg2))/length(p2(p2>=pavg2));
% final_fr2=floor((ff(p2>=pavg22)));
% %final_fr=floor(F(p1>pavg/1.45));
% final_freq2=final_fr2(final_fr2<700)
% final_fff2=zeros(length(pks));
% for i=1:length(pks)
% final_fff2(i)=ff(y_fft==pks(i));
% end
% j=findforsonglist(final_freq2);

% xlim([0 1000]);
% xlabel('Frequency (Hz)')
% ylabel('Amplitude')
% title('Frequency Response of Tuning Fork A4')
%%
function a=find(f,abc)
[x1,Fs] = audioread(f);
%disp(snr(x1));
%info=audioinfo('./beatles_14.wav');
n=8;
x1 = x1(:,1);
low_freq=15.1/(Fs/2);
highfreq=700/(Fs/2);
[z,p,k]=butter(n,[low_freq,highfreq],'bandpass');
[sos,glp] = zp2sos(z,p,k)  ;       % Obtain second-order sections form
%fvtool(sos) ;
x1out=filtfilt(sos,glp,x1);
%disp(snr(x1out));
[POWER,F]=pwelch(x1out,[],[],[],Fs);%computes power density vs frequency power density=power/frequency
if abc==0
    figure()
plot(pwelch(x1out,[],[],[],Fs));
figure()
spectrogram(x1out);
end
p1=POWER;
pavg1=(max(p1)-min(p1))/2;
%pav=(max(p1)-max(p1(~(p1==max(p1)))))
%pavg1=max(p1)-6*pav;
%pavg=sum(p1(p1>=pavg1))/length(p1(p1>=pavg1));
final_fr=floor((F(p1>=pavg1)));
%final_fr=floor(F(p1>pavg1/1.45));
final_freq=final_fr(final_fr<700);
a=findforsonglist(final_freq);
end
%%
function J=findforsonglist(freq)
J=cell(11,1);
j=0;
y1="qwe";
y="qwe";
for i=1:length(freq)
    k=findifsinglevaluechord(freq(i),y1);
    if (~strcmp(k,y)&& ~strcmp(k,y1))
    y=k;
    j=j+1;
    J{j}=y;
    end
end
J=J(1:j);
%{
for l=1:length(J)
    if(J(l)==0)
        J=J(1:l-1);
    end
end
%}
return;
end
%%
function y=findifsinglevaluechord(x,k)
open_6th_string=[82.4,"E";87.3 "F";92.5 ,"F#";97.0,"G";103.8,"G#";];
open_5th_string=[109.0 ,"A";116.5,"A#"; 123.5, "B";130.8 ,"C";138.6 ,"C#";];
open_4th_string=[146.8 ,"D";155.6 ",D#";164.8 ,"E";174.6 ,"F";185.0 ,"F#";];
open_3rd_string=[195.0 ,"G" ;207.6 ,"G#";220.0, "A";233.1, "A#";];
open_2nd_string=[246.9 ,"B";261.6 ,"C" ;277.2 ,"C#";293.6, "D";311.6, "D#";];
open_1st_string=[329.6 ,"E";349.2,"F";370.0 ,"F#";394.6 ,"G";415.3, "G#";];
fifth_fret_on_1st_string=[440.0, "A";466.1 ,"A#";493.8, "B";522.2 ,"C";554.3 ,"C#";587.3 ,"D"; 622.2 ,"D#";659.2 ,"E";];
%twelweth_fret_on_1st_string=[];
%guitar_list=["C",16.35;"C#",17.32;"D",18.35;"Eb",19.45;"E",20.60;"F",21.83;"F#",22.12;"G",24.50;"G#",25.96;"A",27.50;"Bb",29.14;"B",30.87;];
for i=1:length(open_6th_string)
if(floor(str2double(open_6th_string(i,1)))==x)%||(floor(str2double(open_6th_string(i,1)))<=x+3&&floor(str2double(open_6th_string(i,1)))>=x-3))
    s1="6th_string";
    s2=open_6th_string(i,2);
    m=strcat(s1,s2);
    if(strcmp(m,k)==0);

    y=m;  %disp(m);
    end
    return;
end
end
for i=1:length(open_5th_string)
if(floor(str2double(open_5th_string(i,1)))==x)%||(floor(str2double(open_5th_string(i,1)))<=x+3&&floor(str2double(open_5th_string(i,1)))>=x-3))
    s1="5th_string";
    s2=open_5th_string(i,2);
    m=strcat(s1,s2);
     if(strcmp(m,k)==0)
        y=m;  %disp(m);
    end
    return;
end
end
for i=1:length(open_4th_string)
if(floor(str2double(open_4th_string(i,1)))==x)%||(floor(str2double(open_4th_string(i,1)))<=x+1&&floor(str2double(open_4th_string(i,1)))>=x-1))
    s1="4th_string";
    s2=open_4th_string(i,2);
    m=strcat(s1,s2);
     if(strcmp(m,k)==0)
    y=m;  %disp(m);
    end
    return;
end
end
for i=1:length(open_3rd_string)
if(floor(str2double(open_3rd_string(i,1)))==x)%||(floor(str2double(open_3rd_string(i,1)))<=x+3&&floor(str2double(open_3rd_string(i,1)))>=x-3))
    s1="3rd_string";
    s2=open_3rd_string(i,2);
    m=strcat(s1,s2);
     if(strcmp(m,k)==0)
    y=m;  %disp(m);
    end
    return;
end
end
for i=1:length(open_2nd_string)
if(floor(str2double(open_2nd_string(i,1)))==x)%||(floor(str2double(open_2nd_string(i,1)))<=x+3&&floor(str2double(open_2nd_string(i,1)))>=x-1))
    s1="2nd_string";
    s2=open_2nd_string(i,2);
    m=strcat(s1,s2);
    if(strcmp(m,k)==0)
    y=m;  %disp(m);
    end
    return;
end
end
for i=1:length(open_1st_string)
if(floor(str2double(open_1st_string(i,1)))==x)%||(floor(str2double(open_1st_string(i,1)))<=x+3&&floor(str2double(open_1st_string(i,1)))>=x-3))
    s1="1st_string";
    s2=open_1st_string(i,2);
    m=strcat(s1,s2);
     if(strcmp(m,k)==0)
    y=m;  %disp(m);
    end
    return;
end
end
for i=1:length(fifth_fret_on_1st_string)
if(floor(str2double(fifth_fret_on_1st_string(i,1)))==x)%||(floor(str2double(fifth_fret_on_1st_string(i,1)))<=x+3&&floor(str2double(fifth_fret_on_1st_string(i,1)))>=x-3))
    s1="1st_string";
    s2=fifth_fret_on_1st_string(i,2);
    m=strcat(s1,s2);
    y=m;  %disp(m);
    return;
end
end
%{
for i=1:length(twelweth_fret_on_1st_string)
if(floor(str2double(fifth_fret_on_1st_string(i,1)))==x)%||(floor(str2double(fifth_fret_on_1st_string(i,1)))<=x+3&&floor(str2double(fifth_fret_on_1st_string(i,1)))>=x-3))
    s1="1st_string";
    s2=twelweth_fret_on_1st_string(i,2);
    m=strcat(s1,s2);
     if(strcmp(m,k)==0)
    y=m;
      %disp(m);
    end
    return;
end
end

for i=1:length(guitar_list)
for j=1:8
if(floor(str2double(guitar_list(i,2))*j)==x)%||(floor(str2double(guitar_list(i,2))*j)>=x-3&&floor(str2double(guitar_list(i,2))*j)<=x+3))
    s1=num2str(j);
    s2=guitar_list(i,1);
    m=strcat(s1,s2);
    if(strcmp(m,k)==0)
    y=m;
        %disp(m);
    end
    return;
end
end

end
%}
y="qwe";
return;
end

%%
% F=150;
% t=0:1/F:1;F=150;
% t=0:1/F:1;
% f1=5;
% x=cos(2*pi*t*f1);%+complex(0,1)*sin(2*pi*t*f)+cos(2*pi*400*t);
% nfft=1024;
% X=fft(x,nfft);
% X=X(1:nfft/2);
% mx=abs(X);
% x=cos(2*pi*t*f1);
% nfft=1024;
% X=fft(x,nfft);
% %figure(1); spectrogram(x1,[],[],[],1000);
% [s,f,t,P]=spectrogram(x1,[],[],[],Fs);
% freq_in_khz=f/Fs;
% freq_in_khz(1:200);
%figure(2);

%%
function y = bpfilt(signal, f1, f2, fs, isplot)
% Bandpass filtering
%
% Syntax:
%   y = bpfilt(signal, f1, f2, [options])
%
% Description:
%   This function performs bandpass filtering of a time series 
%   with rectangle window.
%
% Input Arguments:
%   signal 	- a column vector of time series.
%   f1 		- the lower bound of frequencies (in Hz).
%   f2 		- the upper bound of frequencies (in Hz).
%
% Options:
%   fs      - the sampling frequency in Hz. Default is 1 Hz.
%   isplot  - whether to produce plots.
%
% Output Arguments:
%   y 		- the filtered time series.
%
% Examples:
%   fs = 100;
%   t  = 1:1/fs:10;
%   x  = sin(t);
%   y  = bpfilt(x,20,30);
% getting options
if nargin < 4 || isempty(fs)
	fs = 1;
end
if nargin < 5 || isempty(isplot)
	isplot = 1;
end
%define variables
if isrow(signal)
    signal = signal';
end
N  = length(signal);
dF = fs/N;
f  = (-fs/2:dF:fs/2-dF)';
%Band-Pass Filter:
if isempty(f1) || f1==-Inf
    BPF = (abs(f) < f2);
elseif isempty(f2) || f2==Inf
    BPF = (f1 < abs(f));
else
    BPF = ((f1 < abs(f)) & (abs(f) < f2));
end
%{
% Power spectrum of the band-pass filter
if isplot
    figure;
    plot(f,BPF);
    title(sprintf('Power spectrum of the band-pass filter in (%.3f, %.3f) Hz',f1,f2));
end
%}
%Power spectrum of the original signal
signal 	 = signal-mean(signal);
spektrum = fftshift(fft(signal))/N;
if (isplot == 1)
    figure;
    subplot(2,1,1);
    plot(f,abs(spektrum));
    title('Power spectrum of the original signal');
end
%Power spectrum of the band-pass filtered signal
%spektrum = BPF.*spektrum;
if (isplot == 1)
    subplot(2,1,2);
    plot(f,abs(spektrum));
    title(sprintf('Power spectrum of the band-pass filtered signal in (%.3f, %.3f) Hz',f1,f2));
end
%The band-pass filtered time series
y = ifft(ifftshift(spektrum)); %inverse ifft
y = real(y);
if (isplot==1)
    time = 1/fs*(0:N-1)';
	
    figure;
    subplot(2,1,1);
    plot(time,signal);
    title('The original time series');
    subplot(2,1,2);
    plot(time,y);
    title(sprintf('The band-pass filtered time series in (%.3f, %.3f) Hz',f1,f2));
end
end