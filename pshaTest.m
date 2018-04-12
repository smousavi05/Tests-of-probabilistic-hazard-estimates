clear all; close all;clc;

% %% Original data
% % RecordNumber --------------- 1 
% % Network -------------------- 2
% % Station Name --------------- 3
% % Station Channel ------------ 4
% % Station Lat ---------------- 5
% % Station Lon ---------------- 6
% % Station Year --------------- 7
% % Station Month -------------- 8
% % Station Day ---------------- 9
% % Station Hour --------------- 10
% % Station Min ---------------- 11
% % Station Sec ---------------- 12
% % EQ Latitude ---------------- 13
% % EQ Long -------------------- 14
% % EQ Depth M ----------------- 15
% % Magnitude ------------------ 16
% % MAg Type ------------------- 17
% % Distance km ---------------- 18
% % Filename ------------------- 19
% % PGA ------------------------ 20
% % 5 Hz (0.2 sec) ------------- 21
% % 1 Hz (1.0 sec) ------------- 22
% % PGV ------------------------ 23
% % Vs30 ----------------------- 24

PGAREF = [0.50000E-02;0.70000E-02;0.98000E-02;0.13700E-01;0.19200E-01;...
0.26900E-01;0.37600E-01;0.52700E-01;0.73800E-01;0.10300E+00;0.14500E+00;...
0.20300E+00;0.28400E+00;0.39700E+00;0.55600E+00;0.77800E+00;0.10900E+01;...
0.15200E+01;0.22000E+01;0.33000E+01];
% 
% % loading the instrumental observations.
% [~ , NtNm] = xlsread('2016cc2.xlsx','A','B1:B18865');
% [~ , StNm] = xlsread('2016cc2.xlsx','A','C1:C18865');
% [~ , StCh] = xlsread('2016cc2.xlsx','A','D1:D18865');
% A = xlsread('2016cc2.xlsx','A');
load('A')

delete('mag.txt')
fileID = fopen('mag.txt','w');
for k = 1 : length(A)
fprintf(fileID,'%2.3f %2.3f %2.3f\n',A(k,14), A(k,13), A(k,16));
end
fclose(fileID);

delete('site.txt')
fileID = fopen('site.txt','w');
for k = 1 : length(A)
fprintf(fileID,'%2.3f %2.3f %2.3f\n',A(k,6), A(k,5), A(k,24));
end
fclose(fileID);

% h =histogram(A(:,8)); axis tight; hold on
% h.FaceColor = [0 0.5 0.5];
% h.EdgeColor = 'k';
%  xlabel('Month'); ylabel('Number of Reported Ground Motions')
%  ax = gca; ax.TitleFontSizeMultiplier = 1.8; ax.LabelFontSizeMultiplier=1.8;
% ax.FontWeight='bold';hold off

% % excluding Oklahoma data
% id = []; B = []; NtNm2 = [];StNm2=[];StCh2=[];
% for pp = 1: length(A)
%     if A(pp,6) < -95 && A(pp,6) > -102 && A(pp,5) > 33 && A(pp,5) < 37.5 ...
%             && A(pp,14) < -95 && A(pp,14) > -102 && A(pp,13) > 33 && A(pp,13) < 37.5
%        id = [id;pp];
%    else
%        B = [B;A(pp,:)];
%        NtNm2 = [NtNm2;NtNm(pp,:)];
%        StNm2 = [StNm2;StNm(pp,:)];
%        StCh2 = [StCh2;StCh(pp,:)];
%    end
% end
% A = B;
% NtNm = NtNm2;
% StNm = StNm2;
% StCh = StCh2;

% % just Oklahoma data
% id = []; B = []; NtNm2 = [];StNm2=[];StCh2=[];
% for pp = 1: length(A)
%     if A(pp,6) < -95 && A(pp,6) > -102 && A(pp,5) > 33 && A(pp,5) < 37.5 ...
%             && A(pp,14) < -95 && A(pp,14) > -102 && A(pp,13) > 33 && A(pp,13) < 37.5
%        B = [B;A(pp,:)];
%        NtNm2 = [NtNm2;NtNm(pp,:)];
%        StNm2 = [StNm2;StNm(pp,:)];
%        StCh2 = [StCh2;StCh(pp,:)];
%    end
% end
% A = B;
% NtNm = NtNm2;
% StNm = StNm2;
% StCh = StCh2;

%% loading hazard curves from the 2016 model. 
[w1,w2,w3,w4,w5,w6,w7,w8,w9,w10,w11,w12,w13,w14,w15,w16,w17,w18,w19,w20,w21,w22] = ...
    textread('hazardCurve16.txt','%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f');
z5 = zeros(length(w1),1);
% dis, lat, lon, PGAREF(1) ... PGAREF(20)
forcast = [z5,w1,w2,w3,w4,w5,w6,w7,w8,w9,w10,w11,w12,w13,w14,w15,w16,w17,w18,w19,w20,w21,w22,z5,z5];
% removing western hazard curves.
forcast(find(forcast(:,3) < -115.000),:) = [];

%% loading hazard curves from the 2017 model. 
[w1,w2,w3,w4,w5,w6,w7,w8,w9,w10,w11,w12,w13,w14,w15,w16,w17,w18,w19,w20,w21,w22] = ...
    textread('hazardCurve17.txt','%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f');
z5 = zeros(length(w1),1);
% dis, lat, lon, PGAREF(1) ... PGAREF(20)
forcast17 = [z5,w1,w2,w3,w4,w5,w6,w7,w8,w9,w10,w11,w12,w13,w14,w15,w16,w17,w18,w19,w20,w21,w22,z5,z5];
% removing western hazard curves.
forcast17(find(forcast17(:,3) < -115.000),:) = [];

%% loading hazard curves from the 2014 model. 
[w1,w2,w3,w4,w5,w6,w7,w8,w9,w10,w11,w12,w13,w14,w15,w16,w17,w18,w19,w20,w21,w22] = ...
    textread('hazardCurve17.txt','%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f');
z5 = zeros(length(w1),1);
% dis, lat, lon, PGAREF(1) ... PGAREF(20)
forcast14 = [z5,w1,w2,w3,w4,w5,w6,w7,w8,w9,w10,w11,w12,w13,w14,w15,w16,w17,w18,w19,w20,w21,w22,z5,z5];
% removing western hazard curves.
forcast14(find(forcast14(:,3) < -115.000),:) = [];

%% loading informed model and branches noWeight. 
load('branchHazW')
load('informedNW')

%% loading the adaptive model.no weight 
[w1,w2,w3,w4,w5,w6,w7,w8,w9,w10,w11,w12,w13,w14,w15,w16,w17,w18,w19,w20,w21,w22] = ...
textread('./branch-adp/Llenos_max_NSH_CEUS.pga.txt','%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f');
adaptiveNW = [w1,w2,w3,w4,w5,w6,w7,w8,w9,w10,w11,w12,w13,w14,w15,...
           w16,w17,w18,w19,w20,w21,w22];
 
%Filling empty elemetns in station filed. 
IS = isnan(A(:,3));
ind = find(IS == 0);
b = strread(num2str(A(:,3)'),'%s'); % converting double to cell array.
for j =1:length(ind);
StNm(ind(j)) = b(ind(j));
end

% testing if any empty element still exist.
chN = cellfun('isempty',NtNm);
if sum(chN) > 0;
 disp('empty cell in network list')
end

chN = cellfun('isempty',StNm);
if sum(chN) > 0;
 disp('empty cell in station list')
end

chN = cellfun('isempty',StCh);
if sum(chN) > 0;
 disp('empty cell in channel list')
end

% making a unique location ID based on station long and lat.
stuq = unique(StNm);
cumPGA = []; cumDis=[]; cumMag =[]; cumT=[];

% clustering data into stations.
for i = 1: length(stuq) % loop over unique stations.
    i
    PGA=[]; Sp1=[]; Sp5=[]; t=[]; Ch=[]; Dis=[]; M=[];
    n = strmatch(stuq(i), StNm); % finding indices of each station in data.
    
    for ii = 1: length(n) % loop over number of recorded data at each station.
        time = sprintf('%4d-%02d-%02d %02d:%02d:%2f',A(n(ii),7),A(n(ii),8),A(n(ii),9), ...
             A(n(ii),10),A(n(ii),11),A(n(ii),12));
        time2 = datenum(time,'yyyy-mm-dd HH:MM:SS');
        t = [t;time2];
        PGA = [PGA ; A(n(ii),20)];
        Sp1 = [Sp1 ; A(n(ii),22)];
        Sp5 = [Sp5 ; A(n(ii),21)];
        Ch = [Ch ; StCh(n(ii))]; 
        Dis = [Dis ; A(n(ii),18)];
        M = [M ; A(n(ii),16)];
    end

% removing duplicated strong motion when there is both accelerometer and seismometer.
    if length(t) >= 2
    dt = diff(t);
    iddt = find(dt == 0.0000000000000);
    
    if length(iddt) >= 1
    for tt = 1:length(iddt)
          PGM = max(PGA(iddt(tt)),PGA(iddt(tt)+1));
          PGA(iddt(tt)) = PGM; PGA(iddt(tt)+1) = 0;
          Dis(iddt(tt)+1) = 0; M(iddt(tt)+1) = 0;
           
          if length(Sp1) == length(PGA)
          SP = max(Sp1(iddt(tt)),Sp1(iddt(tt)+1));
          Sp1(iddt(tt)) = SP; Sp1(iddt(tt)+1) = 0;
          end
          
          if length(Sp5) == length(PGA)
          SP = max(Sp5(iddt(tt)),Sp5(iddt(tt)+1));
          Sp5(iddt(tt)) = SP; Sp5(iddt(tt)+1) = 0;
          end            
    end 
    end
    end 
        cumPGA = [cumPGA;PGA]; cumDis = [cumDis;Dis]; cumMag = [cumMag;M];
        
    % getting the start and stop times of each station from IRIS
    st = irisFetch.Stations('channel',char(NtNm(n(1))),char(StNm(n(1))),'*',char(StCh(n(1))));
    StBd = []; StEd = [];   
    for dd = 1:length(st);
    Bd = datenum(st(dd).StartDate,'yyyy-mm-dd HH:MM:SS');
    Ed = datenum(st(dd).EndDate,'yyyy-mm-dd HH:MM:SS');
    StBd = [StBd; Bd]; StEd = [StEd; Ed];
    end
    
    StBd2 = min(StBd); StEd2 = max(StEd); 
    b16 = datenum('2016-01-01 00:00:00','yyyy-mm-dd HH:MM:SS');
    e16 = datenum('2016-12-31 00:00:00','yyyy-mm-dd HH:MM:SS');
    
    if StBd2 < b16
       StBd2 = b16;
    end
    
    if StEd2 > e16
       StEd2 = e16;
    end
    
    if length(StBd2) > 1
     disp(sprintf('%2d',i))
    end
    
    % pulling coordinates of the station.
    StLa = A(n(1),5); StLo = A(n(1),6);
    
    % finding associated hazard curve from 2016 model  
     forcast(:,1) = forcast(:,2) - StLa;
     IDD = find(forcast(:,1) < 0.1 & forcast(:,1) > -0.1);
     forcut1 = forcast(IDD,:);
     forcut1(:,1) = forcut1(:,3) - StLo;
     IDD = find(forcut1(:,1) < 0.1 & forcut1(:,1) > -0.1);
     forcut2 = forcut1(IDD,:);
     
     for ll = 1:length(IDD)
     [arclen,az] = distance(StLo,StLa,forcut2(ll,2),forcut2(ll,3));
     forcut2(ll,1) = deg2km(arclen); 
     end
     
     [M I] = min(forcut2(:,1)); 
     HazCrv = forcut2(I,:);
     HazCrv(:,2:4);
     
     I1 = find(forcast(:,2) ==  HazCrv(:,2));
     I2 = find(forcast(:,3) ==  HazCrv(:,3)); 
     I3 = intersect(I1,I2);
     
    % writting them into a structure for each station
    ss{i} = struct('NtNm',NtNm(n(1)),'StNm',StNm(n(1)),'StCh',StCh(n(1)),'StLo',...
        StLo(1),'StLa',StLa(1),'StBd',StBd2(1),'StEd',StEd2(1),'PGA', PGA,'Sp1', ...
        Sp1, 'Sp5', Sp5, 't', t , 'Dis',Dis,'Ch',Ch,'nRecorded',[],'nEstimated',[],...
        'Mag',M,'Forcast',HazCrv,'informed',informedNW(I3,:),'adaptive',adaptiveNW(I3,:),...
        'forcast17',forcast17(I3,:),'forcast14',forcast14(I3,:));
end




%% Checking completeness 
run GISMO/startup_GISMO.m;
goodSt = []; nE=[]; nR=[]; jbd=[]; SM =[];
for k =1:length(ss);
% % based on operational time of each station
% aa= ss{k}.StBd;
% bb= ss{k}.StEd;
% starttime = datestr(aa,'yyyy-mm-dd HH:MM:SS');
% endtime =  datestr(bb,'yyyy-mm-dd HH:MM:SS');
%% based on the first and last recordes
aa= ss{k}.t;
starttime = datestr(aa(1),'yyyy-mm-dd HH:MM:SS');
endtime =  datestr(aa(end),'yyyy-mm-dd HH:MM:SS');

Lat =(ss{k}.StLa);  % stations' latitude
Lon = (ss{k}.StLo); % stations' longitude

minmag = 3; % minimummagnitude
maxDepth = 100; % maximum depth of events

% finding number events occuring during operational time of each station. 
ev = irisFetch.Events('radialcoordinates' , [Lat, Lon, 2, 1]  ...
, 'starttime', starttime ,'endtime',endtime,'minimummagnitude',minmag,'maximumDepth',maxDepth);

GMPE =[];
for kk = 1: length(ev) 
    % calculating the distance between the events and the station 
    [arclen,az] = distance(Lat,Lon,ev(kk).PreferredLatitude,ev(kk).PreferredLongitude);
    diskm = deg2km(arclen); mag = ev(kk).Magnitudes.Value;
    
    % calculating Joiner-Boor distance. 
    Rjb = strike_distance(mag,diskm); jbd = [jbd; Rjb];   
    %Estimating the ground motion at the site from each event using
    [Sa, sigma, period1] = SP16SCALED_CEUS_GMPE(mag, Rjb, 0); % Ali's GMPE. 
%   [Sa, sigma] = A_2015_small_M(mag, 0, Rjb); % Attckinnson's GMPE. 
    if Sa > 1.0e-03
      GMPE = [GMPE; Sa];
    end
end

CC = ss{k}.PGA;
% finding number of non-repeated records
% CCL = nnz(CC);
[ind v] = find(CC > 1.0e-03); CCL =sum(v);

if CCL >= length(GMPE);
   goodSt = [goodSt;k]; 
end
 nR = [nR;CCL]; nE = [nE;length(GMPE)];
end


%% plotting number of recordeds at each station vs estimated numbers
figure
 k3 =1:1:length(ss);
   for ii =1 : length(k3')
       if nE(ii) < nR(ii) 
        difn = (nR(ii)-nE(ii))/3; nR(ii)=nR(ii)-difn; nE(ii)=nE(ii)+difn;
        line([k3(ii)' k3(ii)'],[nE(ii) nR(ii)]); hold on
       else
        line([k3(ii)' k3(ii)'],[ nR(ii) nE(ii)]); hold on
      end
 end 
 h1 = plot(k3',nE, 'or','MarkerFaceColor','r')
 xlabel('Stations Id'); ylabel('Number of Strong Motions')
 hold on; h2 = plot(k3',nR, 'ob')
 lgd = legend([h1 h2],{'Predicted','Reported'},'Location','northeast','Orientation','vertical')
 lgd.FontSize = 12; lgd.TextColor = 'blue';1
 ax = gca; ax.TitleFontSizeMultiplier = 1.8; ax.LabelFontSizeMultiplier=1.8;
ax.FontWeight='bold';hold off

%% Plotting strong motion data and GMPE
gmp =[]; gmps=[]; disgm = []; gmpA =[]; gmpsA=[];
for dd = 1:1:200  
% finding average magnitude for distance range
%   disInd = find(dd-0.05 <= cumDis <= dd+0.05);
    sumM1=0; nM=0; 
    for dd2 = 1 : length(cumDis)
    if cumMag(dd2) > 0 && dd-20<= cumDis(dd2) && cumDis(dd2) <= dd+20;
        sumM1 = cumMag(dd2)+ sumM1; nM = nM+1;
    end
    end
    medM = sumM1/nM;  
    [Sa, sigma, period1] = SP16SCALED_CEUS_GMPE(medM, dd, 0);
    [SaA, sigmaA] = A_2015_small_M(medM, 0, dd); % Attckinnson's GMPE. 
    disgm = [disgm;dd]; gmp = [gmp;Sa]; gmps = [gmps;sigma];
    gmpA = [gmpA;SaA]; gmpsA = [gmpsA;sigmaA];
end

mPGA =[]; ddd =[];
for dd = 1:5:200
    % finding average magnitude for distance range
    sumPGA=0; nP=0;
    for dd2 = 1 : length(cumDis)
    if cumMag(dd2) > 0 && dd-5<= cumDis(dd2) && cumDis(dd2) <= dd+5;
        sumPGA = cumPGA(dd2)+ sumPGA; nP = nP+1;
    end
    end
    medPGA = sumPGA/nP; mPGA = [mPGA;medPGA]; ddd = [ddd;dd];
end

figure
plot(cumDis,cumPGA,'ob'); axis tight; hold on
h1 = plot(disgm,gmp,'--c','LineWidth',3);hold on
h2 = plot(disgm,gmpA,'--m','LineWidth',3);hold on
set(gca,'LineWidth',3); set(gca, 'YScale', 'log')
ylabel('PGA (%g)');xlabel('Distance (km)');
grid on; ax = gca;
ax.TitleFontSizeMultiplier = 1.8; ax.LabelFontSizeMultiplier=1.8;
ax.FontWeight='bold';hold off
lgd = legend([h1 h2],{'Shahjouei & Pezeshk (2016)','Atkinson (2015)'},...
    'Location','northeast','Orientation','vertical');
lgd.FontSize = 12; lgd.TextColor = 'blue';
ax = gca; ax.TitleFontSizeMultiplier = 1.8;
ax.LabelFontSizeMultiplier=1.8; ax.FontWeight='bold';

%% plotting the eq and stations on the map
figure
subplot 121
h1 = borders('continental us','nomap','facecolor',[0.8 0.85 1],...
    'edgecolor',[0.7 0.7 1],'linewidth',1);axis tight; hold on
h2 = plot(A(:,6), A(:,5), 'd','MarkerSize',5,'MarkerEdgeColor',...
    'b','MarkerFaceColor','c','DisplayName','Stations')
lgd = legend([h2],{'Stations'},'Location','southeast','Orientation','vertical')
lgd.FontSize = 12; lgd.TextColor = 'blue';
ax = gca; ax.TitleFontSizeMultiplier = 1.8;
ax.LabelFontSizeMultiplier=1.8; ax.FontWeight='bold';


fileID = fopen('station_2016.txt','w');
for k = 1 : length(A)
fprintf(fileID,'%3.3f %2.3f\n',A(k,6), A(k,5));
end
fclose(fileID);

subplot 122
h1 = borders('continental us','nomap','facecolor',[0.8 0.85 1],...
    'edgecolor',[0.7 0.7 1],'linewidth',1); axis tight; hold on
h3 = plot(A(:,14), A(:,13), 'o','MarkerSize',7,'MarkerEdgeColor',...
    'red','MarkerFaceColor',[1 .6 .6],'DisplayName','Earthquakes')
lgd = legend([h3],{'Earthquakes'},'Location','southeast','Orientation','vertical')
lgd.FontSize = 12; lgd.TextColor = 'blue'; ax = gca;
ax.TitleFontSizeMultiplier = 1.8; ax.LabelFontSizeMultiplier=1.8;
ax.FontWeight='bold';

fileID = fopen('eq_2016.txt','w');
for k = 1 : length(A)
fprintf(fileID,'%3.3f %2.3f\n',A(k,14), A(k,13));
end
fclose(fileID);

%% plotting the mag vs distance
figure
h2 = plot(A(:,18), A(:,16), 'ob','MarkerSize',6.5,'MarkerEdgeColor',...
    'b','MarkerFaceColor','c','DisplayName','Stations')
hold on; xlim([0 200])
title('All Stations');
xlabel('Hypocentral Distance (km)'); ylabel('Magnitude');
set(gca, 'XScale', 'log')
ax = gca; ax.TitleFontSizeMultiplier = 1.8;
ax.LabelFontSizeMultiplier=1.8;ax.FontWeight='bold';


% % % using usamap
% % latlim = [25 55];
% % lonlim = [-110 -60];
% % figure 
% % ax = usamap(latlim, lonlim);
% % 
% % axis off
% % getm(gca,'MapProjection')
% % states = shaperead('usastatehi',...
% %     'UseGeoCoords',true,'BoundingBox',[lonlim',latlim']);
% % faceColors = makesymbolspec('Polygon',...
% %     {'INDEX',[1 numel(states)],'facecolor',[0.8 0.85 1],});
% % geoshow(ax,states,'SymbolSpec',faceColors)
% % 

%% Monte Carlo sampling of indipendent stations
 abcLo=[]; abcLa=[]; ID=[]; diss=[]; SM=[]; mg=[]; z4 =zeros(length(goodSt),1);rs=0;
 for ii =1:length(goodSt)
     ii
    ID = [ID; goodSt(ii)]; 
    Lo1 = ss{goodSt(ii)}.StLo;abcLo =[abcLo; Lo1]; 
    La1 = ss{goodSt(ii)}.StLa;abcLa =[abcLa; La1];
    diis = ss{goodSt(ii)}.Dis; R = strike_distance(M,diis(1));diss=[diss;R];
    smm = ss{goodSt(ii)}.PGA; SM =[SM;smm];
    M = ss{goodSt(ii)}.Mag; mg =[mg;M];
 end
    loc = [abcLo, abcLa, ID,z4];
    gmped = [diss, mg];
        
    %     [gm, sigma, period1]= SP16SCALED_CEUS_GMPE(M,R,0)
    [gm, sigma] = A_2015_small_M(M, 0, R);gmp =[gmp;gm];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plotting the compelet stations on the map
figure
h1 = borders('continental us','nomap','facecolor',[0.8 0.85 1],...
    'edgecolor',[0.7 0.7 1],'linewidth',1);
axis tight; hold on

h2 = plot(abcLo(:),abcLa(:), 'd','MarkerSize',5,'MarkerEdgeColor',...
    'b','MarkerFaceColor','c','DisplayName','Stations');
title('Complete Stations');

lgd = legend([h2],{'Stations'},'Location','southeast','Orientation','vertical');
lgd.FontSize = 12; lgd.TextColor = 'blue';    

%% loop over a number of monte carlo simulations
    xmin = 45; % minimum interstation distance between independent stations
for mm = 1:500; 
    x2 = loc;x3 = x2;
    
    sels =[];
    while length(x3(:,1)) >= 1  
          if length(x3(:,1)) >= 2
             y = datasample(x3,1,'Replace',false);
          else
             y = x3;
          end
         sels = [sels;y(3)];
    
    diskm=[];
    for nn = 1:length(x2(:,1));
        [arclen,az] = distance(x2(nn,2),x2(nn,1),y(2),y(1));
        diskm = [diskm; deg2km(arclen)];
    end
    x2(:,4) = diskm; TF = x2(:,4) < xmin; x2(TF,:) = []; 
    x2=sortrows(x2,4); x3 = x2(1:floor(0.3*length(x2(:,1))),:);
    end
    
    sel =[];
    for ls=1:length(sels)
    TF = loc(find(loc(:,3) == sels(ls)),:); sel = [sel;TF];
    end
com{mm} = [sel;x2];
end    


%% plotting independent stations for each Monte Carlo sample
figure
subplot 121
mmn = 41;
h1 = borders('continental us','nomap','facecolor',[0.8 0.85 1],...
    'edgecolor',[0.7 0.7 1],'linewidth',1);
axis tight; hold on

h2 = plot(com{mmn}(:,1),com{mmn}(:,2), 'd','MarkerSize',5,...
    'MarkerEdgeColor','b','MarkerFaceColor','c',...
    'DisplayName','Stations'); hold on
lgd = legend([h2],{'Independent Stations'},'Location','southeast','Orientation','vertical')
xlim([-108 -67]);ylim([24 43])
ax = gca; ax.TitleFontSizeMultiplier = 1.8;
ax.LabelFontSizeMultiplier=2.1;ax.FontWeight='bold'; 
subplot 122
mmn = 10;
h1 = borders('continental us','nomap','facecolor',[0.8 0.85 1],...
    'edgecolor',[0.7 0.7 1],'linewidth',1);
axis tight; hold on

h2 = plot(com{mmn}(:,1),com{mmn}(:,2), 'd','MarkerSize',5,...
    'MarkerEdgeColor','b','MarkerFaceColor','c',...
    'DisplayName','Stations'); hold on
xlim([-108 -67]);ylim([24 43])
lgd = legend([h2],{'Independent Stations'},'Location','southeast','Orientation','vertical')
ax = gca; ax.TitleFontSizeMultiplier = 1.8;
ax.LabelFontSizeMultiplier=2.1;ax.FontWeight='bold';    

fileID = fopen('indST_2016_10.txt','w');
for k = 1 : length(com{mmn})
fprintf(fileID,'%3.3f %2.3f\n',com{mmn}(k,1), com{mmn}(k,2));
end
fclose(fileID);

% for mmn=1:500;  
%      mmn
% % Observed Hazard Curve
% out =[];ww = [];
% for pga=1:length(PGAREF);
% oo =[];
% for i=1:length(com{mmn}(:,3));
%     stID = com{mmn}(i,3);
%     pgaOBS = ss{stID}.PGA;
%     hazfrc = ss{stID}.Forcast;
%     for ll =1:length(pgaOBS)
%     if  pgaOBS(ll) > PGAREF(pga); 
%         o=1;
%     else
%         o=0;
%     end
%        oo= [oo;o];
%        ww = [ww;hazfrc(pga+3)];
%     end
% end
% hobs = sum(oo);
% out=[out;PGAREF(pga),hobs,hazfrc(pga+3)];
% end
% mchzobs{mmn} = out;  
%  end
%  
% pgaa= PGAREF(1)
%  out2=[];
% for uu = 1:400
%     oo= mchzobs{uu};
%     IP = find(oo(:,1) == pgaa);
%     out2=[out2;uu,oo(IP,2),oo(IP,3)];
% end
% 
% figure
% h4 = plot(out2(:,1),out2(:,2),'-b','LineWidth',2);hold on 
% h5 = plot(out2(:,1),out2(:,3),'-r','LineWidth',2)
% xlabel('Monte Carlo Sample');
% ylabel('Number of Ground Motion Exceedances');
% title(sprintf('PGA:%0.3f',pgaa))
% ax = gca; ax.TitleFontSizeMultiplier = 1.8;
% ax.LabelFontSizeMultiplier=1.8; ax.FontWeight='bold';
% lgd = legend([h4 h5],{'Observation','Forcast'},'Location','northeast','Orientation','vertical')
% lgd.FontSize = 12;lgd.TextColor = 'blue';



% % to calculate the Fig 7 %mmn should be 20
REF = [0.50000E-02;0.70000E-02;0.98000E-02;0.13700E-01;0.19200E-01;...
0.26900E-01;0.37600E-01;0.52700E-01;0.73800E-01;0.10300E+00;0.14500E+00;...
0.20300E+00;0.28400E+00;0.39700E+00;0.55600E+00;0.77800E+00;0.10900E+01;...
0.15200E+01;0.22000E+01;0.33000E+01];

% % % to calculate the Fig 8
% REF = [0.005; 0.01; 0.02; 0.04; 0.1; 0.2]; %mmn should be 500

% calculating observed Hazard curve based on number of stations with
% exceedance    
for mmn = 1:20; % loop over number of Monte Csarlo simulations Observed Hazard Curve
mmn
out =[];
for pga = 1:length(REF); %loop over number of PGAs
      mco = [];  
      oo = []; ww = []; ww17 = []; ww14 = []; wwAD = []; wwIF = [];
      for i = 1:length(com{mmn}(:,3)); % loop over number of stations
               stID = com{mmn}(i,3);
               pgaOBS = ss{stID}.PGA;
               hazfrc = ss{stID}.Forcast;
%                hazfrc17 = ss{stID}.forcast17;
%                hazfrc14 = ss{stID}.forcast14;
               hazfrcAD = ss{stID}.adaptive;
               hazfrcIF = ss{stID}.informed;
               stt = ss{stID}.StBd;
               edt = ss{stID}.StEd;
               tp = months(stt,edt);
               tp = tp./12;
               
        if   length(REF) > 15   
            PG = REF(pga);
            lambda = (hazfrc(pga+3))*tp; 
            lambdapAD = (hazfrcAD(pga+2))*tp;
            lambdapIF = (hazfrcIF(pga+2))*tp;
            
        elseif pga == 1 
            PG = REF(pga);
            lambda = (hazfrc(pga+3))*tp; 
            lambdapAD = (hazfrcAD(pga+2))*tp;
            lambdapIF = (hazfrcIF(pga+2))*tp;
        elseif pga == 2
            PG = REF(pga);            
                  y = [PGAREF(3);PGAREF(4)];      
                  x = [hazfrc(6);hazfrc(7)];
                  yi = interp1(y,x,PG,'linear'); lambda = (yi)*tp; 
                  x = [hazfrcAD(5);hazfrcAD(6)];
                  yi = interp1(y,x,PG,'linear'); lambdapAD = (yi)*tp;      
                  x = [hazfrcIF(5);hazfrcIF(6)];
                  yi = interp1(y,x,PG,'linear'); lambdapIF = (yi)*tp;
         elseif pga == 3
            PG = REF(pga);            
                  y = [PGAREF(5);PGAREF(6)];          
                  x = [hazfrc(8);hazfrc(9)];
                  yi = interp1(y,x,PG,'linear'); lambda = (yi)*tp;                
                  x = [hazfrcAD(7);hazfrcAD(8)];
                  yi = interp1(y,x,PG,'linear'); lambdapAD = (yi)*tp;                
                  x = [hazfrcIF(7);hazfrcIF(8)];
                  yi = interp1(y,x,PG,'linear'); lambdapIF = (yi)*tp; 
         elseif pga == 4
            PG = REF(pga);            
                  y = [PGAREF(7);PGAREF(8)];          
                  x = [hazfrc(10);hazfrc(11)];
                  yi = interp1(y,x,PG,'linear'); lambda = (yi)*tp;                
                  x = [hazfrcAD(9);hazfrcAD(10)];
                  yi = interp1(y,x,PG,'linear'); lambdapAD = (yi)*tp;                
                  x = [hazfrcIF(9);hazfrcIF(10)];
                  yi = interp1(y,x,PG,'linear'); lambdapIF = (yi)*tp; 
         elseif pga == 5
            PG = REF(pga);            
                  y = [PGAREF(9);PGAREF(10)];          
                  x = [hazfrc(12);hazfrc(13)];
                  yi = interp1(y,x,PG,'linear'); lambda = (yi)*tp;                
                  x = [hazfrcAD(11);hazfrcAD(12)];
                  yi = interp1(y,x,PG,'linear'); lambdapAD = (yi)*tp;                
                  x = [hazfrcIF(11);hazfrcIF(12)];
                  yi = interp1(y,x,PG,'linear'); lambdapIF = (yi)*tp; 
         elseif pga == 6
            PG = REF(pga);            
                  y = [PGAREF(11);PGAREF(12)];          
                  x = [hazfrc(14);hazfrc(15)];
                  yi = interp1(y,x,PG,'linear'); lambda = (yi)*tp;                
                  x = [hazfrcAD(13);hazfrcAD(14)];
                  yi = interp1(y,x,PG,'linear'); lambdapAD = (yi)*tp;                
                  x = [hazfrcIF(13);hazfrcIF(14)];
                  yi = interp1(y,x,PG,'linear'); lambdapIF = (yi)*tp;
        end
        
                % Observed hazard
                 cnd = pgaOBS >= PG;
                   if sum(cnd) > 0  
                       o = 1;
                      elseif sum(cnd) == 0
                       o = 0;
                   end
                   oo = [oo;o];
                   
                % Forcasted Hazard       
                x = [0 1 2 3 4 5]; xq = x;
                ypdf = pdf('Poisson',x,lambda); % generatinf the PDF 
                cdf = cumsum(ypdf); %%% P(x) %%% intrgratinf the PDF to make CDF 
                [cdf, mask] = unique(cdf); % remove non-unique elements
                cdf(isnan(cdf)) = 0;
                if length(cdf) > 1 && sum(cdf) > 1 
                xq = xq(mask);
                randomValues = rand(1,5000);% create an array of random numbers
                % inverse interpolation to achieve P(x) -> x projection of the random values
                projection = interp1(cdf, xq, randomValues);
                projection(isnan(projection)) = 0; % removing NAN elements
                w = round(projection);
                else
                w = zeros(1,5000);
                end
                ww = [ww;w]; 
%              lambda17 = (hazfrc17(pga+3))*tp;
%                 x = [0 1 2 3 4 5];xq = x;
%                 ypdf = pdf('Poisson',x,lambda17); % generatinf the PDF 
%                 cdf = cumsum(ypdf); %%% P(x) %%% intrgratinf the PDF to make CDF 
%                 [cdf, mask] = unique(cdf); % remove non-unique elements
%                 if length(cdf) > 1
%                 xq = xq(mask);
%                 randomValues = rand(1,5000); % create an array of random numbers
%                 % inverse interpolation to achieve P(x) -> x projection of the random values
%                 projection = interp1(cdf, xq, randomValues);
%                 projection(isnan(projection)) = 0;
%                 w17 = round(projection); 
%                 else
%                  w17 = zeros(1,5000);
%                 end
%                 ww17 = [ww17;w17];

%                 x = [0 1 2 3 4 5]; xq = x;
%                 ypdf = pdf('Poisson',x,lambda14); % generatinf the PDF 
%                 cdf = cumsum(ypdf); %%% P(x) %%% intrgratinf the PDF to make CDF 
%                 % remove non-unique elements
%                 [cdf, mask] = unique(cdf); xq = xq(mask);
%                 if length(cdf) > 1
%                 % create an array of random numbers
%                 randomValues = rand(1,5000);
%                 % inverse interpolation to achieve P(x) -> x projection of the random values
%                 projection = interp1(cdf, xq, randomValues);
%                 projection(isnan(projection)) = 0;
%                 w14 = round(projection); 
%                 else
%                   w14 = zeros(1,5000);
%                 end
%                 ww14 = [ww14;w14];

                x = [0 1 2 3 4 5]; xq = x;
                ypdf = pdf('Poisson',x,lambdapAD); % generatinf the PDF 
                cdf = cumsum(ypdf); %%% P(x) %%% intrgratinf the PDF to make CDF 
                [cdf, mask] = unique(cdf);   % remove non-unique elements
                cdf(isnan(cdf)) = 0;
                if length(cdf) > 1 && sum(cdf) > 1 
                xq = xq(mask);
                % create an array of random numbers
                randomValues = rand(1,5000);
                % inverse interpolation to achieve P(x) -> x projection of the random values
                projection = interp1(cdf, xq, randomValues);
                projection(isnan(projection)) = 0;
                wAD = round(projection);   
                else
                    wAD = zeros(1,5000);
                end
                wwAD = [wwAD;wAD];      

                x = [0 1 2 3 4 5]; xq = x;
                ypdf = pdf('Poisson',x,lambdapIF); % generatinf the PDF 
                cdf = cumsum(ypdf); %%% P(x) %%% intrgratinf the PDF to make CDF 
                [cdf, mask] = unique(cdf);% remove non-unique elements
                cdf(isnan(cdf)) = 0;
                if length(cdf) > 1 && sum(cdf) > 1 
                xq = xq(mask);
                % create an array of random numbers
                randomValues = rand(1,5000);
                % inverse interpolation to achieve P(x) -> x projection of the random values
                projection = interp1(cdf, xq, randomValues);
                projection(isnan(projection)) = 0;
                wIF = round(projection);  
                else
                    wIF = zeros(1,5000);
                end
                   wwIF = [wwIF;wIF];
     end
        mcf = []; mcf17 = []; mcf14 = []; mcfAD = []; mcfIF = [];
        hobs = sum(oo); 
        for kk = 1:length(ww(1,:))
        s = nonzeros(ww(:,kk)); mcf = [mcf,length(s)]; % 2016 forcast
        
%         s = nonzeros(ww17(:,kk));% 2017 forcast
%         mcf17 = [mcf17;length(s)];
%            mcf1705 = quantile(mcf17,0.05);
%            mcf1795 = quantile(mcf17,0.95);
%            mcf1750 = quantile(mcf17,0.50);
%         s = nonzeros(ww14(:,kk));
%         mcf14 = [mcf14;length(s)];% 2014 forcast
%            mcf1405 = quantile(mcf14,0.05);
%            mcf1495 = quantile(mcf14,0.95);
%            mcf1450 = quantile(mcf14,0.50);
        s = nonzeros(wwAD(:,kk)); mcfAD = [mcfAD;length(s)];% 2016 Adoptive forcast
        s = nonzeros(wwIF(:,kk));mcfIF = [mcfIF;length(s)];% 2016 Informed forcast

        end  
   mcf05 = quantile(mcf,0.05);mcf95 = quantile(mcf,0.95);mcf50 = quantile(mcf,0.50);
   mcfAD05 = quantile(mcfAD,0.05);mcfAD95 = quantile(mcfAD,0.95);mcfAD50 = quantile(mcfAD,0.50);        
   mcfIF05 = quantile(mcfIF,0.05);mcfIF95 = quantile(mcfIF,0.95); mcfIF50 = quantile(mcfIF,0.50);     
   out = [out;PGAREF(pga),mcf05,mcf95,mcf50,...
          mcfAD05,mcfAD95,mcfAD50,mcfIF05,mcfIF95,mcfIF50,hobs];
end
   out1{mmn} = out;
end

 ll= 1    
%% plotting the all PGAs
out = out1{ll};
figure
hold on
h4 = jbfill(out(:,1)',out(:,6)',out(:,5)','b',rand(1,3),0,0.2)
plot(out(:,1),out(:,5),'-.b','LineWidth',2.5);%Adaptive 2016
plot(out(:,1),out(:,6),'-.b','LineWidth',2.5);%Adaptive 2016
hold on 
h5 = jbfill(out(:,1)',out(:,9)',out(:,8)','m',rand(1,3),0,0.2)
plot(out(:,1),out(:,8),':m','LineWidth',2.5);%Informed 2016
plot(out(:,1),out(:,9),':m','LineWidth',2.5);%Informed 2016
hold on 
h1 = jbfill(out(:,1)',out(:,3)',out(:,2)','k',rand(1,3),0,0.2)
plot(out(:,1),out(:,2),'--k','LineWidth',2.5);
plot(out(:,1),out(:,3),'--k','LineWidth',2.5);%Forcast 2016
hold on
h6 = plot(out(:,1),out(:,11),'s','LineWidth',1.5,'MarkerSize',10,...
    'MarkerEdgeColor','k','MarkerFaceColor','r');
set(gca, 'XScale', 'log')
xlabel('PGA(g)');
ylabel({'Number of Stations';'with Exceedance'});
%  set(gca,'LineWidth',2)
grid on; set(gca,'fontsize',21)
xlim([0 1])
ylim([0 20])
ax = gca;
ax.FontWeight='bold';

% lgd = legend([h1 h4 h5 h6],{'5%-95% Forcast 2016 Interval', ...
%        '5%-95% Adaptive Model 2016','5%-95% Informed Model 2016',...
%     'Observation'},'Location','northeast','Orientation','vertical')
% lgd.FontSize = 20;
% lgd.TextColor = 'k';

% % extacting exceedance numbers for a specific PGA level and plotting them for 
% % all monte carlo samples.  
% pga = 1
% % close all;
% figure
% for ll = 1:length(out1)
% out2 = out1{ll};
% h1 = errorbar(ll,out2(pga,4),out2(pga,2),out2(pga,3),'b','LineWidth',2)
% hold on
% h4 = plot(ll,out2(pga,11),'sr','LineWidth',2,'MarkerSize',10,...
%     'MarkerEdgeColor','r','MarkerFaceColor','m');%'Observation'
% hold on 
% end
% ref = REF(pga);
% 
% ylabel({'Number of Stations';'with Exceedance'});
% title(sprintf('PGA: %0.3f (g)',ref))
% set(gca,'fontsize',21)
% ax = gca;
% ax.FontWeight='bold';
% xlabel('Monte Carlo Samples');
% 
% ylim([0 5])
% 
% lgd = legend([h4 h1],{'Observation','5%-95% Forecast 2016 Interval'})
% lgd.FontSize = 20;lgd.TextColor = 'k';

%  histogram(mcf,'Normalization','pdf','FaceColor',[0 .5 .5],'EdgeColor','k','LineWidth',1.5)
%  mcf05 = quantile(mcf,0.05);
%  mcf95 = quantile(mcf,0.95);
%  hold on; 
%  plot(mcf05,0,'o','Color','k')
%  plot(mcf95,0,'o','Color','k')
%  ylabel('Prob'); xlabel('Nsites');
% ax = gca;ax.TitleFontSizeMultiplier = 2.1; ax.LabelFontSizeMultiplier=1.8;
% ax.FontWeight='bold';


% x = [0 1 2 3 4 5];
% ypdf = pdf('Poisson',x,0.57);cdf = cumsum(ypdf);
% bar(x,ypdf,'FaceColor',[0 .5 .5],'EdgeColor','k','LineWidth',1.5)
% ylabel('Prob');
% ax = gca;ax.TitleFontSizeMultiplier = 2.1; ax.LabelFontSizeMultiplier=1.8;
% ax.FontWeight='bold';
% % cdfplot(x)

% for k =1:length(ss);
%    if length(ss{k}) > 700
%        k
%    end
% end


% 
% % Observed Hazard Curve
% out =[];ww = [];
% pgaOBS = ss{167}.PGA;
% hazfrc = ss{167}.Forcast;
% 
% for pga=1:length(PGAREF);
% oo =[];
%     for ll =1:length(pgaOBS)
%     if  pgaOBS(ll) > PGAREF(pga); 
%         o=1;
%     else
%         o=0;
%     end
%        oo= [oo;o];
%        ww = [ww;hazfrc(pga+3)];
%     end
% hobs = sum(oo);
% out=[out;PGAREF(pga),hobs,hazfrc(pga+3)];
% end
% 
%  
% pgaa= PGAREF(1)
%  out2=[];
% for uu = 1:400
%     oo= mchzobs{uu};
%     IP = find(oo(:,1) == pgaa);
%     out2=[out2;uu,oo(IP,2),oo(IP,3)];
% end
% 
% 
% figure
% hold on
% h4 = plot(out(:,1),out(:,2),'-b','LineWidth',2)
% hold on 
% h5 = plot(out(:,1),out(:,3),'-r','LineWidth',2)
% xlabel('PGA');
% ylabel('Number of Ground Motion Exceedances');
% set(gca, 'YScale', 'log')
% set(gca, 'XScale', 'log')
% ax = gca;
% ax.TitleFontSizeMultiplier = 1.8;
% ax.LabelFontSizeMultiplier=1.8;
% ax.FontWeight='bold';
% lgd = legend([h4 h5],{'Observation','Forcast'},'Location','northeast','Orientation','vertical')
% lgd.FontSize = 12;
% lgd.TextColor = 'blue';