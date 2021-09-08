
count1_array=[];
count2_array=[];
for abc=1:1
filename = 'C:\Users\vmcx3\Desktop\Data set\Electricity data\TOSG_cont\2016.csv';
delimiter = ',';
startRow = 2;

formatSpec = '%C%f%f%f%f%[^\n\r]';

fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string', 'EmptyValue', NaN, 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');

fclose(fileID);

%% Post processing for unimportable data.
% No unimportable data rules were applied during the import, so no post
% processing code is included. To generate code which works for
% unimportable data, select unimportable cells in a file and regenerate the
% script.

%% Create output variable
File_data = table(dataArray{1:end-1}, 'VariableNames', {'localminute','dataid','use','gen','grid'});

%% Clear temporary variables
clearvars filename delimiter startRow formatSpec fileID dataArray ans;
% File_data = csvread('2014.csv',1,0,[1,0,1511481,3]);

I1 = unique(File_data(:,2));
T1= unique(File_data(:,1));
I=table2array(I1);
T=table2array(T1);
[height, width] = size(File_data);
File_dataT1 = File_data(:,1);
File_dataT2 = File_data(:,2);
File_dataT3 = File_data(:,3);
File_dataT4 = File_data(:,4);
File_dataT5 = File_data(:,5);

File_dataP1 = table2array(File_dataT1);
File_dataS1 = cellstr(File_dataP1);
File_data1 = datetime(File_dataS1,'InputFormat','yyyy-MM-dd HH:mm:ss');
File_data2 = table2array(File_dataT2);
File_data3 = table2array(File_dataT3);
File_data4 = table2array(File_dataT4);
File_data5 = table2array(File_dataT5);


N = length(I);
M=ceil(3*N/10);
SW = 100;                                                                   % SW is the species width in watts.
R = 0;                                                                      % Initialising richness (Number of unique species) to 0.
q = 3.5;
count = 0;
Meter_Readings = abs(File_data3(:)*1000);                                       

for i=1:size(Meter_Readings)
    if Meter_Readings(i) > 3000
        count = count + 1;
        Meter_Readings(i) = 3000;
    end
    if Meter_Readings(i) == 0
        Meter_Readings(i) = 1;
    end
end
Maximum_Reading = max(Meter_Readings);                                      % Retrieves maximum reading to find richness.
Minimum_Reading = min(Meter_Readings);                                      % Retrieves minimum reading to find richness.
R = ceil((Maximum_Reading-Minimum_Reading)/SW);                             % Calculates the value R (Number of unique species)

% max_day = 302;
% min_day = min(floor(File_data(:,2)/100));
frame_size = 15;
% max_day = max_day - rem((max_day - min_day),frame_size) - 1;
number_of_frames = ceil(datenum(File_data1(end)-File_data1(1))/frame_size);

P_matrix = zeros(10, N, R);                                                         % Declaring matrix to store probability of R species of N meters
for z=1:size(File_data3)/4
    frame_number = ceil((datenum(File_data1(z) - File_data1(1)) + 1) / frame_size);
            meter_number = find(I==File_data2(z));
            species_number = ceil(abs(Meter_Readings(z))/SW);
            P_matrix(frame_number,meter_number,species_number)= P_matrix(frame_number,meter_number,species_number) + 1;
end


disp(1);


P = zeros(10, N, R);                                                         % Declaring matrix to store probability of R species of N meters
for x=1:1
for y=1:N
for z=1:R
P(x,y,z) = (P_matrix(x,y,z)+1)/(sum(P_matrix(x,y,:))+R);
end
end
end
% Z = zeros(N, number_of_frames, R, R);
% for v = 1:N
%     for w = 1:number_of_frames
%         for x = 1:R
%             for y=1:R
%                 Z(v,w,y,x) = P(w,v,y);
%             end
%         end
%     end
% end
Z1_final{N} = [];
for w=1:N
for x=1:R
Z1_final{w}(x) =  0;
for y=1:1
Z1_final{w}(x) = Z1_final{w}(x) + abs(P(y,w,x)-P(y,w,x));
end
Z1_final{w}(x) = Z1_final{w}(x)/number_of_frames;
end
end
% P_final{N}=[];
% for w = 1:N
%     for x=1:R
%         P_final{w}(x,1) = (P(1,w,x)+P(number_of_frames,w,x))/2;
%     end
% end
Div=zeros(N,1);
for w=1:N
Div(w,1) = sum(Z1_final{w}(1:R));
end


count = 0;
Meter_Readings = abs(File_data3(:)*1000);                                       % Retrieves maximum reading to find richness.

for i=1:size(Meter_Readings)
    if Meter_Readings(i) > 3000
        count = count + 1;
        Meter_Readings(i) = 3000;
    end
    if Meter_Readings(i) == 0
        Meter_Readings(i) = 1;
    end
end
Maximum_Reading = max(Meter_Readings);                                      % Retrieves maximum reading to find richness.
Minimum_Reading = min(Meter_Readings);                                      % Retrieves minimum reading to find richness.
R = ceil((Maximum_Reading-Minimum_Reading)/SW);                             % Calculates the value R (Number of unique species)


number_of_frames = ceil(datenum(File_data1(end)-File_data1(1))/frame_size);


% % % P_matrix_a = zeros(1, N, R);                                                         % Declaring matrix to store probability of R species of N meters
% % % 
% % % for z=1:size(File_data3)
% % %     meter_number = find(I==File_data2(z));
% % %             frame_number = ceil((datenum(File_data1(z) - File_data1(1)) + 1) / frame_size);
% % %             species_number = ceil(abs(Meter_Readings(z))/SW);
% % %         if frame_number==2
% % %             P_matrix_a(frame_number-1,meter_number,species_number)= P_matrix_a(frame_number-1,meter_number,species_number) + 1;
% % %         end
% % % end
% % % 
% % % P_a = zeros(1, N, R);                                                         % Declaring matrix to store probability of R species of N meters
% % % for x=1:1
% % %     for y=1:N
% % %         for z=1:R
% % %             P_a(x,y,z) = (P_matrix_a(x,y,z)+1)/(sum(P_matrix_a(x,y,:))+R);
% % %         end
% % %     end
% % % end
% % % 
% % % Z1_final_a{N} = [];
% % % for w=1:N
% % %     for x=1:R
% % %         for y=1:R
% % %         Z1_final_a{w}(x,y) =  abs(P_a(1,w,x)-P(1,w,y));
% % %         end
% % %     end
% % % end
% % % 
% % % ordinariness1{N}=[];
% % % for w=1:N
% % % for i=1:R
% % %     ordinariness1{w}(i)=0;
% % %     for j=1:R
% % %         if i==j
% % %         ordinariness1{w}(i) = ordinariness1{w}(i) + (Z1_final_a{w}(i,j)*P(1,w,i))^(0.5);
% % %         end
% % %     end
% % % end
% % % end
% % % 
% % % Div_a=zeros(N,1);
% % % 
% % % for w=1:N
% % %         Div_a(w,1) = Div_a(w,1) + sum(ordinariness1{w});
% % % end
% % % disp(Div_a);



    
% for z=1:size(File_data3)
%         meter_number = find(I_test==File_data2(z));
%         if ismember(I_test(meter_number),I)
%             meter_number1 = find(I==I_test(meter_number));
%             frame_number = ceil((datenum(File_data1(z) - File_data1(1)) + 1) / frame_size);
%             species_number = ceil(abs(Meter_Readings(z))/SW);
%             P_matrix_a(frame_number,meter_number1,species_number)= P_matrix_a(frame_number,meter_number1,species_number) + 1;
%         end
% end
% 
% P_a = zeros(number_of_frames, N, R);                                                         % Declaring matrix to store probability of R species of N meters
% for x=1:number_of_frames
%     for y=1:N
%         for z=1:R
%             P_a(x,y,z) = (P_matrix_a(x,y,z)+1)/(sum(P_matrix_a(x,y,:))+R);
%         end
%     end
% end
% 
% % Z_a = zeros(N, number_of_frames, R, R);
% % for v = 1:N
% %     for w = 1:number_of_frames
% %         for x = 1:R
% %             for y=1:R
% %                 Z_a(v,w,y,x) = P_a(w,v,y);
% %             end
% %         end
% %     end
% % end
% 
% Z1_final_a{N} = [];
% for w=1:N
% for x=1:R
%         Z1_final_a{w}(x) =  0;
%         for y=1:number_of_frames
%             Z1_final_a{w}(x) = Z1_final_a{w}(x) + abs(P_a(y,w,x)-P(y,w,x));
%         end
%         Z1_final_a{w}(x) = Z1_final_a{w}(x)/number_of_frames;
% end
% end
% 
% % P_final_a{N}=[];
% % for w = 1:N    
% %     for x=1:R
% %         P_final_a{w}(x,1) = (P_a(1,w,x)+P_a(number_of_frames,w,x))/2;
% %     end
% % end
% 
% Div_a=zeros(N,1);
% 
% for w=1:N
%     Div_a(w,1) = sum(Z1_final_a{w}(1:R));
% end
% disp(Div_a);

affected = randperm(M,M);
affected=affected;%Select M random numbers for selecting columns of attacked meters
affected_meters_array=sort(affected);                                                                %Sorts the vector in previous step and loads to new vector affected_meters
affected_meters=zeros(length(affected_meters_array),1);
for x=1:length(affected_meters_array)
    affected_meters(x) = I(affected_meters_array(x));
end

dmin=80;                                                                                       %Loading Delta min
dmax=120;                                                                                      %Loading Delta max

House_partition{N}=[];
House_partition1{N}=[];
for x=1:size(File_data1)                                                                        % value j Loops the first column of File_data for the house in the previous loop
        meter_number = find(I==File_data2(x));
            House_partition{meter_number} = [House_partition{meter_number} abs(File_data3(x)*1000)];
end

Meterwise_unaffected{N}=[];
Meterwise_affected{N}=[];

for x=1:size(I)
    Meterwise_unaffected{x}=House_partition{x};
end

count=1;

for x=1:size(I)
        if x==find(I==affected_meters(count))
            d=dmin+(dmax-dmin).*rand(length(House_partition{x}),1);       
            for y=1:length(House_partition{x})
                House_partition{x}(y)=House_partition{x}(y)+d(y,1);
            end
            count=count+1;
            if count>M
                break;
            end
        end
end

for x=1:size(I)
    Meterwise_affected{x}=House_partition{x};
end

File_data_affected1 = File_data1;
File_data_affected2 = File_data2;
File_data_affected3 = File_data3;
File_data_affected4 = File_data4;
File_data_affected5 = File_data5;

for x=1:size(File_data1)
    meter_number = find(I==File_data2(x));
            File_data_affected3(x) = House_partition{meter_number}(1);
            House_partition{meter_number}=House_partition{meter_number}(2:end);
end

count = 0;
Meter_Readings_affected = File_data_affected3(:);                                       % Retrieves maximum reading to find richness.
for i=1:size(Meter_Readings_affected)
    if Meter_Readings_affected(i) > 3000
        count = count + 1;
        Meter_Readings_affected(i) = 3000;
    end
    if Meter_Readings_affected(i) == 0
        Meter_Readings_affected(i) = 1;
    end
end
Maximum_Reading = max(Meter_Readings_affected);                                      % Retrieves maximum reading to find richness.
Minimum_Reading = min(Meter_Readings_affected);                                      % Retrieves minimum reading to find richness.
R = ceil((Maximum_Reading-Minimum_Reading)/SW);                             % Calculates the value R (Number of unique species)


P_matrix_affected = zeros(24, N, R);                                                         % Declaring matrix to store probability of R species of N meters

for z=1:size(File_data3)/4
        frame_number = ceil((datenum(File_data1(z) - File_data1(1)) + 1) / frame_size);
            meter_number = find(I==File_data2(z));
            species_number = ceil(abs(Meter_Readings_affected(z))/SW);
            P_matrix_affected(frame_number,meter_number,species_number)= P_matrix_affected(frame_number,meter_number,species_number) + 1;
end

P_affected = zeros(24, N, R);                                                         % Declaring matrix to store probability of R species of N meters
for x=1:24
    for y=1:N
        for z=1:R
            P_affected(x,y,z) = (P_matrix_affected(x,y,z)+1)/(sum(P_matrix_affected(x,y,:))+R);
        end
    end
end
disp(3);

% Z_affected = zeros(N, number_of_frames, R, R);
% for v = 1:N
%     for w = 1:number_of_frames
%         for x = 1:R
%             for y=1:R
%                 Z_affected(v,w,y,x) = P_affected(w,v,y);
%             end
%         end
%     end
% end

for w=1:N
for x=1:R
for y=1:1
wei(y,w,x)=((P(y,w,x)-min(P(y,w,:)))/(max(P(y,w,:))-min(P(y,w,:))));
end
end
end


Z1_final_affected{N} = [];
for w=1:N
    for x=1:R
        for y=1:R
        Z1_final_affected{w}(x,y) =  abs(P_affected(2,w,x)-P(1,w,y));
        end
    end
end


% Z1_final_affected{N} = [];
% for w=1:N
%     for x=1:R
%         Z1_final_affected{w}(x) =  0;
%         for y=1:1
% %             [sortedP, sortedInds] = sort(P(y,w,:),'descend');
% %             place = find(sortedP==P(y,w,x));
% %             weight = 2/(1.2^(place(1)-1));
%             Z1_final_affected{w}(x) = Z1_final_affected{w}(x) + abs(P_affected(y,w,x)-P(y,w,x));
%         end
%         %Z1_final_affected{w}(x) = Z1_final_affected{w}(x);
%     end
% end

% P1_affected_final{N}=[];
% P1_affected_first{N}=[];
% for w = 1:N    
%     for x=1:R
%         P1_affected_final{w}(x,1) = (P_affected(number_of_frames,w,x)+P_affected(1,w,x))/2;
%     end
% end

ordinariness{N}=[];
for w=1:N
for i=1:R
    ordinariness{w}(i)=0;
    for j=1:R
        if i==j
        ordinariness{w}(i) = ordinariness{w}(i) + (Z1_final_affected{w}(i,j)*P(1,w,i))^0.5;
        end
    end
end
end

Div_affected=zeros(N,1);
for w=1:N
for i=1:R
Div_affected(w,1) = Div_affected(w,1) + (ordinariness{w}(i)*(1-P(1,w,i)));
end
end

% disp(Div_affected);

% % % for w=1:N
% % %     Div_affectedM(w,1) = Div_affected(w,1)-Div_a(w,1);
% % % end


% idx=kmeans(Div_affected,2);
% hold on;
% for i=1:size(idx)
% if i>60
%     plot(i,Div_affected(i),'bo','MarkerSize',10,'LineWidth',2);
% else
%     plot(i,Div_affected(i),'rx','MarkerSize',10,'LineWidth',2);
% end
% end



hold on;
for i=1:N
if i>M
LH(1)=plot(i,Div_affected(i),'bo','MarkerSize',10,'LineWidth',2);
else
LH(2)=plot(i,Div_affected(i),'rx','MarkerSize',10,'LineWidth',2);
end
end

count1=0;
count2=0;
for i=1:M
if Div_affected(i)<0.55
count1 = count1+1;
end
end
for i=M+1:N
if Div_affected(i)>0.55
count2 = count2+1;
end
end
disp(count1/M);
disp(count2/(N-M));





% dmin=120;                                                                                       %Loading Delta min
% dmax=140;                                                                                      %Loading Delta max
% 
% House_partition{N}=[];
% House_partition1{N}=[];
% for x=1:size(File_data1)                                                                        % value j Loops the first column of File_data for the house in the previous loop
%         meter_number = find(I==File_data2(x));
%             House_partition{meter_number} = [House_partition{meter_number} abs(File_data3(x)*1000)];
% end
% 
% Meterwise_unaffected{N}=[];
% Meterwise_affected{N}=[];
% 
% for x=1:size(I)
%     Meterwise_unaffected{x}=House_partition{x};
% end
% 
% count=1;
% 
% for x=1:size(I)
%         if x==find(I==affected_meters(count))
%             d=dmin+(dmax-dmin).*rand(length(House_partition{x}),1);       
%             for y=1:length(House_partition{x})
%                 House_partition{x}(y)=House_partition{x}(y)+d(y,1);
%             end
%             count=count+1;
%             if count>M
%                 break;
%             end
%         end
% end
% 
% for x=1:size(I)
%     Meterwise_affected{x}=House_partition{x};
% end
% 
% File_data_affected1 = File_data1;
% File_data_affected2 = File_data2;
% File_data_affected3 = File_data3;
% File_data_affected4 = File_data4;
% File_data_affected5 = File_data5;
% 
% for x=1:size(File_data1)
%     meter_number = find(I==File_data2(x));
%             File_data_affected3(x) = House_partition{meter_number}(1);
%             House_partition{meter_number}=House_partition{meter_number}(2:end);
% end
% 
% count = 0;
% Meter_Readings_affected = File_data_affected3(:);                                       % Retrieves maximum reading to find richness.
% for i=1:size(Meter_Readings_affected)
%     if Meter_Readings_affected(i) > 3000
%         count = count + 1;
%         Meter_Readings_affected(i) = 3000;
%     end
%     if Meter_Readings_affected(i) == 0
%         Meter_Readings_affected(i) = 1;
%     end
% end
% Maximum_Reading = max(Meter_Readings_affected);                                      % Retrieves maximum reading to find richness.
% Minimum_Reading = min(Meter_Readings_affected);                                      % Retrieves minimum reading to find richness.
% R = ceil((Maximum_Reading-Minimum_Reading)/SW);                             % Calculates the value R (Number of unique species)
% 
% 
% P_matrix_affected = zeros(24, N, R);                                                         % Declaring matrix to store probability of R species of N meters
% 
% for z=1:size(File_data3)/4
%         frame_number = ceil((datenum(File_data1(z) - File_data1(1)) + 1) / frame_size);
%             meter_number = find(I==File_data2(z));
%             species_number = ceil(abs(Meter_Readings_affected(z))/SW);
%             P_matrix_affected(frame_number,meter_number,species_number)= P_matrix_affected(frame_number,meter_number,species_number) + 1;
% end
% 
% P_affected = zeros(24, N, R);                                                         % Declaring matrix to store probability of R species of N meters
% for x=1:24
%     for y=1:N
%         for z=1:R
%             P_affected(x,y,z) = (P_matrix_affected(x,y,z)+1)/(sum(P_matrix_affected(x,y,:))+R);
%         end
%     end
% end
% ab=0.29;
% bb=0.08;
% nu=0.05;
% Z1_final_affected{N} = [];
% temp1{N}=[];
% temp2{N}=[];
% for w=1:N
% for x=1:R
% for y=1:R
% temp1{w}(x,y) =  abs(P_affected(2,w,x)-P(1,w,y));
% temp2{w}(x,y) = 1/((1+ab*exp(-1*bb*100*temp1{w}(x,y)))^(1/nu));
% Z1_final_affected{w}(x,y) = temp2{w}(x,y);
% end
% end
% end
% ordinariness{N}=[];
% for w=1:N
% for i=1:R
% ordinariness{w}(i)=0;
% for j=1:R
% if i==j
% ordinariness{w}(i) = ordinariness{w}(i) + (Z1_final_affected{w}(i,j)*P(1,w,i))^0.5;
% end
% end
% end
% end
% Div_affected=zeros(N,1);
% for w=1:N
% for i=1:R
% Div_affected(w,1) = Div_affected(w,1) + (ordinariness{w}(i)*(1-P(1,w,i)));
% end
% end
% hold on;
% for i=1:N
% if i>M
% LH(1)=plot(i,Div_affected(i),'bo','MarkerSize',10,'LineWidth',2);
% else
% LH(2)=plot(i,Div_affected(i),'rx','MarkerSize',10,'LineWidth',2);
% end
% end
% count1=0;
% count2=0;
% for i=1:M
% if Div_affected(i)<0.55
% count1 = count1+1;
% end
% end
% for i=M+1:N
% if Div_affected(i)>0.55
% count2 = count2+1;
% end
% end
% disp(count1/M);
% disp(count2/(N-M));






% ab=0.5;
% bb=0.07;
% nu=0.05;
% hold on;
% for i=1:100
% res(i) = 1/((1+ab*exp(-1*bb*i))^(1/nu));
% plot(i,res(i),'ro');
% end


% 
% for w=1:N
% for x=1:R
% for y=1:1
% wei(y,w,x)=((P(y,w,x)-min(P(y,w,:)))/(max(P(y,w,:))-min(P(y,w,:))));
% end
% end
% end
% 
% Z1_final_affected{N} = [];
% for w=1:N
% for x=1:R
% Z1_final_affected{w}(x) =  0;
% for y=1:1
% Z1_final_affected{w}(x) = Z1_final_affected{w}(x) + (abs(P_affected(y,w,x)-P(y,w,x)))*(wei(y,w,x));
% end
% end
% end
% 
% Div_affected=zeros(N,1);
% for w=1:N
% Div_affected(w,1) = sum(Z1_final_affected{w}(1:R));
% end
% idx=kmeans(Div_affected,2);
% hold on
% for i=1:size(idx)
% if idx(i)==1
% plot(i,Div_affected(i),'bo');
% else
% plot(i,Div_affected(i),'ro');
% end
% end







% 
% final=[];
% count=1;   
% for x=1:N
%     if x==find(I==affected_meters(count))
%        final = [final Div(x)];
%        count=count+1;
%        if count>M
%            break;
%        end
%     end
% end
% 
% count=1;   
% for x=1:N
%     if x==find(I==affected_meters(count))
%        final = [final Div_affected(x)];
%        count=count+1;
%        if count>M
%            break;
%        end
%     end
% end
% 
% count=1;   
% for x=1:N
%     if x==find(I==affected_meters(count))
%        final = [final Div_a(x)];
%        count=count+1;
%        if count>M
%            break;
%        end
%     end
% end
% 
% final_ratio=[];
% count=1;   
% for x=1:N
%     if x==find(I==affected_meters(count))
%        final_ratio = [final_ratio Div_affected(x)/Div(x)];
%        count=count+1;
%        if count>M
%            break;
%        end
%     end
% end
% 
% count=1;   
% for x=1:N
%     if x==find(I==affected_meters(count))
%        final_ratio = [final_ratio Div_a(x)/Div(x)];
%        count=count+1;
%        if count>M
%            break;
%        end
%     end
% end
% 
% count=1;   
% for x=1:N
%     if x==find(I==affected_meters(count))
%        final_ratio = [final_ratio Div_affected(x)/Div_a(x)];
%        count=count+1;
%        if count>M
%            break;
%        end
%     end
% end

% 
% unaffected_meters_array=[];
% count=1;   
% for x=1:N
%     if x~=find(I==affected_meters_array(count))
%        unaffected_meters_array = [unaffected_meters_array x];
%        count=count+1;
%     end
% end
% 
% final_ratio_unaffected=[];
% count=1;   
% for x=1:N
%     if x~=find(I==affected_meters(count))
%        final_ratio_unaffected = [final_ratio_unaffected Div_a(x)/Div(x)];
%        count=count+1;
%     end
% end
        
% for w=1:N
%     res(w)=sum(Z1_final_affected{w}(1:R));
% end
%                
% idx = kmeans(res',2);
% 
% for i=1:N
%     if idx(i)==1
%         plot(i,res(i),'b.');
%         hold on;
%     end
%     if idx(i)==2
%         plot(i,res(i),'r.');
%         hold on;
%     end
% end
% 
% count1=0;
% count2=0;
% for i=1:300
%     if idx(i)==1
%         count1=count1+1;
%     end
% end
% for i=301:1275
%     if idx(i)==2
%         count2=count2+1;
%     end
% end
% disp(count1);
% disp(count2);
% count1_array(abc)=count1;
% count2_array(abc)=count2;
% disp(count1_array);
end


% for w=1:N
% for x=1:R
% for y=1:R
% Z1_final_affected_test{w}(x,y) = (abs(Z_affected(w,2,x,y)-Z(w,2,x,y))+abs(Z_affected(w,1,x,y)-Z(w,1,x,y)))/2;
% end
% end
% end
% idx=[];
% for w=1:N
% res_test(w)=sum(Z1_final_affected_test{w}(1:R));
% end
% idx = kmeans(res_test',2);
% count1=0;
% count2=0;
% for i=1:300
% if idx(i)==1
% count1=count1+1;
% end
% end
% for i=301:1275
% if idx(i)==2
% count2=count2+1;
% end
% end
% disp(count1);
% disp(count2);


% idx = kmeans(Div_affected,2);
% 
% for i=1:N
%     if idx(i)==1
%         plot(i,Div_affected(i),'b.');
%         hold on;
%     end
%     if idx(i)==2
%         plot(i,Div_affected(i),'r.');
%         hold on;
%     end
% end
%  
% count1=0;
% count2=0;
% for i=1:300
%     if idx(i)==1
%         count1=count1+1;
%     end
% end
% for i=301:1275
%     if idx(i)==2
%         count2=count2+1;
%     end
% end
% disp(count1);
% disp(count2);
% count1_array(abc)=count1;
% count2_array(abc)=count2;