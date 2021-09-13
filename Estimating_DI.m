
count1_array=[];
count2_array=[];
Div_affected=zeros(192,2);
Div_affected_est=zeros(192,2);
Div=zeros(192,2);
for abc=1:2
filename = 'C:\Users\vmcx3\Desktop\Data sets\Electricity data\TOSG_cont\2016.csv';
delimiter = ',';
startRow = 2;

formatSpec = '%C%f%f%f%f%[^\n\r]';

fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string', 'EmptyValue', NaN, 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');

fclose(fileID);

% Post processing for unimportable data.
% No unimportable data rules were applied during the import, so no post
% processing code is included. To generate code which works for
% unimportable data, select unimportable cells in a file and regenerate the
% script.

% Create output variable
File_data = table(dataArray{1:end-1}, 'VariableNames', {'localminute','dataid','use','gen','grid'});

% Clear temporary variables
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
M=60;
SW = 100;                                                                   % SW is the species width in watts.
R = 0;                                                                      % Initialising richness (Number of unique species) to 0.
q = 0.5;
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
frame_size = 20;
% max_day = max_day - rem((max_day - min_day),frame_size) - 1;
number_of_frames = ceil(datenum(File_data1(end)-File_data1(1))/frame_size);

P_matrix = zeros(1, N, R);                                                         % Declaring matrix to store probability of R species of N meters
for z=1:size(File_data3)*abc/12
    frame_number = ceil((datenum(File_data1(z) - File_data1(1)) + 1) / frame_size);
     if frame_number==abc
            meter_number = find(I==File_data2(z));
            species_number = ceil(abs(Meter_Readings(z))/SW);
            P_matrix(1,meter_number,species_number)= P_matrix(1,meter_number,species_number) + 1;
    end
end

P1 = zeros(1, N, R); 

P = zeros(1, N, R);                                                         % Declaring matrix to store probability of R species of N meters
for x=1:1
for y=1:N
for z=1:R
P(x,y,z) = (P_matrix(x,y,z)+1)/(sum(P_matrix(x,y,:))+R);
end
end
end




affected = randperm(M,M);
affected=affected;                                                                                             %Select M random numbers for selecting columns of attacked meters
affected_meters_array=sort(affected);                                                                %Sorts the vector in previous step and loads to new vector affected_meters
affected_meters=zeros(length(affected_meters_array),1);
for x=1:length(affected_meters_array)
    affected_meters(x) = I(affected_meters_array(x));
end

dmin=350;                                                                                       %Loading Delta min
dmax=450;                                                                                      %Loading Delta max

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
            [House_partition{x},o]=sort(House_partition{x});
            d=sort(d,'descend');            
            for y=1:length(House_partition{x})
                House_partition{x}(y)=House_partition{x}(y)+d(y,1);
            end
            for i=1:length(House_partition{x})
                House_partition1{x}(o(i))=House_partition{x}(i);
            end
            House_partition{x}=House_partition1{x};
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


P_matrix_affected = zeros(1, N, R);                                                         % Declaring matrix to store probability of R species of N meters

for z=1:size(File_data3)*abc/12
        frame_number = ceil((datenum(File_data1(z) - File_data1(1)) + 1) / frame_size);
        if frame_number==abc+1
            meter_number = find(I==File_data2(z));
            species_number = ceil(abs(Meter_Readings_affected(z))/SW);
            P_matrix_affected(frame_number-abc,meter_number,species_number)= P_matrix_affected(frame_number-abc,meter_number,species_number) + 1;
        end
end

P_affected = zeros(1, N, R);                                                         % Declaring matrix to store probability of R species of N meters
for x=1:1
    for y=1:N
        for z=1:R
            P_affected(x,y,z) = (P_matrix_affected(x,y,z)+1)/(sum(P_matrix_affected(x,y,:))+R);
        end
    end
end



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
        Z1_final_affected{w}(x,y) =  abs(P_affected(1,w,x)-P(1,w,y));
        end
    end
end


ordinariness{N}=[];
for w=1:N
for i=1:R
    ordinariness{w}(i)=0;
    for j=1:R
        if i==j
        ordinariness{w}(i) = ordinariness{w}(i) + (Z1_final_affected{w}(i,j)*P(1,w,i))^q;
        end
    end
end
end


for w=1:N
for i=1:R
Div_affected(w,abc) = Div_affected(w,abc) + (ordinariness{w}(i)*(1-P(1,w,i)));
end
end




%..................................................................Estimation......................................................

shift = ((dmax+dmin)/2)/SW;

shift_flr = floor(shift);
shidt_rem = shift - floor(shift);
 

P_matrix_affected_est = zeros(1, N, R);                                                         % Declaring matrix to store estimated probability of R species of N meters



for z = 1:N
        for y = 1:R
            if y <= shift_flr
                P_matrix_affected_est(1,z,y) = 0;
            else
                P_matrix_affected_est(1,z,y)= (1-shidt_rem)*P_matrix(1,z,y-shift_flr) + shidt_rem*P_matrix(1,z,y-shift_flr+1);
            end
            if y == R
                P_matrix_affected_est(1,z,y) = P_matrix_affected_est(1,z,y) + sum(P_matrix(1,z,y-shift_flr+1:y));
            end
        end
end

P_affected_est = zeros(1, N, R);                                                         % Declaring matrix to store probability of R species of N meters
for x=1:1
    for y=1:N
        for z=1:R
            P_affected_est(x,y,z) = (P_matrix_affected_est(x,y,z)+1)/(sum(P_matrix_affected_est(x,y,:))+R);
        end
    end
end



for w=1:N
for x=1:R
for y=1:1
wei(y,w,x)=((P(y,w,x)-min(P(y,w,:)))/(max(P(y,w,:))-min(P(y,w,:))));
end
end
end


Z1_final_affected_est{N} = [];
for w=1:N
    for x=1:R
        for y=1:R
        Z1_final_affected_est{w}(x,y) =  abs(P_affected_est(1,w,x)-P(1,w,y));
        end
    end
end


ordinariness_est{N}=[];
for w=1:N
for i=1:R
    ordinariness_est{w}(i)=0;
    for j=1:R
        if i==j
        ordinariness_est{w}(i) = ordinariness_est{w}(i) + (Z1_final_affected_est{w}(i,j)*P(1,w,i))^q;
        end
    end
end
end


for w=1:N
for i=1:R
Div_affected_est(w,abc) = Div_affected_est(w,abc) + (ordinariness_est{w}(i)*(1-P(1,w,i)));
end
end


hold on;
for i=1:N
if Div_affected(i)<0.55
    plot(i,Div_affected_est(i),'bo','MarkerSize',10,'LineWidth',2);
else
    plot(i,Div_affected_est(i),'rx','MarkerSize',10,'LineWidth',2);
end
end



L(1) = "Real DI Score";
L(2) = "Estimated DI Score";

hold on;

for i=1:M
    if i<M
        plot(i,Div_affected_est(i),'rx','MarkerSize',10,'LineWidth',2);
    else
        LH(1) = plot(i,Div_affected_est(i),'rx','MarkerSize',10,'LineWidth',2);
    end
end

for i=1:M
    if i<M
        plot(i,Div_affected(i),'bo','MarkerSize',10,'LineWidth',2);
    else
        LH(2) = plot(i,Div_affected(i),'bo','MarkerSize',10,'LineWidth',2);
    end
end

legend(LH,L);
ylabel('Diversity Index');
xlabel('Smart Meter ID');
title('Estimating Diversity Index')


end

