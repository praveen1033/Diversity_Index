Contains MATLAB code for detecting data falsification in AMI / Smartgrid using the model from following paper.

PAPER NAME: A Diversity Index based Scoring Framework for Identifying Smart Meters Launching Stealthy Data Falsification Attacks.

Abstract: A challenging problem in Advanced Metering Infrastructure (AMI) of smart grids is the identification of smart meters under the control of a stealthy adversary, that inject very low margins of stealthy data falsification. The problem is challenging due to wide legitimate variation in both individual and aggregate trends in real world power consumption data, making such stealthy attacks unrecognizable by existing approaches. In this paper, via proposed modified diversity index scoring metric, we propose a novel information-theory inspired data driven device anomaly classification framework to identify compromised meters launching low margins of stealthy data falsification attacks. Specifically, we draw a parallelism between the effects of data falsification attacks and ecological balance disruptions and identify required mathematical modifications in existing Renyi Entropy and Hill’s Diversity Entropy measures. These modifications such as expected self-similarity with weighted abundance shifts across various temporal scales, and diversity order are appropriately embedded in our resulting framework. The resulting diversity index score is used to classify smart meters launching additive, deductive, and alternating switching attack types with high sensitivity (as low as 100W) compared to the existing works that perform poorly at margins of false data below 400W. Our proposed theory is validated with two different real smart meter datasets from USA and Ireland. Experimental results demonstrate successful detection sensitivity from very low to high margins of false data, thus reducing undetectable strategy space of attacks in AMI for adversary having complete knowledge of our method

DATASET SOURCE: 

DATA DESCRIPTION: Pecan Street Texas Dataset
  The smart grid data is available in time intervals of 1 hour. 
  The power consumption in the data is in kilowatts (kW) and 1 kilowatt = 1000 watts.
  A sample of the smart meter data set is shown below.
  
| time                | id   | usage  | generation | grid    |
|---------------------|------|--------|------------|---------|
| 01/01/2016 00:00:00 | 26   | 0.5271 | -0.00126   | 0.52586 |
| 01/01/2016 00:00:00 | 9052 | 0.3711 | -0.00723   | 0.53711 |
| 01/01/2016 01:00:00 | 26   | 0.4631 | -0.00006   | 0.46231 |

We need the first 3 columns for this paper. 
The first column shows the time. The sample shows we have the data readings per each hour.
The second column shows the smart meter id number. This could be used for meter by meter seperation of data.
The third column is the usage recorded in kilowatts. We will convert it to watts for our project.

We will start with how to load the datafile. First, we will see how to load the file using MATLAB.

To load the csv file to matrix (we named File_data), use the following command
File_data = csvread('File6new.csv',1,1);                                      % Loads file to matrix called File_data

To extract the vector of list of smart meters, which is the list of smart meter ids
I = unique(File_data(:,2));                                                   % I includes all smart meter ids.

To get the number of smart meters in the dataset, we have to find the length of the vector for unique smart meter ids
N = length(I);                                                                % N is the number of smart meters.

To convert the data into watts, we have to multiply each reading by 1000
Meter_Readings = File_data(:,3)*1000;                                         % Converts the usage in kilowatts to watts.

Now, we will see how to clean the data.

WHY FILTERING NECESSARY ?


FILTERING LOGIC? 



