# Matthew Barrett 20 July 2020
# The following code does a number of things
# 1. It calculates the length A and A*sqrt(2) to manufacture RCS target of a known size
# 2. It checks the size A and A*sqrt(2) is correctly calculated
# 3. It calculates the actual RCS of the targets already made that were incorrectly calculated to begin with
# 4. it process the data obtained by importing each CSV file, taking the average of the intensity and saving this
#     intensity to array
# 5. It plots the intensity data saved against distance.
import interpolate as interpolate
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
#######################################################################################################################
###### Calculate RCS Target size
#The following code will calculate the distance A and A2 for the radar cross section triangular targets based on the
#RCS specified. The quantities A and A2 are in mm
# num_Targets = 7
# RCS = [400, 200, 120, 60, 30, 10, 5]
# RCS_lengthA_old = [78.3597e-3, 55.4074e-3, 42.9184e-3, 30.3479e-3, 0, 0, 0]
# RCSold = [78.3597e-3, 55.4074e-3, 42.9184e-3, 30.3479e-3]
# Triangular_Dis = np.zeros(num_Targets)
# Triangular_Dis_A = np.zeros(num_Targets)
# Triangular_Dis_A2 = np.zeros(num_Targets)
# RCS_check = np.zeros(num_Targets)
# RCS_check2 = np.zeros(num_Targets)
# alpha = 0.0039
# for iterate in range(0, num_Targets):
#     Triangular_Dis[iterate] = (RCS[iterate] * 3*pow(alpha, 2))/(4*np.pi)
#     Triangular_Dis_A[iterate] = Triangular_Dis[iterate]**(1./4.)
#     Triangular_Dis_A2[iterate] = Triangular_Dis_A[iterate] * (2**(1/2))
# print("RCS Triangular Target Length A")
# print(Triangular_Dis_A*1000)
# print("RCS Triangular Target Length 2*sqrt(A)")
# print(Triangular_Dis_A2*1000)
# ##### Back calculate to check calculation
# for it in range(0, num_Targets):
#     RCS_check[it] = (4*np.pi*pow(Triangular_Dis_A[it], 4)) / (3*pow(alpha, 2))
#     RCS_check2[it] = (4 * np.pi * pow(RCS_lengthA_old[it], 4)) / (3 * pow(alpha, 2))
# print("RCS Check of New Triangular calculations")
# print(RCS_check)
# print("Evaluate the old RCS measurements in m^2 for the targets originally incorrectly calculated")
# print(RCS_check2)
# #######################################################################################################################

#The below code processes and plots the first set of experiments
#The designation Values T1, T2, T3 and T4 represent the triangular Radar cross section tagets. T4 being the largest and
#T1 being the smallest
#set the number of files here
number_of_files_ute = 12
number_of_files_T6 = 6
number_of_files_T5 = 13
number_of_files_T4 = 18             #The number of files taken for the Target T4
number_of_files_T3 = 5              # The number of files taken for the Target T3
number_of_files_T2 = 10              # The number of files taken for the Target T2
number_of_files_T1 = 12             # The number of files taken for the Target T1

Pr_W = np.zeros(number_of_files_T1)
Pr_dB = np.zeros(number_of_files_T1)

Notarget = [-10, -10]                 #The sets the minimum intensity when no target is down range. ie sensor noise
Dis = [0, 32]                         #A 2 variable array to plot a dotted line for sensor noise


#######################################################################################################################
#######################################################################################################################
# Create empty arrays in order to find the average intensity for all files

length_ute = number_of_files_ute                            #Total number of files to be processed
int_len_ute = np.zeros(length_ute)                          #Pre-allocating array for intesity - Length, Sum and Average
sum_int_ute = np.zeros(length_ute)
ave_int_ute = np.zeros(length_ute)
Distance_ute = np.zeros(length_ute)

length_T6 = number_of_files_T6                            #Total number of files to be processed
int_len_T6 = np.zeros(length_T6)                          #Pre-allocating array for intesity - Length, Sum and Average
sum_int_T6 = np.zeros(length_T6)
ave_int_T6 = np.zeros(length_T6)
Distance_T6 = np.zeros(length_T6)

length_T5 = number_of_files_T5                            #Total number of files to be processed
int_len_T5 = np.zeros(length_T5)                          #Pre-allocating array for intesity - Length, Sum and Average
sum_int_T5 = np.zeros(length_T5)
ave_int_T5 = np.zeros(length_T5)
Distance_T5 = np.zeros(length_T5)

length_T4 = number_of_files_T4                            #Total number of files to be processed
int_len_T4 = np.zeros(length_T4)                          #Pre-allocating array for intesity - Length, Sum and Average
sum_int_T4 = np.zeros(length_T4)
ave_int_T4 = np.zeros(length_T4)
Distance_T4 = np.zeros(length_T4)

length_T3 = number_of_files_T3                            #Total number of files to be processed
int_len_T3 = np.zeros(length_T3)                          #Pre-allocating array for intesity - Length, Sum and Average
sum_int_T3 = np.zeros(length_T3)
ave_int_T3 = np.zeros(length_T3)
Distance_T3 = np.zeros(length_T3)

length_T2 = number_of_files_T2                            #Total number of files to be processed
int_len_T2 = np.zeros(length_T2)                          #Pre-allocating array for intesity - Length, Sum and Average
sum_int_T2 = np.zeros(length_T2)
ave_int_T2 = np.zeros(length_T2)
Distance_T2 = np.zeros(length_T2)

length_T1 = number_of_files_T1                            #Total number of files to be processed
int_len_T1 = np.zeros(length_T1)                          #Pre-allocating array for intesity - Length, Sum and Average
sum_int_T1 = np.zeros(length_T1)
ave_int_T1 = np.zeros(length_T1)
Distance_T1 = np.zeros(length_T1)

#Extraxct intensity from each range of CSV file for each target and then calculate average intensity
#### Please note: The name of the files increment by 1, but the distace taken of each file increments by 2.5.
# The name of the file was changed to increment by 1 instead of 2.5 to make the code easier.
for i in range(5, number_of_files_T1+5):                 #NOTE CHANGE: for i in range(1, number_of_files+1):
    df = pd.read_csv("T1_{}.csv".format(i))              #Read in the first file and subsequent files after looping
    int_T1 = df.intensity                                #Read the Intensity column of the CSV file
    int_len_T1[i-5] = len(int_T1)                        #Find the length of the column
    sum_int_T1[i-5] = sum(int_T1)                        #Find the sum of the column
    ave_int_T1[i-5] = sum_int_T1[i-5] / int_len_T1[i-5]  #Find the average of the column
    ave_int_T1[i-5] = ave_int_T1[i-5] -30               #convert from dbm to db
    if (i == 5):
        Distance_T1[0] = 5
    else:
        Distance_T1[i-5] = Distance_T1[i-6] + 2.50  # Set the distance to start from 0 and increment up

#### Please note: The name of the files increment by 1, but the distace taken of each file increments by 2.5.
# The name of the file was changed to increment by 1 instead of 2.5 to make the code easier.
for i in range(5, number_of_files_T2+5):                 #NOTE CHANGE: for i in range(1, number_of_files+1):
    df = pd.read_csv("T2_{}.csv".format(i))              #Read in the first file and subsequent files after looping
    int_T2 = df.intensity                                #Read the Intensity column of the CSV file
    int_len_T2[i-5] = len(int_T2)                        #Find the length of the column
    sum_int_T2[i-5] = sum(int_T2)                        #Find the sum of the column
    ave_int_T2[i-5] = sum_int_T2[i-5] / int_len_T2[i-5]  #Find the average of the column
    ave_int_T2[i-5] = ave_int_T2[i-5] -30               #convert from dbm to db
    if (i == 5):
        Distance_T2[0] = 5
    else:
        Distance_T2[i-5] = Distance_T2[i-6] + 2.5  # Set the distance to start from 0 and increment up

#### Please note: The name of the files increment by 1, but the distace taken of each file increments by 2.5.
# The name of the file was changed to increment by 1 instead of 2.5 to make the code easier.
for i in range(5, number_of_files_T3+5):                 #NOTE CHANGE: for i in range(1, number_of_files+1):
    df = pd.read_csv("T3_{}.csv".format(i))              #Read in the first file and subsequent files after looping
    int_T3 = df.intensity                                #Read the Intensity column of the CSV file
    int_len_T3[i-5] = len(int_T3)                        #Find the length of the column
    sum_int_T3[i-5] = sum(int_T3)                        #Find the sum of the column
    ave_int_T3[i-5] = sum_int_T3[i-5] / int_len_T3[i-5]  #Find the average of the column
    ave_int_T3[i-5] = ave_int_T3[i-5] -30               #convert from dbm to db
    if (i == 5):
        Distance_T3[0] = 5
    else:
        Distance_T3[i-5] = Distance_T3[i-6] + 2.50  # Set the distance to start from 0 and increment up

for i in range(5, number_of_files_T4+5):                 #NOTE CHANGE: for i in range(5, number_of_files+6):
    df = pd.read_csv("T4_{}.csv".format(i))              #Read in the first file and subsequent files after looping
    int_T4 = df.intensity                                #Read the Intensity column of the CSV file
    int_len_T4[i-5] = len(int_T4)                        #Find the length of the column
    sum_int_T4[i-5] = sum(int_T4)                        #Find the sum of the column
    ave_int_T4[i-5] = sum_int_T4[i-5] / int_len_T4[i-5]  #Find the average of the column
    ave_int_T4[i-5] = ave_int_T4[i-5] -30               #convert from dbm to db
    if (i == 5):
        Distance_T4[0]=5
    else:
        Distance_T4[i-5] = Distance_T4[i-6] + 1.0            #Set the distance to start from 0 and increment up

for i in range(5, number_of_files_T5+5):                 #NOTE CHANGE: for i in range(1, number_of_files+1):
    df = pd.read_csv("T5_{}.csv".format(i))              #Read in the first file and subsequent files after looping
    int_T5 = df.intensity                                #Read the Intensity column of the CSV file
    int_len_T5[i-5] = len(int_T5)                        #Find the length of the column
    sum_int_T5[i-5] = sum(int_T5)                        #Find the sum of the column
    ave_int_T5[i-5] = sum_int_T5[i-5] / int_len_T5[i-5]  #Find the average of the column
    ave_int_T5[i-5] = ave_int_T5[i-5] -30               #convert from dbm to db
    if (i == 5):
        Distance_T5[0] = 5
    else:
        Distance_T5[i-5] = Distance_T5[i-6] + 1.0            #Set the distance to start from 0 and increment up

for i in range(5, number_of_files_T6+5):                 #NOTE CHANGE: for i in range(1, number_of_files+1):
    df = pd.read_csv("T6_{}.csv".format(i))              #Read in the first file and subsequent files after looping
    int_T6 = df.intensity                                #Read the Intensity column of the CSV file
    int_len_T6[i-5] = len(int_T6)                        #Find the length of the column
    sum_int_T6[i-5] = sum(int_T6)                        #Find the sum of the column
    ave_int_T6[i-5] = sum_int_T6[i-5] / int_len_T6[i-5]  #Find the average of the column
    ave_int_T6[i-5] = ave_int_T6[i-5] -30               #convert from dbm to db
    if (i == 5):
        Distance_T6[0] = 5
    else:
        Distance_T6[i-5] = Distance_T6[i-6] + 1.0            #Set the distance to start from 0 and increment up

#### Please note: The name of the files increment by 1, but the distace taken of each file increments by 2.5.
# The name of the file was changed to increment by 1 instead of 2.5 to make the code easier.
for i in range(5, number_of_files_ute+5):                 #NOTE CHANGE: for i in range(1, number_of_files+1):
    df = pd.read_csv("Ute_{}.csv".format(i))              #Read in the first file and subsequent files after looping
    int_ute = df.intensity                                #Read the Intensity column of the CSV file
    int_len_ute[i-5] = len(int_ute)                        #Find the length of the column
    sum_int_ute[i-5] = sum(int_ute)                        #Find the sum of the column
    ave_int_ute[i-5] = sum_int_ute[i-5] / int_len_ute[i-5]  #Find the average of the column
    ave_int_ute[i-5] = ave_int_ute[i-5] -30               #convert from dbm to db
    if (i == 5):
        Distance_ute[0] = 5
    else:
        Distance_ute[i-5] = Distance_ute[i-6] + 2.50            #Set the distance to start from 0 and increment up


# #Interpolate data so that a straigt line fit can be made from each data set.
mT1, bT1 = np.polyfit(Distance_T1, ave_int_T1, 1)

mT2, bT2 = np.polyfit(Distance_T2, ave_int_T2, 1)

mT3, bT3 = np.polyfit(Distance_T3, ave_int_T3, 1)
mT4, bT4 = np.polyfit(Distance_T4, ave_int_T4, 1)
mT5, bT5 = np.polyfit(Distance_T5, ave_int_T5, 1)
mT6, bT6 = np.polyfit(Distance_T6, ave_int_T6, 1)
mTute, bTute = np.polyfit(Distance_ute, ave_int_ute, 1)

# #######################################################################################################################
# #The code below will plot the results

# #A plot of the theoretical return intensity vs range and the Actual intensity vs range
# #
# fig=plt.figure(7)
# plt.plot(Distance_T1, Pr_dB, 'r', label="Theoretical return intensity")
# plt.plot(Distance_T1, ave_int_T1, 'b', label="Experimental return intensity")
# plt.plot(Distance_T1, Pr_dB, 'ko')
# plt.plot(Distance_T1, ave_int_T1, 'k*')
#
# plt.legend()
# plt.xlabel('Distance (m)')
# plt.ylabel('Intensity (dB)')
# plt.title('Theoretical and Experimental Return Intenisity of Target T4', fontsize=14)
# plt.grid()
#

fig = plt.figure()
plt.plot(Distance_T5, ave_int_T5, 'ro', label="Best fit")
plt.plot(Distance_T5, mT5*Distance_T5+bT5, 'k', label="$RCS=30m^2$")
plt.legend()
plt.xlabel('Distance (m)')
plt.ylabel('Intensity (dB)')
plt.title('Radar Calibration Using RCS Triangular Targets ')
plt.savefig('RCS30m^2')

fig = plt.figure()
plt.plot(Distance_T6, ave_int_T6, 'ro', label="Best fit")
plt.plot(Distance_T6, mT6*Distance_T6+bT6, 'k', label="$RCS=10m^2$")
plt.legend()
plt.xlabel('Distance (m)')
plt.ylabel('Intensity (dB)')
plt.title('Radar Calibration Using RCS Triangular Targets ')
plt.savefig('RCS10m^2')

fig = plt.figure()
plt.plot(Distance_ute, ave_int_ute, 'ro', label ="Best fit")
plt.plot(Distance_ute, mTute*Distance_ute+bTute, 'k', label="Utility Vehicle")
plt.legend()
plt.xlabel('Distance (m)')
plt.ylabel('Intensity (dB)')
plt.title('Radar Calibration of Average Utility Vehicle ')
plt.savefig('Ute_Distance')

fig=plt.figure()
plt.plot(Distance_T1, mT1*Distance_T1+bT1, 'k', label="RCS=$400m^2$")
plt.plot(Distance_T2, mT2*Distance_T2+bT2, 'r', label="RCS=$200m^2$")
plt.plot(Distance_T3, mT3*Distance_T3+bT3, 'b', label="RCS=$120m^2$")
plt.plot(Distance_T4, mT4*Distance_T4+bT4, 'g', label="RCS=$60m^2$")
plt.plot(Dis,Notarget, 'm--', label="No Target")
plt.legend()
plt.xlabel('Distance (m)')
plt.ylabel('Intensity (dB)')
plt.title('Radar Calibration Using RCS Triangular Targets', fontsize=14)
plt.grid()
plt.savefig('BigTargets_SamePlot.png')
#
# A sub plot of all intensity vs range for each RCS target

fig = plt.figure(figsize=(2, 3))

sub1 = fig.add_subplot(221)
sub1.set_title('$400 m^2$ RCS Triangular Target')
sub1.plot(Distance_T1, ave_int_T1, 'ro', label="RCS=$400m^2$ raw data")
sub1.plot(Distance_T1, mT1*Distance_T1+bT1, 'k', label='RCS=$400m^2$ best fit')
plt.xlabel('Distance (m)')
plt.ylabel('Intensity (dB)')
plt.legend()
plt.grid()

sub2 = fig.add_subplot(222)
sub2.set_title('$200 m^2$ RCS Triangular Target')
sub2.plot(Distance_T2, ave_int_T2, 'ro', label="RCS=$200m^2$ raw Data")
sub2.plot(Distance_T2, mT2*Distance_T2+bT2, 'k', label='RCS=$200m^2$ best fit')
plt.xlabel('Distance (m)')
plt.ylabel('Intensity (dB)')
plt.legend()
plt.grid()

sub3 = fig.add_subplot(223)
sub3.set_title('$120 m^2$ RCS Triangular Target')
sub3.plot(Distance_T3, ave_int_T3, 'ro', label="RCS=$120m^2$ raw data")
sub3.plot(Distance_T3, mT3*Distance_T3+bT3, 'k', label="RCS=$120m^2$ best fit")
plt.xlabel('Distance (m)')
plt.ylabel('Intensity (dB)')
plt.legend()
plt.grid()

sub4 = fig.add_subplot(224)
sub4.set_title('$60 m^2$ RCS Triangular Target')
sub4.plot(Distance_T4, ave_int_T4, 'ro', label="RCS=$60m^2$ raw data")
sub4.plot(Distance_T4, mT4*Distance_T4+bT4, 'k', label='RCS=$60m^2$ best fit')
plt.xlabel('Distance (m)')
plt.ylabel('Intensity (dB)')
plt.legend()
plt.grid()
fig.suptitle('Maximum Detectable Range of Various Triangular RCS Targets')
# plt.tight_layout()
# plt.savefig('subPlot.png')

####################################################################################################################

#A plot of the theoretical return intensity vs range and the Actual intensity vs range
# Theoretical Radar range calculation - This code will calculate the recieved power a the same distances used in the
# experiments for trinagular target T4. The results are then plotted using code further down.
# The maximum returned distance from the
Pt = (pow(10, (26/10)))/1000 #Convert Output power from dBm to W
Gt = pow(10, (10.5/10)) # transmitter gain convert from dbi to a ratio
Gr = pow(10, (48/10)) #receiver gain agin convet to ratio
alpha = 0.0039
RCS= [400, 200, 120, 60, 30, 10]

for it in range(0, len(Distance_T1)):
    Pr_W[it] = (Pt * Gt * Gr * (pow(alpha, 2)) * RCS[0]) / ((pow((4*np.pi), 3)) * (pow(Distance_T1[it], 4)))
    # mW[it] = Pr_W[it]*1000
    Pr_dB[it] = (10 * np.log10(Pr_W[it])) #convert power in W to power in dbmW - this is for theoretical curve
    # ave_int_T1[it] = ave_int_T1[it]-30
    print(Distance_T1)
################################################
fig=plt.figure()
plt.plot(Distance_T1, Pr_dB, 'r', label="Theoretical return intensity")
plt.plot(Distance_T1, ave_int_T1, 'b', label="Experimental return intensity")
plt.plot(Distance_T1, Pr_dB, 'ko')
plt.plot(Distance_T1, ave_int_T1, 'k*')
plt.legend()
plt.xlabel('Distance (m)')
plt.ylabel('Intensity (dB)')
plt.title('Theoretical and Experimental Return Intenisity of Target $RCS = 400m^2$', fontsize=14)
plt.grid()
plt.savefig("Theoretical_and_Experimental.png")
plt.show()
