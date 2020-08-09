
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

number_of_files_T4 = 18
number_of_files_ute = 10

length_T4 = number_of_files_T4                          #Total number of files to be processed
data_point= [5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22]
range_len_T4 = np.zeros(length_T4)                          #Pre-allocating array for intesity - Length, Sum and Average
sum_range_T4 = np.zeros(length_T4)
ave_range_T4 = np.zeros(length_T4)
Distance_T4 = np.zeros(length_T4)

length_ute = number_of_files_ute                            #Total number of files to be processed
range_len_ute = np.zeros(length_ute)                          #Pre-allocating array for intesity - Length, Sum and Average
sum_range_ute = np.zeros(length_ute)
ave_range_ute = np.zeros(length_ute)
Distance_ute = np.zeros(length_ute)

for i in range(5, number_of_files_T4+5):                 #NOTE CHANGE: for i in range(5, number_of_files+6):
    df = pd.read_csv("T4_{}.csv".format(i))              #Read in the first file and subsequent files after looping
    range_T4 = df.range                                #Read the Intensity column of the CSV file
    range_len_T4[i-5] = len(range_T4)                        #Find the length of the column
    sum_range_T4[i-5] = sum(range_T4)                        #Find the sum of the column
    ave_range_T4[i-5] = sum_range_T4[i-5] / range_len_T4[i-5]  #Find the average of the column

    # ave_range_T4[i-5] = ave_range_T4[i-5] -30               #convert from dbm to db
    if (i == 5):
        Distance_T4[0]=5
    else:
        Distance_T4[i-5] = Distance_T4[i-6] + 1.0            #Set the distance to start from 0 and increment up
#### Please note: The name of the files increment by 1, but the distace taken of each file increments by 2.5.
# The name of the file was changed to increment by 1 instead of 2.5 to make the code easier.
for i in range(5, number_of_files_ute+5):                 #NOTE CHANGE: for i in range(1, number_of_files+1):
    df = pd.read_csv("Ute_{}.csv".format(i))              #Read in the first file and subsequent files after looping
    range_ute = df.range                                #Read the Intensity column of the CSV file
    range_len_ute[i-5] = len(range_ute)                        #Find the length of the column
    sum_range_ute[i-5] = sum(range_ute)                        #Find the sum of the column
    ave_range_ute[i-5] = sum_range_ute[i-5] / range_len_ute[i-5]  #Find the average of the column

    if (i == 5):
        Distance_ute[0] = 5
    else:
        Distance_ute[i-5] = Distance_ute[i-6] + 2.50            #Set the distance to start from 0 and increment up


fig = plt.figure()
plt.plot(Distance_T4, data_point,'k--', label ="Calibration Line - Manually Measured Range")
# plt.plot(Distance_T4, data_point,'k*')
plt.plot(Distance_T4, ave_range_T4,'r', linewidth=2, label ="Radar Range Measurement - Triangular Target T4")
# plt.plot(Distance_T4, ave_range_T4,'ro')
plt.plot(Distance_ute, ave_range_ute, 'b-', linewidth=2, label='Radar Range Measurement - Utility Vehicle')
plt.legend()
plt.grid()
plt.xlabel('Measured Distance (m)')
plt.ylabel('Radar Distance (m)')
plt.title('Evaluation of Range Accuracy')
plt.savefig('RangeComparisson')
plt.show()
