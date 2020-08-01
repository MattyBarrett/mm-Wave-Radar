import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

number_of_files_angle = 91

length_angle = number_of_files_angle                            #Total number of files to be processed
int_len_angle = np.zeros(length_angle)                          #Pre-allocating array for intesity - Length, Sum and Average
sum_int_angle = np.zeros(length_angle)
ave_int_angle = np.zeros(length_angle)
Distance_angle = np.zeros(length_angle)

for i in range(46, number_of_files_angle+46):                 #NOTE CHANGE: for i in range(1, number_of_files+1):
    df = pd.read_csv("A{}.csv".format(i))              #Read in the first file and subsequent files after looping
    int_angle = df.intensity                                #Read the Intensity column of the CSV file
    int_len_angle[i-46] = len(int_angle)                        #Find the length of the column
    sum_int_angle[i-46] = sum(int_angle)                        #Find the sum of the column
    ave_int_angle[i-46] = sum_int_angle[i-46] / int_len_angle[i-46]  #Find the average of the column
    ave_int_angle[i - 46] = ave_int_angle[i - 46] - 30                  # convert from dbm to db
    if (i == 46):
        Distance_angle[0] = 46
    else:
        Distance_angle[i-46] = Distance_angle[i-47] + 1.0  # Set the distance to start from 0 and increment up

print(Distance_angle)
fig=plt.figure()
plt.plot(Distance_angle, ave_int_angle, 'ro', label="Data Point")
plt.plot(Distance_angle, ave_int_angle, 'k', label="Antenna Pattern 76-77 GHz")
plt.legend()
plt.xlabel('Angle (Degrees)')
plt.ylabel('Intensity (dB)')
plt.title('Radar ', fontsize=14)
plt.grid()
plt.savefig('Antenna_Pattern.png')
plt.show()