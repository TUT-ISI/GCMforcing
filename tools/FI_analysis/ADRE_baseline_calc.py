import numpy as np

my_values_ADRE = np.array((-1.4,-1.5,-1.8,-1.5,-1.6,-1.6,-1.6,-1.7,-1.6,-1.4,-1.3,-1.4))

my_values_AOD = np.array((0.113, 0.119, 0.130, 0.117, 0.127, 0.123, 0.125, 0.125, 0.121, 0.113, 0.105, 0.110))

mean_value_ADRE = np.mean(my_values_ADRE)
mean_value_AOD = np.mean(my_values_AOD)

print(mean_value_ADRE)
print(mean_value_AOD)
