# run stochastic fits
import pandas as pd
import numpy as np
from helpers import lupus_system_stoch, fit_stochastic, residual_stoch, g_stoch


new_data = pd.read_csv("data/upc_data.csv")
# I'm willing to live without negative days
new_data = new_data[new_data['day'] >= 0]
new_data = new_data[new_data['day'] <= 500]
data_list = []
for eye_d in new_data['ID'].unique():
    data_list.append(new_data[new_data['ID'] == eye_d].sort_values("day"))

stochastic_dict = dict.fromkeys(new_data['ID'].unique())

for dataset in data_list:
    fig, ax = plt.subplots(1, 1)
    res = fit_stochastic(dataset, ax)
    stochastic_dict[dataset['ID'].to_list()[0]] = (res, fig, ax)
    print(stochastic_dict)


np.save('stochastic_dict.npy', stochastic_dict) 

# Load
#read_dictionary = np.load('my_file.npy',allow_pickle='TRUE').item()
#print(read_dictionary['hello']) # displays "world"

