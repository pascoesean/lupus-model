# run non-stochastic fits
# this is the template; can be run on its own or can be split into data chunks and 
# each can be run separately (embarassingly parallel)

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from helpers import  fit_fixed
from tqdm import tqdm

new_data = pd.read_csv("data/upc_data.csv")
# I'm willing to live without negative days
new_data = new_data[new_data['day'] >= 0]
new_data = new_data[new_data['day'] <= 500]
data_list = []

eye_d_list=  new_data['ID'].unique()[0:5]

for eye_d in eye_d_list:
    data_list.append(new_data[new_data['ID'] == eye_d].sort_values("day"))


fixed_dict = dict.fromkeys(eye_d_list)

for dataset in tqdm(data_list):
    fig, ax = plt.subplots(1, 1)
    res = fit_fixed(dataset, ax)
    fixed_dict[dataset['ID'].to_list()[0]] = (res, fig, ax)
    print(fixed_dict)


np.save('data/saved_dictionaries/fixed_dict.npy', fixed_dict) 

# Load
#read_dictionary = np.load('my_file.npy',allow_pickle='TRUE').item()
#print(read_dictionary['hello']) # displays "world"

