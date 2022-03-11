import bmxobs as bmx
import matplotlib.pyplot as plt 
import numpy as np
import scipy.optimize as opt
import scipy.signal as sgn 
import os 

def linear(x,a,b):
    return a*x+b


################denoising###################
def denoise(dat_file,channel, n_pix, Time = [0,1]):
    
    print("denoising")
    #Determine the fit constants a,b in af+b for every time step
    data = np.transpose(dat_file[int(channel)])
    freq = dat_file.freq[int(channel[2])]
    x_fit_top = freq[freq.size-n_pix:freq.size]
    x_fit_bot = freq[0:n_pix]

    size = data[0].size
    fit_array = np.zeros(shape = (4,size))
    for iter in np.arange(size):
        y_fit_top = np.flip(np.transpose(data[0:n_pix, iter]))
        y_fit_bot = np.flip(np.transpose(data[np.shape(data)[0]-n_pix:np.shape(data)[0], iter]))
        const_out_top, junk = opt.curve_fit(linear, x_fit_top, y_fit_top)
        const_out_bot, junk = opt.curve_fit(linear,x_fit_bot,y_fit_bot)

        const_out_top = np.transpose(const_out_top)
        const_out_bot = np.transpose(const_out_bot)
        fit_array[:,iter] = np.concatenate((const_out_top,const_out_bot), axis = 0)

        os.system('clear')
        per = 100*float(iter)/float(size)
        print("%.5f" % per, "% of denoising done")
        

    avg_fit_array = np.empty(shape = (2,data[0].size))
    avg_fit_array[0,:] = (fit_array[0,:]+fit_array[2,:])/2.0
    avg_fit_array[1,:] = (fit_array[1,:]+fit_array[3,:])/2.0

    #Smoothen fit of (a,b) to determine gain
    Times = np.linspace(Time[0], Time[1], num = data[0,:].size, endpoint = True)

    for i in np.arange(0,2):
        avg_fit_array[i,:] = sgn.savgol_filter(avg_fit_array[i,:], 1001, 3)
        

    gain_array = np.zeros(shape = np.shape(data))

    for iter in np.arange(data[0].size):
        gain_array[:,iter] = np.flip(linear(freq, avg_fit_array[0,iter], avg_fit_array[1,iter]))
    
    ones = np.ones(shape = np.shape(data))

    denoised_out = np.divide(data,gain_array) - ones 
    return denoised_out