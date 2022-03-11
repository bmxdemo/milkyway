import os
import bmxobs as bmx
import matplotlib.pyplot as plt 
import numpy as np
import math
import statistics as stats
from vlsr import vlsr 

def correct_doppler(bmx_obj, channel, sub_div, orig_size = True, data_from_source = True, data_array = None):

    #Define function to calculate doppler factors
    def doppler_fact(rel_vel):
        beta = rel_vel/299792.458
        return np.sqrt((1+beta)/(1-beta))

    def correct_ss(sub_div, doppler_facts): #Correct using subsampling

        def bin_freq(ubin_data, data_shape, ubin_freqsize, bin_freqsize):
            #Calculate number of data points to be averaged
            n_avg = (int(ubin_freqsize)-1)/(int(bin_freqsize)-1)

            data_orig = np.zeros(shape = data_shape)

            for time in np.arange(data_shape[1]):
                for index in np.arange(bin_freqsize):
                    if index == 0:
                        data_orig[index,time] = stats.mean(ubin_data[0:int((n_avg+1)/2),time])
                    
                    if index == bin_freqsize-1:
                        data_orig[index, time] = stats.mean(ubin_data[int(sub_div*index - (n_avg-1)/2):int(sub_div*index + 1),time])
                    
                    if index != 0 and index != bin_freqsize-1:
                        data_orig[index,time] = stats.mean(ubin_data[int(sub_div*index-(n_avg-1)/2):int(sub_div*index+(n_avg+1)/2), time])
                
                os.system('clear')
                print(100.0*float(time)/float(data_shape[1]),"% doppler correction post-processing done")
            return data_orig

        #Extract data array from bmx object
        if data_from_source:
            data = np.transpose(bmx_obj[int(channel)])
        
        if not data_from_source:
            data = data_array

        data_shape = np.shape(data)
        freq = np.flip(np.transpose(bmx_obj.freq[int(channel[2])]))
        freq_shape = len(freq)
        freq_ssc = np.transpose(np.flip(np.linspace(freq[-1],freq[0],sub_div*(freq_shape-1)+1)))
        freq_ss = np.zeros(shape = (sub_div*(freq_shape-1)+1,1))
        delt_f = (freq[-1]-freq[0])/(sub_div*(freq_shape-1))

        #Create new data array 
        data_new = np.full((sub_div*(freq_shape-1)+1,data_shape[1]), np.NaN)

        #Create the frequency array at the source
        for time in np.arange(data_shape[1]):
            freq_ss = doppler_facts[time]*freq_ssc
            J = np.zeros(shape = (sub_div*(freq_shape-1)+1,1))
            #Algorithm to find nearest sub-sampled data point
            for i in np.arange(sub_div*(freq_shape-1)+1):
                j_plus = ((freq_ss[i]-freq[0])/delt_f)+0.5
                J[i] = int(math.modf(j_plus)[1])
              
            freq_ss = J

            #Rearrange data elements according to location given by freq_ss
            for element_numb in np.arange(sub_div*(freq_shape-1)+1):
                ss_count = element_numb%sub_div

                if freq_ss[element_numb] >= 0  and freq_ss[element_numb] <= sub_div*(freq_shape-1):
                    data_new[int(freq_ss[element_numb]),time] = data[int((element_numb - ss_count)/sub_div),time]

            #Fill zero data entries in data_new
            val = 0
            for element_numb in np.arange(sub_div*(freq_shape-1)+1):

                if np.isnan(data_new[element_numb,time]):
                    data_new[element_numb,time] = val
                    
                else:
                    val = data_new[element_numb,time]
            
            os.system('clear')
            per = 100.0*float(time)/float(data_shape[1])
            print("%.5f" % per,"% doppler correction done")

        if orig_size == True:
            data_final = bin_freq(data_new, data_shape, sub_div*(freq_shape-1)+1, freq_shape)

        return data_final if orig_size == True else data_new 
    
    #Extract velocity data
    vel = vlsr(bmx_obj.ra*180/np.pi, 40.8, bmx_obj.mjd)
    v_rad = vel[0]+vel[1]+vel[2]

    #Calculate the doppler factors for all times
    doppler_facts = doppler_fact(v_rad)

    return correct_ss(sub_div, doppler_facts)
                


                    



            








