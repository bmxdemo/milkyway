import bmxobs as bmx
import matplotlib.pyplot as plt 
import numpy as np
import statistics as stats

def bin_data(bmx_obj, channel, bin_fact, axis = "Time"):
    data = np.transpose(bmx_obj[int(channel)])
    data_shape = np.shape(data)
    
    if axis == "Time":
        ubin_len = float(data_shape[1])
    
    if axis == "Frequency":
        ubin_len = float(data_shape[0])
    
    if axis != "Time" and axis != "Frequency": 
        print("ERROR: The allowed options for 'axis' are 'Time' and 'Frequency'")
        quit()

    bin_len = ubin_len/float(bin_fact)

    rem = (ubin_len-1.0)%(bin_len-1.0)
    n = (ubin_len-1.0-rem)/(bin_len-1.0)
    bin_len = int(((ubin_len-1.0)/n)+1)
    
    new_dat_list = np.zeros(bin_len)
    

    if axis == "Time":
        data_out = np.zeros(shape = (data_shape[0], bin_len))
        for row_index in np.arange(data_shape[0]):
            for element in np.arange(bin_len):
                if element == 0:
                    new_dat_list[element] = stats.mean(data[row_index, 0:int((n+1)/2)])
                    
                if element == bin_len-1:
                    new_dat_list[element] = stats.mean(data[row_index, int(n*element - (n-1)/2):int(n*element + 1)])
                    
                if element != 0 and element != bin_len-1:
                    new_dat_list[element] = stats.mean(data[row_index, int(n*element-(n-1)/2):int(n*element+(n+1)/2)])

            data_out[row_index, :] = new_dat_list

            
    
    if axis == "Frequency":
        data_out = np.zeros(shape = (bin_len, data_shape[1]))
        for col_index in np.arange(data_shape[1]):
            for element in np.arange(bin_len):
                if element == 0:
                    new_dat_list[element] = stats.mean(data[0:int((n+1)/2),col_index])
                    
                if element == bin_len-1:
                    new_dat_list[element] = stats.mean(data[int(n*element - (n-1)/2):int(n*element + 1),col_index])
                    
                if element != 0 and element != bin_len-1:
                    new_dat_list[element] = stats.mean(data[int(n*element-(n-1)/2):int(n*element+(n+1)/2),col_index])        

            data_out[:, col_index] = new_dat_list 

    
    return data_out

data = bmx.BMXObs("/home/chandrahas/Desktop/BMX/pas/210903_0000", channels="111")
data_out = bin_data(data,"111",40,"Frequency")
print(np.shape(data_out))

################### PLOT ##########################
extent = [0,1,data.freq[1][0],data.freq[1][-1]]
im = plt.imshow(data_out, cmap = "Blues", extent = extent, aspect = 'auto', interpolation = None)
plt.colorbar()
plt.show()                    







    


