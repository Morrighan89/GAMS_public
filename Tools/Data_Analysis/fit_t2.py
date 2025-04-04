import numpy as np
import matplotlib.pyplot as plt
import pydicom
from scipy.optimize import curve_fit
from scipy.io import loadmat
import os
from tqdm import tqdm

class RelaxometryFitting:
    def __init__(self, filepaths, matlab_mask_filepath,TEs, variable_name='ROIMask'):
        self.images = self.load_dicom_images(filepaths)
        self.mask = self.read_mask_from_matlab(matlab_mask_filepath, variable_name)
        self.TEs=TEs


    @staticmethod
    def load_dicom_images(filepaths):
        images = [pydicom.dcmread(filepath).pixel_array for filepath in filepaths]
        return np.array(images)

    @staticmethod
    def read_mask_from_matlab(file_path, variable_name='ROIMask'):
        mask_data = loadmat(file_path)[variable_name]
        return np.array(mask_data)

    @staticmethod
    def apply_extraction_mask(images, mask,mask_ind=1):
        # Flatten the images and mask to 1D arrays
        flat_images = images.reshape(images.shape[0],images.shape[1]*images.shape[2])
        flat_mask = mask.flatten()
        unique, counts = np.unique(flat_mask, return_counts=True)
        extracted_pixels=np.empty((images.shape[0],counts[mask_ind]))
        
        for i in range(images.shape[0]):
            extracted_pixels[i,:]=flat_images[i,flat_mask == mask_ind]
        return extracted_pixels

    @staticmethod
    def fit_function(x, *params):
        # Define a decaying exponential function
        return params[0] * np.exp(-x/params[1]) + params[2]

    def fit_pixel_across_images(self,mask_ind=1):
        num_images = len(self.images)
        pixel_values = self.apply_extraction_mask(self.images, self.mask,mask_ind)

        x_data = np.array(self.TEs)
        params=[]
        for j in range(pixel_values.shape[1]):
            # Initial guess for the parameters
            initial_params = [pixel_values[0,j], 0.1, 0]

            # Fit the decaying exponential function to the pixel values
            param, _ = curve_fit(self.fit_function, x_data, pixel_values
            [0:6,j], p0=initial_params)
            params.append(param)

        params=np.array(params)
        avg_params=np.average(params,0)
        # Plot the fitted function
        plt.errorbar(x_data, np.average(pixel_values,1),yerr=np.std(pixel_values,1), color='r', capsize=4, marker='o',linestyle='none',label='Pixel Values')
        x_forplot=np.linspace(np.min(x_data),np.max(x_data),100)
        plt.plot(x_forplot, self.fit_function(x_forplot, *avg_params), 'b-', label=f'Fitted Curve {avg_params[1]*1000:.2f}')
        plt.xlabel('Pixel Index')
        plt.ylabel('Pixel Value')
        plt.legend()
        plt.show()
        return avg_params[1]*1000

