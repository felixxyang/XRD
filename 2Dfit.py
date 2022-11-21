"""
data fitting on 2Theta v.s. Intensity curve from Rigaku XRD
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy import optimize
import scipy.integrate as integrate
import math
import click
import sys
import json
import statistics

def _1Voigt(x, ampG1, cenG1, sigmaG1, cenL1, widL1):
    return (ampG1*(1/(sigmaG1*(np.sqrt(2*np.pi))))*(np.exp(-((x-cenG1)**2)/(2*(sigmaG1)**2)))) +\
              (((1-ampG1)*widL1**2/(math.pi*((x-cenL1)**2+widL1**2))) )  

def linearFit(x, k, b):
    return k * x + b

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]

def FWHM(X,Y):
    half_max = np.max(Y) / 2
    max_index = np.argmax(Y)
    left_nearest = find_nearest(Y[:max_index], half_max)
    right_nearest = find_nearest(Y[max_index:], half_max)
    left_index = np.where(Y[:max_index] == left_nearest)
    right_index = np.where(Y[max_index:] == right_nearest) + max_index
    return X[right_index[0][0]]-X[left_index[0][0]]

def Rsquared(Y,Yfit):
    ss_res = np.sum((Y - Yfit)**2)
    ss_tot = np.sum((Y - np.mean(Y))**2)
    R2 = 1 - (ss_res / ss_tot)
    return R2

def AdjRsquared(Y,Yfit,pars):
    ss_res = np.sum((Y - Yfit)**2)
    ss_tot = np.sum((Y - np.mean(Y))**2)
    n = len(Y)
    K = len(pars)
    AR2 = 1 - ((ss_res/ (n - K)  ) / (ss_tot/ (n - 1)  ))
    return AR2

def integration(df, pars, fwhm):

    angle_arr = np.asarray(df['angle'])
    int_arr = np.asarray(df['intensity'])
    peak = find_nearest( angle_arr , pars[3])
    peak_ind = np.where(angle_arr == peak)[0][0]
    left_cutoff = find_nearest( angle_arr[:peak_ind] , peak - 3*fwhm)
    right_cutoff = find_nearest( angle_arr[peak_ind:] , peak + 3*fwhm)
    lc_index = np.where(angle_arr[:peak_ind] == left_cutoff)[0][0]
    rc_index = np.where(angle_arr[peak_ind:] == right_cutoff)[0][0] + peak_ind
    
    pars_l, cov_l = optimize.curve_fit(f=linearFit, xdata=angle_arr[:lc_index], ydata=int_arr[:lc_index])
    pars_r, cov_r = optimize.curve_fit(f=linearFit, xdata=angle_arr[rc_index:], ydata=int_arr[rc_index:])
    ave_k = (pars_l[0] + pars_r[0]) / 2
    ave_b = (pars_l[1] + pars_r[1]) / 2
    
    for i in range(len(angle_arr)):
        int_arr[i] -= linearFit(angle_arr[i], ave_k, ave_b)  
        
    area = np.trapz(int_arr, angle_arr)
    
    return area, peak_ind, lc_index, rc_index, angle_arr, int_arr
    #result = integrate.quad(lambda x:_1Voigt(x, *pars), df['angle'].iat[0], df['angle'].iat[-1])
    #return result[0]

#spicify local xlsx file (sample growth temperatures) location here 
def findTemp(sample):
    df_temp = pd.read_excel("C:\\Users\\felix\\OneDrive\\Desktop\\Temperatures.xlsx")
    for i in range(df_temp.shape[0]):
        if df_temp['sample'][i] == sample:
            temp = df_temp['temp'][i]
    return temp

@click.command()
@click.option(
    "--input-file",
    "-input",
    prompt=True,
    help=(
        "Intensity vs 2Theta txt file from Rigaku"
    )
)
@click.option(
    "--output-file",
    "-output",
    prompt=True,
    help=(
        "Output file for data of this measurement"
    )
)

def handle_input(input_file, output_file):
    
    if sys.platform == "win32":
        sampleName = input_file.split("\\")[-1].split("_")[1]
        planeName = input_file.split("\\")[-1].split("_")[2].split(".")[0]
    else:
        sampleName = input_file.split("/")[-1].split("_")[1]
        planeName = input_file.split("/")[-1].split("_")[2].split(".")[0]
    
    df = pd.read_csv(input_file, delim_whitespace=True, header=0, names=['angle','intensity'], skiprows = 777)
    fig1, ax1 = plt.subplots()
    ax1.plot(df['angle'],df['intensity'],label="measurement")
    ax1.set_xlabel('2$\Theta$')
    ax1.set_ylabel('intensity')
    ax1.set_title(sampleName + '_' + planeName)
    #fit
    max_position = df['angle'][df['intensity'][df['intensity'] == np.max(df['intensity'])].index[0]]
    wid = FWHM(df['angle'], df['intensity'])
    #sigma = FWHM(df['angle'], df['intensity'])/2.355 
    #Guassian is turned off in Voigt to get better fit
    guess_prms = [0, 1, 1, max_position, wid]
    pars, cov = optimize.curve_fit(f=_1Voigt, xdata=df['angle'], ydata=df['intensity'],p0=guess_prms)
    ax1.plot(df['angle'], _1Voigt(df['angle'], *pars), label="fitted curve")
    #ax1.legend()
    #plt.show()
    
    #residula
    #fig2, ax2 = plt.subplots()
    #ax2.plot(df['angle'], df['intensity']-_1Voigt(df['angle'], *pars))
    #ax2.set_title('residual')
    #ax2.set_xlabel('2$\Theta$')
    #ax2.set_ylabel('difference in intensity')
    #plt.show()
    
    
    fwhm = FWHM(df['angle'], _1Voigt(df['angle'], *pars))
    print("The FWHM is", fwhm)
    
    area, p, lci, rci, new_angle, new_int = integration(df, pars, fwhm)
    print("The arae under the peak is", area)
    ax1.axvline(df['angle'][lci], color='r')
    ax1.axvline(df['angle'][rci], color='r')
    ax1.axhline(0, color='r')
    #ax1.axvline(df['angle'][p], color='orange')
    ax1.plot(new_angle, new_int, label='background cut-out')
    ax1.legend()
    plt.show()
    
    
    R2 = Rsquared( df['intensity'], _1Voigt(df['angle'], *pars))
    print("The coefficient of determination, R^2, is", R2)
    AdjustedR2 = AdjRsquared( df['intensity'], _1Voigt(df['angle'], *pars), pars)
    print("The adjusted R^2 is", AdjustedR2)
    theta = pars[3] / 2
    print("Peak postition is located at", 2 * theta)
    #Rigaku uses Cu as X-ray source
    if planeName == "CFS004":
        a = (1.540593 / (2 * np.sin(theta * np.pi / 180))) * 4
        print("Out of plane lattice Constant is", a)
    elif planeName == "CFS002":
        a = (1.540593 / (2 * np.sin(theta * np.pi / 180))) * 2
        print("Out of plane lattice Constant is", a)
    
    #finding sample's growth temperature
    temp = int(findTemp(sampleName))
    
    #create dictionary
    data_dict = {}
    data_dict[sampleName+"_"+planeName] = {
        "temp" : temp,
        "a" : a,
        "area" : area,
        "fwhm" : fwhm,
        "R2" : R2,
        "AdjustedR2": AdjustedR2
        }
    
    with open(output_file, "w") as ofs:
        json.dump(data_dict, ofs, indent=2)
       
if __name__ == "__main__":
    handle_input()
    
    

    