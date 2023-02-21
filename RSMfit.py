"""
data fitting on RSM from Rigaku XRD
"""

import numpy as np
from scipy.optimize import curve_fit
from scipy.signal import convolve
import matplotlib.pyplot as plt
import click

def gaussian(x, y, x0, y0, xalpha, yalpha, A):
    return A * np.exp( -((x-x0)/xalpha)**2 -((y-y0)/yalpha)**2)

def lorentzian(x, y, x0, y0, w, A):
    return A * (w**2) / ((x-x0)**2 + (y-y0)**2 + (w)**2)

def voigt(x, y, x0, y0, xalpha, yalpha, Ag, w, Al):
    return gaussian(x, y, x0, y0, xalpha, yalpha, Ag) + lorentzian(x, y, x0, y0, w, Al)

def _voigt(M, *args):
    x, y = M
    arr = np.zeros(x.shape)
    arr += voigt(x, y, *args[0:7])
    return arr
 
def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

def gauss_kern(size, sizey=None):
    """ Returns a normalized 2D gauss kernel array for convolutions """
    size = int(size)
    if not sizey:
        sizey = size
    else:
        sizey = int(sizey)
    x, y = np.mgrid[-size:size+1, -sizey:sizey+1]
    g = np.exp(-(x**2/float(size)+y**2/float(sizey)))
    return g / g.sum()

def blur_image(im, n, ny=None) :
    """ blurs the image by convolving with a gaussian kernel of typical
        size n. The optional keyword argument ny allows for a different
        size in the y direction.
    """
    g = gauss_kern(n, sizey=ny)
    improc = convolve(im,g, mode='same')
    return(improc)

@click.command()
@click.option(
    "--input-file",
    "-input",
    prompt=True,
    help=(
        "RMS txt file from Rigaku"
    )
)
def handle_input(input_file):
    #create empty list for Omega, 2theta and Intensity
    w = []
    t2 = []
    I = []

    # Build array of lines from file, strip newlines
    mylines = []                                                    # Declare an empty list.
    with open (input_file, 'rt') as myfile:    # Open txt doc for reading text.
        for myline in myfile:                                       # For each line in the file,
            mylines.append(myline.rstrip('\n'))                     # strip newline and add to list.

    #Extract w
    for line in mylines:
        if line[:26] == "*MEAS_COND_AXIS_POSITION-0":
            w.append(float(line[28:-1]))

    #Locate all lines' indexes above all the t2 and I data sections
    lineAbove = []
    index = 0
    for line in mylines:
        if line == "#Attenuator_coefficient=1.0000":
            lineAbove.append(index)
        index += 1

    #Find the index of the line below the first t2 and I data section and use that to calculate the length of each data section
    firstLineBelow = mylines[1:].index("*FILE_COMMENT \"\"")
    t2Length = firstLineBelow - lineAbove[0] 

    #Extract t2
    for line in mylines[lineAbove[0]+1 : firstLineBelow+1]:
        i = line.index(" ")
        t2.append(float(line[:i]))
        
    #Extract I - 2D Array
    for lineIndex in lineAbove:
        Ii = []
        for line in mylines[lineIndex + 1 : lineIndex + t2Length + 1]:
            i = line.index(" ")
            Ii.append(float(line[i+1:]))
        I.append(Ii)

    #Reorder w and I, because w are not organized in order in the text file 
    w, I = zip(*sorted(zip(w, I)))

    #Convert the data lists into arrays
    y = np.asarray(w)
    x = np.asarray(t2)
    Z = np.asarray(I)

    # The two-dimensional domain of the fit.
    X, Y = np.meshgrid(x, y)
    
    # Plot the 3D figure of the raw data.
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    ax.plot_surface(X, Y, Z, cmap='plasma')
    ax.set_zlim(0,np.max(Z)+2)
    plt.show()
    
    #smoothing the raw data and plot
    blur_Z = blur_image(Z, 3)
    """
    print(Z)
    print(type(Z))
    print(len(Z))
    print(len(Z[0]))
    print(blur_Z)
    print(type(blur_Z))
    print(len(blur_Z))
    print(len(blur_Z[0]))
    """
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    ax.plot_surface(X, Y, blur_Z, cmap='plasma')
    ax.set_zlim(0,np.max(blur_Z)+2)
    plt.show()
    
    # Initial guesses to the fit parameters. 
    xFWHM = abs(x[np.where(Z == np.max(Z))[1][0]] - 
                x[find_nearest(Z[np.where(Z == np.max(Z))[0][0]], np.max(Z)/2)])
    xSigma = xFWHM / 2.355

    curve_y = []

    for i in Z:
        curve_y.append(i[np.where(Z == np.max(Z))[1][0]])

    yFWHM = abs(y[np.where(Z == np.max(Z))[0][0]] - 
                y[find_nearest(curve_y, np.max(Z)/2)])
    ySigma = yFWHM / 2.355

    guess_prms = [x[int(len(x)/2)],y[int(len(y)/2)],xSigma,ySigma,np.max(Z)/2,((xFWHM+yFWHM)/2),np.max(Z)/2]

    # We need to ravel the meshgrids of X, Y points to a pair of 1-D arrays.
    xdata = np.vstack((X.ravel(), Y.ravel()))

    # Do the fit, using our custom _gaussian function which understands our
    # flattened (ravelled) ordering of the data points.
    popt, pcov = curve_fit(_voigt, xdata, Z.ravel(), p0=guess_prms, maxfev=100000)
    fit = np.zeros(Z.shape)
    fit += voigt(X, Y, *popt[0:7])
    print('Fitted parameters:')
    print(popt)
    
    #RMS residual
    rms = np.sqrt(np.mean((Z - fit)**2))
    print('RMS residual =', rms)

    # Plot the 3D figure of the fitted function.
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    ax.plot_surface(X, Y, fit, cmap='plasma')
    cset = ax.contourf(X, Y, Z-fit, zdir='z', offset=-4, cmap='plasma')
    ax.set_zlim(-4,np.max(fit))
    plt.show()

    # Plot the 3D figure of the residuals
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    ax.plot_surface(X, Y, Z-fit, cmap='plasma')
    plt.show()

    # Plot the test data as a 2D image and the fit as overlaid contours.
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.imshow(Z, origin='lower', cmap='plasma',
              extent=(x.min(), x.max(), y.min(), y.max()))
    ax.contour(X, Y, fit, colors='w')
    plt.show()

    print('Peak location is at 2theta = ', popt[0] , 'and w = ' , popt[1] )
    Q_o = (np.sin(popt[1]) + np.sin(popt[0] - popt[1] )) / (1.540593)
    Q_i = (np.cos(popt[1]) - np.cos(popt[0] - popt[1] )) / (1.540593)
    print('Out of plane Q is', Q_o , 'In plane Q is', Q_i)

    ###Fitting data after smoothing
    # Initial guesses to the fit parameters. 
    xFWHM = abs(x[np.where(blur_Z == np.max(blur_Z))[1][0]] - 
                x[find_nearest(blur_Z[np.where(blur_Z == np.max(blur_Z))[0][0]], np.max(blur_Z)/2)])
    xSigma = xFWHM / 2.355

    curve_y = []

    for i in blur_Z:
        curve_y.append(i[np.where(blur_Z == np.max(blur_Z))[1][0]])

    yFWHM = abs(y[np.where(blur_Z == np.max(blur_Z))[0][0]] - 
                y[find_nearest(curve_y, np.max(blur_Z)/2)])
    ySigma = yFWHM / 2.355

    guess_prms = [x[int(len(x)/2)],y[int(len(y)/2)],xSigma,ySigma,np.max(blur_Z)/2,((xFWHM+yFWHM)/2),np.max(blur_Z)/2]

    # We need to ravel the meshgrids of X, Y points to a pair of 1-D arrays.
    xdata = np.vstack((X.ravel(), Y.ravel()))

    # Do the fit, using our custom _gaussian function which understands our
    # flattened (ravelled) ordering of the data points.
    popt, pcov = curve_fit(_voigt, xdata, blur_Z.ravel(), p0=guess_prms, maxfev=100000)
    fit = np.zeros(blur_Z.shape)
    fit += voigt(X, Y, *popt[0:7])
    print('Fitted parameters:')
    print(popt)
    
    #RMS residual
    rms = np.sqrt(np.mean((blur_Z - fit)**2))
    print('RMS residual =', rms)

    # Plot the 3D figure of the fitted function.
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    ax.plot_surface(X, Y, fit, cmap='plasma')
    cset = ax.contourf(X, Y, blur_Z-fit, zdir='z', offset=-4, cmap='plasma')
    ax.set_zlim(-4,np.max(fit))
    plt.show()

    # Plot the 3D figure of the residuals
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    ax.plot_surface(X, Y, blur_Z-fit, cmap='plasma')
    plt.show()

    # Plot the test data as a 2D image and the fit as overlaid contours.
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.imshow(blur_Z, origin='lower', cmap='plasma',
              extent=(x.min(), x.max(), y.min(), y.max()))
    ax.contour(X, Y, fit, colors='w')
    plt.show()

    print('Peak location is at 2theta = ', popt[0] , 'and w = ' , popt[1] )
    Q_o = (np.sin(popt[1]) + np.sin(popt[0] - popt[1] )) / (1.540593)
    Q_i = (np.cos(popt[1]) - np.cos(popt[0] - popt[1] )) / (1.540593)
    print('Out of plane Q is', Q_o , 'In plane Q is', Q_i)
   
    ###END fiting data after smoothing 
    
if __name__ == "__main__":
    handle_input()







