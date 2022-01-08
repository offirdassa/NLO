from eddington import fitting_function
import numpy as np

@fitting_function(n=4,
                  name="NLO")
def super(a, x):
    ##constats:
    x=x*np.pi/180
    wavelength=0.805 ##micrometers
    n_0_w = 1.6603
    n_e_w = 1.5454
    n_0_2w=1.6924
    n_e_2w=1.5676

    theta_0=np.arcsin(np.sqrt((1/n_0_w**2-1/n_0_2w**2)/(1/n_e_2w**2-1/n_0_2w**2))) ## theta(0)
    ##
    theta=theta_0+np.arcsin(np.sin(x)/n_0_w) ## theta(alpha)
    l_alpha=1/(np.sqrt(1-(np.sin(x))**2/(n_0_w**2))) ## L(alpha)
    dk_alpha=(2/wavelength)*(n_0_w - np.sqrt((n_e_2w**2 * n_0_2w**2)/(n_0_2w**2 * (np.sin(theta))**2 + n_e_2w**2 * (np.cos(theta)**2)))) ## the pi is included inside of np.sinc
    return a[0] + a[3]*(np.sinc(a[1]*dk_alpha*l_alpha+a[2]))**2 ## y axis

