import inspect
import numpy as np

def list_variables(scope):
    # Get the current function name
    function_name = inspect.currentframe().f_back.f_code.co_name

    # Print the function name
    print(f"Python Function: {function_name}")

    # Print the header
    print(f"{'Variable':<15} | {'Value'}")

    # Print a separator
    print("-" * 30)

    # Iterate and print each key-value pair
    for key, value in scope.items():
        # Check if the value is a list or numpy array and has more than 5 elements
        if isinstance(value, (list, np.ndarray)) and len(value) > 5:
            # Print only the first and last few elements
            display_value = f"[{value[0]}, {value[1]}, ..., {value[-2]}, {value[-1]}]"
        else:
            # Otherwise, print the value as it is
            display_value = value

        print(f"{key:<15} | {display_value}")

import csv
import numpy as np
from datetime import datetime

def save_variables_to_csv(scope, excluded_vars=['params', 'axs'], file_path=None):
    # Create a filename with current datetime
    dt = datetime.now().strftime('%Y%m%d%H%M%S')
    if not file_path:
        file_path = f'pythonmainvariables{dt}.csv'

    with open(file_path, mode='w', newline='') as file:
        writer = csv.writer(file)

        # Write the header
        writer.writerow(['Variable', 'Value'])

        # Iterate and write each key-value pair
        for key, value in scope.items():
            if key not in excluded_vars:
                # Handle array-like structures
                if isinstance(value, (list, np.ndarray)):
                    # Convert the entire array to a string representation
                    value_str = ', '.join(map(str, value))
                else:
                    # Convert other types of values to string
                    value_str = str(value)

                writer.writerow([key, value_str])

# Usage
# Inside your function or scope where you want to debug
# save_variables_to_csv(locals())


        
import numpy as np

def rayleighdebye(diameter, lambda_, n_background, n_sphere, theta, R):
    """
    Calculate the intensity of scattered light using the Rayleigh-Debye scattering model.

    This function computes the light intensity scattered by a sphere with a given diameter,
    refractive index, and at a specific angle theta using the Rayleigh-Debye scattering theory.

    Parameters:
    diameter (float): Diameter of the scattering sphere.
    lambda_ (float): Wavelength of the incident light.
    n_background (float): Refractive index of the background medium.
    n_sphere (float): Refractive index of the sphere.
    theta (float): Scattering angle (in radians).
    R (float): Distance from the scatterer to the observation point.

    Returns:
    float: Intensity of the scattered light at the given angle.

    Notes:
    The formula used in this function is derived from the Rayleigh-Debye scattering theory,
    which is suitable for particles that are small compared to the wavelength of light.
    """

    radius = diameter / 2
    k = 2 * np.pi * n_background / lambda_
    x = k * radius
    n_rel = n_sphere / n_background  # real refractive index
    u = 2 * k * radius * np.sin(theta / 2)
    f_theta = 1 + (np.cos(theta))**2  # the base Rayleigh form factor.
    P_theta = 3 * (np.sin(u) - u * np.cos(u)) / u**3  # the form factor P(theta)
    Rr = P_theta * f_theta
    I = abs(Rr * (abs((n_rel**2 - 1) / (n_rel**2 + 2))**2) * (1 / (2 * R**2)) * (((2 * np.pi) / lambda_)**4) * (radius**6))
    # Combining equations to get the intensity.

    return I

import numpy as np


def SEALS(d, D, a, dcorr, P, NA, lamvec):
    """
    SEALS Function Description:
    Here you can include a detailed description of the function, its parameters, and what it computes.

    Parameters:
    - d: [Description of d]
    - D: [Description of D]
    - a: [Description of a]
    - dcorr: [Description of dcorr]
    - P: [Description of P]
    - NA: [Description of NA]
    - lamvec: [Description of lamvec]

    Returns:
    - y: [Description of y]
    - theta: [Description of theta]
    """

    # The computations from the MATLAB function converted to Python syntax
    y = 1/6 * D * np.tan(a - np.arcsin(lamvec / d - np.sin(a))) / (1 + np.tan(a - np.arcsin(lamvec / d - np.sin(a))) * np.tan(a))
    y = y - y[-1]
    ycenter = (y[0] - y[-1]) / 2
    theta = np.degrees(np.arctan(2 / P * (y - ycenter + dcorr) * np.tan(np.arcsin(NA))))

    # Return the computed values as a tuple
    return y, theta

import numpy as np
from scipy.special import jv, yv

def mie(npar, nmed, dia, lambda_, angl, r):
    """
    Calculate the Mie scattering parameters.

    Parameters:
    npar (float): Refractive index of the particle.
    nmed (float): Refractive index of the medium.
    dia (float): Diameter of the particle.
    lambda_ (float): Wavelength of light.
    angl (array): Array of angles for scattering.
    r (float): Distance parameter.

    Returns:
    tuple: A tuple containing the scattering parameters:
           - sigma_s (float): Scattering cross-section.
           - I_p (array): Intensity of p-polarized light.
           - I_s (array): Intensity of s-polarized light.
           - an (array): Mie coefficient an.
           - bn (array): Mie coefficient bn.
           - T_p (array): Phase of p-polarized light.
           - T_s (array): Phase of s-polarized light.
    """

    pi = np.pi
    x = pi * dia / (lambda_ / nmed)
    phi = pi  # Scattering plane
    k = 2 * pi / (lambda_ / nmed)
    A = pi * dia**2 / 4
    nmax = round(2 + x + 4 * x**(1/3))
    m = npar / nmed
    n = np.arange(1, nmax + 1)
    nu = n + 0.5
    z = m * x
    m2 = m**2

    # Spherical Bessel functions
    jx = jv(nu, x) * np.sqrt(0.5 * pi / x)
    jz = jv(nu, z) * np.sqrt(0.5 * pi / z)
    yx = yv(nu, x) * np.sqrt(0.5 * pi / x)
    hx = jx + 1j * yx

    # Lower order calculations
    j1x = np.hstack((np.sin(x) / x, jx[:-1]))
    j1z = np.hstack((np.sin(z) / z, jz[:-1]))
    y1x = np.hstack((-np.cos(x) / x, yx[:-1]))
    h1x = j1x + 1j * y1x
    ax = x * j1x - n * jx
    az = z * j1z - n * jz
    ahx = x * h1x - n * hx

    # Mie coefficients
    an = (m2 * jz * ax - jx * az) / (m2 * jz * ahx - hx * az)
    bn = (jz * ax - jx * az) / (jz * ahx - hx * az)
    print(f"an[0]: {an[0]}, bn[0]: {bn[0]}")  # Print first element for checking

    # Scattering calculations
    x2 = x * x
    anp = np.real(an)
    anpp = np.imag(an)
    bnp = np.real(bn)
    bnpp = np.imag(bn)
    en = (2 * n + 1) * (anp**2 + anpp**2 + bnp**2 + bnpp**2)
    q = np.sum(en)
    qsca = 2 * q / x2

    # Electric field calculations
    E_theta = np.zeros(500)
    E_phi = np.zeros(500)
    for a in range(500):
        #print(f"500: Loop iteration (Python): {a}")
        theta = angl[a]
        u = np.cos(theta)
        p = np.zeros(nmax)
        t = np.zeros(nmax)
        p[0] = 1  # Calculation of spherical harmonic functions
        t[0] = u
        p[1] = 3 * u
        t[1] = 3 * np.cos(2 * np.arccos(u))
        for j in range(3, nmax):
            #print(f"Current value of j: {j}")
            p1 = (2 * j - 1) / (j - 1) * p[j - 1] * u
            p2 = j / (j - 1) * p[j - 2]
            p[j] = p1 - p2

            t1 = j * u * p[j]
            t2 = (j + 1) * p[j - 1]
            t[j] = t1 - t2

        n2 = (2 * n + 1) / (n * (n + 1))
        pin = n2 * p[:nmax] # pin = n2.*p;
        tin = n2 * t[:nmax] # tin = n2.*t;
        S1 = an * pin + bn * tin # calculation process for the electric field
        S2 = an * tin + bn * pin
        S1_c = np.sum(S1) # scattering amplitude
        S2_c = np.sum(S2)
        E_theta[a] = np.exp(1j * k * r) / (-1j * k * r) * np.cos(phi) * S2_c # parallel
        E_phi[a] = np.exp(1j * k * r) / (1j * k * r) * np.cos(phi) * S1_c # perpendicular

    T_p = np.angle(E_theta) # finding phase
    T_s = np.angle(E_phi)
    sigma_s = qsca * A
    I_p = np.abs(E_theta)**2 # finding intensity
    I_s = np.abs(E_phi)**2
    I_tot = I_p + I_s
    print("Finished Mie function in Python")
    return sigma_s, I_p, I_s, an, bn, T_p, T_s

import numpy as np
import matplotlib.pyplot as plt

def main():
    """
    Main function to handle the flow of calculations and plotting for SEALS, Mie, and Rayleigh-Debye models.
    It uses default values for various parameters, computes scattering profiles, and generates plots.
    """
    # Define default parameters
    params = {
        "dia": 9940e-9,  # diameter in meters
        "npar": 1.39,    # refractive index of sphere
        "nmed": 1,       # refractive index of medium
        "r": 0.1,        # distance of measurement device from particle
        "d": 9.0909e-7,  # groove spacing
        "D": 0.065,      # distance between gratings
        "a": 0.9023,     # grating tilt angle
        "dcorr": -4.2448e-4, # lens correctional term
        "P": 0.0058,     # diameter of lens
        "NA": 0.7,       # numerical aperture of lens <gpt> what is a numerical aperature? <gpt>
        "mangle": 20,    # measurement angle
        "lam1": 1580e-9, # SEALS laser lower wavelength limit
        "lam2": 1600e-9, # SEALS laser upper wavelength limit
        # lambda_ and lamvec will be calculated below
    }

    # Calculate lamvec based on lam1 and lam2 before calling SEALS
    params["lambda_"] = (params["lam1"] + params["lam2"]) / 2
    params["lamvec"] = np.linspace(params["lam1"], params["lam2"], 500)



    # Now you can call SEALS with the parameters from the dictionary
    y, theta = SEALS(params["d"], params["D"], params["a"], params["dcorr"], params["P"], params["NA"], params["lamvec"])

    # Call SEALS function to get y and theta
    y, theta = SEALS(params["d"], params["D"], params["a"], params["dcorr"], params["P"], params["NA"], params["lamvec"])

    # Call the mie function to get relevant variables
    sigma_s, I_p, I_s, an, bn, T_p, T_s = mie(params["npar"], params["nmed"], params["dia"], params["lambda_"], theta, params["r"])

    # Similarly, Call the rayleighdebye function to get I
    I = rayleighdebye(params["dia"], params["lambda_"], params["nmed"], params["npar"], theta, params["r"])

    # Define matplotlib plot size
    figx, figy = 8,8

    # # Plotting Section
    # # Plot 1: Vertical Displacement vs. Wavelength
    # plt.figure(1,figsize=(figx, figy))
    # plt.plot(params["lamvec"], y)
    # plt.xlabel('Wavelength (m)')
    # plt.ylabel('Vertical Displacement (m)')
    # plt.title('Vertical Displacement vs. Wavelength')
    # plt.show()

    # # Plot 2: Scattering Angle vs. Wavelength
    # plt.figure(2,figsize=(figx, figy))
    # theta_adj = theta + params["mangle"]  # adjust theta by the measurement angle
    # plt.plot(params["lamvec"], theta_adj)
    # plt.title('Scattering Angle vs. Wavelength')
    # plt.xlabel('Wavelength (m)')
    # plt.ylabel('Scattering Angle (deg)')
    # plt.show()

    fig, axs = plt.subplots(1, 2, figsize=(10, 5))  # Creating a 1x2 grid

    # Plot 1: Vertical Displacement vs. Wavelength
    axs[0].plot(params["lamvec"], y)
    axs[0].set_title('Vertical Displacement vs. Wavelength')
    axs[0].set_xlabel('Wavelength (m)')
    axs[0].set_ylabel('Vertical Displacement (m)')

    # Plot 2: Scattering Angle vs. Wavelength
    theta_adj = theta + params["mangle"]  # adjust theta by the measurement angle
    axs[1].plot(params["lamvec"], theta_adj)
    axs[1].set_title('Scattering Angle vs. Wavelength')
    axs[1].set_xlabel('Wavelength (m)')
    axs[1].set_ylabel('Scattering Angle (deg)')

    plt.tight_layout()
    plt.show()

    band = 20 * 10**-9
    c = 3 * 10**8
    #lambda_ = (1580 * 10**-9 + 1600 * 10**-9) / 2  # assuming lambda is the average of lam1 and lam2

    vband = c / (params["lambda_"]**2) * band
    #vvec = np.linspace(c / (1600 * 10**-9), c / (1580 * 10**-9), 500)  # converting lam1 and lam2
    vvec = np.linspace(c/params["lam2"], c/params["lam1"],500)
    lineshape = vband / (2 * np.pi * ((vvec - c / params["lambda_"])**2 + (vband / 2)**2))
    lineshape = lineshape / np.max(lineshape)  # normalizing i.e taking into account the laser intensity distribution.

    # # Mie Scattering Plots  3 - 6
    # # Plot 3: Intensity vs. Scattering Angle for Mie Scattering
    # plt.figure(3,figsize=(figx, figy))
    # I_tot = (I_p + I_s)
    # plt.plot(theta, 10 * np.log10(I_tot))  ## I_tot = I_p + I_s
    # plt.title('Intensity vs. Scattering Angle for Mie Scattering')
    # plt.xlabel('Scattering Angle (deg)')
    # plt.ylabel('Intensity (dB)')
    # plt.grid(True)
    # plt.show()

    # # Plot 4: Intensity (log scale) vs. Scattering Angle for Mie Scattering
    # plt.figure(4,figsize=(figx, figy))
    # ax = plt.subplot(111, polar=True)
    # ax.plot(np.deg2rad(theta), np.log10(I_p + I_s) - np.min(np.log10(I_p + I_s)))
    # plt.title('Intensity (log scale) vs. Scattering Angle for Mie Scattering')
    # plt.show()

    # # Plot 5: Intensity vs. Wavelength for Mie Scattering
    # plt.figure(5,figsize=(figx, figy))
    # Im = (I_p + I_s) * lineshape  # Assuming lineshape is defined
    # plt.plot(params["lamvec"], 10 * np.log10(Im))
    # plt.title('Intensity vs. Wavelength for Mie Scattering')
    # plt.xlabel('Wavelength (m)')
    # plt.ylabel('Intensity (dB)')
    # plt.grid(True)
    # plt.show()

    # # Plot 6: Phase vs. Scattering Angle for Mie Scattering
    # plt.figure(6,figsize=(figx, figy))
    # plt.plot(params["lamvec"], np.rad2deg(T_p), label='P-polarized')
    # plt.plot(params["lamvec"], np.rad2deg(T_s), label='S-polarized')
    # plt.title('Phase vs. Scattering Angle for Mie Scattering')
    # plt.xlabel('Wavelength (m)')
    # plt.ylabel('Phase (deg)')
    # plt.legend()
    # plt.grid(True)
    # plt.show()


    fig, axs = plt.subplots(1, 4, figsize=(20, 5))  # Creating a 1x4 grid

    # Plot 3: Intensity vs. Scattering Angle for Mie Scattering
    I_tot = (I_p + I_s)
    axs[0].plot(theta, 10 * np.log10(I_tot))
    axs[0].set_title('Intensity vs. Scattering Angle for Mie Scattering')
    axs[0].set_xlabel('Scattering Angle (deg)')
    axs[0].set_ylabel('Intensity (dB)')

    # Plot 4: Intensity (log scale) vs. Scattering Angle for Mie Scattering (Polar plot)
    axs[1] = plt.subplot(1, 4, 2, polar=True)
    axs[1].plot(np.deg2rad(theta), np.log10(I_p + I_s) - np.min(np.log10(I_p + I_s)))
    axs[1].set_title('Intensity (log scale) vs. Scattering Angle for Mie Scattering')

    # Plot 5: Intensity vs. Wavelength for Mie Scattering
    Im = (I_p + I_s) * lineshape  # Assuming lineshape is defined
    axs[2].plot(params["lamvec"], 10 * np.log10(Im))
    axs[2].set_title('Intensity vs. Wavelength for Mie Scattering')
    axs[2].set_xlabel('Wavelength (m)')
    axs[2].set_ylabel('Intensity (dB)')

    # Plot 6: Phase vs. Scattering Angle for Mie Scattering
    axs[3].plot(params["lamvec"], np.rad2deg(T_p), label='P-polarized')
    axs[3].plot(params["lamvec"], np.rad2deg(T_s), label='S-polarized')
    axs[3].set_title('Phase vs. Scattering Angle for Mie Scattering')
    axs[3].set_xlabel('Wavelength (m)')
    axs[3].set_ylabel('Phase (deg)')
    axs[3].legend()

    plt.tight_layout()
    plt.show()


    # # Rayleigh-Debye Scattering Plots:

    # # Plot 7: Intensity vs. Scattering Angle for Rayleigh Scattering
    # plt.figure(7,figsize=(figx, figy))
    # plt.plot(theta, 10 * np.log10(I))
    # plt.title('Intensity vs. Scattering Angle for Rayleigh Scattering')
    # plt.xlabel('Scattering Angle (deg)')
    # plt.ylabel('Intensity (dB)')
    # plt.grid(True)
    # plt.show()

    # # Plot 8: Intensity (log scale) vs. Scattering Angle for Rayleigh Scattering
    # plt.figure(8,figsize=(figx, figy))
    # ax = plt.subplot(111, polar=True)
    # ax.plot(np.deg2rad(theta), np.log10(I) - np.min(np.log10(I)))
    # plt.title('Intensity (log scale) vs. Scattering Angle for Rayleigh Scattering')
    # plt.show()

    # # Plot 9: Intensity vs. Wavelength for Rayleigh Scattering
    # plt.figure(9,figsize=(figx, figy))
    # plt.plot(params["lamvec"], 10 * np.log10(I))
    # plt.title('Intensity vs. Wavelength for Rayleigh Scattering')
    # plt.xlabel('Wavelength(m)')
    # plt.ylabel('Intensity (dB)')
    # plt.grid(True)
    # plt.show()

    fig, axs = plt.subplots(1, 3, figsize=(15, 5))  # Creating a 1x3 grid

    # Plot 7: Intensity vs. Scattering Angle for Rayleigh Scattering
    axs[0].plot(theta, 10 * np.log10(I))
    axs[0].set_title('Intensity vs. Scattering Angle for Rayleigh Scattering')
    axs[0].set_xlabel('Scattering Angle (deg)')
    axs[0].set_ylabel('Intensity (dB)')

    # Plot 8: Intensity (log scale) vs. Scattering Angle for Rayleigh Scattering (Polar plot)
    axs[1] = plt.subplot(1, 3, 2, polar=True)
    axs[1].plot(np.deg2rad(theta), np.log10(I) - np.min(np.log10(I)))
    axs[1].set_title('Intensity (log scale) vs. Scattering Angle for Rayleigh Scattering')

    # Plot 9: Intensity vs. Wavelength for Rayleigh Scattering
    axs[2].plot(params["lamvec"], 10 * np.log10(I))
    axs[2].set_title('Intensity vs. Wavelength for Rayleigh Scattering')
    axs[2].set_xlabel('Wavelength (m)')
    axs[2].set_ylabel('Intensity (dB)')

    plt.tight_layout()
    plt.show()
    file_directory = '/Users/colincasey/UCLA/jalalilab/jalali-lab-first-assignment/consoleTXToutputs'
    file_name = 'pythonmainvariables.csv'  # You can choose any file name you like
    full_file_path = f"{file_directory}/{file_name}"

    save_variables_to_csv(locals(), file_path=full_file_path)
    list_variables(locals())
    
    # save_variables_to_csv(locals(), file_path='/Users/colincasey/UCLA/jalalilab/jalali-lab-first-assignment/consoleTXToutputs')
    list_variables(locals())


if __name__ == "__main__":
    main()