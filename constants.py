import numpy as np
from dataclasses import dataclass


# Constants (all in SI units)
k_B = 1.380649e-23  # Boltzmann constant in J/K
hbar = 1.054572e-34  # reduced Planck's constant
h = hbar * 2 * np.pi  # Planck constant
amu = 1.660539066e-27  # 1 amu in kg
c = 299792458  # speed of light in m/s
mu_B = 9.2740100657e-24  # Bohr magneton in J/T
a_0 = 5.29177210544e-11  # Bohr radius in m
q_e = 1.602176634-19  # electron charge
mu_0 = 1.25663706127e-6  # vacuum permeability
epsilon_0 = 8.8541878188e-12  # vacuum permittivity

# atomic properties


class Cs133:
    """TODO!
    This class stores atomic properties of cesium-133."""
    mass = 132.905451931 * amu  # atomic mass in kg
    I = 7/2  # nuclear spin
    g_I = -0.00039885395  # nuclear g-factor

    class GS:
        g_J = 2.00254032  # Lande g-factor
        a_HF = 2298157.9425  # hyperfine constant in Hz

    class D2:
        omega = 2 * np.pi * 351.72571850e12  # angular frequency in Hz
        Gamma = 2 * np.pi * 5.234e6  # decay rate in Hz
        wavelength = 852.34727582e-9  # wavelength in m
        a_HF = np.nan  # hyperfine constant in Hz

    class D1:
        omega = 2 * np.pi * 335.116048807e12  # angular frequency in Hz
        Gamma = 2 * np.pi * 4.575e6  # decay rate in Hz
        wavelength = 894.59295987e-9  # wavelength in m


class Li6:
    """TODO!
    This class stores atomic properties of lithium-6."""
    # mass = 132.905451931 * amu  # atomic mass in kg
    # I = 7/2  # nuclear spin
    # g_I = -0.00039885395  # nuclear g-factor

    # class GS:
    #     g_J = 2.00254032  # Lande g-factor
    #     a_HF = 2298157.9425  # hyperfine constant in Hz

    # class D2:
    #     omega = 2 * np.pi * 351.72571850e12  # angular frequency in Hz
    #     Gamma = 2 * np.pi * 5.234e6 # decay rate in Hz
    #     wavelength = 852.34727582e-9  # wavelength in m
    #     a_HF = np.nan # hyperfine constant in Hz

    # class D1:
    #     omega = 2 * np.pi * 335.116048807e12 # angular frequency in Hz
    #     Gamma = 2 * np.pi * 4.575e6 # decay rate in Hz
    #     wavelength = 894.59295987e-9 # wavelength in m


def breit_rabi(field, atom="cs"):
    """TODO!
    This function uses the Breit-Rabi formula to get the energy difference between the highest and lowest 
    hyperfine levels of an atom in a magnetic field. 

    Parameters:
    field (float or array-like): Magnetic field values in Gauss.
    atom_type (int): The type of atom (6 for lithium, 133 for cesium).

    Returns:
    float or numpy.ndarray: The energy difference in MHz.
    """

    h = 6.62606957E-34
    u_B_si = 9.27400915E-24

    # in MHz/G
    u_B = u_B_si * 1E-4 * 1E-6 / h
    out_energy = u_B * field

    if atom_type == 6:  # lithium
        g_J = 2.0023010
        g_I = -0.0004476540
        a_HF = 152.1368407
        I = 1
    elif atom_type == 133:  # cesium
        g_J = 2.00254032
        g_I = -0.00039885395
        a_HF = 2298.1579425
        I = 7 / 2
    elif atom_type == 40:  # potassium
        g_J = 2.00229421
        g_I = 0.000176490
        a_HF = -285.7308
        I = 4
    else:
        raise ValueError("Invalid atom type. Must be 6, 133, or 40.")

    x = (g_J - g_I) * u_B / (a_HF * (I + 1 / 2)) * field

    temp1 = a_HF * (I + 1 / 2) / 2
    temp2 = 1 + 4 * m_f * x / (2 * I + 1) + x ** 2

    branch = branch * np.ones_like(x)
    if m_f == -(I + 1 / 2) * np.sign(a_HF):
        funny_sign = x * np.sign(a_HF) > 1
        branch[funny_sign] = branch[funny_sign] * -1

    out_energy = -a_HF / 4 + g_I * u_B * m_f * \
        field + temp1 * branch * np.sqrt(temp2)

    return out_energy




def aCsCs(B):
    """Cs-Cs scattering length as a function of magnetic field.
    
    This function calculates the scattering length of two Cs atoms in the 
    |3, 3> state as a function of magnetic field. This function interpolates 
    data from Berninger et al. "Feshbach resonances, weakly bound molecular 
    states, and coupled-channel potentials for cesium at high magnetic fields" 
    Phys. Rev. A 87, 032517 (2013). The data is stored in 'Cs_a_vs_B.txt'.
    
    Parameters
    ----------
    B : float or array-like
        Magnetic field values in Gauss.

    Returns
    -------
    a : float or array-like, same shape as B
        Scattering lengths in units of Bohr radii corresponding to the input 
        magnetic field values.
    """

    # Load the data from the paper
    Cs_a_vs_B = np.loadtxt('Cs_a_vs_B.txt')

    # Interpolate the data
    a = np.interp(B, Cs_a_vs_B[:, 0], Cs_a_vs_B[:, 1], left=np.nan, 
                  right=np.nan)
    return a


def aLiCs(B, species='a') -> float | np.ndarray:
    """Li-Cs scattering length as a function of magnetic field.

    This function calculates the Li-Cs scattering length at a given magnetic 
    field B, for Li atoms in the state given by species and Cs atoms in the
    lowest hyperfine state. The scattering lengths are calculated using the 
    resonance positions and widths from Jacob's Nature Physics paper.
    
    Parameters
    ----------
    B : float or array-like
        Magnetic field values in Gauss.
    species : str, default: 'a'
        The state of the Li atoms ('a', 'b', or 'c').

    Returns
    -------
    a : float or array-like, same shape as B
        Scattering lengths in units of Bohr radii corresponding to the input 
        magnetic field values.
    """

    if species == 'a':
        return -29.4 * (1 - (-58.21) / (B - 842.829) - (-4.55) / (B - 892.648))
    elif species == 'b':
        return -29.6 * (1 - (-0.37) / (B - 816.113) - (-57.45) / (B - 888.577) - (-4.22) / (B - 943.033))
    else:
        raise ValueError("Invalid species. Must be 'a' or 'b'.")


def aLiLi(B: float | np.ndarray, species: str) -> float | np.ndarray:
    """Li-Li scattering length as a function of magnetic field.
    
    Calculate the Li-Li scattering length at a given magnetic field B, for Li atoms in the state given by species.
    Values from Table 2 in "Precise Characterization of $^{6}\mathrm{Li}$ Feshbach Resonances Using 
    Trap-Sideband-Resolved RF Spectroscopy of Weakly Bound Molecules" by Zurn, et. al. PRL 110 13 135301 (2013).
    Parameters:
    B (float or array-like): Magnetic field values in Gauss.
    species (str): The state of the Li atoms ('ab', 'ac', or 'bc').
    Returns:
    numpy.ndarray or float: Scattering lengths in Bohr radii corresponding to the input magnetic field values.
    """

    # Set parameters depending on species combination
    if species == 'ab':
        abg = -1582
        delta = -262.3
        B0 = 832.18
    elif species == 'ac':
        abg = -1770
        delta = -166.6
        B0 = 689.68
    elif species == 'bc':
        abg = -1490
        delta = -222.3
        B0 = 809.76
    else:
        raise ValueError("Invalid species. Must be 'ab', 'ac', or 'bc'.")

    # Calculate scattering length
    a = abg * (1 - delta / (B - B0))

    return a


def aLiCs_molscat(B: float | np.ndarray, species: str) -> float | np.ndarray:
    """
    Calculate the Li-Cs scattering length from molscat calculations stored in 'LiCs_a_vs_B.txt'. This will likely be 
    less accurate than the simple model, but on the other hand provides support for Li-c species. 
    Parameters:
    B (float or array-like): Magnetic field values in Gauss. Only values between 800 and 980 are valid currently.
    species (str): The state of the Li atoms ('a', 'b', or 'c').
    Returns:
    numpy.ndarray or float: Scattering lengths in Bohr radii corresponding to the input magnetic field values.
    """

    # Load the data from the file
    M = np.loadtxt('LiCs_a_vs_B.txt', skiprows=1, delimiter=',')
    Bload = M[:, 0]

    # Select the correct species
    if species == 'a':
        aload = M[:, 1]
    elif species == 'b':
        aload = M[:, 2]
    elif species == 'c':
        aload = M[:, 3]
    else:
        raise ValueError("Invalid species. Must be 'a', 'b', or 'c'.")

    # Interpolate the data
    return np.interp(B, Bload, aload, left=np.nan, right=np.nan)


if __name__ == '__main__':

    # Check constants
    print(Cs133.D2.Gamma)

    # Test aCsCs function
    print('Testing aCsCs function')
    # Cs_a_vs_B = np.loadtxt('Cs_a_vs_B.txt')
    # print(Cs_a_vs_B)
    B = np.linspace(880, 895, 15)
    a = aCsCs(B)
    print(a)

    # Test aLiCs function
    print('Testing aLiCs function')
    B = np.array([890, 892, 893, 895])
    a = aLiCs(B, 'a')
    print(a)

    # Test aLiLi function
    print('Testing aLiLi function')
    B = np.linspace(700, 1000, 5)
    a = aLiLi(B, 'ab')
    print(a)

    # Test aLiCs_molscat function
    print('Testing aLiCs_molscat function')
    # LiCs_a_vs_B = np.loadtxt('LiCs_a_vs_B.txt', skiprows=1, delimiter=',')
    # print(LiCs_a_vs_B)
    B = np.array([890, 892, 893, 895])
    a = aLiCs_molscat(B, 'a')
    print(a)
