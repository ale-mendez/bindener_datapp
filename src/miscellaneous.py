"""

Module with miscellaneous functions used in experimental and perturbative modules

@author: Ale Mendez
"""
import os
import numpy as np
import re
import pandas as pd


def print_title(title):
    print('#'*(len(title)+4))
    print('#',title,'#')
    print('#'*(len(title)+4))
    

def isNaN(value):
    '''
    Check if value (int, float, list, array or serie) is/contains NaN value. If so, returns True 
    '''
    if isinstance(value, (int, float)):
        return value != value
    elif isinstance(value, list):
        return np.isnan(value).any()
    elif isinstance(value, pd.Series):
        return value.isnull().values.any()

def column_name(units):
    '''
    Write column name for energy using units given
    '''
    return 'Energy('+units+')'

def determine_energy_units(columns):
    '''
    Function to determine the energy units of a Dataframe with "columns"
    '''
    unit = [name for col in columns if "energy" in col or "Energy" in col for name in re.split(r'[()]', col) if name in ['eV', 'Rydberg', 'Hartree']]
    if len(unit) > 0:
        return unit[0]

def convert_energy_units(df, input_units, units):
    '''
    Convert binding energy values from default to given by user
    '''
    # convert data from eV to units
    def_colname = column_name(input_units)
    colname = column_name(units)
    ener = df[def_colname]
    if input_units == 'eV':
        df[colname] = convert_energy_from_eV(ener, units)
    # convert data from Rydbergs to units
    elif input_units == 'Rydberg':
        df[colname] = convert_energy_from_Rydberg(ener, units)
    # convert data from Hartree to units
    elif input_units == 'Hartree':
        df[colname] = convert_energy_from_Hartree(ener, units)
    return df

# energy conversion 

def physical_constants(constant_name):
    from scipy.constants import physical_constants
    return physical_constants[constant_name]

def convert_energy_from_eV(ener, units):
    ''' Convert energy units from eV to Hartree or Rydberg '''
    conv_hartree_to_eV = physical_constants('Hartree energy in eV')[0]
    conv = None
    conv_ener = None
    if units=='eV': conv = 1.0
    if units=='Hartree': conv = conv_hartree_to_eV
    if units=='Rydberg': conv = conv_hartree_to_eV/2.0
    if not isNaN(ener): conv_ener = ener/conv # esto solo funciona para float, cuando uso una serie, tengo que usar .any()
    return conv_ener


def convert_energy_from_Rydberg(ener, units):
    '''
    Convert energy units from Rydberg to Hartree or eV
    '''
    conv_hartree_to_eV = physical_constants('Hartree energy in eV')[0]
    conv = None
    conv_ener = None
    if units=='eV': conv = 0.50*conv_hartree_to_eV
    if units=='Hartree': conv = 0.50
    if units=='Rydberg': conv = 1.0
    if not isNaN(ener): conv_ener = ener*conv
    return conv_ener


def convert_energy_from_Hartree(ener, units):
    '''
    Convert energy units from Hartree to Rydberg or eV
    '''
    conv_hartree_to_eV = physical_constants('Hartree energy in eV')[0]
    conv = None
    conv_ener = None
    if units=='eV': conv = conv_hartree_to_eV
    if units=='Hartree': conv = 1.0
    if units=='Rydberg': conv = 2.0
    if not isNaN(ener): conv_ener = ener*conv
    return conv_ener


def shorten_units(units):
    if units == 'Rydberg': 
        return 'Ryd'
    elif units == 'Hartree': 
        return 'Eh'
    if units == 'eV': 
        return units


# compute error stuff

def ratio(experiment, calculation):
    return experiment/calculation

def relative_error(experiment, calculation):
    return ratio(experiment, calculation)-1

def abs_relative_error(experiment, calculation):
    return abs(ratio(experiment, calculation)-1)

def percent_relative_error(experiment, calculation):
    return abs_relative_error(experiment, calculation)*100

def mean_value(df_column):
    x = df_column.dropna()
    return sum(abs(x))/len(x)

# check stuff 

def check_energy_units(units):
    if units not in ['eV','Rydberg','Hartree']: 
        raise ValueError('variable units must be eV, Rydberg or Hartree.')
    return units

def check_file_exists(fpath):
    return os.path.isfile(fpath)

def check_folder_exists(folder):
    return os.path.isdir(folder) 

# other functions

def periodic_table(element_str):
    ''' 
    Inherits periodic table properties given by periodictable module 

    Example: 
    >> periodic_table('He').name
    'helium'
    >> periodic_table('He').number
    2
    
    '''
    import periodictable
    for el in periodictable.elements:
        if el.symbol == element_str: 
            return el

def FEG_params(ne, at_density, at_weight):
    ''' 
    Compute FEG parameters 
    
    Input:
            ne         -- number of electrons in FEG
            at_density -- atomic density (g/cm^3)
            at_weight  -- atomic number (g)
    
    Output:
            rs -- Wigner-Seitz radii (a.u.)
            EF -- Fermi energy (Hartree)

    '''
    no_au = (8.916E-2*ne*at_density/at_weight)
    rs = (3.0/(4.0*np.pi*no_au))**(1/3)
    EF = 1.837/rs**2
    return rs, EF

# plotting functions

def plot_fermi_energy(axs, FEG_electrons, element):
    FEGparams = FEG_params(FEG_electrons,element.density,element.mass)
    axs.axhline(FEGparams[1], 0, 1, ls='dashed', color='tab:gray')


def plot_gridlines(axs):
    axs.set_axisbelow(True)
    axs.yaxis.grid(color='gray', linestyle='dotted')
    axs.xaxis.grid(color='gray', linestyle='dotted')

