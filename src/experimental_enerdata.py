"""

Module for processing and reading experimental binding energies database

@author: Ale Mendez
"""
import pandas as pd
import os 
import src.miscellaneous as misc

class experimentalData:


    def __init__(self, folder, units):

        assert units == 'eV' or 'Rydberg' or 'Hartree', 'units should be eV, Rydberg or Hartree'

        self.filename = 'ElectronBindingEnergies'
        self.file_ext = '.tsv'
        self.file_proc_key = ['dat','ref']
        self.folder = folder
        self.units = units
        self.raw_file_path = None
        self.raw_table = None
        self.idx = None
        self.orbs = None
        self.elements = None
        self.bindener = None
        self.ref_table = None
        self.dat_table = None
        self.load_database()


    def load_database(self):
        '''
        Read all the experimental binding energy data available
        '''
        proc_raw = self.check_database_files(key='raw')
        if proc_raw:
            self.load_raw_table()
            self.define_global_variables()

        proc_other = self.check_database_files()
        if proc_other:
            self.load_energy_table()
            self.load_reference_table()
            print('proc')
        else:
            self.proc_raw_table()
            self.write_processed_tables()
            print('raw')


    def check_database_files(self, key=None):
        # check if file with raw experimental data exists
        if key == 'raw': 
            raw_file = ''.join([self.filename, self.file_ext])
            self.raw_file_path = os.path.join(self.folder, raw_file)
            if not os.path.isfile(self.raw_file_path):
                raise OSError(f'{self.raw_file_path} does not exist.')

        # check if processed files exist
        else:
            for key in self.file_proc_key:
                proc_file_path = self.processed_filepath(key)
                if not os.path.isfile(proc_file_path): 
                    return False
        return True


    def load_raw_table(self):
        self.raw_table = pd.read_csv(self.raw_file_path, sep='\t', index_col=0)


    def processed_filepath(self, key):
        proc_file = '_'.join([self.filename, key])
        proc_file = ''.join([proc_file, self.file_ext])
        proc_filepath = os.path.join(self.folder, proc_file)
        return proc_filepath


    def load_energy_table(self):
        datpath = self.processed_filepath(self.file_proc_key[0])
        self.dat_table = pd.read_csv(datpath, sep='\t', index_col=0)


    def load_reference_table(self):
        refpath = self.processed_filepath(self.file_proc_key[1])
        self.ref_table = pd.read_csv(refpath, sep='\t', index_col=0)


    def define_global_variables(self):
        self.idx = self.raw_table.index
        self.orbs = self.raw_table.columns[1:]
        self.elements = self.raw_table['Element'].tolist()


    def proc_raw_table(self):
        self.proc_references()
        self.proc_data()


    def proc_references(self):
        '''
        Function to convert references flags in raw data to numeric references
        '''
        self.ref_table = self.raw_table.copy()
        for i in self.idx: 
            for o in self.orbs:
                val = self.raw_table.loc[i][o]
                ref = val
                if not misc.isNaN(val): 
                     ref = [1]
                     if '*' in val: ref = [2]
                     if '+' in val: ref = [3]
                     if 'a' in val: ref.append('a')
                     if 'b' in val: ref.append('b')
                self.ref_table.at[i,o] = ref


    def proc_data(self):
        '''
        Function to process data table: convert string values to float values
        '''
        self.dat_table = self.raw_table.copy()
        for i in self.idx: 
            for o in self.orbs:
                val = self.raw_table.loc[i][o]
                if not misc.isNaN(val): 
                    for char in '*+ab': val = val.replace(char,'')
                    if 'g' in val: val = val.replace('g','9') # fix on Williams compilation pdf
                    val = float(val)
                self.dat_table.at[i,o] = val


    def element_binding_energies(self, element_str, bprint=False):
        '''
        Selects element binding energy data from table and prints it in output file

            element: (str)  element symbol, e.g. 'He'
            units: (str)  units for converting binding energies
            print: (bool) print output file with element data

        '''
        # define an element object (inherits methods from periodictable.elements module)
        element = misc.periodic_table(element_str) 
        self.check_element_data(element.symbol)

        self.bindener = self.extract_element_data(element.number)
        if bprint: self.print_element_data(element.symbol)
        return self.bindener


    def check_element_data(self, element_symbol):
        ''' 
        Check if element_symbol (str) has any data associated 
        '''
        if element_symbol not in self.elements:
            raise ValueError(f'No data found for {element_symbol}.')


    def extract_element_data(self, element_number):
        '''
        Extracts binding energy data from table according to element selected
        '''
        defcolname = misc.column_name('eV')
        colname = misc.column_name(self.units)
        ener = self.dat_table.loc[element_number][1:].tolist()
        bindener = pd.DataFrame(index=self.orbs)
        bindener.index.name = 'Orbital'
        bindener[defcolname] = ener
        bindener[colname] = [misc.convert_energy_from_eV(e, self.units) for e in ener]
        bindener['Reference'] = self.ref_table.loc[element_number][1:].tolist()
        bindener = bindener.dropna()
        return bindener


    def print_element_data(self, element_symbol):
        '''
        Function to print binding energy data for selected element in units defined
        '''
        fout = os.path.join(self.folder, element_symbol+'_experiment.dat')
        colname = misc.column_name(self.units)
        print_bindener = self.bindener.dropna()
        dictdata = dict(zip(print_bindener.index, print_bindener[colname]))
        with open(fout, 'w') as f:
            print("{}\t{}\t{}".format('Orb', colname, 'Reference'), file=f)
            for orb, ener_val in dictdata.items():
                if not misc.isNaN(ener_val) and 'eV' not in self.units: 
                    ref = self.bindener.loc[orb]['Reference']
                    ener_eV = self.bindener.loc[orb][misc.column_name('eV')]
                    ener_val = self.significant_figures(ener_eV, ener_val)
                print("{}\t{}\t{}".format(orb, ener_val, ref),file=f)
            self.print_footnote(f)


    def significant_figures(self, ener, conv):
        '''
        Function to truncate the converted energy values with the same number of
        significate figures as experimental data
        '''
        from math import log10
        if ener == 0:
            return str(ener)
        # count the number of significate figures (icsig) in experimental value (eV)
        for i in range(6):
            a = ener*(10**i)
            Da = a-round(a)
            if Da == 0: 
                iord = i
                icsig = int(log10(abs(ener*10**iord)))+1
                break
        # truncate converted values using icsing 
        aux = 1
        if conv<1: aux += 1 # add one more integer for "0."
        if conv<0.1: aux += 1 # add one more integer for "0.0"
        if conv<0.01: aux += 1 # add one more integer for "0.00"
        if conv<0.001: aux += 1 # add one more integer for "0.000"
        conv_string = str(conv)[:icsig+aux] # valores truncados (no redondeados) <=== CAMBIAR a redondear
        return conv_string


    def print_footnote(self, f):
        '''
        Function to print footnote with references
        '''
        print('\n# References:',file=f)
        print('# Experimental values compiled by Williams, G. (1995)',file=f)
        print('# [1] J. A. Bearden and A. F. Burr, "Reevaluation of X­Ray Atomic Energy Levels," Rev. Mod. Phys. 39, (1967) p.125\n'+
              '# [2] M. Cardona and L. Ley, Eds., Photoemission in Solids I: General Principles (Springer­Verlag, Berlin, 1978) with additional corrections\n'+
              '# [3] J. C. Fuggle and N. Mårtensson, "Core­Level Binding Energies in Metals," J. Electron Spectrosc. Relat. Phenom. 21, (1980) p.275',file=f)


    def write_processed_tables(self):
        # print reference table to file
        refpath = self.processed_filepath(self.file_proc_key[1])
        self.ref_table.to_csv(refpath, sep='\t')
        # print energy table to file
        datpath = self.processed_filepath(self.file_proc_key[0])
        self.dat_table.to_csv(datpath, sep='\t')