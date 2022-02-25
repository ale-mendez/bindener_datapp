import src.miscellaneous as misc
import pandas as pd
import os


def read_diracfock_bindener(folder):
    filename = 'ElectronBindingEnergies.tsv'
    fpath = os.path.join(folder, filename)
    bindener = pd.read_csv(fpath, sep='\t', header='infer', index_col=0)

    bindener_dict = dict()
    for atom in bindener.index:
        df = bindener.loc[atom:atom]
        df = df.transpose()
        df = df.rename(columns={atom: 'Energy(Hartree)'})
        df = df.dropna()
        bindener_dict[atom] = df

    return bindener_dict


class theoreticalData:

    def __init__(self, folder, units):

        assert units == 'eV' or 'Rydberg' or 'Hartree', 'units should be eV, Rydberg or Hartree'
        self.folder = folder
        self.ener_filename = 'bindener.dat'
        self.input_units = None
        self.units = units
        self.bindener_data = self.load_database()
        self.bindener = None
        

    def load_database(self):
        '''
        Read all the theoretical binding energy data available in folder
        '''
        if not misc.check_folder_exists(self.folder):
            raise IOError(f'{self.folder} does not exists.')

        # dirac-fock data is given with a different format (tsv)
        if 'dirac-fock' in self.folder:
            bindener_dict = read_diracfock_bindener(self.folder)
        else:
            # list atoms with theoretical data
            for root, dirs, files in os.walk(self.folder):
                if root == self.folder:
                    atoms = [atom for atom in dirs if misc.periodic_table(atom)]
                    break

            # load binding energy theoretical data
            bindener_dict = dict()
            for atom in atoms:
                fpath = os.path.join(self.folder, atom, self.ener_filename)
                try:
                    df = pd.read_csv(fpath, sep='\t', comment='#', index_col=0)
                    # check units and convert energy
                    self.input_units = misc.determine_energy_units(df.columns)
                    if self.input_units != self.units:
                        df = misc.convert_energy_units(df, self.input_units, self.units)
                except:
                    df = None
                bindener_dict[atom] = df

        return bindener_dict


    def element_binding_energies(self, element):
        """
        Output perturbative data for element given
        """
        return self.bindener_data[element]
