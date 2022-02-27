import pandas as pd
import src.miscellaneous as misc
import plotly.express as px
import plotly.graph_objects as go
import numpy as np
import src.miscellaneous as misc
import os
import src.experimental_enerdata as expapp
import src.theoretical_enerdata as theoapp


def pull_bindener_data(atom, units, data_folder):

    try:
        pathdir = os.path.join(self.main_folder, data_folder)
        if data_folder == 'experimental':
            atom_df = expapp.experimentalData(pathdir, units).element_binding_energies(atom)
        else:
            atom_df = theoapp.theoreticalData(pathdir, units).element_binding_energies(atom)
    except:
        atom_df = None

    return atom_df


class bindingEnergies:

    def __init__(self, main_folder=None, atom_symbol=None, units=None):
        self.atom_symbol = atom_symbol
        self.atom = misc.periodic_table(self.atom_symbol)
        self.units = units
        self.main_folder = main_folder
        self.experiment = pull_bindener_data(atom_symbol, units, 'experimental')
        self.relativistic = pull_bindener_data(atom_symbol, units, 'perturbative')
        self.diracfock = pull_bindener_data(atom_symbol, units, 'dirac-fock')
        # self.semirelat = pull_bindener_data(atom_symbol, units, 'semi-relativistic')
        self.hartreefock = pull_bindener_data(atom_symbol, units, 'hartree-fock')
        self.bindener = self.arrange_data_to_dataframe(units)
        self.orbitals = self.bindener.index
        self.methods = self.bindener.columns
        self.fermi_energy = None
        self.bindener_error = self.compute_relative_errors()


    def arrange_data_to_dict(self):
        input = {
            'Experimental': self.experiment,
            'Relativistic': self.relativistic,
            'Dirac-Fock': self.diracfock,
            # 'Semi-relativistic': self.semirelat,
            'Hartree-Fock': self.hartreefock}
        return {key: value for key, value in input.items() if value is not None}


    def get_orbitals(self, bindener_dict):
        orbs = None
        methods = list(bindener_dict.keys())
        if len(methods) == 1 and 'Experimental' in methods:
            orbs = bindener_dict['Experimental'].index
        if 'Relativistic' in methods:
            orbs = bindener_dict['Relativistic'].index
        elif 'Dirac-Fock' in methods:
            orbs = bindener_dict['Dirac-Fock'].index
        return orbs


    def arrange_data_to_dataframe(self, units):
        bindener_dict = self.arrange_data_to_dict()
        orbs = self.get_orbitals(bindener_dict)
        print
        bindener = pd.DataFrame(index=orbs)
        for method in bindener_dict.keys():
            data = bindener_dict[method]
            col_ener = [col for col in data.columns if units in col][0]
            if method in ['Experimental', 'Relativistic']:
                data = pd.DataFrame({method: data[col_ener]}, index=data.index)
            else:
                data = self.arrange_nonrelat_energies(data, col_ener, orbs, method)
            bindener = bindener.join(data)
        bindener = bindener.dropna(axis=0, how='all')
        return bindener


    def arrange_nonrelat_energies(self, df, col, orbs, method):
        nonrelat_orbs = df.index
        ener_arr = dict()
        for nlm in orbs:
            match = [nl for nl in nonrelat_orbs if (nlm[:-1] == nl) or (nlm == nl)]
            if match: 
                ener_arr[nlm] = df[col][match[0]]
        df_arr = pd.DataFrame.from_dict(ener_arr, orient='index', columns=[method])
        return df_arr


    def compute_relative_errors(self):
        methods = self.methods
        relat_err = pd.DataFrame(index=self.orbitals)
        if 'Experimental' in methods:
            exp = self.bindener['Experimental']
            methods = methods.drop('Experimental')
            for method in methods:
                relat_err[method] = (exp - self.bindener[method]) / exp
        return relat_err


    def compute_FEG_parameters(self, ne):
        ne = float(ne)
        C = 0.5 * (9 * np.pi / 4) ** (2/3)
        density = self.atom.density # units: g/cm^3
        mass = self.atom.mass # units (g)
        # atomic density in atomic units
        atomic_den = (8.916E-2 * ne * density / mass)
        # wigner-seitz radii
        rs = (3 / (4 * np.pi * atomic_den)) ** (1/3)
        # fermi energy
        Ef_hartree = C / rs ** 2
        # convert units
        Ef_units = misc.convert_energy_from_Hartree(Ef_hartree, self.units)
        self.fermi_energy = Ef_units
        return rs, Ef_units


    def log_minor_and_major_ticks(self, df):

        vmin = min([min(df[col]) for col in df.columns])
        imin = np.log10(vmin).round()
        imin = imin if 1 * 10 ** imin < vmin else imin - 1

        vmax = max([max(df[col]) for col in df.columns])
        imax = np.log10(vmax).round()
        imax = imax if 1 * 10 ** imax > vmax else imax + 1

        nexp = np.arange(imin, imax + 1, 1)
        og_ticksvals = [1 * 10 ** i for i in nexp]
        nw_ticksvals = [el for i in nexp for el in np.arange(1 * 10 ** i, 1 * 10 ** (i + 1), 1 * 10 ** i)]
        tickstext = [f"{val:.0e}" if val in og_ticksvals else "" for val in nw_ticksvals]

        return nw_ticksvals, tickstext


    def binding_energies_graph(self, orbitals=None, methods=None):

        if orbitals is None: orbitals = self.orbitals
        if methods is None: methods = self.methods

        markers = ['circle-open', 'square-open', 'diamond-open', 'triangle-up-open', 'cross-open']
        dmarkers = dict(zip(methods, markers))

        fig = px.scatter(
            self.bindener, 
            x = orbitals, 
            y = methods, 
            log_y = True, 
            template = 'simple_white')

        if self.fermi_energy:
            fig.add_trace(
                go.Scatter(
                    x = [orbitals[0], orbitals[-1]],
                    y = [self.fermi_energy, self.fermi_energy],
                    mode = "lines",
                    line = dict(color='grey', width=1.5, dash='dash'),
                    showlegend = False)
            )

        for i, marker in enumerate(dmarkers.values()):
            fig['data'][i]['marker']['symbol'] = marker

        fig.update_layout(
            title = f"Binding energies for {self.atom_symbol}",
            xaxis_title = "Orbitals",
            yaxis_title = f"Binding energies ({misc.shorten_units(self.units)})",
            legend_title = f"",
            hovermode = "x unified",
            plot_bgcolor = 'white',
            legend = dict(
                yanchor = "top",
                y = 0.99,
                xanchor = "left",
                x = 0.8
                )
            )

        fig.update_xaxes(
            # tickvals = [0] + list(orbitals) + [1],
            # ticktext = [''] + list(orbitals) + ['']
            ticks="inside",
            showgrid=True
            )

        nw_ticksvals, tickstext = self.log_minor_and_major_ticks(self.bindener)
        fig.update_yaxes(
            type = 'log',
            tickvals = nw_ticksvals,
            ticktext = tickstext,
            automargin = False,
            ticks = "inside",
            showgrid=True)

        fig.update_traces(
            marker = dict(line=dict(width=1.5)),
            hovertemplate = '%{y:.3f} ' + misc.shorten_units(self.units))

        return fig


    def plot_binding_energies(self):
        fig = self.binding_energies_graph()
        fig.show()


    def relative_error_graph(self, error=False):
        pass

