
import dash
import dash_core_components as dcc
import dash_html_components as html
import dash_bootstrap_components as dbc
from dash.dependencies import Input, Output, State
import periodictable
import src.structure_data as struc
import src.miscellaneous as misc
from src.app_styling import *

atoms = {el.symbol: el.number for i, el in enumerate(periodictable.elements) if i > 0}
atoms_options = [
    {
        'label': symbol + f'({number})',
        'value': symbol
    }
    for symbol, number in atoms.items()
]

app = dash.Dash(external_stylesheets=[dbc.themes.BOOTSTRAP])

atoms_dropdown = dbc.Card(
    [
        dbc.Label('Atoms', html_for="dropdown"),
        dcc.Dropdown(
            options=atoms_options,
            value=None,
            placeholder="Select an atom",
            id='input-atoms'
        ),
    ],
    className="mb-3",
)

methods_dropdown = dbc.Card(
    [
        dbc.Label('Methods', html_for="dropdown"),
        dcc.Dropdown(
            multi=True,
            id='dropdown-methods'
        ),
    ],
    className="mb-3",
)

units = 'Hartree'
units_items = dbc.Card(
    [
        dbc.Label('Energy units'),
        dbc.RadioItems(
            id='input-units',
            options=[
                {
                    'label': 'Hartree',
                    'value': 'Hartree'
                },
                {
                    'label': 'Rydberg',
                    'value': 'Rydberg'
                },
                {
                    'label': 'eV',
                    'value': 'eV'
                }
            ],
            value=units)
    ],
    className="mb-3",
    style={
        'margin': 'auto'
    }
)

feg_input = dbc.Card(
    [
        html.Label('Electrons in FEG: '),
        dcc.Input(
            id='feg-params',
            value=0,
            type='number',
            min=0,
            # max=bindener.atom.number
        ),
        html.Br(),
        html.Table([
            html.Tr([html.Td(['Wigner-Seitz radii:']), html.Td(id='output-rs')]),
            html.Tr([html.Td(['Fermi energy:']), html.Td(id='output-EF')])
        ]),
    ]
)

param_sidebar = html.Div(
    [
        html.H2('Parameters', style=TEXT_STYLE),
        html.Hr(),
        dbc.Form([atoms_dropdown, methods_dropdown, units_items, feg_input])
    ],
    style=SIDEBAR_STYLE
)

config = {
    'displayModeBar': True,
    'toImageButtonOptions': {
        # 'filename': f'bindener_{bindener.atom_symbol}',
        'filename': f'bindener',
        'format': 'svg',
        'scale': 1
    },
    'edits': {
        'legendPosition': True,
        'legendText': True,
        'titleText': True,
        'axisTitleText': True}
}

bindener_graph = dcc.Graph(
    id='output-bindener',
    config=config
    )

download_fig = html.Div(
    [
        dbc.Button("Download Image", id="btn_image", className="me-md-2"),
        # html.Button("Download Image", id="btn_image"),
        dcc.Download(id="download-image")
    ],
)

download_csv = html.Div(
    [
        dbc.Button("Download Data", id="btn_data", className="me-md-2"),
        # html.Button("Download Data", id="btn_data"),
        dcc.Download(id="download-data")
    ],
)

download_buttons = html.Div(
    [
        download_fig, 
        download_csv
    ],
    className="d-grid gap-2 d-md-flex justify-content-md-end",
)

content = html.Div(
    [
        html.H2('Binding Energy Dashboard', style=TEXT_STYLE),
        html.Hr(),
        bindener_graph,
        download_buttons,
        # orbitals_slider
    ],
    style=CONTENT_STYLE
)

app.layout = html.Div([
    param_sidebar,
    content          
])

@app.callback(
    Output(component_id="dropdown-methods", component_property="options"),
    Output(component_id="dropdown-methods", component_property="value"),
    Output(component_id="feg-params", component_property="max"),
    Output(component_id="output-rs", component_property="children"),
    Output(component_id="output-EF", component_property="children"),
    Output(component_id="output-bindener", component_property="figure"),
    # Output(component_id="output-bindener", component_property="config"),
    Input(component_id="input-atoms", component_property="value"),
    Input(component_id='input-units', component_property='value'),
    Input(component_id="dropdown-methods", component_property="value"),
    Input(component_id='feg-params', component_property='value'),
)
def update_methods(input_atoms, input_units, input_methods, input_nFEG=0):

    methods = []
    out_methods = []
    fig = {}
    max_nFEG = None
    rs_string = ''
    Ef_string = ''
    filename = 'bindener'
    download = None

    if input_atoms:

        # create object with bindener data
        bindener = struc.bindingEnergies(atom_symbol=input_atoms, units=input_units)
        filename = f'bindener_{bindener.atom_symbol}'
        
        # make list of methods with data
        methods = list(bindener.methods)
        orbitals = list(bindener.orbitals)

        out_methods = methods
        ibool = set(methods) == set(input_methods)
        if not ibool and input_methods:
            out_methods = input_methods

        # compute FEG parameters
        max_nFEG = bindener.atom.number
        if input_nFEG is not None and input_nFEG > 0:
            units_short = misc.shorten_units(input_units)
            rs, Ef = bindener.compute_FEG_parameters(input_nFEG)
            rs_string += f'{rs:.2f} a.u.'
            Ef_string += f'{Ef:.2f} {units_short}'

        # create figure
        fig = bindener.binding_energies_graph(orbitals=orbitals, methods=out_methods)

    # return methods, methods, max_nFEG, rs_string, Ef_string, fig, filename
    # return methods, methods, max_nFEG, rs_string, Ef_string, fig, download
    return methods, out_methods, max_nFEG, rs_string, Ef_string, fig

# @app.callback(
#     Output("download-image", "data"),
#     Input("btn_image", "n_clicks"),
#     prevent_initial_call=True,
# )
# def func(n_clicks):

#     if n_clicks:
#         download = dcc.send_file('...')
#     return download


if __name__ == '__main__':
    app.run_server(debug=True, use_reloader=True)
