import base64
import dash
from dash import dcc, html
import dash_bootstrap_components as dbc
from rnapolis import annotator, parser
from io import BytesIO, StringIO
from waitress import serve

app = dash.Dash(__name__, use_pages=True, external_stylesheets=[dbc.themes.BOOTSTRAP], suppress_callback_exceptions=True,  meta_tags=[
        {"name": "viewport", "content": "width=device-width, initial-scale=1"}
    ])

app.layout = html.Div(
    [
        html.Div(
            [
                dcc.Location(id='url', refresh=False),
                html.H1("RNA Graph", style={'textAlign': 'left'}),
                html.Div([
                    dcc.Upload(
                    id='upload-data',
                    className = 'upload-data',
                    children=html.Div(['Drag and Drop or ', html.A('Select a PDB or CIF File')]),
                    multiple=False,
                    style = {'font-size': '16px'}
                    ),
                    html.Div(id='upload-message', className='upload-message'),
                ]),
                dcc.Store(id="store"),
                dcc.Store(id='processed-data'),
                dbc.NavbarSimple(
                    children=[
                        dbc.NavItem(
                            dcc.Link(page['name'], href=page['path'], id=f"{page['name'].lower()}-link", className="nav-link"),  
                        )
                        for page in dash.page_registry.values()
                    ],
                    color = '#021425',
                )
            ],
            className="header", 
        ),
        html.Div(dash.page_container, className="content"),
        html.Div(
            html.H1("Welcome to RNA Graph application", className='h1', id='welcome-text',),
        ),
    ],
    className="app-container",
)

rna_nucleotides = ['A', 'C', 'G', 'U', 'I']
dna_nucleotides = ['DA', 'DC', 'DG', 'DU', 'DI', 'DT']
def check_nucleotide_type(decoded_data):
    pdb_content = decoded_data.decode('utf-8')
    
    for line in pdb_content.splitlines():
        if line.startswith('ATOM') or line.startswith('HETATM'):
            residue_name = line[17:20].strip()  
            if residue_name in rna_nucleotides:
                return "RNA"
            elif residue_name in dna_nucleotides:
                return "DNA"
    return "Other"
@app.callback(
    [
        dash.dependencies.Output('upload-message', 'children'),
        dash.dependencies.Output('store', 'data'),
        dash.dependencies.Output('processed-data', 'data'),
        dash.dependencies.Output('welcome-text', 'style'),
    ] + [
        dash.dependencies.Output(f"{page['name'].lower()}-link", 'className') for page in dash.page_registry.values()
    ],
    [dash.dependencies.Input('upload-data', 'contents'), dash.dependencies.Input('url', 'pathname')],
    dash.dependencies.State('upload-data', 'filename'),
)
def update_active_link(contents, pathname, filename):
    molviewer_class = 'nav-link'
    RNAgraph_class = 'nav-link'

    if pathname == '/':
        molviewer_class = 'nav-link active'
    elif pathname == '/page-2':
        RNAgraph_class = 'nav-link active'

    if contents is None:
        return [None, None, None, {'display': 'block'}, molviewer_class, RNAgraph_class]

    try:
        content_type, content_string = contents.split(',')
        decoded = base64.b64decode(content_string)
        file_ext = filename.split('.')[-1].lower()

        if file_ext in ['pdb', 'cif']:
            file_base64 = base64.b64encode(decoded).decode()
            file_data_url = f"data:{content_type};base64,{file_base64}"

            structure_type = check_nucleotide_type(decoded)
            if structure_type == "Other":
                return [html.Div('*Invalid file format. Please upload a PDB or CIF file.'), None, None, {'display': 'block'}, molviewer_class, RNAgraph_class]
            
            interactions = calculate_interactions(decoded)
            return [None, {'url': file_data_url, 'ext': file_ext}, interactions, {'display': 'none'}, molviewer_class, RNAgraph_class]

        return [html.Div('*Invalid file format. Please upload a PDB or CIF file.'), None, None, {'display': 'none'}, molviewer_class, RNAgraph_class]

    except Exception as e:
        return [html.Div('*There was an error processing this file.'), None, {'display': 'none'}, molviewer_class, RNAgraph_class]

def calculate_interactions(decoded_data):
    structure = parser.read_3d_structure(StringIO(decoded_data.decode('utf-8')))
    
    available_interactions = {
        'phosphodiester': annotator.extract_base_interactions(structure).basePhosphateInteractions,
        'c_base_base': [pair for pair in annotator.extract_base_interactions(structure).basePairs if pair.lw.name == 'cWW'],
        'nc_base_base': [pair for pair in annotator.extract_base_interactions(structure).basePairs if pair.lw.name != 'cWW'],
        'stacking': annotator.extract_base_interactions(structure).stackings
    }

    return available_interactions

if __name__ == "__main__":
    #app.run_server(debug=True)
    serve(app.server, host="0.0.0.0", port=8050)