import base64
import dash
from dash import dcc, html
import dash_bootstrap_components as dbc
from rnapolis import annotator, parser
from io import BytesIO, StringIO
import tempfile

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
                    style = {'fontSize': '16px'}
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
            [
                html.H1("Welcome to RNA Graph application", className='h1', id='welcome-text',),
                html.P("Upload your RNA structure files (PDB or CIF format) to visualize 3D structures and analyze key interactions.", className='p1'),
                html.Div(
                [
                    html.Div(
                        [
                            html.Div(
                                [
                                    html.Img(
                                        src="/assets/download.png", 
                                        style={"width": "48px", "height": "48px"}  
                                    ),
                                    html.H2("Upload and Visualize RNA Structures", className="step-title"),
                                ],
                                style={"display": "flex", "alignItems": "center", "gap": "12px", "justifyContent": "center"}  
                            ),
                            html.P(
                                "Easily upload your RNA structure files in PDB or CIF format to visualize them in 3D. Our tool supports interactive exploration, allowing you to manipulate the structure, rotate it, zoom in, and more.",
                                className="p2",
                            ),
                        ],
                        className="step-box",
                    ),
                    html.Div(
                        [
                            html.Div(
                                [
                                    html.Img(
                                        src="/assets/image.png", 
                                        style={"width": "48px", "height": "48px"}  
                                    ),
                                    html.H2("Switch Between Views", className="step-title"),
                                ],
                                style={"display": "flex", "alignItems": "center", "gap": "12px"}  
                            ),
                            html.P(
                                "Choose from different views to analyze your RNA structures. You can switch between a detailed atomic view for atom-level precision and an interaction network view to focus on key molecular interactions.",
                                className="p2",
                            ),
                        ],
                        className="step-box",
                    ),
                    html.Div(
                        [
                            html.Div(
                                [
                                    html.Img(
                                        src="/assets/filters.png", 
                                        style={"width": "48px", "height": "48px"}  
                                    ),
                                    html.H2("Customize Interactions and Colors", className="step-title"),
                                ],
                                style={"display": "flex", "alignItems": "center", "gap": "12px", "justifyContent": "center"}  
                            ),
                            html.P(
                                "You can display various types of interactions such as phosphodiester bonds, canonical and non-canonical base-base interactions, and stacking interactions. Plus, adjust the colors of nucleotides and interaction lines to suit your preferences, making your visualizations both more informative and visually appealing.",
                                className="p2",
                            ),
                        ],
                        className="step-box",
                    ),
                ],
                className="steps-container",
                ),
            ],
            id='title-container',
            className='title-container',
        ),
        
    ],
    id = "app-container",
    className="app-container",
)

rna_nucleotides = ['A', 'C', 'G', 'U', 'I']
dna_nucleotides = ['DA', 'DC', 'DG', 'DU', 'DI', 'DT']
def check_nucleotide_type_and_completeness(decoded_data, ext):
    content = decoded_data.decode('utf-8')
    nucleotides_found = False
    is_complete = True
    issues = ""

    if ext == 'pdb':
        coordinates_present = False

        for line in content.splitlines():
            if line.startswith('ATOM') or line.startswith('HETATM'):
                if line.startswith('ATOM'):
                    residue_name = line[17:20].strip() 

                    if residue_name in rna_nucleotides + dna_nucleotides:
                        nucleotides_found = True
                       
                try:
                    x = float(line[30:38].strip())
                    y = float(line[38:46].strip())
                    z = float(line[46:54].strip())
                    coordinates_present = True
                except ValueError:
                    coordinates_present = False
        if not coordinates_present:
            issues = "*Missing atomic coordinates."
            return "Other", is_complete, issues
        
    elif ext == 'cif':
        coordinates_present = False

        for line in content.splitlines():
            if line.startswith('_atom_site'):
                continue
            elif line.startswith('_'):
                continue
            elif line.startswith('ATOM'):
                parts = line.split()
                if len(parts) > 6:
                    if parts[5] in rna_nucleotides + dna_nucleotides:
                        nucleotides_found = True

                    try:
                        x = float(parts[10])
                        y = float(parts[11])
                        z = float(parts[12])
                        coordinates_present = True
                    except (IndexError, ValueError):
                        coordinates_present = False
        if not nucleotides_found:
            issues = "*File does not contain any RNA or DNA structure."
            return "Other", is_complete, issues

        if not coordinates_present:
            issues = "*Missing atomic coordinates."
            return "Other", is_complete, issues


    return "RNA", is_complete, issues

def extract_structure_name(decoded_data, ext):
    content = decoded_data.decode('utf-8')
    pdb_id = "Unknown PDB ID"

    if ext == 'pdb':
        for line in content.splitlines():
            if line.startswith("HEADER"):
                pdb_id = line[62:66].strip() 
                break

    elif ext == 'cif':
        for line in content.splitlines():
            if line.startswith("_entry.id"):
                parts = line.split(maxsplit=1)
                if len(parts) > 1:
                    pdb_id = parts[1].strip().strip('"')  
                break

    return pdb_id


@app.callback(
    [
        dash.dependencies.Output('upload-message', 'children'),
        dash.dependencies.Output('store', 'data'),
        dash.dependencies.Output('processed-data', 'data'),
        dash.dependencies.Output('title-container', 'style'),
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
        return [None, None, None, {'display': 'flex'}, molviewer_class, RNAgraph_class]

    try:
        content_type, content_string = contents.split(',')
        decoded = base64.b64decode(content_string)
        file_ext = filename.split('.')[-1].lower()

        if file_ext in ['pdb', 'cif']:
            file_base64 = base64.b64encode(decoded).decode()
            file_data_url = f"data:{content_type};base64,{file_base64}"

            structure_type, is_complete, issues = check_nucleotide_type_and_completeness(decoded, file_ext)
            if structure_type == "Other" or not is_complete:
                return [html.Div(issues), None, None, {'display': 'flex'}, molviewer_class, RNAgraph_class]
            
            if pathname == '/':
                return [None, {'url': file_data_url, 'ext': file_ext, 'name': None}, None, {'display': 'none'}, molviewer_class, RNAgraph_class]
            
            structure_name = extract_structure_name(decoded, file_ext)
            interactions = calculate_interactions(decoded, file_ext)
            return [None, {'url': file_data_url, 'ext': file_ext, 'name': structure_name}, interactions, {'display': 'none'}, molviewer_class, RNAgraph_class]

        return [html.Div('*Invalid file format. Please upload a PDB or CIF file.'), None, None, {'display': 'flex'}, molviewer_class, RNAgraph_class]

    except Exception as e:
        return [html.Div('*There was an error processing this file.'), None, None, {'display': 'flex'}, molviewer_class, RNAgraph_class]

def calculate_interactions(decoded_data, ext):
    try:
        
        if ext == 'pdb':
            decoded_string = decoded_data.decode('utf-8')
            structure = parser.read_3d_structure(StringIO(decoded_string))
        elif ext == 'cif':
            decoded_string = decoded_data.decode('utf-8')
            with tempfile.NamedTemporaryFile(delete=False, suffix='.cif') as temp_file:
                temp_file.write(decoded_data)  
                temp_file.flush()  
                temp_file.seek(0)  
                try:
                    with open(temp_file.name, 'r') as read_file:
                        structure = parser.read_3d_structure(read_file)  
                except Exception as parse_error:
                    print(f"Error while parsing CIF file: {parse_error}")
                    return None

        available_interactions = {
            'phosphodiester': annotator.extract_base_interactions(structure).basePhosphateInteractions,
            'c_base_base': [pair for pair in annotator.extract_base_interactions(structure).basePairs if pair.lw.name == 'cWW'],
            'nc_base_base': [pair for pair in annotator.extract_base_interactions(structure).basePairs if pair.lw.name != 'cWW'],
            'stacking': annotator.extract_base_interactions(structure).stackings
        }

        return available_interactions

    except Exception as e:
        print(f"Error in calculate_interactions: {e}")
        return None

if __name__ == "__main__":
    #app.run_server(debug=True, dev_tools_ui=False)
    app.run_server(debug=True)
    #serve(app.server, host="0.0.0.0", port=8050)