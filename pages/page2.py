import dash
from dash import dcc, html, callback, Output, Input, State
import base64
import numpy as np
import plotly.graph_objects as go
from Bio.PDB import PDBParser, MMCIFParser
from io import BytesIO, StringIO

rna_nucleotides = ['A', 'C', 'G', 'U', 'I']
dna_nucleotides = ['DA', 'DC', 'DG', 'DU', 'DI', 'DT']

colors = []
color_map = {
        'A': 'green',
        'C': 'blue',
        'G': 'red',
        'U': 'yellow',
        'I': 'gray',
        'DA': 'green',
        'DC': 'blue',
        'DG': 'red',
        'DU': 'yellow',
        'DI': 'gray',
        'DT': 'purple'
    }

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

dash.register_page(__name__, path='/page-2', name="RNA Graph")

layout = html.Div( 
    [
        # sel-nucleotide-info container
        html.Div(
            id='sel-nucleotide-info',
            className='sel-nucleotide-info',
            children=[
                html.P("Nucleotide Information", className='title'),
                html.Div(
                    className='label-info',
                    children=[
                        html.Label('Nucleotide'),
                        html.Label(id='nucleotide-input')
                    ]
                ),
                html.Div(
                    className='label-info',
                    children=[
                        html.Label('Chain'),
                        html.Label(id='chain-input')
                    ]
                ),
                html.Div(
                    children=[
                        html.Div(
                            className='label-info',
                            children=[
                                html.Label("Coordinates"),
                                html.Button(
                                    id='toggle-coordinates-btn',
                                    children='▼',  
                                    style={
                                        'background': 'none', 
                                        'border': 'none', 
                                        'font-size': '18px', 
                                        'cursor': 'pointer'
                                    }
                                ),
                            ],
                            style={'display': 'flex', 'align-items': 'start'}
                        ),
                        html.Div(
                            id='coordinates-container',
                            children=[
                                html.P("x: ", id='coord-x', style={'padding-left': '16px'}),
                                html.P("y: ", id='coord-y', style={'padding-left': '16px'}),
                                html.P("z: ", id='coord-z', style={'padding-left': '16px'}),
                            ],
                            style={'display': 'none'}
                        )
                    ]
                ),
                html.Button('Clear', id='clear-button', className='clear-button')  
            ],
            style={'display': 'none', 'flex-grow': '1'} 
        ),
        # Main container with the graph
        html.Div(
            [
                html.H1(
                    "Welcome to RNA Graph application",
                    id='title',
                    className='h1',
                    style={'text-align': 'center', 'margin' : 'auto', 'position' : 'absolute'}  
                ),
                dcc.Graph(
                    id='rna-graph',
                    style={'display': 'none'},  
                    className='graph'
                ),
                html.Div(id='upload-message', style={'color': 'red'}),
            ],
            className='rna-graph-container',
            style={'flex-grow': '2'}  
        ),  

        # Sidebar
        html.Div(
            id='side-bar',
            className='side-bar',
            style={'flex-grow': '1', 'display': 'none', 'border-width': '1px', 'border-style': 'solid', 'border-color': '#19191a'}, 
            children=[
                html.P(
                    "Nucleotide Interactions",
                    id='title2',
                    className='title'
                ),
                dcc.Checklist(
                    id='interaction-type',
                    className='checklist-label',
                    options=[
                        {'label': 'Wiazania fosfodiestrowe', 'value': 'phosphodiester'},
                        {'label': 'Canonical interactions', 'value': 'c_base_base'},
                        {'label': 'Non-canonical interactions', 'value': 'nc_base_base'},
                        {'label': 'Stacking interactions', 'value': 'stacking'}
                    ],
                    value=['phosphodiester'],  
                    inputClassName='checklist-input'
                )
            ]
        )   
    ],
    style={'display': 'flex', 'align-items': 'start', 'width': '100%', 'overflow': 'hidden', 'background-color': '#fafafb', 'height': '610px'},
)


@callback(
    Output('rna-graph', 'figure', allow_duplicate=True),
    Output('rna-graph', 'style'),
    Output('title', 'style'),
    Output('side-bar', 'style'),
    Output('sel-nucleotide-info', 'style', allow_duplicate=True),
    Input('store', 'data'),
    State('upload-data', 'filename'),
    prevent_initial_call = True
)
def update_rna_graph(data, filename):
    if data is None or filename is None:
        return go.Figure(), {'display' : 'none'}, {'display' : 'block'}, {'display' : 'none'}, {'display' : 'none'}

    file_ext = filename.split('.')[-1].lower()
    if file_ext not in ['pdb', 'cif']:
        return go.Figure(), {'display' : 'none'}, {'display' : 'block'}, {'display' : 'none'}, {'display' : 'none'}

    content_type, content_string = data.get('url').split(',')
    decoded = base64.b64decode(content_string)
    file_base64 = base64.b64encode(decoded).decode()
    file_data_url = f"data:{content_type};base64,{file_base64}"

    if file_ext == 'pdb':
        parser = PDBParser()
        structure = parser.get_structure(id=filename.split('.')[0], file=StringIO(decoded.decode('utf-8')))
    elif file_ext == 'cif':
        parser = MMCIFParser()
        structure = parser.get_structure(structure_id=filename.split('.')[0], filename=StringIO(decoded.decode('utf-8')))

    structure_type = check_nucleotide_type(decoded)
    if structure_type == "Other":
        return go.Figure(), {'display' : 'none'}, {'display' : 'block'}, {'display' : 'none'}, {'display' : 'none'}

    points = []
    nucleotide_info = []

    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.resname in rna_nucleotides or residue.resname in dna_nucleotides:
                    coord = []
                    for atom in residue:
                        coord.append(atom.get_coord())
                    coord_array = np.array(coord)
                    center = np.mean(coord_array, axis=0) 
                    points.append(center)                 
                    colors.append(color_map.get(residue.resname))
                    nucleotide_info.append( {
                        "Nucleotide" : residue.resname,
                        "Chain_id" : chain.id,
                        "Coordinate" : center,
                        "Color" : color_map.get(residue.resname),
                        "Nucleotide_id" : residue.id[1]
                    })

    points_array = np.array(points)

    fig = go.Figure(data=[go.Scatter3d(
        x=points_array[:, 0],  
        y=points_array[:, 1],  
        z=points_array[:, 2],  
        mode='markers',
        marker=dict(
            size=6,
            color=colors,
            opacity=1.0,
            line=dict(
                color='white',
                width=0.5
            )
        ),
        hoverinfo='text',
        hovertemplate="Nucleotide: %{customdata[0]} %{customdata[4]}<br>Chain: %{customdata[1]}<br>Coordinates:<br> x: %{customdata[2][0]:.2f}<br> y: %{customdata[2][1]:.2f}<br> z: %{customdata[2][2]:.2f}<extra></extra>",
        customdata=[
            [nucleotide['Nucleotide'], nucleotide['Chain_id'], nucleotide['Coordinate'], nucleotide['Color'], nucleotide['Nucleotide_id']] 
            for nucleotide in nucleotide_info
        ],
    )])
    fig.update_scenes(xaxis_showspikes=False, yaxis_showspikes=False, zaxis_showspikes=False)
    fig.update_layout(
        autosize=True,
        scene=dict(
            xaxis=dict(visible=False),
            yaxis=dict(visible=False),
            zaxis=dict(visible=False)
        ),
        height=610,
        width=1000,
        margin=dict(l=0, r=0, t=0, b=0),
        paper_bgcolor='#fafafb',
        clickmode='event+select',
        dragmode="select",
        newselection_mode="gradual",
    )

    return fig, {'display' : 'block'}, {'display': 'none'}, {'display' : 'block'}, {'display' : 'none'}

@callback(
    [
        Output('nucleotide-input', 'children'), 
        Output('chain-input', 'children'),       
        Output('sel-nucleotide-info', 'style', allow_duplicate=True),    
        Output('rna-graph', 'figure', allow_duplicate=True)              
    ],
    Input('rna-graph', 'clickData'),
    State('rna-graph', 'figure'),
    State('rna-graph', 'relayoutData'),
    prevent_initial_call=True
)
def display_selected_info(clickData, current_figure, relayoutData): 
    fig = go.Figure(current_figure)

    if clickData:
        point_data = clickData['points'][0]['customdata']
        if point_data:
            nucleotide_value = f"{point_data[0]}{point_data[4]}"
            chain_value = f"{point_data[1]}"

            updated_colors = []
            for point in fig.data[0].customdata:
                original_color = point[3]
                if point in [point['customdata'] for point in clickData['points']]:
                    updated_colors.append(original_color)
                    fig.data[0].marker.size = [12 if p in [point['customdata'] for point in clickData['points']] else 10 for p in fig.data[0].customdata]
                else:
                    if original_color == 'green':
                        updated_colors.append('rgba(0, 100, 0, 0.6)')  
                    elif original_color == 'blue':
                        updated_colors.append('rgba(0, 0, 255, 0.6)')  
                    elif original_color == 'red':
                        updated_colors.append('rgba(255, 0, 0, 0.6)')  
                    elif original_color == 'yellow':
                        updated_colors.append('rgba(255, 255, 0, 0.6)')
                    elif original_color == 'purple':
                        updated_colors.append('rgba(155, 39, 176, 0.6)')  
                    else:
                        updated_colors.append('rgba(128, 128, 128, 0.8)')

            fig.data[0].marker.color = updated_colors
            if relayoutData and 'scene.camera' in relayoutData:
                fig.update_layout(scene_camera=relayoutData['scene.camera'])
            return nucleotide_value, chain_value, {'display': 'block'}, fig        
        else:
            return '', '', {'display': 'none'}, fig
    else:
        return '', '', {'display': 'none'}, fig

@callback(
    [
        Output('coordinates-container', 'style'), 
        Output('coord-x', 'children'),             
        Output('coord-y', 'children'),         
        Output('coord-z', 'children'),            
        Output('toggle-coordinates-btn', 'children') 
    ],
    [
        Input('toggle-coordinates-btn', 'n_clicks'), 
        Input('rna-graph', 'clickData')             
    ],
    [
        State('coordinates-container', 'style'), 
        State('coord-x', 'children'),
        State('coord-y', 'children'),
        State('coord-z', 'children')
    ],
    prevent_initial_call=True
)
def toggle_coordinates_visibility(n_clicks, clickData, current_style, current_x, current_y, current_z):
    if n_clicks is None or current_style['display'] == 'block':
        return {'display': 'none'}, current_x, current_y, current_z, '▼'
    elif current_style['display'] == 'none':
        point_data = clickData['points'][0]['customdata']

        coord_x = f"x: {point_data[2][0]:.2f}"
        coord_y = f"y: {point_data[2][1]:.2f}"
        coord_z = f"z: {point_data[2][2]:.2f}"
        return {'display': 'block'}, coord_x, coord_y, coord_z, '▲' 

@callback(
    Output('sel-nucleotide-info', 'style', allow_duplicate=True),
    Output('rna-graph', 'figure', allow_duplicate=True),
    Input('clear-button', 'n_clicks'),
    State('rna-graph', 'figure'),
    prevent_initial_call=True
)
def clear_selection(n_clicks, current_figure):
    if n_clicks is not None:
        fig = go.Figure(current_figure)
        fig.data[0].marker.color = colors  
        fig.data[0].marker.size = 6
        return {'display': 'none'}, fig  
    else:
        return {'display': 'block'}, current_figure

