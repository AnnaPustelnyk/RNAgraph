import dash
from dash import dcc, html, callback, Output, Input, State
import base64
import numpy as np
import plotly.graph_objects as go
from Bio.PDB import PDBParser, MMCIFParser
from io import BytesIO, StringIO

colors = []
color_map = {
        'C': 'green',
        'N': 'blue',
        'O': 'red',
        'S': 'gray'
    }

dash.register_page(__name__, path='/page-2', name="RNA Graph")

layout = html.Div(
    style={'display': 'flex'},  
    children=[
        html.Div(
            id='sel-atom-info',
            style={'margin': '20px', 'display': 'none', 'font-size' : '12px'}
        ),
        html.Div(
            [
                html.H1("Welcome to RNA Graph application", id='title', className='h1', style={'textAlign': 'center'}),
                dcc.Graph(id='rna-graph', style={'display': 'none'}, className='graph'),
                html.Div(id='upload-message', style={'color': 'red'}),
            ],
            style={'flex-grow': '1'}  
        ),
    ]
)

@callback(
    Output('rna-graph', 'figure', allow_duplicate=True),
    Output('rna-graph', 'style'),
    Output('title', 'style'),
    Input('store', 'data'),
    State('upload-data', 'filename'),
    prevent_initial_call = True
)
def update_rna_graph(data, filename):
    if data is None or filename is None:
        return go.Figure(), {'display' : 'none'}, {'display': 'block'}

    file_ext = filename.split('.')[-1].lower()
    if file_ext not in ['pdb', 'cif']:
        return dash.no_update, "Error: Unsupported file format. Please upload a PDB or CIF file."

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

    points = []
    atom_info = []

    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    points.append(atom.get_coord())
                    colors.append(color_map.get(atom.element, 'gray'))
                    atom_info.append({
                        "Residue_name" : residue.resname,
                        "Atom_name" : atom.name,
                        "Chain_id" : chain.id,
                        "Coordinates" : atom.get_coord(),
                        "Color" : color_map.get(atom.element, 'gray')
                    })

    points_array = np.array(points)

    fig = go.Figure(data=[go.Scatter3d(
        x=points_array[:, 0],  
        y=points_array[:, 1],  
        z=points_array[:, 2],  
        mode='markers',
        marker=dict(
            size=4,
            color=colors,
            opacity=1.0,
        ),
        hoverinfo='text',
        hovertemplate="Residue: %{customdata[0]}<br>Atom: %{customdata[1]}<br>Chain: %{customdata[2]}<br>Coordinates:<br> x: %{customdata[3][0]:.2f}<br> y: %{customdata[3][1]:.2f}<br> z: %{customdata[3][2]:.2f}<extra></extra>",
        customdata=[
            [atom['Residue_name'], atom['Atom_name'], atom['Chain_id'], atom['Coordinates'], atom['Color']] 
            for atom in atom_info
        ],
    )])
    fig.update_scenes(xaxis_showspikes=False, yaxis_showspikes=False, zaxis_showspikes=False)
    fig.update_layout(
        scene=dict(
            xaxis=dict(visible=False),
            yaxis=dict(visible=False),
            zaxis=dict(visible=False)
        ),
        height=600,
        width=1200,
        margin=dict(l=0, r=0, t=0, b=0),
        clickmode='event+select',
        dragmode="select",
        newselection_mode="gradual",
    )

    return fig, {'display' : 'block'}, {'display': 'none'}

@callback(
    Output('sel-atom-info', 'children'),
    Output('sel-atom-info', 'style'),
    Output('rna-graph', 'figure', allow_duplicate=True),
    Input('rna-graph', 'clickData'),
    State('rna-graph', 'figure'),
    prevent_initial_call=True
)
def display_selected_info(clickData, current_figure): 
    fig = go.Figure(current_figure)

    if clickData:
        point_data = clickData['points'][0]['customdata']
        if point_data:
            info = [
                html.P("Selected Atom Information:", style={'margin': '0', 'font-weight': 'bold'}),
                html.P(f"Residue: {point_data[0]}", style={'margin': '0'}),
                html.P(f"Atom: {point_data[1]}", style={'margin': '0'}),
                html.P(f"Chain: {point_data[2]}", style={'margin': '0'}),
                html.P(f"Coordinates:", style={'margin': '0'}),
                html.P(f"x: {point_data[3][0]:.2f}", style={'margin': '0'}),
                html.P(f"y: {point_data[3][1]:.2f}", style={'margin': '0'}),
                html.P(f"z: {point_data[3][2]:.2f}", style={'margin': '0'}),
                html.Button('Clear', id='clear-button', className='clear-button', style={'display': 'block'})  # Include button here
            ]

            updated_colors = []
            for point in fig.data[0].customdata:
                original_color = point[4]
                if point in [point['customdata'] for point in clickData['points']]:
                    updated_colors.append(original_color)
                else:
                    if original_color == 'green':
                        updated_colors.append('rgba(0, 169, 0, 0.6)')  # Green with lower opacity
                    elif original_color == 'blue':
                        updated_colors.append('rgba(0, 0, 255, 0.6)')  # Blue with lower opacity
                    elif original_color == 'red':
                        updated_colors.append('rgba(255, 0, 0, 0.6)')  # Red with lower opacity
                    elif original_color == 'yellow':
                        updated_colors.append('rgba(255, 255, 0, 0.6)')  # Yellow with lower opacity
                    else:
                        updated_colors.append('rgba(128, 128, 128, 0.8)')

            fig.data[0].marker.color = updated_colors
            
            return info, {'display': 'block'}, fig
        else:
            return [], {'display': 'none'}, fig
    else:
        return [], {'display': 'none'}, fig
    
@callback(
    Output('sel-atom-info', 'style', allow_duplicate=True),
    Output('rna-graph', 'figure', allow_duplicate=True),
    Input('clear-button', 'n_clicks'),
    State('rna-graph', 'figure'),
    prevent_initial_call=True
)
def clear_selection(n_clicks, current_figure):
    if n_clicks is not None:
        fig = go.Figure(current_figure)
        fig.data[0].marker.color = colors  
        return {'display': 'none'}, fig  
    else:
        return {'display': 'block'}, current_figure

