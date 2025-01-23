import dash
from dash import dcc, html, callback, Output, Input, State, Patch
import base64
import numpy as np
import plotly.graph_objects as go
import rnapolis
from rnapolis import annotator, parser
from Bio.PDB import PDBParser, MMCIFParser
from io import BytesIO, StringIO
from dash.exceptions import PreventUpdate
import dash_bootstrap_components as dbc
import dash_daq as daq
import os

rna_nucleotides = ['A', 'C', 'G', 'U', 'I']
dna_nucleotides = ['DA', 'DC', 'DG', 'DU', 'DI', 'DT']

colors = []
color_map = {
        'A': 'rgb(225, 246, 0)',
        'C': 'rgb(35, 120, 65)',
        'G': 'rgb(228, 34, 23)',
        'U': 'rgb(65, 105, 225)',
        'I': 'rgb(127, 127, 127)',
        'DA': 'rgb(225, 246, 0)',
        'DC': 'rgb(35, 120, 65)',
        'DG': 'rgb(228, 34, 23)',
        'DU': 'rgb(65, 105, 225)',
        'DI': 'rgb(127, 127, 127)',
        'DT': 'rgb(128, 0, 128)'
    }

try:
    dash.register_page(__name__, path='/page-2', name="RNA Graph")
except dash.exceptions.PageError:
    pass

layout = html.Div( 
    [
        # sel-nucleotide-info container
        html.Div(
            id='sel-nucleotide-info',
            className='sel-nucleotide-info',
            children=[
                html.Div(
                    [
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
                                                'fontSize': '18px', 
                                                'cursor': 'pointer'
                                            }
                                        ),
                                    ],
                                    style={'display': 'flex', 'alignItems': 'start'}
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
                    ],
                    className='top-section'
                ), 
                html.Button('Clear', id='clear-button', className='clear-button'), 
            ], 
        ),
        # Main container with the graph
        html.Div(
            [
                dcc.Loading(
                    id = 'loading-1',
                    type = 'circle',
                    overlay_style={"visibility":"visible"},
                    children = dcc.Graph(
                                id='rna-graph',
                                style={'display': 'block'},  
                                className='graph'
                            )
                ),
                dbc.DropdownMenu(
                    label="Colors",  
                    children=[
                        dbc.DropdownMenuItem("Sequence", id="seq"),
                        dbc.DropdownMenuItem("Optional", id="opt"),
                    ],
                    id = "dropdown-menu",
                    color = "primary",
                    direction = "up",
                    toggle_style = {
                        "position": "absolute",
                        "left": "0px",
                        "bottom": "0px",
                        "width": "216px", 
                        "display": "flex",
                        "justify-content": "space-between",
                        "alignItems": "center",
                        "fontFamily": "'roboto', sans-serif",
                        "borderRadius": "12px",
                        "borderColor": "#dadadd",
                        "backgroundColor": "#fafafb",
                        "color": "#1f1f20",
                    },
                ),
                html.Div(
                    [
                        daq.ColorPicker(
                            id='color-picker',
                             
                        ),
                        html.Button('Close', id='close-button', className='clear-button', style={'marginTop': '8px'}),
                    ],
                    id = 'color-picker-container',
                    className = 'color-picker-container',
                    style = {'display': 'none'}
                ),
                html.Div(
                    id = 'structure-name',
                    className = 'structure-name',
                    style = {'display': 'none'}
                )  
            ],
            className='rna-graph-container',
            id = 'rna-graph-container' 
        ),  
        html.Div(
            id='side-bar',
            className='side-bar',
            children=[
                html.Div(
                    [    
                        html.P(
                            "Nucleotide Interactions",
                            id='title2',
                            className='title'
                        ),
                        dcc.Checklist(
                            id='interaction-type',
                            className='checklist-label',
                            options=[
                                {'label': 'Phosphodiester interactions', 'value': 'phosphodiester'},
                                {'label': 'Canonical interactions', 'value': 'c_base_base'},
                                {'label': 'Non-canonical interactions', 'value': 'nc_base_base'},
                                {'label': 'Stacking interactions', 'value': 'stacking'}
                            ],
                            value=[],  
                            inputClassName='checklist-input'
                        ),
                        dcc.Checklist(
                            id='heteroatoms-show',
                            className='checklist-label',
                            options=[
                                {'label': 'Show heteroatoms', 'value': 'heteroatoms'},
                            ], 
                            value=[],
                            inputClassName='checklist2-input',
                        ),
                    ],
                    className='top-section'
                ),
                html.Div(
                    [
                        html.Div(
                            [
                                html.Div(
                                    className='interaction-description',
                                    children = [
                                        html.Hr(id = 'phosphodiester-hr', style={'borderWidth': '2px', 'width': '44px', 'borderColor': 'green', 'opacity': 'unset', 'borderStyle': 'dashed'}),
                                        html.Label('Phosphodiester interactions', style = {'fontSize': '14px'})
                                    ]
                                ),
                                html.Div(
                                    [
                                        dcc.Dropdown(
                                            placeholder="Color",  
                                            options=[
                                                {'label': 'Red', 'value': 'red'},
                                                {'label': 'Green', 'value': 'green'},
                                                {'label': 'Blue', 'value': 'blue'},
                                                {'label': 'Orange', 'value': 'orange'},
                                                {'label': 'Yellow', 'value': 'yellow'},
                                                {'label': 'Violet', 'value': 'violet'},
                                                {'label': 'Gray', 'value': 'gray'},
                                                {'label': 'Black', 'value': 'black'}
                                            ],
                                            value = 'green',
                                            className = 'dropUp',
                                            id = 'phosphodiester-color',
                                            searchable = False,
                                            clearable= False,
                                        ),
                                        dcc.Dropdown(
                                            placeholder="Style",  
                                            options=[
                                                {'label': 'Solid', 'value': 'solid'},
                                                {'label': 'Dashed', 'value': 'dash'},
                                                {'label': 'Dashed Dot', 'value': 'longdash'}
                                            ],
                                            value = 'longdash',
                                            className = 'dropUp',
                                            id = 'phosphodiester-style',
                                            searchable = False,
                                            clearable= False,
                                        ),
                                    ],
                                    className='interaction-dropdown-container'  
                                ),
                            ],
                            className='interaction-container',
                            id = 'phosphodiester-style-container'
                        ),
                        html.Div(
                            [
                                html.Div(
                                    className='interaction-description',
                                    children = [
                                        html.Hr(id = 'canonical-hr', style={'borderWidth': '2px', 'width': '44px', 'borderColor': 'blue', 'opacity': 'unset'}),
                                        html.Label('Canonical interactions', style = {'fontSize': '14px'})
                                    ]
                                ),
                                html.Div(
                                    [
                                        dcc.Dropdown(
                                            placeholder="Color",  
                                            options=[
                                                {'label': 'Red', 'value': 'red'},
                                                {'label': 'Green', 'value': 'green'},
                                                {'label': 'Blue', 'value': 'blue'},
                                                {'label': 'Orange', 'value': 'orange'},
                                                {'label': 'Yellow', 'value': 'yellow'},
                                                {'label': 'Violet', 'value': 'violet'},
                                                {'label': 'Gray', 'value': 'gray'},
                                                {'label': 'Black', 'value': 'black'}
                                            ],
                                            value = 'blue',
                                            className = 'dropUp',
                                            id = 'canonical-color',
                                            searchable = False,
                                            clearable= False,
                                        ),
                                        dcc.Dropdown(
                                            placeholder="Style",  
                                            options=[
                                                {'label': 'Solid', 'value': 'solid'},
                                                {'label': 'Dashed', 'value': 'dash'},
                                                {'label': 'Dashed Dot', 'value': 'longdash'}
                                            ],
                                            value = 'solid',
                                            className = 'dropUp',
                                            id = 'canonical-style',
                                            searchable = False,
                                            clearable= False,
                                        ),
                                    ],
                                    className='interaction-dropdown-container'  
                                ),
                            ],
                            className='interaction-container',
                            id = 'canonical-style-container'
                        ),
                        html.Div(
                            [
                                html.Div(
                                    className='interaction-description',
                                    children = [
                                        html.Hr(id = 'non-canonical-hr', style={'borderWidth': '2px', 'width': '44px', 'borderColor': 'black', 'opacity': 'unset', 'borderStyle': 'dashdot'}),
                                        html.Label('Non-canonical interactions', style = {'fontSize': '14px'})
                                    ]
                                ),
                                html.Div(
                                    [
                                        dcc.Dropdown(
                                            placeholder="Color",  
                                            options=[
                                                {'label': 'Red', 'value': 'red'},
                                                {'label': 'Green', 'value': 'green'},
                                                {'label': 'Blue', 'value': 'blue'},
                                                {'label': 'Orange', 'value': 'orange'},
                                                {'label': 'Yellow', 'value': 'yellow'},
                                                {'label': 'Violet', 'value': 'violet'},
                                                {'label': 'Gray', 'value': 'gray'},
                                                {'label': 'Black', 'value': 'black'}
                                            ],
                                            value = 'black',
                                            className = 'dropUp',
                                            id = 'non-canonical-color',
                                            searchable = False,
                                            clearable= False,
                                        ),
                                        dcc.Dropdown(
                                            placeholder="Style",  
                                            options=[
                                                {'label': 'Solid', 'value': 'solid'},
                                                {'label': 'Dashed', 'value': 'dash'},
                                                {'label': 'Dashed Dot', 'value': 'longdash'}
                                            ],
                                            value = 'solid',
                                            className = 'dropUp',
                                            id = 'non-canonical-style',
                                            searchable = False,
                                            clearable= False,
                                        ),
                                    ],
                                    className='interaction-dropdown-container'  
                                ),
                            ],
                            className='interaction-container',
                            id = 'non-canonical-style-container'
                        ),
                        html.Div(
                            [
                                html.Div(
                                    className='interaction-description',
                                    children = [
                                        html.Hr(id = 'stacking-hr', style={'borderWidth': '2px', 'width': '44px', 'borderColor': 'orange', 'opacity': 'unset', 'borderStyle': 'dashed'}),
                                        html.Label('Stacking interactions', style = {'fontSize': '14px'})
                                    ]
                                ),
                                html.Div(
                                    [
                                        dcc.Dropdown(
                                            placeholder="Color",  
                                            options=[
                                                {'label': 'Red', 'value': 'red'},
                                                {'label': 'Green', 'value': 'green'},
                                                {'label': 'Blue', 'value': 'blue'},
                                                {'label': 'Orange', 'value': 'orange'},
                                                {'label': 'Yellow', 'value': 'yellow'},
                                                {'label': 'Violet', 'value': 'violet'},
                                                {'label': 'Gray', 'value': 'gray'},
                                                {'label': 'Black', 'value': 'black'}
                                            ],
                                            value = 'orange',
                                            className = 'dropUp',
                                            id = 'stacking-color',
                                            searchable = False,
                                            clearable= False,
                                        ),
                                        dcc.Dropdown(
                                            placeholder="Style",  
                                            options=[
                                                {'label': 'Solid', 'value': 'solid'},
                                                {'label': 'Dashed', 'value': 'dash'},
                                                {'label': 'Dashed Dot', 'value': 'longdash'}
                                            ],
                                            value = 'longdash',
                                            className = 'dropUp',
                                            id = 'stacking-style',
                                            searchable = False,
                                            clearable= False,
                                        ),
                                    ],
                                    className='interaction-dropdown-container'  
                                ),
                            ],
                            className='interaction-container',
                            id = 'stacking-style-container'
                        ),
                    ],
                    className='bottom-section'
                ),
            ]
        )   
    ],
    id='contant',
    className='contant'
)


@callback(
    Output('rna-graph', 'figure', allow_duplicate=True),
    Output('rna-graph-container', 'style'),
    Output('structure-name', 'children'),
    Output('structure-name', 'style'),
    Output('contant', 'style'),
    Output('heteroatoms-show', 'options'),
    Input('store', 'data'),
    State('upload-data', 'filename'),
    prevent_initial_call = True
)
def update_rna_graph(data, filename):
    if data is None or filename is None:
        return go.Figure(), None, None, {'display' : 'none'}, {'display' : 'none'}, dash.no_update
    
    colors.clear()
    option = [{'label': 'Show heteroatoms', 'value': 'heteroatoms', 'disabled': False}]

    content_type, content_string = data.get('url').split(',')
    decoded = base64.b64decode(content_string)

    if 'pdb' in filename:
        parser = PDBParser()
        structure = parser.get_structure(id=filename.split('.')[0], file=StringIO(decoded.decode('utf-8')))
    elif 'cif' in filename:
        parser = MMCIFParser()
        structure = parser.get_structure(structure_id=filename.split('.')[0], filename=StringIO(decoded.decode('utf-8')))

    points = []
    heteroatoms =  []
    nucleotide_info = []
    heteroatom_info = []

    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.resname in rna_nucleotides or residue.resname in dna_nucleotides:
                    coord = []
                    heterocoord = []

                    for atom in residue:
                        coord.append(atom.get_coord())
                    
                    coord_array = np.array(coord)
                    center = np.mean(coord_array, axis=0) 
                    points.append(center)  
                    colors.append(color_map.get(residue.resname, 'rgb(16, 16, 16)'))

                    nucleotide_info.append( {
                        "Nucleotide" : residue.resname,
                        "Chain_id" : chain.id,
                        "Coordinate" : center,
                        "Color" : color_map.get(residue.resname, 'rgb(16, 16, 16)'),
                        "Nucleotide_id" : residue.id[1]
                    })
                if 'H_' in residue.id[0]:
                    heterocoord = []
                    for atom in residue:
                        heterocoord.append(atom.get_coord())
                    
                    heteroatom_coord_array = np.array(heterocoord)
                    center = np.mean(heteroatom_coord_array, axis=0) 
                    heteroatoms.append(center)

                    heteroatom_info.append( {
                        "Nucleotide" : residue.resname,
                        "Chain_id" : chain.id,
                        "Coordinate" : center,
                        "Color" : 'black',
                        "Nucleotide_id" : residue.id[1]
                    })
        break
    points_array = np.array(points)
    if points_array.size == 0:
        return go.Figure(), dash.no_update, {'display': 'none'}, {'display' : 'none'}, {'display': 'none'}, dash.no_update
    
    if len(heteroatoms) > 0:
        heteroatoms_array = np.array(heteroatoms)
    else: 
        heteroatoms_array = None
        option = [{'label': 'Show heteroatoms', 'value': 'heteroatoms', 'disabled': True}]

    fig = go.Figure(data=[go.Scatter3d(
        name = 'nucleotides',
        x=points_array[:, 0],  
        y=points_array[:, 1],  
        z=points_array[:, 2],  
        mode='markers',
        marker=dict(
            size=8,
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

    if heteroatoms_array is not None:
        fig.add_trace(go.Scatter3d(
            name='heteroatoms',
            x=heteroatoms_array[:, 0],
            y=heteroatoms_array[:, 1],
            z=heteroatoms_array[:, 2],
            mode='markers',
            marker=dict(
                size=6,
                color='black',
                opacity=0.8,
                line=dict(
                    color='black',
                    width=0.5
                )
            ),
            hoverinfo='text',
            hovertemplate="Heteroatom: %{customdata[0]}<br>Chain: %{customdata[1]}<br>Coordinates:<br>x: %{customdata[2][0]:.2f}<br>y: %{customdata[2][1]:.2f}<br>z: %{customdata[2][2]:.2f}<extra></extra>",
            customdata=[
                [heteroatom['Nucleotide'], heteroatom['Chain_id'], heteroatom['Coordinate'], heteroatom['Color'], heteroatom['Nucleotide_id']]
                for heteroatom in heteroatom_info
            ],
            visible=False
        ))
    fig.update_scenes(xaxis_showspikes=False, yaxis_showspikes=False, zaxis_showspikes=False)
    fig.update_layout(
        autosize=True,
        scene=dict(
            xaxis=dict(visible=False),
            yaxis=dict(visible=False),
            zaxis=dict(visible=False)
        ),
        margin=dict(l=0, r=0, t=0, b=0),
        height = 610,
        paper_bgcolor='#fafafb',
        clickmode='event+select',
        dragmode="select",
        newselection_mode="gradual",
        showlegend=False
    )

    structure_name = data.get('name')
    return fig, {'display' : 'block'}, structure_name, {'display' : 'block'}, {'display' : 'flex'}, option

@callback(
    [
        Output('nucleotide-input', 'children'), 
        Output('chain-input', 'children'),          
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
            updated_colors = [
                point[3] if point in [pt['customdata'] for pt in clickData['points']] else f'rgba({c[4:-1]}, 0.6)'
                    for point, c in zip(fig.data[0].customdata, colors)
            ]
            fig.data[0].marker.color = updated_colors
            fig.data[0].marker.size = [16 if p in [point['customdata'] for point in clickData['points']] else 12 for p in fig.data[0].customdata]

            if len(fig.data) > 1 and fig.data[1].name == 'heteroatoms':
                fig.data[1].marker.color = ['black' if point in [pt['customdata'] for pt in clickData['points']] else f'rgba(0, 0, 0, 0.6)' for point in fig.data[1].customdata]
                fig.data[1].marker.size = [16 if p in [point['customdata'] for point in clickData['points']] else 12 for p in fig.data[1].customdata]
            
            if relayoutData and 'scene.camera' in relayoutData:
                fig.update_layout(scene_camera=relayoutData['scene.camera'])
            return nucleotide_value, chain_value, fig        
        else:
            return PreventUpdate
    else:
        return PreventUpdate

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
        Input('rna-graph', 'clickData'),             
    ],
    [
        State('coordinates-container', 'style'), 
        State('coord-x', 'children'),
        State('coord-y', 'children'),
        State('coord-z', 'children'),
        State('nucleotide-input', 'children')
    ],
    prevent_initial_call=True
)
def toggle_coordinates_visibility(n_clicks, clickData, current_style, current_x, current_y, current_z, nucleotide_input):
        changed_id = [p['prop_id'] for p in dash.callback_context.triggered][0]
        if 'toggle-coordinates-btn' in changed_id and current_style['display'] == 'none':
            if clickData is not None and nucleotide_input != "":
                
                point_data = clickData['points'][0]['customdata']
                coord_x = f"x: {point_data[2][0]:.2f}"
                coord_y = f"y: {point_data[2][1]:.2f}"
                coord_z = f"z: {point_data[2][2]:.2f}"
                
                return {'display': 'block', 'borderStyle': 'solid', 'borderTopWidth': '0px', 'borderRightWidth': '0px', 'borderLeftWidth': '0px', 'borderBottomWidth': '1px', 'BorderColor': '#a6a6a6', 'paddingBottom': '8px'}, coord_x, coord_y, coord_z, '▲' 
            else:
                return {'display': 'none'}, "x: ", "y: ", "z: ", '▼'
        else:
            return {'display': 'none'}, current_x, current_y, current_z, '▼'

@callback(
    Output('rna-graph', 'figure', allow_duplicate=True),
    Output('nucleotide-input', 'children', allow_duplicate=True),
    Output('chain-input', 'children', allow_duplicate=True),
    Output('coord-x', 'children', allow_duplicate=True),
    Output('coord-y', 'children', allow_duplicate=True),
    Output('coord-z', 'children', allow_duplicate=True),
    Output('coordinates-container', 'style', allow_duplicate=True),
    Output('toggle-coordinates-btn', 'children', allow_duplicate=True),
    Input('clear-button', 'n_clicks'),
    State('rna-graph', 'figure'),
    State('rna-graph', 'relayoutData'),
    prevent_initial_call=True
)
def clear_selection(n_clicks, current_figure, relayoutData):
    if n_clicks is not None:
        fig = go.Figure(current_figure)
        fig.data[0].marker.color = colors 
        fig.data[0].marker.size = 6

        if len(fig.data) > 1 and fig.data[1].name == 'heteroatoms':
            fig.data[1].marker.color = 'black'
            fig.data[1].marker.size = 6

        if relayoutData and 'scene.camera' in relayoutData:
                fig.update_layout(scene_camera=relayoutData['scene.camera'])

        return fig, "", "", "x: ", "y: ", "z: ", {'display': 'none'}, '▼'
    else:
        return PreventUpdate
@callback(
    Output('rna-graph', 'figure', allow_duplicate=True),
    Output('heteroatoms-show', 'value'),
    Input('heteroatoms-show', 'value'),
    State('rna-graph', 'figure'),
    State('rna-graph', 'relayoutData'),
    prevent_initial_call=True
)
def show_heteroatoms(values, current_figure, relayoutData):
    current_figure = go.Figure(current_figure)
    
    if len(current_figure.data) > 1  and current_figure.data[1].name == 'heteroatoms':
        if 'heteroatoms' in values:
            current_figure.update_traces(visible=True, selector=dict(name='heteroatoms'))
            if relayoutData and 'scene.camera' in relayoutData:
                current_figure.update_layout(scene_camera=relayoutData['scene.camera'])
            return current_figure, ['heteroatoms']
        else:
            current_figure.update_traces(visible=False, selector=dict(name='heteroatoms'))
            if relayoutData and 'scene.camera' in relayoutData:
                current_figure.update_layout(scene_camera=relayoutData['scene.camera'])
            return current_figure, []
    else:
        return current_figure, []

@callback(
    Output('rna-graph', 'figure', allow_duplicate=True),
    Output('interaction-type', 'options'),
    Input('interaction-type', 'value'),
    Input('store', 'data'),
    Input('rna-graph', 'figure'),
    Input('processed-data', 'data'),
    State('rna-graph', 'relayoutData'),
    prevent_initial_call=True
)
def update_interaction_info(selected_interactions, data, current_figure, interactions, relayoutData):
    available_interactions = interactions if interactions else {
        'phosphodiester': [],
        'c_base_base': [],
        'nc_base_base': [],
        'stacking': []
    }

    interaction_options = [
        {'label': 'Phosphodiester interactions', 'value': 'phosphodiester', 'disabled': True},
        {'label': 'Canonical interactions', 'value': 'c_base_base', 'disabled': True},
        {'label': 'Non-canonical interactions', 'value': 'nc_base_base', 'disabled': True},
        {'label': 'Stacking interactions', 'value': 'stacking', 'disabled': True},
    ]

    current_figure = go.Figure(current_figure)

    if data is not None:   
        interaction_options = [
            {'label': 'Phosphodiester interactions', 'value': 'phosphodiester', 'disabled': not available_interactions['phosphodiester']},
            {'label': 'Canonical interactions', 'value': 'c_base_base', 'disabled': not available_interactions['c_base_base']},
            {'label': 'Non-canonical interactions', 'value': 'nc_base_base', 'disabled': not available_interactions['nc_base_base']},
            {'label': 'Stacking interactions', 'value': 'stacking', 'disabled': not available_interactions['stacking']},
        ]

        if selected_interactions:

            for i, trace in enumerate(current_figure['data']):
                if trace['name'] not in selected_interactions and trace['name'] != 'nucleotides' and trace['name'] != 'heteroatoms':
                    current_figure.update_traces(visible=False, selector=dict(name=trace['name']))
                elif trace['name'] in selected_interactions:
                    current_figure.update_traces(visible=True, selector=dict(name=trace['name']))

            for interaction_type in selected_interactions:
                existing_traces = [trace['name'] for trace in current_figure.data if trace['name'] == interaction_type]
                if interaction_type not in set(existing_traces):
                    interactions = available_interactions.get(interaction_type, [])
                
                    if interactions:
                        if len(current_figure.data) > 1 and current_figure.data[1].name == 'heteroatoms':
                            interaction_lines = create_interaction_lines(interactions, current_figure.data[0].customdata, current_figure.data[1].customdata, interaction_type)
                        else:
                            interaction_lines = create_interaction_lines(interactions, current_figure.data[0].customdata, None, interaction_type)
                        current_figure.add_traces(interaction_lines)
                        current_figure.update_layout(
                            scene=dict(
                                xaxis=dict(visible=False),
                                yaxis=dict(visible=False),
                                zaxis=dict(visible=False)
                            ),
                            showlegend=False
                        )
                        current_figure.update_traces(
                                name=interaction_type, 
                                selector=dict(name=interaction_type)
                        )

            if relayoutData and 'scene.camera' in relayoutData:
                current_figure.update_layout(scene_camera=relayoutData['scene.camera'])

            return current_figure, interaction_options
        
        else:
            if relayoutData and 'scene.camera' in relayoutData:
                current_figure.update_layout(scene_camera=relayoutData['scene.camera'])

            for i, trace in enumerate(current_figure['data']):
                if trace['name'] != 'nucleotides' and trace['name'] != 'heteroatoms':
                    current_figure.update_traces(visible = False, selector=dict(name=trace['name']))

            return current_figure, interaction_options
    else:
        if relayoutData and 'scene.camera' in relayoutData:
                current_figure.update_layout(scene_camera=relayoutData['scene.camera'])

        return current_figure, interaction_options

def create_interaction_lines(interaction_list, nucleotide_info = None, heteroatom_info = None, interaction_type = None):    
    lines_styles = {
        'nc_base_base': {'color': 'black', 'width': 2, 'dash': None},
        'c_base_base': {'color': 'blue', 'width': 2, 'dash': None}, 
        'phosphodiester': {'color': 'green', 'width': 6, 'dash': 'longdash'},
        'stacking': {'color': 'orange', 'width': 6, 'dash': 'longdash'}
    }

    if interaction_list:
        lines = []
        
        for pair in interaction_list:
            nt1 = pair['nt1']
            nt2 = pair['nt2']

            try:
                if nt1['auth']['name'] not in rna_nucleotides and nt1['auth']['name'] not in dna_nucleotides:
                    if heteroatom_info is not None:
                        nt1_info = next(nucleotide for nucleotide in heteroatom_info if nucleotide[4] == nt1['auth']['number'] and nucleotide[1] == nt1['auth']['chain'] and nucleotide[0] == nt1['auth']['name'])
                    else:
                        continue
                else:
                    nt1_info = next(nucleotide for nucleotide in nucleotide_info if nucleotide[4] == nt1['auth']['number'] and nucleotide[1] == nt1['auth']['chain'] and nucleotide[0] == nt1['auth']['name'])
            except StopIteration:
                continue
            
            try:
                if nt2['auth']['name'] not in rna_nucleotides and nt2['auth']['name'] not in dna_nucleotides:
                    if heteroatom_info is not None:
                        nt2_info = next(nucleotide for nucleotide in heteroatom_info if nucleotide[4] == nt2['auth']['number'] and nucleotide[1] == nt2['auth']['chain'] and nucleotide[0] == nt2['auth']['name'])
                    else:
                        continue
                else:
                    nt2_info = next(nucleotide for nucleotide in nucleotide_info if nucleotide[4] == nt2['auth']['number'] and nucleotide[1] == nt2['auth']['chain'] and nucleotide[0] == nt2['auth']['name'])
            except StopIteration:
                continue

            nt1_position = nt1_info[2]
            nt2_position = nt2_info[2]

            style = lines_styles.get(interaction_type, {'color': 'black', 'width': 1, 'dash': None} )

            lines.append(go.Scatter3d(
                x=[nt1_position[0], nt2_position[0]],
                y=[nt1_position[1], nt2_position[1]],
                z=[nt1_position[2], nt2_position[2]],
                mode='lines',
                hoverinfo='none',
                showlegend=False,
                name = interaction_type,
                line=dict(color=style['color'], width=style['width'], dash=style['dash']),
            ))
        return lines if lines else None
    return None

@callback(
    Output('rna-graph', 'figure', allow_duplicate=True),
    Output('color-picker-container', 'style', allow_duplicate = True),              
    Input('seq', 'n_clicks'),
    Input('opt', 'n_clicks'),
    State('rna-graph', 'figure'),
    State('rna-graph', 'relayoutData'),
    prevent_initial_call=True
)
def colors_change(seq_click, opt_click, current_figure,relayoutData): 
    button_id = [p['prop_id'] for p in dash.callback_context.triggered][0]
    style = {'display': 'none'}
    current_figure = go.Figure(current_figure)
    fig = Patch()
    if fig:
        if button_id == 'seq.n_clicks':
            fig.data[0].marker.color = colors
            if len(current_figure.data) > 1 and current_figure.data[1].name == 'heteroatoms':
                fig.data[1].marker.color = 'black'
        elif button_id == 'opt.n_clicks':
            style = {'display': 'flex', 'position' : 'absolute', 'bottom' : '48px', 'left' : '12px'}
        else:
            return dash.no_update, dash.no_update
        
        if relayoutData and 'scene.camera' in relayoutData:
            fig.layout.scene.camera = relayoutData['scene.camera']

        return fig, style
    else:
        return dash.no_update, dash.no_update

@callback(
    Output('rna-graph', 'figure', allow_duplicate=True),
    Input('color-picker', 'value'),
    State('rna-graph', 'figure'),
    State('rna-graph', 'relayoutData'),
    prevent_initial_call='initial_duplicate'
)        
def color_picker_output(color, current_figure, relayoutData):
    fig = Patch()
    current_figure = go.Figure(current_figure)
    if fig and color:
        hex_color = color['hex']
        fig.data[0].marker.color = hex_color

        if len(current_figure.data) > 1 and current_figure.data[1].name == 'heteroatoms':
            fig.data[1].marker.color = hex_color

        if relayoutData and 'scene.camera' in relayoutData:
            fig.layout.scene.camera = relayoutData['scene.camera']

        return fig
    else:
        return dash.no_update

@callback(
    Output('color-picker-container', 'style', allow_duplicate = True),
    Input('close-button', 'n_clicks'),
    prevent_initial_call='initial_duplicate'
)
def close_picker(n_clicks):
    if n_clicks is not None:
        return {'display': 'none'}
    else:
        return dash.no_update


@callback(
    Output('phosphodiester-style-container', 'style'),
    Output('canonical-style-container', 'style'),
    Output('non-canonical-style-container', 'style'),
    Output('stacking-style-container', 'style'),
    Input('interaction-type', 'value'),
)
def interactions_style(value):
    phodphodiester_dis = {'display' : 'none'}
    stacking_dis = {'display' : 'none'}
    canonical_dis = {'display' : 'none'}
    noncanonical_dis = {'display' : 'none'}
    stacking_dis = {'display' : 'none'}
    
    if 'phosphodiester' in value:
        phodphodiester_dis = {'display' : 'block'}
    if 'c_base_base' in value:
        canonical_dis = {'display' : 'block'}
    if 'nc_base_base' in value:
        noncanonical_dis = {'display' : 'block'}
    if 'stacking' in value:
        stacking_dis = {'display' : 'block'}

    return phodphodiester_dis, canonical_dis, noncanonical_dis, stacking_dis

@callback(
    Output('rna-graph', 'figure', allow_duplicate=True),
    Output('phosphodiester-hr', 'style'),
    Input('phosphodiester-color', 'value'),
    Input('phosphodiester-style', 'value'),
    State('rna-graph', 'figure'),
    State('rna-graph', 'relayoutData'),
    prevent_initial_call='initial_duplicate'
)        
def phosphodiester_style(color, style, figure, relayoutData):
    if not figure or 'data' not in figure or not figure['data']:
        return dash.no_update, dash.no_update

    fig = Patch()
    updated = False

    for i in range(len(figure['data'])):
        if figure['data'][i]['name'] == 'phosphodiester': 
            fig.data[i].line.color = color
            fig.data[i].line.dash = style
            fig.data[i].line.width = 6
            updated = True
    
    if relayoutData and 'scene.camera' in relayoutData:
        fig.layout.scene.camera = relayoutData['scene.camera']

    if updated:
        return fig, {'borderWidth': '2px', 'width': '44px', 'borderColor': color, 'opacity': 'unset', 'borderStyle': style}
    else:
        return dash.no_update, dash.no_update

@callback(
    Output('rna-graph', 'figure', allow_duplicate=True),
    Output('canonical-hr', 'style'),
    Input('canonical-color', 'value'),
    Input('canonical-style', 'value'),
    State('rna-graph', 'figure'),
    State('rna-graph', 'relayoutData'),
    prevent_initial_call='initial_duplicate'
)        
def canonical_style(color, style, figure, relayoutData):
    if not figure or 'data' not in figure or not figure['data']:
        return dash.no_update, dash.no_update

    fig = Patch()
    updated = False

    for i in range(len(figure['data'])):
        if figure['data'][i]['name'] == 'c_base_base': 
            fig.data[i].line.color = color
            fig.data[i].line.dash = style
            fig.data[i].line.width = 4
            updated = True
    
    if relayoutData and 'scene.camera' in relayoutData:
        fig.layout.scene.camera = relayoutData['scene.camera']

    if updated:
        return fig, {'borderWidth': '2px', 'width': '44px', 'borderColor': color, 'opacity': 'unset', 'borderStyle': style}
    else:
        return dash.no_update, dash.no_update

@callback(
    Output('rna-graph', 'figure', allow_duplicate=True),
    Output('non-canonical-hr', 'style'),
    Input('non-canonical-color', 'value'),
    Input('non-canonical-style', 'value'),
    State('rna-graph', 'figure'),
    State('rna-graph', 'relayoutData'),
    prevent_initial_call='initial_duplicate'
)        
def non_canonical_style(color, style, figure, relayoutData):
    if not figure or 'data' not in figure or not figure['data']:
        return dash.no_update, dash.no_update

    fig = Patch()
    updated = False

    for i in range(len(figure['data'])):
        if figure['data'][i]['name'] == 'nc_base_base': 
            fig.data[i].line.color = color
            fig.data[i].line.dash = style
            fig.data[i].line.width = 4
            updated = True
    
    if relayoutData and 'scene.camera' in relayoutData:
        fig.layout.scene.camera = relayoutData['scene.camera']

    if updated:
        return fig, {'borderWidth': '2px', 'width': '44px', 'borderColor': color, 'opacity': 'unset', 'borderStyle': style}
    else:
        return dash.no_update, dash.no_update

@callback(
    Output('rna-graph', 'figure', allow_duplicate=True),
    Output('stacking-hr', 'style'),
    Input('stacking-color', 'value'),
    Input('stacking-style', 'value'),
    State('rna-graph', 'figure'),
    State('rna-graph', 'relayoutData'),
    prevent_initial_call='initial_duplicate'
)        
def stacking_style(color, style, figure, relayoutData):
    if not figure or 'data' not in figure or not figure['data']:
        return dash.no_update, dash.no_update

    fig = Patch()
    updated = False

    for i in range(len(figure['data'])):
        if figure['data'][i]['name'] == 'stacking': 
            fig.data[i].line.color = color
            fig.data[i].line.dash = style
            fig.data[i].line.width = 6
            updated = True
    
    if relayoutData and 'scene.camera' in relayoutData:
        fig.layout.scene.camera = relayoutData['scene.camera']

    if updated:
        return fig, {'borderWidth': '2px', 'width': '44px', 'borderColor': color, 'opacity': 'unset', 'borderStyle': style}
    else:
        return dash.no_update, dash.no_update

