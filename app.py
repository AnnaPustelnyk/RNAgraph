import base64
import dash
from dash import dcc, html
import dash_bootstrap_components as dbc

app = dash.Dash(__name__, use_pages=True, external_stylesheets=[dbc.themes.BOOTSTRAP], suppress_callback_exceptions=True)

app.layout = html.Div(
    [
        html.Div(
            [
                dcc.Location(id='url', refresh=False),
                html.H1("RNA Graph", style={'textAlign': 'left'}),
                dcc.Upload(
                    id='upload-data',
                    className = 'upload-data',
                    children=html.Div(['Drag and Drop or ', html.A('Select a PDB or CIF File')]),
                    multiple=False,
                    style = {'font-size': '16px'}
                ),
                dcc.Store(id="store"),
                dbc.NavbarSimple(
                    children=[
                        dbc.NavItem(
                            dcc.Link(page['name'], href=page['path'], id=f"{page['name'].lower()}-link", className="nav-link"),  
                        )
                        for page in dash.page_registry.values()
                    ],
                    color = '#1e3d5c',
                )
            ],
            className="header", 
        ),
        html.Div(dash.page_container, className="content"),
        html.Div(id='upload-message', style={'color': 'red'}),
    ],
    className="app-container",
)

@app.callback(
    [dash.dependencies.Output(f"{page['name'].lower()}-link", 'className') for page in dash.page_registry.values()],
    [dash.dependencies.Input('upload-data', 'contents'), dash.dependencies.Input('url', 'pathname')]
)
def update_active_link(contents, pathname):
    molviewer_class = 'nav-link'
    RNAgraph_class = 'nav-link'

    if pathname == '/':
        molviewer_class = 'nav-link active'

    if pathname == '/page-2':
        RNAgraph_class = 'nav-link active'

    return molviewer_class, RNAgraph_class

if __name__ == "__main__":
    app.run_server(debug=True)