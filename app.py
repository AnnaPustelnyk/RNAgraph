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
                    children=html.Div(['Drag and Drop or ', html.A('Select a PDB or CIF File')]),
                    style={
                        'height': '40px', 
                        'lineHeight': '24px',
                        'borderWidth': '1px', 
                        'borderStyle': 'dashed',
                        'borderRadius': '5px', 
                        'textAlign': 'center', 
                        'margin': '16px',
                        'padding': '10px'
                    },
                    multiple=False
                ),
                dcc.Store(id="store"),
                dbc.NavbarSimple(
                    children=[
                        dbc.NavItem(dcc.Link(page['name'], href=page['path'], id=f"{page['name'].lower()}-link", className="nav-link"))
                        for page in dash.page_registry.values()
                    ],
                    className="navbar-right"
                ),
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

    if contents:
        if pathname == '/':
            molviewer_class = 'nav-link active'

    if pathname == '/page-2':
        RNAgraph_class = 'nav-link active'

    return molviewer_class, RNAgraph_class

if __name__ == "__main__":
    app.run(debug=True)
