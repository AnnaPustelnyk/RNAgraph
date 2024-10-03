import dash
from dash import dcc, html, callback, Output, Input, State
import base64

# Register the page
dash.register_page(__name__, path='/', name="Mol* Viewer")

# Layout of the page
layout = html.Div(
    [
        html.H1("Welcome to RNA Graph application", className='h1', style={'textAlign': 'center'}),
        html.Div(id='molstar-viewer-container'),
        html.Div(id='upload-message', style={'color': 'red'}),
    ]
)

@callback(
    Output('store', 'data'),
    Output('upload-message', 'children'),
    Input('upload-data', 'contents'),
    State('upload-data', 'filename')
)
def update_molstar_view(contents, filename):
    if contents is None or filename is None:
        return dash.no_update, dash.no_update

    content_type, content_string = contents.split(',')
    decoded = base64.b64decode(content_string)
    file_ext = filename.split('.')[-1].lower()

    if file_ext not in ['pdb', 'cif']:
        return dash.no_update, "Error: Unsupported file format. Please upload a PDB or CIF file."
    
    file_base64 = base64.b64encode(decoded).decode()
    file_data_url = f"data:{content_type};base64,{file_base64}"

    return {'url': file_data_url, 'ext': file_ext}, dash.no_update

dash.clientside_callback(
    """
    function(data) {
        console.log("File data received:", data);
        if (!data) {
            console.error("No file data provided.");
            return;
        }

        // Initialize Mol* Viewer
        var viewerElement = document.getElementById('molstar-viewer-container');
        viewerElement.innerHTML = '';  // Clear previous viewer content

        var file_url = data.url;
        var extension = data.ext;
        var structureFormat;
        if (extension === 'cif') {
            structureFormat = 'mmcif';  // For CIF files
        } else if (extension === 'pdb') {
            structureFormat = 'pdb';  // For PDB files
        } else {
            console.error("Unsupported file format:", extension);
            return;  // Exit if format is unsupported
        }

        molstar.Viewer.create(viewerElement, {
            layoutShowControls: false,
            viewportShowExpand: false,
            collapseLeftPanel: false,
            pdbProvider: 'pdbe',
            emdbProvider: 'pdbe'
        }).then(viewerInstance => {
            viewerInstance.loadStructureFromUrl(file_url, structureFormat)
                .catch(error => console.error("Error loading structure:", error));
        });
    }
    """,
    Output('molstar-viewer-container', 'children'),
    Input('store', 'data')
)
