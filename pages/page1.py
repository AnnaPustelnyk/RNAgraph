import dash
from dash import dcc, html, callback, Output, Input, State
import base64

dash.register_page(__name__, path='/', name="Mol* Viewer")

# Layout of the page
layout = html.Div(
    [
        html.Div(id='molstar-viewer-container'),
    ],
    style={'display': 'flex','justifyContent': 'center', 'alignItems': 'center', 'flexDirection' : 'column', 'height': '610px'}
)
dash.clientside_callback(
    """
    function(data) {
        console.log("File data received:", data);
        if (!data) {
            console.error("No file data provided.");
            return;
        }

        var viewerElement = document.getElementById('molstar-viewer-container');
        viewerElement.innerHTML = ''; 

        var file_url = data.url;
        var extension = data.ext;
        var structureFormat;
        if (extension === 'cif') {
            structureFormat = 'mmcif';  
        } else if (extension === 'pdb') {
            structureFormat = 'pdb';  
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
