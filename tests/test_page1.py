import base64
import pytest
from dash import Dash
from dash.testing.application_runners import import_app
from dash.testing.wait import until
import time
import os
from selenium import webdriver
from selenium.webdriver.chrome.options import Options
from selenium.webdriver.chrome.service import Service

# Set Chrome options
options = Options()
options.add_argument('--headless')
options.add_argument('--no-sandbox')
options.add_argument('--disable-dev-shm-usage')

# Specify the path to chromedriver
chromedriver_path = "/usr/local/bin/chromedriver"  # Adjust this path as needed

# Create a Service object with the chromedriver path
service = Service(chromedriver_path)

# Initialize WebDriver with the Service object
driver = webdriver.Chrome(service=service, options=options)

# Paths to your test files
pdb_file_path = 'tests/sample.pdb'
cif_file_path = 'tests/sample.cif'

# Sample file contents
def load_file_content(filename):
    with open(filename, 'rb') as f:
        return f.read()

sample_pdb_content = load_file_content(pdb_file_path)
sample_cif_content = load_file_content(cif_file_path)

@pytest.fixture
def app():
    # Initialize the Dash app for testing
    app = import_app('app')  # Ensure the app is imported correctly
    return app

def test_page1_layout(app):
    # Test if the page layout renders correctly
    client = app.test_client()
    response = client.get('/')

    assert b"molstar-viewer-container" in response.data
    assert b"Welcome to RNA Graph application" in response.data

def test_molstar_callback(app):
    # Test if the Mol* Viewer callback works correctly with a PDB file
    client = app.test_client()

    # Simulate uploading a valid PDB file
    valid_pdb_data = f"data:text/plain;base64,{base64.b64encode(sample_pdb_content).decode()}"
    response = client.post('/', json={'upload-data': {'contents': valid_pdb_data, 'filename': 'sample.pdb'}})
    
    # Wait for callback to complete and check if the viewer is populated
    time.sleep(2)  # Sleep to allow time for callback execution

    # Verify that the Mol* Viewer container contains the structure (check for structure loading)
    response = client.get('/')
    assert b"Mol* Viewer" in response.data  # Make sure the viewer content is visible (or other indicators of success)

def test_invalid_file_format(app):
    # Test the behavior when an unsupported file type is uploaded
    client = app.test_client()

    # Simulate uploading an invalid file
    invalid_file_data = "data:text/plain;base64,InvalidData"
    response = client.post('/', json={'upload-data': {'contents': invalid_file_data, 'filename': 'invalid.txt'}})

    # Verify that the error message is displayed
    assert b"Unsupported file format" in response.data  # Customize according to your app's error message
    assert b"molstar-viewer-container" not in response.data  # Ensure Mol* Viewer is not loaded with invalid data

def test_cif_file_upload(app):
    # Test the behavior when a valid CIF file is uploaded
    client = app.test_client()

    valid_cif_data = f"data:text/plain;base64,{base64.b64encode(sample_cif_content).decode()}"
    response = client.post('/', json={'upload-data': {'contents': valid_cif_data, 'filename': 'sample.cif'}})

    # Wait for callback to complete and check if the viewer is populated
    time.sleep(2)  # Sleep to allow time for callback execution

    # Verify that the Mol* Viewer container contains the structure (check for structure loading)
    response = client.get('/')
    assert b"Mol* Viewer" in response.data  # Ensure structure is visible in viewer

@pytest.mark.parametrize(
    "file_data, expected_error",
    [
        ("data:text/plain;base64,InvalidData", "Invalid file format"),
        ("data:text/plain;base64,", "No file data provided")
    ]
)
def test_invalid_file_uploads(app, file_data, expected_error):
    # Test invalid file uploads (both wrong format and empty file)
    client = app.test_client()

    response = client.post('/', json={'upload-data': {'contents': file_data, 'filename': 'test.txt'}})

    # Wait for the callback to execute
    time.sleep(2)

    # Check if the expected error message appears
    assert expected_error.encode() in response.data

def test_molstar_structure_load_failure(app):
    # Simulate a failure in loading the structure due to incorrect data
    invalid_data = "data:text/plain;base64,InvalidStructureData"
    client = app.test_client()

    response = client.post('/', json={'upload-data': {'contents': invalid_data, 'filename': 'faulty.pdb'}})
    
    # Wait for callback to complete
    time.sleep(2)
    
    # Check if the error message for loading failure is displayed
    assert b"Error loading structure" in response.data  # Customize the error message as needed