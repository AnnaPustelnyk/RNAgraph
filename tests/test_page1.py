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

options = Options()
options.add_argument('--headless')
options.add_argument('--no-sandbox')
options.add_argument('--disable-dev-shm-usage')

chromedriver_path = "/usr/local/bin/chromedriver"  

service = Service(chromedriver_path)

driver = webdriver.Chrome(service=service, options=options)

pdb_file_path = 'tests/sample.pdb'
cif_file_path = 'tests/sample.cif'

def load_file_content(filename):
    with open(filename, 'rb') as f:
        return f.read()

sample_pdb_content = load_file_content(pdb_file_path)
sample_cif_content = load_file_content(cif_file_path)

@pytest.fixture
def app():
    app = import_app('app') 
    return app

def test_page1_layout(app):
    client = app.test_client()
    response = client.get('/')

    assert b"molstar-viewer-container" in response.data
    assert b"Welcome to RNA Graph application" in response.data

def test_molstar_callback(app):
    client = app.test_client()

    valid_pdb_data = f"data:text/plain;base64,{base64.b64encode(sample_pdb_content).decode()}"
    response = client.post('/', json={'upload-data': {'contents': valid_pdb_data, 'filename': 'sample.pdb'}})
    
    time.sleep(2) 

    response = client.get('/')
    assert b"Mol* Viewer" in response.data  

def test_invalid_file_format(app):
    client = app.test_client()

    invalid_file_data = "data:text/plain;base64,InvalidData"
    response = client.post('/', json={'upload-data': {'contents': invalid_file_data, 'filename': 'invalid.txt'}})

    assert b"Unsupported file format" in response.data
    assert b"molstar-viewer-container" not in response.data  

def test_cif_file_upload(app):
    client = app.test_client()

    valid_cif_data = f"data:text/plain;base64,{base64.b64encode(sample_cif_content).decode()}"
    response = client.post('/', json={'upload-data': {'contents': valid_cif_data, 'filename': 'sample.cif'}})

    time.sleep(2) 

    response = client.get('/')
    assert b"Mol* Viewer" in response.data 

@pytest.mark.parametrize(
    "file_data, expected_error",
    [
        ("data:text/plain;base64,InvalidData", "Invalid file format"),
        ("data:text/plain;base64,", "No file data provided")
    ]
)
def test_invalid_file_uploads(app, file_data, expected_error):
    client = app.test_client()

    response = client.post('/', json={'upload-data': {'contents': file_data, 'filename': 'test.txt'}})

    time.sleep(2)

    assert expected_error.encode() in response.data

def test_molstar_structure_load_failure(app):
    invalid_data = "data:text/plain;base64,InvalidStructureData"
    client = app.test_client()

    response = client.post('/', json={'upload-data': {'contents': invalid_data, 'filename': 'faulty.pdb'}})
    
    time.sleep(2)
    
    assert b"Error loading structure" in response.data  