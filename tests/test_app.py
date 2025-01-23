import base64
import pytest
import os
from app import update_active_link, check_nucleotide_type, calculate_interactions

# Load the sample PDB and CIF file contents for testing
def load_file_content(filename):
    with open(filename, 'rb') as f:
        return f.read()

# Paths to your files
pdb_file_path = 'tests/sample.pdb'  # Update this path
cif_file_path = 'tests/sample.cif'  # Update this path

# Sample file contents
sample_pdb_content = load_file_content(pdb_file_path)
sample_cif_content = load_file_content(cif_file_path)

def test_check_nucleotide_type():
    # Test PDB file for RNA
    pdb_data = base64.b64encode(sample_pdb_content).decode()
    assert check_nucleotide_type(base64.b64decode(pdb_data), 'pdb') == "RNA"

    # Test CIF file for RNA
    cif_data = base64.b64encode(sample_cif_content).decode()
    assert check_nucleotide_type(base64.b64decode(cif_data), 'cif') == "RNA"

    # Test invalid nucleotide type
    invalid_data = base64.b64encode(b"Invalid data").decode()
    assert check_nucleotide_type(base64.b64decode(invalid_data), 'pdb') == "Other"

def test_update_active_link():
    # Test with no file uploaded
    output = update_active_link(None, '/', None)
    assert output[0] is None  # upload-message
    assert output[1] is None  # store
    assert output[2] is None  # processed-data
    assert output[3] == {'display': 'block'}  # welcome-text style

    # Test with a valid PDB file    
    valid_pdb_data = f"data:text/plain;base64,{base64.b64encode(sample_pdb_content).decode()}"
    output = update_active_link(valid_pdb_data, '/', '1my9.pdb')
    assert output[0] is None  # upload-message
    assert output[1] is not None  # store should have data
    assert output[2] is not None  # processed-data should have interactions
    assert output[3] == {'display': 'none'}  # welcome-text style

    # Test with a valid CIF file upload
    valid_cif_data = f"data:text/plain;base64,{base64.b64encode(sample_cif_content).decode()}"
    output = update_active_link(valid_cif_data, '/', '1my9.cif')
    assert output[0] is None  # upload-message
    assert output[1] is not None  # store should have data
    assert output[2] is not None  # processed-data should have interactions
    assert output[3] == {'display': 'none'}  # welcome-text style

    # Test with an invalid file upload    
    invalid_file_data = "data:text/plain;base64,InvalidData"
    output = update_active_link(invalid_file_data, '/', 'test.txt')
    assert output[0] is not None  # upload-message should have an error
    assert output[1] is None  # store should be None
    assert output[2] is None  # processed-data should be None
    assert output[3] == {'display': 'none'}  # welcome-text style



def test_calculate_interactions():
    # Test the calculate_interactions function with a valid PDB file
    interactions = calculate_interactions(sample_pdb_content, 'pdb')
    assert interactions is not None  # Ensure interactions are returned
    assert 'phosphodiester' in interactions
    assert 'c_base_base' in interactions
    assert 'nc_base_base' in interactions
    assert 'stacking' in interactions

    # Test the calculate_interactions function with a valid CIF file
    interactions = calculate_interactions(sample_cif_content, 'cif')
    assert interactions is not None  # Ensure interactions are returned
    assert 'phosphodiester' in interactions
    assert 'c_base_base' in interactions
    assert 'nc_base_base' in interactions
    assert 'stacking' in interactions

    # Test the calculate_interactions function with invalid data
    invalid_data = b"Invalid content"
    interactions = calculate_interactions(invalid_data, 'pdb')
    assert interactions is None  # Should return None for invalid data

def test_edge_cases():
    # Test with an empty file
    empty_data = b""
    interactions = calculate_interactions(empty_data, 'pdb')
    assert interactions is None  # Should handle empty input gracefully

    # Test with a very large file (simulate large input)
    large_data = sample_pdb_content * 1000  # Simulating a large file
    interactions = calculate_interactions(large_data, 'pdb')
    assert interactions is not None  # Should handle large input without crashing

    # Test with a file having no RNA
    no_rna_data = base64.b64encode(b"Protein only content").decode()
    assert check_nucleotide_type(base64.b64decode(no_rna_data), 'pdb') == "Other"