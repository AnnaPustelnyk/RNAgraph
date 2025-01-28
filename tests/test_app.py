import base64
import pytest
import os
from app import update_active_link, check_nucleotide_type, calculate_interactions

def load_file_content(filename):
    with open(filename, 'rb') as f:
        return f.read()

pdb_file_path = 'tests/sample.pdb'  
cif_file_path = 'tests/sample.cif' 

sample_pdb_content = load_file_content(pdb_file_path)
sample_cif_content = load_file_content(cif_file_path)

def test_check_nucleotide_type():
    pdb_data = base64.b64encode(sample_pdb_content).decode()
    assert check_nucleotide_type(base64.b64decode(pdb_data), 'pdb') == "RNA"

    cif_data = base64.b64encode(sample_cif_content).decode()
    assert check_nucleotide_type(base64.b64decode(cif_data), 'cif') == "RNA"

    invalid_data = base64.b64encode(b"Invalid data").decode()
    assert check_nucleotide_type(base64.b64decode(invalid_data), 'pdb') == "Other"

def test_update_active_link():
    output = update_active_link(None, '/', None)
    assert output[0] is None 
    assert output[1] is None 
    assert output[2] is None 
    assert output[3] == {'display': 'block'}  
    
    valid_pdb_data = f"data:text/plain;base64,{base64.b64encode(sample_pdb_content).decode()}"
    output = update_active_link(valid_pdb_data, '/', '1my9.pdb')
    assert output[0] is None  
    assert output[1] is not None  
    assert output[2] is not None  
    assert output[3] == {'display': 'none'} 

    valid_cif_data = f"data:text/plain;base64,{base64.b64encode(sample_cif_content).decode()}"
    output = update_active_link(valid_cif_data, '/', '1my9.cif')
    assert output[0] is None 
    assert output[1] is not None  
    assert output[2] is not None  
    assert output[3] == {'display': 'none'}  
   
    invalid_file_data = "data:text/plain;base64,InvalidData"
    output = update_active_link(invalid_file_data, '/', 'test.txt')
    assert output[0] is not None 
    assert output[1] is None  
    assert output[2] is None  
    assert output[3] == {'display': 'none'}  



def test_calculate_interactions():
    interactions = calculate_interactions(sample_pdb_content, 'pdb')
    assert interactions is not None 
    assert 'phosphodiester' in interactions
    assert 'c_base_base' in interactions
    assert 'nc_base_base' in interactions
    assert 'stacking' in interactions

    interactions = calculate_interactions(sample_cif_content, 'cif')
    assert interactions is not None  
    assert 'phosphodiester' in interactions
    assert 'c_base_base' in interactions
    assert 'nc_base_base' in interactions
    assert 'stacking' in interactions

    invalid_data = b"Invalid content"
    interactions = calculate_interactions(invalid_data, 'pdb')
    assert interactions is None  

def test_edge_cases():
    empty_data = b""
    interactions = calculate_interactions(empty_data, 'pdb')
    assert interactions is None 

    large_data = sample_pdb_content * 1000 
    interactions = calculate_interactions(large_data, 'pdb')
    assert interactions is not None 

    no_rna_data = base64.b64encode(b"Protein only content").decode()
    assert check_nucleotide_type(base64.b64decode(no_rna_data), 'pdb') == "Other"