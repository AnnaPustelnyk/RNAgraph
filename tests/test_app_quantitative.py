import time
import base64
from app import check_nucleotide_type, calculate_interactions

# Load the sample PDB and CIF file contents for testing
def load_file_content(filename):
    with open(filename, 'rb') as f:
        return f.read()

# Paths to your test files
pdb_file_path = 'tests/sample.pdb'  # Update this path
cif_file_path = 'tests/sample.cif'  # Update this path

# Sample file contents
sample_pdb_content = load_file_content(pdb_file_path)
sample_cif_content = load_file_content(cif_file_path)

def measure_execution_time(function, *args):
    """Utility to measure execution time of a function."""
    start_time = time.time()
    result = function(*args)
    end_time = time.time()
    return result, end_time - start_time

def test_check_nucleotide_type_performance():
    """Test execution time for check_nucleotide_type with varying sizes."""
    small_data = sample_pdb_content
    large_data = sample_pdb_content * 1000  # Simulate a large file
    
    # Small file
    _, time_small = measure_execution_time(check_nucleotide_type, small_data, 'pdb')
    print(f"Execution time for small PDB file: {time_small:.6f} seconds")
    
    # Large file
    _, time_large = measure_execution_time(check_nucleotide_type, large_data, 'pdb')
    print(f"Execution time for large PDB file: {time_large:.6f} seconds")

def test_calculate_interactions_performance():
    """Test execution time for calculate_interactions with varying sizes."""
    small_data = sample_pdb_content
    large_data = sample_pdb_content * 1000  # Simulate a large file
    
    # Small file
    _, time_small = measure_execution_time(calculate_interactions, small_data, 'pdb')
    print(f"Execution time for small PDB file: {time_small:.6f} seconds")
    
    # Large file
    _, time_large = measure_execution_time(calculate_interactions, large_data, 'pdb')
    print(f"Execution time for large PDB file: {time_large:.6f} seconds")

def test_scalability_check_nucleotide_type():
    """Test scalability of check_nucleotide_type with progressively larger files."""
    scales = [1, 10, 100, 500, 1000]  # Multipliers for file size
    for scale in scales:
        scaled_data = sample_pdb_content * scale
        _, exec_time = measure_execution_time(check_nucleotide_type, scaled_data, 'pdb')
        print(f"check_nucleotide_type execution time for scale {scale}: {exec_time:.6f} seconds")

def test_scalability_calculate_interactions():
    """Test scalability of calculate_interactions with progressively larger files."""
    scales = [1, 10, 100, 500, 1000]  # Multipliers for file size
    for scale in scales:
        scaled_data = sample_pdb_content * scale
        _, exec_time = measure_execution_time(calculate_interactions, scaled_data, 'pdb')
        print(f"calculate_interactions execution time for scale {scale}: {exec_time:.6f} seconds")

def test_robustness_large_inputs():
    """Stress test with extremely large inputs."""
    # Simulate a very large file
    very_large_data = sample_pdb_content * 5000
    try:
        interactions, exec_time = measure_execution_time(calculate_interactions, very_large_data, 'pdb')
        print(f"Execution time for very large input: {exec_time:.6f} seconds")
        assert interactions is not None, "Interactions should not be None for valid data"
    except Exception as e:
        print(f"Failed to process very large input: {e}")


