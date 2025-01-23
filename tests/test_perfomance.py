import time
import base64
import pytest
from unittest.mock import MagicMock, patch
from unittest import mock
import multiprocessing
from pages.page2 import update_rna_graph  # Zmienna, którą musisz zaimportować z Twojego projektu
from app import update_active_link

# Test wydajnościowy - Czas ładowania pliku PDB
import time
import base64

def handle_user_request(request_id):
    # Simulate some work being done
    import time
    time.sleep(2)  # This simulates processing time

def test_load_pdb_file_large():
    file_path = "tests/large_file.pdb"  # Example path for a large PDB file

    # Read the file content and encode it in base64 format
    with open(file_path, "rb") as file:
        file_content = file.read()
    encoded_content = base64.b64encode(file_content).decode()
    contents = f"data:file/pdb;base64,{encoded_content}"
    filename = "large_file.pdb"
    pathname = "/"  # Example of the path

    start_time = time.time()
    
    # Call the update_active_link function
    result = update_active_link(contents, pathname, filename)
    
    end_time = time.time()
    load_time = end_time - start_time
    print(f"Czas ładowania pliku large_file.pdb: {load_time:.2f} sekund")
    
    # You can assert based on what the function returns, e.g., no error message
    assert load_time < 5, f"Expected loading time to be less than 5 seconds, but got {load_time:.2f} seconds"
    assert result[0] is None, "Expected no error message for valid file"


# Test obciążeniowy - Wydajność przy dużej liczbie użytkowników
def test_performance_with_multiple_users_10():
    start_time = time.time()
    processes = []
    for i in range(10):
        p = multiprocessing.Process(target=handle_user_request, args=(i+1,))
        processes.append(p)
        p.start()

    for p in processes:
        p.join()

    end_time = time.time()
    total_time = end_time - start_time
    print(f"Czas wykonania dla 10 użytkowników: {total_time:.2f} sekund")
    assert total_time < 30  # Oczekiwany czas poniżej 30 sekund

def test_performance_with_multiple_users_50():
    start_time = time.time()
    processes = []
    for i in range(50):
        p = multiprocessing.Process(target=handle_user_request, args=(i+1,))
        processes.append(p)
        p.start()

    for p in processes:
        p.join()

    end_time = time.time()
    total_time = end_time - start_time
    print(f"Czas wykonania dla 50 użytkowników: {total_time:.2f} sekund")
    assert total_time < 100  # Oczekiwany czas poniżej 100 sekund


# Test obciążeniowy - Skalowanie przy rosnącej liczbie użytkowników
def test_scaling_with_users():
    for num_users in range(10, 101, 10):
        start_time = time.time()
        processes = []

        for i in range(num_users):
            p = multiprocessing.Process(target=handle_user_request, args=(i+1,))
            processes.append(p)
            p.start()

        for p in processes:
            p.join()

        end_time = time.time()
        total_time = end_time - start_time
        print(f"Czas wykonania dla {num_users} użytkowników: {total_time:.2f} sekund")
        assert total_time < 300  # Oczekiwany czas poniżej 300 sekund


# Test obciążeniowy - Wydajność serwera przy dużym obciążeniu
def test_server_performance():
    start_time = time.time()
    processes = []

    for i in range(100):  # Test dla 100 użytkowników
        p = multiprocessing.Process(target=handle_user_request, args=(i+1,))
        processes.append(p)
        p.start()

    for p in processes:
        p.join()

    end_time = time.time()
    total_time = end_time - start_time
    print(f"Czas wykonania przy dużym obciążeniu (100 użytkowników): {total_time:.2f} sekund")
    assert total_time < 300  # Oczekiwany czas poniżej 300 sekund

