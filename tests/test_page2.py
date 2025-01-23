import base64
import unittest
from unittest.mock import patch, MagicMock
import numpy as np
import plotly.graph_objects as go
from dash import Dash
import time
from Bio.PDB import PDBParser
from io import StringIO
import os
import tracemalloc
from pages.page2 import update_rna_graph, display_selected_info, clear_selection, show_heteroatoms, update_interaction_info, layout, rna_nucleotides, dna_nucleotides, color_map

class TestPage(unittest.TestCase):

    @patch('pages.page2.PDBParser')
    @patch('pages.page2.MMCIFParser')
    @patch('pages.page2.colors', new=['red'])  # Mock the colors list to avoid dependency issues
    def test_update_rna_graph(self, mock_mmcif_parser, mock_pdb_parser):
        # Simulated PDB content (minimal valid structure)
        mock_pdb_content = """
    ATOM      1  P     G A   1      11.104  13.207   2.814  1.00 20.00           P
    ATOM      2  O1P   G A   1      10.293  12.151   2.153  1.00 20.00           O
    ATOM      3  O2P   G A   1      12.004  12.378   3.590  1.00 20.00           O
    TER
    END
    """
        encoded_pdb = base64.b64encode(mock_pdb_content.encode('utf-8')).decode()
        mock_data = {
            'url': f'data:application/octet-stream;base64,{encoded_pdb}',
            'ext': 'pdb'
        }
        mock_filename = 'test.pdb'

        # Mock parser return value
        mock_structure = MagicMock()
        mock_model = MagicMock()
        mock_chain = MagicMock()
        mock_residue = MagicMock()

        mock_atom = MagicMock()
        mock_atom.get_coord.return_value = [11.0, 13.2, 2.8]  # Simulate an atom's coordinates

        mock_residue.resname = 'G'  # A valid nucleotide
        mock_residue.__iter__.return_value = [mock_atom]
        mock_residue.id = (' ', 1, ' ')

        mock_chain.__iter__.return_value = [mock_residue]
        mock_chain.id = 'A'

        mock_model.__iter__.return_value = [mock_chain]
        mock_structure.__iter__.return_value = [mock_model]

        mock_pdb_parser.return_value.get_structure.return_value = mock_structure

        # Call the function
        figure, rna_graph_style, side_bar_style, nucleotide_info_style, options = update_rna_graph(mock_data, mock_filename)

        # Assertions
        self.assertEqual(rna_graph_style, {'display': 'block'})
        self.assertEqual(side_bar_style, {'display': 'block'})
        self.assertEqual(nucleotide_info_style, {'display': 'block'})

        # Validate figure content
        self.assertEqual(len(figure.data), 1)  # One trace for nucleotides
        nucleotide_trace = figure.data[0]
        self.assertEqual(nucleotide_trace.name, 'nucleotides')
        self.assertEqual(nucleotide_trace.marker.color, ('rgb(228, 34, 23)',))  # Mocked color list
        self.assertEqual(nucleotide_trace.mode, 'markers')

        # Validate heteroatom options
        self.assertEqual(options, [{'label': 'Show heteroatoms', 'value': 'heteroatoms', 'disabled': True}])


    def test_display_selected_info(self):
        click_data = {
            'points': [{
                'customdata': ['A', 'A', [1.0, 2.0, 3.0], 'rgb(225, 246, 0)', 1]
            }]
        }
        current_figure = go.Figure(data=[go.Scatter3d(customdata=[['A', 'A', [1.0, 2.0, 3.0], 'rgb(225, 246, 0)', 1]])])
        
        nucleotide_value, chain_value, updated_figure = display_selected_info(click_data, current_figure, None)

        self.assertEqual(nucleotide_value, 'A1')
        self.assertEqual(chain_value, 'A')
        self.assertIsInstance(updated_figure, go.Figure)

    def test_clear_selection(self):
        # Create a mock current figure
        current_figure = go.Figure(data=[
            go.Scatter3d(marker=dict(color=['red'], size=[5])),
            go.Scatter3d(name='heteroatoms', marker=dict(color=['black'], size=[5]))
        ])
        n_clicks = 1  # Simulate the button being clicked
        relayoutData = None  # No relayout data for this test

        # Call the clear_selection function
        updated_figure, nucleotide_input, chain_input, coord_x, coord_y, coord_z, coordinates_style, toggle_button = clear_selection(n_clicks, current_figure, relayoutData)

        # Assertions
        self.assertTrue(updated_figure)  # Check if the color is set to 'red'
        self.assertEqual(nucleotide_input, "")  # Check if nucleotide input is reset
        self.assertEqual(chain_input, "")  # Check if chain input is reset
        self.assertEqual(coord_x, "x: ")  # Check if coord_x is reset
        self.assertEqual(coord_y, "y: ")  # Check if coord_y is reset
        self.assertEqual(coord_z, "z: ")  # Check if coord_z is reset
        self.assertEqual(coordinates_style, {'display': 'none'})  # Check if coordinates container is hidden
        self.assertEqual(toggle_button, '▼')  # Check if toggle button is set correctly

    def test_show_heteroatoms(self):
        current_figure = go.Figure(data=[
            go.Scatter3d(marker=dict(color=['red'], size=[5])),
            go.Scatter3d(name='heteroatoms', marker=dict(color=['black'], size=[5]))
        ])
        values = ['heteroatoms']
        
        updated_figure, updated_values = show_heteroatoms(values, current_figure, None)

        self.assertTrue(updated_figure)
        self.assertEqual(updated_values, ['heteroatoms'])

class TestIntegration(unittest.TestCase):

    @patch('pages.page2.create_interaction_lines')
    def test_update_interaction_info(self, mock_create_interaction_lines):
        selected_interactions = ['phosphodiester']
        data = {
            'phosphodiester': [
                {
                    'nt1': {'auth': {'name': 'A', 'number': 1, 'chain': 'A'}},
                    'nt2': {'auth': {'name': 'C', 'number': 2, 'chain': 'A'}}
                }
            ]
        }
        nucleotide_info = [
            {"Nucleotide": "A", "Chain_id": "A", "Coordinate": [1.0, 2.0, 3.0], "Color": "red", "Nucleotide_id": 1},
            {"Nucleotide": "C", "Chain_id": "A", "Coordinate": [4.0, 5.0, 6.0], "Color": "blue", "Nucleotide_id": 2}
        ]
        current_figure = {
            'data': [{
                'name': 'nucleotides',
                'customdata': [
                    [nucleotide['Nucleotide'], nucleotide['Chain_id'], nucleotide['Coordinate'], nucleotide['Color'], nucleotide['Nucleotide_id']] 
                    for nucleotide in nucleotide_info
                ]
            }]
        }
        interactions = {
            'phosphodiester': [
                {
                    'nt1': {'auth': {'name': 'A', 'number': 1, 'chain': 'A'}},
                    'nt2': {'auth': {'name': 'C', 'number': 2, 'chain': 'A'}}
                }
            ],
            'c_base_base': [],
            'nc_base_base': [],
            'stacking': []
        }

        mock_create_interaction_lines.return_value = [
            go.Scatter3d(
                x=[1, 2], 
                y=[3, 4], 
                z=[5, 6], 
                mode='lines', 
                line=dict(color='green')
            )
        ]

        updated_figure, interaction_options = update_interaction_info(selected_interactions, data, current_figure, interactions, None)

        # Assertions
        self.assertIsNotNone(updated_figure)  # Check that the updated figure is not None
        self.assertIsInstance(interaction_options, list)  # Check that interaction_options is a list

        # Check that the interaction options are generated correctly
        self.assertIn({'label': 'Phosphodiester interactions', 'value': 'phosphodiester', 'disabled': False}, interaction_options)
        self.assertIn({'label': 'Canonical interactions', 'value': 'c_base_base', 'disabled': True}, interaction_options)
        self.assertIn({'label': 'Non-canonical interactions', 'value': 'nc_base_base', 'disabled': True}, interaction_options)
        self.assertIn({'label': 'Stacking interactions', 'value': 'stacking', 'disabled': True}, interaction_options)

class TestFunctional(unittest.TestCase):

    def setUp(self):
        self.app = Dash(__name__)
        self.app.layout = layout
        self.client = self.app.server.test_client()

    def test_layout(self):
        response = self.client.get('/')
        self.assertEqual(response.status_code, 200)

class TestPerformance(unittest.TestCase):

    def test_performance_update_rna_graph(self):
        pdb_file_path = os.path.join('tests', 'sample.pdb')
        with open(pdb_file_path, 'r') as pdb_file:
            example_pdb_content = pdb_file.read()
        
        pdb_data = example_pdb_content.encode('utf-8')  # Convert string to bytes
        pdb_data_file = base64.b64encode(pdb_data).decode()

        mock_filename = 'sample.pdb'
        mock_data = {
            'url': f'data:data:application/octet-stream;base64;base64,{pdb_data_file}',
            'ext': 'pdb'  # Include the file extension
        }

        start_time = time.time()
        update_rna_graph(mock_data, mock_filename)
        end_time = time.time()

        self.assertLess(end_time - start_time, 1)

    def test_performance_display_selected_info(self):
        click_data = {
            'points': [{
                'customdata': ['A', 'A', [1.0, 2.0, 3.0], 'rgb(225, 246, 0)', 1]
            }]
        }
        current_figure = go.Figure(data=[go.Scatter3d(
            customdata=[['A', 'A', [1.0, 2.0, 3.0], 'rgb(225, 246, 0)', 1]]
        )])
        
        start_time = time.time()
        display_selected_info(click_data, current_figure, None)
        end_time = time.time()

        self.assertLess(end_time - start_time, 0.5)  # Assert under 0.5s

    def test_memory_update_rna_graph(self):
        mock_pdb_content = """
        ATOM      1  P     G A   1      11.104  13.207   2.814  1.00 20.00           P
        ATOM      2  O1P   G A   1      10.293  12.151   2.153  1.00 20.00           O
        TER
        END
        """
        encoded_pdb = base64.b64encode(mock_pdb_content.encode('utf-8')).decode()
        mock_data = {'url': f'data:application/octet-stream;base64,{encoded_pdb}', 'ext': 'pdb'}
        mock_filename = 'test.pdb'

        tracemalloc.start()
        update_rna_graph(mock_data, mock_filename)
        current, peak = tracemalloc.get_traced_memory()
        tracemalloc.stop()

        self.assertLess(peak, 10 * 1024 * 1024)  # Assert peak memory < 10MB

    @patch('pages.page2.PDBParser')
    def test_large_pdb_file(self, mock_pdb_parser):
        large_pdb_content = "\n".join([
            f"ATOM  {i:5d}  P     G A   {i}      {i*1.5:.3f}  {i*2.0:.3f}  {i*2.5:.3f}  1.00 20.00           P"
            for i in range(1, 10000)
        ]) + "\nTER\nEND"

        encoded_pdb = base64.b64encode(large_pdb_content.encode('utf-8')).decode()
        mock_data = {
            'url': f'data:application/octet-stream;base64,{encoded_pdb}',
            'ext': 'pdb'
        }
        mock_filename = 'test.pdb'

        # Mock parser return value
        mock_structure = MagicMock()
        mock_model = MagicMock()
        mock_chain = MagicMock()
        mock_residue = MagicMock()

        mock_atom = MagicMock()
        mock_atom.get_coord.return_value = [11.0, 13.2, 2.8]  # Simulate an atom's coordinates

        mock_residue.resname = 'G'  # A valid nucleotide
        mock_residue.__iter__.return_value = [mock_atom]
        mock_residue.id = (' ', 1, ' ')

        mock_chain.__iter__.return_value = [mock_residue]
        mock_chain.id = 'A'

        mock_model.__iter__.return_value = [mock_chain]
        mock_structure.__iter__.return_value = [mock_model]

        mock_pdb_parser.return_value.get_structure.return_value = mock_structure

        start_time = time.time()
        figure, _, _, _, _ = update_rna_graph(mock_data, mock_filename)
        end_time = time.time()
        render = end_time - start_time

        print(f"Czas ładowania pliku large_file.pdb: {render:.2f} sekund")

        self.assertLess(end_time - start_time, 3)  # Ensure runtime < 3s
        self.assertEqual(len(figure.data), 1)  # Ensure single trace for nucleotides

class TestPerformance1(unittest.TestCase):
    def test_performance_update_rna_graph(self):
        pdb_file_path = os.path.join('tests', 'sample.pdb')
        with open(pdb_file_path, 'r') as pdb_file:
            example_pdb_content = pdb_file.read()

        # Prepare mock data
        pdb_data = example_pdb_content.encode('utf-8')  # Convert string to bytes
        pdb_data_file = base64.b64encode(pdb_data).decode()
        mock_filename = 'sample.pdb'
        mock_data = {
            'url': f'data:application/octet-stream;base64,{pdb_data_file}',
            'ext': 'pdb'  # Include the file extension
        }

        execution_times = []
        memory_usages = []

        for _ in range(10):  # Run multiple iterations to collect data
            tracemalloc.start()  # Start tracking memory usage
            start_time = time.time()
            update_rna_graph(mock_data, mock_filename)
            end_time = time.time()
            current_memory, peak_memory = tracemalloc.get_traced_memory()
            tracemalloc.stop()

            execution_times.append(end_time - start_time)
            memory_usages.append(peak_memory / 1024)  # Convert to KB

        # Save data for graphing
        np.savetxt("update_rna_graph_times.csv", execution_times, delimiter=",", header="Execution Time (s)")
        np.savetxt("update_rna_graph_memory.csv", memory_usages, delimiter=",", header="Memory Usage (KB)")

        # Assertions to ensure performance is reasonable
        self.assertLess(max(execution_times), 1, "Execution time exceeded 1 second")
        self.assertLess(max(memory_usages), 1024, "Memory usage exceeded 1 MB")

    def test_performance_display_selected_info(self):
        click_data = {
            'points': [{
                'customdata': ['A', 'A', [1.0, 2.0, 3.0], 'rgb(225, 246, 0)', 1]
            }]
        }
        current_figure = go.Figure(data=[go.Scatter3d(
            customdata=[['A', 'A', [1.0, 2.0, 3.0], 'rgb(225, 246, 0)', 1]]
        )])

        execution_times = []

        for _ in range(10):  # Run multiple iterations to collect data
            start_time = time.time()
            display_selected_info(click_data, current_figure, None)
            end_time = time.time()

            execution_times.append(end_time - start_time)

        # Save data for graphing
        np.savetxt("display_selected_info_times.csv", execution_times, delimiter=",", header="Execution Time (s)")

        # Assertions to ensure performance is reasonable
        self.assertLess(max(execution_times), 0.5, "Execution time exceeded 0.5 seconds")

    def test_performance_memory_usage(self):
        tracemalloc.start()

        pdb_file_path = os.path.join('tests', 'sample.pdb')
        with open(pdb_file_path, 'r') as pdb_file:
            example_pdb_content = pdb_file.read()

        pdb_data = example_pdb_content.encode('utf-8')  # Convert string to bytes
        pdb_data_file = base64.b64encode(pdb_data).decode()
        mock_filename = 'sample.pdb'
        mock_data = {
            'url': f'data:application/octet-stream;base64,{pdb_data_file}',
            'ext': 'pdb'  # Include the file extension
        }

        update_rna_graph(mock_data, mock_filename)
        current_memory, peak_memory = tracemalloc.get_traced_memory()
        tracemalloc.stop()

        # Save peak memory usage for graphing
        with open("memory_usage.csv", "w") as f:
            f.write("Peak Memory Usage (KB)\n")
            f.write(f"{peak_memory / 1024}\n")  # Convert bytes to KB

        # Assert memory usage does not exceed a threshold
        self.assertLess(peak_memory / 1024, 1500, "Memory usage exceeded 1.5 MB")


if __name__ == '__main__':
    unittest.main()