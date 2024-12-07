import plotly.graph_objects as go
import numpy as np
from Bio.PDB import PDBParser, MMCIFParser
from io import StringIO
import rnapolis
from rnapolis import annotator, parser

rna_nucleotides = ['A', 'C', 'G', 'U', 'I']
dna_nucleotides = ['DA', 'DC', 'DG', 'DU', 'DI', 'DT']

colors = []
color_map = {
    'A': 'green',
    'C': 'blue',
    'G': 'red',
    'U': 'yellow',
    'I': 'gray',
    'DA': 'green',
    'DC': 'blue',
    'DG': 'red',
    'DU': 'yellow',
    'DI': 'gray',
    'DT': 'purple'
}

def check_nucleotide_type(decoded_data):
    pdb_content = decoded_data.decode('utf-8')
    for line in pdb_content.splitlines():
        if line.startswith('ATOM') or line.startswith('HETATM'):
            residue_name = line[17:20].strip()  
            if residue_name in rna_nucleotides:
                return "RNA"
            elif residue_name in dna_nucleotides:
                return "DNA"
    return "Other"

filename = '3oxd.pdb'

parserpdb = PDBParser()
structure = parserpdb.get_structure(id='3oxd', file='3oxd.pdb')

points = []
nucleotide_info = []
colors = []
noncanonical_interactions = []
canonical_interactions = []

# Extract nucleotide coordinates and information
for model in structure:
    for chain in model:
        for residue in chain:
            if residue.resname in rna_nucleotides or residue.resname in dna_nucleotides or 'H_' in residue.id[0]:
                coord = []
                for atom in residue:
                    coord.append(atom.get_coord())
                coord_array = np.array(coord)
                center = np.mean(coord_array, axis=0)  
                points.append(center)
                
                color = 'black' if 'H_' in residue.id[0] else color_map.get(residue.resname)
                colors.append(color)

                nucleotide_info.append({
                    "Nucleotide": residue.resname,
                    "Chain_id": chain.id,
                    "Coordinate": center,
                    "Color": color,
                    "Nucleotide_id": residue.id[1]
                })

# Load the 3D structure
with open('3oxd.pdb', 'r') as pdb_file:
    structure1 = parser.read_3d_structure(pdb_file)

# Extract interactions
basepair_interactions = annotator.extract_base_interactions(structure1).basePairs
for pair in basepair_interactions:
    if pair.lw.name != 'cWW':
        noncanonical_interactions.append(pair)
    else:
        canonical_interactions.append(pair)
print(noncanonical_interactions, '\n')
print('\n KANONICZNE')
print(canonical_interactions, '\n')
phosphata_interactions = annotator.extract_base_interactions(structure1).basePhosphateInteractions
ribose_interactions = annotator.extract_base_interactions(structure1).baseRiboseInteractions
stackings_interactions = annotator.extract_base_interactions(structure1).stackings

# Store the lines representing interactions
line_coordinates = []

# Function to handle line creation for different types of interactions
def create_interaction_lines(interaction_list, color, width, dash=None):
    if interaction_list:
        lines = []
        
        for pair in interaction_list:
            nt1 = pair.nt1
            nt2 = pair.nt2

            # Find nucleotide information
            nt1_info = next(nucleotide for nucleotide in nucleotide_info if nucleotide['Nucleotide_id'] == nt1.auth.number and nucleotide['Chain_id'] == nt1.auth.chain and nucleotide['Nucleotide'] == nt1.auth.name)
            nt2_info = next(nucleotide for nucleotide in nucleotide_info if nucleotide['Nucleotide_id'] == nt2.auth.number and nucleotide['Chain_id'] == nt2.auth.chain and nucleotide['Nucleotide'] == nt2.auth.name)

            # Get the coordinates
            nt1_position = nt1_info['Coordinate']
            nt2_position = nt2_info['Coordinate']

            # Create a line for this interaction
            lines.append(go.Scatter3d(
                x=[nt1_position[0], nt2_position[0]],
                y=[nt1_position[1], nt2_position[1]],
                z=[nt1_position[2], nt2_position[2]],
                mode='lines',
                hoverinfo='none',
                line=dict(color=color, width=width, dash=dash),
            ))
        return lines
    return None

# Create the lines for different interactions with different styles
noncanonical_lines = create_interaction_lines(noncanonical_interactions, color='black', width=4, dash='dashdot', )
canonical_lines = create_interaction_lines(canonical_interactions, color='blue', width=6)
phosphate_lines = create_interaction_lines(phosphata_interactions, color='green', width=6, dash='dash')
ribose_lines = create_interaction_lines(ribose_interactions, color='red', width=6, dash='dot')
stacking_lines = create_interaction_lines(stackings_interactions, color='orange', width=6, dash='longdash')

# Prepare the nucleotides for plotting
points_array = np.array(points)

# Create the 3D scatter plot for the nucleotides
scatter = go.Scatter3d(
    x=points_array[:, 0],  
    y=points_array[:, 1],  
    z=points_array[:, 2],  
    mode='markers',
    marker=dict(
        size=6,
        color=colors,
        opacity=1.0,
        line=dict(
            color='white',
            width=0.5
        )
    ),
    hoverinfo='text',
    hovertemplate="Nucleotide: %{customdata[0]} %{customdata[4]}<br>Chain: %{customdata[1]}<br>Coordinates:<br> x: %{customdata[2][0]:.2f}<br> y: %{customdata[2][1]:.2f}<br> z: %{customdata[2][2]:.2f}<extra></extra>",
    customdata=[ 
        [nucleotide['Nucleotide'], nucleotide['Chain_id'], nucleotide['Coordinate'], nucleotide['Color'], nucleotide['Nucleotide_id']] 
        for nucleotide in nucleotide_info
    ],
)

# Initialize figure with nucleotides
fig = go.Figure(data=[scatter])

# Add interaction lines for each type
fig.add_traces(noncanonical_lines)
fig.add_traces(canonical_lines)
fig.add_traces(phosphate_lines)
fig.add_traces(ribose_lines)
fig.add_traces(stacking_lines)

# Update plot layout
fig.update_scenes(xaxis_showspikes=False, yaxis_showspikes=False, zaxis_showspikes=False)
fig.update_layout(
    autosize=True,
    scene=dict(
        xaxis=dict(visible=False),
        yaxis=dict(visible=False),
        zaxis=dict(visible=False)
    ),
    height=610,
    width=1000,
    margin=dict(l=0, r=0, t=0, b=0),
    paper_bgcolor='#fafafb',
    clickmode='event+select',
    dragmode="select",
    newselection_mode="gradual",
)

# Show the figure
fig.show()
