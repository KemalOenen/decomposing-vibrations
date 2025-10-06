import pandas as pd
import plotly.graph_objects as go 


def plot_contribution_matrix(contribution_matrix,normal_coord_harmonic_frequencies,all_internals_string):
    """ 
    Function to plot the contribution matrix using plotly

    Parameters:
    ----------

    contribution_heatmap --> list
        the contribution matrix as a list

    normal_coord_harmonic_frequencies --> list
        list of the harmonic frequencies, these will be the columns of the heatmap

    all_interals_string --> list
        list of all internal coordinates a string     
    """
    
    
    # Build the columns for the matrix (harmonic frequencies)
    columns = {}
    for i,freq in enumerate(normal_coord_harmonic_frequencies):
        columns[i] = freq

    # Prepare Dataframe
    heatmap_df = pd.DataFrame(contribution_matrix)
    heatmap_df = heatmap_df.rename(columns = columns)
    heatmap_df.index = all_internals_string

    # Format values for annotation
    heatmap_df_formatted = heatmap_df.copy()
    heatmap_df_formatted = heatmap_df_formatted.round(1).astype(str)

    # Get user input for customization
    title_input = input("Add Title to Figure? (y/n): ").strip().lower()
    if title_input =='y':
        plot_title = input("Enter the title: ")
    else:
        plot_title = ""

    width_input = input("Enter plot width (press Enter for default 800): ").strip()
    height_input = input("Enter plot height (press Enter for default 800): ").strip()
    fontsize_input = input("Enter annotation font size (press Enter for default 10): ").strip()

    plot_width = int(width_input) if width_input else 800
    plot_height = int(height_input) if height_input else 800
    font_size = int(fontsize_input) if fontsize_input else 10

    fig = go.Figure(data=go.Heatmap(
        z = heatmap_df.values,
        x = [str(element) for element in heatmap_df.columns],
        y = heatmap_df.index,
        colorscale="Blues",
        colorbar = dict(title="Contribution %"),
        zmin=0,
        zmax=100, 
        text=heatmap_df_formatted.values,
        texttemplate="%{text}",
        textfont={"size":font_size}
    )
    )

    fig.update_layout(
        title=plot_title if plot_title else None,
        yaxis={'autorange':"reversed"},
        width=plot_width,
        height=plot_height       
    )


    fig.write_html("heatmap_contribution_table.html")

def plot_sankey_diagram(ContributionTable,bonds,angles,linear_angles,dihedrals,out_of_plane):
    """ 
    Function to generate a sankey diagramm of the nomodeco contributions
    """

    # User Input --> Minimal contribution for plotting 
    min_contribution = input("Enter minimum number of contribution for plotting (integer): ")
    min_contribution = int(min_contribution)

    def get_coordinate_type(coord):
        """ 
        Helper Function to detemrine the type of coordinate
        """
        # Check if coord is tuple
        if not isinstance(coord, tuple):
            coord_tuple = tuple(coord.strip("()").split(","))
            coord_tuple= tuple(element.strip() for element in coord_tuple)

        if coord_tuple in bonds:
            return "bond"
        elif coord_tuple in angles:
            return "angle"
        elif coord_tuple in linear_angles:
            return "linear_angle"
        elif coord_tuple in dihedrals:
            return "dihedral"
        elif coord_tuple in out_of_plane:
            return "out-of-plane"


            
        
        
    # Prepare for Sankey Diagramm
    labels = []
    source = []
    target = []
    value = []
    link_labels = []
    link_colors = []
    coord_types = []

    # Color mapping for coordinate types:

    color_map = {
        "bond": "rgba(255, 99, 132, 0.8)",  # Red
        "angle": "rgba(54, 162, 235, 0.8)",  # Blue
        "dihedral": "rgba(75, 192, 192, 0.8)",  # Green
        "out-of-plane": "rgba(255, 206, 86, 0.8)",  # Yellow
        "linear_angle": "rgba(153, 102, 255, 0.8)",  # Purple
    }

    

    # Create Nodes
    int_coords_label = ContributionTable["Internal Coordinate"].tolist()
    modes = [str(col) for col in ContributionTable.columns if col not in ["Internal Coordinate", "Intrinsic Frequencies"]]
    labels = int_coords_label + modes

    # Create Mappings from names to indices
    coord_indices = {coord: idx for idx,coord in enumerate(int_coords_label)}
    mode_indices = {mode: idx + len(int_coords_label) for idx,mode in enumerate(modes)}

    # Build up the links
    for _, row in ContributionTable.iterrows():
        coord = row["Internal Coordinate"]
        
        # Get coordinate type
        coord_type = get_coordinate_type(coord)
        coord_types.append(coord_type)
        for mode in modes:
            contribution = float(row[mode])
            if contribution > min_contribution:
                source.append(coord_indices[coord])
                target.append(mode_indices[mode])
                value.append(contribution)
                link_labels.append(f"{coord} to {mode}: {contribution:.1f}%")
                link_colors.append(color_map[coord_type]) # Assign color based on type

    legend_trace = []
    for coord_type, color in color_map.items():
        legend_trace.append(
            go.Scatter(
                x=[None],y=[None],
                mode= "markers",
                marker=dict(color=color, size=10),
                name=coord_type,
                hoverinfo="none"
            )
        )
    
    # Create Sankey Diagram
    fig = go.Figure(
        data=[go.Sankey(
            node=dict(
                pad=15,
                thickness=20,
                line=dict(color="black",width=0.5),
                label = labels,
                color="blue"
            ),
            link=dict(
                source=source,
                target=target,
                value=value,
                label=link_labels,
                color=link_colors,
                hovertemplate="%{label}<extra></extra>",
                
            )
        ),
            *legend_trace
        ]
    )

    # Change general layout
    fig.update_layout(
        title_text = "Internal Coordinate Contributions",
        font_size = 10,
        height = 800,
        showlegend=True,
        legend=dict(
            orientation="h",
            yanchor="bottom",
            y=1.02,
            xanchor="right",
            x=1
        )
    )
    fig.write_html("sankey_diagram.html")
    
    


