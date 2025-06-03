import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

def prep_df_all_for_plot(df_consolidated: pd.DataFrame):

    # df_consolidated has the following columns:
    # Project, Case_ID, Sample_Name, Sample_Type, Raw_Count, Filtered_Count

    #create a violin plot for the counts within "Primary Tumor" and "Blood Derived Normal" for each project_id
    df_plot = df_consolidated[df_consolidated['Sample_Type'].isin(['Primary Tumor', 'Blood Derived Normal'])]

    #if within the same Case_ID there are more than one of the same Sample_Type (multiple files), then take only the one with largest filtered counts
    df_plot = df_plot.loc[df_plot.groupby(['Case_ID', 'Project', 'Sample_Type'])['Filtered_Count'].idxmax()]

    return df_plot

def plot_violin_for_all_TCGA_projects(df_plot, data_type):
    """
    Generates and saves a violin plot comparing counts across 'Project'
    for 'Primary Tumor' and 'Blood Derived Normal' sample types,
    with vertically stacked median labels near the top.

    Args:
        df_plot (pd.DataFrame): DataFrame containing the data with columns
                                 'Project', 'Case_ID', 'Sample_Name', 'Sample_Type',
                                 'Raw_Count', 'Filtered_Count'.
        data_type (str): Column name for the data to plot (e.g., 'Raw_Count' or 'Filtered_Count').
    """

    # Ensure 'Sample_Type' is categorical for consistent ordering and coloring
    df_plot = df_plot.copy()
    df_plot['Sample_Type'] = pd.Categorical(df_plot['Sample_Type'], categories=["Primary Tumor", "Blood Derived Normal"], ordered=True)

    # Ensure the data_type column exists
    if data_type not in df_plot.columns:
        print(f"Error: Data column '{data_type}' not found in DataFrame.")
        return
    # Ensure the data column is numeric and handle potential non-finite values
    df_plot[data_type] = pd.to_numeric(df_plot[data_type], errors='coerce')
    df_plot.dropna(subset=[data_type], inplace=True) # Drop rows where data is NaN
    if not np.isfinite(df_plot[data_type]).all():
         print(f"Warning: Non-finite values found in '{data_type}'. Plot might be affected.")
         # Optionally filter out non-finite values if needed:
         # df_plot = df_plot[np.isfinite(df_plot[data_type])]

    # Set up the plot figure size
    plt.figure(figsize=(24, 12))

    # Define the color palette for violins explicitly
    violin_palette = {"Primary Tumor": "tab:blue", "Blood Derived Normal": "tab:red"}

    # Get the sorted unique projects for consistent plotting order
    if 'Project' not in df_plot.columns:
        print("Error: 'Project' column not found in DataFrame.")
        return
    project_ids_ordered = sorted(df_plot['Project'].dropna().unique())

    # Create the violin plot using the sorted order
    ax = sns.violinplot(
        x='Project',
        y=data_type,
        hue='Sample_Type',
        order=project_ids_ordered, # Ensure seaborn uses the same order
        hue_order=["Primary Tumor", "Blood Derived Normal"],
        data=df_plot,
        split=True,
        inner="quartile",
        palette=violin_palette,
        linewidth=2,
        density_norm='width' # Scale the area of the violins to be equal
    )

    #overlay a stripplot to show the individual data points (also split by Sample_Type)
    # Use a smaller size for the stripplot to avoid cluttering
    sns.stripplot(
        x='Project',
        y=data_type,
        hue='Sample_Type',
        order=project_ids_ordered, # Ensure seaborn uses the same order
        hue_order=["Primary Tumor", "Blood Derived Normal"],
        data=df_plot,
        dodge=True,
        palette=violin_palette,
        size=4, # Adjusted size for better visibility
        alpha=0.5, # Slight transparency for better visibility of overlapping points
        linewidth=1, # Line width for the points
    )

    # --- Median Label Adjustments ---
    # Get the y-axis limits *after* plotting the violins
    y_min, y_max = ax.get_ylim()

    # Set the final y-axis limit *before* calculating label positions
    # This ensures label positions are relative to the final axis height
    final_y_max = y_max * 1.05
    ax.set_ylim(bottom=0, top=final_y_max)

    if data_type == "Raw_Count":
        #add grid lines to the plot, every 1000 units
        ax.yaxis.set_major_locator(plt.MultipleLocator(1000))
        ax.yaxis.set_minor_locator(plt.MultipleLocator(500))
    elif data_type == "Filtered_Count":
        #add grid lines to the plot, every 100 units
        ax.yaxis.set_major_locator(plt.MultipleLocator(100))
        ax.yaxis.set_minor_locator(plt.MultipleLocator(50))
    
    ax.yaxis.grid(which='major', linestyle='--', linewidth=0.5, color='gray')
    ax.yaxis.grid(which='minor', linestyle=':', linewidth=0.5, color='lightgray')
    ax.xaxis.grid(False)  # Disable grid on the x-axis

    # Calculate base position and offset based on the *final* y-axis height
    base_label_y_position = final_y_max * 0.96 # Target 90% of the final height
    vertical_offset = final_y_max * 0.015 # Adjust this multiplier for spacing

    # Calculate and display median values for each sample type
    for sample_type in df_plot['Sample_Type'].cat.categories:
        # Calculate medians and averages for the current sample type, ensuring we group by the ordered projects
        medians = df_plot[df_plot['Sample_Type'] == sample_type].groupby('Project')[data_type].median()
        averages = df_plot[df_plot['Sample_Type'] == sample_type].groupby('Project')[data_type].mean()
        # Reindex medians and averages to match the plot order, filling missing values if necessary
        medians = medians.reindex(project_ids_ordered)
        averages = averages.reindex(project_ids_ordered)

        if sample_type == 'Primary Tumor':
            text_color = 'tab:blue'
            label_y = base_label_y_position + vertical_offset # Place above base (closer to top)
        elif sample_type == 'Blood Derived Normal':
            text_color = 'tab:red' 
            label_y = base_label_y_position - vertical_offset # Place below base
        else:
            text_color = 'black'
            label_y = base_label_y_position # Default position

        # Add text labels for each project using the index from the ordered list
        for i, project in enumerate(project_ids_ordered):
            # Check if median exists and is finite for this project/sample_type
            if project in medians.index and pd.notna(medians[project]) and np.isfinite(medians[project]):
                median_value = medians[project]
                average_value = averages[project] if project in averages.index else np.nan
                ax.text(
                    x=i,                             # X position: use index from the ordered list
                    y=label_y,                       # Y position: adjusted vertically around 90%
                    s=f'A:{average_value:.0f}, M:{median_value:.0f}',         # Text: Formatted median value
                    ha='center',                     # Horizontal alignment
                    va='center',                     # Vertical alignment
                    color=text_color,                # Text color (desaturated)
                    fontsize=16,                     # Font size
                    fontweight='bold'                # Bold font weight
                )
            # else: # Optional: handle cases where a project might miss a sample type or median is NaN/inf
            #    print(f"Notice: No valid median found or plotted for {project} - {sample_type}")


    # --- Plot Formatting ---
    # Rotate x-axis labels for better readability
    plt.xticks(rotation=90, fontsize=16, fontweight='bold')
    plt.yticks(fontsize=16, fontweight='bold')

    # Set plot title and labels
    plot_title = f"MMBIR {data_type.replace('_', ' ')} of Primary Tumor and Blood Derived Normal by Project"
    plt.title(plot_title, fontsize=24, pad=20, fontweight='bold') # Add padding to title
    plt.xlabel("Project ID", fontsize=16, fontweight='bold')
    plt.ylabel(f"Number of MMBIR signatures ({data_type.replace('_', ' ')})", fontsize=16, fontweight='bold')

    # Add a legend, put inside the plot area to the right
    plt.legend(title="Sample Type", loc='upper right', fontsize=16, title_fontsize=18, markerscale=1.5, frameon=True, fancybox=True)
    plt.gca().get_legend().set_bbox_to_anchor((1.0, 0.85)) # Adjust legend position
    

    # Improve layout to prevent labels overlapping figure boundaries
    plt.tight_layout(rect=[0, 0, 0.93, 1]) # Adjust layout rectangle slightly more for legend

    # Save the plot
    file_name = f"violin_plot_{data_type}.png"
    plt.savefig(file_name, dpi=300, bbox_inches='tight')
    plt.close()  # Close the plot to free up memory

    print(f"Violin plot saved as {file_name}")

def print_project_medians(df_plot, data_type):
    """
    Prints the median counts for each Project and Sample_Type in the DataFrame.

    Args:
        df_plot (pd.DataFrame): DataFrame containing the data with columns
        'Project', 'Case_ID', 'Sample_Name', 'Sample_Type', 'Raw_Count', 'Filtered_Count'.
        data_type (str): Type of data to calculate medians for ('raw' or 'filtered').
    """
    print("Medians for each Project and Sample_Type:")
    for project_id in df_plot['Project'].unique():
        for sample_type in df_plot['Sample_Type'].unique():
            median_value = df_plot[(df_plot['Project'] == project_id) & (df_plot['Sample_Type'] == sample_type)][data_type].median()
            print(f"Project ID: {project_id}, Sample Type: {sample_type}, Median: {median_value:.2f}")

df_consolidated = pd.read_csv("all_mmbir_calls_consolidated_0.5.tsv", sep="\t")
df_plot = prep_df_all_for_plot(df_consolidated)


data_type = "Raw_Count" # "Raw_Count" or "Filtered_Count"
print_project_medians(df_plot, data_type)
plot_violin_for_all_TCGA_projects(df_plot=df_plot, data_type=data_type)

data_type = "Filtered_Count" # "Raw_Count" or "Filtered_Count"
print_project_medians(df_plot, data_type)
plot_violin_for_all_TCGA_projects(df_plot=df_plot, data_type=data_type)

