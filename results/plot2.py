import matplotlib.pyplot as plt
import os




# "system.ruby.l1_cntrl0.L1Dcache.m_demand_misses",  # L1D for core 0
# "system.ruby.l1_cntrl0.L1Icache.m_demand_misses",  # L1I for core 0
# "system.ruby.l1_cntrl1.L1Dcache.m_demand_misses",  # L1D for core 1
# "system.ruby.l1_cntrl1.L1Icache.m_demand_misses" 

# Function to count the number of stat files in the directory
def count_stats(folder_path, stat_extensions='.txt'):
    files = os.listdir(folder_path)
    stat_count = sum(1 for file in files if file.lower().endswith(stat_extensions))
    return stat_count

# Function to extract the value from the file for a specific metric
def extract_value(filename, find):
    with open(filename, 'r') as file:
        for line in file:
            if find in line:
                value = float(line.split()[1].strip())
                return value
    return None

# Main function to plot L1 I and D cache hits for both cores
def main(x_label, plot_label, convert):
    src = "newprot/"
    metrics = [
        "system.ruby.l1_cntrl0.L1Dcache.m_demand_hits",  # L1D for core 0
        "system.ruby.l1_cntrl0.L1Icache.m_demand_hits",  # L1I for core 0
        "system.ruby.l1_cntrl1.L1Dcache.m_demand_hits",  # L1D for core 1
        "system.ruby.l1_cntrl1.L1Icache.m_demand_hits"   # L1I for core 1
    ]

    file_labels = ["old.txt", "new.txt"]  # Assume you have old.txt and new.txt
    data = {label: [] for label in file_labels}  # Dictionary to store values from each file
    
    # Extract data for each file and each metric
    for file_label in file_labels:
        file_path = os.path.join(src, file_label)
        for metric in metrics:
            value = extract_value(file_path, metric)
            if value is not None:
                data[file_label].append(value / convert)
            else:
                data[file_label].append(0)

    # Create bar plot
    fig, ax = plt.subplots()
    bar_width = 0.35  # Width of the bars
    index = range(len(metrics))  # Number of bars per group
    
    # Plot for old.txt
    bars_old = ax.bar([i - bar_width/2 for i in index], data["old.txt"], bar_width, label='Default', color='tab:blue')
    
    # Plot for new.txt
    bars_new = ax.bar([i + bar_width/2 for i in index], data["new.txt"], bar_width, label='Modified', color='tab:orange')

    # Set labels, title, and ticks
    ax.set_ylabel(x_label)
    ax.set_title(plot_label)
    ax.set_xticks(index)
    ax.set_xticklabels([
        'L1_D 0 Core 0', 'L1_I 1 Core 0', 
        'L1_D 0 Core 1', 'L1_I 1 Core 1'
    ])
    ax.legend()

    # Add exact values on top of each bar, rotated vertically
    for bars in [bars_old, bars_new]:
        for bar in bars:
            height = bar.get_height()
            ax.text(
                bar.get_x() + bar.get_width() / 2,  # X position
                height/2,  # Y position
                f'{height:.2f}',  # Value (formatted to 2 decimal places)
                ha='center',  # Horizontal alignment
                va='bottom',  # Vertical alignment
                rotation='vertical'  # Rotate text vertically
            )

    plt.tight_layout()
    plt.show()

# Run the main function
if __name__ == "__main__":
    x_label = "Hits"
    plot_label = "L1D and L1I hits"
    convert = 1.0  # No conversion in this case

    main(x_label, plot_label, convert)
