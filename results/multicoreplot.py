import matplotlib.pyplot as plt
import os

# Function to extract the value from the file for a specific metric
def extract_value(filename, find):
    with open(filename, 'r') as file:
        for line in file:
            if find in line:
                value = float(line.split()[1].strip())
                return value
    return None

# Main function to plot DRAM power consumption
def main(x_label, plot_label, convert):
    src = "newprot/"
    metric = "system.ruby.network.average_flit_latency"  # Single metric to extract
    file_labels = ["old.txt", "new.txt"]  # Assume you have old.txt and new.txt
    data = []  # List to store values from each file

    # Extract data for each file
    for file_label in file_labels:
        file_path = os.path.join(src, file_label)
        value = extract_value(file_path, metric)
        if value is not None:
            data.append(value / convert)
        else:
            data.append(0)  # In case the value is not found, use 0 or a default value

    # Create bar plot
    fig, ax = plt.subplots()
    bar_width = 0.35  # Width of the bars
    index = range(len(file_labels))  # Number of bars (one for each file)
    
    # Plot the bars
    bars = ax.bar(index, data, bar_width, color=['tab:blue', 'tab:orange'], tick_label=['dafault', 'Modified'])

    # Set labels, title, and ticks
    ax.set_ylabel(x_label)
    ax.set_title(plot_label)

    # Add exact values on top of each bar, rotated vertically
    for bar in bars:
        height = bar.get_height()
        ax.text(
            bar.get_x() + bar.get_width() / 2,  # X position
            height / 2,  # Y position
            f'{height:.2f}',  # Value (formatted to 2 decimal places)
            ha='center',  # Horizontal alignment
            va='bottom',  # Vertical alignment
            rotation='vertical'  # Rotate text vertically
        )

    plt.tight_layout()
    plt.show()

# Run the main function
if __name__ == "__main__":
    x_label = "Latency"
    plot_label = "Avg flit latency"
    convert = 1  # No conversion in this case

    main(x_label, plot_label, convert)
