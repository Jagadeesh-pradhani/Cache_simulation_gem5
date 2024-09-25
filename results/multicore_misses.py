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

# Main function to plot L1 I and D cache misses for all 16 cores and print comparison results
def main(x_label, plot_label, convert):
    src = "newprot/"
    metrics = []
    # Generate the metric list for 16 cores
    for core_id in range(16):
        metrics.append(f"system.ruby.l1_cntrl{core_id}.L1Dcache.m_demand_misses")  # L1D misses for each core
        metrics.append(f"system.ruby.l1_cntrl{core_id}.L1Icache.m_demand_misses")  # L1I misses for each core

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
    fig, ax = plt.subplots(figsize=(15, 8))  # Adjusted size for better visibility
    bar_width = 0.35  # Width of the bars
    index = range(len(metrics))  # Number of bars per group
    
    # Plot for old.txt
    bars_old = ax.bar([i - bar_width/2 for i in index], data["old.txt"], bar_width, label='Default', color='tab:blue')
    
    # Plot for new.txt
    bars_new = ax.bar([i + bar_width/2 for i in index], data["new.txt"], bar_width, label='Modified', color='tab:orange')

    # Set labels, title, and ticks
    ax.set_ylabel(x_label)
    ax.set_title(plot_label)
    ax.set_xticks(range(0, len(metrics), 2))  # Show ticks every 2 bars (one for each core's L1D and L1I)
    ax.set_xticklabels([f'Core {i // 2}' for i in range(0, len(metrics), 2)], rotation=90)

    ax.legend()

    # Add exact values on top of each bar, rotated vertically
    for bars in [bars_old, bars_new]:
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

    # Print comparison results where new misses are lesser than old
    print("\nComparison of 'new' vs 'old' values (Lesser in new):")
    for i, metric in enumerate(metrics):
        old_value = data["old.txt"][i]
        new_value = data["new.txt"][i]
        if new_value < old_value:
            percentage_decrease = ((old_value - new_value) / old_value) * 100 if old_value != 0 else 0
            print(f"{metric}: New value is lesser by {percentage_decrease:.2f}% (Old: {old_value}, New: {new_value})")

# Run the main function
if __name__ == "__main__":
    x_label = "Misses"
    plot_label = "L1D and L1I misses for 16 cores"
    convert = 1.0  # No conversion in this case

    main(x_label, plot_label, convert)
