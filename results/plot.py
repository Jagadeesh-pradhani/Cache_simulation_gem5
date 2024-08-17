import matplotlib.pyplot as plt
import os
import argparse

# Function to count the number of stat files in the directory
def count_stats(folder_path, stat_extensions='.txt'):
    files = os.listdir(folder_path)
    stat_count = sum(1 for file in files if file.lower().endswith(stat_extensions))
    return stat_count

# Function to extract "system.ruby.network.average_flit_latency" from the file
def extract_flit_latency(filename, find):
    with open(filename, 'r') as file:
        for line in file:
            if find in line:
                latency_value = float(line.split()[1].strip())
                return latency_value
    return None

# Main function
def main(x_label, find, plot_label):
    src = "stats/"
    values = []
    counts = []
    x_axis = []
    bar_colors = []
    Files = []
    path = []

    filename = os.listdir(src)

    for i in range(count_stats(src)):
        file = os.path.join(src, filename[i])
        path.append(file)
        Files.append(filename[i])

    for i in path:
        value = extract_flit_latency(i, find)
        counts.append(value)

    for i in Files:
        x_axis_name = i.split('.')[0].strip()
        x_axis.append(x_axis_name)

    for i in Files:
        color = 'tab:blue'
        bar_colors.append(color)

    fig, ax = plt.subplots()
    ax.bar(x_axis, counts, color=bar_colors)

    ax.set_ylabel(x_label)
    ax.set_title(plot_label)

    plt.show()

# Argument parser setup
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot from stats files.")
    parser.add_argument("--x_label", type=str, help="Name for the x-axis", default="X-axis")
    parser.add_argument("--Y_value", type=str, help="String to search for in the stats files.")
    parser.add_argument("--plot_label", type=str, help="Name for the plot.", default='Plots')

    args = parser.parse_args()
    main(args.x_label, args.Y_value, args.plot_label)
