import os
import matplotlib.pyplot as plt

# Function to extract the value from the file for a specific metric
def extract_value(filename, find):
    with open(filename, 'r') as file:
        for line in file:
            if find in line:
                value = float(line.split()[1].strip())
                return value
    return None

# Function to calculate AMAT for L1
def calculate_amat(hits, misses, l1_request_latency, l1_response_latency, l1_to_l2_latency, l2_request_latency, l2_response_latency):
    total_accesses = hits + misses
    if total_accesses == 0:
        return float('inf')  # To handle division by zero
    miss_rate = misses / total_accesses
    # AMAT calculation considering latencies
    amat = (
        (hits * (l1_request_latency + l1_response_latency)) + 
        (misses * (l1_request_latency + l1_response_latency + l1_to_l2_latency + l2_request_latency + l2_response_latency))
    ) / total_accesses
    return amat

# Main function to calculate and compare AMAT for L1 caches
def main(latencies):
    src = "newprot/"
    metrics = []

    # Generate the metric list for 16 cores (L1D and L1I for each core)
    for core_id in range(16):
        metrics.append(f"system.ruby.l1_cntrl{core_id}.L1Dcache")  # L1D for each core
        metrics.append(f"system.ruby.l1_cntrl{core_id}.L1Icache")  # L1I for each core

    file_labels = ["old.txt", "new.txt"]  # Assume you have old.txt and new.txt
    results = {label: [] for label in file_labels}  # Dictionary to store AMAT results

    # Extract hits and misses for each metric and calculate AMAT
    for file_label in file_labels:
        file_path = os.path.join(src, file_label)
        for metric in metrics:
            hits = extract_value(file_path, f"{metric}.m_demand_hits")
            misses = extract_value(file_path, f"{metric}.m_demand_misses")
            if hits is not None and misses is not None:
                # Calculate AMAT for L1 cache
                amat_l1 = calculate_amat(hits, misses, *latencies)
                results[file_label].append(amat_l1)  # Store AMAT value
            else:
                results[file_label].append(float('inf'))  # If data is missing

    # Print comparison results
    print(f"{'Metric':<50} {'Old AMAT L1':<15} {'New AMAT L1':<15} {'% Improvement L1'}")

    # Prepare data for plotting
    labels = []
    old_amat_values = []
    new_amat_values = []

    for i in range(0, len(metrics), 2):  # Process pairs of L1D and L1I for each core
        metric_d = metrics[i]    # L1D
        metric_i = metrics[i + 1]  # L1I

        old_amat_l1_d = results["old.txt"][i]
        new_amat_l1_d = results["new.txt"][i]
        old_amat_l1_i = results["old.txt"][i + 1]
        new_amat_l1_i = results["new.txt"][i + 1]

        # Calculate percentage improvements for L1D
        improvement_l1_d = ((old_amat_l1_d - new_amat_l1_d) / old_amat_l1_d * 100) if old_amat_l1_d != float('inf') else float('nan')
        # Calculate percentage improvements for L1I
        improvement_l1_i = ((old_amat_l1_i - new_amat_l1_i) / old_amat_l1_i * 100) if old_amat_l1_i != float('inf') else float('nan')

        print(f"{metric_d:<50} {old_amat_l1_d:<15.2f} {new_amat_l1_d:<15.2f} {improvement_l1_d:.2f}%")
        print(f"{metric_i:<50} {old_amat_l1_i:<15.2f} {new_amat_l1_i:<15.2f} {improvement_l1_i:.2f}%")

        # Prepare labels and values for plotting
        labels.append(f'L1D {i//2} Core')
        labels.append(f'L1I {i//2} Core')
        old_amat_values.append(old_amat_l1_d)
        old_amat_values.append(old_amat_l1_i)
        new_amat_values.append(new_amat_l1_d)
        new_amat_values.append(new_amat_l1_i)

    # Plotting AMAT values
    bar_width = 0.35
    index = range(len(labels))

    fig, ax = plt.subplots(figsize=(15, 8))
    bars_old = ax.bar([i - bar_width/2 for i in index], old_amat_values, bar_width, label='Old AMAT', color='tab:blue')
    bars_new = ax.bar([i + bar_width/2 for i in index], new_amat_values, bar_width, label='New AMAT', color='tab:orange')

    # Set labels, title, and ticks
    ax.set_ylabel('AMAT (Cycles)')
    ax.set_title('AMAT Comparison for L1 Caches')
    ax.set_xticks(index)
    ax.set_xticklabels(labels,rotation=90)
    ax.legend()

    # Add exact values on top of each bar
    for bars in [bars_old, bars_new]:
        for bar in bars:
            height = bar.get_height()
            ax.text(
                bar.get_x() + bar.get_width() / 2,  # X position
                height/2,  # Y position (slightly above the bar)
                f'{height:.2f}',  # Value (formatted to 2 decimal places)
                ha='center',  # Horizontal alignment
                va='bottom',   # Vertical alignment
                rotation='vertical' 
            )

    plt.tight_layout()
    plt.show()

# Run the main function
if __name__ == "__main__":
    # Define the latency parameters
    latencies = [
        2,  # l1_request_latency
        2,  # l1_response_latency
        1,  # l1_to_l2_latency
        2,  # l2_request_latency
        2   # l2_response_latency
    ]
    
    main(latencies)
