import os

# Function to extract the value from the file for a specific metric
def extract_value(filename, find):
    with open(filename, 'r') as file:
        for line in file:
            if find in line:
                value = float(line.split()[1].strip())
                return value
    return None

# Function to calculate AMAT
def calculate_amat(hits, misses, hit_time, miss_penalty):
    total_accesses = hits + misses
    if total_accesses == 0:
        return float('inf')  # To handle division by zero
    miss_rate = misses / total_accesses
    amat = hit_time + (miss_rate * miss_penalty)
    return amat

# Main function to calculate and compare AMAT for old and new stats
def main(hit_time, miss_penalty):
    src = "newprot/"
    metrics = []
    
    # Generate the metric list for 16 cores
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
                amat = calculate_amat(hits, misses, hit_time, miss_penalty)
                results[file_label].append(amat)
            else:
                results[file_label].append(float('inf'))  # If data is missing, assume worst-case

    # Print comparison results
    print(f"{'Metric':<50} {'Old AMAT':<15} {'New AMAT':<15} {'% Improvement'}")
    for i, metric in enumerate(metrics):
        old_amat = results["old.txt"][i]
        new_amat = results["new.txt"][i]
        if old_amat != float('inf'):
            improvement = ((old_amat - new_amat) / old_amat) * 100
            print(f"{metric:<50} {old_amat:<15.4f} {new_amat:<15.4f} {improvement:.2f}%")
        else:
            print(f"{metric:<50} {'N/A':<15} {'N/A':<15} {'N/A'}")

# Run the main function
if __name__ == "__main__":
    hit_time = 2.0  # Assume hit time is 1 cycle for L1 cache
    miss_penalty = 3.0  # Assume miss penalty is 50 cycles for accessing L2 or memory
    
    main(hit_time, miss_penalty)
