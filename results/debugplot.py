import re
import matplotlib.pyplot as plt

# Initialize lists to store time and occupied values
times = []
occupied_values = []

# Open the file and read line by line
with open('output.out', 'r') as file:
    for line in file:
        # Search for the relevant lines using a regular expression
        match = re.search(r'(\d+): system\.ruby\.l2_cntrl0\.L2cache: Hello World! Occupied value : (\d+)\.', line)
        if match:
            # Extract the time and occupied value from the line
            time = int(match.group(1))
            occupied_value = int(match.group(2))

            # Append the values to the lists
            times.append(time)
            occupied_values.append(occupied_value)

# Plot the data
plt.figure(figsize=(10, 6))
plt.plot(times, occupied_values, marker='o', linestyle='-')
plt.title('Occupied Value over Time for L2 Cache')
plt.xlabel('Time')
plt.ylabel('Occupied Value')
plt.grid(True)
plt.show()
