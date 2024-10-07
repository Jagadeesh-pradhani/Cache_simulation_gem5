# Open the input and output files
with open('output.out', 'r') as infile, open('extracted_output.out', 'w') as outfile:
    # Loop through each line in the input file
    for line in infile:
        # Check if the line contains the specific string
        if 'system.ruby.l1_cntrl8: MESI_Two_Level-L1cache.sm:' in line:
            # Write the matching line to the output file
            if 'system.ruby.l1_cntrl8: MESI_Two_Level-L1cache.sm:290: NotPresent' not in line:
                outfile.write(line)

print("Lines containing 'system.ruby.l1_cntrl8: MESI_Two_Level-L1cache.sm:' have been extracted and saved.")
