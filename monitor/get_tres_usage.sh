#!/bin/bash

# Define the output file to store results
output_file="tres_usage_results.txt"

# Clear the output file if it exists
> "$output_file"

# Iterate over all ecsbatch.log files
for log_file in ecsbatch.log.*; do
    # Extract job ID from the filename (assumes the format ecsbatch.log.$job-ID)
    job_id=$(echo "$log_file" | awk -F. '{print $3}')
    
    # Run sacct command and get TRESUsageInMax
    tres_usage=$(sacct -j "$job_id" --format="TRESUsageInMax%80" --noheader | tail -n 1)
    
    # Extract only the values for disk, mem, pages, and vmem
    disk=$(echo "$tres_usage" | awk -F'fs/disk=' '{print $2}' | awk -F',' '{print $1}')
    
    # Extract mem and normalize (convert M to K if necessary, ensure integer)
    mem=$(echo "$tres_usage" | awk -F'mem=' '{print $2}' | awk -F',' '{print $1}' | xargs) # Trim whitespace
    if [[ "$mem" == *M ]]; then
        mem_value=$(echo "$mem" | sed 's/M//') # Remove M
        mem_K=$(echo "$mem_value * 1024" | bc) # Multiply by 1024 to convert M to K
        mem=$(echo "$mem_K" | awk '{printf "%d", $0}') # Ensure integer format
    else
        mem_value=$(echo "$mem" | sed 's/K//') # Just remove K
        mem=$(echo "$mem_value" | awk '{printf "%d", $0}') # Ensure integer format
    fi
    
    # Extract pages    
    pages=$(echo "$tres_usage" | awk -F'pages=' '{print $2}' | awk -F',' '{print $1}' | xargs) # Trim whitespace
    
    # Extract vmem and normalize (convert M to K if necessary, ensure integer)
    vmem=$(echo "$tres_usage" | awk -F'vmem=' '{print $2}' | awk -F',' '{print $1}' | xargs) # Trim whitespace
    if [[ "$vmem" == *M ]]; then
        vmem_value=$(echo "$vmem" | sed 's/M//') # Remove M
        vmem_K=$(echo "$vmem_value * 1024" | bc) # Multiply by 1024 to convert M to K
        vmem=$(echo "$vmem_K" | awk '{printf "%d", $0}') # Ensure integer format
    else
        vmem_value=$(echo "$vmem" | sed 's/K//') # Just remove K
        vmem=$(echo "$vmem_value" | awk '{printf "%d", $0}') # Ensure integer format
    fi

    # Append the result to the output file in the format: job_id disk mem pages vmem
    echo "$job_id $disk $mem $pages $vmem" >> "$output_file"
done

echo "TRES usage values (job_id, disk, mem, pages, vmem) saved to $output_file"
