#!/bin/bash

# Define the output file to store results
output_file="rss_usage_results.txt"

# Clear the output file if it exists
> "$output_file"

# Iterate over all ecsbatch.log files
for log_file in ecsbatch.log.*; do
    # Extract job ID from the filename (assumes the format ecsbatch.log.$job-ID)
    job_id=$(echo "$log_file" | awk -F. '{print $3}')
    
    # Run sacct command and get MaxRSS and AveRSS
    rss_usage=$(sacct -j "$job_id" --format="MaxRSS,AveRSS" --noheader | tail -n 1)
    
    # Extract MaxRSS and AveRSS
    max_rss=$(echo "$rss_usage" | awk '{print $1}')
    ave_rss=$(echo "$rss_usage" | awk '{print $2}')
    
    # Normalize MaxRSS and AveRSS to MB (convert K to MB if necessary)
    if [[ "$max_rss" == *K ]]; then
        max_rss_value=$(echo "$max_rss" | sed 's/K//') # Remove K
        max_rss_MB=$(echo "scale=2; $max_rss_value / 1024" | bc) # Divide by 1024 to convert K to MB
        max_rss=$(echo "$max_rss_MB" | awk '{printf "%.2f", $0}') # Ensure two decimal format
    elif [[ "$max_rss" == *M ]]; then
        max_rss=$(echo "$max_rss" | sed 's/M//') # Just remove M
    fi

    if [[ "$ave_rss" == *K ]]; then
        ave_rss_value=$(echo "$ave_rss" | sed 's/K//') # Remove K
        ave_rss_MB=$(echo "scale=2; $ave_rss_value / 1024" | bc) # Divide by 1024 to convert K to MB
        ave_rss=$(echo "$ave_rss_MB" | awk '{printf "%.2f", $0}') # Ensure two decimal format
    elif [[ "$ave_rss" == *M ]]; then
        ave_rss=$(echo "$ave_rss" | sed 's/M//') # Just remove M
    elif [[ "$ave_rss" =~ ^[0-9]+$ ]]; then
        # No suffix, assume the number is in KB and convert to MB
        ave_rss_value=$(echo "$ave_rss")
        ave_rss_MB=$(echo "scale=2; $ave_rss_value / (1024*1024)" | bc) # Divide by 1024^2 to convert byte to MB
        ave_rss=$(echo "$ave_rss_MB" | awk '{printf "%.2f", $0}') # Ensure two decimal format
    fi

    # Append the result to the output file in the format: job_id max_rss ave_rss
    echo "$job_id $max_rss $ave_rss" >> "$output_file"
done

echo "RSS usage values (job_id, max_rss, ave_rss) saved to $output_file"
