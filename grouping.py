import os

def combine_files(i, num_procs=4):
    combined_data = []
    
    for j in range(num_procs):
        filename = f"./sol.{i}_proc_{j}.dat"
        if os.path.exists(filename):
            with open(filename, 'r') as file:
                combined_data.extend(file.readlines())
        else:
            print(f"Warning: File {filename} does not exist and will be skipped.")
    
    output_filename = f"./sol.{i}_.dat"
    with open(output_filename, 'w') as outfile:
        outfile.writelines(combined_data)
    print(f"Combined data for i={i} written to {output_filename}")

def main(end_i, num_procs=4):
    start_i = 0
    for i in range(start_i, end_i + 1):
        combine_files(i, num_procs)

if __name__ == "__main__":
    # Set the end_i value here
    end_i = 200# Adjust as needed// our cas 2

    main(end_i)
