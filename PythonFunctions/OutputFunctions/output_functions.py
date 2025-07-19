##############################################################################
##
## Author:      Jamie E. Simon (modified)
##
## Description: 
##
##############################################################################

## Import modulues
import numpy as np

def numpysavetxt(filename, data, headers, decimals=12, delimiter=','):
    """
    Saves numerical data to a text file with specified decimal places and delimiter.

    Parameters:
        filename (str): Name of the output file.
        data (numpy.ndarray): Data to be saved.
        headers (list): List of column headers.
        decimals (int): Number of decimal places (default is 6).
        delimiter (str): Delimiter for separating columns (default is ',').
    """
    # Determine the maximum width required
    all_values = data.flatten()
    max_value_length = max(len(f"{val:.{decimals}e}") for val in all_values)  # Adjust width based on decimals
    max_header_length = max(len(name) for name in headers)  # Find longest header
    max_width = max(max_value_length, max_header_length)  # Use the largest width

    # Define format string dynamically based on user input
    fmt = f'%-{max_width}.{decimals}e'  # Left-aligned with user-defined precision

    # Apply the same format to all columns
    fmt_list = [fmt] * data.shape[1]

    # Generate a properly spaced header
    header_str = ' '.join(f"{name:<{max_width}}" for name in headers)

    # Save the file with aligned columns
    np.savetxt(filename, data, header=header_str, fmt=fmt_list, delimiter=delimiter, comments='')

    print(f"File '{filename}' saved successfully with {decimals} decimal places!")
    
    return