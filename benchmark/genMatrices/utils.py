import numpy as np


def SaveMatrix_csv (matrix, output_file_path):

    with open(output_file_path, 'w') as f:

        # write the size
        size = np.shape(matrix)[0]
        f.write(f"{size}\n")
        
        # write an empty line
        f.write("\n")
        
        # write the matrix
        for row in matrix:
            line = "\t".join(map(str, row))
            f.write(f"{line}\n")


def LoadMatrix_csv (input_file_path):

    with open(input_file_path, 'r') as f:
        lines = f.read().splitlines()
    
    # first line contains the size of the matrix
    size = int(lines[0].strip())
    
    # remaining lines contain the matrix (skipping the blank line at index 1)
    matrix_data = [list(map(int, line.split('\t'))) for line in lines[2:2+size]]
    
    # convert to numpy array
    matrix = np.array(matrix_data)

    return matrix
    