import numpy as np

# Function to read the file and store data into a numpy array
def read_file_to_array(file_name):
    with open(file_name, 'r') as file:
        # Read the first line and extract nx, ny, nz
        first_line = file.readline().strip()
        nx, ny, nz = map(int, first_line.split(','))

        # Read the second line and extract the data
        second_line = file.readline().strip()
        data = list(map(int, second_line.split(',')))

        # Convert the data into a numpy array with shape (nx, ny, nz)
        array = np.array(data).reshape((nx, ny, nz))

    return array

# Example usage
file_name = 'data.txt'
array = read_file_to_array(file_name)
print(array)
