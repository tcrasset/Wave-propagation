import numpy as np
import matplotlib.pyplot as plt
import struct
import pylab

def contourInput(filename):
    # Extracting data from mapped data
    bytes_data = []
    with open(filename, "rb") as binary_file:
        data = binary_file.read()
        a = struct.unpack_from('d',data,0)[0] # double
        b = struct.unpack_from('d',data,8)[0]
        X = struct.unpack_from('I',data,16)[0] # int
        Y = struct.unpack_from('I',data,20)[0]

        print(a)
        print(b)
        print(X)
        print(Y)

        it= struct.iter_unpack("d", data)
        for value in it:
            bytes_data.append(value[0])

    # Ignore first 3 values (24 bytes = 8 + 8 + 4 +4)
    depth = np.array(bytes_data[3:])

    # Reshape and reverse to have the right orientation
    depth = depth.reshape(Y,X)[::-1]

    myContourPlot(depth, 'Contour plot input')


def contourOutput(filename):
    # Extracting data from mapped data
    bytes_data = []
    with open(filename, "rb") as binary_file:
        data = binary_file.read()
        N = struct.unpack_from('I',data,0)[0] # int
        M = struct.unpack_from('I',data,4)[0]

        print(N)
        print(M)

        it= struct.iter_unpack("d", data)
        for value in it:
            bytes_data.append(value[0])

    # Ignore first  value (8 bytes = 4 +4)
    depth = np.array(bytes_data[1:])

    # Reshape and reverse to have the right orientation
    depth = depth.reshape(M,N)[::-1]

    myContourPlot(depth, 'Contour plot output')

def myContourPlot(matrix, title):
    # Plotting
    fig, ax = plt.subplots(figsize=(7,6))

    ax.set_title(title)
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    # ax.xaxis.tick_top()
    # ax.xaxis.set_label_position('top') 
    # plt.gca().invert_yaxis()

    plt.legend()
    # Create contour lines or level curves using matplotlib.pyplot module
    contours = ax.contourf(matrix)
    artists, labels = contours.legend_elements()
    ax.legend(handles=artists, labels=labels,
                loc='upper center', bbox_to_anchor=(0.5, -0.08),
            ncol=3, fancybox=True, shadow=True)
    fig.tight_layout()
    plt.show()



if __name__=='__main__':
    contourInput('test_map.dat')
    contourOutput('eta_test.dat')