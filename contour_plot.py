import numpy as np
import matplotlib.pyplot as plt
import struct
import pylab


# Extracting data from mapped data
bytes_data = []
with open("test_map.dat", "rb") as binary_file:
    data = binary_file.read()
    a = struct.unpack_from('d',data,0)[0]
    b = struct.unpack_from('d',data,8)[0]
    X = struct.unpack_from('q',data,16)[0]
    Y = struct.unpack_from('q',data,24)[0]
    print(a)
    print(b)
    print(X)
    print(Y)

    it= struct.iter_unpack("d", data)
    for value in it:
        bytes_data.append(value[0])

depth = np.array(bytes_data[4:]).reshape(X,Y)
# print(depth)


# Plotting
fig, ax = plt.subplots(figsize=(6,5))
plt.subplots_adjust(top=0.8)

ax.set_title('Contour plot')
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.xaxis.tick_top()
ax.xaxis.set_label_position('top') 
plt.gca().invert_yaxis()


# Create contour lines or level curves using matplotlib.pyplot module
contours = ax.contourf(depth)

# Display the contour plot
plt.show()