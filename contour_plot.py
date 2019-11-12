import numpy as np
import matplotlib.pyplot as plt
import struct
import pylab


# Extracting data from mapped data
bytes_data = []
with open("serverFiles/sriLanka.dat", "rb") as binary_file:
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


# Plotting
fig, ax = plt.subplots(figsize=(6,5))
plt.subplots_adjust(top=0.8)

ax.set_title('Contour plot')
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.xaxis.tick_top()
ax.xaxis.set_label_position('top') 
plt.gca().invert_yaxis()

plt.legend()


# Create contour lines or level curves using matplotlib.pyplot module
contours = ax.contourf(depth)

# Display the contour plot
plt.show()