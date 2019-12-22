import numpy as np
import matplotlib.pyplot as plt
import struct
import pylab
import os

def contourInput(filename, output_filename):
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

	myContourPlot(depth, 'Contour plot input', output_filename)


def contourOutput(filename, output_filename):
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

	myContourPlot(depth, 'Contour plot output', output_filename)

def myContourPlot(matrix, title, output_filename):
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
	# ax.legend(handles=artists, labels=labels,
	#             loc='upper center', bbox_to_anchor=(0.5, -0.08),
	#         ncol=3, fancybox=True, shadow=True)
	plt.colorbar(contours, ax=ax)
	fig.tight_layout()
	fig.savefig(output_filename)
	plt.show()


if __name__=='__main__':
	# contourInput('test_map.dat')
	directory = "/home/tom/Documents/Uliege/Master2/HPC/Project2/Results"

	# foldername = "matrices_of_mod_space_250_4_4" # Exploded
	# filename = "eta_40000"	

	# foldername = "matrices_of_mod_space_500_4_4" # Exploded
	# filename = "eta_40000"

	# foldername = "matrices_of_mod_space_1000_4_4"
	# filename = "eta_40000"

	# foldername = "matrices_of_mod_space_2000_4_4"
	# filename = "eta_40000"
	
	# foldername = "matrices_of_mod_space_10000_4_4"
	# filename = "eta_40000"

	# foldername = "matrices_of_mod_space_25000_4_4"
	# filename = "eta_40000"

	# foldername = "matrices_of_mod_time_10_4_4" # Exploded
	# filename = "eta_1000"

	# foldername = "matrices_of_mod_time_1_4_4" # Exploded
	# filename = "eta_100000"

	# foldername = "matrices_of_mod_time_0,1_4_4"
	# filename = "eta_1000"

	# foldername = "matrices_of_mod_time_0,01_4_4"
	# filename = "eta_200000"

	# foldername = "matrices_of_mod_time_0,05_4_4"
	# filename = "eta_40000"

	# foldername = "matrices_of_mod_time_0,001_4_4"
	# filename = "eta_1500000"


	foldername = "matrices_of_sriLanka_implicit_25000_2_1"
	filename = "eta_100"


	dat = ".dat"
	svg = ".svg"
	# output_filename = "{}/{}/{}.svg".format(directory, foldername, filename)
	output_filename = "{}/{}/{}.png".format(directory, foldername, filename)
	input_filename = "{}/{}/{}.dat".format(directory, foldername, filename)
	contourOutput(input_filename, output_filename)




