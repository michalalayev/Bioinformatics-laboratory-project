# This script can be run from the terminal inside IDE, since plt.show() does not work from the regular terminal
# Argument to pass: the pkl file name of the plot to show

import pickle
import matplotlib.pyplot as plt
import sys

plot_file_name = sys.argv[1]

plot_file = open(plot_file_name, 'rb')
ax = pickle.load(plot_file)
plt.show()
print("done")