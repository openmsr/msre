import matplotlib.pyplot as plt
from parameters import *

# 1: radiatior coolant inlet temp
# 2: radiator coolant outlet temp
# 3: radiator air outlet temp
# 4: heat exchanger fuel inlet temp
# 5: heat exchanger fuel node 1 temp
# 6: heat exchanger fuel node 2 temp
# 7: heat exchanger fuel node 3 temp
# 8: heat exchanger fuel node 4 temp
# 9: heat exchanger tube node 1 temp
# 10: heat exchanger tube node 2 temp
# 11: heat exchanger coolant inlet temp
# 12: heat exchanger coolant node 1
# 13: heat exchanger coolant node 2
# 14: heat exchanger coolant node 3
# 15: heat exchanger coolant node 4
# 16: core fuel inlet temp
# 17: k
# 18: C1
# 19: C2
# 20: C3
# 21: C4
# 22: C5
# 23: C6
# 24: core graphite temp
# 25: core fuel node 1 temp
# 26: core fuel node 2 temp

filename1 = "sim_out_1000.0_1.txt"
filename5 = "sim_out_1000.0_5.txt"
filename8 = "sim_out_1000.0_8.txt"

sol1 = []
k_file = open(filename1, 'r')
for k in k_file.readlines():
    sol1.append(k.split())
k_file.close()

sol5 = []
k_file = open(filename5, 'r')
for k in k_file.readlines():
    sol5.append(k.split())
k_file.close()

sol8 = []
k_file = open(filename8, 'r')
for k in k_file.readlines():
    sol8.append(k.split())
k_file.close()

sol1 = [[float(j) for j in s] for s in sol1]
sol5 = [[float(j) for j in s] for s in sol5]
sol8 = [[float(j) for j in s] for s in sol8]

k = 10000
test_pow1 = [(1*s[13]-1) for s in sol1]
test_pow5 = [(5*s[13]-5) for s in sol5]
test_pow8 = [(8*s[13]-8) for s in sol8]
tidx1 = [t[0] for t in enumerate(sol1) if (t[1][0] >= 500.00 and t[1][0] <= 800.00)]
tidx5 = [t[0] for t in enumerate(sol5) if (t[1][0] >= 500.00 and t[1][0] <= 800.00)]
tidx8 = [t[0] for t in enumerate(sol8) if (t[1][0] >= 500.00 and t[1][0] <= 800.00)]
t1 = [s[0] for s in sol1[tidx1[0]:tidx1[-1]]]
t5 = [s[0] for s in sol5[tidx5[0]:tidx5[-1]]]
t8 = [s[0] for s in sol8[tidx8[0]:tidx8[-1]]]

# Create a figure and a 3x1 grid of subplots
fig, axs = plt.subplots(3, 1, figsize=(6, 12))

# Plot data on the first subplot
axs[0].plot(t1, test_pow1[tidx1[0]:tidx1[-1]])
axs[0].set_title('1 MW')
#axs[0].set_xlabel('x')
axs[0].set_ylabel('dP')
axs[0].set_xticklabels([])
axs[0].set_xticks([])

# Plot data on the second subplot
axs[1].plot(t5, test_pow5[tidx5[0]:tidx5[-1]])
axs[1].set_title('5 MW')
#axs[1].set_xlabel('x')
axs[1].set_ylabel('dP')
axs[1].set_xticklabels([])
axs[1].set_xticks([])

# Plot data on the third subplot
axs[2].plot(t8, test_pow8[tidx8[0]:tidx8[-1]])
axs[2].set_title('8 MW')
axs[2].set_xlabel('x')
axs[2].set_ylabel('dP')

plt.tight_layout()
plt.show()