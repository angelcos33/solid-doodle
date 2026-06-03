from pathlib import Path
import random
import matplotlib.pyplot as plt
import numpy as np

def plot_data(data):
    fig, ax = plt.subplots()
    ax.plot(data[:, 0], data[:, 1], lw=2)
    ax.plot(data[:, 0], data[:, 2], lw=2)
    ax.plot(data[:, 0], data[:, 3], lw=2)

datadir = Path("./examples/data")
datafiles = [f for f in datadir.glob("*.txt")]
I_p = []

for file in datafiles:
    print(file.stem)
    data = np.loadtxt(file,
                    delimiter='\t',
                    skiprows=3)
    
    plot_data(data)
    print(data[1, :])
    I_p.append(data[:, 3])
print(I_p)

# Convert degrees into radians
data[:, 0] = np.deg2rad(data[:, 0])
I_p_isotropic = I_p[0] * 0.34 + I_p[1] * 0.67

#data[0][:, 3] = data[0][:, 3] * 0.34 + data[1][:, 3] * 0.67
res = np.transpose(np.vstack((data[:, 0], I_p_isotropic))) 
print(res)

np.savetxt(datafiles[0].with_stem("setfos_simple_emitter_isotropic"),
            res,
            delimiter='\t')

plt.show()