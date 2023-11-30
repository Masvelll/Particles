from matplotlib.animation import PillowWriter
from matplotlib import pyplot as plt
import moldynamics as md
import numpy as np
import time

size = 7
B = md.Box(size)
B.GenerateMols(50, 1)
Mol_arr = B.GetMols()


fig = plt.figure(figsize=(12, 6))
ax = fig.add_subplot(121, projection='3d')
ax_en = fig.add_subplot(122)
Ke = []
Pe = []
Time = [0.005 * 5]
#ax.scatter3D(X, Y, Z, s=10, alpha=1)

metadata = dict(title='Movie Test', artist='Matplotlib',
                comment='Movie support!')
writer = PillowWriter(fps=120, metadata=metadata)

def plot_update():
    st = time.time()
    s = 200
    for i in range(s): B.update()
    et = time.time()
    print(f"It took {(et - st) * 1000} miliseconds to update {s} times! ", end='')
    Pe.append(B.get_penergy())
    Ke.append(B.get_kenergy())
    ax.clear()
    ax_en.clear()
    ax.set_xlim3d([0, size])
    ax.set_ylim3d([0, size])
    ax.set_zlim3d([0, size])
    ax_en.plot(Time, Pe, label='potential')
    ax_en.plot(Time, Ke, label='kinetic')
    ax_en.plot(Time, np.array(Pe) + np.array(Ke), label='full')
    ax_en.legend()
    Time.append(Time[-1] + Time[0])
    img = []

    st = time.time()
    particles = B.GetMols()
    particles = np.array(B.GetMols())
    X = particles[:, 0]
    Y = particles[:, 1]
    Z = particles[:, 2]
    
    img.append(ax.scatter3D(X, Y, Z, c="red", s=10, alpha=1))
    et = time.time()
    print("Plotted!")
    print(f"It took {(et - st) * 1000} milliseconds to plot!")
    
    return img

def create_gif():
    with writer.saving(fig, "c_ani.gif", 150):
        frame_number = 400
        for i in range(frame_number):
            # print('\n'*10)
            print(f"========== FRAME {i+1}/{frame_number} ==========")
            plot_update()
            writer.grab_frame()
            

    print("file created")
    
create_gif()