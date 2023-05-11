import random
import time
import numpy as np
import seaborn as sns
import matplotlib
#matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.animation import PillowWriter
#посчитать коэф диффузии
#понижение температуры --somewhat ok
#скоростной верле -- ok

epsilon = 1.0
sigma = 1.0
cut_off = sigma*10
dt = 0.005
number_of_particles = 50

t = 0
potential_energy = []
kinetic_energy = []
Time = []

size = 7
sum_dr = 0

class Particle:
    
    def __init__(self, pos):
        self.m = 1
        
        self.x = pos[0]
        self.vx = 0
        self.y = pos[1]
        self.vy = 0
        self.z = pos[2]
        self.vz = 0
        self.force = [0, 0, 0]
        self.force_prev = [0, 0, 0]

def calc_distance(p1, p2):
    rx = p1.x - p2.x
    ry = p1.y - p2.y
    rz = p1.z - p2.z

    
    if abs(rx) > size / 2:
        rx -= size * (round(rx / size))
    if abs(ry) > size / 2:
        ry -= size * (round(ry / size))
    if abs(rz) > size / 2:
        rz -= size * (round(rz / size))

    r2 = rx ** 2 + ry ** 2 + rz ** 2

    
    return r2

def calc_forces(p1, p2):
    rx = p1.x - p2.x
    ry = p1.y - p2.y
    rz = p1.z - p2.z

    r2 = calc_distance(p1, p2)

    if r2 < cut_off ** 2:
        if abs(rx) > size / 2:
            rx -= size * (round(rx / size))
        if abs(ry) > size / 2:
            ry -= size * (round(ry / size))
        if abs(rz) > size / 2:
            rz -= size * (round(rz / size))
        r2 = rx ** 2 + ry ** 2 + rz ** 2
            
        if r2 != 0:
            force_r = - 24 * epsilon * ((sigma**6 / r2**3 - 2 * (sigma**12 / r2**6)))
        else:
            return 0, 0, 0

        
        fx = force_r * rx / r2
        fy = force_r * ry / r2
        fz = force_r * rz / r2
        return fx, fy, fz
    return 0, 0, 0

def calc_energy(particles):
    kinetic_energy = 0
    potential_energy = 0
    n = len(particles)
    for i in range(n-1):
        p1 = particles[i]
        for j in range(i, n):
            p2 = particles[j]
            if p1 != p2:
                r2 = calc_distance(p1, p2)
                potential_energy += 4 * epsilon * (sigma ** 12 / r2 **6 \
                                                   - sigma ** 6 / r2 ** 3)


        kinetic_energy += p1.m / 2 * (p1.vx ** 2 + p1.vy ** 2 + p1.vz ** 2)

    #print("kinetic >> ", kinetic_energy)
    #print("potential >> ", potential_energy)
    return kinetic_energy, potential_energy

def change_pos(x):
    
    if x > size:
        #print("bom")
        x -= size
    if x < 0:
        #print("bom")
        x += size
    return x

def update_force(particles, p1):
    force = [0, 0, 0]
    for p2 in particles:
        if p1 != p2:
            forces = calc_forces(p1, p2)
            force[0] += forces[0]
            force[1] += forces[1]
            force[2] += forces[2]
    p1.force = force
def update_forces(particles):
    n = len(particles)
    for i in range(n):
        force = [0, 0, 0]
        for j in range(n):
            if i != j:
                forces = calc_forces(particles[i], particles[j])
                force[0] += forces[0]
                force[1] += forces[1]
                force[2] += forces[2]
        particles[i].force_prev = particles[i].force
        particles[i].force = force

def update_position(particles, dt):
    energy_to_down = False
    for p in particles:
        fx_prev, fy_prev, fz_prev = p.force_prev    
        
        p.x = p.x + p.vx * dt + 1 / 2 * fx_prev / p.m * dt ** 2
        p.y = p.y + p.vy * dt + 1 / 2 * fy_prev / p.m * dt ** 2
        p.z = p.z + p.vz * dt + 1 / 2 * fz_prev / p.m * dt ** 2
        
        p.force_prev = p.force

        p.x = change_pos(p.x)
        p.y = change_pos(p.y)
        p.z = change_pos(p.z)

    for p in particles:
        update_force(particles, p)
        fx, fy, fz = p.force
        fx_prev, fy_prev, fz_prev = p.force_prev
        p.vx = p.vx + 1 / 2 * (fx_prev / p.m + fx / p.m) * dt
        p.vy = p.vy + 1 / 2 * (fy_prev / p.m + fy / p.m) * dt
        p.vz = p.vz + 1 / 2 * (fz_prev / p.m + fz / p.m) * dt

        if p.vx ** 2 + p.vy ** 2 + p.vz ** 2 > 100:
            energy_to_down = True
            
            p.vx *= 0.9
            p.vy *= 0.9
            p.vz *= 0.9
    if energy_to_down: print("Energy decreased! ", end='')

        
                
    

        

def update_position_s(particles, dt):
    for p in particles:
        fx, fy, fz = p.force
        x_prev_next = p.x
        y_prev_next = p.y
        z_prev_next = p.z
        p.x = -p.x_prev + 2 * p.x + fx / p.m * dt ** 2
        p.y = -p.y_prev + 2 * p.y + fy / p.m * dt ** 2
        p.z = -p.z_prev + 2 * p.z + fz / p.m * dt ** 2
        p.x_prev = x_prev_next
        p.y_prev = y_prev_next
        p.z_prev = z_prev_next

        p.x, p.x_prev = change_pos(p.x, p.x_prev)
        p.y, p.y_prev = change_pos(p.y, p.y_prev)
        p.z, p.z_prev = change_pos(p.z, p.z_prev)

        p.vx = (p.x - p.x_prev) / dt
        p.vy = (p.y - p.y_prev) / dt
        p.vz = (p.z - p.z_prev) / dt



def update(particles, dt, n):
    global t, potential_energy, kinetic_energy
    
    for i in range(n):
        t += dt
        Time.append(t)

        #update_forces(particles)
        st = time.time()
        update_position(particles, dt)
        
        k, p = calc_energy(particles)
        potential_energy.append(p)
        kinetic_energy.append(k)
        et = time.time()
        print(f"Took {(et - st) * 1000} ms to update! ")

def count_speed(particles):
    distr = []
    for p in particles:
        v2 = p.vx ** 2 + p.vy ** 2 + p.vz ** 2
        distr.append(v2)
    return distr

def plot_update(particles):
    st = time.time()
    update(particles, dt , 10)
    et = time.time()
    print(f"It took {(et - st) * 1000} miliseconds to update 10 times! ", end='')
    ax.clear()
    ax2.clear()
    ax3.clear()
    ax.set_xlim3d([0, size])
    ax.set_ylim3d([0, size])
    ax.set_zlim3d([0, size])
    img = []

    ax.set_title("Animation")
    ax2.set_title("Energy")
    
    ax2.plot(Time, kinetic_energy, label="kinetic_energy")
    ax2.plot(Time, potential_energy, label="potential_energy")
    ax2.plot(Time, np.array(potential_energy) + np.array(kinetic_energy), label="full")
    ax2.legend()
    #print(count_speed(particles))
    sns.kdeplot(count_speed(particles), ax=ax3)

    st = time.time()
    X, Y, Z = [], [], []
    for p in particles:
        X.append(p.x)
        Y.append(p.y)
        Z.append(p.z)
    
    img.append(ax.scatter3D(X, Y, Z, c="red", s=10, alpha=1))
    et = time.time()
    print("Plotted!")
    print(f"It took {(et - st) * 1000} milliseconds to plot!")
    
    return img

def generate_positions(n):
    pos = []
    l = int(n ** (1 / 3)) + 1
    bound = size * 0.9
    for x in np.linspace(size - bound, bound, l):
        for y in np.linspace(size - bound, bound, l):
            for z in np.linspace(size - bound, bound, l):
                pos.append([x, y, z])
    pos = np.array(pos[:n])
    ran = 0
    #ran = np.random.rand(n) / 50
    x = np.linspace(0, size - 0.5, n)
    #np.random.shuffle(x)
    x += ran
    
    #ran = np.random.rand(n) / 50
    y = np.linspace(0, size - 0.5, n)
    #np.random.shuffle(y)
    y += ran

    #ran = np.random.rand(n) / 50
    z = np.linspace(0, size - 0.5, n)
    #np.random.shuffle(z)
    z += ran

    return pos
    #return np.transpose(np.vstack([x, y, z]))




positions = generate_positions(number_of_particles)
particles = []
#p1 = Particle([1, 1, 1])
#p2 = Particle([4, 1, 1])
#particles.append(p1)
#particles.append(p2)


def create_part():
    for i in range(number_of_particles):
        pos = positions[i]
        p = Particle(pos)
        particles.append(p)
create_part()   

metadata = dict(title='Movie Test', artist='Matplotlib',
                comment='Movie support!')
writer = PillowWriter(fps=120, metadata=metadata)

fig = plt.figure(figsize=(12, 6))

ax = fig.add_subplot(121, projection='3d')
ax2 = fig.add_subplot(222)
ax3 = fig.add_subplot(224)



def do_ani():

    line_ani = animation.FuncAnimation(fig, plot_update, [particles], interval=15, blit=False, repeat=True)
    plt.show()

    plt.plot(Time, kinetic_energy, label="kinetic")
    plt.plot(Time, potential_energy, label="potential")
    plt.plot(Time, np.array(potential_energy) + np.array(kinetic_energy), label="full")
    plt.legend()
    plt.show()

def create_gif():
    with writer.saving(fig, "ani.gif", 150):
        frame_number = 200
        for i in range(frame_number):
            print('\n'*10)
            print(f"========== FRAME {i+1}/{frame_number} ==========")
            plot_update(particles)
            writer.grab_frame()
            

    print("file created")

#for i in range(300):
#    print("frame >> ", i)
#    update(particles, dt , 100)

#do_ani()
create_gif()
#with open("data8.txt", "w") as file:
#    file.write("\n".join(map(str, count_speed(particles))))
#    print(count_speed(particles))




