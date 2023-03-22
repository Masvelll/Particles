import random
import numpy as np
import seaborn as sns
import matplotlib
#matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.animation import PillowWriter


epsilon = 1.0
sigma = 1.0
cut_off = sigma*10
dt = 0.01

t = 0
potential_energy = []
kinetic_energy = []
time = []

size = 5

class Particle:
    
    def __init__(self, pos):
        self.m = 1
        self.x = pos[0]
        self.x_prev = pos[0]
        self.y = pos[1]
        self.y_prev = pos[1]
        self.z = pos[2]
        self.z_prev = pos[2]
        self.force = [0, 0, 0]

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

        vx = (p1.x - p1.x_prev) / dt
        vy = (p1.y - p1.y_prev) / dt
        vz = (p1.z - p1.z_prev) / dt

        kinetic_energy += p1.m / 2 * (vx ** 2 + vy ** 2 + vz ** 2)

    #print("kinetic >> ", kinetic_energy)
    #print("potential >> ", potential_energy)
    return kinetic_energy, potential_energy

def change_pos(x, x_prev):
    
    if x > size:
        #print("bom")
        x -= size
        x_prev -= size
    if x < 0:
        #print("bom")
        x += size
        x_prev += size
    return x, x_prev

def update_position(particles, dt):
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
        particles[i].force = force

def update(particles, dt, n):
    global t, potential_energy, kinetic_energy
    for i in range(n):
        t += dt
        time.append(t)
        
        update_forces(particles)
        update_position(particles, dt)
        
        k, p = calc_energy(particles)
        potential_energy.append(p)
        kinetic_energy.append(k)

def count_speed(particles):
    distr = []
    for p in particles:
        vx = (p.x - p.x_prev) / dt
        vy = (p.y - p.y_prev) / dt
        vz = (p.z - p.z_prev) / dt
        v2 = vx ** 2 + vy ** 2 + vz ** 2
        distr.append(v2)
    return distr

def plot_update(particles):
    update(particles, dt , 1)
    ax.clear()
    ax2.clear()
    ax3.clear()
    ax.set_xlim3d([0, size])
    ax.set_ylim3d([0, size])
    ax.set_zlim3d([0, size])
    img = []

    ax.set_title("Animation")
    ax2.set_title("Energy")
    
    ax2.plot(time, kinetic_energy, label="kinetic_energy")
    ax2.plot(time, potential_energy, label="potential_energy")
    ax2.plot(time, np.array(potential_energy) + np.array(kinetic_energy), label="full")
    ax2.legend()
    #print(count_speed(particles))
    sns.kdeplot(count_speed(particles), ax=ax3)
    
    for p in particles:
        #print(">> ", p.x, p.y, p.z)
        #print(p.force)
        img.append(ax.scatter3D(p.x, p.y, p.z, c="red", s=10))
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



number_of_particles = 60
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

fig = plt.figure("Lennard-Jones", figsize=(12, 6))

ax = fig.add_subplot(121, projection='3d')
ax2 = fig.add_subplot(222)
ax3 = fig.add_subplot(224)



def do_ani():

    anim = animation.FuncAnimation(fig, plot_update, [particles], interval=15, blit=False, repeat=True)
    plt.show()

    fig2 = plt.figure("Energy graphs")
    plt.plot(time, kinetic_energy, label="kinetic")
    plt.plot(time, potential_energy, label="potential")
    plt.plot(time, np.array(potential_energy) + np.array(kinetic_energy), label="full")
    plt.legend()
    plt.show()

def create_gif():
    with writer.saving(fig, "Best_ani_ever3.gif", 100):
        for i in range(500):
            print("frame {} done".format(i))
            plot_update(particles)
            writer.grab_frame()

    print("file created")

do_ani()
#create_gif()



