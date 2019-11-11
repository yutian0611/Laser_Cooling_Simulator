import random
import time
from tkinter import *

import numpy as np

initial_v = 5
photon_niu = 1
h = 1
number_of_atoms=30
radius=25

height = 800
width = 1500

master = Tk()
master.geometry("1500x800")
w = Canvas(master, width=1500, height=800)
w.pack()
t = 0
dt = 0.02


class Atom:

    def __init__(self, Theta_i, x_i, y_i):
        self.v = np.array([initial_v * np.cos(Theta_i), initial_v * np.sin(Theta_i)])
        self.pos = np.array([x_i, y_i])
        self.r = radius
        self.cross_section = 200
        self.mass = 1
        self.p = self.v * self.mass

    def set_v_from_p(self):
        self.v = self.p / self.mass


def vcollision(b1, b2):
    pos_diff= b1.pos-b2.pos
    v1prime = b1.v - (b1.pos - b2.pos) * np.dot(b1.v - b2.v, b1.pos - b2.pos) / (pos_diff[0]**2+pos_diff[1]**2)
    v2prime = b2.v - (b2.pos - b1.pos) * np.dot(b2.v - b1.v, b2.pos - b1.pos) / (pos_diff[0]**2+pos_diff[1]**2)
    return v1prime, v2prime

def temperature():
    sum = 0
    for a in all_atoms:
        sum += a.v[0] ** 2 + a.v[1] ** 2
    ms = sum / len(all_atoms)
    rms = np.sqrt(ms)
    return rms

all_atoms = []
for i in range(0, number_of_atoms):
    temp_atom = Atom(random.uniform(0, 2 * np.pi), random.uniform(25, 1475), random.uniform(25, 775))
    all_atoms.append(temp_atom)

all_circles = []
for i in range(0, number_of_atoms):
    temp_circle = w.create_oval(all_atoms[i].pos[0] - all_atoms[i].r, all_atoms[i].pos[1] - all_atoms[i].r, all_atoms[i].pos[0] + all_atoms[i].r, all_atoms[i].pos[1] + all_atoms[i].r, fill="green")
    all_circles.append(temp_circle)
    w.update()
'''        
    def single_atom_amplitude(self, any_normal_unit_vector, theta):
        # Calculate atomic form factor first, then use it into the calculation of amplitude
        AFF = self.element_object.c
        l = global_constants.lambda_XR_in_Angstrom
        q = 4 * math.pi * math.sin(theta) / l
        n = any_normal_unit_vector
        Qr = q * (n.dot(self.R))
        for i in range(0, 4):
            AFF = AFF + self.element_object.a[i] * math.exp(-self.element_object.b[i] * (q / (4 * math.pi)) ** 2)

        amplitude = AFF * cmath.exp(global_constants.unit_i * Qr)
        # Uncommentalize the following line and related lines in general_testing.py when a check of Atomic Form Factor figure is wanted.
        # return AFF
        return amplitude
'''

while (1):

    for i in range(0, number_of_atoms):
        if (all_atoms[i].pos[0] < all_atoms[i].r and all_atoms[i].v[0] < 0):
            all_atoms[i].v[0] = - all_atoms[i].v[0]
        if (all_atoms[i].pos[0] > width - all_atoms[i].r and all_atoms[i].v[0] > 0):
            all_atoms[i].v[0] = - all_atoms[i].v[0]
        if (all_atoms[i].pos[1] < all_atoms[i].r and all_atoms[i].v[1] < 0):
            all_atoms[i].v[1] = - all_atoms[i].v[1]
        if (all_atoms[i].pos[1] > height - all_atoms[i].r and all_atoms[i].v[1] > 0):
            all_atoms[i].v[1] = - all_atoms[i].v[1]

        for j in range(i + 1, number_of_atoms):
            pos_diff=all_atoms[i].pos - all_atoms[j].pos
            if (np.sqrt(pos_diff[0]**2+pos_diff[1]**2) <= 2 * all_atoms[i].r and (all_atoms[i].pos - all_atoms[j].pos) @ (all_atoms[i].v - all_atoms[j].v) < 0):
                all_atoms[i].v, all_atoms[j].v = vcollision(all_atoms[i], all_atoms[j])


        all_atoms[i].pos = all_atoms[i].pos + all_atoms[i].v
        w.move(all_circles[i], all_atoms[i].v[0], all_atoms[i].v[1])
    temperature_indicator = w.create_text(200, 20, fill="black", font="Times 20 italic bold", text="The temperature is " + str(round(temperature(), 2)))
    w.update()
    w.delete(temperature_indicator)
    time.sleep(dt)
    t = t + dt
    # print(t)
mainloop()
