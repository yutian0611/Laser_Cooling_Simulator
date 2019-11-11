import random
import time
from tkinter import *

import numpy as np

speed_up=1
initial_v = 5
speed_of_light = 20
h = 1
number_of_atoms = 1
radius = 200
intensity=1
mass_of_atoms=0.5

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
        self.mass = mass_of_atoms
        self.p = self.v * self.mass

    def set_v_from_p(self):
        self.v = self.p / self.mass
    def set_p_from_v(self):
        self.p = self.v * self.mass


class Photon:

    def __init__(self, lamb, Theta_i, x_i, y_i):
        self.pos = np.array([x_i, y_i])
        self.v = np.array([speed_of_light * np.cos(Theta_i), speed_of_light * np.sin(Theta_i)])
        self.wavelength = lamb
        self.p = np.array([(1 / self.wavelength) * np.cos(Theta_i), (1 / self.wavelength) * np.sin(Theta_i)])


all_incident_photons = []
all_emitted = []


def vcollision(b1, b2):
    pos_diff = b1.pos - b2.pos
    v1prime = b1.v - (b1.pos - b2.pos) * np.dot(b1.v - b2.v, b1.pos - b2.pos) / (pos_diff[0] ** 2 + pos_diff[1] ** 2)
    v2prime = b2.v - (b2.pos - b1.pos) * np.dot(b2.v - b1.v, b2.pos - b1.pos) / (pos_diff[0] ** 2 + pos_diff[1] ** 2)
    return v1prime, v2prime


all_atoms = []
for i in range(0, number_of_atoms):
    temp_atom = Atom(0, 300, 400)
    all_atoms.append(temp_atom)

all_circles = []
for i in range(0, number_of_atoms):
    temp_circle = w.create_oval(all_atoms[i].pos[0] - all_atoms[i].r, all_atoms[i].pos[1] - all_atoms[i].r, all_atoms[i].pos[0] + all_atoms[i].r, all_atoms[i].pos[1] + all_atoms[i].r, fill="green")
    all_circles.append(temp_circle)
    w.update()

all_incident_dots = []
all_emitted_dots = []

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
    # a source that generates photons at a specific rate
    for i in range(0, intensity):
        temp_rand = int(np.random.normal(loc=25.0, scale=5.0, size=None))
        if (random.uniform(0, 1) > 0.5 and temp_rand > 20 and temp_rand < 30):
            temp_photon = Photon(temp_rand, np.pi, width, random.uniform(height / 2 - 100, height / 2 + 100))
            all_incident_photons.append(temp_photon)
            temp_dot = w.create_oval(temp_photon.pos[0] - 10, temp_photon.pos[1] - 10, temp_photon.pos[0] + 10, temp_photon.pos[1] + 10, fill="#" + str(temp_photon.wavelength - 20) + "0" + str(30 - temp_photon.wavelength))
            all_incident_dots.append(temp_dot)

    # behavior of incident photons
    for i in all_incident_photons:
        if (i.pos[0] < -100 or i.pos[0] > width + 100 or i.pos[1] < -100 or i.pos[1] > height + 100):
            all_incident_dots.remove(all_incident_dots[all_incident_photons.index(i)])
            all_incident_photons.remove(i)

        else:
            i.pos = i.pos + i.v
            w.move(all_incident_dots[all_incident_photons.index(i)], i.v[0], i.v[1])
    # behavior of atoms
    for i in range(0, number_of_atoms):
        # collision at walls
        if (all_atoms[i].pos[0] < radius and all_atoms[i].v[0] < 0):
            all_atoms[i].v[0] *= -1
            all_atoms[i].set_p_from_v()
        if (all_atoms[i].pos[0] > width - radius and all_atoms[i].v[0] > 0):
            all_atoms[i].v[0] *= -1
            all_atoms[i].set_p_from_v()
        if (all_atoms[i].pos[1] < radius and all_atoms[i].v[1] < 0):
            all_atoms[i].v[1] *= -1
            all_atoms[i].set_p_from_v()
        if (all_atoms[i].pos[1] > height - radius and all_atoms[i].v[1] > 0):
            all_atoms[i].v[1] *= -1
            all_atoms[i].set_p_from_v()

        # collision between atoms
        for j in range(i + 1, number_of_atoms):
            pos_diff = all_atoms[i].pos - all_atoms[j].pos
            if (np.sqrt(pos_diff[0] ** 2 + pos_diff[1] ** 2) <= 2 * all_atoms[i].r and (all_atoms[i].pos - all_atoms[j].pos) @ (all_atoms[i].v - all_atoms[j].v) < 0):
                all_atoms[i].v, all_atoms[j].v = vcollision(all_atoms[i], all_atoms[j])
                all_atoms[i].set_p_from_v()
                all_atoms[j].set_p_from_v()
        # absorption of photons
        for k in all_incident_photons:
            pos_diff = all_atoms[i].pos - k.pos
            distance = np.sqrt(pos_diff[0] ** 2 + pos_diff[1] ** 2)
            if (distance <= radius and (all_atoms[i].v[0] - k.v[0]) / k.wavelength > 0.95 and (all_atoms[i].v[0] - k.v[0]) / k.wavelength < 1.05):
                #print('p_i=', all_atoms[i].p)
                #print('v_i=', all_atoms[i].v)
                all_atoms[i].p = all_atoms[i].p + k.p
                all_atoms[i].set_v_from_p()
                #print('p_f=', all_atoms[i].p)
                #print('v_f=', all_atoms[i].v)
                w.delete(all_incident_dots[all_incident_photons.index(k)])
                all_incident_dots.remove(all_incident_dots[all_incident_photons.index(k)])
                all_incident_photons.remove(k)

                # an absorption is accompanied by an emission of photons
                temp_emitted = Photon(20, random.uniform(0, 2 * np.pi), all_atoms[i].pos[0], all_atoms[i].pos[1])
                temp_dot = w.create_oval(temp_emitted.pos[0] - 10, temp_emitted.pos[1] - 10, temp_emitted.pos[0] + 10, temp_emitted.pos[1] + 10 , fill="black")
                all_emitted.append(temp_emitted)
                #print(len(all_emitted))
                all_emitted_dots.append(temp_dot)
                #print(len(all_emitted_dots))

                all_atoms[i].p = all_atoms[i].p - temp_emitted.p
                all_atoms[i].set_v_from_p()

        all_atoms[i].pos = all_atoms[i].pos + all_atoms[i].v
        w.move(all_circles[i], all_atoms[i].v[0], all_atoms[i].v[1])

    for e in all_emitted:
        if (e.pos[0] < -100 or e.pos[0] > width + 100 or e.pos[1] < -100 or e.pos[1] > height + 100):
            all_emitted_dots.remove(all_emitted_dots[all_emitted.index(e)])
            all_emitted.remove(e)
        else:
            e.pos = e.pos + e.v
            w.move(all_emitted_dots[all_emitted.index(e)], e.v[0], e.v[1])
    w.update()
    time.sleep(dt/speed_up)
    print(all_atoms[0].v)
    t = t + dt
    # print(t)
mainloop()
