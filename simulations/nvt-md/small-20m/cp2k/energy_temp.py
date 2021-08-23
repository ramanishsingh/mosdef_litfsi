import os
import numpy as np
import mdtraj as md
import matplotlib.pyplot as plt
import signac
from subprocess import call
import os

project = signac.get_project()


def temp(job):

    with job:
        print(job.ws)
        os.chdir(job.ws)
        call("awk '{print $4}' litfsi-1.ener > temp.txt", shell=True)

        fig, ax = plt.subplots()
        f = open("temp.txt", "r")
        temp = []
        a = f.readline()

        while a is not "":
            a = f.readline()
            a.strip()
            if a == "":
                break

            temp.append(float(a))
        f.close()

        call("awk '{print $2}' litfsi-1.ener > time.txt", shell=True)

        f = open("time.txt", "r")
        time = []
        a = f.readline()

        while a is not "":
            a = f.readline()
            a.strip()
            if a == "":
                break

            time.append(float(a) / 1000)

        f.close()

        print(len(time))
        print(len(temp))
        print(time[0:5])
        print(temp[0:5])
        print(time[500])
        print(temp[500])
        plt.plot(time, temp)
        plt.xlabel("Time (ps)")
        plt.ylabel("Temperature (K)")

        np.savetxt(
            "temp",
            np.transpose(np.vstack([time, temp])),
            header="time (ps)\tTemp (K)\t",
        )
        os.remove("time.txt")
        os.remove("temp.txt")
        plt.savefig("temperature.pdf")
        plt.close()


def pot(job):
    with job:
        os.chdir(job.ws)
        print(job.ws)
        call("awk '{print $5}' litfsi-1.ener > pot.txt", shell=True)

        fig, ax = plt.subplots()
        f = open("pot.txt", "r")
        pot = []
        a = f.readline()

        while a is not "":
            a = f.readline()
            a.strip()
            if a == "":
                break

            pot.append(float(a))
        f.close()

        call("awk '{print $2}' litfsi-1.ener > time.txt", shell=True)

        f = open("time.txt", "r")
        time = []
        a = f.readline()

        while a is not "":
            a = f.readline()
            a.strip()
            if a == "":
                break

            time.append(float(a) / 1000)
        f.close()
        print(len(time))
        print(len(pot))
        plt.plot(time, pot)
        plt.xlabel("Time (ps)")
        plt.ylabel("Potential energy (Ha)")
        plt.ticklabel_format(useOffset=False)
        plt.tight_layout()
        np.savetxt(
            "pot.txt",
            np.transpose(np.vstack([time, pot])),
            header="time (ps)\tPotential energy (Ha)\t",
        )
        os.remove("time.txt")
        os.remove("pot.txt")

        plt.savefig("pot.pdf")
        plt.close()


if __name__ == "__main__":
    for job in project.find_jobs():
        temp(job)
        pot(job)
