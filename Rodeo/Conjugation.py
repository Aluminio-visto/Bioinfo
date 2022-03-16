#!/usr/bin/env python

import os
import argparse
import numpy as np
from matplotlib import pyplot as plt

def get_arguments():

    parser = argparse.ArgumentParser(prog='ODE.py', description='ODE solver for SIR conjugative models')
    
    input_group = parser.add_argument_group('Input', 'Input parameters')

    input_group.add_argument('-m', '--model', dest="model", metavar="Model used",
                             type=str, required=False, default='MOB', help='Model used for numerical simulation, whether mobilization (-m MOB) or conjugation (-m CONJ)')

    input_group.add_argument('-f', '--fi', dest="fi", metavar="growth rate",
                             type=float, required=False, default=0.0066, help='Growth rate per second, should be around 0.0006 for E.coli')
    input_group.add_argument('-kenc', '--kenc', metavar="encounter rate",
                             type=float, required=False, default=0.06, help='Area (in microns) covered by a Donor cell per unit time, should be in the 0.1-0.01s')
    input_group.add_argument('-kdis', '--kdis', metavar="failure rate",
                             type=float, required=False, default=0.01, help='Non-succesful contact events per unit time')

    input_group.add_argument('-t', '--tau', metavar="tau", type=float, required=False, default=1200,
                             help='Handling time -from complex D-R to D+Transconjugant- may take hundreds to thousands of seconds')
    input_group.add_argument('-R', '--Recipients', type=float, required=False, default=0.1,
                             help='Initial concentration of recipient cells')
    input_group.add_argument('-D', '--Donors', type=float, required=False, default=0.01,
                             help='Initial concentration of donor cells')
    input_group.add_argument('-T', '--Transconjugants', type=float, required=False, default=0,
                             help='Initial concentration of transconjugant cells')
    arguments = parser.parse_args()

    return arguments

args=get_arguments()

class ODESolver:

    def __init__(self, f):
        self.f=f
    
    def advance(self):
        raise NotImplementedError
    
    def set_initial_conditions(self, U0):     # U0 son las condiciones iniciales
        if isinstance(U0, (int, float)):      # si solo hay un único número, es que solo hay una ecuación (ODE escalar)
            self.number_of_equations=1
            U0 = float(U0)
        else:
            U0=np.asarray(U0)                 # si tenemos un sistema de ecuaciones, U0 no será un nº sino un array
            self.number_of_equations=U0.size  #el nº de ecuaciones coincide con el tamaño del array U0
        self.U0=U0

    def solve(self, time_points):
        self.t       = np.asarray(time_points)                 # Resuelve a lo largo de t puntos
        n            = self.t.size                             # Este es el nº de puntos en el tiempo 
        self.u       = np.zeros((n, self.number_of_equations)) # Hacemos un array vacío para las soluciones (n eqs)x(t puntos)
        self.u[0, :] = self.U0                                 # El punto t0 corresponde a la condicion inicial U0

        for i in range (n-1):
            self.i      = i
            self.u[i+1] = self.advance()
        return self.u[:i+2], self.t[:i+2]

class ForwardEuler(ODESolver):
    def advance(self):
        u, f, i, t = self.u, self.f, self.i, self.t
        dt         = t[i+1]-t[i]
        return u[i, :] + dt*f(u[i, :], t[i])

class MOB:
    def __init__(self, fi, kenc, kdis, kon, D0, R0, T0, DR0):   # conjugation model requires these parameters and initial values

        if isinstance(fi, (float, int)):
            self.fi=lambda t: fi 
        elif callable:
            self.fi=fi

        if isinstance(kenc, (float, int)):
            self.kenc=lambda t: kenc 
        elif callable:
            self.kenc=kenc

        if isinstance(kdis, (float, int)):
            self.kdis=lambda t: kdis 
        elif callable:
            self.kdis=kdis

        if isinstance(kon, (float, int)):
            self.kon=lambda t: kon 
        elif callable:
            self.kon=kon

        self.initial_conditions = [D0, R0, T0, DR0]

    def __call__(self, u, t):

        D, R, T, DR = u

        return np.asarray([
            -self.kenc(t)*D*R + self.kdis(t)*DR + self.kon(t)*DR + self.fi(t)*D,
            -self.kenc(t)*D*R + self.kdis(t)*DR                  + self.fi(t)*R,
                                                  self.kon(t)*DR + self.fi(t)*T,
             self.kenc(t)*D*R - self.kdis(t)*DR - self.kon(t)*DR 
        ])

class CONJ:
    def __init__(self, fi, kenc, kdis, kon, D0, R0, T0, DR0, TR0):   # COnjugation model requires the same parameters and initial values plus T-R intermediate

        if isinstance(fi, (float, int)):
            self.fi=lambda t: fi 
        elif callable:
            self.fi=fi

        if isinstance(kenc, (float, int)):
            self.kenc=lambda t: kenc 
        elif callable:
            self.kenc=kenc

        if isinstance(kdis, (float, int)):
            self.kdis=lambda t: kdis 
        elif callable:
            self.kdis=kdis

        if isinstance(kon, (float, int)):
            self.kon=lambda t: kon 
        elif callable:
            self.kon=kon

        self.initial_conditions = [D0, R0, T0, DR0, TR0]

    def __call__(self, u, t):

        D, R, T, DR, TR = u

        return np.asarray([
            -self.kenc(t)*D*R + self.kdis(t)*DR + self.kon(t)*DR + self.fi(t)*D,
            -self.kenc(t)*D*R + self.kdis(t)*DR                  + self.fi(t)*R -self.kenc(t)*T*R + self.kdis(t)*TR,
                                                  self.kon(t)*DR + self.fi(t)*T -self.kenc(t)*T*R + self.kdis(t)*TR + 2*self.kon(t)*TR,
             self.kenc(t)*D*R - self.kdis(t)*DR - self.kon(t)*DR,
             self.kenc(t)*T*R - self.kdis(t)*TR - 2*self.kon(t)*TR
        ])

if __name__ == "__main__":

    fi   = args.fi
    kenc = args.kenc
    kdis = args.kdis
    tau  = args.tau
    kon  = 1/tau
    D0   = args.Donors
    R0   = args.Recipients
    DR0  = 0
    T0   = 0
    TR0  = 0
    model= args.model

    conjugation = model(fi, kenc, kdis, kon, D0, R0, T0, DR0, TR0)
    solver = ForwardEuler(conjugation)

    solver.set_initial_conditions(conjugation.initial_conditions)

    time_steps = np.linspace(0, 120, 7200)

    u, t = solver.solve(time_steps)

    fig, (ax1,ax2) = plt.subplots(1, 2)
    ax1.plot(t, (u[:,0] + u[:,3]), label="Don")
    ax1.plot(t, u[:,1], label="Rec")
    ax1.plot(t, u[:,2], label="Tc")
    #ax1.plot(t, u[:,3], label="D-R")
    #ax1.plot(t, u[:,4], label="T-R")
    ax1.set_ylabel(r'[Cells/$\mu m² $]')
    ax1.set_xlabel("Time (min.)")
    ax1.set_ylim(10E-07, 1)
    ax1.set_title("Cell concentration through time")
    ax1.legend()
    ax1.set_yscale("log")

    Don   =   u[:,3] + u[:,0]
    ratio =   u[:,2]/(Don)
    ax2.plot(t, ratio, label="Tc/Don")
    ax2.set_yscale('log')
    ax2.set_ylim(10E-07, 1)
    ax2.set_title("Ratio Tc : D through time")
    ax2.set_xlabel("Time (min.)")
    ax2.legend()

    fig.suptitle('Plasmid spreading by conjugation')
    plt.show()



