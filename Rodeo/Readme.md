# Modelling plasmid spreading with Python and R

_In this project we are going to build a mathematical model simulating how plasmids spread in a bacterial population, and we will simultae it using two different approaches and languages: Python and R and its library rodeo._

## Introduction 🚀

Plasmids are extra-chromosomal DNA molecules coding antibiotic resistance and virulence genes that may be expressed by their hosts, usually microorganisms. These plasmids can be transferred among hosts (even among really different species) in a process called "conjugation". This process is directional, and it involves the use of a molecular "wire" called Type IV Secretion System, energy spending and closeness between donor and recipient cells.

![tumblr_lmeun1D3ZL1qfdkp5o1_400](https://user-images.githubusercontent.com/77884314/124498511-c003d980-ddbc-11eb-9d1c-8c53e6bc5b21.gif)

This poses a huge problem for infectious disease control, as pathogenic bacteria tend to collect plasmids to resist every antibiotic we know. In the lab where I made my PhD we used to measure the parameters of this process of conjugation in order to make accurate models of plasmid (and their antibiotic resistance genes) spreading.


### The parameters 📋

As conjugation happens in the same timescale than bacterial replication, these models differ from usual human epidemiologics where the infection itself takes a negligible fraction of our life. Hence, the most relevant parameters for the model are:
- Bacterial growth speed, psi (ψ).
- Food or substrate availability ([S]).
- Total bacterial concentration ([N]).
- Initial Donor (plasmid carriers) to Recipient cells ratio ([D]/[R]).
- Conjugation efficiency, gamma (γ).

Besides, conjugation resembles a sexually-transmitted disease in the spatial and temporal constraints limiting the infected (here, Donor) ability to transmit the infectious particle (i.e., plasmid); then, there are two parameters that drive the plasmid spreading and that are of special interest to us:
- Number of infections per hour or per cell cycle caused by one Donor.
- Donor-Recipient distance distribution.

### Building the model 🔧

There are a bunch of mathematical models that have been applied to bacterial conjugation. Maybe the most common is an adaptation (published by Levin in the 70s) of the ubiquituous SIR model in which we exchange Susceptible individuals by recipient (R) cells and infected individuals by donor cells (D). Some models mantain the Recovered population (as some unstable plasmids may be lost after some generations in some carrier cells), but as long as our model runs for short times (i.e., 1-3 hours), plasmids lost by instability may be neglected.

The general adaptation of this model to Donor, Recipient and Transconjugant cells follows the system:

<img src="https://render.githubusercontent.com/render/math?math=\frac{d[D]}{dt} = \psi[D],">

<img src="https://render.githubusercontent.com/render/math?math=\frac{d[R]}{dt} = \psi[R] - \gamma[D][R],">

<img src="https://render.githubusercontent.com/render/math?math=\frac{d[T]}{dt} = \psi[T] + \gamma[D][R],">


To this very basic model, many additions have been included throughout the last 50 years to refine predictions and adjust the outcome to experimental results. In our case example we'll define γ as:

<img src="https://render.githubusercontent.com/render/math?math=\gamma = \frac{k[R]}{1 + \tau[R]},">


where k is the searching rate (the area a single cell covers by unit time) and tau the time it takes for a donor cell to pass the plasmid to a recipient one. 
With this system of ODEs we are prepared to build the computational model and numerically solve it in order to make further predictions. As you will see, the real implemented model is a little bit more complicated than this one, but in essence, they are the same.



## Computational modelling ⚙️



### Building an ODE solver in Python 🔩

First, we'll build up a class that will accept our ODE system and will let us tweak the parameters, add new equations to the system, etc.

```python
class ODESolver:

    def __init__(self, f):
        self.f=f
    
    def advance(self):
        raise NotImplementedError
    
    def set_initial_conditions(self, U0):     
        if isinstance(U0, (int, float)):     
            self.number_of_equations=1
            U0 = float(U0)
        else:
            U0=np.asarray(U0)                
            self.number_of_equations=U0.size 
        self.U0=U0

    def solve(self, time_points):
        self.t       = np.asarray(time_points)                
        n            = self.t.size                            
        self.u       = np.zeros((n, self.number_of_equations))
        self.u[0, :] = self.U0                                

        for i in range (n-1):
            self.i      = i
            self.u[i+1] = self.advance()
        return self.u[:i+2], self.t[:i+2]

class ForwardEuler(ODESolver):
    def advance(self):
        u, f, i, t = self.u, self.f, self.i, self.t
        dt         = t[i+1]-t[i]
        return u[i, :] + dt*f(u[i, :], t[i])
```
And, after the ODE solver, we may build up the model itself for mobilization (plamids that once entered a recipient do not conjugate anymore) and for conjugation (plasmids that are able to conjugate once again from their new hosts): 

```python
class MOB:
    def __init__(self, fi, kenc, kdis, kon, D0, R0, T0, DR0):   
        # Mobilization model requires these parameters and initial values

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
    def __init__(self, fi, kenc, kdis, kon, D0, R0, T0, DR0, TR0):   
        # Conjugation models require the same parameters and initial values plus a T-R intermediate

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
```

### Running the script ⌨️

Now that we have build our model and solver, we may build the rest of the program: this first chunk will serve us to take in different parameters every time we run the script:

```python
def get_arguments():

    parser = argparse.ArgumentParser(prog='ODE.py', description='ODE solver for SIR conjugative models')
    
    input_group = parser.add_argument_group('Input', 'Input parameters')

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
```
This will allow us to tweak the initial concentration of cells, the time it takes a cell to fully accept a plasmid, etcetera. Now we may pass to the last step, solving the system and plotting the output using matplotlib:

```python

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

    conjugation = CONJ(fi, kenc, kdis, kon, D0, R0, T0, DR0, TR0)
    solver = ForwardEuler(conjugation)

    solver.set_initial_conditions(conjugation.initial_conditions)

    time_steps = np.linspace(0, 120, 7200)

    u, t = solver.solve(time_steps)

    fig, (ax1,ax2) = plt.subplots(1, 2)
    ax1.plot(t, (u[:,0] + u[:,3]), label="Don")
    ax1.plot(t, u[:,1], label="Rec")
    ax1.plot(t, u[:,2], label="Tc")
    ax1.plot(t, u[:,3], label="D-R")
    ax1.plot(t, u[:,4], label="T-R")
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


    plt.show()
```

## Output 📦

Thus, with this simple model we may predict how will the system evolve for the first couple hours:

### Mobilization
![Mobilization(tau100)](https://user-images.githubusercontent.com/77884314/158076267-0d3b7bba-f2b5-4aad-911e-50772966d469.png)

### Conjugation
![Conjugation(tau10)](https://user-images.githubusercontent.com/77884314/158076275-98437e3b-f261-48ef-9b91-19ae3a53c40e.png)



## Built with 🛠️

* [Visual CODE Studio](https://code.visualstudio.com/) - My preferred MarkDown editor, and everything else.
* [Python 3.9](https://www.python.org) - Programming language
* [Matplotlib](https://matplotlib.org) - The graphic Python library
* [rodeo](https://cran.r-project.org/rodeoVignette) - A numeric simulator from R



## Author ✒️


* **Jorge Rodríguez Grande** - *Part of my PhD work* - [aluminio](https://github.com/Aluminio-visto)

