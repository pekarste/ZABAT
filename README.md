# ZABAT - Thermodynamic model

ZABAT is a project about Zinc-air batteries, and this part of the project is about creating a *thermodynamic model*. In this model there are two main parts:
1. Pourbaix diagram of Zn with temperature dependence
2. Concentration profile for Zn in contact with additives (inorganic and organic)
These models are purely based on thermodynamics, and does not say anything about the kinetics of the processes. 

The projecct is written in [Python] and the last deliverable is an app made in [Streamlit] combining these two parts into one.
This repo holds the modelling work associated with m-era.net ZABAT WP 3: Dielectric.

## Pourbaix diagram
A Pourbaix diagram in general is a diagram based purely on thermodynamic properties and equilibria conditions. It describes what compound/species is most stable for different potentials and pH. The diagram is made up of intersecting lines representing different equilibria creating different domains. The lines can be either:

1. Horisontal: An equilibrium only dependent on the potential. Purely electrochemical.
2. Sloped line: An equilibrium depending on both the potential and the pH. Both electrochemical and chemical
3. Vertical line: An equilibrium only dependent on pH. Purely chemical

Usually the temperature and pressure is set/fixed so that the pH and the potentials are the only variables changing, also making a 2D representation possible. In addition, there is usually assumed a total activity of dissolved species, this also influences the equilibria occuring. This activity is often assumed to be equal to the concentration. It is also customary to add the lines for the *OER* and *HER* in the diagram, since this have practical influences for the metal in contact with water. 

The thermodynamic data used in this project is gathered from BEVERSKOG et al. [source] and SI Chemical Data [source], and assumed to be able to describe the equilibria in the temperature range 25-100 degrees celsius. 
[Add picture here]

## Concentration profile
The concentration profile is just a concentration vs pH diagram showing the amount of dissolved species in the solution based on a start concentration which is assumed. It shows how the concentration of species varies for different pH. These concentrations are governed by a set of chemical equilibria described by the law of mass action, giving rise to a set of equations which must be solved simultaneously. It is purely based on thermodynamics. The equilibrium constant for the different reactions is gathered from  ... [source]. Similar projects have been done by ...
The concentration profile of different Zn-species is studied both with inorganic additives $KOH$ (which effectively just decides the pH and gives more $OH^{-}$) and $NH_{3}$, but also organic additives like...

# Running the program
In order to run the program you need a working directory with Python, and the following packages if they are not installed:
1. Numpy
2. Matplotlib
3. Scipy
Which you can use pip to install or conda if you use anaconda.

```bash
pip install numpy, matplotlib, scipy
```

# Supports
The main contributors to the *Thermodynamic model* are
* Pål Emil England Karstensen [england1501@gmail.com](mailto:england1501@gmail.com) and [pal.karstensen@sintef.no](mailto:pal.karstensen@sintef.no). Mainly the Pourbaix diagram
* Sidsel Meli Hanetho [SidselMeli.Hanetho@sintef.no](mailto:SidselMeli.Hanetho@sintef.no). Mainly the concentration profile
Both have contributed to 
