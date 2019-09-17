# Simulating realistic quantum states

> Simulating non-Gaussian quantum state preparation with loss via the [fastest known classical method to calculate loop hafnians](https://github.com/XanaduAI/thewalrus).

This repository contains the source code used to produce the results presented in
*"Simulating realistic non-Gaussian state preparation"* [Phys. Rev. A *100*, 022341 (2019)](https://journals.aps.org/pra/abstract/10.1103/PhysRevA.100.022341)

## Contents

The included script `cubic_phase.py` uses the loop hafnian algorithm to explore the loss
parameter space of a heralded weak cubic phase state, as given in
[Phys. Rev. A *100*, 022341 (2019)](https://journals.aps.org/pra/abstract/10.1103/PhysRevA.100.022341). The script generates
contour plots of the resulting fidelity, Wigner Log Negativity, and probability as a function
of heralding and heralded loss.

## Requirements

To perform the quantum state preparation and parameter search, the script uses the
Gaussian backend of [Strawberry Fields](https://github.com/XanaduAI/strawberryfields), 
 the [The Walrus](https://github.com/XanaduAI/thewalrus) library is required for loop hafnian evaluation,
and [matplotlib](https://matplotlib.org/) for graph generation.

All of these prerequisites can be installed via `pip`:

```bash
pip install strawberryfields thewalrus matplotlib
```

## Authors

Nicolás Quesada, Luke Helt, Josh Izaac, Juan Miguel Arrazola, Rahaneh Shahrokhshahi,
Casey Myers, and Krishna Kumar Sabapathy.

If you are doing any research using this source code, the Hafnian library, and
Strawberry Fields, please cite the following three papers:

> Nicolás Quesada, Luke Helt, Josh Izaac, Juan Miguel Arrazola, Rahaneh Shahrokhshahi,
Casey Myers, and Krishna Kumar Sabapathy. "Simulating realistic non-Gaussian state preparation",
[Phys. Rev. A *100*, 022341 (2019)](https://journals.aps.org/pra/abstract/10.1103/PhysRevA.100.022341). The script generates
> 

> Andreas Björklund, Brajesh Gupt, and Nicolás Quesada. "A faster hafnian formula
> for complex matrices and its benchmarking on the Titan supercomputer", [Journal of Experimental Algorithmics (JEA) 24.1 (2019): 11](https://doi.org/10.1145/3325111).

> Nathan Killoran, Josh Izaac, Nicolás Quesada, Ville Bergholm, Matthew Amy, and
> Christian Weedbrook. "Strawberry Fields: A Software Platform for Photonic Quantum Computing",
> [Quantum, 3, 129 (2019)](https://quantum-journal.org/papers/q-2019-03-11-129/).

## License

This source code is free and open source, released under the Apache License, Version 2.0.
