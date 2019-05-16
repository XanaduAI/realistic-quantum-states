# Simulating realistic quantum states

#### Simulating non-Gaussian quantum state preparation with loss via the
[fastest known classical method to calculate loop hafnians](https://github.com/XanaduAI/hafnian).

This repository contains the source code used to produce the results presented in
*"Simulating realistic non-Gaussian state preparation"* [arXiv:todo](https://arxiv.org/abs/1809.04680).

## Contents

The following three scripts use the loop hafnian algorithm to explore the loss
parameter space for three specific non-Gaussian state preparation schemes. Each script generates
contour plots of the resulting fidelity, Wigner Log Negativity, and probability as a function
of heralding and heralded loss.

* `cubic_phase.py`: explores the loss parameter space of a heralded weak cubic phase state,
  as first determined in [arXiv:1809.04680](https://arxiv.org/abs/1809.04680).

* `fock.py`: explores the loss parameter space of a heralded single photon state.

* `cat_state.py`: explores the loss parameter space of a heralded cat state.


## Requirements

To perform the quantum state preparation and parameter search, these scripts use the
Gaussian backend of [Strawberry Fields](https://github.com/XanaduAI/strawberryfields).
In addition, the `hafnian` library is required for efficient computation of loop hafnians,
and `matplotlib` for graph generation.

All of these prerequisites can be installed via `pip`:

```bash
pip install strawberryfields hafnian matplotlib
```

## Authors

Nicolás Quesada, Luke Helt, Josh Izaac, Juan Miguel Arrazola, Rahaneh Shahrokhshahi,
Casey Myers, and Krishna Kumar Sabapathy.

If you are doing any research using this source code, the Hafnian library, and
Strawberry Fields, please cite the following three papers:

> Nicolás Quesada, Luke Helt, Josh Izaac, Juan Miguel Arrazola, Rahaneh Shahrokhshahi,
Casey Myers, and Krishna Kumar Sabapathy. "Simulating realistic non-Gaussian state preparation",
> arXiv, 2019. arxiv:1805.12498

> Andreas Björklund, Brajesh Gupt, and Nicolás Quesada. "A faster hafnian formula
> for complex matrices and its benchmarking on the Titan supercomputer", arXiv, 2018.
> arxiv:1805.12498

> Nathan Killoran, Josh Izaac, Nicolás Quesada, Ville Bergholm, Matthew Amy, and
> Christian Weedbrook. "Strawberry Fields: A Software Platform for Photonic Quantum Computing",
> Quantum, 3, 129 (2019).

## License

This source code is free and open source, released under the Apache License, Version 2.0.
