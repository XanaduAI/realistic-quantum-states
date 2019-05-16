# Copyright 2019 Xanadu Quantum Technologies Inc.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
"""explores the loss parameter space of a heralded weak cubic phase state,
as first determined in https://arxiv.org/abs/1809.04680."""
import copy

import numpy as np
from scipy.integrate import simps

from matplotlib import pyplot as plt

import strawberryfields as sf
from strawberryfields import ops

from hafnian.quantum import density_matrix


def cubic_phase_state(a, cutoff):
    """The weak cubic phase state.

    |psi> = |0> + i*sqrt(3/2)*a|1> + a*i|3>

    Args:
        a (float): the ON state parameter
        cutoff (int): the Fock basis truncation

    Returns:
        array: the density matrix rho=|psi><psi|
    """
    ket = np.zeros([cutoff], dtype=np.complex128)
    ket[0] = 1.0
    ket[1] = 1j * np.sqrt(3 / 2) * a
    ket[3] = 1j * a
    ket = ket / np.linalg.norm(ket)
    dm = np.einsum("i,j->ij", ket, np.conj(ket))
    return dm


def wigner_log_negativity(rho, l, cutoff):
    r"""Calculates and returns the Wigner Log Negativity
    of the quantum state:

    .. math:: WLN = \log(\int |W| dx dp)

    Args:
        rho (array): density matrix of the quantum state
        l (float): maximum grid size -x_max < x < x_max
            over which to compute the Wigner function
        cutoff (int): the Fock basis truncation

    Returns:
        float: the Wigner log negativity
    """
    x = np.linspace(-l, l, 100)
    p = np.linspace(-l, l, 100)

    Q, P = np.meshgrid(x, p)
    A = (Q + P * 1.0j) / (2 * np.sqrt(2 / 2))

    Wlist = np.array([np.zeros(np.shape(A), dtype=complex) for k in range(cutoff)])

    # Wigner function for |0><0|
    Wlist[0] = np.exp(-2.0 * np.abs(A) ** 2) / np.pi

    # W = rho(0,0)W(|0><0|)
    W = np.real(rho[0, 0]) * np.real(Wlist[0])

    for n in range(1, cutoff):
        Wlist[n] = (2.0 * A * Wlist[n - 1]) / np.sqrt(n)
        W += 2 * np.real(rho[0, n] * Wlist[n])

    for m in range(1, cutoff):
        temp = copy.copy(Wlist[m])
        # Wlist[m] = Wigner function for |m><m|
        Wlist[m] = (2 * np.conj(A) * temp - np.sqrt(m) * Wlist[m - 1]) / np.sqrt(m)

        # W += rho(m,m)W(|m><m|)
        W += np.real(rho[m, m] * Wlist[m])

        for n in range(m + 1, cutoff):
            temp2 = (2 * A * Wlist[n - 1] - np.sqrt(m) * temp) / np.sqrt(n)
            temp = copy.copy(Wlist[n])
            # Wlist[n] = Wigner function for |m><n|
            Wlist[n] = temp2

            # W += rho(m,n)W(|m><n|) + rho(n,m)W(|n><m|)
            W += 2 * np.real(rho[m, n] * Wlist[n])

    return np.log(simps(simps(np.abs(W / 2), p), x))


def circuit(cutoff, l1=0.85, l2=1):
    """Runs the heralded circuit with specified parameters,
    returning the output fidelity to the requested weak cubic phase state,
    the post-selection probability, and the Wigner log negativity.

    Args:
        cutoff (int): the Fock basis truncation
        l1 (float): squeeze cavity loss
        l2 (float): PNR loss

    Returns:
        tuple: a tuple containing the output fidelity to the target ON state,
            the probability of post-selection, the state norm before entering the beamsplitter,
            the state norm after exiting the beamsplitter, and the density matrix of the output state.
    """
    # weak cubic phase state parameter
    a = 0.53
    # the Fock state measurement of mode 0 to be post-selected
    m1 = 1
    # the Fock state measurement of mode 1 to be post-selected
    m2 = 2

    # define target state
    target_state = cubic_phase_state(a, cutoff)

    # gate parameters for the heralded quantum circuit.
    # squeezing magnitudes
    sq_r = [0.71, 0.67, -0.42]
    # squeezing phase
    sq_phi = [-2.07, 0.06, -3.79]
    # displacement magnitudes
    d_r = [-0.02, 0.34, 0.02]
    # beamsplitter theta
    bs_theta1, bs_theta2, bs_theta3 = [-1.57, 0.68, 2.5]
    # beamsplitter phi
    bs_phi1, bs_phi2, bs_phi3 = [0.53, -4.51, 0.72]

    # quantum circuit prior to entering the beamsplitter
    eng, q = sf.Engine(3)
    with eng:
        for k in range(3):
            ops.Sgate(sq_r[k], sq_phi[k]) | q[k]
            ops.Dgate(d_r[k]) | q[k]
            ops.LossChannel(l1) | q[k]

        ops.BSgate(bs_theta1, bs_phi1) | (q[0], q[1])
        ops.BSgate(bs_theta2, bs_phi2) | (q[1], q[2])
        ops.BSgate(bs_theta3, bs_phi3) | (q[0], q[1])
        ops.LossChannel(l2) | q[0]
        ops.LossChannel(l2) | q[1]

    state = eng.run("gaussian", cutoff_dim=cutoff)
    mu = state.means()
    cov = state.cov()

    rho = density_matrix(mu, cov, post_select={0: m1, 1: m2}, cutoff=cutoff, hbar=2)

    # probability of measuring m1 and m2
    prob = np.abs(np.trace(rho))

    # output state
    if prob != 0:
        rho = rho / prob

    # fidelity with the target state
    fidelity = np.abs(np.trace(np.einsum("ij,jk->ik", rho, target_state)))
    return fidelity, prob, rho


def loss_parameter_search(cutoff, filename):
    """Explores the loss parameter space of the prepared weak
    cubic phase states, by varying the values of squeeze cavity
    and PNR loss.

    The results are saved to a NumPy npz file, with arrays
    named ``fid``, ``prob``, and ``wln``.

    Args:
        cutoff (int): the Fock basis truncation
        filename (str): the name of the ``npz`` file
            to save the parameter search results to

    Returns:
        tuple[array]: a tuple of three arrays, containing the
        fidelity, probability, and Wigner log negativity respectively
    """
    dl = 0.05
    l1_list = np.arange(dl, 1 + dl, dl)
    l2_list = np.arange(dl, 1 + dl, dl)

    fid = []
    prob = []
    wln = []

    for l1 in l1_list:
        for l2 in l2_list:
            print("\nRunning simulation with l1 = ", l1, " l2 = ", l2)
            f, p, rho = circuit(cutoff, l1=l1, l2=l2)

            print("fidelity: {}".format(f))
            print("probability: {}".format(p))
            print("rho: {}".format(rho[:4, :4]))

            fid.append([l1, l2, f])
            prob.append([l1, l2, p])

            w = wigner_log_negativity(rho, l=10, cutoff=cutoff)
            wln.append([l1, l2, w])

    fid = np.array(fid)
    prob = np.array(prob)
    wln = np.array(wln)
    wln[:, 2][np.abs(wln[:, 2]) < 1e-5] = 0
    np.savez(filename, fid=fid, prob=prob, wln=wln)

    return fid, prob, wln


def contour_plot(data, ax=None, fig=None, title="None", cmap="Greens"):
    """Generates a contour plot of the loss parameter space
    for the specified quantum state data.

    Args:
        data (array): array of the form `[[eta1, eta2, value], ...]``
        ax (axis): matplotlib axis to draw the contour plot on.
            If none is provided, an axis is created dynamically.
        fig (figure): matplotlib figure to draw the contour plot axis on.
            If none is provided, a figure is created dynamically.
        title (str): the contour plot title
        cmap (str): a valid matplotlib color map

    Returns:
        tuple[figure, axis]: a tuple containing the matplotlib figure and axis
    """
    L1, L1_idx = np.unique(data.T[0], return_inverse=True)
    L2, L2_idx = np.unique(data.T[1], return_inverse=True)
    Z = np.empty(L1.shape + L2.shape)
    Z.fill(np.nan)
    Z[L1_idx, L2_idx] = data.T[2]

    if ax is None and fig is None:
        fig, ax = plt.subplots(1, 1, figsize=(12, 7))
        ax.set_xlabel(r"$\eta_1$")
        ax.set_ylabel(r"$\eta_2$")

    max_abs = np.max(np.abs(Z))
    ax.imshow(
        Z,
        extent=[0, 1, 0, 1],
        origin="lower",
        cmap=cmap,
        alpha=0.5,
        interpolation="hermite",
        aspect="equal",
        vmin=0,
        vmax=max_abs,
    )
    contours = ax.contour(L1, L2, Z, 10, colors="black")
    ax.clabel(contours, inline=True, fontsize=10)

    ax.set_xlim(0.5, 1)
    ax.set_ylim(0.5, 1)

    for item in (
        [ax.title, ax.xaxis.label, ax.yaxis.label] + ax.get_xticklabels() + ax.get_yticklabels()
    ):
        item.set_fontsize(16)

    ax.set_title("{}".format(title))

    if ax is None and fig is None:
        ax.set_xlabel(r"$\eta_1$")
        ax.set_ylabel(r"$\eta_2$")

    return fig, ax


def plot(fid, prob, wln, filename, cmap="Greens"):
    """Generates a 1x3 row of contour plots containing loss parameter space
    plots of state fidelity, probability, and Wigner log negativity.

    Args:
        fid (array): fidelity data of the form `[[eta1, eta2, value], ...]``
        prob (array): probability data of the form `[[eta1, eta2, value], ...]``
        wln (array): Wigner log negativity data of the form `[[eta1, eta2, value], ...]``
        cmap (str): a valid matplotlib color map
        filename (str): the name of the image file to save the plot to
    """
    fig, ax = plt.subplots(1, 3, sharex=False, sharey=False, figsize=(12, 6))
    contour_plot(fid, title="Fidelity", ax=ax[0], fig=fig, cmap=cmap)
    contour_plot(prob, title="Probability", ax=ax[2], fig=fig, cmap=cmap)
    contour_plot(wln, title="Wigner Log Negativity", ax=ax[1], fig=fig, cmap=cmap)

    ax[0].set_ylabel(r"$\eta_1$")
    ax[0].set_xlabel(r"$\eta_2$")
    ax[1].set_xlabel(r"$\eta_2$")
    ax[2].set_xlabel(r"$\eta_2$")

    plt.tight_layout(rect=[0, 0.03, 1, 0.8])
    plt.savefig(filename)


if __name__ == "__main__":
    cutoff = 15
    fid, prob, wln = loss_parameter_search(cutoff, "cubic_phase_a0.53.npz")
    plot(fid, prob, wln, "cubic_phase_a0.53.pdf")
