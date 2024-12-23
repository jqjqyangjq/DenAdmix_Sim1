from collections import defaultdict
import pandas as pd
import numpy as np
from numba import njit
import json
from math import exp, ceil


class HMMParam:
    def __init__(self, state_names, starting_probabilities, transitions, emissions):
        self.state_names = np.array(state_names)
        self.starting_probabilities = np.array(starting_probabilities)
        self.transitions = np.array(transitions)
        self.emissions = np.array(emissions)

    def __str__(self):
        out = f"> state_names = {self.state_names.tolist()}\n"
        out += f"> starting_probabilities = {np.matrix.round(self.starting_probabilities, 4).tolist()}\n"
        out += f"> transitions = {np.matrix.round(self.transitions, 8).tolist()}\n"
        out += f"> emissions = {np.matrix.round(self.emissions, 8).tolist()}"
        return out

    def __repr__(self):
        return f"{self.__class__.__name__}({self.state_names}, {self.starting_probabilities}, {self.transitions}, {self.emissions})"


# Read HMM parameters from a json file
def read_HMM_parameters_from_file(filename):

    if filename is None:
        return get_default_HMM_parameters()

    with open(filename) as json_file:
        data = json.load(json_file)

    return HMMParam(
        state_names=data["state_names"],
        starting_probabilities=data["starting_probabilities"],
        transitions=data["transitions"],
        emissions=data["emissions"],
    )


# Set default parameters
def get_default_HMM_parameters():
    return HMMParam(
        state_names=["Human", "Archaic"],
        starting_probabilities=[0.98, 0.02],
        transitions=[[0.9999, 0.0001], [0.02, 0.98]],
        #emissions = [0.00006, 0.00006]
        emissions=[0.00001896, 0.0002],
    )


@njit
def fwd_step(alpha_prev, E, trans_mat):
    alpha_new = (alpha_prev @ trans_mat) * E
    n = np.sum(alpha_new)
    return alpha_new / n, n


@njit
def forward(probabilities, transitions, init_start):
    n = len(probabilities)
    forwards_in = np.zeros((n, len(init_start)))
    scale_param = np.ones(n)
    for t in range(n):    
        if t == 0:
            forwards_in[t, :] = init_start * probabilities[t, :]
            scale_param[t] = np.sum(forwards_in[t, :])
            forwards_in[t, :] = forwards_in[t, :] / scale_param[t]
        else:
            forwards_in[t, :], scale_param[t] = fwd_step(
                forwards_in[t - 1, :], probabilities[t, :], transitions
            )
    return forwards_in, scale_param


@njit
def bwd_step(beta_next, E, trans_mat, n):
    beta = (trans_mat * E) @ beta_next
    return beta / n


@njit
def backward(emissions, transitions, scales):
    n, n_states = emissions.shape
    beta = np.ones((n, n_states))
    for i in range(n - 1, 0, -1):
        beta[i - 1, :] = bwd_step(beta[i, :], emissions[i, :], transitions, scales[i])
    return beta


def logoutput(pars, loglikelihood, iteration, log_file):
    n_states = len(pars.emissions)
    with open(log_file, "a") as f:
    # Make header
        if iteration == 0:
            print_emissions = "\t".join(["emis{0}".format(x + 1) for x in range(n_states)])
            print_starting_probabilities = "\t".join(
                ["start{0}".format(x + 1) for x in range(n_states)]
            )
            print_transitions = "\t".join(
                ["trans{0}_{0}".format(x + 1) for x in range(n_states)]
            )
            print(
                "it",
                "ll",
                print_starting_probabilities,
                print_emissions,
                print_transitions,
                sep="\t",
                file = f
            )
            print(
                "it",
                "ll",
                print_starting_probabilities,
                print_emissions,
                print_transitions,
                sep="\t"            )

        # Print parameters
        print_emissions = "\t".join([str(x) for x in np.matrix.round(pars.emissions, 8)])
        print_starting_probabilities = "\t".join(
            [str(x) for x in np.matrix.round(pars.starting_probabilities, 5)]
        )
        print_transitions = "\t".join(
            [str(x) for x in np.matrix.round(pars.transitions, 6).diagonal()]
        )
        print(
            iteration,
            round(loglikelihood, 6),
            print_starting_probabilities,
            print_emissions,
            print_transitions,
            sep="\t",
            file = f
        )
        print(
            iteration,
            round(loglikelihood, 6),
            print_starting_probabilities,
            print_emissions,
            print_transitions,
            sep="\t"
        )
def write_post_to_file(Z, chr_index, w_index, outfile):
    Z_df = pd.DataFrame()
    for i in range(len(Z)):
        phase = len(Z[i])
    # Create temporary dataframes for each element
        for phase_ in range(phase):
            temp_df = pd.DataFrame()
            temp_df['Chr'] = [chr_index[i]] * len(w_index[i])  # Assign Chr value
            temp_df['phase'] = phase_
            temp_df['Window'] = [item for item in w_index[i]]  # Extract the values from W[i]
            temp_df['state1'] = [item for item in Z[i][phase_][:,0]]  # Flatten and assign values from Z[i]
            temp_df['state2'] = [item for item in Z[i][phase_][:,1]]
        # Concatenate temporary dataframe to the main dataframe
            Z_df = pd.concat([Z_df, temp_df], ignore_index=True)
    Z_df.to_csv(outfile, index=False)
    
def decode_from_params(raw_obs, chr_index, w, pars, post_file, m_rates, window_size = 1000,):
    n_chr = len(chr_index)
    n_states = len(pars.starting_probabilities)
    n_windows = np.ones(n_chr)
    GLL = []
    SNP2BIN = []
    w_index = []
    Z = []
    E = []
    n_gt = 2
    previous_ll = -np.Inf
    SNP = []
    PG = []
    S = []
    fwd = []
    bwd = []
    scales = []
    for chr in range(n_chr):
        n_windows[chr] = w[chr][-1] - w[chr][0] + 1
        GLL_, SNP2BIN_ = linearize_obs(raw_obs[chr], w[chr])
        w_start = w[chr][0]
        w_index_ = np.arange(w_start, w[chr][-1] + 1)    # including missing windows in between
        SNP2BIN_ -= w_start
        n_snp = len(GLL_)
        n_windows_ = round(n_windows[chr])
        SNP_ = np.zeros((n_snp, n_states, n_gt))
        Z_ = np.zeros((n_windows_, n_states))  # P(Z | O)
        E_ = np.ones((n_windows_, n_states))  # P(O | Z)
        GLL.append(GLL_)
        SNP2BIN.append(SNP2BIN_)
        w_index.append(w_index_)
        Z.append(Z_)
        E.append(E_)
        SNP.append(SNP_)
        PG_ = np.zeros((n_snp, n_states, n_gt))
        PG.append(PG_)
    p = pars.emissions
    for chr in range(n_chr):
        n_windows_ = round(n_windows[chr])
        S = np.zeros(n_windows_)
        update_emissions_scale(E[chr], SNP[chr], GLL[chr], p, SNP2BIN[chr], S, m_rates[chr])
        fwd, scales = forward(E[chr], pars.transitions, pars.starting_probabilities)
        bwd = backward(E[chr], pars.transitions, scales)
        Z[chr][:] = fwd * bwd
    write_post_to_file(Z, chr_index, w_index, post_file)

def write_HMM_to_file(hmmparam, outfile):
    data = {key: value.tolist() for key, value in vars(hmmparam).items()}
    json_string = json.dumps(data, indent = 2)
    with open(outfile, 'w') as out:
        out.write(json_string)
