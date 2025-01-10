import numpy as np
from tqdm import tqdm
import json
import time
import datetime
from ldpc import BpOsdDecoder
from bposd.css import css_code
import scipy


class css_decode_sim:
    """
    A class for simulating BP+OSD decoding of CSS codes

    Note
    ....
    The input parameters can be entered directly or as a dictionary.

    Parameters
    ----------

    hx: numpy.ndarray
        The hx matrix of the CSS code.
    hz: numpy.ndarray
        The hz matrix of the CSS code.
    error_rate: float
        The physical error rate on each qubit.
    xyz_error_bias: list of ints
        The relative bias for X, Y and Z errors.
    seed: int
        The random number generator seed.
    target_runs: int
        The number of runs you wish to simulate.
    bp_method: string
        The BP method. Choose either: 1) "minimum_sum"; 2) "product_sum".
    ms_scaling_factor: float
        The minimum sum scaling factor (if applicable)
    max_iter: int
        The maximum number of iterations for BP.
    osd_method: string
        The OSD method. Choose from: 1) "osd_cs"; 2) "osd_e"; 3) "osd0".
    channel_update: string
        The channel update method. Choose form: 1) None; 2) "x->z"; 3) "z->x".
    output_file: string
        The output file to write to.
    save_interval: int
        The time in interval (in seconds) between writing to the output file.
    check_code: bool
        Check whether the CSS code is valid.
    tqdm_disable: bool
        Enable/disable the tqdm progress bar. If you are running this script on a HPC
        cluster, it is recommend to disable tqdm.
    run_sim: bool
        If enabled (default), the simulation will start automatically.
    hadamard_rotate: bool
        Toggle Hadamard rotate. ON: 1; OFF; 0
    hadamard_rotate_sector1_length: int
        Specifies the number of qubits in sector 1 for the Hadamard rotation.
    error_bar_precision_cutoff: float
        The simulation will stop after this precision is reached.
    """

    def __init__(self, hx=None, hz=None, **input_dict):
        # default input values
        default_input = {
            "error_rate": None,
            "xyz_error_bias": [1, 1, 1],
            "target_runs": 100,
            "seed": 0,
            "bp_method": "minimum_sum",
            "ms_scaling_factor": 0.625,
            "max_iter": 0,
            "osd_method": "osd_cs",
            "osd_order": 2,
            "save_interval": 2,
            "output_file": None,
            "check_code": 1,
            "tqdm_disable": 0,
            "run_sim": 1,
            "channel_update": "x->z",
            "hadamard_rotate": 0,
            "hadamard_rotate_sector1_length": 0,
            "error_bar_precision_cutoff": 1e-3,
        }

        # apply defaults for keys not passed to the class
        for key in input_dict.keys():
            self.__dict__[key] = input_dict[key]
        for key in default_input.keys():
            if key not in input_dict:
                self.__dict__[key] = default_input[key]

        # output variables
        output_values = {
            "K": None,
            "N": None,
            "start_date": None,
            "runtime": 0.0,
            "runtime_readable": None,
            "run_count": 0,
            "bp_converge_count_x": 0,
            "bp_converge_count_z": 0,
            "bp_success_count": 0,
            "bp_logical_error_rate": 0,
            "bp_logical_error_rate_eb": 0,
            "osd0_success_count": 0,
            "osd0_logical_error_rate": 0.0,
            "osd0_logical_error_rate_eb": 0.0,
            "osdw_success_count": 0,
            "osdw_logical_error_rate": 0.0,
            "osdw_logical_error_rate_eb": 0.0,
            "osdw_word_error_rate": 0.0,
            "osdw_word_error_rate_eb": 0.0,
            "min_logical_weight": 1e9,
        }

        for key in output_values.keys():  # copies initial values for output attributes
            if key not in self.__dict__:
                self.__dict__[key] = output_values[key]

        # the attributes we wish to save to file
        temp = []
        for key in self.__dict__.keys():
            if key not in [
                "channel_probs_x",
                "channel_probs_z",
                "channel_probs_y",
                "hx",
                "hz",
            ]:
                temp.append(key)
        self.output_keys = temp

        # random number generator setup
        if self.seed == 0 or self.run_count != 0:
            self.seed = np.random.randint(low=1, high=2**32 - 1)
        np.random.seed(self.seed)
        print(f"RNG Seed: {self.seed}")

        # the hx and hx matrices
        self.hx = scipy.sparse.csr_matrix(hx).astype(np.uint8)
        self.hz = scipy.sparse.csr_matrix(hz).astype(np.uint8)
        self.N = self.hz.shape[1]  # the block length
        if (
            self.min_logical_weight == 1e9
        ):  # the minimum observed weight of a logical operator
            self.min_logical_weight = self.N
        self.error_x = np.zeros(self.N).astype(np.uint8)  # x_component error vector
        self.error_z = np.zeros(self.N).astype(np.uint8)  # z_component error vector

        # construct the CSS code from hx and hz
        self._construct_code()

        # setup the error channel
        self._error_channel_setup()

        # setup the BP+OSD decoders
        self._decoder_setup()

        if self.run_sim:
            self.run_decode_sim()

    def _single_run(self):
        """
        The main simulation procedure
        """

        # randomly generate the error
        self.error_x, self.error_z = self._generate_error()

        if self.channel_update is None:
            # decode z
            synd_z = self.hx @ self.error_z % 2
            self.bpd_z.decode(synd_z)

            # decode x
            synd_x = self.hz @ self.error_x % 2
            self.bpd_x.decode(synd_x)

        elif self.channel_update == "z->x":
            # decode z
            synd_z = self.hx @ self.error_z % 2
            self.bpd_z.decode(synd_z)

            # update the channel probability
            self._channel_update(self.channel_update)

            # decode x
            synd_x = self.hz @ self.error_x % 2
            self.bpd_x.decode(synd_x)

        elif self.channel_update == "x->z":
            # decode x
            synd_x = self.hz @ self.error_x % 2
            self.bpd_x.decode(synd_x)

            # update the channel probability
            self._channel_update(self.channel_update)

            # decode z
            synd_z = self.hx @ self.error_z % 2
            self.bpd_z.decode(synd_z)

        # compute the logical and word error rates
        self._encoded_error_rates()

    def _channel_update(self, update_direction):
        """
        Function updates the channel probability vector for the second decoding component
        based on the first. The channel probability updates can be derived from Bayes' rule.
        """

        # x component first, then z component
        if update_direction == "x->z":
            decoder_probs = np.zeros(self.N)
            for i in range(self.N):
                if self.bpd_x.osdw_decoding[i] == 1:
                    if (self.channel_probs_x[i] + self.channel_probs_y[i]) == 0:
                        decoder_probs[i] = 0
                    else:
                        decoder_probs[i] = self.channel_probs_y[i] / (
                            self.channel_probs_x[i] + self.channel_probs_y[i]
                        )
                elif self.bpd_x.osdw_decoding[i] == 0:
                    decoder_probs[i] = self.channel_probs_z[i] / (
                        1 - self.channel_probs_x[i] - self.channel_probs_y[i]
                    )

            self.bpd_z.update_channel_probs(decoder_probs)

        # z component first, then x component
        elif update_direction == "z->x":
            self.bpd_z.osdw_decoding
            decoder_probs = np.zeros(self.N)
            for i in range(self.N):
                if self.bpd_z.osdw_decoding[i] == 1:
                    if (self.channel_probs_z[i] + self.channel_probs_y[i]) == 0:
                        decoder_probs[i] = 0
                    else:
                        decoder_probs[i] = self.channel_probs_y[i] / (
                            self.channel_probs_z[i] + self.channel_probs_y[i]
                        )
                elif self.bpd_z.osdw_decoding[i] == 0:
                    decoder_probs[i] = self.channel_probs_x[i] / (
                        1 - self.channel_probs_z[i] - self.channel_probs_y[i]
                    )

            self.bpd_x.update_channel_probs(decoder_probs)

    def _encoded_error_rates(self):
        """
        Updates the logical and word error rates for OSDW, OSD0 and BP (before post-processing)
        """

        # OSDW Logical error rate
        # calculate the residual error
        residual_x = (self.error_x + self.bpd_x.osdw_decoding) % 2
        residual_z = (self.error_z + self.bpd_z.osdw_decoding) % 2

        # check for logical X-error
        if (self.lz @ residual_x % 2).any():
            logical_weight = np.sum(residual_x)
            if logical_weight < self.min_logical_weight:
                self.min_logical_weight = int(logical_weight)

        # check for logical Z-error
        elif (self.lx @ residual_z % 2).any():
            logical_weight = np.sum(residual_z)
            if logical_weight < self.min_logical_weight:
                self.min_logical_weight = int(logical_weight)
        else:
            self.osdw_success_count += 1

        # compute logical error rate
        self.osdw_logical_error_rate = 1 - self.osdw_success_count / self.run_count
        self.osdw_logical_error_rate_eb = np.sqrt(
            (1 - self.osdw_logical_error_rate)
            * self.osdw_logical_error_rate
            / self.run_count
        )

        # compute word error rate
        self.osdw_word_error_rate = 1.0 - (1 - self.osdw_logical_error_rate) ** (
            1 / self.K
        )
        self.osdw_word_error_rate_eb = (
            self.osdw_logical_error_rate_eb
            * ((1 - self.osdw_logical_error_rate_eb) ** (1 / self.K - 1))
            / self.K
        )

        # OSD0 logical error rate
        # calculate the residual error
        residual_x = (self.error_x + self.bpd_x.osd0_decoding) % 2
        residual_z = (self.error_z + self.bpd_z.osd0_decoding) % 2

        # check for logical X-error
        if (self.lz @ residual_x % 2).any():
            logical_weight = np.sum(residual_x)
            if logical_weight < self.min_logical_weight:
                self.min_logical_weight = int(logical_weight)

        # check for logical Z-error
        elif (self.lx @ residual_z % 2).any():
            logical_weight = np.sum(residual_z)
            if logical_weight < self.min_logical_weight:
                self.min_logical_weight = int(logical_weight)
        else:
            self.osd0_success_count += 1

        # compute logical error rate
        self.osd0_logical_error_rate = 1 - self.osd0_success_count / self.run_count
        self.osd0_logical_error_rate_eb = np.sqrt(
            (1 - self.osd0_logical_error_rate)
            * self.osd0_logical_error_rate
            / self.run_count
        )

        # compute word error rate
        self.osd0_word_error_rate = 1.0 - (1 - self.osd0_logical_error_rate) ** (
            1 / self.K
        )
        self.osd0_word_error_rate_eb = (
            self.osd0_logical_error_rate_eb
            * ((1 - self.osd0_logical_error_rate_eb) ** (1 / self.K - 1))
            / self.K
        )

        # BP Logical error rate
        # check for convergence
        if self.bpd_z.converge:
            self.bp_converge_count_z += 1
        if self.bpd_x.converge:
            self.bp_converge_count_x += 1

        if self.bpd_z.converge and self.bpd_x.converge:
            # calculate the residual error
            residual_x = (self.error_x + self.bpd_x.bp_decoding) % 2
            residual_z = (self.error_z + self.bpd_z.bp_decoding) % 2

            # check for logical X-error
            if (self.lz @ residual_x % 2).any():
                pass

            # check for logical Z-error
            elif (self.lx @ residual_z % 2).any():
                pass
            else:
                self.bp_success_count += 1

        # compute logical error rate
        self.bp_logical_error_rate = 1 - self.bp_success_count / self.run_count
        self.bp_logical_error_rate_eb = np.sqrt(
            (1 - self.bp_logical_error_rate)
            * self.bp_logical_error_rate
            / self.run_count
        )

        # compute word error rate
        self.bp_word_error_rate = 1.0 - (1 - self.bp_logical_error_rate) ** (1 / self.K)
        self.bp_word_error_rate_eb = (
            self.bp_logical_error_rate_eb
            * ((1 - self.bp_logical_error_rate_eb) ** (1 / self.K - 1))
            / self.K
        )

    def _construct_code(self):
        """
        Constructs the CSS code from the hx and hz stabilizer matrices.
        """

        print("Constructing CSS code from hx and hz matrices...")
        if isinstance(self.hx, (np.ndarray, scipy.sparse.spmatrix)) and isinstance(
            self.hz, (np.ndarray, scipy.sparse.spmatrix)
        ):
            qcode = css_code(self.hx, self.hz)
            self.lx = qcode.lx
            self.lz = qcode.lz
            self.K = qcode.K
            self.N = qcode.N
            print("Checking the CSS code is valid...")
            if self.check_code and not qcode.test():
                raise Exception(
                    "Error: invalid CSS code. Check the form of your hx and hz matrices!"
                )
        else:
            raise Exception("Invalid object type for the hx/hz matrices")
        return None

    def _error_channel_setup(self):
        """
        Sets up the error channels from the error rate and error bias input parameters
        """

        xyz_error_bias = np.array(self.xyz_error_bias)
        if xyz_error_bias[0] == np.inf:
            self.px = self.error_rate
            self.py = 0
            self.pz = 0
        elif xyz_error_bias[1] == np.inf:
            self.px = 0
            self.py = self.error_rate
            self.pz = 0
        elif xyz_error_bias[2] == np.inf:
            self.px = 0
            self.py = 0
            self.pz = self.error_rate
        else:
            self.px, self.py, self.pz = (
                self.error_rate * xyz_error_bias / np.sum(xyz_error_bias)
            )

        if self.hadamard_rotate == 0:
            self.channel_probs_x = np.ones(self.N) * (self.px)
            self.channel_probs_z = np.ones(self.N) * (self.pz)
            self.channel_probs_y = np.ones(self.N) * (self.py)

        elif self.hadamard_rotate == 1:
            n1 = self.hadamard_rotate_sector1_length
            self.channel_probs_x = np.hstack(
                [np.ones(n1) * (self.px), np.ones(self.N - n1) * (self.pz)]
            )
            self.channel_probs_z = np.hstack(
                [np.ones(n1) * (self.pz), np.ones(self.N - n1) * (self.px)]
            )
            self.channel_probs_y = np.ones(self.N) * (self.py)
        else:
            raise ValueError(
                f"The hadamard rotate attribute should be set to 0 or 1. Not '{self.hadamard_rotate}"
            )

        self.channel_probs_x.setflags(write=False)
        self.channel_probs_y.setflags(write=False)
        self.channel_probs_z.setflags(write=False)

    def _decoder_setup(self):
        """
        Setup for the BP+OSD decoders
        """

        self.ms_scaling_factor = float(self.ms_scaling_factor)

        # decoder for Z errors
        self.bpd_z = BpOsdDecoder(
            self.hx,
            channel_probs=self.channel_probs_z + self.channel_probs_y,
            max_iter=self.max_iter,
            bp_method=self.bp_method,
            ms_scaling_factor=self.ms_scaling_factor,
            osd_method=self.osd_method,
            osd_order=self.osd_order,
        )

        # decoder for X-errors
        self.bpd_x = BpOsdDecoder(
            self.hz,
            channel_probs=self.channel_probs_x + self.channel_probs_y,
            max_iter=self.max_iter,
            bp_method=self.bp_method,
            ms_scaling_factor=self.ms_scaling_factor,
            osd_method=self.osd_method,
            osd_order=self.osd_order,
        )

    def _generate_error(self):
        """
        Generates a random error on both the X and Z components of the code
        distributed according to the channel probability vectors.
        """

        for i in range(self.N):
            rand = np.random.random()
            if rand < self.channel_probs_z[i]:
                self.error_z[i] = 1
                self.error_x[i] = 0
            elif (
                self.channel_probs_z[i]
                <= rand
                < (self.channel_probs_z[i] + self.channel_probs_x[i])
            ):
                self.error_z[i] = 0
                self.error_x[i] = 1
            elif (
                (self.channel_probs_z[i] + self.channel_probs_x[i])
                <= rand
                < (
                    self.channel_probs_x[i]
                    + self.channel_probs_y[i]
                    + self.channel_probs_z[i]
                )
            ):
                self.error_z[i] = 1
                self.error_x[i] = 1
            else:
                self.error_z[i] = 0
                self.error_x[i] = 0

        return self.error_x, self.error_z

    def run_decode_sim(self):
        """
        This function contains the main simulation loop and controls the output.
        """

        # save start date
        self.start_date = datetime.datetime.fromtimestamp(time.time()).strftime(
            "%A, %B %d, %Y %H:%M:%S"
        )

        pbar = tqdm(
            range(self.run_count + 1, self.target_runs + 1),
            disable=self.tqdm_disable,
            ncols=0,
        )

        start_time = time.time()
        save_time = start_time

        for self.run_count in pbar:
            self._single_run()

            pbar.set_description(
                f"d_max: {self.min_logical_weight}; OSDW_WER: {self.osdw_word_error_rate*100:.3g}±{self.osdw_word_error_rate_eb*100:.2g}%; OSDW: {self.osdw_logical_error_rate*100:.3g}±{self.osdw_logical_error_rate_eb*100:.2g}%; OSD0: {self.osd0_logical_error_rate*100:.3g}±{self.osd0_logical_error_rate_eb*100:.2g}%;"
            )

            current_time = time.time()
            save_loop = current_time - save_time

            if (
                int(save_loop) > self.save_interval
                or self.run_count == self.target_runs
            ):
                save_time = time.time()
                self.runtime = save_loop + self.runtime

                self.runtime_readable = time.strftime(
                    "%H:%M:%S", time.gmtime(self.runtime)
                )

                if self.output_file is not None:
                    f = open(self.output_file, "w+")
                    print(self.output_dict(), file=f)
                    f.close()

                if (
                    self.osdw_logical_error_rate_eb > 0
                    and self.osdw_logical_error_rate_eb / self.osdw_logical_error_rate
                    < self.error_bar_precision_cutoff
                ):
                    print(
                        "\nTarget error bar precision reached. Stopping simulation..."
                    )
                    break

        return json.dumps(self.output_dict(), sort_keys=True, indent=4)

    def output_dict(self):
        """
        Function for formatting the output
        """

        output_dict = {}
        for key, value in self.__dict__.items():
            if key in self.output_keys:
                output_dict[key] = value
        # return output_dict
        return json.dumps(output_dict, sort_keys=True, indent=4)
