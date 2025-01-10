import numpy as np
from ldpc import mod2
from tqdm import tqdm
import scipy


def gf2_to_gf4(bin):
    n = int(len(bin) / 2)
    gf4 = np.zeros(n).astype(int)
    for i in range(n):
        if bin[i] == 1 and bin[i + n] == 0:
            gf4[i] = 1
        elif bin[i] == 0 and bin[i + n] == 1:
            gf4[i] = 3
        elif bin[i] == 1 and bin[i + n] == 1:
            gf4[i] = 2
        else:
            gf4[i] = 0
    return gf4


class stab_code:
    def __init__(self, hx, hz, name=None):
        if not scipy.sparse.issparse(hx) and not scipy.sparse.issparse(hz):
            hx = scipy.sparse.csr_matrix(hx).astype(np.uint8)
            hz = scipy.sparse.csr_matrix(hz).astype(np.uint8)

        if name is None:
            self.name = "<Unamed stabiliser code>"
        else:
            self.name = name

        self.hx = hx.astype(np.uint8)
        self.hz = hz.astype(np.uint8)
        self.init_code()

        self.h = scipy.sparse.hstack([self.hx, self.hz])
        self.l = scipy.sparse.hstack([self.lx, self.lz])

    def init_code(self):
        self.h = scipy.sparse.hstack([self.hx, self.hz])
        self.N = self.hx.shape[1]
        self.K = self.N - mod2.rank(self.h)
        self.compute_logical_operators()
        self.D = np.nan

    def compute_logical_operators(self):
        # compute logical operators
        # Kernel H

        ker_H = mod2.kernel(scipy.sparse.hstack([self.hz, self.hx]))

        rankH = mod2.rank(self.h)

        log_stack = scipy.sparse.vstack([self.h, ker_H])
        pivots = mod2.pivot_rows(log_stack)[rankH:]
        self.l = log_stack[pivots]
        self.lx = self.l[:, 0 : self.N]
        self.lz = self.l[:, self.N : 2 * self.N]

        self.K = int(self.l.shape[0] / 2)

    def compute_code_distance(self, return_logicals=False):
        if self.N > 10:
            print(
                "Warning: computing a code distance of codes with N>10 will take a long time."
            )

        re, r, _, _ = mod2.row_echelon(self.h)
        stab_basis = re[0:r]
        logical_stack = scipy.sparse.vstack([stab_basis, self.l])
        all_logicals = mod2.row_span(logical_stack)[1:].toarray()
        # np.argmin(np.sum(all_logicals,axis=1))

        d_min = self.N
        min_indices = []
        min_logicals = []
        for i in tqdm(range(len(all_logicals))):
            logical = all_logicals[i]
            logical = gf2_to_gf4(logical)
            temp = np.count_nonzero(logical)
            if temp < d_min:
                d_min = temp
                min_indices = [i]
                min_logicals = [logical]
            elif temp == d_min:
                min_indices.append(i)
                min_logicals.append(logical)

        # d_min=np.min( np.sum(all_logicals,axis=1) )
        self.D = d_min

        # print(all_logicals)

        if return_logicals:
            return np.array(min_logicals)

        return d_min

    def test(self, show_tests=True):
        valid_code = True

        # if self.K==np.nan: self.compute_dimension()

        code_label = f"{self.code_params}"

        if show_tests:
            print(f"{self.name}, {code_label}")

        try:
            assert self.N == self.hz.shape[1] == self.lz.shape[1] == self.lx.shape[1]
            assert self.K == self.lz.shape[0] // 2 == self.lx.shape[0] // 2
            if show_tests:
                print(" -Block dimensions: Pass")
        except AssertionError:
            valid_code = False
            print(" -Block dimensions incorrect")

        try:
            print(type(self.hx))
            print(type(self.hz))
            hx_hz = (self.hz @ self.hx.T) + (self.hx @ self.hz.T)
            hx_hz = hx_hz.data % 2
            assert not np.any(hx_hz.data)
            if show_tests:
                print(" -PCMs commute hz@hx.T==0: Pass")
        except AssertionError:
            valid_code = False
            print(" -PCMs commute hz@hx.T==0: Fail")

        # if show_tests and valid_code: print("\t-PCMs commute hx@hz.T == hz@hx.T ==0: Pass")

        try:
            assert ((self.hx @ self.lz.T + self.hz @ self.lx.T).data % 2).any() == 0
        except AssertionError:
            valid_code = False
            print(" -lx \in ker{hz} AND lz \in ker{hx}: Fail")

        if show_tests and valid_code:
            print("\t-lx \in ker{hz} AND lz \in ker{hx}: Pass")

        try:
            lx_lz = self.lx @ self.lz.T + self.lz @ self.lx.T
            lx_lz.data = lx_lz.data % 2
            print(type(self.lz))

            print(self.lz.toarray())

            assert mod2.rank(lx_lz) == self.l.shape[0]
            if show_tests:
                print(" -lx and lz anticommute: Pass")
        except AssertionError:
            valid_code = False
            print(" -lx and lz anitcommute: Fail")

        # if show_tests and valid_code: print("\t- lx and lz anitcommute: Pass")

        if show_tests and valid_code:
            print(f"{self.name} is a valid stabiliser code w/ params {code_label}")

        return valid_code

    @property
    def code_params(self):
        return f"[[{self.N},{self.K},{self.D}]]"
