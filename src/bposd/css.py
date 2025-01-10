import numpy as np
from ldpc import mod2
from bposd import stab
import scipy


class css_code:
    def __init__(
        self,
        hx=np.array([[]]),
        hz=np.array([[]]),
        code_distance=np.nan,
        name="<Unnamed CSS code>",
    ):
        if not scipy.sparse.issparse(hx) and not scipy.sparse.issparse(hz):
            hx = scipy.sparse.csr_matrix(hx).astype(np.uint8)
            hz = scipy.sparse.csr_matrix(hz).astype(np.uint8)

        self.hx = hx.astype(np.uint8)  # hx pcm
        self.hz = hz.astype(np.uint8)  # hz pcm

        self.lx = scipy.sparse.csr_matrix([]).astype(np.uint8)  # x logicals
        self.lz = scipy.sparse.csr_matrix([]).astype(np.uint8)  # z logicals

        self.N = np.nan  # block length
        self.K = np.nan  # code dimension
        self.D = code_distance  # code distance
        self.L = np.nan  # max column weight
        self.Q = np.nan  # max row weight

        _, nx = self.hx.shape
        _, nz = self.hz.shape
        try:
            assert nx == nz
        except AssertionError:
            raise Exception(
                "Error: hx and hz matrices must have equal numbers of columns!"
            )

        if nx != 0:
            self.compute_dimension()
            self.compute_logicals()

        self.name = name

    def compute_dimension(self):
        self.N = self.hx.shape[1]
        assert self.N == self.hz.shape[1], "Code block length (N) inconsistent!"

        self.K = self.N - mod2.rank(self.hx) - mod2.rank(self.hz)
        return self.K

    def to_stab_code(self):
        hx = scipy.sparse.vstack([np.zeros(self.hz.shape, dtype=np.uint8), self.hx])
        hz = scipy.sparse.vstack([self.hz, np.zeros(self.hx.shape, dtype=np.uint8)])
        return stab.stab_code(hx, hz)

    @property
    def h(self):
        hx = scipy.sparse.vstack([np.zeros(self.hz.shape, dtype=np.uint8), self.hx])
        hz = scipy.sparse.vstack([self.hz, np.zeros(self.hx.shape, dtype=np.uint8)])
        return scipy.sparse.hstack([hx, hz])

    @property
    def l(self):
        lx = scipy.sparse.vstack([np.zeros(self.lz.shape, dtype=np.uint8), self.lx])
        lz = scipy.sparse.vstack([self.lz, np.zeros(self.lx.shape, dtype=np.uint8)])
        return scipy.sparse.hstack([lx, lz])

    def compute_code_distance(self):
        temp = self.to_stab_code()
        self.D = temp.compute_code_distance()
        return self.D

    def compute_logicals(self):
        def compute_lz(hx, hz):
            # lz logical operators
            # lz\in ker{hx} AND \notin Im(Hz.T)

            ker_hx = mod2.nullspace(hx)  # compute the kernel basis of hx
            # in the below we row reduce to find vectors in kx that are not in the image of hz.T.
            log_stack = scipy.sparse.vstack([hz, ker_hx])

            rank_hz = mod2.rank(hz)

            pivots = mod2.pivot_rows(log_stack)[rank_hz:]
            log_ops = log_stack[pivots]
            return log_ops

        if self.K == np.nan:
            self.compute_dimension()
        self.lx = compute_lz(self.hz, self.hx)
        self.lz = compute_lz(self.hx, self.hz)

        return self.lx, self.lz

    @property
    def code_params(self):
        try:
            self.N
        except AttributeError:
            self.N = np.nan
        try:
            self.K
        except AttributeError:
            self.K = np.nan
        try:
            self.D
        except AttributeError:
            self.D = np.nan
        try:
            self.L
        except AttributeError:
            self.L = np.nan
        try:
            self.Q
        except AttributeError:
            self.Q = np.nan

        return f"({self.L},{self.Q})-[[{self.N},{self.K},{self.D}]]"

    def test(self, show_tests=True):
        valid_code = True

        if self.K == np.nan:
            self.compute_dimension()

        if show_tests:
            print(f"{self.name}")

        try:
            assert self.N == self.hz.shape[1] == self.lz.shape[1] == self.lx.shape[1]
            assert self.K == self.lz.shape[0] == self.lx.shape[0]
            if show_tests:
                print(" -Block dimensions: Pass")
        except AssertionError:
            valid_code = False
            print(" -Block dimensions incorrect")

        try:
            assert not np.any((self.hz @ self.hx.T).data % 2)
            if show_tests:
                print(" -PCMs commute hz@hx.T==0: Pass")
        except AssertionError:
            valid_code = False
            print(" -PCMs commute hz@hx.T==0: Fail")

        try:
            assert not np.any((self.hx @ self.hz.T).data % 2)
            if show_tests:
                print(" -PCMs commute hx@hz.T==0: Pass")
        except AssertionError:
            valid_code = False
            print(" -PCMs commute hx@hz.T==0: Fail")

        # if show_tests and valid_code: print("\t-PCMs commute hx@hz.T == hz@hx.T ==0: Pass")

        try:
            assert not np.any((self.hz @ self.lx.T).data % 2)
        except AssertionError:
            valid_code = False
            print(" -lx \in ker{hz} AND lz \in ker{hx}: Fail")

        try:
            assert not np.any((self.hx @ self.lz.T).data % 2)
            if show_tests:
                print(" -lx \in ker{hz} AND lz \in ker{hx}: Pass")
        except AssertionError:
            valid_code = False
            print(" -lx \in ker{hz} AND lz \in ker{hx}: Fail")

        # if show_tests and valid_code: print("\t-lx \in ker{hz} AND lz \in ker{hx}: Pass")

        try:
            lx_lz = self.lx @ self.lz.T
            lx_lz.data = lx_lz.data % 2
            assert mod2.rank(lx_lz) == self.K
            if show_tests:
                print(" -lx and lz anticommute: Pass")
        except AssertionError:
            valid_code = False
            print(" -lx and lz anticommute: Fail")

        # if show_tests and valid_code: print("\t- lx and lz anitcommute: Pass")

        if show_tests and valid_code:
            print(
                f" -{self.name} is a valid CSS code w/ params [{self.N},{self.K},{self.D}]"
            )

        return valid_code
