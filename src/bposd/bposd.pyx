#cython: language_level=3, boundscheck=False, wraparound=False, initializedcheck=False, cdivision=True

import numpy as np
from scipy.special import comb as nCr

cdef class bposd_decoder:

    '''
    A class implementing a belief propagation + ordered statistics
    decoding for LDPC codes

    Parameters
    ----------
    mat: numpy.ndarray
        The parity check matrix of the binary code in numpy.ndarray format.
    error_rate: float64, optional
        The bit error rate.
    max_iter: int, optional
        The maximum number of iterations for the BP decoder. If max_iter==0, the BP algorithm
        will iterate n times, where n is the block length of the code.
    bp_method: str OR int, optional
        The BP method. Currently three methods are implemented: 1) "ps": product sum updates;
        2) "ms": min-sum updates; 3) "msl": min-sum log updates
    ms_scaling_factor: float64, optional
        Sets the min-sum scaling factor for the min-sum BP method
    channel_probs: list, optional
        This parameter can be used to set the initial error channel across all bits.
    osd_order: str or int, optional
        Sets the OSD order.
    osd_method: str or int, optional
        The OSD method. Currently three methods are availbe: 1) "osd_0": Zero-oder OSD; 2) "osd_e": exhaustive OSD;
        3) "osd_cs": combination-sweep OSD.

    '''


    def __cinit__(self,mat, error_rate=None, max_iter=0, bp_method=0,ms_scaling_factor=0.625,channel_probs=[None],osd_order=-1, osd_method=1):

        self.MEM_ALLOCATED=False

        cdef i,j

        self.m=mat.shape[0]
        self.n=mat.shape[1]

        #setup BP decoder
        #nb. couldn't get cython class inheritance to work so implemented it manually.
        self.bpd=bp_decoder(mat,error_rate=error_rate,\
        max_iter=max_iter,\
        ms_scaling_factor=ms_scaling_factor,\
        channel_probs=channel_probs, bp_method=bp_method)
        self.H=self.bpd.H
        self.max_iter=self.bpd.max_iter
        self.bp_method=self.bpd.bp_method
        self.ms_scaling_factor=self.bpd.ms_scaling_factor
        assert self.n==self.H.n_cols #validate number of bits in mod2sparse format
        assert self.m==self.H.n_rows #validate number of checks in mod2sparse format
        self.error=self.bpd.error
        self.synd=self.bpd.synd
        self.bp_decoding_synd=self.bpd.bp_decoding_synd
        self.bp_decoding=self.bpd.bp_decoding
        self.channel_probs=self.bpd.channel_probs
        self.log_prob_ratios=self.bpd.log_prob_ratios

        #memory allocation for OSD specific attributes
        self.osd0_decoding=<char*>calloc(self.n,sizeof(char)) #the OSD_0 decoding
        self.osdw_decoding=<char*>calloc(self.n,sizeof(char)) #the osd_w decoding

        #OSD setup

        # OSD method
        if str(osd_method).lower() in ['OSD_0','osd_0','0','osd0']:
            osd_method=0
            osd_order=0
        elif str(osd_method).lower() in ['osd_e','1','osde','exhaustive','e']:
            osd_method=1
            if osd_order>15:
                print("WARNING: Running the 'OSD_E' (Exhaustive method) with search depth greater than 15 is not recommended. Use the 'osd_cs' method instead.")
        elif str(osd_method).lower() in ['osd_cs','2','osdcs','combination_sweep','combination_sweep','cs']:
            osd_method=2
        else:
            raise ValueError(f"ERROR: OSD method '{osd_method}' invalid. Please choose from the following methods: 'OSD_0', 'OSD_E' or 'OSD_CS'.")

        self.osd_order=int(osd_order)
        self.osd_method=int(osd_method)

        self.encoding_input_count=0
        
        if self.osd_order>-1:
            self.rank=mod2sparse_rank(self.H)
            try:
                assert self.osd_order<=(self.n - self.rank)
            except AssertionError:
                self.osd_order=-1
                raise ValueError(f"For this code, the OSD order should be set in the range 0<=osd_oder<={self.n - self.rank}.")
            self.cols=<int*>calloc(self.n,sizeof(int)) 
            self.orig_cols=<int*>calloc(self.n,sizeof(int))
            self.rows=<int*>calloc(self.m,sizeof(int))
            self.k=self.n-self.rank

        if self.osd_order>0:
            self.y=<char*>calloc(self.n,sizeof(char))
            self.g=<char*>calloc(self.m,sizeof(char))
            self.Htx=<char*>calloc(self.m,sizeof(char))
            self.Ht_cols=<int*>calloc(self.k,sizeof(int)) 

        if osd_order==0: pass
        elif self.osd_order>0 and self.osd_method==1: self.osd_e_setup()
        elif self.osd_order>0 and self.osd_method==2: self.osd_cs_setup()
        elif self.osd_order==-1: pass
        else: raise Exception(f"ERROR: OSD method '{osd_method}' invalid")

        self.MEM_ALLOCATED=True

    cdef void osd_e_setup(self):

        self.encoding_input_count=2**self.osd_order
        self.osdw_encoding_inputs=<char**>calloc(self.encoding_input_count,sizeof(char*))
        for i in range(self.encoding_input_count):
            self.osdw_encoding_inputs[i] = decimal_to_binary_reverse(i, self.n - self.rank)

    cdef void osd_cs_setup(self):

        cdef int kset_size=self.n-self.rank

        assert self.osd_order<=kset_size


        self.encoding_input_count=kset_size+nCr(self.osd_order,2)
        
        self.osdw_encoding_inputs=<char**>calloc(self.encoding_input_count,sizeof(char*))
        cdef int total_count=0
        for i in range(kset_size):
            self.osdw_encoding_inputs[total_count] = <char*>calloc(kset_size,sizeof(char))
            self.osdw_encoding_inputs[total_count][i]=1
            total_count+=1

        for i in range(self.osd_order):
            for j in range(self.osd_order):
                if i<j:
                    self.osdw_encoding_inputs[total_count] = <char*>calloc(kset_size,sizeof(char))
                    self.osdw_encoding_inputs[total_count][i]=1
                    self.osdw_encoding_inputs[total_count][j]=1
                    total_count+=1

        assert total_count==self.encoding_input_count


    cdef char* decode_cy(self, char* syndrome):

        self.bpd.synd=syndrome
        self.bpd.bp_decode_cy()

        # print(double2numpy(self.bpd.log_prob_ratios,self.n))


        if self.osd_order==-1: return self.bp_decoding

        #if BP has converged, return the BP solution
        if self.bpd.converge==1:
            for j in range(self.n): self.osd0_decoding[j]=self.osdw_decoding[j]=self.bp_decoding[j]
            return self.osd0_decoding

        #if BP doesn't converge, run OSD post-processing
        self.osd()

        if self.osd_order==0:
            for j in range(self.n): self.osdw_decoding[j]=self.osd0_decoding[j]
            return self.osd0_decoding
        else:
            return self.osdw_decoding

    cpdef np.ndarray[np.int_t, ndim=1] decode(self, np.ndarray[np.int_t, ndim=1] syndrome):
        """
        Runs the BP+OSD decoder for a given syndrome.

        Parameters
        ----------
        syndrome: numpy.ndarray
            The syndrome to be decoded.

        Returns
        -------
        numpy.ndarray
            The BP+OSD decoding in numpy.ndarray format.

        """
        self.synd=numpy2char(syndrome,self.synd)
        self.decode_cy(self.synd)
        if self.osd_order==-1: return char2numpy(self.bp_decoding,self.n)
        else: return char2numpy(self.osdw_decoding,self.n)

    #OSD Post-processing
    cdef int osd(self):
        cdef int i, j
        cdef long int l
        cdef mod2sparse *L
        cdef mod2sparse *U

        #allocating L and U matrices 
        L=mod2sparse_allocate(self.m,self.rank)
        U=mod2sparse_allocate(self.rank,self.n)

        #sort the columns on the basis of the soft decisions
        soft_decision_col_sort(self.log_prob_ratios,self.cols, self.n)

        #save the original sorted column order
        for i in range(self.n):
            self.orig_cols[i]=self.cols[i]

        #find the LU decomposition of the ordered matrix
        mod2sparse_decomp_osd(
            self.H,
            self.rank,
            L,
            U,
            self.rows,
            self.cols)


        #solve the syndrome equation with most probable full-rank submatrix
        LU_forward_backward_solve(
            L,
            U,
            self.rows,
            self.cols,
            self.synd,
            self.osd0_decoding)


        if self.osd_order==0:
            mod2sparse_free(U)
            mod2sparse_free(L)
            return 1


        #return the columns outside of the information set to their orginal ordering (the LU decomp scrambles them)
        cdef int check, counter, in_pivot
        cdef mod2sparse* Ht=mod2sparse_allocate(self.m,self.k)

        counter=0

        for i in range(self.n):
            check=self.orig_cols[i]
            in_pivot=0
            for j in range(self.rank):
                if self.cols[j]==check:
                    in_pivot=1
                    break
            
            if in_pivot==0:
                self.cols[counter+self.rank]=check
                counter+=1

        #create the HT matrix
        for i in range(self.k):
            self.Ht_cols[i]=self.cols[i+self.rank]

        mod2sparse_copycols(self.H,Ht,self.Ht_cols)



        # cdef osd_0_weight=bin_char_weight(self.osd0_decoding,self.n)

        cdef double osd_min_weight=0
        for i in range(self.n):
            # osd_min_weight+=self.osd0_decoding[i]*log(1/self.channel_probs[i])
            # osd_min_weight+=self.osd0_decoding[i]
            if self.osd0_decoding[i]==1:
                osd_min_weight+=log(1/self.channel_probs[i])


        for i in range(self.n):
            self.osdw_decoding[i]=self.osd0_decoding[i]

        cdef double solution_weight
        cdef char *x



        for l in range(self.encoding_input_count):
            x=self.osdw_encoding_inputs[l]
            mod2sparse_mulvec(Ht,x,self.Htx)
            for i in range(self.m):
                self.g[i]=self.synd[i]^self.Htx[i]

            LU_forward_backward_solve(
                L,
                U,
                self.rows,
                self.cols,
                self.g,
                self.y)

            for i in range(self.k):
                self.y[self.Ht_cols[i]]=x[i]

            solution_weight=0.0
            for i in range(self.n):
                # solution_weight+=self.y[i]*log(1/self.channel_probs[i])
                # solution_weight+=self.y[i]
                if self.y[i]==1:
                    solution_weight+=log(1/self.channel_probs[i])

            if solution_weight<osd_min_weight:
                osd_min_weight=solution_weight
                for i in range(self.n):
                    self.osdw_decoding[i]=self.y[i]

        mod2sparse_free(Ht)
        mod2sparse_free(U)
        mod2sparse_free(L)
        return 1

    def update_channel_probs(self,channel):
        """
        Function updates the channel probabilities for each bit in the BP decoder.

        Parameters
        ----------
        channel: numpy.ndarray
            A list of the channel probabilities for each bit

        Returns
        -------
        NoneType
        """
        cdef j
        for j in range(self.n): self.channel_probs[j]=channel[j]  

    @property
    def bp_method(self):
        """
        Getter for the BP method
        
        Returns
        -------
        str
        """
        if self.bp_method==0: return "product_sum"
        elif self.bp_method==1: return "minimum_sum"
        elif self.bp_method==2: return "product_sum_log"
        elif self.bp_method==3: return "minimum_sum_log"

    @property
    def iter(self):
        """
        Getter. Returns the number of iterations in the last round of BP decoding.
        
        Returns
        -------
        numpy.ndarray
        """
        return self.bpd.iter

    @property
    def ms_scaling_factor(self):
        """
        Getter. Returns the min-sum scaling factor.

        Returns
        -------
        float64
        """
        return self.ms_scaling_factor

    @property
    def max_iter(self):
        """
        Getter. Returns the maximum interation depth for the BP decoder.
        
        Returns
        -------
        int
        """
        return self.max_iter

    @property
    def converge(self):
        """
        Getter. Returns `1' if the last round of BP succeeded (converged) and `0' if it failed.
        
        Returns
        -------
        int
        """
        return self.bpd.converge

    @property
    def bp_decoding(self):
        """
        Getter. Returns the recovery vector from the last round of BP decoding.
        
        Returns
        -------
        numpy.ndarray
        """
        return char2numpy(self.bp_decoding,self.n)

    @property
    def log_prob_ratios(self):
        """
        Getter. Returns the soft-decision log probability ratios from the last round of BP
        decoding.
        
        Returns
        -------
        numpy.ndarray
        """
        return double2numpy(self.log_prob_ratios,self.n)

    @property
    def osd_method(self):
        """
        Getter. Returns the OSD method.
        
        Returns
        -------
        str
        """
        if self.osd_order==-1: return None
        if self.osd_method==0: return "osd_0"
        if self.osd_method==1: return "osd_e"
        if self.osd_method==2: return "osd_cs"
    
    @property
    def osd_order(self):
        """
        Getter. Returns the OSD order.
        
        Returns
        -------
        int
        """
        return self.osd_order

    @property
    def osdw_decoding(self):
        """
        Getter. Returns the recovery vector from the last round of BP+OSDW decoding.
        
        Returns
        -------
        numpy.ndarray
        """        
        return char2numpy(self.osdw_decoding,self.n)

    @property
    def osd0_decoding(self):
        """
        Getter. Returns the recovery vector from the last round of BP+OSD0 decoding.
        
        Returns
        -------
        numpy.ndarray
        """        
        return char2numpy(self.osd0_decoding,self.n)


    def __dealloc__(self):
        
        if self.MEM_ALLOCATED==True:
    
            free(self.osd0_decoding)
            free(self.osdw_decoding)

            if self.osd_order>-1:
                free(self.cols)
                free(self.rows)
                free(self.orig_cols)

            if self.osd_order>0:
                free(self.Htx)
                free(self.g)
                free(self.y)
                free(self.Ht_cols)

            if self.encoding_input_count!=0:
                for i in range(self.encoding_input_count):
                    free(self.osdw_encoding_inputs[i])

