import  numpy as np
cimport numpy as np
import scipy.linalg as spl
import pyross.utils
from libc.math cimport exp, pow, sqrt
cimport cython
import warnings
from types import ModuleType
import os


DTYPE   = np.float
ctypedef np.float_t DTYPE_t
@cython.wraparound(False)
@cython.boundscheck(False)
@cython.cdivision(True)
@cython.nonecheck(False)
cdef class ContactMatrixFunction:
    """Generates a time dependent contact matrix

    Parameters
    ----------
    CH: 2d np.array
        Contact matrix at home
    CW: 2d np.array
        Contact matrix at work
    CS: 2d np.array
        Contact matrix at school
    CO: 2d np.array
        Contact matrix at other locations
    """

    cdef:
        np.ndarray CH, CW, CS, CO
        double aW, aS, aO
        Protocol protocol

    def __init__(self, CH, CW, CS, CO):
        self.CH, self.CW, self.CS, self.CO = CH, CW, CS, CO

    def contactMatrix(self, t, **kwargs):
        cdef np.ndarray[DTYPE_t, ndim=2] C
        aW, aS, aO = self.protocol(t, **kwargs)
        C = self.CH + aW*self.CW + aS*self.CS + aO*self.CO
        return C

    cpdef get_individual_contactMatrices(self):
        return self.CH, self.CW, self.CS, self.CO

    def constant_contactMatrix(self, aW=1, aS=1, aO=1):
        '''Constant contact matrix

        Parameters
        ----------
        aW: float, optional
            Scale factor for work contacts.
        aS: float, optional
            Scale factor for school contacts.
        aO: float, optional
            Scale factor for other contacts.

        Returns
        -------
        contactMatrix: callable
            A function that takes t as an argument and outputs the contact matrix
        '''

        self.protocol=ConstantProtocol(aW, aS, aO)
        return self.contactMatrix

    def interventions_temporal(self,times,interventions):
        '''Temporal interventions

        Parameters
        ----------
        time: np.array
            Ordered array with temporal boundaries between the different interventions.
        interventions: np.array
            Ordered matrix with prefactors of CW, CS, CO matrices
            during the different time intervals.
            Note that len(interventions) = len(times) + 1

        Returns
        -------
        contactMatrix: callable
            A function that takes t as an argument and outputs the contact matrix
        '''

        self.protocol = TemporalProtocol(np.array(times), np.array(interventions))
        return self.contactMatrix

    def interventions_threshold(self,thresholds,interventions):
        '''Temporal interventions

        Parameters
        ----------
        threshold: np.array
            Ordered array with temporal boundaries between the different interventions.
        interventions: np.array
            Array of shape [K+1,3] with prefactors during different phases of intervention
            The current state of the intervention is defined by
            the largest integer "index" such that state[j] >= thresholds[index,j] for all j.

        Returns
        -------
        contactMatrix: callable
            A function that takes t as an argument and outputs the contact matrix
        '''

        self.protocol = ThresholdProtocol(np.array(thresholds), np.array(interventions))
        return self.contactMatrix

    def intervention_custom_temporal(self, intervention_func, **kwargs):
        '''Custom temporal interventions

        Parameters
        ----------
        intervention_func: callable
            The calling signature is `intervention_func(t, **kwargs)`,
            where t is time and kwargs are other keyword arguments for the function.
            The function must return (aW, aS, aO).
        kwargs: dict
            Keyword arguments for the function.

        Returns
        -------
        contactMatrix: callable
            A function that takes t as an argument and outputs the contact matrix.

        Examples
        --------
        An example for an custom temporal intervetion that
        allows for some anticipation and reaction time

        >>> def fun(t, width=1, loc=0) # using keyword arguments for parameters of the intervention
                a = (1-np.tanh((t-loc)/width))/2
                return a, a, a
        >>> contactMatrix = generator.intervention_custom_temporal(fun, width=5, loc=10)
        '''

        self.protocol = CustomTemporalProtocol(intervention_func, **kwargs)
        return self.contactMatrix




cdef class Protocol:

    def __init__(self):
        pass

    def __call__(self, t):
        pass

cdef class ConstantProtocol(Protocol):
    cdef:
        double aW, aS, aO

    def __init__(self, aW=1, aS=1, aO=1):
        self.aW = aW
        self.aS = aS
        self.aO = aO

    def __call__(self, double t):
        return self.aW, self.aS, self.aO

cdef class TemporalProtocol(Protocol):
    cdef:
        np.ndarray times, interventions

    def __init__(self, np.ndarray times, np.ndarray interventions):
        self.times = times
        self.interventions = interventions

    def __call__(self, double t):
        cdef:
            Py_ssize_t index
            np.ndarray t_arr=self.times, prefac_arr=self.interventions
        index = np.argmin( t_arr < t)
        if index == 0:
            if t >= t_arr[len(t_arr)-1]:
                index = -1
        return prefac_arr[index,0], prefac_arr[index,1], prefac_arr[index,2]




cdef class ThresholdProtocol(Protocol):
    cdef:
        np.ndarray thresholds, interventions

    def __init__(self, np.ndarray thresholds, np.ndarray interventions):
        self.thresholds = thresholds
        self.interventions = interventions

    def __call__(self, double t, S=None, Ia=None, Is=None):
        cdef:
            np.ndarray[DTYPE_t, ndim=1] state
            np.ndarray thresholds=self.thresholds, prefac_arr=self.interventions
            Py_ssize_t index
        state = np.concatenate((S, Ia, Is))
        index = np.argmin((thresholds <= state ).all(axis=1))
        if index == 0:
            N = len(thresholds)
            if (thresholds[N-1] <= state ).all():
                index = N
        return prefac_arr[index,0], prefac_arr[index,1], prefac_arr[index,2]




cdef class CustomTemporalProtocol(Protocol):
    cdef:
        object intervention_func
        dict kwargs

    def __init__(self, intervention_func, **kwargs):
        self.intervention_func = intervention_func
        self.kwargs = kwargs

    def __call__(self, double t):
        return self.intervention_func(t, **self.kwargs)




cdef class SpatialContactMatrix:
    '''A class for generating a spatial compartmental model with commute data

    Let :math:`\\mu, \\nu` denote spatial index and i, j denote age group index.

    .. math::
        C^{\\mu\\nu}_{ij} = \\frac{1}{N^\\mu_i} \widetilde{C}^{\\mu \\nu}_{ij}

    Parameters
    ----------
    b: float
        Parameter b in the above equation
    populations: np.array(n_loc, M)
        Populations of regions by age groups. Here n_loc is the number of regions
        and M is the number of age groups.
    areas: np.array(n_loc)
        Areas of the geographical regions.
    commutes: np.array(n_loc, n_loc, M)
        Each entry commute[mu, nu, i] needs to be the number of people of age
        group i commuting from :math:`\\mu` to :math:`\\nu`.
        Entries with :math:`\\mu = \\nu` are ignored.
    '''
    cdef:
        readonly np.ndarray Ni, pops, commute_fraction,
        readonly np.ndarray density_factor, norm_factor, local_contacts
        readonly Py_ssize_t n_loc, M
        readonly double work_ratio

    def __init__(self, double b, double work_ratio, np.ndarray populations,
                    np.ndarray areas, np.ndarray commutes):
        cdef:
            double [:, :] densities
        self.n_loc = populations.shape[0]
        if self.n_loc != commutes.shape[0] or self.n_loc !=commutes.shape[1]:
            raise Exception('The first two dimensions of `commutes` must be \
                             the number of regions')
        if self.n_loc != areas.shape[0] :
            raise Exception('The first dimension of populations \
                             areas must be equal')

        self.M = populations.shape[1]
        if self.M != commutes.shape[2]:
            raise Exception('The last dimension of `commutes` must be \
                             the number of age groups')

        self.Ni = np.sum(populations, axis=0)
        self.work_ratio = work_ratio

        # for pops, etc, pops[0] gives the values without commutors
        # and pops[1] gives the values with commutors
        self.pops = np.empty((2, self.n_loc, self.M), dtype=DTYPE)
        self.commute_fraction = np.empty((self.n_loc, self.n_loc, self.M), dtype=DTYPE)
        self.density_factor = np.empty((2, self.n_loc, self.M, self.M), dtype=DTYPE)
        self.norm_factor = np.zeros((2, self.M, self.M), dtype=DTYPE)
        self.local_contacts = np.empty((2, self.n_loc, self.M, self.M), dtype=DTYPE)

        self.pops[0] = populations
        self._process_commute(populations.astype('float'), commutes.astype('float'))
        rho = np.sum(populations)/np.sum(areas)
        self._compute_density_factor(b, rho, areas.astype('float'))

    def spatial_contact_matrix(self, np.ndarray CM):
        self._compute_local_contacts(CM.astype('float'))
        cdef:
            Py_ssize_t mu, nu, i, j, M=self.M, n_loc=self.n_loc
            double [:, :, :] f=self.commute_fraction,
            double [:, :, :] C=self.local_contacts[0], CC=self.local_contacts[1]
            double [:, :] pop=self.pops[0], commute_time_pop=self.pops[1]
            np.ndarray spatial_CM
            double p, cc, work_ratio=self.work_ratio
        spatial_CM = np.zeros((n_loc, M, n_loc, M))
        for i in range(M):
            for j in range(M):
                for mu in range(n_loc):
                    p = pop[mu, i]
                    spatial_CM[mu, i, mu, j] = C[mu, i, j]*(1-work_ratio)
                    for nu in range(n_loc):
                        for eta in range(n_loc):
                            cc = CC[eta, i, j]
                            f1 = f[mu, eta, i]
                            f2 = f[nu, eta, j]
                            spatial_CM[mu, i, nu, j] += f1*f2*cc*work_ratio
                        spatial_CM[mu, i, nu, j] /= p
        return spatial_CM

    cdef _process_commute(self, double [:, :] pop, double [:, :, :] commutes):
        cdef:
            Py_ssize_t i, mu, nu
            double influx, outflux, p, f
            double [:, :, :] pops=self.pops
            double [:, :, :] commute_fraction=self.commute_fraction
        for i in range(self.M):
            for mu in range(self.n_loc):
                influx = np.sum(commutes[:, mu, i])
                outflux = np.sum(commutes[mu, :, i])
                p = pops[0, mu, i] + influx - outflux
                pops[1, mu, i] = p
                commute_fraction[mu, mu, i] = 1
                for nu in range(self.n_loc):
                    if nu != mu:
                        f = commutes[nu, mu, i]/p
                        commute_fraction[nu, mu, i] = f
                        commute_fraction[mu, mu, i] -= f


    cdef _compute_density_factor(self, double b, double rho, double [:] areas):
        cdef:
            Py_ssize_t mu, i, j, a
            double rhoi, rhoj
            double [:, :, :] pops=self.pops
            double [:, :, :, :] density_factor=self.density_factor
            double [:, :, :] norm_factor=self.norm_factor
        for mu in range(self.n_loc):
            for i in range(self.M):
                for j in range(self.M):
                    for a in range(2):
                        rhoi = pops[a, mu, i]/areas[mu]
                        rhoj = pops[a, mu, j]/areas[mu]
                        density_factor[a, mu, i, j] = pow(rhoi*rhoj/rho**2, b)
                        norm_factor[a, i, j] += density_factor[a, mu, i, j]


    cdef _compute_local_contacts(self, double [:, :] CM):
        cdef:
            Py_ssize_t mu, i, j, M=self.M, n_loc=self.n_loc
            double [:] Ni = self.Ni
            double [:, :, :, :] local_contacts=self.local_contacts
            double [:, :, :] norm=self.norm_factor
            double [:, :, :, :] density_factor=self.density_factor
            double c, d
        for mu in range(n_loc):
            for i in range(M):
                for j in range(M):
                    c = CM[i, j] * Ni[i]
                    for a in range(2):
                        d = density_factor[a, mu, i, j]
                        local_contacts[a, mu, i, j] = c * d / norm[a, i, j]




cdef class MinimalSpatialContactMatrix:
    '''A class for generating a minimal spatial compartmental model

    Let :math:`\\mu, \\nu` denote spatial index and i, j denote age group index,

    .. math::
        C^{\\mu\\nu}_{ij} = \\frac{1}{a_{ij}} k^{\\mu\\nu} g^{\\mu\\nu}_{ij} f^{\\mu\\nu}_{ij} C_{ij}

    where
    .. math::
        a_{ij} = \sum_{\mu\nu} k^{\mu\nu}_{ij} g^{\mu\nu}_{ij} \sqrt{ \frac{ N_i^\mu N_j^\nu} {N_i N_j}}
        k^{\mu \nu} = \exp ( - b_{\mu\nu} | r^\mu - r^\nu |) \\
        b_{\mu, \nv} = (\rho^\mu \rho^\nu)^{-c} \\
        g^{\mu\nu}_{ij} = (\rho^\mu_i \rho^\nu_j)^b \\
        f^{\mu\nu}_{ij} = \sqrt{ \frac{N_i N_j^\nu}{N_i^\mu N_j}}

    where :math:`C_{ij}` is the non-spatial age structured contact matrix,
    :math:`\rho^\mu_i` is the population density of age group i
    in spatial location :math:`\mu`,
    :math:`\rho_i` is the total population density of age group i,
    :math:`r^\mu` is the spatial position of the :math:`\mu`th location,
    :math:`N_i^\mu` is the population of the age group i at :math:`\mu`,
    :math:`N_i` is the total population of the age group i,
    and b, c are free parameters

    Parameters
    ----------
    b: float
        Parameter b in the above equation.
    c: float
        Paramter c in the above equation.
    populations: np.array(n_loc, M)
        The population of each age group (total number of age group = M)
        in each location (total number of location = n_loc ).
    areas: np.array(n_loc)
        The area of each geographical region.
    coordinates: np.array(n_loc, 2)
        The GPS coordinates of each geographical region.
    '''

    cdef:
        Py_ssize_t n_loc, M
        readonly np.ndarray Ni, spatial_kernel, density_factor, rescale_factor, normalisation_factor

    def __init__(self, double b, double c, np.ndarray populations,
                        np.ndarray areas, np.ndarray coordinates):
        cdef:
            double [:, :] densities
        self.n_loc = populations.shape[0]
        if self.n_loc != areas.shape[0] or self.n_loc != coordinates.shape[0]:
            raise Exception('The first dimension of populations \
                             areas and coordinates must be equal')
        self.M = populations.shape[1]
        self.Ni = np.sum(populations, axis=0)
        self.spatial_kernel = np.empty((self.n_loc, self.n_loc))
        self.density_factor = np.empty((self.n_loc, self.M, self.n_loc, self.M))
        self.rescale_factor = np.empty((self.n_loc, self.M, self.n_loc, self.M))
        self.normalisation_factor = np.zeros((self.M, self.M))

        densities = populations/areas[:, np.newaxis]
        self._compute_spatial_kernel(densities, coordinates, c)
        self._compute_density_factor(densities, b)
        self._compute_rescale_factor(populations)
        self._compute_normalisation_factor(populations)

    def spatial_contact_matrix(self, double [:, :] CM):
        cdef:
            Py_ssize_t n_loc=self.n_loc, M=self.M, i, j, m, n
            double [:] Ni=self.Ni
            np.ndarray spatial_CM=np.empty((n_loc, M, n_loc, M),
                                                        dtype=DTYPE)
            double [:, :] a=self.normalisation_factor, k=self.spatial_kernel
            double [:, :, :, :] g=self.density_factor, f=self.rescale_factor

        for i in range(n_loc):
            for j in range(n_loc):
                for m in range(M):
                    for n in range(M):
                        spatial_CM[i, m, j, n]=k[i, j]*g[i, m, j, n]*f[i, m, j, n]/a[m, n]*CM[m, n]
        return spatial_CM

    cdef _compute_spatial_kernel(self, double [:, :] densities,
                                        double [:, :] coordinates, double c):
        cdef:
            Py_ssize_t n_loc=self.n_loc, i, j
            double [:] densities_sum_over_age
            double [:, :] spatial_kernel=self.spatial_kernel
            double d, k, rhoi, rhoj
        densities_sum_over_age = np.sum(densities, axis=1)
        for i in range(n_loc):
            spatial_kernel[i, i] = 1
            rho_i = densities_sum_over_age[i]
            for j in range(i+1, n_loc):
                rho_j = densities_sum_over_age[j]
                d = pyross.utils.distance_on_Earth(coordinates[i], coordinates[j])
                k = exp(-d*pow(rho_i*rho_j, -c))
                spatial_kernel[i, j] = k
                spatial_kernel[j, i] = k

    cdef _compute_density_factor(self, double [:, :] densities, double b):
        cdef:
            Py_ssize_t n_loc=self.n_loc, M=self.M, i, j, m, n
            double [:, :, :, :] density_factor=self.density_factor
            double rho_im, rho_jn, f
        for i in range(n_loc):
            for m in range(M):
                rho_im = densities[i, m]
                for j in range(n_loc):
                    for n in range(M):
                        rho_jn = densities[j, n]
                        f = pow(rho_im*rho_jn, b)
                        density_factor[i, m, j, n] = f

    cdef _compute_rescale_factor(self, double [:, :] populations):
        cdef:
            Py_ssize_t n_loc=self.n_loc, M=self.M, i, j, m, n
            double [:] Ni=self.Ni
            double [:, :, :, :] rescale_factor=self.rescale_factor
        for i in range(n_loc):
            for j in range(n_loc):
                for m in range(M):
                    for n in range(M):
                        rescale_factor[i, m, j, n] = sqrt(Ni[m]/Ni[n]*populations[j, n]/populations[i, m])

    cdef _compute_normalisation_factor(self, double [:, :] populations):
        cdef:
            Py_ssize_t n_loc=self.n_loc, M=self.M, i, j, m, n
            double [:] Ni=self.Ni
            double [:, :] norm_fac=self.normalisation_factor, k=self.spatial_kernel
            double [:, :, :, :] g=self.density_factor, f=self.rescale_factor
            double a
        for m in range(M):
            for n in range(M):
                for i in range(n_loc):
                    for j in range(n_loc):
                        a = sqrt(populations[i, m]*populations[j, n]/(Ni[m]*Ni[n]))
                        norm_fac[m, n] += k[i, j]*g[i, m, j, n]*a

@cython.wraparound(False)
@cython.boundscheck(False)
@cython.cdivision(True)
@cython.nonecheck(False)
cdef class SIR(ContactMatrixFunction):


    def basicReproductiveRatio(self, data, state='constant'):
        C = self.CH + self.CW + self.CS + self.CO
        # matrix for linearised dynamics
        alpha = data['alpha']                         # infection rate
        beta  = data['beta']                         # infection rate
        gIa   = data['gIa']                          # removal rate of Ia
        gIs   = data['gIs']                          # removal rate of Is
        fsa   = data['fsa']                          # the self-isolation parameters
        M     = data['M']
        Ni    = data['Ni']

        L0 = np.zeros((M, M))
        L  = np.zeros((2*M, 2*M))

        if state=='constant':
            for i in range(M):
                for j in range(M):
                    L0[i,j]=C[i,j]*Ni[i]/Ni[j]

            L[0:M, 0:M]     = alpha*beta/gIa*L0
            L[0:M, M:2*M]   = fsa*alpha*beta/gIs*L0
            L[M:2*M, 0:M]   = ((1-alpha)*beta/gIa)*L0
            L[M:2*M, M:2*M] = fsa*((1-alpha)*beta/gIs)*L0
            r0 = np.max(np.linalg.eigvals(L))

        else:
            t = data.get('t');
            Nt = t.size
            r0 = np.zeros((Nt))

            for tt in range(Nt):
                S = np.array((data['X'][tt,0:M]))
                for i in range(M):
                    for j in range(M):
                        L0[i,j]=C[i,j]*S[i]/Ni[j]

                L[0:M, 0:M]     = alpha*beta/gIa*L0
                L[0:M, M:2*M]   = fsa*alpha*beta/gIs*L0
                L[M:2*M, 0:M]   = ((1-alpha)*beta/gIa)*L0
                L[M:2*M, M:2*M] = fsa*((1-alpha)*beta/gIs)*L0
                r0[tt] = np.real(np.max(np.linalg.eigvals(L)))

        return r0





"""
---
KreissPy
https://gitlab.com/AustenBolitho/kreisspy
sublibrary dedicated to calculating the transient effects of non-normal
contact matricies
---
"""
def _epsilon_eval(z, A, ord=2):
    """
    Finds the value of \epsilon for a given complex number and matrix.
    Uses the first definition of the pseudospectrum in Trfethen & Embree
    ord="svd" uses fourth definition (may be faster)


    inputs:
    z: length 2 array representing a complex number
    A: an MxM matrix
    order: order of the matrix norm given from associated vector norm
    default is regular L2 norm -> returns maximum singular value.
    accepted inputs are any in spl.norm or "svd"
    """
    z=np.array(z)
    A=np.array(A)
    zc = complex(z[0], z[1])
    try :
        ep = 1/spl.norm(spl.inv(zc*np.eye(*A.shape)-A),ord=ord)
        # ep = spl.norm(zc*np.eye(*A.shape)-A,ord=ord)
    except TypeError:
        if ord=="svd":
            ep = np.min(spl.svdvals(zc*np.eye(*A.shape)-A))
        else: raise Exception("invalid method")
    return ep


def _inv_epsilon_eval(z, A, ord=2):
    """
    Finds the value of 1/\epsilon for a given complex number and matrix.
    Uses the first definition of the pseudospectrum in Trfethen & Embree
    ord="svd" uses fourth definition (may be faster)


    inputs:
    z: length 2 array representing a complex number
    A: an MxM matrix
    order: order of the matrix norm given from associated vector norm
    default is regular L2 norm -> returns maximum singular value.
    accepted inputs are any in spl.norm or "svd"
    """
    z=np.array(z)
    A=np.array(A)
    zc = complex(z[0], z[1])
    try :
        iep = spl.norm(spl.inv(zc*np.eye(*A.shape)-A),ord=ord)
    except TypeError:
        if ord=="svd":
            iep = 1/np.min(spl.svdvals(zc*np.eye(*A.shape)-A))
        else: raise Exception("invalid method")
    return iep


def _kreiss_eval(z, A, theta=0, ord=2):
    """
    Kreiss constant guess for a matrix and pseudo-eigenvalue.

    inputs:
    z: length 2 array representing a complex number
    A: an MxM matrix
    theta: normalizing factor found in Townley et al 2007, default 0
    ord: default 2, order of matrix norm
    """
    z=np.array(z)
    A=np.array(A)
    kg = (z[0]-theta)*_inv_epsilon_eval(z, A, ord=ord)
    return kg


def _inv_kreiss_eval(z, A, theta=0, ord=2):
    """
    1/Kreiss constant guess for a matrix and pseudo-eigenvalue.
    for minimizer

    inputs:
    z: length 2 array representing a complex number
    A: an MxM matrix
    theta: normalizing factor found in Townley et al 2007, default 0
    ord: default 2, order of matrix norm
    """
    z=np.array(z)
    A=np.array(A)
    ikg = _epsilon_eval(z, A, ord=ord)/np.real(z[0]-theta) if z[0]-theta > 0 else np.inf
    # print(z[0]-theta)
    return ikg


def _transient_properties(guess, A, theta=0, ord=2):
    """
    returns the maximal eigenvalue (spectral abcissa),
    initial groth rate (numerical abcissa),
    the Kreiss constant (minimum bound of transient)
    and time of transient growth

    inputs:
    A: an MxM matrix
    guess: initial guess for the minimizer
    theta: normalizing factor found in Townley et al 2007, default 0
    ord: default 2, order of matrix norm

    returns: [spectral abcissa, numerical abcissa, Kreiss constant ,
              duration of transient, henrici's departure from normalcy']
    """
    from scipy.optimize import minimize
    A = np.array(A)
    if np.array_equal(A@A.T, A.T@A):
        warnings.warn("The input matrix is normal")
        # print("The input matrix is normal")
    evals = spl.eigvals(A)
    sa = evals[
        np.where(np.real(evals) == np.amax(np.real(evals)))[0]
    ]
    na = np.real(np.max(spl.eigvals((A+A.T)/2)))
    m = minimize(_inv_kreiss_eval, guess, args=(A, theta, ord),
                 bounds=((0, None), (None, None)))
    K = 1/m.fun
    tau = np.log(1/m.fun)/m.x[0]
    evals2 = np.dot(evals,np.conj(evals))
    frobNorm = spl.norm(A,ord='fro')
    henrici = np.sqrt(frobNorm**2-evals2)#/frobNorm
    return np.array([sa, na, K, tau, henrici],dtype=np.complex64)


def _first_estimate( A, tol=0.001):
    """
    Takes the eigenvalue with the largest real part

    returns a first guess of the
    maximal pseudoeigenvalue in the complex plane
    """
    evals = spl.eigvals(A)
    revals = np.real(evals)
    idxs = np.where(revals == np.amax(revals))[0]
    mevals = evals[idxs]
    iguesses = []
    for evl in mevals:
        guess = []
        a, b = np.real(evl), np.imag(evl)
        if a > 0:
            guess = [a+tol, b]
        else:
            guess = [tol, b]
        iguesses.append(guess)
    return iguesses


def characterise_transient(A, tol=0.001, theta=0, ord=2):
    """
    The maximal eigenvalue (spectral abcissa),
    initial groth rate (numerical abcissa),
    the Kreiss constant (minimum bound of transient)
    and time of transient growth

    Parameters 
    -----------
    A    : an MxM matrix
    tol  : Used to find a first estimate of the pseudospectrum
    theta: normalizing factor found in Townley et al 2007, default 0
    ord  : default 2, order of matrix norm

    Returns 
    ---------
    [spectral abcissa, numerical abcissa, Kreiss constant, 
    duration of transient, henrici's departure from normalcy]

    """
    guesses = _first_estimate(A, tol)
    transient_properties = [1, 0, 0]
    for guess in guesses:
        tp = _transient_properties(guess, A, theta, ord)
        if tp[2] > transient_properties[2]:
            transient_properties = tp
        else:
            pass
    return transient_properties









'''
---
Contact structure: Projecting social contact matrices in 152 countries using
contact surveys and demographic data,
Kiesha Prem, Alex R. Cook, Mark Jit, PLOS Computational Biology, (2017)

---
Below we provide the contact matrix for some of the countries using
data from above
'''


def getCM(country='India', sheet=1):
    """
    Method to compute contact matrices of a given country
    
    The data is read from sheets at: 

    https://github.com/rajeshrinet/pyross/tree/master/examples/data/contact_matrices_152_countries

    Parameters
    ----------
    country: string 
        Default is India
    sheet: int 
        Default is 1
        sheet takes value 1 and 2

    Returns
    ----------
    four np.arrays: CH, CW, CS, CO of the given country 

    CH - home, 

    CW - work, 

    CS - school, 

    CO - other locations 

    """

    u1 ='https://raw.githubusercontent.com/rajeshrinet/pyross/master/examples'
    u2 ='/data/contact_matrices_152_countries/MUestimates_'

    if sheet==1:
        uH = u1 + u2 + 'home_1.xlsx'
        uW = u1 + u2 + 'work_1.xlsx'
        uS = u1 + u2 + 'school_1.xlsx'
        uO = u1 + u2 + 'other_locations_1.xlsx'
    elif sheet==2:
        uH = u1 + u2 + 'home_2.xlsx'
        uW = u1 + u2 + 'work_2.xlsx'
        uS = u1 + u2 + 'school_2.xlsx'
        uO = u1 + u2 + 'other_locations_2.xlsx'
    else:
        raise Exception('There are only two sheets, choose 1 or 2')

    import pandas as pd
    CH = np.array(pd.read_excel(uH,  sheet_name=country))
    CW = np.array(pd.read_excel(uW,  sheet_name=country))
    CS = np.array(pd.read_excel(uS,  sheet_name=country))
    CO = np.array(pd.read_excel(uO,  sheet_name=country))
    return CH, CW, CS, CO

