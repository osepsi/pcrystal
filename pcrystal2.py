import numpy as np
from scipy.spatial import distance
import matplotlib.pyplot as plt
import math
import os
import fcrystal


# from IPython.core.debugger import Tracer

class PCrystal:

    def __init__(self, Nx, Ny, Nz, Lx, Ly, Lz, dt, tmax,
                 gr2_max=3, L_addx=0, L_addy=0, L_addz=0,
                 img_folder="."):

        self.Lx = Lx
        self.Ly = Ly
        self.Lz = Lz

        self.L_addx = L_addx
        self.L_addy = L_addy
        self.L_addz = L_addz

        self.V_tot = (Lx + L_addx) * (Ly + L_addy) * (Lz + L_addz)

        self.Nx = Nx
        self.Ny = Ny
        self.Nz = Nz
        self.N3 = Nx * Ny * Nz

        self.dx = float(Lx) / Nx
        self.dy = float(Ly) / Ny
        self.dz = float(Lz) / Nz

        self.t_series = np.arange(0, tmax, dt)
        self.dt = self.t_series[1] - self.t_series[0]
        self.Nt = len(self.t_series)

        if not os.path.isdir(img_folder):
            os.makedirs(img_folder)

        self.img_folder = img_folder
        self.gr2_max = gr2_max

        self.x = np.linspace(-Lx / 2. + self.dx * 0.5, Lx / 2. - self.dx * 0.5, num=Nx)
        self.y = np.linspace(-Ly / 2. + self.dy * 0.5, Ly / 2. - self.dy * 0.5, num=Ny)
        self.z = np.linspace(-Lz / 2. + self.dz * 0.5, Lz / 2. - self.dz * 0.5, num=Nz)

    def _init_run(self):

        self.aaa = np.zeros((self.Nx, self.Ny, self.Nz), dtype=np.uint16)

        self.o = np.zeros((0, 3))
        self.r = np.zeros((0))
        self.V = np.zeros((0), dtype=np.uint16)  # volume of the crystalline phase in units of pixels
        self.Nsp = 0

        self.to_grow = np.array([], dtype=np.uint16)
        self.gr2 = np.array([], dtype=np.uint8)
        self.k_current = np.array([0])

    def grow(self, n_init, G, F, random_starting_radius=0, savefig=0, plot=0):

        self._init_run()

        plt.ioff()

        dr = G * self.dt
        self._add_n_new(n_init, dr, random_starting_radius)

        for it, t in enumerate(self.t_series):
            self._grow_all_1step(dr)

            n_new_sph = np.random.poisson(F * self.dt * (1 - self.k_current[-1]) * self.V_tot)
            self._add_n_new(n_new_sph, dr, 1)

            if savefig:
                plt.figure()
                plt.imshow(self.aaa[:, :, self.Nz / 2], cmap=plt.get_cmap('Paired'))
                plt.savefig(self.img_folder + "/nov" + str(it) + ".png", bbox_inches='tight')
                plt.close()

            print('-----------------------------------')
            print(it, "of", self.Nt, "| k = ", self.k_current[-1] * 100, "%", "| growing", len(self.to_grow), "of",
                  self.Nsp)
            if self.k_current[-1] >= 0.8 and len(self.to_grow) == 0:
                print("volume full.")
                break

        # delete sphs on edge???

        plt.ion()

        if plot:
            rad = (self.V * self.dx * self.dy * self.dz / np.pi * 0.75) ** (1. / 3.)
            plt.figure()
            plt.hist(rad)
            plt.ylabel('freq')
            plt.xlabel('equiv. radius')

            plt.figure()
            # plt.plot(t_series,V_tbl,label ='illesztes')
            plt.plot(self.t_series[0:it + 1], self.k_current[:-1], label='cdf')
            plt.xlabel('t')
            plt.title('cdf')

        return self.t_series[0:it + 1], self.k_current[:-1], self.V * self.dx * self.dy * self.dz

    # def spontaneous_growth(self, G, F, savefig = 0, plot = 0):
    #
    #     self._init_run()
    #
    #     plt.ioff()
    #
    #     dr = G*self.dt
    #     for it, t in enumerate(self.t_series):
    #         self._grow_all_1step(dr)
    #
    #         n_new_sph = np.random.poisson(F*self.dt*(1-self.k_current[-1])*self.V_tot)
    #         self._add_n_new(n_new_sph, dr, 1)
    #
    #         if savefig:
    #             plt.figure()
    #             plt.imshow(self.aaa[:,:,self.Nz/2],cmap=plt.get_cmap('Paired'))
    #             plt.savefig(self.img_folder+"/nov"+str(it)+".png",bbox_inches='tight')
    #             plt.close()
    #
    #         print("-----------------------------------")
    #         print(it, "of", self.Nt, "| k = ", self.k_current[-1]*100, "%", "| growing", len(self.to_grow), "of", self.Nsp)
    #         if self.k_current[-1] >= 0.8 and len(self.to_grow) == 0:
    #             print("volume full.")
    #             break
    #
    #     # delete sphs on edge???
    #
    #     plt.ion()
    #
    #     if plot:
    #         rad = (self.V*self.dx*self.dy*self.dz/np.pi*0.75)**(1./3.)
    #         plt.figure()
    #         plt.hist(rad)
    #         plt.ylabel('freq')
    #         plt.xlabel('equiv. radius')
    #
    #         plt.figure()
    #         # plt.plot(t_series,V_tbl,label ='illesztes')
    #         plt.plot(self.t_series[0:it+1], self.k_current[:-1], label = 'cdf')
    #
    #     return (self.t_series[0:it+1],self.k_current[:-1],self.V*self.dx*self.dy*self.dz)

    def _add_n_new(self, nsp, r0, random_starting_radius=0):
        for isp in range(nsp):
            self._add_new(r0, random_starting_radius)

    def _add_new(self, r0, random_starting_radius):

        # get new radius
        if random_starting_radius:
            r_start = np.random.rand() * r0
        else:
            r_start = r0

        # get origin
        do_it_again = 1
        while do_it_again:
            new_o = (np.random.rand(1, 3) - 0.5) * np.array(
                [[self.Lx + self.L_addx, self.Ly + self.L_addy, self.Lz + self.L_addz]])
            do_it_again = np.any(np.sum((self.o - new_o) ** 2, axis=1) < (self.r + r_start) ** 2)

        self.o = np.append(self.o, new_o, axis=0)
        self.r = np.append(self.r, 0)
        self.V = np.append(self.V, 0)
        self.Nsp = self.Nsp + 1
        self.to_grow = np.append(self.to_grow, self.Nsp - 1)
        self.gr2 = np.append(self.gr2, self.gr2_max)
        self._grow_isp(self.Nsp - 1, r_start)

    def _grow_isp(self, i_sp, dr):

        self.r[i_sp] = self.r[i_sp] + dr

        # if self.r[i_sp] < self.dx:
        #     num_r = self.dx
        # else:
        #     num_r = self.r[i_sp]

        num_r = self.r[i_sp]

        x_min_ind = max(0, int(math.ceil((-num_r - self.x[0] + self.o[i_sp][0]) / self.dx)))
        x_max_ind = min(int(math.floor((num_r - self.x[0] + self.o[i_sp][0]) / self.dx)), self.Nx - 1)

        y_min_ind = max(0, int(math.ceil((-num_r - self.y[0] + self.o[i_sp][1]) / self.dy)))
        y_max_ind = min(int(math.floor((num_r - self.y[0] + self.o[i_sp][1]) / self.dy)), self.Ny - 1)

        z_min_ind = max(0, int(math.ceil((-num_r - self.z[0] + self.o[i_sp][2]) / self.dz)))
        z_max_ind = min(int(math.floor((num_r - self.z[0] + self.o[i_sp][2]) / self.dz)), self.Nz - 1)

        # find enclosing cube indicies
        # x_tmp_ind = np.nonzero(np.logical_and(self.x>=self.o[i_sp][0]-num_r, self.x<=self.o[i_sp][0]+num_r))
        # y_tmp_ind = np.nonzero(np.logical_and(self.y>=self.o[i_sp][1]-num_r, self.y<=self.o[i_sp][1]+num_r))
        # z_tmp_ind = np.nonzero(np.logical_and(self.z>=self.o[i_sp][2]-num_r, self.z<=self.o[i_sp][2]+num_r))

        # if(x_max_ind-x_tmp_ind[0][-1] != 0):
        #     print x_tmp_ind
        #     print x_max_ind

        # if np.any(x_tmp_ind) and np.any(y_tmp_ind) and np.any(z_tmp_ind):
        if x_min_ind <= x_max_ind:
            sph2 = ((self.x[x_min_ind:x_max_ind + 1, np.newaxis, np.newaxis] - self.o[i_sp][0]) ** 2 +
                    (self.y[np.newaxis, y_min_ind:y_max_ind + 1, np.newaxis] - self.o[i_sp][1]) ** 2 +
                    (self.z[np.newaxis, np.newaxis, z_min_ind:z_max_ind + 1] - self.o[i_sp][2]) ** 2)
            # sph = sph2 <= self.r[i_sp]**2
            sph = np.logical_and(sph2 <= self.r[i_sp] ** 2, sph2 >= (self.r[i_sp] - dr) ** 2)

            parta = self.aaa[x_min_ind:x_max_ind + 1, y_min_ind:y_max_ind + 1, z_min_ind:z_max_ind + 1][sph]

            zind = parta == 0
            parta[zind] = i_sp + 1

            self.aaa[x_min_ind:x_max_ind + 1, y_min_ind:y_max_ind + 1, z_min_ind:z_max_ind + 1][sph] = parta

            dV = np.sum(zind)

            ###########################################
            # parta = self.aaa[x_min_ind:x_max_ind+1, y_min_ind:y_max_ind+1, z_min_ind:z_max_ind+1][sph]
            # inva = np.logical_not(parta)

            # self.aaa[x_min_ind:x_max_ind+1, y_min_ind:y_max_ind+1, z_min_ind:z_max_ind+1][sph] = inva.astype(np.uint16)*(i_sp+1)+parta
            # dV = np.sum(inva)

            ############################################
            # self.aaa[x_min_ind:x_max_ind+1, y_min_ind:y_max_ind+1, z_min_ind:z_max_ind+1][sph] += inva.astype(np.uint16)*(i_sp+1)

            ############################################
            # tmp2 = np.nonzero(np.logical_and(np.logical_not(self.aaa[x_min_ind:x_max_ind+1,
            #                                                          y_min_ind:y_max_ind+1,
            #                                                          z_min_ind:z_max_ind+1]),
            #                                  sph))
            # tmp2 = (tmp2[0]+x_min_ind,
            #         tmp2[1]+y_min_ind,
            #         tmp2[2]+z_min_ind)

            ## tmp3 = np.argmin((np.array(o)[:,0,np.newaxis]-x[tmp2[0]][np.newaxis,:])**2+
            ##                  (np.array(o)[:,1,np.newaxis]-y[tmp2[1]][np.newaxis,:])**2+
            ##                  (np.array(o)[:,2,np.newaxis]-z[tmp2[2]][np.newaxis,:])**2,axis=0) == i_sp

            ## aaa[(tmp2[0][tmp3],tmp2[1][tmp3],tmp2[2][tmp3])] = i_sp+1

            # self.aaa[tmp2] = i_sp+1
            # dV = len(self.aaa[tmp2])

            self.V[i_sp] = self.V[i_sp] + dV

            if dV == 0:
                self.gr2[i_sp] = self.gr2[i_sp] - 1
                if self.gr2[i_sp] == 0:
                    # print "finishing", i_sp
                    self.to_grow = np.delete(self.to_grow, np.nonzero(self.to_grow == i_sp)[0][0])
            else:
                self.gr2[i_sp] = self.gr2_max

    def _grow_all_1step(self, dr):

        for i_sp in np.random.permutation(self.to_grow):  # kicsit megbolonditjuk a sorrendet :)
            self._grow_isp(i_sp, dr)
        self.k_current = np.append(self.k_current, np.sum(self.V) / float(self.N3))


def mc_cdf(N_seeds, LL, n_rnd_pos, G, mode="i", F=0, n_rnd_loop=1):
    """
    Calculates the conversion curve of instantaneous nucleation with calculating at random test positions the crystallization time.
    Growth speed of G. Initial number of seeds is N_seeds.
    Length of computational cube is LL. Number of random test positions is n_rnd_pos.
    
    typ: "i" instantaneous
         "s" spontaneous
         
    If you want to increase timespan, increase number of seeds (N_seeds) or decrease LL
    """

    os = np.random.rand(N_seeds, 3) * LL

    if mode == "i":
        t0 = np.zeros(N_seeds)
    elif mode == "s":
        tmax = N_seeds / F / LL ** 3
        t0 = np.random.rand(N_seeds) * tmax
        # check whether seeds will be born

        d0 = distance.cdist(os, os)
        np.fill_diagonal(d0, np.inf)  # eliminate problem of self killing
        t_kill = np.min(d0 / G + t0[:, np.newaxis], axis=0)  # the first time a seed would be eaten if not yet born
        os = os[t_kill > t0]
        t0 = t0[t_kill > t0]

    ts = []
    for i_loop in range(n_rnd_loop):
        # get n random positions
        test_pos = np.random.rand(n_rnd_pos, 3) * LL

        # find distances between points and centers

        dd = distance.cdist(os, test_pos)

        ts.append(np.min(dd / G + t0[:, np.newaxis], axis=0))

    ts = np.array(ts).flatten()

    if mode == "i":
        t_sorted = np.sort(ts)
    elif mode == "s":
        t_sorted = np.sort(ts[ts <= tmax])

    print(len(t0), "particles started growing")

    return (t_sorted, np.arange(len(t_sorted)) / (n_rnd_pos * n_rnd_loop - 1.))


def fcrystal_init(rho_inst, nucl_rate, growth_rate, n_test0=1000, n_sphs=10000, verbose=False):
    """

    :param rho_inst:
    :param nucl_rate:
    :param growth_rate:
    :param n_test0:
    :param n_sphs:
    :param verbose:
    :return: sph_centers, sph_birth_times, cube_edge_length, tmax, n_inst
    """

    n_thrm0 = 0
    n_inst = n_sphs - n_thrm0
    if rho_inst == 0:
        comp_vol = 1.
    else:
        comp_vol = n_inst / rho_inst  # computational volume
    LL = comp_vol ** (1. / 3.)
    rho_thrm0 = 0.01 / comp_vol  # N_thrm / comp_vol, with N_thrm << 1
    tmax = rho_thrm0 / nucl_rate

    if verbose:
        print("running initial simulation...")
        print("n_thrm0   n_inst    LL        tmax")
    while 1:
        # adjust tmax
        #        LL
        #        n_inst
        #
        if verbose:
            print("{:<10.4g}{:<10.4g}{:<10.4g}{:<10.4g}".format(n_thrm0, n_inst, LL, tmax))

        t0s = np.concatenate((np.zeros(n_inst), np.sort(np.random.rand(n_thrm0) * tmax)))
        centers = np.random.rand(n_thrm0 + n_inst, 3) * LL

        (new_n, new_c, new_t0) = fcrystal.sp_seeds(centers, t0s, growth_rate, 0)
        # contains really born instant and thermal seeds
        # N_thrm_real = new_n - n_inst
        # rho_thrm_real = N_thrm_real/comp_vol
        # print(new_n, "born of", N_thrm+n_inst, "seeds")

        (t_conv_arr, vols_arr) = fcrystal.sr_cdf(new_c[0:new_n], new_t0[0:new_n], LL, growth_rate, n_test0)
        ts3 = np.sort(t_conv_arr[t_conv_arr < tmax], axis=None)
        ks3 = np.arange(len(ts3)) / (n_test0 - 1.)

        # if ks3.size == 0:
        #     print(0, N_thrm)
        # else:
        #     print(ks3[-1], N_thrm)

        if ks3.size == 0:
            tmax = tmax * 10
        elif ks3[-1] < 0.1:
            tmax = tmax * 5
        elif ks3[-1] < 0.5:
            tmax = tmax * 3
        elif ks3[-1] < 0.99:
            tmax = tmax * 2
        else:
            break

        # defining variables for the next iteration
        rho_thrm0 = tmax * nucl_rate
        comp_vol = n_sphs / (rho_inst + rho_thrm0)
        LL = comp_vol ** (1. / 3.)
        n_thrm0 = int(round(comp_vol * rho_thrm0))
        n_inst = n_sphs - n_thrm0

    return new_c[0:new_n], new_t0[0:new_n], LL, tmax, n_inst


def fcrystal_main(c_sphs_arr, t0_sphs_arr, LL, G, tmax, n_test=1000000):
    """

    :param c_sphs_arr: center of spherulites, list of 3length lists
    :param t0_sphs_arr: birth time of nuclei, list of floats
    :param LL: edge length of computational cube
    :param G: growth rate
    :param tmax: max time of simulation
    :param n_test: number of test points
    :return: ts, ks, vols
    """

    (t_conv_arr, vols_arr) = fcrystal.sr_cdf(c_sphs_arr, t0_sphs_arr, LL, G, n_test)
    ts3 = np.sort(t_conv_arr[t_conv_arr < tmax], axis=None)
    ks3 = np.arange(len(ts3)) / (n_test - 1.)

    return ts3, ks3, vols_arr


def fcrystal_cdf(rho_inst, F, G, n_test=1000000, n_sphs=10000, n_test0=1000, verbose=False):
    """
    atu (arbitrary time unit)
    alu (arbitrary length unit)


    :param rho_inst: density of instantly born sphs
    :param F: birth rate in 1/atu/alu^3
    :param G: growth rate in alu/atu
    :param n_test: number of test points
    :param n_sphs: number of spherulites
    :param n_test0: number of test points for initial simulation
    :param verbose:
    :return: (t [array], k [array], spherulite_volumes [array], N_thrm, N_inst, total_volume)
    """

    c_sphs_arr, t0_sphs_arr, LL, tmax, N_inst = fcrystal_init(rho_inst, F, G, n_test0, n_sphs, verbose)
    ts3, ks3, vols_arr = fcrystal_main(c_sphs_arr, t0_sphs_arr, LL, G, tmax, n_test)

    nn = len(c_sphs_arr)
    N_thrm = nn - N_inst

    if verbose:
        # print("kmax: {:.5g}".format(ks3[-1]))
        print("f_inst    f_therm")
        print("{:<10.2%}{:<10.2%}".format(float(N_inst) / nn, float(N_thrm) / nn))

    return ts3, ks3, vols_arr, N_thrm, N_inst, LL**3
