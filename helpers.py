# helper functions
import numpy as np
import pandas as pd
rng = np.random.default_rng()
import pandas as pd
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from lmfit import minimize, Parameters

######################################################################################
# FUNCTIONS FOR PLOTTING AND INVESTIGATING PARAMETERS
######################################################################################

def plot_ts4(sim, title = ""):
    Its = list_plot(sim[::,(0,1)], plotjoined=True, color="grey", legend_label="I")
    Pts = list_plot(sim[::,(0,2)], plotjoined=True, color="blue", legend_label="P")
    Dts = list_plot(sim[::,(0,3)], plotjoined=True, color="red", legend_label="D")
    Ats = list_plot(sim[::,(0,4)], plotjoined=True, color="purple", legend_label="A", title = title)
    return Its+Pts+Dts+Ats

## ~~ ORIGINAL MODEL FROM HANNAH ~~ ##
# Define Variables and parameters


def lupus_flare_fixed(init = (0.1, 0.4, 1.7, 0.1),
                     tmax = 500,
                     kid = 1,
                     kip = 0.025, #rate of immune complex removal from system
                     kpi = 0.13, #rate of mediator activation and recruitment
                     kpp = 0.02, #immune response amplified by existing inflammatory response (kpi)
                     kpd = 0.001, #rate of activation for pro-inflammatory agents as a result of cytokine release or induced by damaged tissue
                     mup = 0.06, #decay of pro-inflammatory mediators
                     kdip = 0.025, #rate of phagocytosis of immune complexes by immune cells
                     kdp = 0.27, #rate at which collateral damage is produced by pro-inflammatory mediators 
                     mud = 0.04, #decay rate of damage
                     kap = 0.022, #intrarenal production of anti-inflammatory mediators
                     kad = 0.22, #intrarenal rate of tissue damage ??
                     mua = 2.2, #rate of anti-inflammatory agent degradation
                     Ainf = 0.45):

    t_range = srange(0, tmax, 1)

    def fsid_t(t): #time dependency Sid
        sid = 0.005
        if 0 <= t <= 60: 
            sid = 0.005 # 6m-4m
        elif 60 < t <= 120:
            sid = 0.003 # 4m-2m
        elif 120 < t <= 165:
            sid = 0.015 # 2m-2w
        elif 165 < t <= 180:
            sid = 0.015 # 2w-flare
        elif 180 < t <= 195:
            sid = 0.002 # flare-2w
        elif 195 < t <= 225:
            sid = 0.002 # 2w-6w
        elif 225 < t <= 240:
            sid = 0.002 # 6w-2m
        else:
            sid = 0.012 #2m-4m
        return sid 

    def fsi_t(t): #time dependency Si
        si = 0.002
        if 0 <= t <= 60: 
            si = 0.002
        elif 60 < t <= 120:
            si = 0.001
        elif 120 < t <= 165:
            si = 0.005
        elif 165 < t <= 180:
            si = 0.005
        elif 180 < t <= 195:
            si = 0.001
        elif 195 < t <= 225:
            si = 0.001
        elif 225 < t <= 240:
            si = 0.001
        else:
            si = 0.005
        return si 
    
    def fsa_t(t): #time dependency Sa
        sa = 0.05
        if 0 <= t <= 60: 
            sa = 0.05
        elif 60 < t <= 120:
            sa = 0.1
        elif 120 < t <= 165:
            sa = 0.05
        elif 165 < t <= 180:
            sa = 0.05
        elif 180 < t <= 195:
            sa = 0.3
        elif 195 < t <= 225:
            sa = 0.3
        elif 225 < t <= 240:
            sa = 0.3
        else:
            sa = 0.1
        return sa 

    def systems (IPDA, t, si_t,sid_t,sa_t):
        I, P, D, A = IPDA 
        def f(x):
            return x/((1+A/Ainf)**2) 

        Idot = f(si_t(t)) + f(sid_t(t))*(D**2/(kid**2 + D**2)) - kip*f(P)*I
        Pdot = f(kpi*I + kpp*P) + f(kpd*D) - mup*P
        Ddot = kdip*f(P)*I + kdp*f(P) - mud*D
        Adot = sa_t(t) + f(kap*P + kad*D) - mua*A

        LupusSystem = (Idot, Pdot, Ddot, Adot)
        return LupusSystem

    IPDAsim = odeint(systems, init, t_range, args=(fsi_t, fsid_t, fsa_t))
    IPDAsim = np.insert(IPDAsim, 0, t_range, axis=1)

    return IPDAsim

def investigate_fixed(ax, **kwargs):
    simulation = lupus_flare_fixed(**kwargs)
    title_string = "Original Model"
    for key, value in kwargs.items():
        addition = ", %s = %s" % (key, round(value, 4))
        title_string += addition

    plot_ts4(simulation).matplotlib(sub = ax)
    ax.set_title(title_string)
    ax.set_xlabel("day")
    ax.set_ylabel("state value")
    ax.legend()

def generate_ts_params(tmax = 600, flare_odds = 0.05, flare_length = 5, si_0 = 0.002, sid_0 = 0.015, sa_0 = 0.5): # DEFAULT IS 0.2, changed to: 0.001, 0.07, 0.5, 1
    
    # set some values to keep track of if we are on a flare
    flare_ts = []
    flare_tracker = False
    flare_day = 0
    for i in range(tmax):
        event = rng.binomial(n = 1, p = flare_odds)
        if event == 1:
            flare_tracker = True
        if (flare_tracker) & (flare_day < flare_length):
            flare_ts.append("flare")
            flare_day += 1
        else:
            flare_ts.append("no flare")
        if flare_day == flare_length:
            flare_tracker = False
            flare_day = 0

    by_time_params = pd.DataFrame(srange(0, tmax, 1), columns = ["time"], dtype = "float")

    by_time_params['flare_status'] = np.array(flare_ts)

    by_time_params['si'] = np.where(by_time_params['flare_status'] == "flare",
                                si_0 *5,
                                si_0)
    by_time_params['sid'] = np.where(by_time_params['flare_status'] == "flare",
                                sid_0 *5,
                                sid_0)
    by_time_params['sa'] = np.where(by_time_params['flare_status'] == "flare",
                                sa_0/4,
                                sa_0)
    return by_time_params

def lupus_flare_tsinput(by_time_params = generate_ts_params(tmax = 600).to_dict(orient = "index"),
                         init = (0.1, 0.4, 1.7, 0.1), # tuple: (I0, P0, D0, A0)
                         tmax=550, # number of days to simulate 
                         kid = 1,
                         kip = 0.025, #rate of immune complex removal from system
                         kpi = 0.13, #rate of mediator activation and recruitment
                         kpp = 0.02, #immune response amplified by existing inflammatory response (kpi)
                         kpd = 0.001, #rate of activation for pro-inflammatory agents as a result of cytokine release or induced by damaged tissue
                         mup = 0.06, #decay of pro-inflammatory mediators
                         kdip = 0.025, #rate of phagocytosis of immune complexes by immune cells
                         kdp = 0.27, #rate at which collateral damage is produced by pro-inflammatory mediators 
                         mud = 0.04, #decay rate of damage
                         #sa = 0., #addition of anti-inflammatory drugs
                         kap = 0.022, #intrarenal production of anti-inflammatory mediators
                         kad = 0.22, #intrarenal rate of tissue damage ??
                         mua = 2.2, #rate of anti-inflammatory agent degradation
                         Ainf = 0.45
                        ):

    Lupusvars = list(var("I", "P", "D", "A"))
    t_range = srange(0, tmax, 1)

    def fsid_t(t): #time dependency Sid
        integer = round(t)
        return by_time_params.get(integer).get("sid")

    def fsi_t(t): #time dependency Si
        integer = round(t)
        return by_time_params.get(integer).get("si")

    def fsa_t(t): #time dependency Sa
        integer = round(t)
        return by_time_params.get(integer).get("sa")


    def systems (IPDA, t, si_t,sid_t,sa_t):
        I, P, D, A = IPDA 
        def f(x):
            return x/((1+A/Ainf)**2) 

        Idot = f(si_t(t)) + f(sid_t(t))*(D**2/(kid**2 + D**2)) - kip*f(P)*I
        Pdot = f(kpi*I + kpp*P) + f(kpd*D) - mup*P
        Ddot = kdip*f(P)*I + kdp*f(P) - mud*D
        Adot = sa_t(t) + f(kap*P + kad*D) - mua*A

        LupusSystem = (Idot, Pdot, Ddot, Adot)
        return LupusSystem

    IPDAsim = odeint(systems, init, t_range, args=(fsi_t, fsid_t, fsa_t))
    IPDAsim = np.insert(IPDAsim, 0, t_range, axis=1)

    return IPDAsim

def investigate_stochastic(ax, **kwargs):
    simulation = lupus_flare_tsinput(**kwargs)
    title_string = "Stochastic Model"
    for key, value in kwargs.items():
        addition = ", %s = %s" % (key, round(value, 4))
        title_string += addition

    plot_ts4(simulation).matplotlib(sub = ax)
    ax.set_title(title_string)
    ax.set_xlabel("day")
    ax.set_ylabel("state value")
    ax.legend()



######################################################################################
# FUNCTIONS FOR DATA FITTING
######################################################################################

# this was generated randomly; then saved. could be regenerated
by_time_params = pd.read_csv("data/by_time_params.csv")
by_time_params = by_time_params[['time', 'flare_status', 'si', 'sid', 'sa']].to_dict(orient = "index")


def lupus_system(init, t, paras):
        I, P, D, A = init 

        # get params
        si = paras['si'].value #rate that immune complexes deposit in the kidneys
        sa = paras['sa'].value 
        sid = paras['sid'].value #immune response to accumulation of damaged cells
        kip = paras['kip'].value #rate of immune complex removal from system
        kpi = paras['kpi'].value #rate of mediator activation and recruitment
        kpp = paras['kpp'].value #immune response amplified by existing inflammatory response (kpi)
        kpd = paras['kpd'].value #rate of activation for pro-inflammatory agents as a result of cytokine release or induced by damaged tissue
        mup = paras['mup'].value #decay of pro-inflammatory mediators
        kdip = paras['kdip'].value #rate of phagocytosis of immune complexes by immune cells
        kdp = paras['kdp'].value #rate at which collateral damage is produced by pro-inflammatory mediators 
        mud = paras['mud'].value #decay rate of damage
        kap = paras['kap'].value #intrarenal production of anti-inflammatory mediators
        kad = paras['kad'].value #intrarenal rate of tissue damage ??
        mua = paras['mua'].value #rate of anti-inflammatory agent degradation
        Ainf = 0.45
        def f(x):
            return x/((1+(A/Ainf)**2))
            #return x
        
        Idot = f(si) + f(sid)*(D**2/(1 + D**2)) - kip*f(P)*I
        Pdot = f(kpi*I + kpp*P) + f(kpd*D) - mup*P
        Ddot = kdip*f(P)*I + kdp*f(P) - mud*D
        Adot = sa + f(kap*P + kad*D) - mua*A

        LupusSystem = [Idot, Pdot, Ddot, Adot]
        return LupusSystem

def g(t_range, init, paras):
    """
    Solution to the ODE x'(t) = f(t,x,k) with initial condition x(0) = x0
    """
    IPDAsim = odeint(lupus_system, init, t_range, args=(paras,))
    IPDAsim = np.insert(IPDAsim, 0, t_range, axis=1) # add time values
    return IPDAsim


def residual(paras, t_range, data):

    """
    compute the residual between actual data and fitted data
    """
    
    init = paras['I_0'].value, paras['P_0'].value, paras['D_0'].value, paras['A_0'].value
    model = g(t_range, init, paras)

    # you only have data for damage as function of time
    damage_model = model[:, 3]

    return (damage_model - data).ravel()


def fit_fixed(data, ax):
    t_measured = np.array(data['day'])
    uPCR_measured = np.array(data['uPCR'])
    #plt.figure()
    ax.scatter(t_measured, uPCR_measured, marker='o', color='b', label='measured data', s=5)

    init = (0.1, 0.5, uPCR_measured[0], 0.1)
    #t_range = range(0, t_measured.max() + 10, 1)

    # set parameters including bounds; you can also fix parameters (use ,vary=False)
    params = Parameters()

    params.add('sa', value=1, min=0.0001, max=0.5, vary = False)
    params.add('si', value=0.1, min=0.0001, max=0.5, vary = False)
    params.add('sid', value=0.25, min=0.0001, max=0.5, vary = False)
    
    params.add('I_0', value=0.1, vary = False)
    params.add('P_0', value=0.5, vary = False)
    params.add('D_0', value=uPCR_measured[0], vary = False)
    params.add('A_0', value=0.1, vary = False)

    params.add('kip', value=0.027, min=0.0001, max=0.5)
    params.add('kpi', value=0.01, min=0.0001, max=1.)
    params.add('kpp', value=0.018, min=0.0001, max=2., vary = False) #
    params.add('kpd', value=0.001, min=0.0001, max=5., vary = False)#
    params.add('mup', value=0.66, min=0.0001, max=1., vary = False)
    params.add('kdip', value=0.104, min=0.0001, max=0.5)
    params.add('kdp', value=0.001, min=0.0001, max=1., vary = False)
    params.add('mud', value=0.087, min=0.0001, max=0.1, vary = False)
    params.add('kap', value=0.001, min=0.0001, max=0.1)
    params.add('kad', value=0.001, min=0.0001, max=0.5, vary = False)
    params.add('mua', value=2.2, min=0.0001, max=5., vary = False)#

    # fit model
    result = minimize(residual, params, args=(t_measured, uPCR_measured), method='leastsq')  # leastsq method
    # check results of the fit
    data_fitted = g(t_measured, init, result.params)

    # plot fitted data
    ax.plot(t_measured, data_fitted[:, 3], '-', linewidth=2, color='red', label='fitted data')
    ax.set_title("Fixed Fit: " + data['ID'].to_list()[0])
    ax.legend()
    ax.set_xlim([0, max(t_measured) + 20])
    # display fitted statistics

    #plt.show()

    return (result)



def lupus_system_stoch(init, t, paras):
        I, P, D, A = init 

        # get params
        #si = paras['si'].value #rate of immune complex removal from system
        #sid = paras['sid'].value #rate of immune complex removal from system
        #sa = paras['sa'].value #rate of immune complex removal from system
        
        kip = paras['kip'].value #rate of immune complex removal from system
        kpi = paras['kpi'].value #rate of mediator activation and recruitment
        kpp = paras['kpp'].value #immune response amplified by existing inflammatory response (kpi)
        kpd = paras['kpd'].value #rate of activation for pro-inflammatory agents as a result of cytokine release or induced by damaged tissue
        mup = paras['mup'].value #decay of pro-inflammatory mediators
        kdip = paras['kdip'].value #rate of phagocytosis of immune complexes by immune cells
        kdp = paras['kdp'].value #rate at which collateral damage is produced by pro-inflammatory mediators 
        mud = paras['mud'].value #decay rate of damage
        kap = paras['kap'].value #intrarenal production of anti-inflammatory mediators
        kad = paras['kad'].value #intrarenal rate of tissue damage ??
        mua = paras['mua'].value #rate of anti-inflammatory agent degradation

        # get our parameters which were time dependent
        time_int = round(t)
        si = by_time_params.get(time_int).get("si")
        sid = by_time_params.get(time_int).get("sid")
        sa = by_time_params.get(time_int).get("sa")
        
        def f(x):
            return x/((1+(A/0.45)**2)) 
            #return x
        
        Idot = f(si) + f(sid)*(D**2/(1 + D**2)) - kip*f(P)*I
        Pdot = f(kpi*I + kpp*P) + f(kpd*D) - mup*P
        Ddot = kdip*f(P)*I + kdp*f(P) - mud*D
        Adot = sa + f(kap*P + kad*D) - mua*A

        LupusSystem = [Idot, Pdot, Ddot, Adot]
        return LupusSystem

def g_stoch(t_range, init, paras):
    """
    Solution to the ODE x'(t) = f(t,x,k) with initial condition x(0) = x0
    """
    IPDAsim = odeint(lupus_system_stoch, init, t_range, args=(paras,))
    IPDAsim = np.insert(IPDAsim, 0, t_range, axis=1) # add time values
    return IPDAsim


def residual_stoch(paras, t_range, data):

    """
    compute the residual between actual data and fitted data
    """

    init = paras['I_0'].value, paras['P_0'].value, paras['D_0'].value, paras['A_0'].value
    model = g_stoch(t_range, init, paras)

    # you only have data for damage as function of time
    damage_model = model[:, 3]

    return (damage_model - data).ravel()

def fit_stochastic(data, ax):
    t_measured = np.array(data['day'])
    uPCR_measured = np.array(data['uPCR'])
    #plt.figure()
    ax.scatter(t_measured, uPCR_measured, marker='o', color='b', label='measured data', s=5)

    init = (0.1, 0.5, uPCR_measured[0], 0.1)
    #t_range = range(0, t_measured.max() + 10, 1)


    # set parameters including bounds; you can also fix parameters (use ,vary=False)
    params = Parameters()
    
    params.add('I_0', value=0.1, vary = False)
    params.add('P_0', value=0.5, vary = False)
    params.add('D_0', value=uPCR_measured[0], vary = False)
    params.add('A_0', value=0.1, vary = False)

    params.add('kip', value=0.027, min=0.0001, max=0.5)
    params.add('kpi', value=0.01, min=0.0001, max=1.)
    params.add('kpp', value=0.018, min=0.0001, max=2., vary = False) #
    params.add('kpd', value=0.001, min=0.0001, max=5., vary = False)#
    params.add('mup', value=0.66, min=0.0001, max=1., vary = False)
    params.add('kdip', value=0.104, min=0.0001, max=0.5)
    params.add('kdp', value=0.001, min=0.0001, max=1., vary = False)
    params.add('mud', value=0.087, min=0.0001, max=0.1, vary = False)
    params.add('kap', value=0.001, min=0.0001, max=0.1)
    params.add('kad', value=0.001, min=0.0001, max=0.5, vary = False)
    params.add('mua', value=2.2, min=0.0001, max=5., vary = False)#

    # fit model
    result = minimize(residual_stoch, params, args=(t_measured, uPCR_measured), method='leastsq')  # leastsq nelder
    # check results of the fit
    data_fitted = g_stoch(t_measured, init, result.params)

    # plot fitted data
    ax.plot(t_measured, data_fitted[:, 3], '-', linewidth=2, color='red', label='fitted data')
    ax.set_title("Stochastic Fit: " + data['ID'].to_list()[0])
    ax.legend()
    ax.set_xlim([0, max(t_measured) + 20])
    # display fitted statistics

    #plt.show()

    return (result)