# helper functions
import numpy as np
import pandas as pd
#rng = np.random.default_rng()
import pandas as pd
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from lmfit import minimize, Parameters

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