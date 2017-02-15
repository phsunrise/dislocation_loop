import numpy as np

'''
Inputs:
    pos is a N*3 array, or [x,y,z]
    R is loop radius
Output:
    return an array of length N
'''
def dampingfunc(pos, R, form, params=[]):
    if form == 'nodamp':
        if len(pos.shape) == 1:
            return 1.
        else:
            return np.ones(len(pos))

    elif form == 'fermi': # params should have two arguments
        if not len(params) == 2:
            raise ValueError("Length of parameters for fermi function should be 2!")

        if len(pos.shape) == 1: # if pos is just a 1d array
            pos = np.array([pos])

        fr = R*np.sqrt(2.)
        hr = R*params[0]
        sigma = (fr**2.-hr**2.)/np.log(1/params[1]-1.)
        mu = hr**2./sigma
        disp = np.sqrt(pos[:,0]**2+pos[:,1]**2)
        #fer = (1.+np.exp(-mu))/(1.+np.exp(pos[:,2]**2/sigma-mu))*(disp<1.5*R)
        fer = (1.+np.exp(-mu))/(1.+np.exp(np.sum(pos**2, axis=1)/sigma-mu))
        if len(fer) == 1:
            return fer[0]
        else:
            return fer

    elif form == 'gaussian':
        if not len(params) == 1:
            raise ValueError("Length of parameters for gaussian function should be 1 (D/R)!")
        
        if len(pos.shape) == 1: # if pos is just a 1d array
            pos = np.array([pos])
        
        D = R*1.*params[0]
        gaus = np.exp(-0.5*(np.sum(pos**2, axis=1)/D**2))
        if len(gaus) == 1:
            return gaus[0]
        else:
            return gaus 
    


if __name__=='__main__':
    import matplotlib.pyplot as plt
    from W_parameters import R, funcform, funcparams 
    import sys
    print "R=%.1f" % R
    print "%s" % funcform
    print "params:", funcparams

    r = np.linspace(0., 400., 200)
    pos = np.zeros((len(r), 3))
    pos[:, 2] = r
    #plt.plot(r/R, fermi(pos, R, params=[2.5, 0.9]), label='orig.')
    plt.plot(r/R, dampingfunc(pos, R, funcform, funcparams))
    plt.xlabel("r/R")

    plt.savefig("fermifunc.pdf")
    plt.show()
