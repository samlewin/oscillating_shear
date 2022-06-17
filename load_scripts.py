import numpy as np
import h5py
#These are python scripts for quick loading of the data into numpy array format. 
class sim_out:
    """Load data from a given directory rundir. Data is loaded\
    into a dictionary with keys 'mean', 'movie' and 'tke'.
    start and end denote the corresponding time steps of interest."""
    def __init__(self, rundir, start, end):
        try:
            self.data = {'mean': h5py.File(rundir + 'mean.h5','r'),   \
                'movie': h5py.File(rundir + 'movie.h5','r'), \
                'tke': h5py.File(rundir + 'tke.h5', 'r')}#), \
        except:
            self.data = {'mean': h5py.File(rundir + 'mean.h5','r'),   \
                'tke': h5py.File(rundir + 'tke.h5', 'r')}
        self.start = start
        self.end = end
        self.grid = np.array(self.data['mean']['gyf/0001'])
    def integrate_y(self, vec, grid_uniform=False, section=None):
        """Integrate a vertical profile over the y-values in self.grid. If
        grid_uniform is True assumes a uniform grid. A sub-section of a grid
        to integrate over can be specified using: 
        section=np.arange(start_index, end_index)."""
        if grid_uniform:
            grid_size = vec.shape[0] - 1
            int_y = 0
            for i in range(grid_size-1):
                int_y += 0.5*(vec[i] + vec[i+1])/grid_size
            return int_y
        else:
            if section is None:
                gyf = self.grid
            else:
                gyf = self.grid[section]
                vec = vec[section]
            grid_size = gyf.shape[0]
            DYF = gyf[1:] - gyf[:-1]
            LY = gyf[-1] - gyf[0]
            int_y = 0
            for i in range(grid_size-1):
                int_y += 0.5*(vec[i] + vec[i+1])*DYF[i]/LY
            return int_y
        
def to_timestep(i):
    """Convert integer timestep to a string for loading from file"""
    if i<10:
        dname='000'+str(i)
    elif i<100:
        dname='00'+str(i)
    elif i<1000:
        dname='0'+str(i)
    else:
        dname=''+str(i)           
    return dname

def sim_scalar_timeseries(rundir, scalar):
    """Outputs a vector containing the time series evolution
    of a scalar output from diablo.    
    Note scalar must be a key of data.rundir['mean']""" 
    scalar_ts = np.zeros((rundir.end-rundir.start))
    for i in range(rundir.start, rundir.end):
        if i+1<10:
            dname='000'+str(i+1)
        elif i+1<100:
            dname='00'+str(i+1)
        elif i+1<1000:
            dname='0'+str(i+1)
        else:
            dname=''+str(i+1)
        scalar_ts[i-rundir.start] = np.array(rundir.data['mean'][scalar + '/' + dname])
    return scalar_ts

def sim_mean_timeseries(rundir, profile, section=None):
    """Outputs a vector containing the time series evolution
    of the integral mean of a vertical profile output from diablo.
    Integral is taken over a sub-section of the y-grid if specified by
    the section parameter (use np.arange).
    Note profile must be a key of data.rundir['mean']."""
    mean_ts = np.zeros((rundir.end - rundir.start))
    for i in range(rundir.start, rundir.end):
        if i+1<10:
            dname='000'+str(i+1)
        elif i+1<100:
            dname='00'+str(i+1)
        elif i+1<1000:
            dname='0'+str(i+1)
        else:
            dname=''+str(i+1)
        if profile == 'epsilon' or profile == 'epsilon_prime' or profile == 'epsilon_3d':
            mean_ts[i-rundir.start] = rundir.integrate_y(np.array(rundir.data['tke'][profile + '/' + dname]),
                                               section=section)
        else:
            mean_ts[i-rundir.start] = rundir.integrate_y(np.array(rundir.data['mean'][profile + '/' + dname]),
                                               section=section)
    return mean_ts    

def sim_profile_timeseries(rundir, profile, section=None): 
    """Outputs an array containing the evolution of a vertical profile output
    from diablo. Vertical profile is limited to a sub-section of the y-grid if
    specified by the section parameter (use np.arange)."""
    if section is None:
        profile_ts = np.zeros((rundir.end-rundir.start, rundir.grid.shape[0]))
    else:
        profile_ts = np.zeros((rundir.end-rundir.start, section.shape[0]))
    if profile == 'epsilon' or profile == 'epsilon_prime' or profile=='epsilon_3d':
        data_dict = rundir.data['tke']
    else:
        data_dict = rundir.data['mean']
    for i in range(rundir.start, rundir.end):
        if i+1<10:
            dname='000'+str(i+1)
        elif i+1<100:
            dname='00'+str(i+1)
        elif i+1<1000:
            dname='0'+str(i+1)
        else:
            dname=''+str(i+1)
        if section is None:
            profile_ts[i-rundir.start, :] = np.array(data_dict[profile + '/' + dname])
    return profile_ts