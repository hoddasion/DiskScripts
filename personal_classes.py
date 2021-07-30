# seperate file to contain class definitions
# created by Wilhelm Hodder on 12/03/2019

import assist_functions.auxiliary_functions as aux # personal auxiliary functions for data management and more
import numpy as np

class Variable_1D:
    """
    class for one dimensional data variables, e.g. time series of temperature or pressure, etc...
    """
    # constructors
    def __init__(self): # empty instantiation
        self.data = np.array([])
        self.time = np.array([])
        self.legs = []
    def __init__(self, data, time, name): # instantiation of class object with data inside
        # first turn data into numpy array and account for filler values
        data_ = np.float32(np.array(data))
        time_ = np.float32(np.array(time))
        try:
            data_[data_==-9999] = np.nan
            time_[time_==-9999] = np.nan
        except:
            print('Value error for filler to nan conversion.')
        # now assign
        self.data = data_
        self.time = time_
        self.name = name
        self.legs = []
        
    # manual destructor
    def erase(self): # [currently not comprehensive enough, needs changing]
        del self.data; del self.time; del self.name
        
    # modifying functions
    def change_data(self, new_data, new_time): # [erases current data contained and] replaces it with new data
        self.data = np.array(new_data)
        self.time  = np.array(new_time)
    def change_name(self, new_name):
        self.name = new_name
    def add_leg(self, lower, upper, name):
        leg = Variable_leg(self.data, self.time, lower, upper, name)
        self.legs.append(leg)
    
    # stats
    def calc_mean(self):
        self.mean = np.mean(self.data)
    def calc_meadian(self):
        self.median = np.median(self.data)
    def calc_std(self):
        self.standard = np.std(self.data)
        
        
        
class Variable_leg(Variable_1D):
    
    # constructors
    def __init__(self, parent_data, parent_time, lower, upper, name):
        # instantiate directly from base class object; parent must be base class object
        self.data, self.time = aux.extract_chunk(lower, upper, parent_data, parent_time)
        self.name = name
        self.upper = upper
        self.lower = lower
        
    # modifying functions
    def reset_bounds(self, parent, lower, upper): # set new boundaries and change data and time accordingly from parent object
        del self.lower; self.lower = lower
        del self.upper; self.upper = upper
        del self.data; del self.time
        self.data, self.time = aux.extract_chunk(self.lower, self.upper, parent.data, parent.time)
        
    def stats(self): # calculate standard deviation and mean of data
        data_ = self.data[np.where(self.data != np.nan)] # exclude nans from calculation
        self.mean = np.mean(data_)
        self.median = np.median(data_)
        self.standard = np.std(data_)
        
class VarContainer:
    """
    Class to contain and manage large amounts of Variable_1D objects
    """
    # constructor
    def __init__(self, object_list):
        self.variables = object_list