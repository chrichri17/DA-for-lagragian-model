

class Ensemble:
    """ Object to store all the member of our ensemble kalman filter

    Attributes
    ----------
    list (list) :  container storing all the member


    Methods
    -------
    append(member) : add a new member to the ensemble

    propagate(steps) : propagate the lagragian model for
                       each member of the ensemble
    """
    def __init__(self):
        self.list = []
    
    def __iter__(self):
        return iter(self.list)
    
    def __getitem__(self, key):
        return self.list[key]
    
    def append(self, member):
        self.list.append(member)
    
    def propagate(self, steps):
        for member in self.list:
            member.propagate(steps)
