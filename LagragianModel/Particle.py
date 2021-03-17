

class Particle:
    """ Class that model one particle """

    def __init__(self, **kwargs):
        """
            INITIALISE PARTICLE
            
            state  = For each particle, a value 0 or 1 ("cold" or "ignited"). 
                     If "ignited", it moves until it dies according to a decay.
            tmem   = Timescale of the decay (s) for each particle -  for now, this
                     reflects the QMAX
            height = Height of release (m) for each particle - for a 3-D implementation
            x      = x-position of each particle
            y      = y-position of each particle
            
            factor = From 0 to 1, according to a Gaussian, centred at max burning
            ux     = mean wind speed in X-dir (later, it can be f(Z))
            vy     = mean wind speed in Y-dir
            time, with spread so as to account for IGNTIME & BURNDUR.
            In this version, all cells emit equal number of particles, but they
            could have different MEMORY to account for different QMAX
        """
        self._properties = {"state", "height", "x", "y", "tmem", "ux", "vy", "type", "factor"}
        self.values = {props: 0 for props in self._properties}
        self.update(**kwargs)
    
    def update(self, **kwargs):
        for key in kwargs:
            if key in self._properties:
                self.values[key] = kwargs[key]
    
    def update_key(self, key, value):
        """ Increment self.values[key] by value """
        if key in self._properties:
            self.values[key] += value

    def get_all(self):
        """ Returns all properties """
        props = self.values
        return props["state"], props["type"], props["x"], props["y"], \
            props["ux"], props["vy"], props["factor"], props["tmem"], \
            props["height"]

    def __getitem__(self, key):
        if key in self._properties:
            return self.values[key]
        raise KeyError(f"particle doesn't have property {key}.")