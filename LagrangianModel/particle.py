
class Particle:
    """Object to model a convection or a radiation particle.

    Params:
    ------
    state  (float): For each particle, a value 0 or 1 ("cold" or "ignited"). 
                    If "ignited", it moves until it dies according to a decay.
    tmem   (float): Timescale of the decay (s) for each particle - for now, this
                    reflects the QMAX
    height (float): Height of release (m) for each particle - for a 3-D implementation
    x      (float): x-position of each particle
    y      (float): y-position of each particle
    
    factor (float): From 0 to 1, according to a Gaussian, centred at max burning
    ux     (float): mean wind speed in X-dir (later, it can be f(Z))
    vy     (float): mean wind speed in Y-dir


    Attributes:
    -----------
    props_name: all particle's properties name
                READ ONLY
    values: a dict containing all values for each key in props
            READ ONLY


    Methods:
    --------

    update(self, **kwargs): update the particle's value using kwargs dict

    get_from_keys(self, keys): return a tuple with the values of each property
                               in the iterable 'keys' IN THE GIVEN ORDER.

    Example:
    -------

    >>> particle = Particle(state=1, x=20, y=20)
    >>> particle.x                                # same as particle["x"]
    20
    >>> particle.update(x=0.3, vy=6)
    >>> particle
    Particle(state=1, x=0.3, y=20, ux=0, vy=6)
    >>> particle.props_name
    {'x', 'factor', 'ux', 'type', 'y', 'vy', 'tmem', 'state', 'height'}
    >>> particle.get_from_keys(["state", "ux", "vy"])
    (1, 0, 6)


    time, with spread so as to account for IGNTIME & BURNDUR.
    In this version, all cells emit equal number of particles, but they
    could have different MEMORY to account for different QMAX
    """

    def __init__(self, **kwargs):
        self._properties = {"state", "height", "x", "y", "tmem", "ux", "vy", "type", "factor"}
        self._props = {"state", "height", "x", "y", "ux", "vy"}
        self._values = {props: 0 for props in self._properties}
        self.update(**kwargs)
    
    @property
    def values(self):
        """Get a key: value dictionnary containning the properties
        of the particle."""
        return self._values

    @property
    def props_name(self):
        """Get the particle's properties' name."""
        return self._properties

    def update(self, **kwargs):
        for key, value in kwargs.items():
            if key in self._properties:
                self._values[key] = value

    def get_from_keys(self, keys):
        val = ()
        for key in keys:
            val += (self[key], )
        return val

    def __getattr__(self, key):
        if key in self._properties:
            return self._values[key]
        else:
            raise AttributeError(f"particle does not have property '{key}'")

    def __setattr__(self, key, value):
        if key in {"_properties", "_values", "_props"}:
            object.__setattr__(self, key, value)
        elif key in self._properties:
            self._values[key] = value

    def __getitem__(self, key):
        if key in self._properties:
            return self._values[key]
        raise KeyError(f"particle doesn't have property {key}.")
    
    def __add__(self, value):
        """Add 'value' to each property of the particle excluding the state."""
        if isinstance(value, (int, float)):
            state = self._values["state"]
            for key in self._props:
                self._values[key] += value
            self._values["state"] = state
        else:
            raise ValueError("Addition with type {} not supported".format(type(value)))
        return self.__class__(**self._values)

    def __radd__(self, value):
        return self.__add__(value)
    
    # Comparison of two particles or a particle and a float.
    # Here we just consider the state.
    def __eq__(self, particle):
        if isinstance(particle, Particle):
            return self._values["state"] == particle["state"]
        elif isinstance(particle, (int, float)):
            return self._values["state"] == particle
    
    def __ne__(self, particle):
        return not self.__eq__(particle)
    
    def __gt__(self, particle):
        if isinstance(particle, Particle):
            return self._values["state"] > particle["state"]
        elif isinstance(particle, (int, float)):
            return self._values["state"] > particle
    
    def __ge__(self, particle):
        return self.__gt__(particle) or self.__eq__(particle)

    def __lt__(self, particle):
        if isinstance(particle, Particle):
            return self._values["state"] < particle["state"]
        elif isinstance(particle, (int, float)):
            return self._values["state"] < particle

    def __le__(self, particle):
        return self.__lt__(particle) or self.__eq__(particle)

    def __repr__(self):
        properties = "state x y ux vy".split()
        state, x, y, ux, vy = self.get_from_keys(properties)
        return f"Particle({state=}, {x=}, {y=}, {ux=}, {vy=})"

if __name__ == '__main__':
    particle =Particle(state=1, x=20, y=20)
    print(particle)
    particle.update(x=0.3, vy=6)
    print(particle)
    print(particle.props_name)
    print(particle.get_from_keys(["state", "ux", "vy"]))
    print(Particle.__doc__)