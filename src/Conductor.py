class conductor:
    def __init__(self, conductorJson, name):
        self.name = name
        self.circ_mils = conductorJson['circ_mils']
        self.out_diam = conductorJson['out_diam']
        self.gmr = conductorJson['gmr']
        self.r_60 = conductorJson['r_60']
        self.amp_cap = conductorJson['amp_cap']
