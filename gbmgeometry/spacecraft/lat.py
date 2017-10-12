from gbmgeometry.spacecraft.geometry import Volume


class LAT(Volume):

    def __init__(self):
        # cm

        height = 106.
        z_origin = 158.6 + 106. / 2.
        x_origin = 0.
        y_origin = 0.
        x_width = 190.3
        y_width = 190.3

        super(LAT, self).__init__(name='LAT',
                                  x_origin=x_origin,
                                  y_origin=y_origin,
                                  z_origin=z_origin,
                                  height=height,
                                  x_width=x_width,
                                  y_width=y_width,
                                  color='#421352'
                                  )


class LATRadiator(Volume):
    def __init__(self,name,sign):
        # cm

        height = 158.6
        z_origin = 158.6 / 2.
        x_origin = 0.
        y_origin = sign*96.2
        x_width = 190.3
        y_width = 5.2

        super(LATRadiator, self).__init__(name=name,
                                                 x_origin=x_origin,
                                                 y_origin=y_origin,
                                                 z_origin=z_origin,
                                                 height=height,
                                                 x_width=x_width,
                                                 y_width=y_width,
                                                 color='darkgrey'
                                                 )


class LATRadiatorMinus(LATRadiator):
    def __init__(self):


        super(LATRadiatorMinus, self).__init__('LAT Radiator-', -1)

class LATRadiatorPlus(LATRadiator):
    def __init__(self):


        super(LATRadiatorPlus, self).__init__('LAT Radiator+',1)

