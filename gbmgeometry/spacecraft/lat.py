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

        super(LAT, self).__init__(x_origin=x_origin,
                                  y_origin=y_origin,
                                  z_origin=z_origin,
                                  height=height,
                                  x_width=x_width,
                                  y_width=y_width,
                                  color='k'
                                  )


class LAT_radiator_minus(Volume):
    def __init__(self):
        # cm

        height = 158.6
        z_origin = 158.6 / 2.
        x_origin = 0.
        y_origin = -96.2
        x_width = 190.3
        y_width = 5.2

        super(LAT_radiator_minus, self).__init__(x_origin=x_origin,
                                                 y_origin=y_origin,
                                                 z_origin=z_origin,
                                                 height=height,
                                                 x_width=x_width,
                                                 y_width=y_width,
                                                 color='darkgrey'
                                                 )


class LAT_radiator_plus(Volume):
    def __init__(self):
        # cm

        height = 158.6
        z_origin = 158.6 / 2.
        x_origin = 0.
        y_origin = 96.2
        x_width = 190.3
        y_width = 5.2

        super(LAT_radiator_plus, self).__init__(x_origin=x_origin,
                                                y_origin=y_origin,
                                                z_origin=z_origin,
                                                height=height,
                                                x_width=x_width,
                                                y_width=y_width,
                                                color='darkgrey'
                                                )