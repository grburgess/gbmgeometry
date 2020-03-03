from gbmgeometry.spacecraft.geometry import Volume


class LAT(Volume):
    def __init__(self):
        # cm

        height = 106.0
        z_origin = 158.6 + 106.0 / 2.0
        x_origin = 0.0
        y_origin = 0.0
        x_width = 190.3
        y_width = 190.3

        active_surfaces = ["+x", "-x", "+y", "-y", "+z"]

        super(LAT, self).__init__(
            name="LAT",
            x_origin=x_origin,
            y_origin=y_origin,
            z_origin=z_origin,
            height=height,
            x_width=x_width,
            y_width=y_width,
            color="#777777",
            active_surfaces=active_surfaces,
        )


class LATRadiator(Volume):
    def __init__(self, name, sign, active_surfaces):
        # cm

        height = 158.6
        z_origin = 158.6 / 2.0
        x_origin = 0.0
        y_origin = sign * 96.2
        x_width = 190.3
        y_width = 5.2

        super(LATRadiator, self).__init__(
            name=name,
            x_origin=x_origin,
            y_origin=y_origin + sign * y_width / 2,
            z_origin=z_origin,
            height=height,
            x_width=x_width,
            y_width=y_width,
            color="darkgrey",
            active_surfaces=active_surfaces,
        )


class LATRadiatorMinus(LATRadiator):
    def __init__(self):
        active_surfaces = ["+x", "-x", "+y", "-z"]

        super(LATRadiatorMinus, self).__init__("LAT Radiator-", -1, active_surfaces)


class LATRadiatorPlus(LATRadiator):
    def __init__(self):
        active_surfaces = ["+x", "-x", "-y", "-z"]

        super(LATRadiatorPlus, self).__init__("LAT Radiator+", 1, active_surfaces)
