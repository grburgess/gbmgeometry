from gbmgeometry.spacecraft.geometry import Volume


class SolarPanel(Volume):
    def __init__(self, sign, name):
        lenght = 492.3
        width = 151.0

        xtra_y_mount = 150.0

        x_origin = 0.0
        y_origin = sign * (lenght / 2.0 + xtra_y_mount)
        z_origin = 41.5

        y_width = lenght
        x_width = 1.0
        height = width

        super(SolarPanel, self).__init__(
            name=name,
            x_origin=x_origin,
            y_origin=y_origin,
            z_origin=z_origin,
            height=height,
            x_width=x_width,
            y_width=y_width,
            color="#3498DB",
            active_surfaces=["+x", "-x"],
        )


class SolarPanelPlus(SolarPanel):
    def __init__(self):
        super(SolarPanelPlus, self).__init__(1.0, "Solar Panel +")


class SolarPanelMinus(SolarPanel):
    def __init__(self):
        super(SolarPanelMinus, self).__init__(-1.0, "Solar Panel -")
