---
jupyter:
  jupytext:
    formats: ipynb,md
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.2'
      jupytext_version: 1.8.0
  kernelspec:
    display_name: Python 3
    language: python
    name: python3
---

# The PositionInterpolator

Here we cover the detail of the **PositionInterpolator**. This tool allows you to gather lots of information about what is occurring during an orbit or trigger.


## The Constructor
The constructor is build a **PositionInterpolator** from both a Trigdat file and a daily position history file. 

```python
from gbmgeometry import PositionInterpolator
from gbmgeometry.utils.package_utils import get_path_of_data_file

trigdat_file = get_path_of_data_file('glg_trigdat_all_bn080916009_v02.fit')

pi = PositionInterpolator.from_trigdat(trigdat_file)

posthist_file = get_path_of_data_file('glg_poshist_all_151013_v00.fit')

pi = PositionInterpolator.from_poshist(posthist_file)
```

You can even build them from HdF5 file versions of the normal fits files. What? Yes we provide a converter for that as the HDF5 is much smaller and faster to load

```python
from gbmgeometry import convert_poshist2hdf5, convert_trigdat2hdf5

convert_poshist2hdf5(posthist_file, 'my_posthist.h5')

convert_trigdat2hdf5(trigdat_file, 'my_trigdat.h5')
```

and to load them

```python

trigdat_file = get_path_of_data_file('trigdat.h5')

pi = PositionInterpolator.from_trigdat_hdf5(trigdat_file)

posthist_file = get_path_of_data_file('posthist.h5')

pi = PositionInterpolator.from_poshist_hdf5(posthist_file)
```

## Plotting
We will do the fun stuff first. As the **PositionInterpolator** contains a lot of information about where Fermi is, we can use it to plot where the spacecraft is relative to Earth.




### Static

```python
from gbmgeometry import plot_in_space

pi = PositionInterpolator.from_poshist_hdf5(posthist_file)

# get the min and max time store in the file
tmin, tmax = pi.minmax_time()

plot_in_space(pi, tmin+1500);
```

Cool. But where are the detectors pointing?

```python
plot_in_space(pi, tmin+2000, show_detector_pointing=True, show_orbit=False, earth_time='day');
```

Where is the moon? 

**That's no moon.**

```python
plot_in_space(pi, tmin, show_detector_pointing=True,
              show_moon=True, 
              earth_time='midnight',
              show_orbit=False);
```

Yes, yes, you can also show the Sun... but it is not very useful at the moment.


### Animation

Let's go crazy and watch Fermi orbit

```python
from gbmgeometry import animate_in_space
animate_in_space(pi,
                 n_step=200,
                 interval=1100,
                 show_stars=True, show_detector_pointing=True, realistic=True, earth_time='day');
```

## Basic functions

We can check to see if Fermi is active at any time. For example, it is shut off during SAA passage.

```python
pi.is_fermi_active(tmin)
```

```python
import matplotlib.pyplot as plt
%matplotlib inline

fig, ax = plt.subplots()

time = np.linspace(tmin, tmax, 5000)

ax.plot(time-tmin,pi.is_fermi_active(time) )
```

```python

```
