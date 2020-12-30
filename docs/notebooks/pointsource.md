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
    display_name: Python3
    language: python
    name: python3
---

# Plotting Points on the Sky

It is sometimes interesting to see how various point on the sky relate to the position of Fermi. 
Here we will demonstrate but plotting the position of a [BALROG](https://academic.oup.com/mnras/article-abstract/476/2/1427/4670828?redirectedFrom=fulltext) posterior on the sky.


```python
from gbmgeometry import plot_in_space, PositionInterpolator
from gbmgeometry.utils.plotting.sky_point import balrog_to_skypoints
from gbmgeometry.utils.package_utils import get_path_of_data_file

pi = PositionInterpolator.from_trigdat(get_path_of_data_file('balrog_trig.fits'))

skypoints = balrog_to_skypoints(get_path_of_data_file('balrog.fits'), new_nside=2**6, cmap='viridis', as_point=True)
```

We have constructed a series of **SkyPoints** from the BALROG posterior. We have chosen to plot them as points on the sky.

```python
plot_in_space(pi, 0, sky_points=skypoints,
              show_detector_pointing=True, 
              show_moon=True,earth_time='day', 
              show_stars=True);
```

Alternatively, we could plot them as rays from Fermi.

```python
skypoints = balrog_to_skypoints(get_path_of_data_file('balrog.fits'), new_nside=2**4, cmap='winter_r', as_point=False)

plot_in_space(pi, 0,
              sky_points=skypoints,
              earth_time='thats_no_moon', 
              show_stars=True);
```
