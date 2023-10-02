# Subduction zone landscapes libary 

This library is designed to assist you in calculating the surface deformation field resulting from numerous inelastic earthquake cycles.

**Usage:**

Duplicate the library using GitHub and then use the function `funsLongTermUpliftClassV3.ComputeLongTermDisplacement(dipAngle, distanceLocked, distanceSemiLocked, cutoff)`

- `dipAngle`: Dip angle of the forearc [degrees]
- `distanceLocked`: Surface distance locked from the trench [km] (d1)
- `distanceSemiLocked`: Horizontal distance from the trench where the interface is fully creeping [km] (d2)
- `cutoff`: For simplicity, use 0.05

**Geometry:**
```
   surface

t      d1    d2
---------------|
\              |
  \    forearc |    Depth
    \          |
      \        |
        \      |
          \    |
            \  |
              \|

```           
                  
- From the trench to d1 - fully locked (coupling = 1).
- From d1 to d2, coupling transitions linearly from fully locked (1) to fully creeping (0).

**Function returns:**

The function will return three objects:

1. An xarray object with inelastic vertical motion along the domain.
2. Pandas dataframe with a list of earthquakes producing inelastic deformation.
3. Max along-strike standard deviation.

Feel free to contact me for more complicated cases.
