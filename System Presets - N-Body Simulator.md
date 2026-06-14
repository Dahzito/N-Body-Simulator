**- SOLAR SYSTEM PRESET:**



bodies.append(Body(1.989e30, \[0,0,0], \[0, -15.97, -0.42], \[0, 0, 0], "Sun", 6.96e8, 5878, 0.0, 1000, "Star", \[1, 1, 0.7, 1], 0, 1))

\#The RGB values of the Star are automatically determined depending on their surface temperature, also calculated using the Mass-Luminosity Relation, and Stefan-Boltzmann's law.

\# Velocities updated to circular orbital speeds for the Sun mass and the current orbital radii.



bodies.append(Body(3.3011e23, \[0,0,0], \[0, 47526.0, 5835.5], \[5.79e10, 0, 0], "Mercury", 2.4397e6, 440, 0.09, 1000, "Planet", \[156/255, 102/255, 31/255, 0.8], 0, 0.9))



bodies.append(Body(4.8675e24, \[0,0,0], \[0, 34966.0, 2071.2], \[1.082e11, 0, 0], "Venus", 6.0518e6, 737, 0.76, 2000, "Planet", \[235/255, 191/255, 122/255, 0.8], 0.92, 0.02))



bodies.append(Body(5.97\*10\*\*24, \[0,0,0], \[0, 30041.0, 0.0], \[1.471e11, 0, 0], "Earth", 6.371e6, 280, 0.31, 2000, "Planet", \[15/255, 59/255, 141/255, 0.80], 0.46, 0.9))

\# Moon velocity set relative to Earth: Earth orbital speed + Moon orbital speed around Earth

bodies.append(Body(7.35\*10\*\*22, \[0,0,0], \[0, 31061.1, 91.7], \[1.471e11 + 3.633e8, 0, 0], "Moon", 1.737e6, 280, 0.11, 2300, "Moon", \[158/255, 158/255, 158/255, 0.80], 0.0, 0.95))



bodies.append(Body(6.4171e23, \[0,0,0], \[0, 24122.5, 779.2], \[2.279e11, 0, 0], "Mars", 3.3895e6, 210, 0.25, 2000, "Planet", \[178/255, 34/255, 34/255, 0.8], 0, 0.95))



bodies.append(Body(1.898e27, \[0,0,0], \[0, 13055.1, 296.9], \[7.785e11, 0, 0], "Jupiter", 6.9911e7, 120, 0.50, 2000, "Planet", \[205/255, 133/255, 63/255, 0.8], 0, 0.9))



bodies.append(Body(5.683e26, \[0,0,0], \[0, 9615.9, 417.3], \[1.433e12, 0, 0], "Saturn", 5.8232e7, 95, 0.34, 2000, "Planet", \[210/255, 180/255, 140/255, 0.8], 0, 0.9))



bodies.append(Body(8.681e25, \[0,0,0], \[0, 6792.2, 91.5], \[2.877e12, 0, 0], "Uranus", 2.5362e7, 76, 0.30, 2000, "Planet", \[173/255, 216/255, 230/255, 0.8], 0, 0.9))



bodies.append(Body(1.02413e26, \[0,0,0], \[0, 5427.0, 167.6], \[4.503e12, 0, 0], "Neptune", 2.4622e7, 72, 0.29, 2000, "Planet", \[72/255, 61/255, 139/255, 0.8], 0, 0.9))



\#Mass, acceleration, velocity, coordinates, name, radii, temperature, albedo, temp\_threshold, type, RGB colors (0-1), greenhouse effect (0-1), ε - Surface Emissitivity (0-1);







**- TEEGARDENS SYSTEM**



bodies.append(Body(0.0970\*M\_Sun, \[0,0,0], \[0, 0, 0], \[0, 0, 0], "Teegardens Star", 0.107\*R\_Sun, 5878, 0.0, 1000, "Red Dwarf Star", \[1, 1, 0.7, 1], 0, 1))

\#The RGB values of the Star are automatically determined depending on their surface temperature, also calculated using the Mass-Luminosity Relation, and Stefan-Boltzmann's law.



bodies.append(Body(1.16\*M\_Earth, \[0,0,0], \[0, 59406, 0.0], \[0.0259\*au\*(1-0.03), 0, 0], "Teegardens b", 1.16\*R\_Earth, 280, 0.31, 2000, "Planet", \[206/255, 104/255, 45/255, 0.8], 0.2, 0.9))



bodies.append(Body(1.05\*M\_Earth, \[0,0,0], \[0, 45271, 0], \[0.0455\*au\*(1-0.04), 0, 0], "Teegardens c", 1.05\*R\_Earth, 120, 0.50, 2000, "Planet", \[226/255, 83/255, 0/255, 0.8], 0.1, 0.9))



bodies.append(Body(0.82\*M\_Earth, \[0,0,0], \[0, 35384, 0], \[0.0791\*au\*(1-0.07), 0, 0], "Teegardens d", 0.82\*R\_Earth, 120, 0.50, 2000, "Planet", \[190/255, 154/255, 36/255, 0.8], 0.0, 0.9))







**- HD 28109 SYSTEM**



bodies.append(Body(1.227\*M\_Sun, \[0,0,0], \[0, 0, 0], \[0, 0, 0], "HD 28109 \*", 1.425\*R\_Sun, 5878, 0.0, 1000, "Star", \[1, 1, 0.7, 1], 0, 1))

\#The RGB values of the Star are automatically determined depending on their surface temperature, also calculated using the Mass-Luminosity Relation, and Stefan-Boltzmann's law.



bodies.append(Body(6.2\*M\_Earth, \[0,0,0], \[0, 88033, 19154], \[0.1357\*au\*(1-0.00710), 0, 0], "HD 28109 b", 2.69\*R\_Earth, 280, 0.31, 2000, "Planet", \[206/255, 104/255, 45/255, 0.8], 0.2, 0.9))



bodies.append(Body(9.20\*M\_Earth, \[0,0,0], \[0, 58619, 10819], \[0.308\*au\*(1-0.0039), 0, 0], "HD 28109 c", 4.13\*R\_Earth, 120, 0.50, 2000, "Planet", \[226/255, 83/255, 0/255, 0.8], 0.1, 0.9))



bodies.append(Body(5.681\*M\_Earth, \[0,0,0], \[0, 50844, 925.7], \[0.411\*au\*(1-0.0054), 0, 0], "HD 28109 d", 3.25\*R\_Earth, 120, 0.50, 2000, "Planet", \[190/255, 154/255, 36/255, 0.8], 0.0, 0.9))



\#Mass, acceleration, velocity, coordinates, name, radii, temperature, albedo, temp\_threshold, type, RGB colors (0-1), greenhouse effect (0-1), ε - Surface Emissitivity (0-1);

