# N-Body-Simulator
A continuation to my, original, 2D N-Body Simulator project ( https://github.com/Dahzito/N-Body-2D-Simulator ), but this time in 3D.

As more and more physics, and code, gets changed on this simulator, the more my "old" 2D N-Body Simulator gets outdated in comparison with this one.
Therefore, even if the "old" project is still pretty much usable, I will consider it deprecated, specially with the changes I am planning on making on this new and exciting continuation of my project.

----------------------

This program uses Newtonian Mechanics, and uses the Yoshida 4th Order Integrator to compute the velocity and position of every body that is added to the program.

Some features include: automatic star's luminosity and temperature calculation using the Mass-Luminosity relation and the Stefan-Boltzmann's law, force felt by every body in the N-Body system, a simple body vaporation system, flux received by every body in the n-body system, and then calculation of the temperature of every planet in it using the Stefan-Boltzmann's law, the body's albedo, surface emissivity and greenhouse factor, and gravity at the surface of every body to eventually see the Roche limit.

Stars have an automatic RGB value system calculation depending on their effective temperature, all the other bodies need a manual input for their RGBA values.

To add a new body to the program, add the following line of code, and change whatever you want:

bodies.append(Body(1.989e30, [0,0], [0, 0], [0, 0], "Sun", 6.96e8, 5878, 0.0, 1000, "Star", [1, 1, 0.7, 1], 0, 1)) #The RGB values of the Star are automatically determined depending on their surface temperature, also calculated using the Mass-Luminosity Relation, and Stefan-Boltzmann's law.

bodies.append(Body(5.97e24, [0,0], [0, 30290], [1.471e11, 0], "Earth", 6.371e6, 280, 0.31, 2000, "Planet", [15/255, 59/255, 141/255, 0.80], 0.46, 0.9)) bodies.append(Body(7.35e22, [0,0], [0, 31312], [1.471e11 + 3.633e8, 0], "Moon", 1.737e6, 280, 0.11, 2300, "Planet", [158/255, 158/255, 158/255, 0.80], 0.0, 0.95))

bodies.append(Body(1.898e27, [0,0], [0, 13070], [7.785e11, 0], "Jupiter", 6.9911e7, 120, 0.50, 2000, "Planet", [205/255, 133/255, 63/255, 0.8], 0, 0.9)) #Mass, acceleration, velocity, coordinates, name, radii, temperature, albedo, temp_threshold, type, RGB colors (0-1), greenhouse effect (0-1), ε - Surface Emissitivity (0-1);

(Or simply change the parameters of the already existing bodies)

This needs to be added after the piece of code that contains: class Body ; and the array bodies = [], btw

Runs on Python with a tad of C++.

----------------------

Some simulations I was able to run on this program:

<img width="1479" height="1526" alt="N-Body Simulator, C++ Test1, Plot 7 @ 1271 56 days" src="https://github.com/user-attachments/assets/f46eb0f0-54c3-46af-9916-bd84d05a80e0" />
<img width="1479" height="1526" alt="N-Body Simulator, C++ Test1, Plot 8 @ 1454 02 days" src="https://github.com/user-attachments/assets/1dde6877-82e5-425e-b051-0b71207f162c" />

<img width="3212" height="1594" alt="N-Body Simulator, HD 28109 System Plot 5 @ 906 67 days" src="https://github.com/user-attachments/assets/cf53b697-2b29-49d4-a04e-126ae277b20e" />
<img width="3212" height="1594" alt="N-Body Simulator, HD 28109 System Plot 7 @ 1271 56 days" src="https://github.com/user-attachments/assets/b58a3705-9a3d-45f9-9cc0-e229e1be5a0f" />

<img width="1479" height="1526" alt="N-Body Simulator, Orbits Visualization, T3 Plot 4 @ 2914 69 days" src="https://github.com/user-attachments/assets/fcc90244-1c32-4e40-b337-27ad702cc76e" />
<img width="1479" height="1568" alt="N-Body Simulator, Orbits Visualization, T3 Plot 9 @ 6564 12 days" src="https://github.com/user-attachments/assets/f97e6f40-f960-490b-9581-ac63a5321ce4" />
