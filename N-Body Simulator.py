import builtins
from math import*
from sys import*
from os import*
import os
import json
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

import core_calculations

Name_sys:str = "HD 28109 System"  # Name of the System that is going to be simulated, going to be shown on the program and window of the same.
auto_center:bool = False          # This allows the Plot to be automatically centered to the center coordinates of all bodies in the simulation
k_update:float = 4                # Number of times per Δt_reset time in which the orbits will be updated

δt:float = 3 #days
δt_default:float = 2 #default timestep for normal conditions
Δt:float = 365.25*20 #days

Δt_reset:float = 365*2

ε:float = 0 #Gravitational softening
force_temp:float = 0 #Temporary variable for storing gravitational force

G = 6.67430 * 10**-11 #Gravitational constant
σ = 5.67037 * 10**-8  #Stefan-Boltzmann constant
b = 2.89777 * 10**-3 #Wien's displacement constant

M_Earth = 5.972e24 # Mass of Earth in kg
R_Earth = 6.371e6  # Radius

M_Sun = 1.989e30   # Mass of Sun in kg
R_Sun = 6.96e8     # Radius of Sun in meters
L_Sun = 3.828e26   # Luminosity of Sun in watts

M_Jupiter = 1.898e27 # Mass of Jupiter in kg
R_Jupiter = 6.9911e7 # Radius of Jupiter in meters

au = 1.49597871e11 # Astronomical Unit in meters

_t = 24*3600 #Conversion factor from days to seconds

MAX_TRAIL = 1000

class Body:
        def __init__(self, mass, acc, vel, coord, name, radii, temperature, albedo, temp_threshold, type, rgb=[0, 0, 0, 1], gh=0, emissivity=0.9):
            
            #Dynamic Properties
            self.acceleration = list(acc)   # List as [acc_x, acc_y]
            self.velocity = list(vel)       # List as [vel_x, vel_y]
            self.coordinates = list(coord)  # List as [coord_x, coord_y]
            self.prev_acceleration = list(acc)  # Store previous acceleration for velocity Verlet

            self.force = 0

            self.potential_energy = 0; self.kinetic_energy = 0; self.mechanical_energy = 0

            #Physical Properties 
            self.name = name
            self.radii = radii
            self.mass = mass

            self.temperature = temperature        # Kelvin
            self.temp_threshold = temp_threshold  # Arbitrary threshold for temperature-based effects (e.g., radiation pressure)

            self.gravity = G * self.mass / self.radii**2  # Surface gravity

            self.luminosity = σ * 4 * pi * self.radii**2 * self.temperature**4
            self.albedo = albedo  # Albedo for energy absorption calculations
            self.flux = 0         # Flux received from other bodies (for energy balance calculations)

            #Other properties
            self.type = type
            # Type of body (e.g., "Main Sequence Star", "Planet", "Asteroid")
            # Star Tupes: SuperGiant, Main Sequence, Red Giant, Red Dwarf, White Dwarf, Neutron Star, Black Hole (for visualization purposes only, no relativistic effects)

            self.status = "Stable"  # Status can be "Stable", "Torn Apart", for Roche's limit effects

            self.wavelength = 0 #Wavelength of the radiation emitted by the body, using Wien's Law
            self.rgb = rgb

            self.gh = gh   # Greenhouse effect factor for planets, where 0 means no greenhouse effect. 
                           # This is a simplified model and does not account for atmospheric composition or other factors that influence greenhouse effects.
            self.emissivity = emissivity  # Surface emissivity for thermal radiation calculations of planets.
            self.heat_capacity = 1000 * self.mass  # Simplified heat capacity proportional to mass (J/K), for temperature change calculations

# rgba(158, 158, 158, 0.80)
bodies = []

bodies.append(Body(0.0970*M_Sun, [0,0,0], [0, 0, 0], [0, 0, 0], "Teegardens Star", 0.1186490*R_Sun, 5878, 0.0, 1000, "Star", [1, 1, 0.7, 1], 0, 1))
#The RGB values of the Star are automatically determined depending on their surface temperature, also calculated using the Mass-Luminosity Relation, and Stefan-Boltzmann's law.

bodies.append(Body(1.16*M_Earth, [0,0,0], [0, 59406, 0.0], [0.0259*au*(1-0.03), 0, 0], "Teegardens b", 1.16*R_Earth, 280, 0.31, 2000, "Planet", [206/255, 104/255, 45/255, 0.8], 0.2, 0.9))

bodies.append(Body(1.05*M_Earth, [0,0,0], [0, 45271, 0], [0.0455*au*(1-0.04), 0, 0], "Teegardens c", 1.05*R_Earth, 120, 0.50, 2000, "Planet", [226/255, 83/255, 0/255, 0.8], 0.1, 0.9))

bodies.append(Body(0.82*M_Earth, [0,0,0], [0, 35384, 0], [0.0791*au*(1-0.07), 0, 0], "Teegardens d", 0.82*R_Earth, 120, 0.50, 2000, "Planet", [190/255, 154/255, 36/255, 0.8], 0.0, 0.9))

#Mass, acceleration, velocity, coordinates, name, radii, temperature, albedo, temp_threshold, type, RGB colors (0-1), greenhouse effect (0-1), ε - Surface Emissitivity (0-1);

def Update_():
    global δt, ε, force_temp
    radius = []
    dt_seconds = δt * _t  # Convert to seconds

    force_temp_dx = 0; force_temp_dy = 0; force_temp_dz = 0;
    
    ε = min(b.radii for b in bodies) * 1 # Softening length based on smallest radius

    #Calculate the distance of every object on the class Body
    for b in range(len(bodies)):
        bodies[b].acceleration = [0 , 0 , 0]
        bodies[b].flux = 0
        
    body_forces = [0.0 for _ in bodies]

    for i in range(len(bodies)):
        force_temp_dx = 0.0
        force_temp_dy = 0.0
        force_temp_dz = 0.0
        for j in range(len(bodies)):
            if i != j:
                dx = bodies[j].coordinates[0] - bodies[i].coordinates[0]
                dy = bodies[j].coordinates[1] - bodies[i].coordinates[1]
                dz = bodies[j].coordinates[2] - bodies[i].coordinates[2]

                r = sqrt(dx**2 + dy**2 + dz**2)
                radius.append(r)

                r_soft = sqrt(dx**2 + dy**2 + dz**2 + ε**2)
                bodies[i].acceleration[0] += (G * bodies[j].mass * dx) / r_soft**3
                bodies[i].acceleration[1] += (G * bodies[j].mass * dy) / r_soft**3
                bodies[i].acceleration[2] += (G * bodies[j].mass * dz) / r_soft**3

                # Accumulate flux from body j
                bodies[i].flux += bodies[j].luminosity / (4 * pi * r**2)

                force_temp_dx += G * bodies[i].mass * bodies[j].mass * dx / r_soft**3
                force_temp_dy += G * bodies[i].mass * bodies[j].mass * dy / r_soft**3
                force_temp_dz += G * bodies[i].mass * bodies[j].mass * dz / r_soft**3
                
                # Roche/tidal disruption check for smaller bodies near a much larger primary
                if bodies[i].type != "Star" and bodies[j].mass > bodies[i].mass and bodies[i].mass > 0:
                    rho_primary = bodies[j].mass / ((4/3) * pi * bodies[j].radii**3)
                    rho_secondary = bodies[i].mass / ((4/3) * pi * bodies[i].radii**3)
                    if rho_secondary > 0:
                        roche_limit = 2.44 * bodies[j].radii * (rho_primary / rho_secondary)**(1/3)
                        if r < roche_limit:
                            bodies[i].status = "Torn Apart"

        body_forces[i] = sqrt(force_temp_dx**2 + force_temp_dy**2 + force_temp_dz**2)
        bodies[i].force = body_forces[i]
        


    #adaptative timestep based on minimum distance to improve accuracy during close encounters
    try:
        a_max = max(np.linalg.norm(body.acceleration)
            for body in bodies)

        δt = np.clip(
        0.05 / np.sqrt(a_max),
        0.01,
        50
        )
    except ValueError:
        δt = 6*δt_default  # Default large distance if no pairs exist (e.g., single body)
    
    # Calculate temperature and mass loss after all flux has been accumulated
    for i in range(len(bodies)):
        temp_new = (bodies[i].flux * ((1 - bodies[i].albedo) / (4*σ*bodies[i].emissivity)))**0.25 if bodies[i].flux > 0 else 0
        
        k0 = 1e5; T0 = 1000
        k = k0 * exp(temp_new / T0)
        
        # Temperature, luminosity, and mass updates based on flux and temperature thresholds
        # NOTE: disruption/status changes are handled by the Roche/tidal check
        # performed earlier (which sets `bodies[i].status = "Torn Apart"` when
        # appropriate). Do not change status here based on bulk acceleration,
        # since tidal/stress calculations are a better predictor.
        if bodies[i].type in ("Planet", "Moon","Asteroid"):
            bodies[i].temperature = temp_new
            if bodies[i].mass <= k * dt_seconds:
                bodies[i].status = "Extreme Mass Loss"
            
            if bodies[i].type in ("Planet", "Moon") and bodies[i].gh >0:
                bodies[i].temperature = ( (bodies[i].flux * (1 - bodies[i].albedo) * (1 + bodies[i].gh)) / (4 * σ * bodies[i].emissivity) ) ** 0.25

            if temp_new > bodies[i].temp_threshold:
                bodies[i].mass -= k * dt_seconds

            power_in = bodies[i].flux * (1 - bodies[i].albedo) * (pi * bodies[i].radii**2)
            power_out = (σ * (bodies[i].temperature**4) * (4 * pi * bodies[i].radii**2) * bodies[i].emissivity) / (1 + bodies[i].gh)

            dT_dt = (power_in - power_out) / bodies[i].heat_capacity
            bodies[i].temperature += dT_dt * dt_seconds
        
        #Temperature, luminosity, and wavelength calculations
        elif "Star" in bodies[i].type:
            Mass_Relative = bodies[i].mass / M_Sun
            if Mass_Relative < 0.43: bodies[i].luminosity = 0.23*L_Sun*((Mass_Relative)**2.3)
            elif 0.43 < Mass_Relative < 2: bodies[i].luminosity = L_Sun*((Mass_Relative)**4)
            elif 2 < Mass_Relative < 55: bodies[i].luminosity = 1.4*L_Sun*((Mass_Relative)**3.5)
            elif Mass_Relative > 55: bodies[i].luminosity = 32000*L_Sun*(Mass_Relative)

            bodies[i].temperature = (bodies[i].luminosity / (4 * pi * bodies[i].radii**2 * σ))**0.25
            bodies[i].wavelength = b / bodies[i].temperature

            T_scaled = max(10, min(400, bodies[i].temperature / 100))
            # RGB color calculation based on temperature using a simplified blackbody approximation
            if T_scaled <= 66:
                bodies[i].rgb[0] = 1
            else:
                bodies[i].rgb[0] = (329.698727446 * ((T_scaled - 60) ** -0.133205)) / 255

            if T_scaled <= 66:
                bodies[i].rgb[1] = (99.4708025861 * log(T_scaled) - 161.1195681661) / 255
            else:
                bodies[i].rgb[1] = (288.1221695283 * ((T_scaled - 60) ** -0.0755148492)) / 255

            if T_scaled >= 66:
                bodies[i].rgb[2] = 1
            elif T_scaled <= 19:
                bodies[i].rgb[2] = 0
            else:
                bodies[i].rgb[2] = (138.5177312231 * log(T_scaled - 10) - 305.0447927307) / 255

            bodies[i].rgb[0] = max(0, min(1, bodies[i].rgb[0]))
            bodies[i].rgb[1] = max(0, min(1, bodies[i].rgb[1]))
            bodies[i].rgb[2] = max(0, min(1, bodies[i].rgb[2]))

        else: 
            if temp_new > bodies[i].temperature: bodies[i].temperature = temp_new

    # Yoshida 4th Order Integrator: for bigger timesteps and better long-term energy conservation
    
    n_b = len(bodies)
    if n_b > 0:
        pos_x = [body.coordinates[0] for body in bodies]
        pos_y = [body.coordinates[1] for body in bodies]
        pos_z = [body.coordinates[2] for body in bodies]

        vel_x = [body.velocity[0] for body in bodies]
        vel_y = [body.velocity[1] for body in bodies]
        vel_z = [body.velocity[2] for body in bodies]

        masses = [b.mass for b in bodies]

        results = core_calculations.Yoshida4thIntegrator(
                pos_x, pos_y, pos_z,
                vel_x, vel_y, vel_z,
                masses, ε, G, dt_seconds
                )

        pos_x, pos_y, pos_z, vel_x, vel_y, vel_z = results

        for i in range(n_b):
            bodies[i].coordinates[0] = pos_x[i]
            bodies[i].coordinates[1] = pos_y[i]
            bodies[i].coordinates[2] = pos_z[i]
            bodies[i].velocity[0] = vel_x[i]
            bodies[i].velocity[1] = vel_y[i]
            bodies[i].velocity[2] = vel_z[i]

#------ Visualization of the N Body system using Matplotlib ------

def calculate_energy():
    energies = [0] 
    #First value is for the total energy of the system, then each body will have its own energy value

    for i in range(len(bodies)):
        b = bodies[i]
        kinetic_energy = 0.5 * b.mass * (b.velocity[0]**2 + b.velocity[1]**2 + b.velocity[2]**2)
        bodies[i].kinetic_energy = kinetic_energy

        # Potential energy: -G * m1 * m2 / r
        potential_energy = 0
        for j in range(len(bodies)):
            if i != j:
                dx = bodies[j].coordinates[0] - b.coordinates[0]
                dy = bodies[j].coordinates[1] - b.coordinates[1]
                dz = bodies[j].coordinates[2] - b.coordinates[2]
                r = sqrt(dx**2 + dy**2 + dz**2 + ε**2)
                potential_energy -= 0.5 * (G * b.mass * bodies[j].mass) / r

                bodies[i].potential_energy = potential_energy

        # Total energy
        total_energy = kinetic_energy + potential_energy
        bodies[i].mechanical_energy = total_energy
        energies[0] += total_energy
        energies.append(total_energy)

    return energies

def calculate_forces():
    forces = [] 
    #Each body will have its own force value

    for i in range(len(bodies)):
        force_temp_dx = 0
        force_temp_dy = 0
        force_temp_dz = 0
        for j in range(len(bodies)):
            if i != j:
                dx = bodies[j].coordinates[0] - bodies[i].coordinates[0]
                dy = bodies[j].coordinates[1] - bodies[i].coordinates[1]
                dz = bodies[j].coordinates[2] - bodies[i].coordinates[2]
                r = sqrt(dx**2 + dy**2 + dz**2 + ε**2)

                force_temp_dx += G*bodies[i].mass*bodies[j].mass*dx/r**3
                force_temp_dy += G*bodies[i].mass*bodies[j].mass*dy/r**3
                force_temp_dz += G*bodies[i].mass*bodies[j].mass*dz/r**3

        # Total energy
        total_force = sqrt(force_temp_dx**2 + force_temp_dy**2 + force_temp_dz**2)  # Magnitude of the total force
        forces.append(total_force)

    return forces

def Visualize(energy_reset_period=Δt_reset):
    plt.rcParams['font.family'] = 'monospace'
    plt.rcParams['font.weight'] = 'normal'
    plt.style.use('dark_background')
    """Visualize the n-body system with energy tracking.
    
    Args:
        energy_reset_period: Reset the energy graph every N days (default: 365)
        force_reset_period: Reset the force graph every N days (default: 365)
    """
    fig = plt.figure(num="N-Body Simulation - Full Tracking", figsize=(24, 18))
    try:
        fig.canvas.manager.set_window_title("N-Body Simulation - Full Tracking")
    except Exception:
        pass
    save_png_template = os.path.join(
        os.path.expanduser("~"),
        "Downloads",
        "Python Plots - 3D N-Body Simulator",
        "N-Body Simulator - Full Tracking - C++ Test4 - Plots",
        "N-Body Simulator, C++ Test4, Plot {num} @ {current_time:.2f} days.png"
    )
    save_json_template = os.path.join(
        os.path.expanduser("~"),
        "Downloads",
        "Python Plots - 3D N-Body Simulator",
        "N-Body Simulator - Full Tracking - C++ Test4 - Data",
        "N-Body Simulator, C++ Test4, Data {num} @ {current_time:.2f} days.json"
    )
    os.makedirs(os.path.dirname(save_png_template), exist_ok=True)
    os.makedirs(os.path.dirname(save_json_template), exist_ok=True)
    save_count = [1]
    next_save_time = [energy_reset_period - 2 * δt]  # Save just before the reset for better visualization of changes
    next_save_time_orbits = [energy_reset_period/k_update - 2* δt]
    manager = plt.get_current_fig_manager()

    try:
        manager.window.showMaximized()
    except AttributeError:
        try:
            manager.window.state('zoomed')
        except AttributeError:
            manager.full_screen_toggle()

    gs = fig.add_gridspec(3, 3, height_ratios=[2, 2, 1.5], width_ratios=[1.5, 1, 1], hspace=0.4, wspace=0.4)

    ax = fig.add_subplot(gs[0:2, 0], projection="3d")
    ax.set_box_aspect([1, 1, 1])
    fig.patch.set_facecolor("#000000")
    for axis in [ax.xaxis, ax.yaxis, ax.zaxis]:
        axis._axinfo["grid"]['color'] = (1, 1, 1, 0.3)
        axis._axinfo["grid"]['linewidth'] = 0.15

    ax.xaxis.set_pane_color((0.15, 0.15, 0.15, 0.1)) 
    ax.yaxis.set_pane_color((0.15, 0.15, 0.15, 0.1)) 
    ax.zaxis.set_pane_color((0.15, 0.15, 0.15, 0.1))

    ax.xaxis.label.set_color("#e0e0e067")
    ax.yaxis.label.set_color('#e0e0e067')
    ax.zaxis.label.set_color('#e0e0e067')
    ax.tick_params(colors='#b0b0b0', labelsize=8)
    ax.grid(True, which='both', linestyle='--', linewidth=0.15, alpha=0.25)

    ax_energy = fig.add_subplot(gs[0, 1])  # Energy plot
    ax_force = fig.add_subplot(gs[1, 1])   # Force plot
    ax_mass = fig.add_subplot(gs[0, 2])    # Mass plot
    ax_flux = fig.add_subplot(gs[1, 2])    # Flux plot
    ax_temperature = fig.add_subplot(gs[2, 2])  # Temperature plot
    
    for ax_sub in [ax_energy, ax_force, ax_mass, ax_flux, ax_temperature]:
        ax_sub.set_facecolor('#0d1117')
        ax_sub.grid(True, which='both', linestyle=':', linewidth=0.4, alpha=0.3, color='#444444')
        ax_sub.tick_params(colors='#b0b0b0', labelsize=8)
        for spine in ax_sub.spines.values():
            spine.set_edgecolor('#444444')
            spine.set_linewidth(0.8)
        ax_sub.xaxis.label.set_color('#e0e0e0')
        ax_sub.yaxis.label.set_color('#e0e0e0')
        ax_sub.title.set_color('#ffffff')

    ax.set_aspect('equal')
    ax.set_title("N-Body Simulation: Orbital Viewer", fontsize=13, pad=10, fontweight='bold')

    # auto limits based on initial system
    max_r = max(max(abs(x.coordinates[0]), abs(x.coordinates[1]), abs(x.coordinates[2])) for x in bodies)
    limit = max(max_r * 1.5, 1e10)

    ax.set_xlim(-limit, limit)
    ax.set_ylim(-limit, limit)
    ax.set_zlim(-limit, limit)

    ax.set_xlabel("x coordinates (m)")
    ax.set_ylabel("y coordinates (m)")
    ax.set_zlabel("z coordinates (m)")

    points = []
    trails = []
    archived_trails = []

    trails_x = [[] for _ in bodies]
    trails_y = [[] for _ in bodies]
    trails_z = [[] for _ in bodies]

    archived_mass_lines = []
    archived_flux_lines = []
    archived_temperature_lines = []
    archived_force_lines = []
    archived_energy_lines = []

    for i, b in enumerate(bodies):
        p, = ax.plot([], [], [], 'o', color=b.rgb, markersize=7, markeredgecolor='white', markeredgewidth=0.5)
        t, = ax.plot([], [], [], color=b.rgb, linewidth=1.2, alpha=0.7)

        points.append(p)
        trails.append(t)

    # Create an external legend for the orbital plot
    fig.subplots_adjust(right=0.80)
    orbit_legend = fig.legend(
        points,
        [b.name for b in bodies],
        loc='center left',
        bbox_to_anchor=(0.86, 0.5),
        bbox_transform=fig.transFigure,
        frameon=True,
        fontsize=9,
        ncol=1,
        title='Bodies',
        title_fontsize=10,
        borderaxespad=0.5,
        labelspacing=0.3,
        handletextpad=0.4
    )
    legend_artists = [orbit_legend] if orbit_legend is not None else []

    # Helper: archive body display assets while removing it from active simulation
    def remove_body(idx):
        # remove the body marker but keep the trail visible
        try:
            points[idx].remove()
        except Exception:
            pass

        try:
            trails[idx].set_alpha(0.55)
            archived_trails.append(trails[idx])
        except Exception:
            pass

        # archive per-body metric lines so the last values remain visible
        try:
            mass_lines[idx].set_alpha(0.55)
            archived_mass_lines.append(mass_lines[idx])
        except Exception:
            pass
        try:
            flux_lines[idx].set_alpha(0.55)
            archived_flux_lines.append(flux_lines[idx])
        except Exception:
            pass
        try:
            temperature_lines[idx].set_alpha(0.55)
            archived_temperature_lines.append(temperature_lines[idx])
        except Exception:
            pass
        try:
            force_lines[idx].set_alpha(0.55)
            archived_force_lines.append(force_lines[idx])
        except Exception:
            pass

        # archive per-body energy line
        try:
            eidx = idx + 1
            if 0 <= eidx < len(energy_lines):
                energy_lines[eidx].set_alpha(0.55)
                archived_energy_lines.append(energy_lines[eidx])
        except Exception:
            pass

        # remove data arrays and active artist references to keep indexing aligned
        try:
            trails_x.pop(idx)
            trails_y.pop(idx)
            trails_z.pop(idx)
        except Exception:
            pass
        try:
            points.pop(idx)
            trails.pop(idx)
        except Exception:
            pass
        try:
            mass_lines.pop(idx); mass_data.pop(idx)
        except Exception:
            pass
        try:
            flux_lines.pop(idx); flux_data.pop(idx)
        except Exception:
            pass
        try:
            temperature_lines.pop(idx); temperature_data.pop(idx)
        except Exception:
            pass
        try:
            force_lines.pop(idx)
        except Exception:
            pass
        try:
            del energy_lines[idx+1]
            del energy_data[idx+1]
        except Exception:
            pass

    info_text = ax.text2D(
        0.02,
        0.98,
        '',
        transform=ax.transAxes,
        fontsize=11,
        verticalalignment='top',
        bbox=dict(facecolor='#1a1f35', alpha=0.3, boxstyle='round,pad=0.6', edgecolor='#444444', linewidth=1),
        color='#e0e0e0',
        family='monospace'
    )

    # Energy plot setup
    initial_energies = calculate_energy()
    energy_lines = []
    energy_data = [[e] for e in initial_energies]
    time_data = [0.0]

    limit = max(abs(e) for e in initial_energies) * 1.5 if initial_energies else 1.0

    for i in range(len(initial_energies)):
        label = "Entire System" if i == 0 else bodies[i-1].name
        color = 'black' if i == 0 else bodies[i-1].rgb
        line, = ax_energy.plot(time_data, energy_data[i], label=label, color=color)
        energy_lines.append(line)

    ax_energy.set_xlim(0, energy_reset_period)
    ax_energy.set_ylim(-limit, limit)
    ax_energy.set_ylabel("Total Mechanical Energy (J)", fontsize=9)
    ax_energy.set_title("Energy Tracking", fontsize=10, fontweight='bold', color='#e0e0e0')


    # Force plot setup
    initial_forces = calculate_forces()
    force_lines = []
    force_data = [[f] for f in initial_forces]
    time_data_f = [0.0]  # Separate time data for force plot

    limit = max(abs(f) for f in initial_forces) * 1.1 if initial_forces else 1.0

    for i in range(len(initial_forces)):
        label = bodies[i].name
        line, = ax_force.plot(time_data_f, force_data[i], label=label, color=bodies[i].rgb)
        force_lines.append(line)

    ax_force.set_xlim(0, energy_reset_period)
    ax_force.set_ylim(0, limit)
    ax_force.set_xlabel("Time (days)", fontsize=9)
    ax_force.set_ylabel("Gravitational Force (N)", fontsize=9)
    ax_force.set_title("Force Analysis", fontsize=10, fontweight='bold', color='#e0e0e0')

    # Mass plot setup
    mass_lines = []
    mass_data = [[b.mass] for b in bodies]
    time_data_mass = [0.0]

    initial_masses = [b.mass for b in bodies]
    limit = max(initial_masses) * 1.5 if initial_masses else 1.0

    for i, b in enumerate(bodies):
        line, = ax_mass.plot(time_data_mass, mass_data[i], label=b.name, color=b.rgb)
        mass_lines.append(line)

    ax_mass.set_xlim(0, energy_reset_period)
    ax_mass.set_ylim(0, limit)
    ax_mass.set_ylabel("Mass (kg)", fontsize=9)
    ax_mass.set_title("Mass Evolution", fontsize=10, fontweight='bold', color='#e0e0e0')

    # Flux plot setup
    flux_lines = []
    flux_data = [[b.flux] for b in bodies]
    time_data_flux = [0.0]

    initial_fluxes = [b.flux for b in bodies]
    max_flux = max(initial_fluxes) if initial_fluxes else 1.0
    limit = max(max_flux * 1.5, 1.0)  # Set a minimum limit to avoid zero range

    for i, b in enumerate(bodies):
        line, = ax_flux.plot(time_data_flux, flux_data[i], color=b.rgb) #, label=b.name
        flux_lines.append(line)

    ax_flux.set_xlim(0, energy_reset_period)
    ax_flux.set_ylim(0, limit)
    ax_flux.set_ylabel("Flux (W/m²)", fontsize=9)
    ax_flux.set_title("Radiation Flux", fontsize=10, fontweight='bold', color='#e0e0e0')
    #ax_flux.legend()

    # Temperature plot setup
    temperature_lines = []
    temperature_data = [[b.temperature] for b in bodies]
    time_data_temperature = [0.0]

    initial_temperatures = [b.temperature for b in bodies]
    limit = max(initial_temperatures) * 1.5 if initial_temperatures else 1.0

    for i, b in enumerate(bodies):
        line, = ax_temperature.plot(time_data_temperature, temperature_data[i], color=b.rgb) #, label=b.name
        temperature_lines.append(line)

    ax_temperature.set_xlim(0, energy_reset_period)
    ax_temperature.set_ylim(0, limit)
    ax_temperature.set_xlabel("Time (days)", fontsize=9)
    ax_temperature.set_ylabel("Temperature (K)", fontsize=9)
    ax_temperature.set_title("Thermal State", fontsize=10, fontweight='bold', color='#e0e0e0')
    #ax_temperature.legend()

    # Track simulation time for adaptive timestep and plotting
    simulation_time = [0.0]
    last_reset_time_energy = [0.0]
    last_reset_time_force = [0.0]
    last_reset_time_mass = [0.0]
    last_reset_time_flux = [0.0]
    last_reset_time_temperature = [0.0]
    collisions = []  # Track collision events: (frame, time, body1, body2)
    statuses = []  # List to track collision events as tuples: (frame, time, message)

    def frame_generator():
        while simulation_time[0] < Δt:
            yield None

    def init():
        return points + trails + archived_trails + energy_lines + archived_energy_lines + force_lines + archived_force_lines + mass_lines + archived_mass_lines + flux_lines + archived_flux_lines + temperature_lines + archived_temperature_lines + legend_artists

    def update(frame):
        Update_()
        simulation_time[0] += δt
        current_time = simulation_time[0]

        for i, b in enumerate(bodies):

            trails_x[i].append(b.coordinates[0])
            trails_y[i].append(b.coordinates[1])
            trails_z[i].append(b.coordinates[2])

            if len(trails_x[i]) > MAX_TRAIL:
                trails_x[i].pop(0)
                trails_y[i].pop(0)
                trails_z[i].pop(0)

            points[i].set_data_3d([b.coordinates[0]], [b.coordinates[1]], [b.coordinates[2]])
            trails[i].set_data_3d(trails_x[i], trails_y[i], trails_z[i])

        # Update energy data
        energies = calculate_energy()
        forces = calculate_forces()

        # Use the actual adaptive timestep clock
        current_time = simulation_time[0]

        #-------------------------------------------- Reset energy graph if it exceeds the reset period --------------------------------------------
        time_since_reset_energy = current_time - last_reset_time_energy[0]
        if time_since_reset_energy >= energy_reset_period:
            last_reset_time_energy[0] = current_time
            energy_data[:] = [[] for _ in energies]
            time_data[:] = []
            time_since_reset_energy = 0

        #-------------------------------------------- Reset force graph if it exceeds the reset period --------------------------------------------
        time_since_reset_force = current_time - last_reset_time_force[0]
        if time_since_reset_force >= energy_reset_period:
            last_reset_time_force[0] = current_time
            force_data[:] = [[] for _ in forces]
            time_data_f[:] = []
            time_since_reset_force = 0

        time_data.append(time_since_reset_energy); time_data_f.append(time_since_reset_force)
        
        for energy_list, line, energy in zip(energy_data, energy_lines, energies):
            energy_list.append(energy)
            line.set_data(time_data, energy_list)
        
        for force_list, line, force in zip(force_data, force_lines, forces):
            force_list.append(force)
            line.set_data(time_data_f, force_list)

        #-------------------------------------------- Handle mass, flux, and temperature resets and updates --------------------------------------------
        time_since_reset_mass = current_time - last_reset_time_mass[0]
        if time_since_reset_mass >= energy_reset_period:
            last_reset_time_mass[0] = current_time
            mass_data[:] = [[] for _ in bodies]
            time_data_mass[:] = []
            time_since_reset_mass = 0

        time_since_reset_flux = current_time - last_reset_time_flux[0]
        if time_since_reset_flux >= energy_reset_period:
            last_reset_time_flux[0] = current_time
            flux_data[:] = [[] for _ in bodies]
            time_data_flux[:] = []
            time_since_reset_flux = 0

        time_since_reset_temperature = current_time - last_reset_time_temperature[0]
        if time_since_reset_temperature >= energy_reset_period:
            last_reset_time_temperature[0] = current_time
            temperature_data[:] = [[] for _ in bodies]
            time_data_temperature[:] = []
            time_since_reset_temperature = 0

        time_data_mass.append(time_since_reset_mass)
        time_data_flux.append(time_since_reset_flux)
        time_data_temperature.append(time_since_reset_temperature)

        for mass_list, line, b in zip(mass_data, mass_lines, bodies):
            mass_list.append(b.mass)
            line.set_data(time_data_mass, mass_list)

        for flux_list, line, b in zip(flux_data, flux_lines, bodies):
            flux_list.append(b.flux)
            line.set_data(time_data_flux, flux_list)

        for temp_list, line, b in zip(temperature_data, temperature_lines, bodies):
            temp_list.append(b.temperature)
            line.set_data(time_data_temperature, temp_list)

        #-------------------------------------------- Build info text with collision history --------------------------------------------
        info_lines = [f"Time passed: {current_time:.1f} days | {current_time / 365.25:.1f} years", f"Timestep: {round(δt, 4)} days"]
        
        # Show last 3 collisions
        recent_collisions = collisions[-3:] if collisions else []
        if recent_collisions:
            info_lines.append("\nRecent collisions:")
            for _, col_time, col_msg in recent_collisions:
                info_lines.append(f"  {col_msg} @ {col_time:.1f}d")

        # Show last 3 status updates
        recent_statuses = statuses[-3:] if statuses else []
        if recent_statuses:
            info_lines.append("\nRecent status updates:")
            for _, status_time, status_msg in recent_statuses:
                info_lines.append(f"  {status_msg} @ {status_time:.1f}d")
        
        info_text.set_text("\n".join(info_lines))

        ax_energy.set_xlim(0, energy_reset_period)
        ax_force.set_xlim(0, energy_reset_period)
        ax_mass.set_xlim(0, energy_reset_period)
        ax_flux.set_xlim(0, energy_reset_period)
        ax_temperature.set_xlim(0, energy_reset_period)

        #-------------------------------------------- Dynamic y-limit updates based on current data (per-body) --------------------------------------------
        # Choose a body index to focus on (clamp to available bodies)
        body_index = 1 #max(0, min(3, len(bodies) - 1))

        # Mass axis: compare against the focused body's own history
        """
        if mass_data and any(mass_data):
            body_mass_indices = [i for i, b in enumerate(bodies) if b.type in ("Planet", "Moon", "Asteroid")]
            max_mass = max(max(mass_data[i]) for i in body_mass_indices if mass_data[i])
            min_mass = min(min(mass_data[i]) for i in body_mass_indices if mass_data[i])
            ax_mass.set_ylim(min_mass * 1.2, -min_mass * 0.3 + max_mass * 1.1)
        """
        if mass_data and len(mass_data) > body_index:
            body_mass_history = mass_data[body_index]
            focused_mass = bodies[body_index].mass
            if body_mass_history:
                current_body_max = max(body_mass_history)
            else:
                current_body_max = 0
            if focused_mass >= current_body_max:
                ax_mass.set_ylim(bodies[body_index].mass - 10**14, focused_mass * 1.0 + 10**14)
        

        # Energy axis: energy_data[0] is whole-system; per-body energy is at index body_index + 1
        eidx = body_index + 1

        if energy_data and any(energy_data[1:]):
            body_energy_indices = [i + 1 for i, b in enumerate(bodies) if b.type in ("Planet", "Moon", "Asteroid")]
            max_energy = max(max(energy_data[i]) for i in body_energy_indices if energy_data[i])
            min_energy = min(min(energy_data[i]) for i in body_energy_indices if energy_data[i])
            ax_energy.set_ylim(min_energy * 1.2, -min_energy * 0.3 + max_energy * 1.1)
        """
        if energy_data and len(energy_data) > eidx:
            body_energy_history = energy_data[eidx]
            focused_energy = bodies[body_index].mechanical_energy
            if body_energy_history:
                current_body_max_e = max(body_energy_history)
                current_body_min_e = min(body_energy_history)
                current_extreme = max(abs(current_body_min_e), abs(current_body_max_e))
            else:
                current_extreme = 0
            if abs(focused_energy) >= current_extreme:
                max_energy = abs(focused_energy)
                min_energy = -max_energy
                ax_energy.set_ylim(min_energy * 1.2, max_energy * 1.2)
        """

        # Force axis: compare against the focused body's own force history
        if force_data and any(force_data):
            max_force = max(max(i) for i in force_data if i)
            min_force = min(min(i) for i in force_data if i)
            ax_force.set_ylim(0, -min_force * 0.3 + max_force * 1.1)

        """
        if force_data and len(force_data) > body_index:
            body_force_history = force_data[body_index]
            focused_force = bodies[body_index].force
            if body_force_history:
                current_body_max_f = max(body_force_history)
            else:
                current_body_max_f = 0
            if focused_force >= current_body_max_f:
                ax_force.set_ylim(0, focused_force * 1.2)
        """

        # Flux axis
        if flux_data and any(flux_data):
            max_flux = max(max(i) for i in flux_data if i)
            min_flux = min(min(i) for i in flux_data if i)
            ax_flux.set_ylim(0, -min_flux * 0.3 + max_flux * 1.1)
        """
        if flux_data and len(flux_data) > body_index:
            body_flux_history = flux_data[body_index]
            focused_flux = bodies[body_index].flux
            if body_flux_history:
                current_body_max_flux = max(body_flux_history)
            else:
                current_body_max_flux = 0
            if focused_flux >= current_body_max_flux:
                ax_flux.set_ylim(0, focused_flux * 1.2)
        """

        # Temperature axis
        if temperature_data and any(temperature_data):
            body_temp_indices = [i for i, b in enumerate(bodies) if b.type in ("Planet", "Moon", "Asteroid")]
            max_temp = max(max(temperature_data[i]) for i in body_temp_indices if temperature_data[i])
            min_temp = min(min(temperature_data[i]) for i in body_temp_indices if temperature_data[i])
            ax_temperature.set_ylim(min_temp * 1.2, -min_temp * 0.3 + max_temp * 1.1)
        """    
        if temperature_data and any(temperature_data):
            max_temp = max(max(i) for i in temperature_data if i)
            min_temp = min(min(i) for i in temperature_data if i)
            ax_temperature.set_ylim(min_temp * 1.2, max_temp * 1.2)
        """
        """
        if temperature_data and len(temperature_data) > body_index:
            body_temp_history = temperature_data[body_index]
            focused_temp = bodies[body_index].temperature
            if body_temp_history:
                current_body_max_temp = max(body_temp_history)
            else:
                current_body_max_temp = 0
            if focused_temp >= current_body_max_temp:
                ax_temperature.set_ylim(0, focused_temp * 1.2)
        """

        #-------------------------------------------- Save the current figure automatically before each reset period. --------------------------------------------

        # The following piece of code allows the Plot to be automatically centered to the center coordinates of all bodies in the simulation
        if auto_center == True:
            while current_time >= next_save_time_orbits[0] and next_save_time_orbits[0] <= Δt:
                po = np.array([body.coordinates for body in bodies])
                center = np.mean(po, axis=0)
                                    
                max_r = np.max(np.abs(po - center))
                limit = max(max_r * 1.0, 1e10)

                ax.set_xlim(-limit + center[0], limit + center[0])
                ax.set_ylim(-limit + center[1], limit + center[1])
                ax.set_zlim(-limit + center[2], limit + center[2]) 
                next_save_time_orbits[0] += energy_reset_period/k_update

        while current_time >= next_save_time[0] and next_save_time[0] <= Δt:
            save_png_path = save_png_template.format(num=save_count[0], current_time=current_time)
            save_json_path = save_json_template.format(num=save_count[0], current_time=current_time)
            try:
                fig.savefig(save_png_path, dpi=200, bbox_inches='tight')
                print(f"Saved plot automatically to: {save_png_path}")
            except Exception as e:
                print(f"Failed to save pre-reset plot: {e}")

            json_data = {
                "save_index": save_count[0],
                "save_time": current_time,
                "reset_threshold": next_save_time[0],
                "\nSystem's Name:": Name_sys,
                "\nentire_system_energy": {
                    "times": list(time_data),
                    "values": list(energy_data[0])
                },
                "bodies": [
                    {
                        "name": bodies[i].name,
                        "times": list(time_data),
                        "energy": list(energy_data[i + 1]) if i + 1 < len(energy_data) else [],
                        "force": list(force_data[i]) if i < len(force_data) else [],
                        "mass": list(mass_data[i]) if i < len(mass_data) else [],
                        "flux": list(flux_data[i]) if i < len(flux_data) else [],
                        "temperature": list(temperature_data[i]) if i < len(temperature_data) else []
                    }
                    for i in range(len(bodies))
                ]
            }

            try:
                with builtins.open(save_json_path, "w", encoding="utf-8") as json_file:
                    json.dump(json_data, json_file, indent=2)
                print(f"Saved plot data JSON automatically to: {save_json_path}")
            except Exception as e:
                print(f"Failed to save pre-reset JSON data: {e}")

            txt_lines = []
            txt_lines.append(f"Save index: {save_count[0]}")
            txt_lines.append(f"Current time: {current_time:.4f} days")
            txt_lines.append(f"Reset threshold: {next_save_time[0]:.4f} days")
            txt_lines.append("")
            for i, body in enumerate(bodies):
                txt_lines.append(f"Body: {body.name}")
                txt_lines.append("Time (days) | Energy (J)        | Force (N)      | Mass (kg)      | Flux (W/m^2)   | Temperature (K)")
                txt_lines.append("-" * 110)

                body_energy = list(energy_data[i + 1]) if i + 1 < len(energy_data) else []
                body_force = list(force_data[i]) if i < len(force_data) else []
                body_mass = list(mass_data[i]) if i < len(mass_data) else []
                body_flux = list(flux_data[i]) if i < len(flux_data) else []
                body_temp = list(temperature_data[i]) if i < len(temperature_data) else []
                row_count = max(
                    len(time_data),
                    len(body_force),
                    len(body_mass),
                    len(body_flux),
                    len(body_temp)
                )

                for idx in range(row_count):
                    t = f"{time_data[idx]:.4f}" if idx < len(time_data) else ""
                    e = f"{body_energy[idx]:.6e}" if idx < len(body_energy) else ""
                    f_val = f"{body_force[idx]:.6e}" if idx < len(body_force) else ""
                    m = f"{body_mass[idx]:.6e}" if idx < len(body_mass) else ""
                    fl = f"{body_flux[idx]:.6e}" if idx < len(body_flux) else ""
                    temp = f"{body_temp[idx]:.6e}" if idx < len(body_temp) else ""
                    txt_lines.append(f"{t:<12} | {e:<16} | {f_val:<14} | {m:<14} | {fl:<14} | {temp:<14}")
                txt_lines.append("")

            save_txt_path = save_json_path[:-5] + ".txt"
            try:
                with builtins.open(save_txt_path, "w", encoding="utf-8") as txt_file:
                    txt_file.write("\n".join(txt_lines))
                print(f"Saved plot data TXT automatically to: {save_txt_path}")
            except Exception as e:
                print(f"Failed to save pre-reset TXT data: {e}")

            save_count[0] += 1
            next_save_time[0] += energy_reset_period

        #-------------------------------------------- Check for collisions between all body pairs --------------------------------------------
        i = 0
        while i < len(bodies):
            collision_occurred = False
            j = i + 1
            while j < len(bodies):
                dx = bodies[j].coordinates[0] - bodies[i].coordinates[0]
                dy = bodies[j].coordinates[1] - bodies[i].coordinates[1]
                dz = bodies[j].coordinates[2] - bodies[i].coordinates[2]
                r = sqrt(dx**2 + dy**2 + dz**2)
                
                if r < (bodies[i].radii + bodies[j].radii):
                    collision_msg = f"Collision: {bodies[i].name} & {bodies[j].name}"
                    collisions.append((frame, current_time, collision_msg))
                    print(f"Collision detected between {bodies[i].name} and {bodies[j].name} at time {current_time:.2f} days!")
                    # Simple collision response: merge bodies (conservation of mass and momentum)
                    total_mass = bodies[i].mass + bodies[j].mass
                    new_velocity_x = (bodies[i].velocity[0] * bodies[i].mass + bodies[j].velocity[0] * bodies[j].mass) / total_mass
                    new_velocity_y = (bodies[i].velocity[1] * bodies[i].mass + bodies[j].velocity[1] * bodies[j].mass) / total_mass
                    new_velocity_z = (bodies[i].velocity[2] * bodies[i].mass + bodies[j].velocity[2] * bodies[j].mass) / total_mass
                    new_coordinates_x = (bodies[i].coordinates[0] * bodies[i].mass + bodies[j].coordinates[0] * bodies[j].mass) / total_mass
                    new_coordinates_y = (bodies[i].coordinates[1] * bodies[i].mass + bodies[j].coordinates[1] * bodies[j].mass) / total_mass
                    new_coordinates_z = (bodies[i].coordinates[2] * bodies[i].mass + bodies[j].coordinates[2] * bodies[j].mass) / total_mass
                            
                    # Update body i to be the merged body, and remove body j
                    bodies[i].mass = total_mass
                    bodies[i].velocity = [new_velocity_x, new_velocity_y, new_velocity_z]
                    bodies[i].coordinates = [new_coordinates_x, new_coordinates_y, new_coordinates_z]
                    remove_body(j)
                    bodies.pop(j)
                    collision_occurred = True
                    break
                j += 1

            if collision_occurred:
                continue

            if bodies[i].mass <= 0 or bodies[i].status in ["Torn Apart", "Reached Roche Limit", "Extreme Mass Loss"]:
                print(f"{bodies[i].name} was lost at {current_time:.2f} days")
                status_msg = f"{bodies[i].name} lost: {bodies[i].status}"
                statuses.append((frame, current_time, status_msg))
                remove_body(i)
                bodies.pop(i)
                continue

            i += 1

        return points + trails + archived_trails + energy_lines + archived_energy_lines + force_lines + archived_force_lines + mass_lines + archived_mass_lines + flux_lines + archived_flux_lines + temperature_lines + archived_temperature_lines + [info_text] + legend_artists
    

    ani = FuncAnimation(
        fig,
        update,
        frames=frame_generator(),
        init_func=init,
        interval=2,
        blit=False,
        repeat=False
    )

    plt.show()

# ---------------- RUN PROGRAM ----------------

Visualize()