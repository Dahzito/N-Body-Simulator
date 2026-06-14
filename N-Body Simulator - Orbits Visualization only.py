from math import*
from sys import*
from os import*
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

δt:float = 3 #days
δt_default:float = 2 #default timestep for normal conditions
Δt:float = 365.25*10 #days

Δt_reset:float = 365*(1/2)

ε:float = 0 #Gravitational softening
Acc_dx:float = 0 
Acc_dy:float = 0
Acc_dz:float = 0
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
            self.prev_acceleration = list(acc)  # Store previous acceleration for velocity Verlet~

            self.force = 0

            self.potential_energy = 0; self.kinetic_energy = 0; self.mechanical_energy = 0

            #Physical Properties 
            self.name = name
            self.radii = radii
            self.mass = mass

            self.temperature = temperature                 # Kelvin
            self.temp_threshold = temp_threshold  # Arbitrary threshold for temperature-based effects (e.g., radiation pressure)

            self.gravity = G * self.mass / self.radii**2  # Surface gravity

            self.luminosity = σ * 4 * pi * self.radii**2 * self.temperature**4
            self.albedo = albedo  # Albedo for energy absorption calculations
            self.flux = 0         # Flux received from other bodies (for energy balance calculations)

            #Other properties
            self.type = type
            # Type of body (e.g., "Star", "Planet", "Asteroid")

            self.status = "Stable"  # Status can be "Stable", "Torn Apart", for Roche's limit effects

            self.wavelength = 0 #Wavelength of the radiation emitted by the body, using Wien's Law
            self.rgb = rgb

            self.gh = gh   # Greenhouse effect factor for planets, where 0 means no greenhouse effect. 
                           # This is a simplified model and does not account for atmospheric composition or other factors that influence greenhouse effects.
            self.emissivity = emissivity  # Surface emissivity for thermal radiation calculations of planets, defaulting to 1 (blackbody) (ranges from 0 to 1)
            self.heat_capacity = 1000 * self.mass  # Simplified heat capacity proportional to mass (J/K), for temperature change calculations

# rgba(190, 154, 36, 0.80)
bodies = []

bodies.append(Body(1.227*M_Sun, [0,0,0], [0, 0, 0], [0, 0, 0], "HD 28109 *", 1.425*R_Sun, 5878, 0.0, 1000, "Star", [1, 1, 0.7, 1], 0, 1))
#The RGB values of the Star are automatically determined depending on their surface temperature, also calculated using the Mass-Luminosity Relation, and Stefan-Boltzmann's law.

bodies.append(Body(6.2*M_Earth, [0,0,0], [0, 88033, 19154], [0.1357*au*(1-0.00710), 0, 0], "HD 28109 b", 2.69*R_Earth, 280, 0.31, 2000, "Planet", [206/255, 104/255, 45/255, 0.8], 0.2, 0.9))

bodies.append(Body(9.20*M_Earth, [0,0,0], [0, 58619, 10819], [0.308*au*(1-0.0039), 0, 0], "HD 28109 c", 4.13*R_Earth, 120, 0.50, 2000, "Planet", [226/255, 83/255, 0/255, 0.8], 0.1, 0.9))

bodies.append(Body(5.681*M_Earth, [0,0,0], [0, 50844, 925.7], [0.411*au*(1-0.0054), 0, 0], "HD 28109 d", 3.25*R_Earth, 120, 0.50, 2000, "Planet", [190/255, 154/255, 36/255, 0.8], 0.0, 0.9))

#Mass, acceleration, velocity, coordinates, name, radii, temperature, albedo, temp_threshold, type, RGB colors (0-1), greenhouse effect (0-1), ε - Surface Emissitivity (0-1);


def Acceleration(dx, dy, dz, num, pos_x=None, pos_y=None, pos_z=None):
    eps_local = min(b.radii for b in bodies) * 1
    ax = 0.0
    ay = 0.0
    az = 0.0

    if pos_x is not None and pos_y is not None and pos_z is not None:
        for j in range(len(bodies)):
            if j == num:
                continue
            rx = pos_x[j] - dx
            ry = pos_y[j] - dy
            rz = pos_z[j] - dz
            r_soft = sqrt(rx**2 + ry**2 + rz**2 + eps_local**2)
            ax += (G * bodies[j].mass * rx) / r_soft**3
            ay += (G * bodies[j].mass * ry) / r_soft**3
            az += (G * bodies[j].mass * rz) / r_soft**3
    else:
        for j, other in enumerate(bodies):
            if j == num:
                continue
            rx = other.coordinates[0] - dx
            ry = other.coordinates[1] - dy
            rz = other.coordinates[2] - dz
            r_soft = sqrt(rx**2 + ry**2 + rz**2 + eps_local**2)
            ax += (G * other.mass * rx) / r_soft**3
            ay += (G * other.mass * ry) / r_soft**3
            az += (G * other.mass * rz) / r_soft**3
    return ax, ay, az

def Update_():
    global δt, ε, force_temp
    radius = []
    dt_seconds = δt * _t  # Convert to seconds

    force_temp_dx = 0; force_temp_dy = 0;
    
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
        0.25 / np.sqrt(a_max),
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
    # Coefficients for the Yoshida 4th order symplectic integrator
    w0 = -1.702414289193
    w1 = 1.3512071919597
    c1 = c4 = w1 / 2.0
    c2 = c3 = (w0 + w1) / 2.0
    d1 = d3 = w1
    d2 = w0

    n_b = len(bodies)
    if n_b > 0:
        pos_x = [body.coordinates[0] for body in bodies]
        pos_y = [body.coordinates[1] for body in bodies]
        pos_z = [body.coordinates[2] for body in bodies]

        vel_x = [body.velocity[0] for body in bodies]
        vel_y = [body.velocity[1] for body in bodies]
        vel_z = [body.velocity[2] for body in bodies]

        def compute_accelerations(px, py, pz):
            axs = [0.0] * n_b
            ays = [0.0] * n_b
            azs = [0.0] * n_b

            for i in range(n_b):
                axs[i], ays[i], azs[i] = Acceleration(px[i], py[i], pz[i], i, px, py, pz)
            return axs, ays, azs

        # Stage 1
        for i in range(n_b):
            pos_x[i] += c1 * vel_x[i] * dt_seconds
            pos_y[i] += c1 * vel_y[i] * dt_seconds
            pos_z[i] += c1 * vel_z[i] * dt_seconds

        axs, ays, azs = compute_accelerations(pos_x, pos_y, pos_z)
        for i in range(n_b):
            vel_x[i] += d1 * dt_seconds * axs[i]
            vel_y[i] += d1 * dt_seconds * ays[i]
            vel_z[i] += d1 * dt_seconds * azs[i]

        # Stage 2
        for i in range(n_b):
            pos_x[i] += c2 * vel_x[i] * dt_seconds
            pos_y[i] += c2 * vel_y[i] * dt_seconds
            pos_z[i] += c2 * vel_z[i] * dt_seconds
        axs, ays, azs = compute_accelerations(pos_x, pos_y, pos_z)
        for i in range(n_b):
            vel_x[i] += d2 * dt_seconds * axs[i]
            vel_y[i] += d2 * dt_seconds * ays[i]
            vel_z[i] += d2 * dt_seconds * azs[i]

        # Stage 3
        for i in range(n_b):
            pos_x[i] += c3 * vel_x[i] * dt_seconds
            pos_y[i] += c3 * vel_y[i] * dt_seconds
            pos_z[i] += c3 * vel_z[i] * dt_seconds
        axs, ays, azs = compute_accelerations(pos_x, pos_y, pos_z)
        for i in range(n_b):
            vel_x[i] += d3 * dt_seconds * axs[i]
            vel_y[i] += d3 * dt_seconds * ays[i]
            vel_z[i] += d3 * dt_seconds * azs[i]

        # Final position update
        for i in range(n_b):
            pos_x[i] += c4 * vel_x[i] * dt_seconds
            pos_y[i] += c4 * vel_y[i] * dt_seconds
            pos_z[i] += c4 * vel_z[i] * dt_seconds

        for i in range(n_b):
            bodies[i].coordinates[0] = pos_x[i]
            bodies[i].coordinates[1] = pos_y[i]
            bodies[i].coordinates[2] = pos_z[i]
            bodies[i].velocity[0] = vel_x[i]
            bodies[i].velocity[1] = vel_y[i]
            bodies[i].velocity[2] = vel_z[i]

#------ Visualization of the N Body system using Matplotlib ------


def Visualize(frame_save_period=Δt_reset):
    plt.rcParams['font.family'] = 'monospace'
    plt.rcParams['font.weight'] = 'normal'
    plt.style.use('dark_background')
    """Visualize the n-body system with energy tracking.
    
    Args:
        energy_reset_period: Reset the energy graph every N days (default: 365)
        force_reset_period: Reset the force graph every N days (default: 365)
    """
    fig = plt.figure(num="N-Body Simulation - Orbits Tracking Only", figsize=(24, 18))
    try:
        fig.canvas.manager.set_window_title("N-Body Simulation - Orbits Tracking Only")
    except Exception:
        pass
    save_template = os.path.join(
        os.path.expanduser("~"),
        "Downloads",
        "Python Plots - 3D N-Body Simulator",
        "N-Body Simulator - Orbits Visualization only",
        "N-Body Simulator - Orbits Visualization - HD 28109 System - Plots",
        "N-Body Simulator, HD 28109 System, Plot {num} @ {current_time:.2f} days.png"
    )
    os.makedirs(os.path.dirname(save_template), exist_ok=True)
    save_count = [1]
    next_save_time = [frame_save_period - 2 * δt]  # Save just before the reset for better visualization of changes
    manager = plt.get_current_fig_manager()

    try:
        manager.window.showMaximized()
    except AttributeError:
        try:
            manager.window.state('zoomed')
        except AttributeError:
            manager.full_screen_toggle()

    gs = fig.add_gridspec(1, 1, height_ratios=[1], width_ratios=[1], hspace=0.35, wspace=0.35)
    ax = fig.add_subplot(1, 1, 1, projection="3d")
    ax.set_box_aspect([1, 1, 1])
    fig.patch.set_facecolor('black')
    ax.set_facecolor('black')
    for axis in [ax.xaxis, ax.yaxis, ax.zaxis]:
        axis._axinfo["grid"]['color'] = (1, 1, 1, 0.3)
        axis._axinfo["grid"]['linewidth'] = 0.15

    ax.xaxis.set_pane_color((0.15, 0.15, 0.15, 0.1)) 
    ax.yaxis.set_pane_color((0.15, 0.15, 0.15, 0.1)) 
    ax.zaxis.set_pane_color((0.15, 0.15, 0.15, 0.1))

    ax.xaxis.label.set_color("#e0e0e067")
    ax.yaxis.label.set_color('#e0e0e067')
    ax.zaxis.label.set_color('#e0e0e067')
    ax.title.set_color('white')
    ax.tick_params(colors='#b0b0b0', labelsize=8)

    ax.set_title("N-Body Simulation: Orbits and Energy Tracking", fontsize=14)

    ax.grid(True, which='both', linestyle='--', linewidth=0.2, alpha=0.1)

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

    for i, b in enumerate(bodies):
        p, = ax.plot([], [], [], 'o', color=b.rgb, markeredgecolor='white', markeredgewidth=0.5, markersize=6)
        t, = ax.plot([], [], [], color=b.rgb, linewidth=1.3, alpha=0.7)

        points.append(p)
        trails.append(t)

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


    info_text = ax.text2D(
        0.02,
        0.98,
        '',
        transform=ax.transAxes,
        fontsize=10,
        verticalalignment='top',
        bbox=dict(facecolor='black', alpha=0.2, boxstyle='round,pad=0.4'),
        color='white'
    )
    info_text.set_animated(True)

    # Track simulation time for adaptive timestep and plotting
    simulation_time = [0.0]

    collisions = []  # Track collision events: (frame, time, body1, body2)
    statuses = []  # List to track collision events as tuples: (frame, time, message)

    def frame_generator():
        while simulation_time[0] < Δt:
            yield None

    def init():
        return points + trails + archived_trails + [info_text]

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

        # Use the actual adaptive timestep clock
        current_time = simulation_time[0]

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


        #-------------------------------------------- Save the current figure automatically 50 days before each reset period. --------------------------------------------
        while current_time >= next_save_time[0] and next_save_time[0] <= Δt:
            save_path = save_template.format(num=save_count[0], current_time=current_time)
            try:
                fig.savefig(save_path, dpi=200, bbox_inches='tight')
                print(f"Saved plot automatically to: {save_path}")
            except Exception as e:
                print(f"Failed to save pre-reset plot: {e}")
            save_count[0] += 1
            next_save_time[0] += frame_save_period

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

        return points + trails + archived_trails + [info_text]
    

    ani = FuncAnimation(
        fig,
        update,
        frames=frame_generator(),
        init_func=init,
        interval=2,
        blit=True,
        repeat=False
    )

    plt.show()

# ---------------- RUN PROGRAM ----------------

Visualize()