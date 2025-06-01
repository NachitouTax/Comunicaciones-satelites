from datetime import datetime
import pytz
import math

# ===============================
# Constantes
# ===============================
mu = 398600.4418  # km^3/s^2
earth_radius = 6378.137  # km
twopi = 2.0 * math.pi

# ===============================
# Leer archivo de entrada
# ===============================
input_file = "datos.txt"
with open(input_file, "r") as f:
    lines = f.read().strip().splitlines()
    tle_line1 = lines[0].strip()
    tle_line2 = lines[1].strip()
    observer_lat_deg = float(lines[2].strip())
    observer_lon_deg = float(lines[3].strip())

# ===============================
# Obtener hora local
# ===============================
timezone = pytz.timezone('America/Santiago')
local_time = datetime.now(timezone)

# ===============================
# Parseo del TLE
# ===============================
def parse_tle(line1, line2):
    inclination = float(line2[8:16])
    raan = float(line2[17:25])
    eccentricity = float("0." + line2[26:33].strip())
    arg_perigee = float(line2[34:42])
    mean_anomaly = float(line2[43:51])
    mean_motion = float(line2[52:63])
    return {
        "inclination": math.radians(inclination),
        "raan": math.radians(raan),
        "eccentricity": eccentricity,
        "arg_perigee": math.radians(arg_perigee),
        "mean_anomaly": math.radians(mean_anomaly),
        "mean_motion": mean_motion
    }

def kepler_equation(M, e, tol=1e-8):
    E = M if e < 0.8 else math.pi
    while True:
        dE = (E - e * math.sin(E) - M) / (1 - e * math.cos(E))
        E -= dE
        if abs(dE) < tol:
            break
    return E

def orbital_to_eci(elements):
    a = (mu / (twopi * elements["mean_motion"] / 86400)**2)**(1/3)
    M = elements["mean_anomaly"]
    e = elements["eccentricity"]
    E = kepler_equation(M, e)
    nu = 2 * math.atan2(math.sqrt(1 + e) * math.sin(E / 2),
                        math.sqrt(1 - e) * math.cos(E / 2))
    r = a * (1 - e * math.cos(E))
    x_orb = r * math.cos(nu)
    y_orb = r * math.sin(nu)
    cos_O, sin_O = math.cos(elements["raan"]), math.sin(elements["raan"])
    cos_i, sin_i = math.cos(elements["inclination"]), math.sin(elements["inclination"])
    cos_w, sin_w = math.cos(elements["arg_perigee"]), math.sin(elements["arg_perigee"])
    x = (cos_O * cos_w - sin_O * sin_w * cos_i) * x_orb + (-cos_O * sin_w - sin_O * cos_w * cos_i) * y_orb
    y = (sin_O * cos_w + cos_O * sin_w * cos_i) * x_orb + (-sin_O * sin_w + cos_O * cos_w * cos_i) * y_orb
    z = (sin_w * sin_i) * x_orb + (cos_w * sin_i) * y_orb
    return x, y, z

def eci_to_alt_az(eci, observer_lat_deg, observer_lon_deg, dt):
    lat = math.radians(observer_lat_deg)
    lon = math.radians(observer_lon_deg)
    jd = 367 * dt.year - int((7 * (dt.year + int((dt.month + 9) / 12))) / 4) + int((275 * dt.month) / 9) + dt.day + 1721013.5 + (dt.hour + dt.minute / 60 + dt.second / 3600) / 24
    GMST = 280.46061837 + 360.98564736629 * (jd - 2451545.0)
    GMST = math.radians(GMST % 360)
    theta = GMST + lon
    cos_theta, sin_theta = math.cos(theta), math.sin(theta)
    x, y, z = eci
    x_ecef = cos_theta * x + sin_theta * y
    y_ecef = -sin_theta * x + cos_theta * y
    z_ecef = z
    R = earth_radius
    obs_x = R * math.cos(lat) * math.cos(lon)
    obs_y = R * math.cos(lat) * math.sin(lon)
    obs_z = R * math.sin(lat)
    rx = x_ecef - obs_x
    ry = y_ecef - obs_y
    rz = z_ecef - obs_z
    sin_lat = math.sin(lat)
    cos_lat = math.cos(lat)
    top_x = -sin_lat * math.cos(lon) * rx - sin_lat * math.sin(lon) * ry + cos_lat * rz
    top_y = -math.sin(lon) * rx + math.cos(lon) * ry
    top_z = cos_lat * math.cos(lon) * rx + cos_lat * math.sin(lon) * ry + sin_lat * rz
    az = math.atan2(top_y, top_x)
    el = math.asin(top_z / math.sqrt(top_x**2 + top_y**2 + top_z**2))
    return (math.degrees(az) + 360) % 360, math.degrees(el)

# Procesar
orbital_elements = parse_tle(tle_line1, tle_line2)
sat_eci = orbital_to_eci(orbital_elements)
azimut, altitud = eci_to_alt_az(sat_eci, observer_lat_deg, observer_lon_deg, local_time)

# Mostrar resultados
print("\n--- Resultados ---")
print("Hora local:", local_time.strftime("%Y-%m-%d %H:%M:%S"))
print(f"Azimut del satélite: {azimut:.2f}°")
print(f"Altitud del satélite: {altitud:.2f}°")
