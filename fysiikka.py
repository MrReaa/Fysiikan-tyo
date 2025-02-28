import pandas as pd
import numpy as np
from math import radians, cos, sin, asin, sqrt
from scipy.signal import butter, filtfilt, find_peaks
import streamlit as st
import matplotlib.pyplot as plt
import folium


df = pd.read_csv("https://raw.githubusercontent.com/MrReaa/Fysiikan-tyo/refs/heads/main/Linear%20Acceleration.csv")
dm = pd.read_csv('https://raw.githubusercontent.com/MrReaa/Fysiikan-tyo/refs/heads/main/Location.csv')
#Suodatettu data

def butter_lowpass_filter(data, cutoff, fs, nyq, order):
    normal_cutoff = cutoff / nyq
    b, a = butter(order, normal_cutoff, btype='low', analog=False)
    y = filtfilt(b, a, data)
    return y

def butter_highpass_filter(data, cutoff, fs, nyq, order):
    normal_cutoff = cutoff / nyq
    b, a = butter(order, normal_cutoff, btype='high', analog=False)
    y = filtfilt(b, a, data)
    return y


#Filttereiden parametrit
T = df['Time (s)'][len(df['Time (s)'])-1] - df['Time (s)'][0]  # Koko datan pituus
n = len(df['Time (s)'])  # Datapisteiden lukumäärä
fs = n / T  # Näytteenottotaajuus
nyq = fs / 2  # Nyqvistin taajuus
order = 3  # Kertaluku
cutoff = 1 / 0.2  # Cutt-off taajuus


df['filter_a_x'] = butter_lowpass_filter(df['Linear Acceleration x (m/s^2)'], cutoff, fs, nyq, order)
df['filter_a_y'] = butter_lowpass_filter(df['Linear Acceleration y (m/s^2)'], cutoff, fs, nyq, order)
df['filter_a_z'] = butter_lowpass_filter(df['Linear Acceleration z (m/s^2)'], cutoff, fs, nyq, order)


#Askelmäärä laskettuna suodatetusta kiihtyvyysdatasta
df['acc_magnitude'] = np.sqrt(df['filter_a_x'] ** 2 + df['filter_a_y'] ** 2 + df['filter_a_z'] ** 2)
peaks, _ = find_peaks(df['acc_magnitude'], height=0.5, distance=10)
step_count = len(peaks)


# Kuljettu matka (GPS-datasta) ja keskinopeus (GPS-datasta)

def haversine(lon1, lat1, lon2, lat2):
    lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2])
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = sin(dlat / 2) ** 2 + cos(lat1) * cos(lat2) * sin(dlon / 2) ** 2
    c = 2 * asin(sqrt(a))
    r = 6371  # Radius of Earth in km
    return c * r  # Distance in km


dm['dist'] = np.zeros(len(dm))
dm['time_diff'] = np.zeros(len(dm))

for i in range(len(dm) - 1):
    dm.loc[i, 'dist'] = haversine(dm['Longitude (°)'][i], dm['Latitude (°)'][i], dm['Longitude (°)'][i + 1], dm['Latitude (°)'][i + 1]) * 1000  # meters
    dm.loc[i, 'time_diff'] = dm['Time (s)'][i + 1] - dm['Time (s)'][i]

total_distance = dm['dist'].sum()
total_time = dm['time_diff'].sum()
average_speed = total_distance / total_time


# Askelpituus (lasketun askelmäärän ja matkan perusteella)

step_length = total_distance / step_count



st.title("Fysiikan loppuprojekti")

st.write(f"### Puhelimeni tippui vahingossa matkalla jonka takia mittauksissa piikki")


st.write(f"#### Askelmäärä: {step_count}")
st.write(f"#### Kokonaismatka: {total_distance:.2f} meters")
st.write(f"#### Keskinopeus: {average_speed:.2f} m/s")
st.write(f"#### Askelpituus: {step_length:.2f} meters")


# Suodatettu kiihtyvyysdata x,y,z-komponenttit
st.subheader('Suodatettu kiihtyvyysdata x,y,z-komponenttit')

fig, axs = plt.subplots(3, 1, figsize=(12, 14))
axs[0].plot(df['Time (s)'], df['filter_a_x'])
axs[0].set_title('Filtered Acceleration X')
axs[1].plot(df['Time (s)'], df['filter_a_y'])
axs[1].set_title('Filtered Acceleration Y')
axs[2].plot(df['Time (s)'], df['filter_a_z'])
axs[2].set_title('Filtered Acceleration Z')

st.pyplot(fig)

# Reitti kartalla

st.subheader('Route Map')

start_lat = dm['Latitude (°)'].mean()
start_long = dm['Longitude (°)'].mean()
my_map = folium.Map(location=[start_lat, start_long], zoom_start=14)

folium.PolyLine(dm[['Latitude (°)', 'Longitude (°)']], color='red', weight=2.5, opacity=1).add_to(my_map)

my_map.save("fysiikka_lopputyo.html")
st.markdown("### Interactive Map")
st.components.v1.html(open("fysiikka_lopputyo.html", 'r').read(), height=500)
