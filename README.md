# Satellite Dynamics in 6-DOF

This repo contains code for simulating the dynamics of a satellite in 6-DOF.

## Problem Setup

Assuming we are on Earth with the following planet parameters:

- Mass: $M = 5.97219 × 10^{24} kg$
- Radius: $R = 6.3781 × 10^6 m$
- Gravitational constant: $G = 6.67430 × 10^{-11} m^3 kg^{-1} s^{-2}$

The satellite is in a circular orbit with the following parameters:

- Mass: $m = 8kg$
- Altitude: $h = 5 x 10^5 m$
- Inclination: $i = 0$

The satellite is equipped with the following sensors:

- IMU

The satellite is equipped with the following actuators:

- Reaction wheels
- Magnetic torquers
- Thrusters

Forces acting on the satellite:

- Gravitational force: $\vec{F_g}$
- ~~Aerodynamic drag: $\vec{F_d}$~~
- ~~Solar radiation pressure: $\vec{F_s}$~~
- ~~Thrust: $\vec{F_t}$~~

Torques acting on the satellite:

- ~~Magnetic torque: $\vec{T_m}$~~
- ~~Reaction wheel torque: $\vec{T_r}$~~
- ~~Thruster torque: $\vec{T_t}$~~

## Satellite Translational Dynamics

Translational dynamics define the motion of the satellite without considering the rotation. The states considered here are just the position and velocity of the satellite in the Earth Centered Inertial (ECI) frame. Let the vector from the center of the Earth to the satellite's center of mass be $\vec{r} = [x, y, z]^T$ and the velocity of the satellite be $\vec{v} = [u, v, w]^T$. Then,

$$
\dot{\vec{r}} = \vec{v}
$$

The acceleration of the satellite $\vec{a}=\dot{\vec{v}}$ is given by:

$$
\dot{\vec{v}} = -\frac{\vec{F}}{m}
$$

The total force acting on the satellite is given by:

$$
\vec{F} = \vec{F_g} + \vec{F_d} + \vec{F_s} + \vec{F_t}
$$

The gravitational force vector acting on the satellite is given by:

$$
\vec{F_g} = -\frac{GMm}{||\vec{r}||^2}\hat{r}
$$

where $\hat{r}$ is the unit vector in the direction of $\vec{r}$.

The aerodynamic drag force vector acting on the satellite is given by:

$$
\vec{F_d} = -\frac{1}{2}\rho v^2 C_d A \hat{v}
$$

where $\rho$ is the air density, $v$ is the velocity of the satellite, $C_d$ is the drag coefficient, and $A$ is the cross-sectional area of the satellite.

The solar radiation pressure force vector acting on the satellite is given by:

$$
\vec{F_s} = \frac{P_s}{c}A \hat{s}
$$

where $P_s$ is the solar radiation pressure, $c$ is the speed of light, and $\hat{s}$ is the unit vector in the direction of the sun.

To fly in orbit, the satellite needs to have a translational velocity that is perpendicular to the gravitational force acting on it. This implies that the centripetal force acting on the satellite is equal to the gravitational force acting on it. The centripetal force is given by:

$$
\vec{F_c} = \frac{m v^2}{||\vec{r}||}\hat{r}
$$

Equating the gravitational force and the centripetal force, we get:

$$
\begin{aligned}
-\frac{GMm}{||\vec{r}||^2}\hat{r} &= \frac{m v^2}{||\vec{r}||}\hat{r} \\
\frac{GM}{||\vec{r}||} &= v^2 \\
v &= \sqrt{\frac{GM}{||\vec{r}||}}
\end{aligned}
$$

We can use this velocity to initialize the satellite in orbit. The $x, y, z$ components of the velocity are related to the inclination of the satellite orbit. We initialize the satellite in a circular orbit, placing it along the $x$-axis. Hence, it will have 0 velocity in the $x$-direction ($p=0$) and the velocity in the $y$-direction will be given by:

$$
\begin{aligned}
q = vcos(i) \\
r = vsin(i)
\end{aligned}
$$

## Satellite Rotational Dynamics

We use quaternions to represent the angular orientation of the satellite. The quaternion $\vec{q} = [q_0, q_1, q_2, q_3]^T$ represents the rotation of the satellite from the ECI frame to the body frame. The angular velocity of the satellite is given by $\vec{\omega} = [p, q, r]^T$. The quaternion derivative is given by:

$$
\dot{\vec{q}} = \frac{1}{2}\begin{bmatrix}
q_4 & -q_3 & q_2 \\
q_3 & q_4 & -q_1 \\
-q_2 & q_1 & q_4 \\
-q_1 & -q_2 & -q_3
\end{bmatrix}\begin{bmatrix}
p \\
q \\
r
\end{bmatrix}
$$

The total qtorque acting on the satellite is given by:

$$
\vec{T} = \vec{T_m} + \vec{T_r} + \vec{T_t}
$$
