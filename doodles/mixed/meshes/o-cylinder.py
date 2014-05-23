import matplotlib.pyplot as plt
from numpy import *

# Basis unit
u = 0.1

# Circle points
O = array([0, 0])
E = 3*u*array([cos(4*pi/4), sin(4*pi/4)])
W = 3*u*array([cos(0*pi/4), sin(0*pi/4)])
N = 3*u*array([cos(2*pi/4), sin(2*pi/4)])
S = 3*u*array([cos(6*pi/4), sin(6*pi/4)])

# Circle arc
x0, y0, r = 0, 0, 3*u
thetas = arange(0, 2*pi, 0.1)
circle_x = x0 + r*cos(thetas)
circle_y = y0 + r*sin(thetas)

# Lower turn points
P = 4*u*array([cos(5*pi/4), sin(5*pi/4)])
L = 4*u*array([cos(6*pi/4), sin(6*pi/4)])
R = 4*u*array([cos(7*pi/4), sin(7*pi/4)])

# Lower arc
thetas = arange(5*pi/4, 7*pi/4, 0.01)
l_arc_x = 4*u*cos(thetas)
l_arc_y = 4*u*sin(thetas)

# Upper turn points
P_ = 4*u*array([cos(3*pi/4), sin(3*pi/4)])
L_ = 4*u*array([cos(2*pi/4), sin(2*pi/4)])
R_ = 4*u*array([cos(1*pi/4), sin(1*pi/4)])

# Upper arc
thetas = arange(1*pi/4, 3*pi/4, 0.01)
u_arc_x = 4*u*cos(thetas)
u_arc_y = 4*u*sin(thetas)

# Lower left arm point
s = 3
Q = P - s*array([-sin(5*pi/4), cos(5*pi/4)])

# Upper left arm point
t = s - 4*u
Q_ = P_ - t*array([-sin(5*pi/4), cos(5*pi/4)])

# Lower right arm point
s = 3
S = R + s*array([-sin(7*pi/4), cos(7*pi/4)])

# Upper right arm point
t = s - 4*u
S_ = R_ + t*array([-sin(7*pi/4), cos(7*pi/4)])

# Plot
plt.figure()

plt.plot(circle_x, circle_y)
plt.plot(l_arc_x, l_arc_y)
plt.plot(u_arc_x, u_arc_y)

plt.plot([Q[0], P[0]], [Q[1], P[1]])
plt.plot([R[0], S[0]], [R[1], S[1]])
plt.plot([Q_[0], P_[0]], [Q_[1], P_[1]])
plt.plot([R_[0], S_[0]], [R_[1], S_[1]])
plt.plot([Q[0], Q_[0]], [Q[1], Q_[1]])
plt.plot([S[0], S_[0]], [S[1], S_[1]])

plt.axis('equal')
plt.show()
