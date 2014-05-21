import matplotlib.pyplot as plt
from numpy import *

# Basis unit
u = 0.1

P = 4*u*array([cos(5*pi/4), sin(5*pi/4)])
B = 4*u*array([cos(6*pi/4), sin(6*pi/4)])
V = 4*u*array([cos(2*pi/4), sin(2*pi/4)])
R = 4*u*array([cos(7*pi/4), sin(7*pi/4)])

thetas = arange(5*pi/4, 7*pi/4, 0.01)
arc_x = 4*u*cos(thetas)
arc_y = 4*u*sin(thetas)

s = 2
Q = P - s*array([-sin(5*pi/4), cos(5*pi/4)])
S = R + s*array([-sin(7*pi/4), cos(7*pi/4)])

t = s - sin(pi/4)*4*u

A = V - t*array([-sin(5*pi/4), cos(5*pi/4)])
B = V + t*array([-sin(7*pi/4), cos(7*pi/4)])

x0, y0, r = 0, 0, 2*u
thetas = arange(0, 2*pi, 0.1)
circle_x = x0 + r*cos(thetas)
circle_y = y0 + r*sin(thetas)

# Plot
plt.figure()

plt.plot(arc_x, arc_y)

plt.plot([Q[0], P[0]], [Q[1], P[1]])
plt.plot([R[0], S[0]], [R[1], S[1]])
plt.plot([A[0], V[0]], [A[1], V[1]])
plt.plot([B[0], V[0]], [B[1], V[1]])

plt.plot([A[0], Q[0]], [A[1], Q[1]])
plt.plot([B[0], S[0]], [B[1], S[1]])

plt.plot(circle_x, circle_y)

plt.axis('equal')
plt.show()

x = (B-V)/dot(B-V, B-V)
y = (S-B)/dot(S-B, S-B)

print sqrt(dot(S-R, S-R))
print sqrt(dot(A-Q, A-Q)), 4*u + 4*u*cos(pi/4)

print dot(x, y), arccos(dot(x, y))*180/pi
