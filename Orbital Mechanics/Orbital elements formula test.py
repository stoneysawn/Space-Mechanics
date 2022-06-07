#Muhammad Yunus bin Shamsul Bahari
#204686
import numpy as np
degree = chr(176)
u = 398600

r = np.array([-6045,-3490,2500])
v = np.array([-3.457,6.618,2.533])
K = np.array([0,0,1])

distance = np.sqrt(np.dot(r,r))
speed = np.sqrt(np.dot(v,v))
rvelocity = np.dot(r,(v/distance))
h = np.cross(r,v)
magh = np.sqrt(np.dot(h,h))
inclination = np.rad2deg(np.arccos(h[2]/magh))
N = np.cross(K,h)
magN = np.sqrt(np.dot(N,N))
RAAN = 0
AOP = 0
true_anomaly = 0

if N[1]>0:
    RAAN = np.arccos(N[0]/magN)
else:
    RAAN = 360-np.rad2deg(np.arccos(N[0]/magN))

e = (1/u)*(((speed**2-u/distance))*r-distance*rvelocity*v)
eccentricity = np.sqrt(np.dot(e,e))

if e[2]>0:
    aop = np.rad2deg(np.arccos(np.dot(N,e)/(magN*eccentricity)))
else:
    aop = 360 - (np.rad2deg(np.arccos(np.dot(N,e)/(magN*eccentricity))))

if rvelocity >= 0:
    true_anomaly = np.rad2deg(np.arccos((np.dot(e,r)/(eccentricity*distance))))
else:
    true_anomaly = 360 - (np.rad2deg(np.arccos((np.dot(e,r)/(eccentricity*distance)))))

print('Distance = {:.4f}km'.format(distance))
print('Speed = {:.4f}km/s '.format(speed))
print('Radial velocity = {:.4f}km/s'.format(rvelocity))
print('\nhx = {:.4f}km^2/s\nhy = {:.4f}km^2/s\nhz = {:.4f}km^2/s'.format(h[0],h[1],h[2]))
print('\nMagnitude h = {:.4f}km^2/s'.format(magh))
print('Inclination = {:.4f}{}'.format(inclination,degree))
print('N = {}km^2/s'.format(N))
print('Magnitude N = {:.4f}'.format(magN))
print('RAAN = {:.4f}'.format(RAAN))
print('e = {}'.format(e))
print('Eccentricity = {:.4f}'.format(eccentricity))
print('Argument of perigee = {:.2f}{}'.format(aop,degree))
print('True anomaly = {:.2f}{}'.format(true_anomaly,degree))

