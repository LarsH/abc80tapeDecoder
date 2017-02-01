import wave
from matplotlib import pyplot as plt
import struct
import math

# Read audio file
f = wave.open('hello.bac.wav')
l = f.getnframes()
data = [None]*l
frames = f.readframes(l)
for i in range(l):
	s = frames[2*i:2*i+2]
	data[i] = struct.unpack('<h',s)[0]

print "Got data"

# Only look at the actual basic file
i = 340000
j = 80000
data = data[i:i+j]

# Normalize level
data = [e + 1700 for e in data]


# Remove ringing
# Ringing has five periods on 38 samples
# amplification is about 1.9
# need zeroes at fs / (38/5)
# n = 1.9 * (cos(pi/7.6) +- i*sin(pi/7.6))
# p(z) = (z-n1)*(z-n2) = z*z -(n1+n2)z + n1*n2
# p(z) = z**2 - 2*1.9**2*cos(pi/7.6)z + 1.9**2

# Compute filter taps
a = 1.0 / 1.9
k1 = 1
k2 = -2*a*a*math.cos(math.pi / 5)
k3  = a*a

# Filter should amplify a constant signal with the factor 1, normalize
A = (k1 + k2 + k3)
k1 /= A
k2 /= A
k3 /= A


filt = [e for e in data]
n = 6
for i in range(n):
	print "Filtering... %u/%u"%(i,n)
	filt = [d1*k1+ d2*k2+ d3*k3 for d1,d2,d3 in zip(filt,filt[1:],filt[2:])]
	# Added mean value filter, to remove high frequency noise
	filt = filt[:1] + [(d1+d2+ d3)/3 for d1,d2,d3 in zip(filt,filt[1:],filt[2:])]


# filt now holds the signal without ringings
#plt.plot(data,'r:')
#plt.plot(filt,'b-')
#plt.show()

print "Starting computations..."
# System naturally halves in 40 samples
# u(t) = A * 0.5 ^ (t/40) = A * exp(t * log(0.5)/40)
# u'(t) =  log(0.5)/40 * u(t)
slope = [0, 0] + [(b - a)/5 for a,b in zip(filt, filt[4:])]
k = math.log(0.5)/40
speed = [d - u*k for d,u in zip(slope, filt)]
force = [0, 0] + [(b - a) for a,b in zip(speed, speed[4:])]
span = [max(force[i:i+10]) - min(force[i:i+10]) for i in range(len(force))]

print "Computations done."

#plt.plot(slope,'r--')
#plt.plot(force,'r:')
#plt.show()


#Hysteresis constants,
thresh = 300
trig = 100

isHigh = False
gotHigh = 0
lastMidpoint = 0

isHalfPeriod = False
bits = []
pos = []
for i,s in enumerate(span):
	if isHigh:
		if s < trig:
			isHigh = False
			midpoint = (gotHigh + i)/2
			pulseLen = midpoint - lastMidpoint
			if pulseLen > 200:
				# start of tape
				pass
			elif pulseLen < 80:
				# short pulse
				if isHalfPeriod:
					isHalfPeriod = False
					bits += [1]
					pos += [i]
					pass
				else:
					isHalfPeriod = True
			else:
				assert not isHalfPeriod, repr(i)
				#long pulse, add bit
				bits += [0]
				pos += [i]

			lastMidpoint = midpoint
			gotLow = i
	else:
		if s > thresh:
			isHigh = True
			gotHigh = i

'''
plt.plot(span,'b-')
plt.plot(pos,[e * max(span) for e in bits],'k+')
plt.show()
'''

# bits now holds the raw bitstream from the audio tape

# bit stream starts with a lot of zeroes, find first 1-bit to synchronize
bits = bits[(bits.index(1)-1)%8:]

# Print out the bits and data
text = ''
for i in range(len(bits)/8):
	l = bits[i*8:i*8+8][::-1]
	c = eval('0b' + ''.join(map(str,l)))
	text += chr(c)
	print l, hex(c), repr(chr(c))
print repr(text)

