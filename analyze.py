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
'''
i = 340000
j = 80000
data = data[i:i+j]
'''

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


# span now holds the signal to extract bits from
# it is a pulse train
# '000' is encoded as '__^___^___^____'  ~120 samples between peaks
# '111' is encoded as '__^_^_^_^_^_^___' ~60 samples between peaks
# '010' is encoded as '__^___^_^_^_____'

#Hysteresis constants
thresh = 300
trig = 150
isHigh = False
# The loop below uses isHigh as a state variable
# to detect peaks with hysterisis

maxPeakIndex = 0
maxPeakValue = 0
lastPeakIndex = 0 # holds the sample index of midpoint of last peak

isHalfPeriod = False
bits = []
pos = []

for i,s in enumerate(span):
	if isHigh:
		# We are on a peak.
		if s > maxPeakValue:
			# New highest value of peak found:
			maxPeakValue = s
			maxPeakIndex = i
		if s < trig:
			# End of peak detected, compute bit value to store
			pulseLen = maxPeakIndex- lastPeakIndex

			if pulseLen > 200:
				# Very long pulse, start of tape. Ignore.
				pass
			elif pulseLen < 90:
				# short pulse, only store bit on every other peak
				if isHalfPeriod:
					# Second short pulse, store bit
					isHalfPeriod = False
					bits += [1]

					pos += [i] # Bit index for debugging
				else:
					# This was the first pulse, wait for the
					# next one to store bit
					isHalfPeriod = True
			else:
				# 90 <= pulselen <= 200
				# This is a long period

				if isHalfPeriod:
					# Safety check, only a half short period
					# was received previously
					print "Warning: out of sync at index:", i
					isHalfPeriod = False
				# long pulse, add bit
				bits += [0]
				pos += [i] # Bit index for debugging

			# Save this peak position for next time
			lastPeakIndex = maxPeakIndex

			# Leave peak state
			isHigh = False
	else:
		# not isHigh, we are between two peaks
		if s > thresh:
			# Peak detected, enter peak state
			maxPeakIndex = i
			maxPeakValue = s
			isHigh = True

# bits now holds the raw bitstream from the audio tape
# pos holds the sample index of the bits

'''
plt.plot(span,'b-')
plt.plot(pos,[e * max(span) for e in bits],'k+')
plt.show()
'''

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

