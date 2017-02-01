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

c = 0
pulses = []
# Simulation of RC filter
for inVolt in data:
	out = inVolt + c
	c -= out*0.1
	pulses += [abs(out)]

pulses = [max(pulses[i:i+20]) for i in range(len(pulses))]
pulses = [min(pulses[i:i+15]) for i in range(len(pulses))]


# pulses now holds the signal to extract bits from
# it is a pulse train
# '000' is encoded as '__^___^___^____'  ~120 samples between peaks
# '111' is encoded as '__^_^_^_^_^_^___' ~60 samples between peaks
# '010' is encoded as '__^___^_^_^_____'

#Hysteresis constants
thresh = 1200
isHigh = False
# The loop below uses isHigh as a state variable
# to detect peaks with hysterisis

maxPeakIndex = 0
peakValue = 0
lastPeakIndex = 0 # holds the sample index of midpoint of last peak

isHalfPeriod = False
bits = []
pos = []

# For debug plotting peak detection
peakPlot = []
peakPos = []

for i,s in enumerate(pulses):
	if isHigh:
		# We are on a peak.
		if s > peakValue:
			# New highest value of peak found:
			peakValue = s
			maxPeakIndex = i
		if s < peakValue - thresh:
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
			peakPlot += [1, 0]
			peakPos += [i, i]
	else:
		# not isHigh, we are between two peaks
		if s < peakValue:
			# New low value found:
			peakValue = s
		if s > thresh + peakValue:
			# Peak detected, enter peak state
			maxPeakIndex = i
			maxPeakValue = s
			isHigh = True
			peakPlot += [0, 1]
			peakPos += [i, i]

# bits now holds the raw bitstream from the audio tape
# pos holds the sample index of the bits


print "Preparing plot..."
m = max(pulses)/3
peakPlot = [m+v*m for v in peakPlot]
bitPlot = [m+v*m for v in bits]
plt.plot(pulses,'g:')
plt.plot(peakPos,peakPlot,'b-')
plt.plot(pos,bitPlot,'k-*')
plt.show()

# bit stream starts with a lot of zeroes, find first 1-bit to synchronize
bits = bits[(bits.index(1)-1)%8:]

# Print out the bits and data
text = ''
for i in range(len(bits)/8):
	l = bits[i*8:i*8+8][::-1]
	c = eval('0b' + ''.join(map(str,l)))
	text += chr(c)
	#print l, hex(c), repr(chr(c))
print repr(text)

