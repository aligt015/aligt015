import numpy as np
import matplotlib.pyplot as plt

#############################################################################
# Question 1
#############################################################################

xb1, xb2 = 1.05, 2.13 #b_1 and b_2 airmass
xb = xb2 - xb1
V = np.array([ 12.01, 12.44, 12.19, 12.89])#v_column
BV= np.array([-0.07, 0.36, 0.69, 1.15])#b-v column

b1 = np.array([ 9.853, 10.693, 10.759, 11.898])#b_1 column
b2 = np.array([ 10.687, 11.479, 11.462, 12.547])#b_2 column

# extinction coeffs
kb = np.polyfit(BV, (b2 - b1)/xb, 1)

# model fuction that will be used to fit the polynomial
def peakfunc(x, m, b):
    return m*x + b

# plot the best-fit and data
plt.figure(figsize=(10,10))
plt.subplot(221)
plt.plot(BV, (b2 - b1)/xb, 'ro', label='B-band')
xdense1 = np.linspace(-0.5,1.5,25)
b_label1 = f'y = {kb[0]:.4f}x + {kb[1]:.4f}'
modelpoints = peakfunc(xdense1, kb[0], kb[1])
plt.plot(xdense1, modelpoints, 'b-', label = b_label1)
plt.legend()

xv1, xv2 = 1.10, 2.48 #v_1 and v_2 airmass
xv = xv2 - xv1

v1 = np.array([8.778, 9.160, 8.873, 9.522])#v_1 column
v2 = np.array([9.427, 9.739, 9.425, 10.001])#v_2 column

# extinction coeffs
kv = np.polyfit(BV, (v2 - v1)/xv, 1)
print('kv is: ', kv)

# plot the best-fit and data
plt.subplot(222)
plt.plot(BV, (v2 - v1)/xv, 'ro', label='V-band')
xdense2 = np.linspace(-0.5, 1.5, 25)
b_label2 = f'y = {kv[0]:.4f}x + {kv[1]:.4f}'
modelpoints = peakfunc(xdense2, kv[0], kv[1])
plt.plot(xdense2, modelpoints, 'b-', label = b_label2)
plt.legend()

#############################################################################
# Question 2
#############################################################################

#extinction corrected magnitudes
v1_c = v1 - (kv[1] + kv[0]*BV)*xv1
v2_c = v2 - (kv[1] + kv[0]*BV)*xv2
b1_c = b1 - (kb[1] + kb[0]*BV)*xb1
b2_c = b2 - (kb[1] + kb[0]*BV)*xb2

st_tr_v = np.stack((v1_c, v2_c))
st_tr_b = np.stack((b1_c, b2_c))
## standard transformation
BV2 = np.array([-0.07, 0.36, 0.69, 1.15, -0.07, 0.36, 0.69, 1.15])
b_st_array = BV + V - st_tr_b
b_st = np.concatenate([b_st_array[0], b_st_array[1]])
v_st_array = V - st_tr_v
v_st = np.concatenate([v_st_array[0], v_st_array[1]])

ab1, ab0 = np.polyfit(BV2, b_st, 1)
av1, av0 = np.polyfit(BV2, v_st, 1)

print('yo is: ', ab1, ab0, av1, av0)

plt.subplot(223)
plt.plot(BV2, b_st, 'ro', label = 'B-band')
b_label3 = f'y = {ab1:.4f}x + {ab0:.4f}'
plt.plot(xdense2, np.polyval([ab1, ab0], xdense2),'b-', label = b_label3)
plt.legend()

plt.subplot(224)
plt.plot(BV2, v_st, 'ro', label = 'V-band')
b_label4 = f'y = {av1:.4f}x + {av0:.4f}'
plt.plot(xdense2, np.polyval([av1, av0], xdense2),'b-', label = b_label4)
plt.legend()

#############################################################################
# Question 3
#############################################################################

#print(kb[0], kb[1], kv[0], kv[1], ab1, ab0, av1, av0)

#I will use the formula given to us in the lecture notes
b = 10.899
v = 9.850
airmass = 1.5

mb = b - (kb[1] + kb[0]*(b - v))*airmass
mv = v - (kv[1] + kv[0]*(b - v))*airmass

b_value = mb + ab0
v_value = mv + av0

#############################################################################
#  Now a little bit of algebra is needed to solve for (B-V)
#  B = b_value + ab1(B-V)
#  V = v_value + av1(B-V)
#  Now subrtracting both sides
#  B-V = b_value + ab1(B-V)-(v_value + av1(B-V))
#  therefore --> B-V = b_value - v_value + (B-V)(ab1 - av1)
#  therefore --> (B-V)(1 - (ab1-av1)) = b_value - v_value
#  finally ---> B-V = -(b_value - v_value)/(1 - (ab1-av1))
#############################################################################

BV = -(b_value - v_value)/(1 - (ab1-av1))
print('B-V for question 3 is: ', BV)

#Now since we know (B-V) we can just plug into our equation and solve for V
V = v - (kv[1] + kv[0]*BV)*airmass
print('V for questin 3 is: ', V)






