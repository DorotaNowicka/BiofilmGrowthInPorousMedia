import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import Normalize, LogNorm
from matplotlib.cm import viridis, summer, cool, autumn, winter


n = 10000
x = np.arange(1, n + 1)
y = np.ones(n)

fig, ax = plt.subplots(figsize=(6, 1.35))

# Create a colormap
cmap = viridis
autumn_reversed = autumn.reversed()
norm_under = Normalize(vmin=1, vmax=4500)
norm_over = Normalize(vmin=4501, vmax=9000)


colors = np.zeros((n, 4))



for i in range(n):
    if 1 <= x[i] <= 4500:
        colors[i] = viridis(norm_under(x[i]))
    elif 4500 < x[i] <= 9000:
        colors[i] = autumn_reversed(norm_over(x[i]))
    elif 9000 < x[i] <= 10000:
        colors[i] = [1, 0, 0, 1]  # Red color in RGBA


# Plot the image
ax.imshow(np.tile(colors, (1, 1, 1)), aspect='auto', extent=[1, n, 0, 1])

ax.text(4500, -0.2, "$F_{max}$", ha='center', va='top', fontsize=16)
ax.text(9000, -0.2, "$2F_{max}$", ha='center', va='top', fontsize=16)
ax.set_xticks([4500, 9000])
ax.tick_params(axis='x', which='major', direction='inout', top=False, bottom=True)
ax.tick_params(axis='x', labelbottom=False)

# Remove axes and labels
# ax.set_xticks([])
ax.set_yticks([])
# ax.set_xlabel('Shear force - color scale', fontsize=14)
# ax.set_xlabel('')

ax.set_ylabel('')
fig.tight_layout()

# Remove the border
for spine in ax.spines.values():
    spine.set_visible(False)

plt.savefig("legenda_shear.png")

n = 100
x = np.arange(1, n + 1)
y = np.ones(n)

fig, ax = plt.subplots(figsize=(6, 1.5))

# Initialize the colors array
colors = np.zeros((n, 4))

# Apply different colormaps for different range

norm_summer = Normalize(vmin=1, vmax=50)
norm_cool = Normalize(vmin=50, vmax=95)

for i in range(n):
    if 1 <= x[i] <= 50:
        colors[i] = summer(norm_summer(x[i]))
    elif 50 < x[i] <= 99:
        colors[i] = cool(norm_cool(x[i]))
    elif 99 < x[i] <= 100:
        colors[i] = [1, 0, 0, 1]  # Red color in RGBA

# Plot the image
ax.imshow(np.tile(colors, (1, 1, 1)), aspect='auto', extent=[1, n, 0, 1])

ax.text(99, 1.1, "$B_{max}$", ha='center', va='bottom', fontsize=16)
ax.text(99, 1.1, "$|$", ha='center', va='top', fontsize=7)

ax.text(50.5, 1.1, "$B_{spread}$", ha='center', va='bottom', fontsize=16)
ax.text(50.5, 1.1, "$|$", ha='center', va='top', fontsize=7)

# Remove axes and labels
ax.set_yticks([])
# ax.set_xlabel('Percent of diameter with bacteria - color scale', fontsize=14)
ax.set_xlabel('')
ax.set_ylabel('')

# Draw dashed vertical line at y=50
# ax.axvline(x=50, linestyle='--', color='black')

# Add text "bacteria spread" to the right of the color bar
# ax.text(102, 0.5, 'bacteria spread', va='center', fontsize=12)

# Remove the border
for spine in ax.spines.values():
    spine.set_visible(False)

fig.tight_layout()


plt.savefig("legenda_bacteria.png")
