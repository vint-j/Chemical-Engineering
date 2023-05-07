import numpy as np
import matplotlib
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt
from numba import njit, prange
from matplotlib.widgets import Slider

@njit
def mandelbrot(c, max_iter):
    z = 0
    n = 0
    while abs(z) <= 2 and n < max_iter:
        z = z*z + c
        n += 1
    return n

@njit(parallel=True)
def mandelbrot_set(xmin, xmax, ymin, ymax, width, height, max_iter, output):
    r1 = np.linspace(xmin, xmax, width)
    r2 = np.linspace(ymin, ymax, height)
    for i in prange(height):
        for j in range(width):
            output[i, j] = mandelbrot(complex(r1[j], r2[i]), max_iter)
    return r1, r2

def tie_dye_color_map():
    return plt.cm.hsv

def display_mandelbrot(xmin, xmax, ymin, ymax, width=800, height=800, max_iter=2000):
    dpi = 80
    img_width = width
    img_height = height

    output = np.empty((height, width), dtype=np.int64)
    r1, r2 = mandelbrot_set(xmin, xmax, ymin, ymax, img_width, img_height, max_iter, output)
    output_log = np.log(1 + output)

    fig, ax = plt.subplots(figsize=(img_width/dpi, img_height/dpi), dpi=dpi)
    plt.subplots_adjust(bottom=0.2)

    ax.set_xticks(np.arange(0, img_width, 3*dpi))
    ax.set_xticklabels(xmin + (xmax - xmin)*np.arange(0, img_width, 3*dpi)/img_width)
    ax.set_yticks(np.arange(0, img_height, 3*dpi))
    ax.set_yticklabels(ymin + (ymax - ymin)*np.arange(0, img_height, 3*dpi)/img_width)

    ax.set_xlabel("Re(c)")
    ax.set_ylabel("Im(c)")
    ax.set_title("Mandelbrot Set with Tie-dye Colors")

    norm = plt.Normalize(output_log.min(), output_log.max())
    colors = tie_dye_color_map()(norm(output_log))
    img = ax.imshow(colors, extent=[xmin, xmax, ymin, ymax])

    ax_color = plt.axes([0.25, 0.1, 0.65, 0.03])
    slider_color = Slider(ax_color, 'Color Shift', 0, 1, valinit=0, valstep=0.01)

    def update(val):
        new_colormap = plt.cm.hsv(np.mod(np.arange(256) + val * 256, 256) / 256)
        new_colormap = plt.cm.colors.ListedColormap(new_colormap)
        new_colors = new_colormap(norm(output_log))
        img.set_data(new_colors)
        fig.canvas.draw()

    slider_color.on_changed(update)



    plt.show()

xmin, xmax, ymin, ymax = -2, 1, -1.5, 1.5
display_mandelbrot(xmin, xmax, ymin, ymax)
