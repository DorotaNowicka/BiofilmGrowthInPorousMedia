import os
import imageio
from config import SimInputData

def make_gif():
    print("I am making gif...")
    dirname = 'continuos_shear/ciekawy_przypadek/24/'
    image_folder = dirname
    gif_path = dirname + 'output.gif'
    duration = 0.5
    print(f"It should be {gif_path}")
    images = []

    # Get the list of image files in the folder
    image_files = sorted([f for f in os.listdir(image_folder) if f.endswith('.png')])

    for image_file in image_files:
        if "par" in str(image_file):
            continue
        if "flow" in str(image_file):
            continue
        print(image_file)
        image_path = os.path.join(image_folder, image_file)
        images.append(imageio.imread(image_path))

    # Save the GIF
    imageio.mimsave(gif_path, images, duration=duration)


make_gif()