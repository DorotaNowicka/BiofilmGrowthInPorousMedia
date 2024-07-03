import os
import imageio
from config import SimInputData
import gc

def make_gif(sid: SimInputData):
    print("I am making gif...")
    image_folder = sid.dirname
    gif_path = sid.dirname + '/output.gif'
    duration = sid.graph_duration
    print(f"It should be {gif_path}")
    images = []

    # Get the list of image files in the folder
    image_files = sorted([f for f in os.listdir(image_folder) if f.endswith('.png')])

    for image_file in image_files:
        image_path = os.path.join(image_folder, image_file)
        images.append(imageio.imread(image_path))
        gc.collect()

    # Save the GIF
    imageio.mimsave(gif_path, images, duration=duration)


