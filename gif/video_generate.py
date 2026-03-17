import numpy as np
# import imageio
import imageio.v2 as imageio


def save_ani(image_list, file_name, fps=5):
    """ Saves the animation created.
    """
    frames = []
    for image_name in image_list:
        frames.append(imageio.imread(image_name))

    imageio.mimsave(file_name+".gif", frames, 'GIF', loop=0, fps=min(50, fps))
    # imageio.mimsave(file_name+".mp4", frames, 'MP4', fps=fps)


def animation(tt, file_name):
    image_list = []
    for t in tt:
        image_list.append("frames/%d.png" % t)

    fps = 144  # 60
    save_ani(image_list, file_name, fps)


tt = np.arange(2000, 7000, 50)
animation(tt, "test")
print("Done!")
