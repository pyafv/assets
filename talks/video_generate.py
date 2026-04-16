import numpy as np
# import imageio
import imageio.v2 as imageio


def save_ani(image_list, file_name, fps=24):
    """ Saves the animation created.
    """
    frames = []
    for image_name in image_list:
        frames.append(imageio.imread(image_name))

    imageio.mimsave(file_name+".gif", frames, 'GIF', loop=0, fps=min(50, fps))
    # imageio.mimsave(file_name+".mp4", frames, 'MP4', fps=fps)


def save_mp4_high_quality(image_list, file_name, fps=60, bitrate="12M", subsampling="420"):
    """
    subsampling: "420" (widely compatible), "444" (sharper color; larger files)
    bitrate: e.g. "8M", "12M", "20M" — raise for larger frames or sharper needs
    """
    # Choose pixel format
    pix_fmt = "yuv420p" if subsampling == "420" else "yuv444p"

    with imageio.get_writer(
        file_name + ".mp4",
        fps=fps,
        codec="libx264",
        # Prevent hidden resizing to macro-block multiples:
        macro_block_size=None,
        # Ask for constant rate factor style quality control OR set bitrate
        # quality=10,             # alt: CRF-like (imageio maps this under the hood)
        bitrate=bitrate,          # explicit bitrate wins if both given
        pixelformat=pix_fmt,
        ffmpeg_params=[
            "-preset", "slow",   # better compression at given quality
            "-crf", "18"         # lower = better; combine with bitrate OR use this alone
        ],
    ) as writer:
        for im_path in image_list:
            frame = imageio.imread(im_path)  # uint8 HxWxC
            writer.append_data(frame)


def animation(tt, file_name):
    image_list = [f"frames/{t}.png" for t in tt]

    fps = 60
    save_ani(image_list, file_name, fps)
    save_mp4_high_quality(image_list, file_name, fps=fps, bitrate="12M", subsampling="420")


tt = np.arange(0, 15_000, 50)
animation(tt, "test")
print("Done!")

