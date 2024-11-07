import imageio
import os

# Directory containing the images.
image_dir = "plots/overnight3d"
# Output video file
output_file = "plots/video_overnight3d_temperature_slices.mp4"

# Get a list of all PNG files in the directory.
images = [img for img in os.listdir(image_dir) if img.endswith("png")]
# Sort the images by filename.
images.sort()

# Create a video writer object
with imageio.get_writer(output_file, fps=10) as writer:
    for filename in images:
        # Read the image
        image = imageio.imread(os.path.join(image_dir, filename))
        # Append the image to the video
        writer.append_data(image)

print(f"Video saved as {output_file}")