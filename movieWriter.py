import os
import pathlib
import cv2
import numpy as np
from PIL import Image
import glob
import time
 
def generatemovie(fps,new_path):

    img_array = []
    #path = os.path.join(r'/gpfs/home/nmanzella/PWA/QuEP')

    old_width = 4800
    old_height = 3000
    dim_scale = 1.0
    new_width = int(old_width * dim_scale)
    new_height = int(old_height * dim_scale)
    print(f"New dimensions for video: {new_width}x{new_height}px")

    #fps = fps
    print(f"Running at {fps} frames per second")

    image_folder = new_path
    dim_scale_ = "{:03.2f}".format(dim_scale).replace(".","-")
    video_name = f'00-progression-movie__{dim_scale_}res__{fps}fps.mp4'
    os.chdir(image_folder)

    images = [img for img in os.listdir(image_folder)
            if img.endswith(".jpg") or
                img.endswith(".jpeg") or
                img.endswith("png")]

    images = sorted(images)
    
    print(f"\nUsing {len(images)} images from: {image_folder}")
    # Array images should only consider
    # the image files ignoring others if any
    #print(f"\nImages included in movie: {images}") 

    img = cv2.imread(os.path.join(image_folder, images[0]))
    frame = cv2.resize(img, dsize=(new_width,new_height), interpolation=cv2.INTER_CUBIC)

    # setting the frame width, height width
    # the width, height of first image
    height, width, layers = frame.shape

    fourcc = cv2.VideoWriter_fourcc(*"h264")
    #0x00000021
    video = cv2.VideoWriter(video_name, 0x31637661, fps, (width, height)) 

    start_time_video = time.time()
    t = time.localtime()
    curr_time = time.strftime("%H:%M:%S", t)
    print("\nMovie generation - START TIME: ", curr_time)

    # Appending the images to the video one by one
    for image in images: 
        img = cv2.imread(os.path.join(image_folder, image))
        frame = cv2.resize(img, dsize=(new_width,new_height), interpolation=cv2.INTER_CUBIC)
        video.write(frame) 

    # Deallocating memories taken for window creation
    cv2.destroyAllWindows() 
    video.release()  # releasing the video generated


    print(f"Movie generated! File stored at {os.path.join(image_folder,video_name)}")
    t_video_end = time.localtime()
    curr_time_video_end = time.strftime("%H:%M:%S", t_video_end)
    print("Movie generation - END TIME: ", curr_time_video_end)
    print("Movie generation - DURATION: ", (time.time() - start_time_video)/60, " min\n")

