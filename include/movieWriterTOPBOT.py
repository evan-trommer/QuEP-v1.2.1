import os
import pathlib
import cv2
import numpy as np
from PIL import Image
import glob
import time
 
def generatemovieTOPBOT(fps,image_folder_1,image_folder_2):

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

    
    dim_scale_ = "{:03.2f}".format(dim_scale).replace(".","-")
    video_name = f'00-progression-movie-TOPBOT__{dim_scale_}res__{fps}fps.mp4'
    
    # GET FIRST IMAGES
    os.chdir(image_folder_1)

    images1 = [img1 for img1 in os.listdir(image_folder_1)
            if img1.endswith(".jpg") or
                img1.endswith(".jpeg") or
                img1.endswith("png")]

    images1 = sorted(images1)
    
    # GET SECOND IMAGES
    os.chdir(image_folder_2)

    images2 = [img2 for img2 in os.listdir(image_folder_2)
            if img2.endswith(".jpg") or
                img2.endswith(".jpeg") or
                img2.endswith("png")]

    images2 = sorted(images2)

    print(f"\nUsing {len(images1)} images from: {image_folder_1}")
    print(f"\nUsing {len(images2)} images from: {image_folder_2}")

    # Array images should only consider
    # the image files ignoring others if any

    img1 = cv2.imread(os.path.join(image_folder_1, images1[0]))
    frame = cv2.resize(img1, dsize=(new_width,new_height), interpolation=cv2.INTER_CUBIC)

    # setting the frame width, height width
    # the width, height of first image
    height, width, layers = img1.shape

    fourcc = cv2.VideoWriter_fourcc(*"h264")
    #0x00000021
    video = cv2.VideoWriter(video_name, 0x31637661, fps, (width, (2*height-250-585))) 

    start_time_video = time.time()
    t = time.localtime()
    curr_time = time.strftime("%H:%M:%S", t)
    print("\nMovie generation - START TIME: ", curr_time)

    # Appending the images to the video one by one
    for i in range(0,len(images1)):
        filenumber = "{:05.1f}".format(i).replace(".","-") 
        img1 = cv2.imread(os.path.join(image_folder_1, f'progression-x-{filenumber}mm.png'))
        img2 = cv2.imread(os.path.join(image_folder_2, f'progression-x-{filenumber}mm.png'))
        cropped_img1 = img1[:2415,:]
        cropped_img2 = img2[250:,:]

        vis = np.concatenate((cropped_img1, cropped_img2), axis=0)
        frame = cv2.resize(vis, dsize=(new_width,new_height), interpolation=cv2.INTER_CUBIC)

        video.write(vis) 

    # Deallocating memories taken for window creation
    cv2.destroyAllWindows() 
    video.release()  # releasing the video generated


    print(f"Movie generated! File stored at {os.path.join(image_folder_2,video_name)}")
    t_video_end = time.localtime()
    curr_time_video_end = time.strftime("%H:%M:%S", t_video_end)
    print("Movie generation - END TIME: ", curr_time_video_end)
    print("Movie generation - END TIME: ", (time.time() - start_time_video)/60, " min\n")

